# 32c. Target-trial subgroup HRs for drug-resistance status, grace-period
# ==============================================================================
# Mirrors 32b but stratifies by `resistance_clean` and applies the 30-day
# grace-period eligibility (as in 30d).
#
# Levels:
#   - Sensitive
#   - Resistant (Any)        — TB-MR/TB-R per cleaning rule
#   - Not Evaluated          — no DST result; programmatic gap
#
# Output: target_trial_resistance_grace_mi.csv
#         (rows can slot into Figure 3C alongside age/sex/HIV/homelessness)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(mice)
  library(broom)
})

.here <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  frames <- sys.frames()
  for (f in rev(frames)) {
    of <- f$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  getwd()
}
source(file.path(.here(), "_paths.R"))

GRACE_DAYS <- 30
GRACE_YRS  <- GRACE_DAYS / 365.25

imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
stopifnot(length(imp_files) > 0)
M <- length(imp_files)
cat(sprintf("[32c] %d imputed datasets; grace = %d d\n", M, GRACE_DAYS))

prepare_imp <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date_start <- as.Date(d$best_start)
  d$date_end   <- as.Date(d$end_date)
  d$tx_duration_yrs <- as.numeric(d$date_end - d$date_start) / 365.25
  d$age_group <- factor(d$age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
  for (v in c("sex", "race_clean", "edu_clean", "hiv_aids", "diabetes",
              "alcohol", "drug_use", "incarcerated", "homelessness",
              "hosp_admission", "clinical_clean", "dot_status",
              "resistance_clean")) {
    if (v %in% names(d)) d[[v]] <- as.factor(d[[v]])
  }
  d
}
imp_list <- lapply(imp_files, prepare_imp)

COVARS <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids",
            "diabetes", "alcohol", "drug_use", "incarcerated",
            "homelessness", "hosp_admission", "clinical_clean",
            "dot_status", "trial_month")

# Build pooled-trials dataset with grace eligibility
build_pooled_grace <- function(d) {
  trial_list <- vector("list", 6)
  for (m in 1:6) {
    start_yrs   <- (m - 1) * 30 / 365.25
    end_yrs     <-  m      * 30 / 365.25
    origin_yrs  <- start_yrs + GRACE_YRS
    tr <- d |>
      dplyr::filter(time_d_tx > origin_yrs) |>
      dplyr::mutate(
        eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
      ) |>
      dplyr::filter(eligible == 1) |>
      dplyr::mutate(
        expose = ifelse(itt_group == "Loss to follow-up" &
                          tx_duration_yrs >= start_yrs &
                          tx_duration_yrs <  end_yrs, 1, 0),
        trial_month = paste0("Month_", m),
        event_d_num = as.numeric(as.character(event_d)),
        time_raw    = time_d_tx - origin_yrs
      )
    trial_list[[m]] <- tr
  }
  pooled <- bind_rows(trial_list)
  pooled$trial_month <- factor(pooled$trial_month, levels = paste0("Month_", 1:6))
  pooled
}

pooled_list <- lapply(imp_list, build_pooled_grace)

# Outcome configurator — same as 32b
apply_outcome <- function(tr, model, cap) {
  t <- tr$time_raw
  e <- tr$event_d_num
  if (model == "overall") {
    tr$time_out  <- pmin(t, cap)
    tr$event_out <- ifelse(t > cap, 0, e)
  } else if (model == "early") {
    tr$time_out  <- pmin(t, 0.5)
    tr$event_out <- ifelse(t <= 0.5 & e == 1, 1, 0)
  } else if (model == "late") {
    tr$time_out  <- pmin(t, cap)
    tr$event_out <- ifelse(t > 0.5 & t <= cap & e == 1, 1, 0)
  }
  tr[tr$time_out > 0, , drop = FALSE]
}

SUBGROUPS <- c("resistance_clean")
CONFIGS <- list(
  list(model = "overall", cap = 2),
  list(model = "early",   cap = 0.5),
  list(model = "late",    cap = 2),
  list(model = "late",    cap = 5)
)

all_rows <- list()
for (sg in SUBGROUPS) {
  cat(sprintf("\n[32c] Subgroup: %s\n", sg))
  base_covs <- setdiff(COVARS, sg)
  lvls <- sort(unique(as.character(pooled_list[[1]][[sg]])))
  lvls <- lvls[!is.na(lvls) & lvls != ""]
  cat(sprintf("  Levels: %s\n", paste(lvls, collapse = ", ")))
  for (lvl in lvls) {
    for (cfg in CONFIGS) {
      fits <- lapply(pooled_list, function(pooled_i) {
        sub_i <- pooled_i[!is.na(pooled_i[[sg]]) & pooled_i[[sg]] == lvl, ]
        sub_i <- apply_outcome(sub_i, cfg$model, cfg$cap)
        if (nrow(sub_i) < 100 || sum(sub_i$event_out) < 5) return(NULL)
        rhs <- paste("expose +", paste(base_covs, collapse = " + "))
        f <- as.formula(paste("Surv(time_out, event_out) ~", rhs))
        tryCatch(coxph(f, data = sub_i, cluster = sinan_clean),
                 error = function(e) NULL)
      })
      fits <- Filter(Negate(is.null), fits)
      if (length(fits) == 0) next
      pooled_fit <- tryCatch(
        summary(pool(as.mira(fits)), exponentiate = TRUE, conf.int = TRUE),
        error = function(e) NULL)
      if (is.null(pooled_fit)) next
      ex <- pooled_fit[pooled_fit$term == "expose", ]
      if (nrow(ex) == 0) next
      all_rows[[length(all_rows) + 1]] <- data.frame(
        Subgroup = sg, Level = lvl,
        model = cfg$model, cap = cfg$cap,
        HR = ex$estimate, CI_L = ex$conf.low, CI_H = ex$conf.high,
        P_Value = ex$p.value, N_imp = length(fits)
      )
      cat(sprintf("  %-17s %-7s cap=%.1f  HR=%.2f (%.2f-%.2f)  p=%.3g  N_imp=%d\n",
                  lvl, cfg$model, cfg$cap, ex$estimate, ex$conf.low, ex$conf.high,
                  ex$p.value, length(fits)))
    }
  }
}

final_df <- bind_rows(all_rows)
out_path <- file.path(ITT_RESULTS_DIR, "target_trial_resistance_grace_mi.csv")
write.csv(final_df, out_path, row.names = FALSE)
cat(sprintf("\n[32c] Wrote %d rows to %s\n", nrow(final_df), out_path))
print(final_df)

# Event counts and N per level for transparency (trial month=2 representative)
cat("\n[32c] Counts by resistance level (pooled trials, after grace eligibility):\n")
for (sg in SUBGROUPS) {
  pooled_i <- pooled_list[[1]]
  lvls <- sort(unique(as.character(pooled_i[[sg]])))
  lvls <- lvls[!is.na(lvls) & lvls != ""]
  for (lvl in lvls) {
    sub <- pooled_i[!is.na(pooled_i[[sg]]) & pooled_i[[sg]] == lvl, ]
    sub_late <- apply_outcome(sub, "late", 2)
    cat(sprintf("  %-17s  N_obs=%d  late events=%d  exposed=%d\n",
                lvl, nrow(sub_late), sum(sub_late$event_out),
                sum(sub_late$expose)))
  }
}
