# 32e. Subgroup target-trial under defn-B + grace
# ==============================================================================
# Mirrors 32d but uses defn-B (LTFU exposure shifted back 30 days).
# Output: target_trial_defnB_subgroups_mi.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(survival); library(mice); library(broom)
})

.here <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  frames <- sys.frames()
  for (f in rev(frames)) { of <- f$ofile; if (!is.null(of)) return(dirname(normalizePath(of))) }
  getwd()
}
source(file.path(.here(), "_paths.R"))

GRACE_YRS <- 30 / 365.25
SHIFT_YRS <- 30 / 365.25

imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
stopifnot(length(imp_files) > 0)
cat(sprintf("[32e] %d imputed datasets; defn-B + grace\n", length(imp_files)))

prepare_imp <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date_start <- as.Date(d$best_start); d$date_end <- as.Date(d$end_date)
  d$tx_duration_yrs <- as.numeric(d$date_end - d$date_start) / 365.25
  d$true_tx_duration_yrs <- pmax(d$tx_duration_yrs - SHIFT_YRS, 1/365.25)
  d$age_group <- factor(d$age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
  for (v in c("sex","race_clean","edu_clean","hiv_aids","diabetes",
              "alcohol","drug_use","incarcerated","homelessness",
              "hosp_admission","clinical_clean","dot_status")) {
    d[[v]] <- as.factor(d[[v]])
  }
  d
}
imp_list <- lapply(imp_files, prepare_imp)

COVARS <- c("age_group","sex","race_clean","edu_clean","hiv_aids",
            "diabetes","alcohol","drug_use","incarcerated",
            "homelessness","hosp_admission","clinical_clean",
            "dot_status","trial_month")

build_pooled <- function(d) {
  trial_list <- vector("list", 6)
  for (m in 1:6) {
    start_yrs <- (m-1) * 30 / 365.25
    end_yrs   <-  m    * 30 / 365.25
    origin_yrs <- start_yrs + GRACE_YRS
    tr <- d |>
      dplyr::filter(time_d_tx > origin_yrs) |>
      dplyr::mutate(eligible = ifelse(itt_group == "Non-LTFU" |
                                        true_tx_duration_yrs >= start_yrs, 1, 0)) |>
      dplyr::filter(eligible == 1) |>
      dplyr::mutate(
        expose = ifelse(itt_group == "Loss to follow-up" &
                          true_tx_duration_yrs >= start_yrs &
                          true_tx_duration_yrs <  end_yrs, 1, 0),
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
pooled_list <- lapply(imp_list, build_pooled)

apply_outcome <- function(tr, model, cap) {
  t <- tr$time_raw; e <- tr$event_d_num
  if (model == "overall") { tr$time_out <- pmin(t,cap); tr$event_out <- ifelse(t>cap,0,e) }
  else if (model == "early") { tr$time_out <- pmin(t,0.5); tr$event_out <- ifelse(t<=0.5 & e==1,1,0) }
  else if (model == "late")  { tr$time_out <- pmin(t,cap); tr$event_out <- ifelse(t>0.5 & t<=cap & e==1,1,0) }
  tr[tr$time_out > 0, , drop = FALSE]
}

SUBGROUPS <- c("age_group","sex","hiv_aids","homelessness")
CONFIGS <- list(
  list(model="overall", cap=2),
  list(model="overall", cap=5),
  list(model="early", cap=0.5),
  list(model="late", cap=2),
  list(model="late", cap=5)
)

all_rows <- list()
for (sg in SUBGROUPS) {
  cat(sprintf("\n[32e] Subgroup: %s\n", sg))
  base_covs <- setdiff(COVARS, sg)
  lvls <- sort(unique(as.character(pooled_list[[1]][[sg]])))
  lvls <- lvls[!is.na(lvls) & lvls != ""]
  for (lvl in lvls) for (cfg in CONFIGS) {
    fits <- lapply(pooled_list, function(p_i) {
      sub_i <- p_i[!is.na(p_i[[sg]]) & p_i[[sg]] == lvl, ]
      sub_i <- apply_outcome(sub_i, cfg$model, cfg$cap)
      if (nrow(sub_i) < 100 || sum(sub_i$event_out) < 5) return(NULL)
      rhs <- paste("expose +", paste(base_covs, collapse=" + "))
      f <- as.formula(paste("Surv(time_out, event_out) ~", rhs))
      tryCatch(coxph(f, data=sub_i, cluster=sinan_clean), error=function(e) NULL)
    })
    fits <- Filter(Negate(is.null), fits)
    if (length(fits) == 0) next
    pf <- tryCatch(summary(pool(as.mira(fits)), exponentiate=TRUE, conf.int=TRUE),
                   error=function(e) NULL)
    if (is.null(pf)) next
    ex <- pf[pf$term == "expose", ]; if (nrow(ex) == 0) next
    all_rows[[length(all_rows)+1]] <- data.frame(
      Subgroup=sg, Level=lvl, model=cfg$model, cap=cfg$cap,
      HR=ex$estimate, CI_L=ex$conf.low, CI_H=ex$conf.high,
      P_Value=ex$p.value, N_imp=length(fits)
    )
    cat(sprintf("  %-13s %-14s %-7s cap=%.1f  HR=%.2f (%.2f-%.2f)\n",
                sg, lvl, cfg$model, cfg$cap, ex$estimate, ex$conf.low, ex$conf.high))
  }
}

final_df <- bind_rows(all_rows)
out_path <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_subgroups_mi.csv")
write.csv(final_df, out_path, row.names = FALSE)
cat(sprintf("\n[32e] Wrote %d rows to %s\n", nrow(final_df), out_path))
