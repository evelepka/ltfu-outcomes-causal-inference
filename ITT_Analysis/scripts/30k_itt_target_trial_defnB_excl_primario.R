# 30k. Defn-B target trial — sensitivity excluding never-treated Primario
# ==============================================================================
# Per Sueli (2026-04-30): Abandono Primario cases without a recorded tx_start
# represent diagnosis-confirmed-but-treatment-never-initiated patients
# (closure rule: >=90 days post-notification). They are coded LTFU under the
# Brazilian programmatic definition but did not receive the intervention,
# so an LTFU-vs-on-treatment causal contrast assigns them an exposure they
# never had.
#
# This sensitivity refits 30h's defn-B grace target trial after excluding
# the 788 sinan_clean IDs in primario_no_txstart_ids.csv (built by
# tmp_primario_sanity.py from the raw TBweb pull).
#
# Output:
#   target_trial_defnB_excl_primario_no_tx.csv   (parallel structure to 30h)
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
SHIFT_YRS  <- 30 / 365.25

excl_path <- file.path(ITT_DATA_DIR, "primario_no_txstart_ids.csv")
stopifnot(file.exists(excl_path))
excl_ids <- read.csv(excl_path, stringsAsFactors = FALSE)$sinan_clean
cat(sprintf("[30k] Exclusion list: %d sinan_clean IDs\n", length(excl_ids)))

imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
stopifnot(length(imp_files) > 0)
M <- length(imp_files)
cat(sprintf("[30k] %d imputed datasets; defn-B shift = %d d; grace = %d d\n",
            M, 30, GRACE_DAYS))

prepare_imp <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  before <- nrow(d)
  d <- d[!d$sinan_clean %in% excl_ids, , drop = FALSE]
  cat(sprintf("  %s: %d -> %d (excluded %d)\n",
              basename(path), before, nrow(d), before - nrow(d)))
  d$date_start <- as.Date(d$best_start)
  d$date_end   <- as.Date(d$end_date)
  d$tx_duration_yrs <- as.numeric(d$date_end - d$date_start) / 365.25
  d$true_tx_duration_yrs <- pmax(d$tx_duration_yrs - SHIFT_YRS, 1/365.25)
  d$age_group <- factor(d$age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
  for (v in c("sex", "race_clean", "edu_clean", "hiv_aids", "diabetes",
              "alcohol", "drug_use", "incarcerated", "homelessness",
              "hosp_admission", "clinical_clean", "dot_status")) {
    d[[v]] <- as.factor(d[[v]])
  }
  d
}
imp_list <- lapply(imp_files, prepare_imp)

COVARS <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids",
            "diabetes", "alcohol", "drug_use", "incarcerated",
            "homelessness", "hosp_admission", "clinical_clean", "dot_status")
RHS <- paste("expose +", paste(COVARS, collapse = " + "))

build_trial_defnB <- function(d, m) {
  start_yrs   <- (m - 1) * 30 / 365.25
  end_yrs     <-  m      * 30 / 365.25
  origin_yrs  <- start_yrs + GRACE_YRS

  d |>
    dplyr::filter(time_d_tx > origin_yrs) |>
    dplyr::mutate(
      eligible = ifelse(itt_group == "Non-LTFU" |
                          true_tx_duration_yrs >= start_yrs, 1, 0)
    ) |>
    dplyr::filter(eligible == 1) |>
    dplyr::mutate(
      expose = ifelse(itt_group == "Loss to follow-up" &
                        true_tx_duration_yrs >= start_yrs &
                        true_tx_duration_yrs <  end_yrs, 1, 0),
      event_d_num = as.numeric(as.character(event_d)),
      time_raw    = time_d_tx - origin_yrs
    )
}

prep_outcome <- function(tr, model, cap) {
  t <- tr$time_raw
  e <- tr$event_d_num
  if (model == "overall") {
    time_out  <- pmin(t, cap)
    event_out <- ifelse(t > cap, 0, e)
  } else if (model == "early") {
    time_out  <- pmin(t, 0.5)
    event_out <- ifelse(t <= 0.5 & e == 1, 1, 0)
  } else if (model == "late") {
    time_out  <- pmin(t, cap)
    event_out <- ifelse(t > 0.5 & t <= cap & e == 1, 1, 0)
  } else stop("unknown model")
  tr$time_out  <- time_out
  tr$event_out <- event_out
  tr[tr$time_out > 0, , drop = FALSE]
}

fit_pooled <- function(m, model, cap) {
  fits <- lapply(imp_list, function(d) {
    tr <- build_trial_defnB(d, m)
    if (nrow(tr) < 50) return(NULL)
    tr <- prep_outcome(tr, model, cap)
    if (sum(tr$event_out) < 5) return(NULL)
    tryCatch(
      coxph(as.formula(paste("Surv(time_out, event_out) ~", RHS)), data = tr),
      error = function(e) NULL
    )
  })
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) return(NULL)
  pooled <- tryCatch(
    summary(pool(as.mira(fits)), exponentiate = TRUE, conf.int = TRUE),
    error = function(e) NULL)
  if (is.null(pooled)) return(NULL)
  ex <- pooled[pooled$term == "expose", ]
  if (nrow(ex) == 0) return(NULL)
  data.frame(
    Trial_Month = paste0("Month_", m),
    model       = model,
    cap         = cap,
    HR          = ex$estimate,
    CI_L        = ex$conf.low,
    CI_H        = ex$conf.high,
    P_Value     = ex$p.value,
    N_imp       = length(fits)
  )
}

CONFIGS <- list(
  list(model = "overall", cap = 2),
  list(model = "overall", cap = 5),
  list(model = "early",   cap = 0.5),
  list(model = "late",    cap = 2),
  list(model = "late",    cap = 5)
)

rows <- list()
for (m in 1:6) {
  for (cfg in CONFIGS) {
    cat(sprintf("  trial=%d model=%-7s cap=%.1f\n", m, cfg$model, cfg$cap))
    r <- fit_pooled(m, cfg$model, cfg$cap)
    if (!is.null(r)) rows[[length(rows) + 1]] <- r
  }
}

final_df <- bind_rows(rows)
out_path <- file.path(ITT_RESULTS_DIR,
                     "target_trial_defnB_excl_primario_no_tx.csv")
write.csv(final_df, out_path, row.names = FALSE)
cat(sprintf("\n[30k] Wrote %d rows to %s\n", nrow(final_df), out_path))
print(final_df)

# Side-by-side vs the primary defn-B (with no exclusion).
prim <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_mi_early_late_array.csv")
if (file.exists(prim)) {
  p <- read.csv(prim, stringsAsFactors = FALSE)
  comp <- merge(
    final_df[, c("Trial_Month", "model", "cap", "HR", "CI_L", "CI_H")],
    p[,        c("Trial_Month", "model", "cap", "HR", "CI_L", "CI_H")],
    by = c("Trial_Month", "model", "cap"),
    suffixes = c("_excl", "_primary")
  )
  comp$pct_diff <- 100 * (comp$HR_excl - comp$HR_primary) / comp$HR_primary
  comp_path <- file.path(ITT_RESULTS_DIR,
                        "target_trial_defnB_excl_vs_primary.csv")
  write.csv(comp, comp_path, row.names = FALSE)
  cat(sprintf("\n[30k] Comparison written to %s\n", comp_path))
  cat("\n[30k] Late, cap=2 yr (the headline window):\n")
  print(comp[comp$model == "late" & comp$cap == 2,
             c("Trial_Month", "HR_primary", "HR_excl", "pct_diff")])
}
