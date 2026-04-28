# 30d. Target-trial HR by month of LTFU, with 30-day grace-period landmark
# ==============================================================================
# Companion to 30c. Addresses classification immortal-time bias from the
# Brazilian programmatic 30-day rule for declaring LTFU:
#
#   * A patient is only classified as LTFU after >=30 days of no contact.
#   * A patient who dies within 30 days of last contact is re-classified
#     as 'Obito on treatment' (ascertained as Non-LTFU in our cohort)
#     rather than as LTFU.
#   * Therefore, the LTFU arm is implicitly conditioned on surviving the
#     30 days following true disengagement; the on-treatment comparator
#     is not.
#
# Fix: apply a symmetric 30-day landmark in BOTH arms.
#   * Eligibility at trial month m: alive at month m + 30 days.
#   * Time-at-risk origin: month m + 30 days.
#   * Comparator handling during grace window: lenient (survival only;
#     LTFU-status switches inside the grace window are not used to
#     exclude — the trial is "started LTFU at month m vs. didn't start
#     LTFU at month m"; switches afterwards are handled implicitly
#     by the cause-specific Cox in the same way as in 30c).
#
# Outputs:
#   target_trial_grace_mi_early_late_array.csv    (parallels 30c)
#
# This script does NOT overwrite the 30/30b/30c output files. It is
# intended as the going-forward primary analysis.
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
cat(sprintf("[30d] %d imputed datasets; grace = %d d (%.4f yr)\n",
            M, GRACE_DAYS, GRACE_YRS))

prepare_imp <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date_start <- as.Date(d$best_start)
  d$date_end   <- as.Date(d$end_date)
  d$tx_duration_yrs <- as.numeric(d$date_end - d$date_start) / 365.25
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

# Build one trial for a given imputation and month, with grace landmark.
# Origin for time-at-risk is month m + grace; events in the grace window
# are excluded by construction (eligibility filter).
build_trial_grace <- function(d, m) {
  start_yrs   <- (m - 1) * 30 / 365.25       # nominal trial start (LTFU window opens)
  end_yrs     <-  m      * 30 / 365.25       # nominal trial start (LTFU window closes)
  origin_yrs  <- start_yrs + GRACE_YRS       # post-grace origin for time-at-risk

  d |>
    # Symmetric landmark: must be alive at origin (month m + grace)
    dplyr::filter(time_d_tx > origin_yrs) |>
    # LTFU eligibility: tx_duration must reach the LTFU window for LTFU arm.
    # Non-LTFU arm: implicitly any tx_duration is fine.
    dplyr::mutate(
      eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
    ) |>
    dplyr::filter(eligible == 1) |>
    dplyr::mutate(
      # Exposed: LTFU happened in the [start, end) window.
      # Note: people whose LTFU date falls in [start, origin) survived the
      # grace and remain exposed; people whose LTFU date falls in [origin, end)
      # also remain exposed.
      expose = ifelse(itt_group == "Loss to follow-up" &
                        tx_duration_yrs >= start_yrs &
                        tx_duration_yrs <  end_yrs, 1, 0),
      event_d_num = as.numeric(as.character(event_d)),
      time_raw = time_d_tx - origin_yrs        # time at risk from grace-shifted origin
    )
}

# Outcome configurator: same family as 30c but computed from grace-shifted time.
# Windows are interpreted relative to the post-grace origin.
prep_outcome <- function(tr, model, cap) {
  t <- tr$time_raw
  e <- tr$event_d_num
  if (model == "overall") {
    time_out  <- pmin(t, cap)
    event_out <- ifelse(t > cap, 0, e)
  } else if (model == "early") {
    # Early window: 6 months from grace-shifted origin
    time_out  <- pmin(t, 0.5)
    event_out <- ifelse(t <= 0.5 & e == 1, 1, 0)
  } else if (model == "late") {
    # Late window: from 6 mo to cap, post-grace origin.
    # No additional landmark; everyone contributes follow-up from origin,
    # but only events in (0.5, cap] count.
    time_out  <- pmin(t, cap)
    event_out <- ifelse(t > 0.5 & t <= cap & e == 1, 1, 0)
  } else stop("unknown model")
  tr$time_out  <- time_out
  tr$event_out <- event_out
  tr[tr$time_out > 0, , drop = FALSE]
}

# Pooled Cox across imputations for a (month, model, cap) combination
fit_pooled <- function(m, model, cap) {
  fits <- lapply(imp_list, function(d) {
    tr <- build_trial_grace(d, m)
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
out_path <- file.path(ITT_RESULTS_DIR, "target_trial_grace_mi_early_late_array.csv")
write.csv(final_df, out_path, row.names = FALSE)
cat(sprintf("\n[30d] Wrote %d rows to %s\n", nrow(final_df), out_path))
print(final_df)

# ----------------------------------------------------------------------------
# Side-by-side comparison vs. the original 30c file (if present)
# ----------------------------------------------------------------------------
orig_path <- file.path(ITT_RESULTS_DIR, "target_trial_mi_early_late_array.csv")
if (file.exists(orig_path)) {
  orig <- read.csv(orig_path, stringsAsFactors = FALSE)
  comp <- merge(
    orig[, c("Trial_Month", "model", "cap", "HR", "CI_L", "CI_H")],
    final_df[, c("Trial_Month", "model", "cap", "HR", "CI_L", "CI_H")],
    by = c("Trial_Month", "model", "cap"),
    suffixes = c("_orig", "_grace")
  )
  comp$delta <- comp$HR_grace - comp$HR_orig
  comp$pct   <- 100 * comp$delta / comp$HR_orig
  cat("\n[30d] Original (no grace) vs grace 30d, side-by-side:\n")
  print(comp[order(comp$model, comp$cap, comp$Trial_Month),
             c("Trial_Month", "model", "cap",
               "HR_orig", "HR_grace", "delta", "pct")])
  comp_path <- file.path(ITT_RESULTS_DIR, "target_trial_grace_vs_original.csv")
  write.csv(comp, comp_path, row.names = FALSE)
  cat(sprintf("[30d] Comparison written to %s\n", comp_path))
}
