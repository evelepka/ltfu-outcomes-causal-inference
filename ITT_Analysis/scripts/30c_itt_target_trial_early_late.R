# 30c. Target-trial HR by month of LTFU, split into Early / Late mortality
# ==============================================================================
# Replaces the single-HR output of 30b with a stacked array covering:
#   - Overall    (event_d, time = time_d_tx - start_yrs, capped at 2 or 5 yr)
#   - Early      (events in first 6 months from trial start; cap = 0.5 yr)
#   - Late       (events after 6 months from trial start; cap = 2 or 5 yr)
# Method: cause-specific Cox. For "Early", late deaths get censored at 0.5.
# For "Late", early deaths get censored at their actual death time (event=0)
# — no landmark restriction; everyone contributes follow-up from trial start.
# Anchor: time_d_tx (years from best_start / treatment initiation).
# Input:  5 miceforest-imputed datasets at ITT_Analysis/data/mi/imp_*.csv
# Output: ITT_Analysis/results/target_trial_mi_early_late_array.csv
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

imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
stopifnot(length(imp_files) > 0)
M <- length(imp_files)
cat(sprintf("[30c] %d imputed datasets\n", M))

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

# Build one trial for a given imputation and month
build_trial <- function(d, m) {
  start_yrs <- (m - 1) * 30 / 365.25
  end_yrs   <-  m      * 30 / 365.25
  d |>
    dplyr::filter(time_d_tx > start_yrs) |>
    dplyr::mutate(eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)) |>
    dplyr::filter(eligible == 1) |>
    dplyr::mutate(
      expose = ifelse(itt_group == "Loss to follow-up" &
                        tx_duration_yrs >= start_yrs &
                        tx_duration_yrs <  end_yrs, 1, 0),
      event_d_num = as.numeric(as.character(event_d)),
      time_raw = time_d_tx - start_yrs
    )
}

# Outcome configurator: for a specified (model, cap), return time + event columns
# model: "overall" | "early" | "late"
# cap: numeric (years)
prep_outcome <- function(tr, model, cap) {
  t <- tr$time_raw
  e <- tr$event_d_num
  if (model == "overall") {
    time_out <- pmin(t, cap)
    event_out <- ifelse(t > cap, 0, e)
  } else if (model == "early") {
    # Early window: (0, 0.5]
    time_out <- pmin(t, 0.5)
    event_out <- ifelse(t <= 0.5 & e == 1, 1, 0)
  } else if (model == "late") {
    # Late window: (0.5, cap]. No landmark: full follow-up, but only count late events.
    time_out <- pmin(t, cap)
    event_out <- ifelse(t > 0.5 & t <= cap & e == 1, 1, 0)
  } else stop("unknown model")
  tr$time_out <- time_out
  tr$event_out <- event_out
  tr[tr$time_out > 0, , drop = FALSE]
}

# Fit pooled Cox across the M imputations for a given trial/model/cap combo
fit_pooled <- function(m, model, cap) {
  fits <- lapply(imp_list, function(d) {
    tr <- build_trial(d, m)
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
    model = model, cap = cap,
    HR = ex$estimate, CI_L = ex$conf.low, CI_H = ex$conf.high,
    P_Value = ex$p.value, N_imp = length(fits)
  )
}

# Loop trials × outcome configurations
rows <- list()
CONFIGS <- list(
  list(model = "overall", cap = 2),
  list(model = "overall", cap = 5),
  list(model = "early",   cap = 0.5),
  list(model = "late",    cap = 2),
  list(model = "late",    cap = 5)
)
for (m in 1:6) {
  for (cfg in CONFIGS) {
    cat(sprintf("  trial=%d model=%-7s cap=%.1f\n", m, cfg$model, cfg$cap))
    r <- fit_pooled(m, cfg$model, cfg$cap)
    if (!is.null(r)) rows[[length(rows) + 1]] <- r
  }
}
final_df <- bind_rows(rows)
out_path <- file.path(ITT_RESULTS_DIR, "target_trial_mi_early_late_array.csv")
write.csv(final_df, out_path, row.names = FALSE)
cat(sprintf("\n[30c] Wrote %d rows to %s\n", nrow(final_df), out_path))
print(final_df)
