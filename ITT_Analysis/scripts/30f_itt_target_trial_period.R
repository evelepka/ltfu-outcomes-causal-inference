# 30f. Target-trial late mortality by calendar period of treatment start
# ==============================================================================
# Tests whether the LTFU effect on mortality differs across calendar time:
#   * Pre-COVID:  best_start in 2013-01-01 .. 2019-12-31
#   * Post-COVID: best_start in 2020-01-01 .. 2023-12-31
#
# Plus a continuous interaction (year of best_start) for trend visualization.
#
# Uses 30-day grace eligibility (mirrors 30d).
# Outputs the late-mortality (cap=2 yr) MI-pooled aHR per period, plus
# pooled across all six monthly trials (pooled trial-month indicator).
#
# Output: target_trial_period_grace_mi.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(mice)
  library(broom)
  library(lubridate)
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
cat(sprintf("[30f] %d imputed datasets; grace = %d d\n", M, GRACE_DAYS))

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
  d$year_start  <- year(d$date_start)
  d$period_2cat <- ifelse(d$year_start <= 2019, "Pre-COVID (2013-2019)",
                          "Post-COVID (2020-2023)")
  d$period_2cat <- factor(d$period_2cat, levels = c("Pre-COVID (2013-2019)",
                                                     "Post-COVID (2020-2023)"))
  d
}
imp_list <- lapply(imp_files, prepare_imp)

COVARS <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids",
            "diabetes", "alcohol", "drug_use", "incarcerated",
            "homelessness", "hosp_admission", "clinical_clean",
            "dot_status", "trial_month")

# ----------------------------------------------------------------------------
# Build pooled trials with grace eligibility (mirrors 32b's pooled-trial design)
# ----------------------------------------------------------------------------
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

apply_outcome_late <- function(tr, cap = 2) {
  t <- tr$time_raw
  e <- tr$event_d_num
  tr$time_out  <- pmin(t, cap)
  tr$event_out <- ifelse(t > 0.5 & t <= cap & e == 1, 1, 0)
  tr[tr$time_out > 0, , drop = FALSE]
}

# ----------------------------------------------------------------------------
# (A) Stratified estimates per period (Pre-/Post-COVID)
# ----------------------------------------------------------------------------
all_rows <- list()
for (lvl in levels(pooled_list[[1]]$period_2cat)) {
  cat(sprintf("\n[30f] Period: %s\n", lvl))
  for (cap in c(2, 5)) {
    fits <- lapply(pooled_list, function(pooled_i) {
      sub_i <- pooled_i[!is.na(pooled_i$period_2cat) & pooled_i$period_2cat == lvl, ]
      sub_i <- apply_outcome_late(sub_i, cap)
      if (nrow(sub_i) < 100 || sum(sub_i$event_out) < 5) return(NULL)
      rhs <- paste("expose +", paste(COVARS, collapse = " + "))
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
      analysis = "stratified",
      period   = lvl,
      cap      = cap,
      term     = "expose",
      HR       = ex$estimate,
      CI_L     = ex$conf.low,
      CI_H     = ex$conf.high,
      P_Value  = ex$p.value,
      N_imp    = length(fits)
    )
    cat(sprintf("  cap=%d  HR=%.2f (%.2f-%.2f)  p=%.3g\n",
                cap, ex$estimate, ex$conf.low, ex$conf.high, ex$p.value))
  }
}

# ----------------------------------------------------------------------------
# (B) Formal interaction test: expose * period_2cat in a single model
# ----------------------------------------------------------------------------
cat("\n[30f] Interaction test (expose * period_2cat) on the late, cap=2 model\n")
for (cap in c(2, 5)) {
  fits <- lapply(pooled_list, function(pooled_i) {
    sub_i <- apply_outcome_late(pooled_i, cap)
    if (nrow(sub_i) < 100) return(NULL)
    rhs_main <- paste("expose * period_2cat +", paste(COVARS, collapse = " + "))
    f <- as.formula(paste("Surv(time_out, event_out) ~", rhs_main))
    tryCatch(coxph(f, data = sub_i, cluster = sinan_clean),
             error = function(e) NULL)
  })
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) next
  pooled_fit <- tryCatch(
    summary(pool(as.mira(fits)), exponentiate = TRUE, conf.int = TRUE),
    error = function(e) NULL)
  if (is.null(pooled_fit)) next
  for (term_name in c("expose", "expose:period_2catPost-COVID (2020-2023)")) {
    ex <- pooled_fit[pooled_fit$term == term_name, ]
    if (nrow(ex) == 0) next
    all_rows[[length(all_rows) + 1]] <- data.frame(
      analysis = "interaction",
      period   = "ALL",
      cap      = cap,
      term     = term_name,
      HR       = ex$estimate,
      CI_L     = ex$conf.low,
      CI_H     = ex$conf.high,
      P_Value  = ex$p.value,
      N_imp    = length(fits)
    )
    cat(sprintf("  cap=%d  term=%-50s  HR=%.2f (%.2f-%.2f)  p=%.3g\n",
                cap, term_name, ex$estimate, ex$conf.low, ex$conf.high,
                ex$p.value))
  }
}

# ----------------------------------------------------------------------------
# (C) Continuous: year_start centered at 2018, plus expose:year_start
# ----------------------------------------------------------------------------
cat("\n[30f] Continuous interaction (expose * year_start_c) on late, cap=2\n")
for (cap in c(2, 5)) {
  fits <- lapply(pooled_list, function(pooled_i) {
    pooled_i$year_start_c <- pooled_i$year_start - 2018
    sub_i <- apply_outcome_late(pooled_i, cap)
    if (nrow(sub_i) < 100) return(NULL)
    rhs_main <- paste("expose * year_start_c +",
                      paste(COVARS, collapse = " + "))
    f <- as.formula(paste("Surv(time_out, event_out) ~", rhs_main))
    tryCatch(coxph(f, data = sub_i, cluster = sinan_clean),
             error = function(e) NULL)
  })
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) next
  pooled_fit <- tryCatch(
    summary(pool(as.mira(fits)), exponentiate = TRUE, conf.int = TRUE),
    error = function(e) NULL)
  if (is.null(pooled_fit)) next
  for (term_name in c("expose", "year_start_c", "expose:year_start_c")) {
    ex <- pooled_fit[pooled_fit$term == term_name, ]
    if (nrow(ex) == 0) next
    all_rows[[length(all_rows) + 1]] <- data.frame(
      analysis = "continuous",
      period   = "ALL",
      cap      = cap,
      term     = term_name,
      HR       = ex$estimate,
      CI_L     = ex$conf.low,
      CI_H     = ex$conf.high,
      P_Value  = ex$p.value,
      N_imp    = length(fits)
    )
    cat(sprintf("  cap=%d  term=%-25s  HR=%.3f (%.3f-%.3f)  p=%.3g\n",
                cap, term_name, ex$estimate, ex$conf.low, ex$conf.high,
                ex$p.value))
  }
}

final_df <- bind_rows(all_rows)
out_path <- file.path(ITT_RESULTS_DIR, "target_trial_period_grace_mi.csv")
write.csv(final_df, out_path, row.names = FALSE)
cat(sprintf("\n[30f] Wrote %d rows to %s\n", nrow(final_df), out_path))
print(final_df)

# ----------------------------------------------------------------------------
# Counts
# ----------------------------------------------------------------------------
cat("\n[30f] Counts by period (after grace eligibility, late events only):\n")
pi <- pooled_list[[1]]
for (lvl in levels(pi$period_2cat)) {
  sub <- pi[!is.na(pi$period_2cat) & pi$period_2cat == lvl, ]
  late <- apply_outcome_late(sub, 2)
  cat(sprintf("  %-25s  N_obs=%d  late events=%d  exposed=%d\n",
              lvl, nrow(late), sum(late$event_out), sum(late$expose)))
}
