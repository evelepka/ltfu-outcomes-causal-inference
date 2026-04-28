# 30h. Target-trial under definition B (LTFU date shifted back 30 days)
# ==============================================================================
# Primary analysis going forward, per Sueli's clarification that TBweb's
# data_de_encerramento for Abandono cases is the case-CLOSURE date (= ~30+
# days after the actual last visit/medication), not the date of last visit.
#
# We therefore shift the LTFU exposure point back by 30 days:
#
#   true_tx_duration_yrs = max(tx_duration_yrs - 30/365.25, 1/365.25)
#
# (Clamped to a minimum of 1 day to handle the ~4% of `Abandono` cases that
# violate strict defn-B with tx_days < 30 — likely data-entry edge cases.)
#
# Trial-month interpretation under defn B is the actual disengagement
# month rather than the classification month (one month earlier than under
# defn A).
#
# Grace-period landmark is still applied (mirrors 30d): both arms must be
# alive at month_m + 30d; time-at-risk origin shifted to month_m + 30d.
#
# Output:
#   target_trial_defnB_mi_early_late_array.csv  (parallels 30d)
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
SHIFT_YRS  <- 30 / 365.25   # defn-B shift back

imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
stopifnot(length(imp_files) > 0)
M <- length(imp_files)
cat(sprintf("[30h] %d imputed datasets; defn-B shift = %d d; grace = %d d\n",
            M, 30, GRACE_DAYS))

prepare_imp <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date_start <- as.Date(d$best_start)
  d$date_end   <- as.Date(d$end_date)
  d$tx_duration_yrs <- as.numeric(d$date_end - d$date_start) / 365.25
  # defn-B: actual disengagement is ~30 days before recorded end_date.
  # Clamp to 1 day to avoid negatives for the ~4% violators.
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

# Build trial under defn B + grace landmark.
build_trial_defnB <- function(d, m) {
  start_yrs   <- (m - 1) * 30 / 365.25      # actual-disengagement window m
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
out_path <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_mi_early_late_array.csv")
write.csv(final_df, out_path, row.names = FALSE)
cat(sprintf("\n[30h] Wrote %d rows to %s\n", nrow(final_df), out_path))
print(final_df)

# Side-by-side vs defn-A grace and original (no-grace) for context.
ga <- file.path(ITT_RESULTS_DIR, "target_trial_grace_mi_early_late_array.csv")
or <- file.path(ITT_RESULTS_DIR, "target_trial_mi_early_late_array.csv")
if (file.exists(ga) && file.exists(or)) {
  defA <- read.csv(ga, stringsAsFactors = FALSE)
  orig <- read.csv(or, stringsAsFactors = FALSE)
  comp <- final_df[, c("Trial_Month", "model", "cap", "HR", "CI_L", "CI_H")]
  names(comp)[4:6] <- c("HR_defnB", "L_defnB", "H_defnB")
  comp <- merge(comp,
                defA[, c("Trial_Month", "model", "cap", "HR", "CI_L", "CI_H")],
                by = c("Trial_Month", "model", "cap"))
  names(comp)[7:9] <- c("HR_defnA_grace", "L_defnA_grace", "H_defnA_grace")
  comp <- merge(comp,
                orig[, c("Trial_Month", "model", "cap", "HR", "CI_L", "CI_H")],
                by = c("Trial_Month", "model", "cap"))
  names(comp)[10:12] <- c("HR_orig", "L_orig", "H_orig")
  comp_path <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_vs_defnA_vs_orig.csv")
  write.csv(comp, comp_path, row.names = FALSE)
  cat(sprintf("\n[30h] Comparison written to %s\n", comp_path))
  cat("\n[30h] Late, cap=2 yr (the headline window):\n")
  print(comp[comp$model == "late" & comp$cap == 2,
             c("Trial_Month", "HR_defnB", "HR_defnA_grace", "HR_orig")])
}
