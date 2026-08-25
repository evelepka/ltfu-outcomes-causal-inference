# 46g. Treatment-start exclusion sensitivity, on the ROLLING LANDMARK
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# The primary cohort requires a recorded treatment start date, which excludes
# 1,394 individuals (789 LTFU, 605 not). The main text promises the reader that
# "this exclusion is described in Appendix 1.1", and Jason's note on the appendix
# asks to keep one sentence WITH THE RESULT. The result that existed came from
# the monthly landmark and was expressed as late-window hazard ratios, so it
# could not be carried into a paper that no longer has windows.
#
# This produces the current-design equivalent: the same rolling landmark fit on
# the primary cohort and on the variant that re-includes the 1,394 with the
# start date imputed from the diagnostic date, or the notification date where
# the diagnostic date is missing (built by 01 with IMPUTE_TX_START=1).
#
# LIKE FOR LIKE
# ---------------------------------------------------------------------------
# Multiple imputation was run on the primary cohort only, so the variant has no
# imputed datasets. Both sides are therefore fitted UNIMPUTED, which makes each
# a complete-case fit. That is deliberate: the question is whether re-including
# the excluded cases moves the estimate, and the only way to isolate that is to
# hold the estimator fixed. The absolute level will match 46c, not the
# MI-pooled primary; the quantity of interest here is the difference.
#
# Usage:  HORIZON_Y=5 Rscript 46g_itt_txstart_sensitivity.R
# Output: ITT_Analysis/results/rolling_txstart_sensitivity_<H>y.csv
# ==============================================================================
suppressPackageStartupMessages({ library(survival); library(dplyr) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

CAP <- HORIZON_Y
VARIANT <- file.path(ITT_DATA_DIR, "itt_cohort_impute_start.csv")
if (!file.exists(VARIANT))
  stop("run  IMPUTE_TX_START=1 python3 01_itt_cohort_selection.py  first")

cat(sprintf("[46g] treatment-start exclusion sensitivity | horizon %g y\n", CAP))

lookup     <- build_cause_lookup(verbose = FALSE)
outcome_lk <- build_outcome_lookup(verbose = FALSE)

fit_one <- function(path, label) {
  d <- prepare_rolling(path, cause_lookup = lookup, outcome_lookup = outcome_lk)
  n_ind <- nrow(d); n_ltfu <- sum(d$is_ltfu)
  tr <- build_rolling(d, comparator = "in_care")
  f  <- fit_rolling(tr, "overall", CAP)
  if (is.null(f)) { cat(sprintf("  %-18s model failed\n", label)); return(NULL) }
  b  <- unname(coef(f)["expose"]); se <- unname(sqrt(diag(vcov(f))["expose"]))
  cat(sprintf("  %-18s N=%s  LTFU=%s  trials=%d  events=%s  aHR %.4f (%.4f-%.4f)\n",
              label, format(n_ind, big.mark = ","), format(n_ltfu, big.mark = ","),
              length(unique(tr$trial_day)), format(f$nevent, big.mark = ","),
              exp(b), exp(b - 1.96 * se), exp(b + 1.96 * se)))
  data.frame(cohort = label, N_individuals = n_ind, N_ltfu = n_ltfu,
             N_trials = length(unique(tr$trial_day)), N_events = f$nevent,
             HR = exp(b), CI_L = exp(b - 1.96 * se), CI_H = exp(b + 1.96 * se),
             stringsAsFactors = FALSE)
}

rows <- Filter(Negate(is.null), list(
  fit_one(COHORT_CSV, "primary"),
  fit_one(VARIANT,    "with imputed start")))
out <- bind_rows(rows)

dest <- file.path(ITT_RESULTS_DIR, sprintf("rolling_txstart_sensitivity_%gy.csv", CAP))
write.csv(out, dest, row.names = FALSE)

if (nrow(out) == 2) {
  a <- out$HR[out$cohort == "primary"]; b <- out$HR[out$cohort == "with imputed start"]
  cat(sprintf("\n  added individuals : %s (%s of them LTFU)\n",
              format(out$N_individuals[2] - out$N_individuals[1], big.mark = ","),
              format(out$N_ltfu[2] - out$N_ltfu[1], big.mark = ",")))
  cat(sprintf("  aHR %.4f -> %.4f   relative change %+.2f%%\n", a, b, 100 * (b - a) / a))
}
cat(sprintf("\n[46g] wrote %s\n", dest))
