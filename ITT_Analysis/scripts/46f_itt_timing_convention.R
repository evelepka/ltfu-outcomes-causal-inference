# 46f. Time-origin convention sensitivity for the ROLLING LANDMARK design
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# Response 3.1, to Reviewer 3's pivotal comment, asserts:
#
#   "Because both conventions place the origin at case closure, the estimates
#    are unchanged to three decimal places; only the interpretation of the
#    exposure timing differs."
#
# Appendix §3.2 is supposed to back that. What it holds instead is Table S3,
# "Late-window adjusted hazard ratios under two conventions" -- monthly-landmark
# late-window values, from a design the paper no longer uses. The claim the
# reviewer will check is therefore unsupported.
#
# THE TWO CONVENTIONS
# ---------------------------------------------------------------------------
#   defnB (primary)  disengagement = recorded closure - 30 d, grace 30 d
#   defnA            disengagement = recorded closure,        grace  0 d
#
# build_rolling() implements both. They place the ORIGIN on the same calendar
# day; they differ only in the x-axis value -- how much therapy the patient is
# recorded as having received before disengaging. So the exposure coefficient
# should be invariant and the timing curve should shift by 30 days. This
# measures that rather than asserting it.
#
# Usage:  HORIZON_Y=5 Rscript 46f_itt_timing_convention.R
# Output: ITT_Analysis/results/rolling_timing_convention_<H>y.csv
# ==============================================================================
suppressPackageStartupMessages({ library(survival); library(dplyr) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

CAP    <- HORIZON_Y
CAUSES <- list(all_cause = NULL, tb_hybrid = "tb_hybrid", nontb_hybrid = "nontb_hybrid")

cat(sprintf("[46f] time-origin convention | horizon %g y\n", CAP))

lookup     <- build_cause_lookup(verbose = FALSE)
outcome_lk <- build_outcome_lookup(verbose = FALSE)

imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]
stopifnot(length(imp_files) > 0)
cat(sprintf("  %d imputation(s)\n", length(imp_files)))

rows <- list()
for (defn in c("defnB", "defnA")) {
  stacks <- lapply(imp_files, function(p)
    build_rolling(prepare_rolling(p, cause_lookup = lookup, outcome_lookup = outcome_lk),
                  comparator = "in_care", carry = c("tb_hybrid", "nontb_hybrid"),
                  defn = defn))
  s1 <- stacks[[1]]
  cat(sprintf("\n  %s: trials=%s  rows=%s  exposed=%s  median trial day=%.0f\n",
              defn, format(length(unique(s1$trial_day)), big.mark = ","),
              format(nrow(s1), big.mark = ","),
              format(sum(s1$expose), big.mark = ","),
              median(s1$trial_day[s1$expose == 1])))

  for (nm in names(CAUSES)) {
    fits <- Filter(Negate(is.null),
                   lapply(stacks, fit_rolling, model = "overall", cap = CAP,
                          cause = CAUSES[[nm]]))
    pl <- pooled_expose(fits)
    if (is.null(pl)) { cat(sprintf("    %-13s pooling failed\n", nm)); next }
    rows[[length(rows) + 1]] <- data.frame(
      convention = defn, cause = nm, model = "overall", cap = CAP,
      HR = pl$hr, CI_L = pl$lo, CI_H = pl$hi, P_Value = pl$p, N_imp = pl$M,
      median_trial_day = median(s1$trial_day[s1$expose == 1]),
      stringsAsFactors = FALSE)
    cat(sprintf("    %-13s aHR %.4f (%.4f-%.4f)\n", nm, pl$hr, pl$lo, pl$hi))
  }
  rm(stacks); invisible(gc())
}

out <- bind_rows(rows)
dest <- file.path(ITT_RESULTS_DIR, sprintf("rolling_timing_convention_%gy.csv", CAP))
write.csv(out, dest, row.names = FALSE)

cat("\n  side by side:\n")
for (nm in names(CAUSES)) {
  a <- out[out$convention == "defnA" & out$cause == nm, ]
  b <- out[out$convention == "defnB" & out$cause == nm, ]
  if (!nrow(a) || !nrow(b)) next
  cat(sprintf("    %-13s defnB %.4f | defnA %.4f | absolute difference %.5f\n",
              nm, b$HR[1], a$HR[1], abs(b$HR[1] - a$HR[1])))
}
cat(sprintf("\n[46f] wrote %s\n", dest))
cat("\n  The claim in response 3.1 holds only if these agree to three decimals.\n")
cat("  The median trial day should differ by ~30 between conventions: that is\n")
cat("  the interpretation shift the response describes, and the only one.\n")
