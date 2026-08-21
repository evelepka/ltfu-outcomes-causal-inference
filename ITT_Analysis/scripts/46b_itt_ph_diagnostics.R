# 46b. Proportional-hazards diagnostics for the rolling landmark Cox models
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# Reviewer 1 comment 17 asks us to check proportional hazards (Schoenfeld) and,
# if violated, to consider stratification, time-varying coefficients or RMST.
# Response 1.17 in the letter currently reads:
#
#   "We have tested this and report Schoenfeld residuals. [Result to be inserted.]"
#
# i.e. it asserts something that has not been done. This does it.
#
# WHICH MODELS
# ---------------------------------------------------------------------------
# The 2026-08-19 decision puts the main effect on the CCW risk-difference scale,
# which carries no PH assumption. What still depends on PH is the ROLLING
# LANDMARK, which now carries the timing analysis and the cause-specific hazard
# ratios. So those are what is tested here:
#
#   overall, whole window
#   late window            <- the one quoted in the manuscript
#   cause-specific late, TB and non-TB
#
# Values live in docs/number-registry.csv, which is untracked: this repo is
# public and the estimates are under review.
#
# Numbering: 49, 50 and 55 already collide between two people's scripts, so this
# takes a free slot adjacent to the rolling family rather than extending the mess.
#
# CAVEAT ON THE TEST
# ---------------------------------------------------------------------------
# `cox.zph` ignores the cluster-robust variance, so its p-values are based on
# the naive information matrix and are ANTI-CONSERVATIVE here: patients appear
# in many trials, so the effective sample size is smaller than the row count and
# a real p-value would be larger. Read a borderline result as "no strong
# evidence", never as "PH confirmed". With ~380k rows, a trivially small
# departure will reach significance regardless — so the SIZE of the Schoenfeld
# correlation matters more than its p-value, and both are reported.
#
# Usage:  Rscript 46b_itt_ph_diagnostics.R            # all imputations
#         N_IMP=1 Rscript 46b_itt_ph_diagnostics.R    # fast
#
# Output: ITT_Analysis/results/rolling_ph_diagnostics.csv
# ==============================================================================
suppressPackageStartupMessages({ library(survival); library(dplyr) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

# The reported timing aHRs use a FIVE-year horizon (handoff 2c), not two. Test
# both: 5 y is the one that answers R1-17 about the reported numbers, 2 y is
# kept because the cause-specific HRs are still quoted at that cap.
CAPS <- as.numeric(strsplit(Sys.getenv('CAPS', unset='2,5'), ',')[[1]])

# model label -> (window, cause)
MODELS <- list(
  list(label = "overall 0-24",      model = "overall", cause = NULL),
  list(label = "late 6-24",         model = "late",    cause = NULL),
  list(label = "late 6-24 TB",      model = "late",    cause = "tb_hybrid"),
  list(label = "late 6-24 non-TB",  model = "late",    cause = "nontb_hybrid")
)

imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]
cat(sprintf("[46b] PH diagnostics | %d imputation(s)\n", length(imp_files)))

lookup     <- build_cause_lookup()
outcome_lk <- build_outcome_lookup()

rows <- list()
for (i in seq_along(imp_files)) {
  cat(sprintf("\n  imputation %d/%d\n", i, length(imp_files)))
  tr <- build_rolling(
    prepare_rolling(imp_files[i], cause_lookup = lookup, outcome_lookup = outcome_lk),
    comparator = "in_care",
    carry = c("tb_hybrid", "nontb_hybrid"))

  for (CAP in CAPS) for (m in MODELS) {
    fit <- fit_rolling(tr, m$model, CAP, cause = m$cause)
    if (is.null(fit)) { cat(sprintf("    %-18s model failed\n", m$label)); next }

    z <- tryCatch(cox.zph(fit, transform = "km"), error = function(e) NULL)
    if (is.null(z)) { cat(sprintf("    %-18s cox.zph failed\n", m$label)); next }

    tab <- as.data.frame(z$table)
    glob <- tab["GLOBAL", ]
    exp_row <- tab[rownames(tab) == "expose", , drop = FALSE]
    if (!nrow(exp_row)) next

    # Sign of the Schoenfeld correlation for `expose`: negative means the
    # exposure effect SHRINKS over follow-up, positive means it grows.
    idx <- which(colnames(z$y) == "expose")
    rho <- if (length(idx)) suppressWarnings(
      cor(z$x, z$y[, idx], method = "spearman", use = "complete.obs")) else NA_real_

    rows[[length(rows) + 1]] <- data.frame(
      imputation = i, cap = CAP, model = m$label,
      expose_chisq = exp_row$chisq, expose_p = exp_row$p,
      global_chisq = glob$chisq, global_p = glob$p,
      expose_rho = as.numeric(rho),
      n_events = fit$nevent, n_rows = fit$n,
      stringsAsFactors = FALSE)

    cat(sprintf("    cap %gy  %-18s expose p=%-9.3g rho=%+.3f   global p=%-9.3g  (%s events)\n",
                CAP, m$label, exp_row$p, as.numeric(rho), glob$p,
                format(fit$nevent, big.mark = ",")))
  }
  rm(tr); invisible(gc())
}

out <- bind_rows(rows)
write.csv(out, file.path(ITT_RESULTS_DIR, "rolling_ph_diagnostics.csv"), row.names = FALSE)

cat("\n\n  SUMMARY across imputations (median):\n")
s <- out |> group_by(cap, model) |>
  summarise(expose_p = median(expose_p), expose_rho = median(expose_rho),
            global_p = median(global_p), events = median(n_events), .groups = "drop")
print(as.data.frame(s), row.names = FALSE, digits = 3)

cat("\n  Reading this:\n")
cat("   * expose_p is the test for the EXPOSURE effect being constant over\n")
cat("     follow-up. That is the assumption the reported hazard ratio needs.\n")
cat("   * expose_rho is the direction and size. Negative = the effect shrinks\n")
cat("     with follow-up time. |rho| below about 0.05 is a trivial departure\n")
cat("     even when p is small, because n is large.\n")
cat("   * cox.zph ignores clustering, so these p-values are anti-conservative.\n")
cat(sprintf("\n[46b] wrote %s\n",
            file.path(ITT_RESULTS_DIR, "rolling_ph_diagnostics.csv")))
