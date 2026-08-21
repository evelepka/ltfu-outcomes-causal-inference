# 44. Rolling landmark: (a) convention-invariance of the estimate,
#     (b) DOT-stratified timing curves.
# ==============================================================================
# (a) CONVENTION-INVARIANCE.  In the rolling design both LTFU-date conventions put
#     the time origin at the recorded CLOSURE date:
#       defnB: disengagement = closure - 30 d, grace 30 d -> origin = closure
#       defnA: disengagement = closure,        grace  0 d -> origin = closure
#     so the arms, the risk sets and therefore the EFFECT ESTIMATE are identical.
#     Only the x-axis value -- how much therapy the patient is credited with
#     receiving -- differs, by exactly 30 days.
#
#     This is a much stronger answer to Reviewer 3 comment 3.1 than a sensitivity
#     table: under the rolling design the choice of convention cannot change the
#     estimate. It changes only the horizontal placement of the timing curve. In
#     the MONTHLY design it did change the estimate (by up to 35%), because there
#     the origin was tied to the bin (m*30) rather than to the individual.
#
# (b) WHAT THE X-AXIS MEANS, AND WHY DOT MATTERS.  The back-shift recovers the last
#     day the patient was SEEN, not the last day they took drugs. Closure dates
#     cluster on a 30-day grid (42.8% within 3 d of a multiple vs 23.3% by chance),
#     and the clustering is stronger for non-DOT (48.3%) than DOT (41.1%).
#     If drugs are dispensed monthly, a self-administering patient last seen on day
#     30 holds medication to ~day 60, so their therapy received runs to the closure
#     date (defnA). Under DOT the last observed dose IS the last visit (defnB).
#     So the timing curve should be less blurred among DOT patients.
#
# Output: ITT_Analysis/results/rolling_landmark_convention.csv
#         ITT_Analysis/results/rolling_landmark_timing_by_dot.csv
# ==============================================================================
suppressPackageStartupMessages({ library(splines) })
.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  getwd()
}
source(file.path(.here(), "_paths.R")); source(file.path(.here(), "_rolling.R"))

outcome_lk <- build_outcome_lookup()
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]
cat(sprintf("[44] %d imputation(s)\n", length(imp_files)))
prepped <- lapply(imp_files, function(p)
  prepare_rolling(p, outcome_lookup = outcome_lk, extra_factors = "dot_status"))

CONFIGS <- list(list(model = "overall", cap = HORIZON_Y),
                list(model = "early",   cap = 0.5),
                list(model = "late",    cap = HORIZON_Y))

# ---- (a) convention invariance -------------------------------------------
cat("\n--- (a) is the estimate invariant to the LTFU-date convention? ---\n")
rows <- list()
for (dn in c("defnB", "defnA")) {
  st <- lapply(prepped, build_rolling, comparator = "in_care",
               carry = "dot_status", defn = dn)
  cat(sprintf("  %s: trials=%d exposed=%s\n", dn,
              length(unique(st[[1]]$trial_day)),
              format(sum(st[[1]]$expose), big.mark = ",")))
  for (cfg in CONFIGS) {
    pl <- pooled_expose(lapply(st, fit_rolling, model = cfg$model, cap = cfg$cap))
    if (is.null(pl)) next
    rows[[length(rows) + 1]] <- data.frame(defn = dn, model = cfg$model,
      HR = pl$hr, CI_L = pl$lo, CI_H = pl$hi, N_imp = pl$M)
    cat(sprintf("    %-8s aHR %.3f (%.3f-%.3f)\n", cfg$model, pl$hr, pl$lo, pl$hi))
  }
}
cv <- bind_rows(rows)
write.csv(cv, file.path(ITT_RESULTS_DIR, "rolling_landmark_convention.csv"),
          row.names = FALSE)
b <- cv$HR[cv$defn == "defnB" & cv$model == "late"]
a <- cv$HR[cv$defn == "defnA" & cv$model == "late"]
if (length(a) && length(b))
  cat(sprintf("  => late aHR defnB %.3f vs defnA %.3f  (%.2f%% apart) -> %s\n",
              b, a, 100 * abs(a - b) / b,
              if (abs(a - b) / b < 0.02) "INVARIANT" else "DIFFERS"))

# ---- (b) DOT-stratified timing -------------------------------------------
cat("\n--- (b) timing curve by DOT status (defnB labelling) ---\n")
stacks <- lapply(prepped, build_rolling, comparator = "in_care",
                 carry = "dot_status", defn = "defnB")
grid <- seq(DAY_MIN, DAY_MAX, by = 5)
trows <- list(); lrows <- list()
for (lv in c("Yes", "No")) {
  subs <- lapply(stacks, function(st)
    st[!is.na(st$dot_status) & as.character(st$dot_status) == lv, , drop = FALSE])
  covs <- setdiff(COVARS, "dot_status")
  pl <- pooled_expose(lapply(subs, fit_rolling, model = "late",
                             cap = HORIZON_Y, covars = covs))
  if (!is.null(pl)) {
    lrows[[length(lrows) + 1]] <- data.frame(dot = lv, HR = pl$hr,
      CI_L = pl$lo, CI_H = pl$hi, N_imp = pl$M)
    cat(sprintf("  DOT=%-4s exposed=%-7s late aHR %.2f (%.2f-%.2f)\n", lv,
                format(sum(subs[[1]]$expose), big.mark = ","),
                pl$hr, pl$lo, pl$hi))
  }
  tf <- Filter(Negate(is.null),
               lapply(subs, fit_rolling, model = "late", cap = HORIZON_Y,
                      covars = covs, timing = TRUE))
  if (!length(tf)) next
  bs <- ns(stacks[[1]]$trial_day, df = 3)
  per <- lapply(tf, function(f) {
    bb <- coef(f); V <- vcov(f); j <- grep("^expose", names(bb))
    X <- cbind(1, predict(bs, grid))
    list(est = as.vector(X %*% bb[j]),
         var = rowSums((X %*% V[j, j, drop = FALSE]) * X))
  })
  E <- do.call(cbind, lapply(per, `[[`, "est"))
  Vr <- do.call(cbind, lapply(per, `[[`, "var"))
  for (r in seq_along(grid)) {
    q <- pool_loghr(E[r, ], sqrt(Vr[r, ])); if (is.null(q)) next
    trows[[length(trows) + 1]] <- data.frame(dot = lv, day = grid[r],
      HR = q$hr, CI_L = q$lo, CI_H = q$hi, N_imp = q$M)
  }
}
td <- bind_rows(trows)
write.csv(td, file.path(ITT_RESULTS_DIR, "rolling_landmark_timing_by_dot.csv"),
          row.names = FALSE)
if (nrow(td)) {
  cat("\n  timing curve, DOT=Yes vs DOT=No:\n")
  for (dd in seq(DAY_MIN, DAY_MAX, by = 25)) {
    y <- td[td$dot == "Yes" & td$day == dd, ]; n <- td[td$dot == "No" & td$day == dd, ]
    if (nrow(y) && nrow(n))
      cat(sprintf("    day %3d: DOT %.2f (%.2f-%.2f)   no-DOT %.2f (%.2f-%.2f)\n",
                  dd, y$HR, y$CI_L, y$CI_H, n$HR, n$CI_L, n$CI_H))
  }
  rng <- function(v) max(v) / min(v)
  cat(sprintf("\n  peak-to-trough ratio of the curve: DOT %.2f   no-DOT %.2f\n",
              rng(td$HR[td$dot == "Yes"]), rng(td$HR[td$dot == "No"])))
  cat("  (a SHARPER curve under DOT is the predicted signature of less\n")
  cat("   measurement error in the inferred disengagement day)\n")
}
cat(sprintf("\n[44] wrote %d convention rows and %d timing rows\n", nrow(cv), nrow(td)))
