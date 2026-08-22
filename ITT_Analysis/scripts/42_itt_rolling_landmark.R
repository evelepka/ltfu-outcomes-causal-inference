# 42. Rolling (continuous-time) landmark analysis
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# The monthly landmark (30d/30h/30i) applies its grace period from the START of
# each disengagement window, not from each patient's own disengagement date:
#
#     start_yrs  <- (m - 1) * 30 / 365.25
#     origin_yrs <- start_yrs + GRACE_YRS      #  = m * 30 days
#
# So the time origin is the END of the disengagement window. Measured on
# imp_01, that gives a mean of only ~18 days of required post-disengagement
# survival instead of 30, and the origin falls BEFORE the LTFU declaration date
# for 90.5% of exposed patients (median 10 days early).
#
# Consequence: a patient who disengages on day 62 and dies on day 85 was never
# declared LTFU (only 23 days absent), so they are recorded as an on-treatment
# death and land in the COMPARATOR arm -- while the exposed arm contains only
# patients who survived to their own declaration. The exposed arm is selected on
# 30-day survival; the comparator only has to reach the window end. That residual
# asymmetry is a mechanism for the sub-unity early-window estimates.
#
# THIS SCRIPT removes that asymmetry by using each patient's ACTUAL
# disengagement day:
#
#   * trial t      = a disengagement day d (continuous time, no monthly bins)
#   * time origin  = d + 30 days = the patient's own declaration date, so
#                    follow-up NEVER begins before the LTFU definition is met
#   * exposed      = patients who disengaged on day d and were alive at d + 30
#   * comparator   = patients alive AND still in care at d + 30, i.e. matched on
#                    time since treatment initiation
#   * estimation   = one stratified Cox model over all stacked trials, so all
#                    ~19k disengagers contribute to a single exposure
#                    coefficient (this is why precision holds, unlike a set of
#                    separate per-month contrasts)
#   * SEs          = cluster-robust on patient id, since a comparator may be
#                    sampled into several trials
#   * timing       = expose x natural spline in d, giving a CONTINUOUS aHR
#                    curve by day of disengagement rather than six bins
#
# WHAT IT DOES NOT FIX: the exposed arm still cannot contain anyone who died
# before declaring, so death still competes with the exposure. That is what the
# clone-censor-weight analysis (CCW_analysis/ccw_v3.py) addresses. Reading the
# two together decomposes the early-window artifact into (a) misaligned origin,
# fixed here, and (b) competing exposure, fixed there.
#
# Window conventions are copied verbatim from 30h so the estimates are directly
# comparable: "late" counts events in (0.5, cap] without left-truncating, so it
# does not condition on surviving to 6 months.
#
# Output: ITT_Analysis/results/rolling_landmark.csv
#         ITT_Analysis/results/rolling_landmark_timing.csv
# ==============================================================================

suppressPackageStartupMessages({ library(splines) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
  getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))


# ---------------------------------------------------------------------------
# Load one imputed dataset and derive day-level clocks
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Build the stacked rolling-trial dataset
# ---------------------------------------------------------------------------

# window conventions copied verbatim from 30h for comparability



# Rubin pooling on the log-HR scale

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
outcome_lk <- build_outcome_lookup()
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
stopifnot(length(imp_files) > 0)
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]
cat(sprintf("[42] rolling landmark | %d imputation(s) | days %d-%d | %d comparators/exposed\n",
            length(imp_files), DAY_MIN, DAY_MAX, K_COMP))

CONFIGS <- list(list(model = "overall", cap = HORIZON_Y),
                list(model = "early",   cap = 0.5),
                list(model = "late",    cap = HORIZON_Y))
COMPARATORS <- c("in_care", "alive_only")

rows <- list(); timing_rows <- list()
for (cmp in COMPARATORS) {
  stacks <- lapply(imp_files, function(p) build_rolling(prepare_rolling(p, outcome_lookup = outcome_lk), comparator = cmp))
  s1 <- stacks[[1]]
  cat(sprintf("  comparator=%-10s trials=%d  rows=%s  exposed=%s  deaths=%s\n",
              cmp, length(unique(s1$trial_day)), format(nrow(s1), big.mark = ","),
              format(sum(s1$expose), big.mark = ","),
              format(sum(s1$event_d_num == 1), big.mark = ",")))

  for (cfg in CONFIGS) {
    fits <- lapply(stacks, fit_rolling, model = cfg$model, cap = cfg$cap)
    fits <- Filter(Negate(is.null), fits)
    if (!length(fits)) next
    est <- vapply(fits, function(f) unname(coef(f)["expose"]), numeric(1))
    se  <- vapply(fits, function(f) sqrt(diag(vcov(f))["expose"]), numeric(1))
    pl  <- pool_loghr(est, se)
    if (is.null(pl)) next
    rows[[length(rows) + 1]] <- data.frame(
      comparator = cmp, model = cfg$model, cap = cfg$cap,
      HR = pl$hr, CI_L = pl$lo, CI_H = pl$hi, P_Value = pl$p, N_imp = pl$M)
    cat(sprintf("    %-10s cap=%.1f  aHR %.2f (%.2f-%.2f)  p=%.3g\n",
                cfg$model, cfg$cap, pl$hr, pl$lo, pl$hi, pl$p))
  }

  # continuous timing curve, late window (the primary contrast)
  tf <- lapply(stacks, fit_rolling, model = "late", cap = HORIZON_Y, timing = TRUE)
  tf <- Filter(Negate(is.null), tf)
  if (length(tf)) {
    grid <- seq(DAY_MIN, DAY_MAX, by = 5)
    per <- lapply(tf, function(f) {
      b <- coef(f); V <- vcov(f)
      nm <- names(b); j <- grep("^expose", nm)
      bs <- ns(stacks[[1]]$trial_day, df = 3)
      X <- cbind(1, predict(bs, grid))
      colnames(X) <- nm[j]
      list(est = as.vector(X %*% b[j]),
           var = rowSums((X %*% V[j, j, drop = FALSE]) * X))
    })
    est <- do.call(cbind, lapply(per, `[[`, "est"))
    vr  <- do.call(cbind, lapply(per, `[[`, "var"))
    for (r in seq_along(grid)) {
      pl <- pool_loghr(est[r, ], sqrt(vr[r, ]))
      if (is.null(pl)) next
      timing_rows[[length(timing_rows) + 1]] <- data.frame(
        comparator = cmp, day = grid[r], HR = pl$hr,
        CI_L = pl$lo, CI_H = pl$hi, N_imp = pl$M)
    }
  }
}

# --- primary abandonment, reported separately (see _rolling.R) -------------
cat("\n--- primary abandonment (Abandono Primario), aligned on declaration day ---\n")
pa_stacks <- lapply(imp_files, function(p) build_primary_abandonment(prepare_rolling(p, outcome_lookup = outcome_lk)))
cat(sprintf("  trials=%d  rows=%s  exposed=%s\n",
            length(unique(pa_stacks[[1]]$trial_day)),
            format(nrow(pa_stacks[[1]]), big.mark = ","),
            format(sum(pa_stacks[[1]]$expose), big.mark = ",")))
for (cfg in CONFIGS) {
  pl <- pooled_expose(lapply(pa_stacks, fit_rolling, model = cfg$model, cap = cfg$cap))
  if (is.null(pl)) { cat(sprintf("    %-10s -- not estimable\n", cfg$model)); next }
  rows[[length(rows) + 1]] <- data.frame(
    comparator = "primary_abandonment", model = cfg$model, cap = cfg$cap,
    HR = pl$hr, CI_L = pl$lo, CI_H = pl$hi, P_Value = pl$p, N_imp = pl$M)
  cat(sprintf("    %-10s cap=%.1f  aHR %.2f (%.2f-%.2f)\n",
              cfg$model, cfg$cap, pl$hr, pl$lo, pl$hi))
}

res <- bind_rows(rows)
write.csv(res, file.path(ITT_RESULTS_DIR, "rolling_landmark.csv"), row.names = FALSE)
tim <- bind_rows(timing_rows)
if (nrow(tim)) write.csv(tim, file.path(ITT_RESULTS_DIR, "rolling_landmark_timing.csv"),
                         row.names = FALSE)
cat(sprintf("\n[42] wrote %d rows to rolling_landmark.csv and %d to rolling_landmark_timing.csv\n",
            nrow(res), nrow(tim)))

if (nrow(tim)) {
  cat("\nContinuous timing curve (late window, in_care comparator):\n")
  sub <- tim[tim$comparator == "in_care" & tim$day %in% seq(DAY_MIN, DAY_MAX, by = 25), ]
  for (i in seq_len(nrow(sub)))
    cat(sprintf("  day %3d: aHR %.2f (%.2f-%.2f)\n",
                sub$day[i], sub$HR[i], sub$CI_L[i], sub$CI_H[i]))
}
