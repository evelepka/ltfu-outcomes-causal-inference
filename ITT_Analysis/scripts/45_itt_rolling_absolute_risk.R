# 45. Rolling landmark: standardized absolute risks and risk differences
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# Scripts 42-44 give hazard ratios only. They stratify on `trial_day`, which
# means 180 separate baseline hazards and therefore NO single baseline to
# standardize against -- so absolute risks are not available from them.
#
# Reviewer 1's pivotal comment pushes toward absolute burden (which populations
# carry the preventable deaths), and Reviewer 3 notes that once the early-death
# bias is handled there is no reason to split the estimand into windows. Both
# point at the same answer: report a RISK DIFFERENCE over the full horizon, which
# needs no proportional-hazards assumption and no window selection.
#
# HOW
#   * one Cox model over the whole stacked rolling dataset, with `trial_day`
#     entering as a SPLINE rather than a stratum, so a single baseline hazard
#     exists. This is the extra assumption absolute risks require: the baseline
#     hazard is common across trials up to a smooth function of trial_day.
#     State it in the manuscript.
#   * g-formula standardization: for every EXPOSED patient, predict cumulative
#     incidence under expose = 1 and expose = 0 from the same fitted model, then
#     average. Each exposed patient appears exactly once in the stack (in their
#     own trial), so there is no re-weighting ambiguity.
#   * the estimand is therefore the effect in the DISENGAGING population (ATT):
#     "what happened to people who disengaged, versus what would have happened
#     had they stayed in care" -- which is what the paper actually claims.
#   * confidence intervals by PATIENT-level bootstrap: resample the cohort,
#     rebuild the stack, refit. Resampling stacked ROWS would be wrong, because
#     a comparator appears in many trials.
#   * multiple imputation pooled on the risk-difference scale (a difference, so
#     natural scale, not log).
#
# Usage:  Rscript 45_itt_rolling_absolute_risk.R          # point estimates
#         B=100 Rscript 45_itt_rolling_absolute_risk.R    # + bootstrap CIs
#         N_IMP=1 ... for a fast check
#
# Output: ITT_Analysis/results/rolling_landmark_absolute.csv
# ==============================================================================
suppressPackageStartupMessages({ library(splines) })
.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R")); source(file.path(.here(), "_rolling.R"))

REPORT_AT <- c(0.5, 1, 2)          # years from the trial origin
SUBGRPS   <- c("age_group", "hiv_aids", "homelessness")

# ---------------------------------------------------------------------------
# Fit one pooled model and standardize to the exposed population
# ---------------------------------------------------------------------------
standardized_risks <- function(tr, covars = COVARS, times = REPORT_AT) {
  tr <- prep_outcome(tr, "overall", max(times))
  if (sum(tr$event_out) < 20) return(NULL)
  cv  <- drop_constant(tr, covars)
  # trial_day as a SPLINE, not a stratum -> a single baseline hazard exists
  rhs <- paste(c("expose", cv, "ns(trial_day, df = 3)"), collapse = " + ")
  fit <- tryCatch(
    coxph(as.formula(paste("Surv(time_out, event_out) ~", rhs)),
          data = tr, cluster = pid, ties = "efron", x = FALSE, model = TRUE),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  bh <- basehaz(fit, centered = FALSE)            # H0(t) at covariates = 0
  H0 <- function(t) {
    i <- findInterval(t, bh$time)
    ifelse(i == 0, 0, bh$hazard[pmax(i, 1)])
  }
  # standardize to the DISENGAGING population: each exposed patient once
  ex <- tr[tr$expose == 1, , drop = FALSE]
  if (!nrow(ex)) return(NULL)
  d1 <- ex; d1$expose <- 1
  d0 <- ex; d0$expose <- 0
  lp1 <- predict(fit, newdata = d1, type = "lp", reference = "zero")
  lp0 <- predict(fit, newdata = d0, type = "lp", reference = "zero")
  out <- lapply(times, function(t) {
    h <- H0(t)
    r1 <- mean(1 - exp(-h * exp(lp1)))
    r0 <- mean(1 - exp(-h * exp(lp0)))
    data.frame(time_y = t, risk1 = 100 * r1, risk0 = 100 * r0,
               rd = 100 * (r1 - r0), rr = r1 / r0,
               hr_unstrat = unname(exp(coef(fit)["expose"])),
               n_exposed = nrow(ex))
  })
  bind_rows(out)
}

pool_rd <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (!length(lst)) return(NULL)
  b <- bind_rows(lst)
  b |> group_by(time_y) |>
    summarise(risk1 = mean(risk1), risk0 = mean(risk0),
              rd = mean(rd), rr = exp(mean(log(rr))),
              hr_unstrat = exp(mean(log(hr_unstrat))),
              n_exposed = round(mean(n_exposed)), .groups = "drop")
}

# ---------------------------------------------------------------------------
# Coherence gate: stratified fit on the SAME window is the reference
# ---------------------------------------------------------------------------
# The g-formula must drop strata(trial_day) for a single baseline hazard to exist.
# To show that substitution is harmless we refit the SAME data WITH strata, on the
# SAME (overall) window, and compare the exposure HRs. A g-formula risk difference
# whose HR disagrees with the stratified HR must not be reported.
TOL <- 0.25
strat_hr_overall <- function(tr, covars = COVARS) {
  f <- fit_rolling(tr, model = "overall", cap = max(REPORT_AT), covars = covars)
  if (is.null(f)) return(NA_real_)
  unname(exp(coef(f)["expose"]))
}
coherent <- function(hr_unstrat, hr_strat) {
  if (!is.finite(hr_unstrat) || !is.finite(hr_strat)) return(NA)
  abs(log(hr_unstrat / hr_strat)) <= log(1 + TOL)
}

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
outcome_lk <- build_outcome_lookup()
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]
B <- as.integer(Sys.getenv("B", unset = 0))
cat(sprintf("[45] standardized absolute risks | %d imputation(s) | bootstrap B=%d\n",
            length(imp_files), B))

prepped <- lapply(imp_files, function(p)
  prepare_rolling(p, outcome_lookup = outcome_lk, extra_factors = "dot_status"))
stacks  <- lapply(prepped, build_rolling, comparator = "in_care")

cat("\n--- overall, standardized to the disengaging population ---\n")
pt <- pool_rd(lapply(stacks, standardized_risks))
if (is.null(pt)) stop("[45] overall model not estimable")
for (i in seq_len(nrow(pt)))
  cat(sprintf("  %4.1f y: disengage %5.2f%%  stay in care %5.2f%%  RD %+5.2f pp  RR %.2f  (HR %.2f)\n",
              pt$time_y[i], pt$risk1[i], pt$risk0[i], pt$rd[i], pt$rr[i],
              pt$hr_unstrat[i]))

sh_all <- mean(vapply(stacks, strat_hr_overall, numeric(1)), na.rm = TRUE)
cat(sprintf("  coherence: unstratified HR %.2f vs stratified HR %.2f -> %s\n",
            pt$hr_unstrat[1], sh_all,
            if (isTRUE(coherent(pt$hr_unstrat[1], sh_all))) "OK" else "CHECK"))
rows <- pt |> mutate(stratum = "overall", level = "all",
                     reportable = isTRUE(coherent(pt$hr_unstrat[1], sh_all)))

# ---- subgroups ------------------------------------------------------------
cat("\n--- by subgroup (standardized within subgroup) ---\n")
for (sg in SUBGRPS) {
  lvls <- sort(unique(as.character(stacks[[1]][[sg]])))
  lvls <- lvls[!is.na(lvls) & nzchar(lvls) & lvls != "NA"]
  covs <- setdiff(COVARS, sg)
  for (lv in lvls) {
    r <- pool_rd(lapply(stacks, function(st) {
      sub <- st[!is.na(st[[sg]]) & as.character(st[[sg]]) == lv, , drop = FALSE]
      if (nrow(sub) < 1000 || sum(sub$expose) < 100) return(NULL)
      standardized_risks(sub, covars = covs)
    }))
    if (is.null(r)) { cat(sprintf("  %-14s %-16s not estimable\n", sg, lv)); next }
    rows <- bind_rows(rows, r |> mutate(stratum = sg, level = lv, reportable = NA))
    r2 <- r[r$time_y == 2, ]
    sh <- mean(vapply(stacks, function(st) {
      sub <- st[!is.na(st[[sg]]) & as.character(st[[sg]]) == lv, , drop = FALSE]
      if (nrow(sub) < 1000 || sum(sub$expose) < 100) return(NA_real_)
      strat_hr_overall(sub, covars = covs)
    }, numeric(1)), na.rm = TRUE)
    ok <- coherent(r2$hr_unstrat, sh)
    flag <- if (isFALSE(ok)) sprintf("  <-- NOT REPORTABLE (stratified %.2f)", sh) else ""
    rows$reportable[nrow(rows) - nrow(r) + seq_len(nrow(r))] <- isTRUE(ok) | is.na(ok)
    cat(sprintf("  %-14s %-16s n_exp=%-6s 2y: %5.2f%% vs %5.2f%%  RD %+5.2f pp  (HR %.2f)%s\n",
                sg, lv, format(r2$n_exposed, big.mark=","), r2$risk1, r2$risk0,
                r2$rd, r2$hr_unstrat, flag))
  }
}

# ---- bootstrap ------------------------------------------------------------
# Cluster bootstrap: resample PATIENTS with replacement (each draw becomes its own
# cluster), rebuild the stack, refit. Resampling stacked ROWS would be wrong,
# because a comparator appears in many trials.
#
# TWO THINGS THAT MADE THE FIRST ATTEMPT WRONG, both fixed here:
#  1. build_rolling() calls set.seed(seed) internally with a FIXED default, so
#     every replicate drew the SAME comparators. That suppressed comparator
#     sampling variability entirely and made the intervals far too narrow.
#     A distinct seed is now passed per replicate.
#  2. Imputations were cycled deterministically. They are now drawn at random per
#     replicate, so the interval carries imputation uncertainty too.
# Diagnostics below compare the bootstrap mean with the point estimate; a large
# gap means the resampling is biased and the interval must not be used.
if (B > 0) {
  cat(sprintf("\n--- cluster bootstrap (%d reps, per-rep seed + random imputation) ---\n", B))
  set.seed(4545)
  imp_pick <- sample.int(length(prepped), B, replace = TRUE)
  bs <- vector("list", B)
  t0 <- Sys.time()
  for (b in seq_len(B)) {
    d  <- prepped[[imp_pick[b]]]
    dd <- d[sample.int(nrow(d), nrow(d), replace = TRUE), , drop = FALSE]
    dd$pid <- seq_len(nrow(dd))
    bs[[b]] <- standardized_risks(
      build_rolling(dd, comparator = "in_care", seed = SEED + 1000L * b))
    if (b %% 25 == 0) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("    rep %d/%d  %.1fs/rep  ETA %.1f min\n",
                  b, B, el / b, el / b * (B - b) / 60))
    }
  }
  bs <- Filter(Negate(is.null), bs)
  cat(sprintf("  %d/%d replicates succeeded\n", length(bs), B))
  if (length(bs) >= 50) {
    bb <- bind_rows(bs)
    ci <- bb |> group_by(time_y) |>
      summarise(rd_lo = quantile(rd, .025), rd_hi = quantile(rd, .975),
                rr_lo = quantile(rr, .025), rr_hi = quantile(rr, .975),
                rd_mean = mean(rd), .groups = "drop")
    # The bootstrap was run for the OVERALL contrast only, so its interval must
    # attach ONLY to the overall rows. Joining on time_y alone broadcast the
    # overall CI onto every subgroup row, making it look as though e.g. the PLHIV
    # risk difference had been bootstrapped.
    rows <- rows |>
      left_join(ci |> select(-rd_mean) |> mutate(stratum = "overall"),
                by = c("time_y", "stratum"))
    cat("\n  overall, with cluster-bootstrap CIs:\n")
    for (i in seq_len(nrow(pt))) {
      c1 <- ci[ci$time_y == pt$time_y[i], ]
      bias <- c1$rd_mean - pt$rd[i]
      cat(sprintf("    %4.1f y: RD %+5.2f pp (%+.2f to %+.2f)   RR %.2f (%.2f-%.2f)\n",
                  pt$time_y[i], pt$rd[i], c1$rd_lo, c1$rd_hi,
                  pt$rr[i], c1$rr_lo, c1$rr_hi))
      cat(sprintf("            diagnostic: bootstrap mean %+.3f vs point %+.3f "
                  , c1$rd_mean, pt$rd[i]))
      cat(sprintf("(bias %+.3f pp, %s)\n", bias,
                  if (abs(bias) < 0.1 * abs(pt$rd[i])) "acceptable" else "CHECK"))
    }
  } else cat("  too few successful replicates for CIs\n")
}

write.csv(rows, file.path(ITT_RESULTS_DIR, "rolling_landmark_absolute.csv"),
          row.names = FALSE)
cat(sprintf("\n[45] wrote %d rows to rolling_landmark_absolute.csv\n", nrow(rows)))
