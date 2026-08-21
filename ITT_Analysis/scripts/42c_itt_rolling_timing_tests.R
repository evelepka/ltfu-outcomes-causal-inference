# 42c. Rolling landmark: does the timing curve actually vary? Two tests.
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# Script 42 emits a smooth timing curve (rolling_landmark_timing.csv) but nothing
# tests whether it is FLAT. Pointwise confidence bands cannot answer that, and
# neither can asking whether two months' intervals overlap.
#
# Two tests, deliberately on the same data:
#
#   A. 3-df spline interaction  expose:ns(trial_day, df = 3)   -> 3 df
#   B. free monthly indicators  expose:factor(disengagement month) -> 5 df
#
# BOTH ARE HERE FOR ONE REASON. A spline is a more powerful test when the truth
# is smooth, so it is not interchangeable with free bins -- and the CCW's timing
# test (CCW_analysis/ccw_timing_heterogeneity.py) uses free bins. Comparing the
# spline result against the CCW would measure the SMOOTHING BASIS, not the
# design. Test B is the matched comparison; test A is the more powerful test
# under an explicit smoothness assumption. Handoff section 7 reports both, and
# they do not agree about significance -- which is the point.
#
# Report the chi-square, not the D1 F: between-imputation variance here is
# negligible, so the F denominator degrees of freedom come out absurdly large.
#
# trial_day IS the disengagement day, so the month index matches ccw_v3's
# m_dis + 1 and the two designs' months are comparable.
#
# Usage:  Rscript 42c_itt_rolling_timing_tests.R
#         N_IMP=1 Rscript 42c_itt_rolling_timing_tests.R      # fast check
# ==============================================================================
suppressPackageStartupMessages({ library(splines); library(survival); library(dplyr) })
.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R")); source(file.path(.here(), "_rolling.R"))

outcome_lk <- build_outcome_lookup()
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]
stacks <- lapply(imp_files, function(p)
  build_rolling(prepare_rolling(p, outcome_lookup = outcome_lk), comparator = "in_care"))
for (i in seq_along(stacks))
  stacks[[i]]$dmon <- factor(pmin(ceiling(stacks[[i]]$trial_day / 30.4), 6))
bs <- ns(stacks[[1]]$trial_day, df = 3)
cat(sprintf("[42c] %d imputations, in_care comparator, days %d-%d\n",
            length(stacks), DAY_MIN, DAY_MAX))
print(table(stacks[[1]]$dmon[stacks[[1]]$expose == 1]))
curves <- list(); bins <- list()

# ===================== TEST A: 3-df spline interaction =====================
# --- Rubin pooling for a vector estimand, and the D1 Wald test -------------
pool_vec <- function(bl, Vl) {
  M <- length(bl); k <- length(bl[[1]])
  qbar <- Reduce(`+`, bl) / M
  Ubar <- Reduce(`+`, Vl) / M
  B <- matrix(0, k, k)
  if (M > 1) for (b in bl) B <- B + tcrossprod(b - qbar)
  if (M > 1) B <- B / (M - 1)
  list(qbar = qbar, Ubar = Ubar, B = B, M = M, k = k,
       Tt = Ubar + (1 + 1 / M) * B)
}
d1_test <- function(pv, idx) {
  q <- pv$qbar[idx]; U <- pv$Ubar[idx, idx, drop = FALSE]
  B <- pv$B[idx, idx, drop = FALSE]; M <- pv$M; k <- length(idx)
  r1 <- (1 + 1 / M) * sum(diag(B %*% solve(U))) / k
  Tt <- (1 + r1) * U
  stat <- as.numeric(t(q) %*% solve(Tt) %*% q) / k
  dfc <- k * (M - 1)
  nu <- if (dfc > 4) 4 + (dfc - 4) * (1 + (1 - 2 / dfc) / r1)^2
        else (k + 1) * (M - 1) * (1 + 1 / r1)^2 / 2
  list(F = stat, df1 = k, df2 = nu, p = pf(stat, k, nu, lower.tail = FALSE),
       chi2 = stat * k, p_chi = pchisq(stat * k, k, lower.tail = FALSE), r1 = r1)
}

for (cap in c(2, 5)) {
  cat(sprintf("\n=================  late window, cap = %g years  =================\n", cap))
  fits <- lapply(stacks, fit_rolling, model = "late", cap = cap, timing = TRUE)
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) { cat("  not estimable\n"); next }
  bl <- lapply(fits, function(f) { b <- coef(f); b[grep("^expose", names(b))] })
  Vl <- lapply(fits, function(f) { j <- grep("^expose", names(coef(f)))
                                   vcov(f)[j, j, drop = FALSE] })
  pv <- pool_vec(bl, Vl)
  cat(sprintf("  events pooled over %d imputations; %d expose terms\n", pv$M, pv$k))

  tt <- d1_test(pv, 2:pv$k)     # drop the main effect: test only the interaction
  cat(sprintf("\n  H0: timing curve FLAT (all %d spline interactions = 0)\n", pv$k - 1))
  cat(sprintf("    D1 F(%d, %.1f) = %.2f   p = %.4f    [chi2(%d) = %.2f, p = %.4f]\n",
              tt$df1, tt$df2, tt$F, tt$p, tt$df1, tt$chi2, tt$p_chi))

  # pointwise curve and specific contrasts
  X <- function(d) { m <- cbind(1, predict(bs, d)); m }
  hr <- function(d) { x <- X(d); e <- as.vector(x %*% pv$qbar)
                      s <- sqrt(rowSums((x %*% pv$Tt) * x))
                      data.frame(day = d, HR = exp(e),
                                 lo = exp(e - 1.96 * s), hi = exp(e + 1.96 * s)) }
  g <- hr(seq(DAY_MIN, DAY_MAX, by = 5))
  pk <- g$day[which.max(g$HR)]
  cat(sprintf("\n  curve peaks at day %d (aHR %.2f)\n", pk, max(g$HR)))
  for (d in c(1, 15, 30, 60, 90, 120, 150, 176)) {
    r <- hr(d); cat(sprintf("    day %3d: aHR %.2f (%.2f-%.2f)\n", d, r$HR, r$lo, r$hi))
  }
  cat("\n  pairwise contrasts on the curve (ratio of aHRs):\n")
  for (p in list(c(pk, 1), c(pk, 176), c(1, 176), c(30, 90))) {
    cv <- as.vector(X(p[1]) - X(p[2]))
    e <- sum(cv * pv$qbar); s <- sqrt(as.numeric(t(cv) %*% pv$Tt %*% cv))
    cat(sprintf("    day %3d vs day %3d: %.2f (%.2f-%.2f)  p = %.4f\n",
                p[1], p[2], exp(e), exp(e - 1.96 * s), exp(e + 1.96 * s),
                2 * pnorm(-abs(e / s))))
  }
  g$cap <- cap; curves[[length(curves) + 1]] <- g
}

# ================== TEST B: free monthly indicators ========================
fit_bins <- function(tr, cap) {
  tr <- prep_outcome(tr, "late", cap, NULL)
  rhs <- paste(c("expose", drop_constant(tr, COVARS)), collapse = " + ")
  f <- as.formula(paste("Surv(time_out, event_out) ~", rhs,
                        "+ expose:dmon + strata(trial_day)"))
  tryCatch(coxph(f, data = tr, cluster = pid, ties = "efron"), error = function(e) NULL)
}
d1 <- function(bl, Vl) {
  M <- length(bl); k <- length(bl[[1]])
  q <- Reduce(`+`, bl) / M; U <- Reduce(`+`, Vl) / M
  B <- matrix(0, k, k); if (M > 1) { for (b in bl) B <- B + tcrossprod(b - q); B <- B/(M-1) }
  r1 <- (1 + 1/M) * sum(diag(B %*% solve(U))) / k
  st <- as.numeric(t(q) %*% solve((1 + r1) * U) %*% q) / k
  dfc <- k * (M - 1)
  nu <- if (dfc > 4) 4 + (dfc - 4) * (1 + (1 - 2/dfc)/r1)^2 else (k+1)*(M-1)*(1+1/r1)^2/2
  list(k = k, F = st, chi2 = st * k, p = pchisq(st * k, k, lower.tail = FALSE),
       nu = nu, r1 = r1, q = q, Tt = (1 + r1) * U)
}

for (cap in c(2, 5)) {
  cat(sprintf("\n=========== late window, cap = %g y, FREE MONTHLY BINS ===========\n", cap))
  fits <- Filter(Negate(is.null), lapply(stacks, fit_bins, cap = cap))
  if (!length(fits)) { cat("  not estimable\n"); next }
  nm <- names(coef(fits[[1]])); j <- grep("^expose:dmon", nm)
  bl <- lapply(fits, function(f) coef(f)[j])
  Vl <- lapply(fits, function(f) vcov(f)[j, j, drop = FALSE])
  tt <- d1(bl, Vl)
  cat(sprintf("  H0: no timing variation (%d free interaction df)\n", tt$k))
  cat(sprintf("    chi2(%d) = %.2f   p = %.4f   (r1 = %.4f)\n",
              tt$k, tt$chi2, tt$p, tt$r1))
  # per-month aHR: month 1 = main effect; others = main + interaction
  bm <- lapply(fits, function(f) { b <- coef(f)
    c(b["expose"], b["expose"] + b[j]) })
  Vm <- lapply(fits, function(f) { V <- vcov(f); k2 <- c(grep("^expose$", nm), j)
    A <- rbind(c(1, rep(0, length(j))), cbind(1, diag(length(j))))
    A %*% V[k2, k2] %*% t(A) })
  pm <- d1(bm, Vm)
  cat("  per-month aHR (disengagement month, 1-6):\n")
  for (i in seq_along(pm$q)) {
    e <- pm$q[i]; s <- sqrt(pm$Tt[i, i])
    cat(sprintf("    month %d: aHR %.2f (%.2f-%.2f)\n", i, exp(e),
                exp(e - 1.96*s), exp(e + 1.96*s)))
    bins[[length(bins) + 1]] <- data.frame(cap = cap, dmon = i, aHR = exp(e),
      CI_L = exp(e - 1.96*s), CI_H = exp(e + 1.96*s))
  }
}

tim <- bind_rows(curves); bn <- bind_rows(bins)
if (nrow(tim)) write.csv(tim, file.path(ITT_RESULTS_DIR,
  "rolling_timing_spline_curve.csv"), row.names = FALSE)
if (nrow(bn)) write.csv(bn, file.path(ITT_RESULTS_DIR,
  "rolling_timing_month_bins.csv"), row.names = FALSE)
cat(sprintf("\n[42c] wrote %d spline-curve rows and %d monthly-bin rows\n",
            nrow(tim), nrow(bn)))
