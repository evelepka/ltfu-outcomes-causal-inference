# 45b. Rolling landmark: standardized risk differences by month of disengagement
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# The 2026-08-19 decision reports TIMING from the rolling landmark, on the risk
# difference scale with a ratio alongside (handoff sections 2c and 2c-bis).
#
# Script 45 cannot supply those numbers. It fits ONE Cox model and standardizes
# through basehaz(), so its risk curves inherit proportional hazards -- and PH is
# what the Schoenfeld test in 46b/46c rejects once follow-up is extended to five
# years. Standardizing a single PH fit would mean reporting risk differences
# BECAUSE PH fails, computed under PH.
#
# So here the exposure effect gets its own coefficient in each of five disjoint
# follow-up periods and the cumulative hazard is assembled piecewise:
#
#   H(t | exposed)   = exp(gX) * sum_k exp(b_k) * [H0(min(t,u_k)) - H0(u_{k-1})]
#   H(t | unexposed) = exp(gX) * H0(t)
#
# The standardized risks therefore do not assume a constant hazard ratio. Both
# the PH-free and the PH-assuming version of the OVERALL risk difference are
# printed, so the cost of the assumption is visible rather than argued about.
#
# WHAT IT DOES NOT DO
# ---------------------------------------------------------------------------
# No confidence intervals. These are point estimates only. A cluster bootstrap on
# patient id -- resampling patients and rebuilding the stack per replicate, as
# script 45 does -- is still needed, and the by-month cells are not interpretable
# without it. Do not put these numbers in the manuscript as they stand.
#
# Within a month trial_day spans ~30 days, so it enters linearly there; across
# all months it is a 3-df spline, matching script 45.
#
# Usage:  Rscript 45b_itt_rolling_rd_by_month.R
#         N_IMP=1 Rscript 45b_itt_rolling_rd_by_month.R      # fast check
# ==============================================================================
suppressPackageStartupMessages({ library(splines); library(survival); library(dplyr) })
.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R")); source(file.path(.here(), "_rolling.R"))

HZ    <- 5                       # horizon, years
CUTS  <- c(0.5, 1, 2, 3)         # -> periods (0,.5] (.5,1] (1,2] (2,3] (3,5]
UPPER <- c(CUTS, HZ); LOWER <- c(0, CUTS); NP <- length(UPPER)
REPORT <- c(2, 5)

fit_flex <- function(tr, covars, within_month) {
  tr <- prep_outcome(tr, "overall", HZ)
  if (sum(tr$event_out) < 20) return(NULL)
  sp <- survSplit(Surv(time_out, event_out) ~ ., data = tr, cut = CUTS,
                  episode = "per")
  for (k in seq_len(NP)) sp[[paste0("ex_p", k)]] <-
    as.numeric(sp$expose == 1 & sp$per == k)
  exv <- paste0("ex_p", seq_len(NP))
  exv <- exv[vapply(exv, function(v) sum(sp[[v]]) > 0, logical(1))]
  cv  <- drop_constant(sp, covars)
  # inside one month trial_day spans ~30 d: linear. across all months: spline.
  tday <- if (within_month) "trial_day" else "ns(trial_day, df = 3)"
  f <- as.formula(paste("Surv(tstart, time_out, event_out) ~",
                        paste(c(exv, cv, tday), collapse = " + ")))
  fit <- tryCatch(coxph(f, data = sp, cluster = pid, ties = "efron",
                        x = FALSE, model = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  list(fit = fit, sp = sp, exv = exv, tr = tr)
}

std_risks <- function(o) {
  if (is.null(o)) return(NULL)
  fit <- o$fit; bh <- basehaz(fit, centered = FALSE)
  H0 <- function(t) { i <- findInterval(t, bh$time)
                      ifelse(i == 0, 0, bh$hazard[pmax(i, 1)]) }
  b <- coef(fit)
  # exposure log-HR in each period (0 if that period had no exposed person-time)
  bk <- vapply(seq_len(NP), function(k) {
    nm <- paste0("ex_p", k); if (nm %in% names(b)) unname(b[nm]) else 0 }, numeric(1))
  ex <- o$sp[o$sp$expose == 1 & o$sp$per == 1, , drop = FALSE]   # one row/patient
  if (!nrow(ex)) return(NULL)
  d0 <- ex; for (v in o$exv) d0[[v]] <- 0
  lp0 <- predict(fit, newdata = d0, type = "lp", reference = "zero")
  bind_rows(lapply(REPORT, function(t) {
    dH <- vapply(seq_len(NP), function(k)
      if (t <= LOWER[k]) 0 else H0(min(t, UPPER[k])) - H0(LOWER[k]), numeric(1))
    H1 <- sum(exp(bk) * dH)
    r1 <- mean(1 - exp(-exp(lp0) * H1))
    r0 <- mean(1 - exp(-exp(lp0) * H0(t)))
    data.frame(time_y = t, risk1 = 100 * r1, risk0 = 100 * r0,
               rd = 100 * (r1 - r0), rr = r1 / r0, n_exposed = nrow(ex))
  }))
}

# PH-assuming comparison: one exposure coefficient over the whole horizon
std_risks_ph <- function(tr, covars, within_month) {
  tr <- prep_outcome(tr, "overall", HZ)
  if (sum(tr$event_out) < 20) return(NULL)
  cv <- drop_constant(tr, covars)
  tday <- if (within_month) "trial_day" else "ns(trial_day, df = 3)"
  f <- as.formula(paste("Surv(time_out, event_out) ~",
                        paste(c("expose", cv, tday), collapse = " + ")))
  fit <- tryCatch(coxph(f, data = tr, cluster = pid, ties = "efron",
                        x = FALSE, model = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  bh <- basehaz(fit, centered = FALSE)
  H0 <- function(t) { i <- findInterval(t, bh$time)
                      ifelse(i == 0, 0, bh$hazard[pmax(i, 1)]) }
  ex <- tr[tr$expose == 1, , drop = FALSE]; if (!nrow(ex)) return(NULL)
  d1 <- ex; d1$expose <- 1; d0 <- ex; d0$expose <- 0
  lp1 <- predict(fit, newdata = d1, type = "lp", reference = "zero")
  lp0 <- predict(fit, newdata = d0, type = "lp", reference = "zero")
  bind_rows(lapply(REPORT, function(t) {
    h <- H0(t)
    data.frame(time_y = t, rd = 100 * (mean(1 - exp(-h * exp(lp1))) -
                                       mean(1 - exp(-h * exp(lp0)))))
  }))
}

pool <- function(lst, cols = NULL) {
  lst <- Filter(Negate(is.null), lst); if (!length(lst)) return(NULL)
  b <- bind_rows(lst)
  out <- b |> group_by(time_y) |>
    summarise(across(any_of(c("risk1", "risk0", "rd")), mean), .groups = "drop")
  if ("rr" %in% names(b))
    out$rr <- (b |> group_by(time_y) |>
      summarise(x = exp(mean(log(rr))), .groups = "drop"))$x
  if ("n_exposed" %in% names(b))
    out$n_exposed <- (b |> group_by(time_y) |>
      summarise(x = round(mean(n_exposed)), .groups = "drop"))$x
  out
}

outcome_lk <- build_outcome_lookup()
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]
stacks <- lapply(imp_files, function(p)
  build_rolling(prepare_rolling(p, outcome_lookup = outcome_lk), comparator = "in_care"))
for (i in seq_along(stacks))
  stacks[[i]]$dmon <- pmin(ceiling(stacks[[i]]$trial_day / 30.4), 6)
cat(sprintf("[rd] %d imputations | horizon %g y | periods %s\n", length(stacks), HZ,
            paste(sprintf("(%g,%g]", LOWER, UPPER), collapse = " ")))

cat("\n=== OVERALL, standardized to the disengaging population ===\n")
ov  <- pool(lapply(stacks, function(s) std_risks(fit_flex(s, COVARS, FALSE))))
ovp <- pool(lapply(stacks, function(s) std_risks_ph(s, COVARS, FALSE)), cols = "rd")
for (i in seq_len(nrow(ov))) {
  php <- ovp$rd[ovp$time_y == ov$time_y[i]]
  cat(sprintf("  %g y: risk %.2f%% vs %.2f%%   RD %+.2f pp   RR %.2f   "
              , ov$time_y[i], ov$risk1[i], ov$risk0[i], ov$rd[i], ov$rr[i]))
  cat(sprintf("[PH-assumed RD %+.2f, diff %+.2f]\n", php, ov$rd[i] - php))
}

cat("\n=== BY MONTH OF DISENGAGEMENT (RD at each horizon) ===\n")
rows <- list()
for (m in 1:6) {
  r <- pool(lapply(stacks, function(s)
    std_risks(fit_flex(s[s$dmon == m, , drop = FALSE], COVARS, TRUE))))
  if (is.null(r)) { cat(sprintf("  month %d: not estimable\n", m)); next }
  for (i in seq_len(nrow(r))) rows[[length(rows) + 1]] <-
    data.frame(dmon = m, time_y = r$time_y[i], risk1 = r$risk1[i],
               risk0 = r$risk0[i], rd = r$rd[i], rr = r$rr[i],
               n_exposed = r$n_exposed[i])
}
res <- bind_rows(rows)
for (t in REPORT) {
  cat(sprintf("\n  --- horizon %g y ---\n", t))
  s <- res[res$time_y == t, ]
  for (i in seq_len(nrow(s)))
    cat(sprintf("  month %d: risk %5.2f%% vs %5.2f%%   RD %+.2f pp   RR %.2f   n=%s\n",
                s$dmon[i], s$risk1[i], s$risk0[i], s$rd[i], s$rr[i],
                format(s$n_exposed[i], big.mark = ",")))
}
ov$dmon <- NA_integer_
allout <- bind_rows(ov[, c("dmon", "time_y", "risk1", "risk0", "rd", "rr",
                          "n_exposed")], res)
write.csv(allout, file.path(ITT_RESULTS_DIR, "rolling_rd_by_month.csv"),
          row.names = FALSE)
cat(sprintf("\n[45b] wrote %d rows to rolling_rd_by_month.csv (dmon=NA is overall)\n",
            nrow(allout)))
