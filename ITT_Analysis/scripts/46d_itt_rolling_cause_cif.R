# 46d. Rolling landmark: cause-specific ABSOLUTE risk via Aalen-Johansen
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# The revised Figure 5 asks for risk differences by cause, and the handoff's
# "considered and rejected" list rules out the way they were previously
# computed:
#
#   "Cause-specific RDs, for now. Needs a cumulative incidence accounting for
#    competing deaths (Aalen-Johansen); the cause-specific hazard model in
#    script 43 is correct for HRs but does not convert into a risk difference.
#    Reporting TB/non-TB RDs computed as 1 - KM would not be coherent with the
#    all-cause RD."
#
# Owner released this on 2026-08-20. This is that estimator.
#
# WHAT WAS WRONG WITH 1 - KM
# ---------------------------------------------------------------------------
# `prep_outcome(cause = ...)` censors a competing-cause death at its own event
# time. That is the right construction for a cause-specific HAZARD ratio, but
# 1 - exp(-H_k) then answers "risk of dying of cause k in a world where no one
# dies of anything else" -- a hypothetical, and the parts do not add up to the
# all-cause risk. Script 49's "83% of the excess is TB" was computed that way
# and is superseded by whatever this produces.
#
# THE CONSTRUCTION
# ---------------------------------------------------------------------------
# Three competing causes, because TB and non-TB do NOT partition the deaths:
# roughly 12% of deaths carry no usable cause, and omitting them would leave the
# survival function too high and every CIF too large.
#
#   k in {tb, nontb, unclassified}
#   H_k(t | X)  = H0_k(t) * exp(lp_k(X))          cause-specific Cox, one per k
#   S(t | X)    = exp( - sum_k H_k(t | X) )        all three hazards
#   CIF_k(t| X) = sum over grid of S(u- | X) * dH_k(u | X)
#
# By construction sum_k CIF_k(t) = 1 - S(t) = the all-cause cumulative
# incidence, so the cause-specific risk differences ADD UP to the all-cause risk
# difference. That coherence is the whole point and is checked at the end.
#
# Standardized to the disengaging population (ATT), as scripts 45 and 45b do:
# every exposed patient contributes once, evaluated under expose = 1 and
# expose = 0 from the same fitted models.
#
# Usage:  Rscript 46d_itt_rolling_cause_cif.R              # point estimates
#         N_IMP=1 Rscript 46d_itt_rolling_cause_cif.R      # fast check
#         B=200 Rscript 46d_itt_rolling_cause_cif.R        # + cluster bootstrap
#
# DO NOT run this while another heavy job is going: a first attempt at two
# concurrent bootstraps drove the machine into swap and the rate fell fourfold.
#
# Output: ITT_Analysis/results/rolling_cause_cif.csv
# ==============================================================================
suppressPackageStartupMessages({ library(splines); library(dplyr) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

HZ      <- as.numeric(Sys.getenv("HZ", unset = "5"))    # horizon, years
# 5 years only: owner decision 2026-08-20, one horizon throughout the paper.
REPORT  <- as.numeric(strsplit(Sys.getenv("REPORT", unset = "5"), ",")[[1]])
NGRID   <- 400                                          # time grid for the integral
CAUSES  <- c("tb", "nontb", "unclass")
BYMONTH <- nzchar(Sys.getenv("BYMONTH", unset = "1"))   # also do months 1-6

# ---------------------------------------------------------------------------
# One cause-specific Cox fit. Competing-cause deaths are censored here, which is
# correct: the CIF assembly below puts them back through S(t).
# ---------------------------------------------------------------------------
fit_cause <- function(tr, which_cause, covars, within_month) {
  d <- tr
  d$time_out  <- pmin(d$time_raw, HZ)
  d$event_out <- as.integer(d$time_raw <= HZ & d[[paste0("ev_", which_cause)]])
  d <- d[d$time_out > 0, , drop = FALSE]
  if (sum(d$event_out) < 15) return(NULL)
  cv   <- drop_constant(d, covars)
  tday <- if (within_month) "trial_day" else "ns(trial_day, df = 3)"
  f <- as.formula(paste("Surv(time_out, event_out) ~",
                        paste(c("expose", cv, tday), collapse = " + ")))
  fit <- tryCatch(coxph(f, data = d, cluster = pid, ties = "efron",
                        x = FALSE, model = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  list(fit = fit, d = d, n_ev = sum(d$event_out))
}

# step function H0 evaluated on a grid
h0_on_grid <- function(fit, grid) {
  bh <- basehaz(fit, centered = FALSE)
  i  <- findInterval(grid, bh$time)
  ifelse(i == 0, 0, bh$hazard[pmax(i, 1)])
}

# ---------------------------------------------------------------------------
# Aalen-Johansen CIF for every cause, standardized to the exposed patients.
# ---------------------------------------------------------------------------
cif_standardized <- function(tr, covars, within_month) {
  fits <- lapply(CAUSES, function(k) fit_cause(tr, k, covars, within_month))
  names(fits) <- CAUSES
  if (any(vapply(fits, is.null, logical(1)))) return(NULL)

  grid <- seq(0, HZ, length.out = NGRID + 1)
  H0   <- lapply(fits, function(o) h0_on_grid(o$fit, grid))

  # one row per exposed patient, from any of the fits (same rows)
  ex <- fits[[1]]$d[fits[[1]]$d$expose == 1, , drop = FALSE]
  if (!nrow(ex)) return(NULL)

  out <- list()
  for (arm in c(1, 0)) {
    dd <- ex; dd$expose <- arm
    lp <- lapply(fits, function(o)
      predict(o$fit, newdata = dd, type = "lp", reference = "zero"))
    e_lp <- lapply(lp, exp)                       # n_exposed vector per cause

    # cumulative hazard for each cause on the grid: outer(n, grid)
    Hk <- lapply(CAUSES, function(k) outer(e_lp[[k]], H0[[k]]))
    names(Hk) <- CAUSES
    Htot <- Reduce(`+`, Hk)
    S    <- exp(-Htot)                            # survival, n x (NGRID+1)
    Slag <- cbind(1, S[, -ncol(S), drop = FALSE]) # S(u-)

    for (k in CAUSES) {
      dH  <- t(diff(t(Hk[[k]])))                  # increments, n x NGRID
      cif <- t(apply(Slag[, -1, drop = FALSE] * dH, 1, cumsum))
      for (tt in REPORT) {
        j <- max(1, which.min(abs(grid[-1] - tt)))
        out[[length(out) + 1]] <- data.frame(
          cause = k, arm = arm, time_y = tt,
          risk = 100 * mean(cif[, j]), stringsAsFactors = FALSE)
      }
    }
    # all-cause from the same S, so the parts must add to it
    for (tt in REPORT) {
      j <- which.min(abs(grid - tt))
      out[[length(out) + 1]] <- data.frame(
        cause = "all", arm = arm, time_y = tt,
        risk = 100 * mean(1 - S[, j]), stringsAsFactors = FALSE)
    }
  }
  w <- bind_rows(out) |>
    tidyr::pivot_wider(names_from = arm, values_from = risk,
                       names_prefix = "risk") |>
    mutate(rd = risk1 - risk0, rr = risk1 / risk0,
           n_exposed = nrow(ex),
           n_events = sum(vapply(fits, function(o) o$n_ev, numeric(1))))
  w
}

# ---------------------------------------------------------------------------
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]
cat(sprintf("[46d] cause-specific CIF (Aalen-Johansen) | %d imputation(s) | horizon %g y\n",
            length(imp_files), HZ))

lookup     <- build_cause_lookup()
outcome_lk <- build_outcome_lookup()
CARRY <- c("tb_hybrid", "nontb_hybrid")

prepped <- lapply(imp_files, prepare_rolling, cause_lookup = lookup,
                  outcome_lookup = outcome_lk)
stacks  <- lapply(prepped, build_rolling, comparator = "in_care", carry = CARRY)

# three mutually exclusive, exhaustive death indicators
for (i in seq_along(stacks)) {
  s <- stacks[[i]]
  died <- s$event_d_num == 1
  s$ev_tb      <- died &  s$tb_hybrid
  s$ev_nontb   <- died & !s$tb_hybrid &  s$nontb_hybrid
  s$ev_unclass <- died & !s$tb_hybrid & !s$nontb_hybrid
  stacks[[i]] <- s
}
s1 <- stacks[[1]]
cat(sprintf("  deaths in the stack: TB %s, non-TB %s, unclassified %s (%.1f%% unclassified)\n",
            format(sum(s1$ev_tb), big.mark = ","),
            format(sum(s1$ev_nontb), big.mark = ","),
            format(sum(s1$ev_unclass), big.mark = ","),
            100 * sum(s1$ev_unclass) / sum(s1$event_d_num == 1)))

pool <- function(lst) {
  lst <- Filter(Negate(is.null), lst); if (!length(lst)) return(NULL)
  bind_rows(lst) |> group_by(cause, time_y) |>
    summarise(risk1 = mean(risk1), risk0 = mean(risk0), rd = mean(rd),
              rr = exp(mean(log(rr))), n_exposed = round(mean(n_exposed)),
              .groups = "drop")
}

cat("\n=== OVERALL, standardized to the disengaging population ===\n")
ov <- pool(lapply(stacks, cif_standardized, covars = COVARS, within_month = FALSE))
ov$dmon <- NA_integer_
for (tt in REPORT) {
  cat(sprintf("\n  --- %g years ---\n", tt))
  s <- ov[ov$time_y == tt, ]
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-8s risk %6.3f%% vs %6.3f%%   RD %+7.3f pp   RR %.2f\n",
                s$cause[i], s$risk1[i], s$risk0[i], s$rd[i], s$rr[i]))
  parts <- sum(s$rd[s$cause != "all"]); whole <- s$rd[s$cause == "all"]
  cat(sprintf("  CHECK  parts %+.3f  vs  all-cause %+.3f  (gap %+.4f)\n",
              parts, whole, parts - whole))
  cat("  (the gap must be ~0: that is what Aalen-Johansen buys over 1 - KM)\n")
}

res <- ov
if (BYMONTH) {
  cat("\n=== BY MONTH OF DISENGAGEMENT ===\n")
  for (i in seq_along(stacks))
    stacks[[i]]$dmon <- pmin(ceiling(stacks[[i]]$trial_day / 30.4), 6)
  rows <- list()
  for (m in 1:6) {
    r <- pool(lapply(stacks, function(s)
      cif_standardized(s[s$dmon == m, , drop = FALSE], COVARS, TRUE)))
    if (is.null(r)) { cat(sprintf("  month %d: not estimable\n", m)); next }
    r$dmon <- m; rows[[length(rows) + 1]] <- r
  }
  bym <- bind_rows(rows)
  for (tt in REPORT) {
    cat(sprintf("\n  --- %g years ---\n", tt))
    for (m in 1:6) {
      s <- bym[bym$dmon == m & bym$time_y == tt, ]
      if (!nrow(s)) next
      g <- function(k) { v <- s$rd[s$cause == k]; if (length(v)) v else NA_real_ }
      cat(sprintf("  month %d: TB %+6.3f   non-TB %+6.3f   unclass %+6.3f   all %+6.3f\n",
                  m, g("tb"), g("nontb"), g("unclass"), g("all")))
    }
  }
  res <- bind_rows(ov, bym)
}

out <- file.path(ITT_RESULTS_DIR, "rolling_cause_cif.csv")
write.csv(res, out, row.names = FALSE)
cat(sprintf("\n[46d] wrote %s\n", out))
cat("  NOTE: point estimates only. Confidence intervals need the same cluster\n")
cat("  bootstrap as script 45b; do not put these in the manuscript without them.\n")
