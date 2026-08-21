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
# PERIOD-FLEXIBLE EXPOSURE (owner decision 2026-08-20, "option 1")
# ---------------------------------------------------------------------------
# The first version of this script gave the exposure a SINGLE coefficient per
# cause over the whole horizon, i.e. it assumed proportional hazards -- while
# 45b, which produces the all-cause number in Figure 4, does not. The two
# therefore disagreed on the same estimand: for month 1 at five years the
# all-cause risk difference from 45b was materially larger than the sum of the
# causes here. A reader adding up Figure 5 and comparing it with Figure 4 panel
# A would find one of them wrong. (Both values are in docs/number-registry.csv,
# which is untracked: this repo is public and the estimates are under review.)
#
# So the exposure now gets its own coefficient in each of five disjoint
# follow-up periods, per cause, exactly as in 45b:
#
#   H_k(t | X, exposed)   = exp(gX) * sum_j exp(b_kj) * [H0_k(min(t,U_j)) - H0_k(L_j)]
#   H_k(t | X, unexposed) = exp(gX) * H0_k(t)
#
# The cause-specific risks are then free of PH in the same way the all-cause one
# is, the two scripts estimate the same thing, and the parts still add to the
# whole because S(t) is assembled from the same three hazards.
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
suppressPackageStartupMessages({ library(splines); library(survival); library(dplyr) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

HZ      <- as.numeric(Sys.getenv("HZ", unset = "5"))    # horizon, years
# 5 years only: owner decision 2026-08-20, one horizon throughout the paper.
# Default to BOTH horizons, always. Evaluating the CIF at an extra time point is
# a vector lookup on a grid that has already been computed; the cost of a
# replicate is the Cox fits, which are identical either way. Run once with
# REPORT=5 on 2026-08-20 and the 2-year panel of Figure 5 then cost a second
# 9.5-hour bootstrap to recover numbers the first run had already computed and
# discarded. Narrow this only with a reason.
REPORT  <- as.numeric(strsplit(Sys.getenv("REPORT", unset = "2,5"), ",")[[1]])
NGRID   <- 400                                          # time grid for the integral
CAUSES  <- c("tb", "nontb", "unclass")
BYMONTH <- nzchar(Sys.getenv("BYMONTH", unset = "1"))   # also do months 1-6

# Follow-up periods for the exposure effect. Identical to 45b's CUTS -- if these
# two ever diverge the scripts stop estimating the same thing again.
CUTS  <- c(0.5, 1, 2, 3)         # -> periods (0,.5] (.5,1] (1,2] (2,3] (3,5]
UPPER <- c(CUTS, HZ); LOWER <- c(0, CUTS); NP <- length(UPPER)

# ---------------------------------------------------------------------------
# Merge any covariate level carrying fewer than MIN_EVENTS_PER_LEVEL events into
# the modal level. Same helper as 50_itt_rolling_subgroup_windows.R.
#
# Needed here and not in 45b because the causes are fitted separately: TB deaths
# are roughly a quarter of all deaths, and once those are split across five
# follow-up periods a rare level (race "Other", the sparse geo4 categories) can
# end up with zero events in a period. coxph then separates and dies with
# "exp overflow due to covariates" -- which is what happened on the first run.
# ---------------------------------------------------------------------------
MIN_EVENTS_PER_LEVEL <- 5

collapse_sparse_levels <- function(sub, vars) {
  for (v in vars) {
    if (!is.factor(sub[[v]]) && !is.character(sub[[v]])) next
    f  <- as.character(sub[[v]])
    ev <- tapply(sub$event_out, f, sum); ev <- ev[!is.na(ev)]
    if (length(ev) < 2 || all(ev >= MIN_EVENTS_PER_LEVEL)) next
    keep <- names(ev)[ev >= MIN_EVENTS_PER_LEVEL]
    ref  <- names(ev)[which.max(ev)]
    if (!length(keep)) next
    sub[[v]] <- relevel(factor(ifelse(f %in% keep, f, ref),
                               levels = union(ref, keep)), ref = ref)
  }
  sub
}

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
  d <- collapse_sparse_levels(d, covars)
  # one exposure coefficient per follow-up period, as in 45b's fit_flex()
  sp <- survSplit(Surv(time_out, event_out) ~ ., data = d, cut = CUTS,
                  episode = "per")
  for (k in seq_len(NP)) sp[[paste0("ex_p", k)]] <-
    as.numeric(sp$expose == 1 & sp$per == k)
  exv <- paste0("ex_p", seq_len(NP))
  exv <- exv[vapply(exv, function(v) sum(sp[[v]]) > 0, logical(1))]
  cv   <- drop_constant(sp, covars)
  tday <- if (within_month) "trial_day" else "ns(trial_day, df = 3)"
  # Two-stage fit. Starting from zero, the TB model diverges: ex_p3 swings from
  # +5.2 to -13.9 between Newton steps and coxph dies with "exp overflow due to
  # covariates". TB is the thinnest of the three causes, and split across five
  # periods and two arms its smallest cell holds only a few dozen events, so the
  # information matrix is poorly conditioned at the origin -- not separated, just
  # badly scaled for a Newton step. The model WITHOUT the time term converges
  # cleanly, so its coefficients are the starting point for the full one. This
  # changes where the optimiser starts, not what it maximises.
  #
  # Only TB actually needs the warm start, and the extra fit is not cheap on the
  # full stack, so try the direct fit first and fall back. Same likelihood either
  # way: a converged fit is the same MLE wherever it started from.
  #
  # A fit that merely stops iterating is NOT usable. Bootstrap replicate 1 of the
  # first trial run returned a TB risk difference far outside the plausible range while
  # printing "ran out of iterations", roughly double the point estimate; left in,
# replicates like that widen the
  # interval with numerical noise rather than sampling variability. So
  # non-convergence is treated exactly like an error: warm-start it, and if that
  # still will not converge, drop the replicate and say so in the count.
  ok <- function(o, imax) !is.null(o) && isTRUE(o$iter < imax)

  f <- as.formula(paste("Surv(tstart, time_out, event_out) ~",
                        paste(c(exv, cv, tday), collapse = " + ")))
  fit <- tryCatch(coxph(f, data = sp, cluster = pid, ties = "efron",
                        x = FALSE, model = TRUE), error = function(e) NULL)
  if (!ok(fit, 20)) {
    f0 <- as.formula(paste("Surv(tstart, time_out, event_out) ~",
                           paste(c(exv, cv), collapse = " + ")))
    fit0 <- tryCatch(coxph(f0, data = sp, ties = "efron"), error = function(e) NULL)
    if (is.null(fit0)) return(NULL)
    nm  <- colnames(model.matrix(f, data = sp))[-1]
    ini <- setNames(rep(0, length(nm)), nm)
    common <- intersect(nm, names(coef(fit0)))
    ini[common] <- coef(fit0)[common]
    fit <- tryCatch(coxph(f, data = sp, cluster = pid, ties = "efron", init = ini,
                          control = coxph.control(iter.max = 100),
                          x = FALSE, model = TRUE), error = function(e) NULL)
    if (!ok(fit, 100)) return(NULL)
  }
  if (is.null(fit)) return(NULL)
  list(fit = fit, sp = sp, exv = exv, n_ev = sum(d$event_out))
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

  # Per cause, two baseline curves on the grid: H0_k(t) for the unexposed arm,
  # and the period-weighted sum_j exp(b_kj) * [H0_k(min(t,U_j)) - H0_k(L_j)] for
  # the exposed one. Everything downstream is common to both arms.
  H0un <- list(); H0ex <- list()
  for (k in CAUSES) {
    bh <- basehaz(fits[[k]]$fit, centered = FALSE)
    H0f <- function(t) { i <- findInterval(t, bh$time)
                         ifelse(i == 0, 0, bh$hazard[pmax(i, 1)]) }
    b  <- coef(fits[[k]]$fit)
    bk <- vapply(seq_len(NP), function(j) {
      nm <- paste0("ex_p", j); if (nm %in% names(b)) unname(b[nm]) else 0 },
      numeric(1))
    H0un[[k]] <- H0f(grid)
    H0ex[[k]] <- vapply(grid, function(t)
      sum(exp(bk) * pmax(0, H0f(pmin(t, UPPER)) - H0f(LOWER))), numeric(1))
  }

  # one row per exposed patient-clone; period 1 rows only, so no one is counted
  # five times. Row sets are identical across causes: the causes differ only in
  # which deaths count as events, never in who is at risk or for how long.
  #
  # Each cause is read off ITS OWN sp, not a shared one. collapse_sparse_levels()
  # runs per cause, so a level surviving in one cause may have been merged away
  # in another; feeding TB's factor coding to the non-TB fit throws "factor geo4
  # has new levels". The row FILTER is identical across causes -- who is at risk
  # and for how long does not depend on which death counts -- so the rows still
  # line up one-to-one and S(t) is assembled over the same people.
  n_ex <- sum(fits[[1]]$sp$expose == 1 & fits[[1]]$sp$per == 1)
  if (!n_ex) return(NULL)

  # linear predictor WITHOUT the exposure terms: the exposure now lives in the
  # baseline assembly above, not in lp.
  e_lp <- list()
  for (k in CAUSES) {
    spk <- fits[[k]]$sp
    d0  <- spk[spk$expose == 1 & spk$per == 1, , drop = FALSE]
    if (nrow(d0) != n_ex) return(NULL)          # row sets must correspond
    for (v in fits[[k]]$exv) d0[[v]] <- 0
    e_lp[[k]] <- exp(predict(fits[[k]]$fit, newdata = d0, type = "lp",
                             reference = "zero"))
  }

  out <- list()
  for (arm in c(1, 0)) {
    H0a <- if (arm == 1) H0ex else H0un
    Hk  <- lapply(CAUSES, function(k) outer(e_lp[[k]], H0a[[k]]))
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
           n_exposed = n_ex,
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

# ---------------------------------------------------------------------------
# Cluster bootstrap. Same construction as 45b: patients are the sampling unit,
# the stack is rebuilt inside each replicate, the imputation is drawn at random
# and the seed differs per replicate.
#
# Draws are appended to disk after EVERY replicate, not at the end. One replicate
# here costs minutes rather than 45b's seconds, so an overnight run that is cut
# short still has to yield usable intervals -- and 45b already lost 300
# replicates once by holding everything in memory until a print statement failed.
# Re-running picks nothing up automatically, but the partial draws file can be
# read straight into the summary below.
# ---------------------------------------------------------------------------
B <- as.integer(Sys.getenv("B", unset = "0"))
if (B > 0) {
  draws_path <- file.path(ITT_RESULTS_DIR, "rolling_cause_cif_draws.csv")
  cat(sprintf("\n=== cluster bootstrap, B=%d ===\n", B))
  cat(sprintf("  draws appended to %s after every replicate\n", basename(draws_path)))
  BOOT_BYMONTH <- nzchar(Sys.getenv("BOOT_BYMONTH", unset = ""))
  cat(sprintf("  scope: overall%s\n",
              if (BOOT_BYMONTH) " + months 1-6" else " only (BOOT_BYMONTH=1 to add months)"))
  set.seed(4646)
  rm(stacks); invisible(gc())          # the stacks are rebuilt per replicate
  imp_pick <- sample.int(length(prepped), B, replace = TRUE)
  if (file.exists(draws_path)) file.remove(draws_path)
  wrote_header <- FALSE
  t0 <- Sys.time(); n_ok <- 0L

  add_ev <- function(s) {
    died <- s$event_d_num == 1
    s$ev_tb      <- died &  s$tb_hybrid
    s$ev_nontb   <- died & !s$tb_hybrid &  s$nontb_hybrid
    s$ev_unclass <- died & !s$tb_hybrid & !s$nontb_hybrid
    s
  }

  for (b in seq_len(B)) {
    d  <- prepped[[imp_pick[b]]]
    dd <- d[sample.int(nrow(d), nrow(d), replace = TRUE), , drop = FALSE]
    dd$pid <- seq_len(nrow(dd))
    sb <- add_ev(build_rolling(dd, comparator = "in_care",
                               carry = CARRY, seed = 7000L + 13L * b))
    rows <- list()
    r <- tryCatch(cif_standardized(sb, COVARS, FALSE), error = function(e) NULL)
    if (!is.null(r)) { r$dmon <- NA_integer_; rows[[length(rows) + 1]] <- r }
    if (BOOT_BYMONTH) {
      sb$dmon <- pmin(ceiling(sb$trial_day / 30.4), 6)
      for (m in 1:6) {
        rm_ <- tryCatch(cif_standardized(sb[sb$dmon == m, , drop = FALSE],
                                         COVARS, TRUE), error = function(e) NULL)
        if (!is.null(rm_)) { rm_$dmon <- m; rows[[length(rows) + 1]] <- rm_ }
      }
    }
    if (length(rows)) {
      dr <- bind_rows(rows); dr$rep <- b
      write.table(dr[, c("rep", "dmon", "cause", "time_y", "rd", "rr")],
                  draws_path, sep = ",", row.names = FALSE,
                  col.names = !wrote_header, append = wrote_header)
      wrote_header <- TRUE; n_ok <- n_ok + 1L
    }
    if (b %% 5 == 0) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("  rep %d/%d  ok=%d  %.0fs/rep  ETA %.1f h\n",
                  b, B, n_ok, el / b, el / b * (B - b) / 3600))
      invisible(gc())
    }
  }

  if (n_ok >= 20) {
    bd <- read.csv(draws_path)
    ci <- bd |> group_by(dmon, cause, time_y) |>
      summarise(rd_lo = quantile(rd, .025, na.rm = TRUE),
                rd_hi = quantile(rd, .975, na.rm = TRUE),
                rr_lo = quantile(rr, .025, na.rm = TRUE),
                rr_hi = quantile(rr, .975, na.rm = TRUE),
                n_reps = sum(!is.na(rd)), .groups = "drop")
    res2 <- left_join(res, ci, by = c("dmon", "cause", "time_y"))
    outb <- file.path(ITT_RESULTS_DIR, "rolling_cause_cif_boot.csv")
    write.csv(res2, outb, row.names = FALSE)      # write BEFORE printing
    cat(sprintf("\n[46d] wrote %s  (%d replicates completed of %d requested)\n",
                outb, n_ok, B))
    for (tt in REPORT) {
      s <- res2[is.na(res2$dmon) & res2$time_y == tt, ]
      cat(sprintf("\n  --- %g years, overall, with %d-replicate CIs ---\n", tt, n_ok))
      for (i in seq_len(nrow(s)))
        cat(sprintf("  %-8s RD %+6.3f (%+6.3f to %+6.3f)%s\n",
                    s$cause[i], s$rd[i], s$rd_lo[i], s$rd_hi[i],
                    if (isTRUE(s$rd_lo[i] < 0 && s$rd_hi[i] > 0)) "  *includes zero*" else ""))
    }
  } else {
    cat(sprintf("\n[46d] only %d replicates completed -- too few for an interval.\n", n_ok))
    cat(sprintf("  Raw draws are in %s.\n", draws_path))
  }
} else {
  cat("  NOTE: point estimates only. Confidence intervals need the cluster\n")
  cat("  bootstrap: B=150 Rscript 46d_itt_rolling_cause_cif.R\n")
}
