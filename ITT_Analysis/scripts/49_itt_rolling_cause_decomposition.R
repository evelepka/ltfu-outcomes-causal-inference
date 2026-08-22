# 49. Rolling landmark: ABSOLUTE decomposition of the excess by cause,
#     plus an external-cause outcome
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# The revision replaces "non-TB mortality is a negative control" with "the
# deaths caused by LTFU concentrate in TB deaths". That is a claim about
# COMPOSITION, and it does not follow from comparing two hazard ratios.
#
# A hazard ratio of 3.87 for TB death and 1.99 for non-TB death is consistent
# with EITHER more excess TB deaths OR more excess non-TB deaths, because non-TB
# deaths are far more numerous at baseline. This is the same relative-versus-
# absolute inversion the paper already reports for subgroups, where PLHIV carry
# the largest absolute risk difference on a mid-range hazard ratio.
#
# So this script answers the question the sentence actually makes: of the excess
# deaths attributable to LTFU, how many are TB and how many are non-TB?
#
#   * standardized 2-year risk under expose=1 and expose=0, per cause, in the
#     disengaging population (ATT) -- same g-formula machinery as script 45
#   * the risk differences are then directly comparable ACROSS causes, because
#     they share a denominator: the same exposed patients
#   * run for BOTH cause classifications side by side (owner's decision,
#     2026-08-19):
#       hybrid  = SIM ICD code, falling back to the TBweb programmatic outcome
#       simonly = SIM ICD code only
#     They disagree about which deaths have a cause at all, and the coverage is
#     differential by arm (see the ascertainment table this script writes), so
#     reporting one without the other would hide that.
#
# EXTERNAL-CAUSE OUTCOME (ICD-10 V01-Y98)
# ---------------------------------------------------------------------------
# Under the "concentrates in TB deaths" framing, external-cause mortality
# (accident, violence) is the one outcome with no plausible causal path from
# treatment interruption -- post-TB lung disease and cardiovascular sequelae can
# produce a "non-TB" death caused by TB, but they cannot produce a road traffic
# death. So an elevated external-cause estimate is evidence of residual
# confounding by social vulnerability rather than of effect, and a null one
# supports the concentration claim.
#
# It is necessarily SIM-only: TBweb records a programmatic outcome, not a
# mechanism of injury, so there is no hybrid equivalent. Interpret against the
# SIM coverage reported below, not against the hybrid column.
#
# The external-cause flag is derived HERE rather than in build_cause_lookup(),
# deliberately: editing `_rolling.R` obligates re-running 42, 43, 45 and 46
# (CLAUDE.md invariant 3), and this is an additive derived column, not a change
# to shared logic. If external cause becomes a reported analysis, move it into
# build_cause_lookup() and re-run the family.
#
# Usage:  Rscript 49_itt_rolling_cause_decomposition.R          # point estimates
#         B=200 Rscript 49_itt_rolling_cause_decomposition.R    # + bootstrap CIs
#         N_IMP=1 Rscript 49_...                                # fast check
#
# Output: ITT_Analysis/results/rolling_cause_decomposition.csv
#         ITT_Analysis/results/rolling_cause_ascertainment.csv
# ==============================================================================
suppressPackageStartupMessages({ library(splines); library(dplyr) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

HORIZON <- 2                       # years from the trial origin
B       <- as.integer(Sys.getenv("B", unset = "0"))

# Causes to decompose. all_cause first: the per-cause risk differences should
# roughly add up to it, and printing them together makes a failure obvious.
CAUSES <- c("all_cause",
            "tb_hybrid", "nontb_hybrid",
            "tb_simonly", "nontb_simonly",
            "ext_simonly")

# ---------------------------------------------------------------------------
# External-cause flag + the ascertainment table, from the same raw extract and
# the same per-patient rule build_cause_lookup() uses (last coded death record
# per patient, NOT the index episode -- ADR-0003).
# ---------------------------------------------------------------------------
add_external_cause <- function(lookup) {
  raw <- read.csv(file.path(DATA_DIR, "Final_table_cleaned.csv"),
                  stringsAsFactors = FALSE)
  raw$dod <- as.Date(raw$dod, format = "%B %d, %Y")
  dr <- raw[!is.na(raw$dod) & !is.na(raw$cause_of_death_code) &
              nzchar(raw$cause_of_death_code), ]
  dr <- dr[order(dr$dod), ]
  dr <- dr[!duplicated(dr$sinan_clean, fromLast = TRUE),
           c("sinan_clean", "cause_of_death_code")]
  dr$cause_of_death_code <- toupper(trimws(dr$cause_of_death_code))
  dr$ext_simonly <- grepl("^[VWXY]", dr$cause_of_death_code)
  dr$has_icd     <- TRUE

  out <- merge(lookup, dr[, c("sinan_clean", "ext_simonly", "has_icd")],
               by = "sinan_clean", all.x = TRUE)
  out$ext_simonly[is.na(out$ext_simonly)] <- FALSE
  out$has_icd[is.na(out$has_icd)]         <- FALSE

  # --- how was the DEATH detected, as opposed to how its cause was coded? ---
  # R4 comment 1 asks for the TBweb-versus-SIM split separately for tuberculosis
  # and non-tuberculosis deaths. Detection and cause-coding are different things
  # and the paper conflates them, so both are tabulated below.
  #   det_sim   : a date of death exists in the mortality registry linkage
  #   det_tbweb : an Obito outcome is recorded in TBweb, in ANY episode (ADR-0003)
  sim_ids <- unique(raw$sinan_clean[!is.na(raw$dod)])
  oc      <- trimws(ifelse(is.na(raw$case_outcome), "", raw$case_outcome))
  tbw_ids <- unique(raw$sinan_clean[oc %in% c("Obito TB", "Obito NTB")])
  out$det_sim   <- out$sinan_clean %in% sim_ids
  out$det_tbweb <- out$sinan_clean %in% tbw_ids

  cat(sprintf("  [cause] external-cause (V01-Y98) deaths with an ICD code: %d\n",
              sum(out$ext_simonly)))
  cat(sprintf("  [cause] death detection: SIM %d, TBweb %d, both %d\n",
              sum(out$det_sim), sum(out$det_tbweb), sum(out$det_sim & out$det_tbweb)))
  out
}

# ---------------------------------------------------------------------------
# Standardized per-cause risk in the disengaging population.
# Identical machinery to script 45, with prep_outcome()'s `cause` argument:
# competing-cause deaths are censored at their own event time.
# ---------------------------------------------------------------------------
cause_risks <- function(tr, cause, covars = COVARS, horizon = HORIZON) {
  cs <- if (identical(cause, "all_cause")) NULL else cause
  tr <- prep_outcome(tr, "overall", horizon, cause = cs)
  if (sum(tr$event_out) < 20) return(NULL)
  cv  <- drop_constant(tr, covars)
  rhs <- paste(c("expose", cv, "ns(trial_day, df = 3)"), collapse = " + ")
  # Capture, rather than swallow, the "coefficient may be infinite" warning:
  # in a sparse cause it means a covariate level has no events, and the affected
  # estimate must not be reported without knowing which level.
  warned <- character(0)
  fit <- withCallingHandlers(
    tryCatch(
      coxph(as.formula(paste("Surv(time_out, event_out) ~", rhs)),
            data = tr, cluster = pid, ties = "efron", x = FALSE, model = TRUE),
      error = function(e) NULL),
    warning = function(w) { warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (is.null(fit)) return(NULL)

  # Name the offending term(s) instead of leaving "variable 20, 21".
  infinite_terms <- character(0)
  if (any(grepl("may be infinite", warned))) {
    nm <- names(coef(fit))
    idx <- unique(unlist(regmatches(warned, gregexpr("[0-9]+", warned))))
    idx <- suppressWarnings(as.integer(idx))
    idx <- idx[!is.na(idx) & idx >= 1 & idx <= length(nm)]
    infinite_terms <- nm[idx]
  }

  bh <- basehaz(fit, centered = FALSE)
  i  <- findInterval(horizon, bh$time)
  H0 <- if (i == 0) 0 else bh$hazard[max(i, 1)]

  ex <- tr[tr$expose == 1, , drop = FALSE]
  if (!nrow(ex)) return(NULL)
  d1 <- ex; d1$expose <- 1
  d0 <- ex; d0$expose <- 0
  r1 <- mean(1 - exp(-H0 * exp(predict(fit, newdata = d1, type = "lp", reference = "zero"))))
  r0 <- mean(1 - exp(-H0 * exp(predict(fit, newdata = d0, type = "lp", reference = "zero"))))

  data.frame(cause = cause, risk1 = 100 * r1, risk0 = 100 * r0,
             rd = 100 * (r1 - r0),
             hr_unstrat = unname(exp(coef(fit)["expose"])),
             n_events = sum(tr$event_out), n_exposed = nrow(ex),
             nonconverged = paste(infinite_terms, collapse = ";"),
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
stopifnot(length(imp_files) > 0)
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]

cat(sprintf("[49] rolling cause decomposition | %d imputation(s) | horizon %g y | B=%d\n",
            length(imp_files), HORIZON, B))

lookup     <- add_external_cause(build_cause_lookup())
outcome_lk <- build_outcome_lookup()

CARRY <- c("tb_hybrid", "nontb_hybrid", "tb_simonly", "nontb_simonly",
           "ext_simonly", "has_icd", "det_sim", "det_tbweb")

cat("  building stacked trials...\n")
prepped <- lapply(imp_files, function(p)
  prepare_rolling(p, cause_lookup = lookup, outcome_lookup = outcome_lk))
stacks <- lapply(prepped, build_rolling, comparator = "in_care", carry = CARRY)

# --- ascertainment table: where does the cause come from, by arm? ------------
s1 <- stacks[[1]]
deaths <- s1[s1$event_d_num == 1, , drop = FALSE]
asc <- deaths |>
  group_by(expose) |>
  summarise(n_deaths   = n(),
            icd_coded  = sum(has_icd),
            pct_icd    = 100 * mean(has_icd),
            hybrid_cls = sum(tb_hybrid | nontb_hybrid),
            pct_hybrid = 100 * mean(tb_hybrid | nontb_hybrid),
            simonly_cls = sum(tb_simonly | nontb_simonly),
            pct_simonly = 100 * mean(tb_simonly | nontb_simonly),
            .groups = "drop")
cat("\n  cause ascertainment among deaths in the stacked data:\n")
print(as.data.frame(asc), row.names = FALSE)
write.csv(asc, file.path(ITT_RESULTS_DIR, "rolling_cause_ascertainment.csv"),
          row.names = FALSE)

# --- R4 comment 1: detection source x cause category, LATE window ------------
# The reviewer asks how many late deaths were detected through TBweb and how
# many through SIM, separately for tuberculosis and non-tuberculosis deaths.
# Detection (was there a death record at all) and cause coding (what it was
# coded as) are different questions; both are reported, because the paper
# currently answers neither.
late <- s1[s1$event_d_num == 1 & s1$time_raw > 0.5 & s1$time_raw <= HORIZON, , drop = FALSE]
late$arm   <- ifelse(late$expose == 1, "LTFU", "in care")
late$cause <- ifelse(late$tb_hybrid, "TB",
              ifelse(late$nontb_hybrid, "non-TB", "unclassified"))
late$detection <- ifelse(late$det_sim & late$det_tbweb, "both",
                  ifelse(late$det_sim, "SIM only",
                  ifelse(late$det_tbweb, "TBweb only", "neither")))
# where the CAUSE came from, which is the narrower question behind R2 comment 3
late$cause_src <- ifelse(!late$tb_hybrid & !late$nontb_hybrid, "none",
                  ifelse(late$has_icd, "SIM ICD code", "TBweb outcome"))

r4 <- late |>
  count(arm, cause, detection, cause_src, name = "n") |>
  group_by(arm) |>
  mutate(pct_of_arm = round(100 * n / sum(n), 1)) |>
  ungroup() |>
  arrange(arm, cause, detection)

cat(sprintf("\n  R4.1 -- late-window deaths (%.1f-%g y from origin), n=%d:\n",
            0.5, HORIZON, nrow(late)))
cat("\n  detection source by cause category:\n")
print(as.data.frame(late |> count(arm, cause, detection, name = "n") |>
                      tidyr::pivot_wider(names_from = detection, values_from = n,
                                         values_fill = 0)), row.names = FALSE)
cat("\n  where the cause code came from:\n")
print(as.data.frame(late |> count(arm, cause, cause_src, name = "n") |>
                      tidyr::pivot_wider(names_from = cause_src, values_from = n,
                                         values_fill = 0)), row.names = FALSE)

write.csv(r4, file.path(ITT_RESULTS_DIR, "rolling_late_death_sources.csv"),
          row.names = FALSE)
cat(sprintf("\n  wrote %s\n",
            file.path(ITT_RESULTS_DIR, "rolling_late_death_sources.csv")))

# Stop here when only the ascertainment tables are wanted, so this can be
# refreshed without re-running the models or the bootstrap.
if (nzchar(Sys.getenv("ASCERTAINMENT_ONLY"))) {
  cat("\n[49] ASCERTAINMENT_ONLY set -- stopping before the models.\n")
  quit(save = "no", status = 0)
}

# --- point estimates ---------------------------------------------------------
cat("\n  fitting per-cause models...\n")
point <- bind_rows(lapply(CAUSES, function(cs) {
  per_imp <- lapply(stacks, function(tr) cause_risks(tr, cs))
  per_imp <- Filter(Negate(is.null), per_imp)
  if (!length(per_imp)) {
    cat(sprintf("    %-14s SKIPPED (too few events)\n", cs)); return(NULL)
  }
  b <- bind_rows(per_imp)
  # Carry the non-convergence flag through the imputation pooling: an estimate
  # whose model had an infinite coefficient must not be reported silently.
  nc <- unique(unlist(strsplit(b$nonconverged[nzchar(b$nonconverged)], ";")))
  r <- data.frame(cause = cs,
                  risk1 = mean(b$risk1), risk0 = mean(b$risk0),
                  rd = mean(b$rd), hr_unstrat = exp(mean(log(b$hr_unstrat))),
                  n_events = round(mean(b$n_events)),
                  n_exposed = round(mean(b$n_exposed)),
                  M = nrow(b),
                  nonconverged = paste(nc, collapse = ";"),
                  stringsAsFactors = FALSE)
  cat(sprintf("    %-14s RD %+6.3f pp   risk1 %5.2f%%  risk0 %5.2f%%  HR %4.2f  (%d events)%s\n",
              cs, r$rd, r$risk1, r$risk0, r$hr_unstrat, r$n_events,
              if (nzchar(r$nonconverged)) paste0("  <== NON-CONVERGED: ", r$nonconverged) else ""))
  r
}))

# --- bootstrap CIs -----------------------------------------------------------
# Cluster bootstrap, copied from script 45 so the intervals in this table are
# built the same way as the ones already in the manuscript. Two mistakes are
# avoided here because 45 already hit them (see its comments):
#   1. resample PATIENTS in the prepped data and REBUILD the stack. Resampling
#      rows of the built stack suppresses comparator sampling variability and
#      makes the interval far too narrow.
#   2. pass a distinct seed per replicate — build_rolling() seeds internally with
#      a fixed default, so every replicate would otherwise draw the same
#      comparators; and draw the imputation at random so the interval carries
#      imputation uncertainty.
if (B > 0) {
  cat(sprintf("\n  cluster bootstrap (%d reps, per-rep seed + random imputation)\n", B))
  set.seed(4949)
  imp_pick <- sample.int(length(prepped), B, replace = TRUE)
  boots <- vector("list", B)
  t0 <- Sys.time()
  for (b in seq_len(B)) {
    d  <- prepped[[imp_pick[b]]]
    dd <- d[sample.int(nrow(d), nrow(d), replace = TRUE), , drop = FALSE]
    dd$pid <- seq_len(nrow(dd))
    tr <- build_rolling(dd, comparator = "in_care", carry = CARRY,
                        seed = SEED + 1000L * b)
    boots[[b]] <- bind_rows(lapply(CAUSES, function(cs) {
      r <- cause_risks(tr, cs)
      if (is.null(r)) return(NULL)
      data.frame(cause = cs, rd = r$rd, stringsAsFactors = FALSE)
    }))
    if (b %% 25 == 0) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("    rep %d/%d  %.1fs/rep  ETA %.1f min\n",
                  b, B, el / b, el / b * (B - b) / 60))
    }
  }
  bd <- bind_rows(Filter(Negate(is.null), boots))
  ci <- bd |> group_by(cause) |>
    summarise(rd_lo = quantile(rd, 0.025, na.rm = TRUE),
              rd_hi = quantile(rd, 0.975, na.rm = TRUE),
              rd_mean = mean(rd, na.rm = TRUE), B_ok = n(), .groups = "drop")
  point <- left_join(point, ci, by = "cause")
  # Diagnostic: a bootstrap mean far from the point estimate means the
  # resampling is biased and the interval must NOT be used (45's convention).
  cat("\n  bootstrap diagnostic (bias = boot mean - point):\n")
  for (i in seq_len(nrow(point))) {
    bias <- point$rd_mean[i] - point$rd[i]
    cat(sprintf("    %-14s RD %+6.3f (%+.3f to %+.3f)  bias %+.3f  %s  [%d reps]\n",
                point$cause[i], point$rd[i], point$rd_lo[i], point$rd_hi[i], bias,
                if (isTRUE(abs(bias) < 0.1 * abs(point$rd[i]))) "acceptable" else "CHECK",
                point$B_ok[i]))
  }
}

out <- file.path(ITT_RESULTS_DIR, "rolling_cause_decomposition.csv")
write.csv(point, out, row.names = FALSE)
cat(sprintf("\n[49] wrote %s\n", out))

# --- does the decomposition add up? -----------------------------------------
if (all(c("all_cause", "tb_hybrid", "nontb_hybrid") %in% point$cause)) {
  g <- function(x) point$rd[point$cause == x]
  cat(sprintf("\n  CHECK  all-cause RD %+.3f pp   vs   TB %+.3f + non-TB %+.3f = %+.3f (hybrid)\n",
              g("all_cause"), g("tb_hybrid"), g("nontb_hybrid"),
              g("tb_hybrid") + g("nontb_hybrid")))
  cat("  (hybrid does not classify every death, so the parts sum to less than the whole;\n")
  cat("   a LARGE gap means unclassified deaths dominate and neither split is reportable.)\n")
}
