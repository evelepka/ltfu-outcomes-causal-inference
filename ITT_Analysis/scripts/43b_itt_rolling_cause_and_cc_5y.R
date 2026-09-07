# 43b. Rolling landmark at FIVE years: cause-specific, and complete-case
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# Appendix sections 3.2 (death-certificate-only cause attribution) and 4.3
# (complete-case sensitivity) report FIVE-YEAR estimates, but the two files they
# were read from --
#     ITT_Analysis/results/rolling_landmark_cause_5y.csv
#     ITT_Analysis/results/rolling_complete_case_5y.csv
# -- had no generator anywhere in the tree. They were orphans: no script wrote
# them, so nothing could refresh them, and they still carried the pre-option-C
# cohort. That was provable rather than suspected -- their all_cause early value
# was 0.889828408587957, digit for digit the value script 44 produced BEFORE the
# 2026-08-24 cohort change, against 0.916363669070408 after it. This script
# closes that hole (CLAUDE.md invariant 8: a reported number needs a source that
# can be re-run).
#
# Script 43 cannot do the job: its cap comes from HORIZON_Y, which defaults to 2,
# and it writes rolling_landmark_cause.csv. Rather than make 43 dual-horizon and
# risk its 2-year output, the five-year variant lives here and passes cap
# explicitly, so the horizon cannot be changed by an environment variable.
#
# COMPLETE-CASE DEFINITION. Stated explicitly because it had to be rediscovered:
#   * the restriction runs over COVARS EXCLUDING geo4, which is complete by
#     construction (an unmatched municipality is assigned the urban reference)
#   * empty strings count as missing, since the CSV encodes missing categoricals
#     both ways
# This reproduces the appendix's reported 102,186 of 171,048 individuals and
# 11,528 of 20,830 lost to follow-up EXACTLY, which is the evidence that it is
# the original rule and not merely a plausible one. The five-year estimates
# change because the cohort changed; the denominators do not.
#
# Output: ITT_Analysis/results/rolling_landmark_cause_5y.csv
#         ITT_Analysis/results/rolling_complete_case_5y.csv
#
# Usage:  Rscript 43b_itt_rolling_cause_and_cc_5y.R
#         N_IMP=1 Rscript 43b_itt_rolling_cause_and_cc_5y.R    # fast check
# ==============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  for (fr in rev(sys.frames())) if (!is.null(fr$ofile)) return(dirname(normalizePath(fr$ofile)))
  getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

CAP5    <- 5
# The appendix reports the complete-case comparison at BOTH horizons, and its
# two-year comparators were stale too (it quoted pooled 1.45/2.05/1.11 against
# 1.49/2.29/1.13 from the refreshed script 43). rolling_complete_case_2y.csv was
# the third orphan file, so both are written here.
CC_CAPS <- c(2, 5)
# Section 1 builds five imputed stacks and takes about 90 minutes. CC_ONLY=1
# skips it and reads the pooled comparators off disk, so the complete-case part
# (one unimputed stack) can be re-run in a couple of minutes.
CC_ONLY <- nzchar(Sys.getenv("CC_ONLY"))
CAUSES  <- c("tb_hybrid", "nontb_hybrid", "tb_broad", "tb_simonly", "nontb_simonly")
CC_CAUSES <- c("all_cause", "tb_hybrid", "nontb_hybrid")
# geo4 is complete by construction, so it cannot drive a complete-case exclusion
CC_COVARS <- setdiff(COVARS, "geo4")

outcome_lk <- build_outcome_lookup()
lookup     <- build_cause_lookup()

imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
stopifnot(length(imp_files) > 0)
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]

cat(sprintf("[43b] five-year cause-specific + complete-case | %d imputation(s)\n",
            length(imp_files)))

# ---------------------------------------------------------------------------
# 1. Cause-specific at five years, all three windows (appendix 3.1 / 3.2)
# ---------------------------------------------------------------------------
if (CC_ONLY) {
  cat("  CC_ONLY: skipping the imputed cause-specific section\n")
  cr <- read.csv(file.path(ITT_RESULTS_DIR, "rolling_landmark_cause_5y.csv"))
} else {
cat("  building stacked trials (imputed)...\n")
stacks <- lapply(imp_files, function(p) {
  d <- prepare_rolling(p, cause_lookup = lookup, outcome_lookup = outcome_lk,
                       extra_factors = c("resistance_clean"))
  build_rolling(d, comparator = "in_care",
                carry = c(CAUSES, "resistance_clean", "period"))
})
s1 <- stacks[[1]]
cat(sprintf("  trials=%d  rows=%s  exposed=%s\n",
            length(unique(s1$trial_day)), format(nrow(s1), big.mark = ","),
            format(sum(s1$expose), big.mark = ",")))

CONFIGS <- list(list(model = "overall", cap = CAP5),
                list(model = "early",   cap = 0.5),
                list(model = "late",    cap = CAP5))

cause_rows <- list()
cat("\n--- cause-specific, cap = 5 y ---\n")
for (cs in c("all_cause", CAUSES)) {
  cause_col <- if (cs == "all_cause") NULL else cs
  for (cfg in CONFIGS) {
    pl <- pooled_expose(lapply(stacks, fit_rolling, model = cfg$model,
                               cap = cfg$cap, cause = cause_col))
    if (is.null(pl)) {
      cat(sprintf("  %-14s %-8s -- not estimable\n", cs, cfg$model)); next
    }
    cause_rows[[length(cause_rows) + 1]] <- data.frame(
      cause = cs, model = cfg$model, cap = cfg$cap,
      HR = pl$hr, CI_L = pl$lo, CI_H = pl$hi, P_Value = pl$p, N_imp = pl$M)
    cat(sprintf("  %-14s %-8s aHR %5.2f (%.2f-%.2f)\n",
                cs, cfg$model, pl$hr, pl$lo, pl$hi))
  }
}
cr <- bind_rows(cause_rows)
write.csv(cr, file.path(ITT_RESULTS_DIR, "rolling_landmark_cause_5y.csv"),
          row.names = FALSE)
cat(sprintf("\n  wrote rolling_landmark_cause_5y.csv (%d rows)\n", nrow(cr)))
}

# ---------------------------------------------------------------------------
# 2. Complete-case at five years (appendix 4.3)
# ---------------------------------------------------------------------------
# Read the PRE-imputation cohort: the imp_*.csv files have no missingness left,
# so the restriction has to be applied before imputation.
cat("\n--- complete-case, cap = 5 y ---\n")
cc <- prepare_rolling(file.path(ITT_DATA_DIR, "itt_cohort.csv"),
                      cause_lookup = lookup, outcome_lookup = outcome_lk,
                      extra_factors = c("resistance_clean"))

blank_to_na <- function(x) {
  if (is.factor(x)) {
    lv <- levels(x)
    if (any(lv == "")) levels(x)[levels(x) == ""] <- NA
    return(x)
  }
  ifelse(!is.na(x) & trimws(as.character(x)) == "", NA, x)
}
for (v in CC_COVARS) cc[[v]] <- blank_to_na(cc[[v]])

complete <- Reduce(`&`, lapply(CC_COVARS, function(v) !is.na(cc[[v]])))
n_all      <- nrow(cc)
n_complete <- sum(complete)
ltfu_all      <- sum(cc$is_ltfu)
ltfu_complete <- sum(complete & cc$is_ltfu)
cat(sprintf("  complete on %d covariates (geo4 excluded, complete by construction):\n",
            length(CC_COVARS)))
cat(sprintf("    %s of %s individuals (%.1f%%)\n",
            format(n_complete, big.mark = ","), format(n_all, big.mark = ","),
            100 * n_complete / n_all))
cat(sprintf("    %s of %s lost to follow-up (%.1f%%)\n",
            format(ltfu_complete, big.mark = ","), format(ltfu_all, big.mark = ","),
            100 * ltfu_complete / ltfu_all))
# The appendix reports 102,186 / 11,528. Assert it rather than print a note: if a
# future change to the covariate set or to the missing-value encoding moves these
# denominators, that is a finding and should stop the run, not scroll past.
stopifnot(n_complete == 102186L, ltfu_complete == 11528L)
cat("  denominators match the appendix (102,186 / 11,528) -- definition confirmed\n")

cc_stack <- build_rolling(cc[complete, , drop = FALSE], comparator = "in_care",
                          carry = c(CAUSES, "resistance_clean", "period"))
cat(sprintf("  stacked complete-case rows=%s  exposed=%s\n",
            format(nrow(cc_stack), big.mark = ","),
            format(sum(cc_stack$expose), big.mark = ",")))

# Pooled comparator for a given horizon. cap 5 comes from this script's own
# output; cap 2 from script 43, which writes rolling_landmark_cause.csv. Reading
# from disk rather than from memory is what lets CC_ONLY work.
mi_ref <- function(cap) {
  f <- if (cap == CAP5) "rolling_landmark_cause_5y.csv" else "rolling_landmark_cause.csv"
  d <- read.csv(file.path(ITT_RESULTS_DIR, f))
  d[d$model == "overall" & d$cap == cap, ]
}

for (cap in CC_CAPS) {
  ref <- mi_ref(cap)
  cat(sprintf("\n  --- complete-case, cap = %g y ---\n", cap))
  cc_rows <- list()
  for (cs in CC_CAUSES) {
    cause_col <- if (cs == "all_cause") NULL else cs
    fit <- fit_rolling(cc_stack, model = "overall", cap = cap, cause = cause_col)
    pl  <- pooled_expose(list(fit))
    if (is.null(pl)) { cat(sprintf("  %-14s -- not estimable\n", cs)); next }
    tr <- prep_outcome(cc_stack, "overall", cap, cause_col)
    r  <- ref[ref$cause == cs, ]
    hm <- if (nrow(r) == 1) r$HR else NA_real_
    cc_rows[[length(cc_rows) + 1]] <- data.frame(
      cause = cs, model = "overall", cap = cap, analysis = "complete_case",
      N_individuals = n_complete, N_events = sum(tr$event_out),
      N_rows = nrow(tr),
      HR = pl$hr, CI_L = pl$lo, CI_H = pl$hi, P_Value = pl$p,
      HR_mi = hm, ratio_cc_mi = pl$hr / hm)
    cat(sprintf("  %-14s complete-case aHR %5.2f (%.2f-%.2f)   pooled %5.2f   ratio %.2f\n",
                cs, pl$hr, pl$lo, pl$hi, hm, pl$hr / hm))
  }
  ccr <- bind_rows(cc_rows)
  out <- sprintf("rolling_complete_case_%gy.csv", cap)
  write.csv(ccr, file.path(ITT_RESULTS_DIR, out), row.names = FALSE)
  cat(sprintf("  wrote %s (%d rows)   ratios %.0f%% to %.0f%% larger\n", out, nrow(ccr),
              100 * (min(ccr$ratio_cc_mi) - 1), 100 * (max(ccr$ratio_cc_mi) - 1)))
}

# Coherence: the ordering the appendix claims is unchanged must actually hold.
ov <- cr[cr$model == "overall", ]
tb_h <- ov$HR[ov$cause == "tb_hybrid"]; nt_h <- ov$HR[ov$cause == "nontb_hybrid"]
tb_s <- ov$HR[ov$cause == "tb_simonly"]; nt_s <- ov$HR[ov$cause == "nontb_simonly"]
cat(sprintf("\nordering, hybrid attribution:      TB %.2f vs non-TB %.2f  -> %s\n",
            tb_h, nt_h, if (tb_h > nt_h) "TB higher" else "CHECK"))
cat(sprintf("ordering, death-certificate only:  TB %.2f vs non-TB %.2f  -> %s\n",
            tb_s, nt_s, if (tb_s > nt_s) "TB higher" else "CHECK"))
