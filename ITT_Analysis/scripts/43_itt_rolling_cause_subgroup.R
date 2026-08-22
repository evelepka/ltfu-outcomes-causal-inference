# 43. Rolling landmark: cause-specific outcomes and subgroups
# ==============================================================================
# Same design as 42 (see _rolling.R): one trial per disengagement DAY, origin at
# the patient's own LTFU declaration date (d + 30), comparator alive and still in
# care at that origin, all trials stacked into one stratified Cox.
#
# Cause attribution uses the 2026-08-16 fix (Obito outcome from ANY episode, not
# just the index one). Cause-specific analyses censor competing-cause deaths at
# their own event time, matching 30e/30i.
#
# Subgroup contrasts restrict the stacked data to a level of the subgroup
# variable, so both arms come from within that subgroup, and drop that variable
# from the covariate set.
#
# Output: ITT_Analysis/results/rolling_landmark_cause.csv
#         ITT_Analysis/results/rolling_landmark_subgroups.csv
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

CAUSES  <- c("tb_hybrid", "nontb_hybrid", "tb_broad", "tb_simonly", "nontb_simonly")
SUBGRPS <- c("age_group", "sex", "hiv_aids", "homelessness",
             "resistance_clean", "period")

outcome_lk <- build_outcome_lookup()
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
stopifnot(length(imp_files) > 0)
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]

cat(sprintf("[43] rolling cause-specific + subgroups | %d imputation(s)\n",
            length(imp_files)))
lookup <- build_cause_lookup()

cat("  building stacked trials...\n")
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
cat("  events available in the stacked data (late window contributions):\n")
for (cs in CAUSES)
  cat(sprintf("    %-14s %s\n", cs, format(sum(s1[[cs]]), big.mark = ",")))

# ---------------------------------------------------------------------------
# 1. Cause-specific, all three windows
# ---------------------------------------------------------------------------
cause_rows <- list()
CONFIGS <- list(list(model = "overall", cap = HORIZON_Y),
                list(model = "early",   cap = 0.5),
                list(model = "late",    cap = HORIZON_Y))
cat("\n--- cause-specific ---\n")
for (cs in c("all_cause", CAUSES)) {
  cause_col <- if (cs == "all_cause") NULL else cs
  for (cfg in CONFIGS) {
    pl <- pooled_expose(lapply(stacks, fit_rolling, model = cfg$model,
                               cap = cfg$cap, cause = cause_col))
    if (is.null(pl)) { cat(sprintf("  %-14s %-8s -- not estimable\n", cs, cfg$model)); next }
    cause_rows[[length(cause_rows) + 1]] <- data.frame(
      cause = cs, model = cfg$model, cap = cfg$cap,
      HR = pl$hr, CI_L = pl$lo, CI_H = pl$hi, P_Value = pl$p, N_imp = pl$M)
    if (cfg$model == "late")
      cat(sprintf("  %-14s late   aHR %5.2f (%.2f-%.2f)  p=%.3g\n",
                  cs, pl$hr, pl$lo, pl$hi, pl$p))
  }
}
cr <- bind_rows(cause_rows)
write.csv(cr, file.path(ITT_RESULTS_DIR, "rolling_landmark_cause.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Subgroups, late window (the primary contrast)
# ---------------------------------------------------------------------------
sub_rows <- list()
cat("\n--- subgroups (late window, all-cause) ---\n")
for (sg in SUBGRPS) {
  if (is.null(s1[[sg]])) { cat(sprintf("  %s: not present, skipped\n", sg)); next }
  lvls <- sort(unique(as.character(s1[[sg]])))
  lvls <- lvls[!is.na(lvls) & nzchar(lvls) & lvls != "NA"]
  covs_sg <- setdiff(COVARS, sg)
  for (lv in lvls) {
    fits <- lapply(stacks, function(st) {
      sub <- st[!is.na(st[[sg]]) & as.character(st[[sg]]) == lv, , drop = FALSE]
      # a subgroup must retain both arms inside at least a few trials
      if (nrow(sub) < 500 || sum(sub$expose) < 50) return(NULL)
      fit_rolling(sub, model = "late", cap = HORIZON_Y, covars = covs_sg)
    })
    pl <- pooled_expose(fits)
    if (is.null(pl)) { cat(sprintf("  %-18s %-16s -- not estimable\n", sg, lv)); next }
    sub_rows[[length(sub_rows) + 1]] <- data.frame(
      Subgroup = sg, Level = lv, model = "late", cap = HORIZON_Y,
      HR = pl$hr, CI_L = pl$lo, CI_H = pl$hi, P_Value = pl$p, N_imp = pl$M)
    cat(sprintf("  %-18s %-16s aHR %5.2f (%.2f-%.2f)\n", sg, lv, pl$hr, pl$lo, pl$hi))
  }
}
sr <- bind_rows(sub_rows)
if (nrow(sr)) write.csv(sr, file.path(ITT_RESULTS_DIR, "rolling_landmark_subgroups.csv"),
                        row.names = FALSE)

cat(sprintf("\n[43] wrote %d cause rows and %d subgroup rows\n", nrow(cr), nrow(sr)))

# coherence check: cause-specific estimates should straddle all-cause
lt <- cr[cr$model == "late", ]
ac <- lt$HR[lt$cause == "all_cause"]
if (length(ac) == 1) {
  tb <- lt$HR[lt$cause == "tb_hybrid"]; nt <- lt$HR[lt$cause == "nontb_hybrid"]
  cat(sprintf("\ncoherence (late): TB %.2f  >  all-cause %.2f  >  non-TB %.2f  -> %s\n",
              tb, ac, nt,
              if (length(tb) && length(nt) && tb > ac && nt < ac) "PASS" else "CHECK"))
}
