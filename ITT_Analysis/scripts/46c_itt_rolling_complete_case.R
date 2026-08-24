# 46c. Complete-case sensitivity for the ROLLING LANDMARK design
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# Methods says: "we repeated the primary rolling landmark cohort analysis in the
# complete-case set, restricted to individuals with non-missing data for all 13
# covariates (Appendix §6.3)". That sentence had nothing behind it.
#
# What Appendix §6.3 actually held was `30m`, run 2026-05/06: a complete-case
# version of the WITHIN-LTFU multivariable models (Cox for mortality, Fine-Gray
# for retreatment). Those models are not the primary analysis, the retreatment
# arm of them was removed from the paper entirely, and 30m predates `geo4`
# joining COVARS on 2026-08-18 -- which is where the "13" came from.
#
# This runs the sensitivity the Methods sentence describes: the same rolling
# landmark, the same covariate set, the same horizon, on individuals with
# complete covariate data and NO imputation.
#
# ON THE COVARIATE COUNT
# ---------------------------------------------------------------------------
# COVARS has 14 entries, but `geo4` cannot exclude anyone: prepare_rolling()
# folds an unmatched municipality into the GEO_REF level rather than leaving it
# NA, by design (custody settings fold in there too). So the complete-case
# restriction bites on 13 covariates whatever the model adjusts for. The script
# measures and prints the geo4 match rate so that claim is checked, not assumed,
# and the manuscript can say "13" honestly.
#
# WHAT COUNTS AS MISSING
# ---------------------------------------------------------------------------
# Mirrors 30m: empty string, "Unknown", "Ignorado", "NA", "Sem Informacao".
# A surveillance form left blank and one marked `Ignorado` are the same event.
#
# Usage:  HORIZON_Y=5 Rscript 46c_itt_rolling_complete_case.R
#         HORIZON_Y=2 Rscript 46c_itt_rolling_complete_case.R   # 2 y kept too
#
# Output: ITT_Analysis/results/rolling_complete_case_<H>y.csv
# ==============================================================================
suppressPackageStartupMessages({ library(survival); library(dplyr) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

CAP <- HORIZON_Y
MISSING_TOKENS <- c("", "Unknown", "Ignorado", "NA", "Sem Informacao")
CAUSES <- list(all_cause = NULL, tb_hybrid = "tb_hybrid", nontb_hybrid = "nontb_hybrid")

cat(sprintf("[46c] rolling landmark, COMPLETE CASE | horizon %g y\n", CAP))

lookup     <- build_cause_lookup(verbose = FALSE)
outcome_lk <- build_outcome_lookup(verbose = FALSE)

# ---------------------------------------------------------------------------
# geo4 match rate, measured BEFORE prepare_rolling() fills it
# ---------------------------------------------------------------------------
raw  <- read.csv(COHORT_CSV, stringsAsFactors = FALSE)["sinan_clean"]
geo  <- merge(raw, build_geo_lookup(verbose = FALSE), by = "sinan_clean", all.x = TRUE)
n_geo_na <- sum(is.na(geo$geo4))
cat(sprintf("  geo4 unmatched municipalities: %s of %s (%.2f%%), folded into '%s'\n",
            format(n_geo_na, big.mark = ","), format(nrow(geo), big.mark = ","),
            100 * n_geo_na / nrow(geo), GEO_REF))

# ---------------------------------------------------------------------------
# Cohort, un-imputed, restricted to complete covariate data
# ---------------------------------------------------------------------------
d <- prepare_rolling(COHORT_CSV, cause_lookup = lookup, outcome_lookup = outcome_lk)
N_total <- nrow(d)

for (v in COVARS) {
  x <- as.character(d[[v]]); x[x %in% MISSING_TOKENS] <- NA
  d[[v]] <- if (v == "age_group")
    factor(x, levels = c("15-24", "25-44", "45-64", "≥65")) else factor(x)
}

miss <- vapply(COVARS, function(v) sum(is.na(d[[v]])), integer(1))
cat("\n  missing per covariate:\n")
for (v in COVARS)
  cat(sprintf("    %-16s %7s (%5.2f%%)%s\n", v, format(miss[[v]], big.mark = ","),
              100 * miss[[v]] / N_total,
              if (v == "geo4") "   <- filled by construction, cannot exclude" else ""))

N_ltfu_total <- sum(d$is_ltfu)
cc <- complete.cases(d[, COVARS])
cat(sprintf("\n  cohort            %s\n  complete case     %s (%.1f%%)\n  dropped           %s (%.1f%%)\n",
            format(N_total, big.mark = ","), format(sum(cc), big.mark = ","),
            100 * sum(cc) / N_total, format(sum(!cc), big.mark = ","),
            100 * sum(!cc) / N_total))
cat(sprintf("  covariates that can exclude: %d of %d\n",
            sum(miss > 0), length(COVARS)))

d <- droplevels(d[cc, , drop = FALSE])
d$geo4 <- relevel(d$geo4, ref = GEO_REF)
d$pid  <- seq_len(nrow(d))          # rebuild: build_rolling indexes on it
cat(sprintf("  LTFU retained     %s of %s (%.1f%%)\n",
            format(sum(d$is_ltfu), big.mark = ","),
            format(N_ltfu_total, big.mark = ","),
            100 * sum(d$is_ltfu) / N_ltfu_total))

# ---------------------------------------------------------------------------
# Same design, same comparator, no imputation
# ---------------------------------------------------------------------------
tr <- build_rolling(d, comparator = "in_care",
                    carry = c("tb_hybrid", "nontb_hybrid"))
cat(sprintf("\n  stacked trials    %s\n  rows              %s\n  exposed rows      %s\n",
            format(length(unique(tr$trial_day)), big.mark = ","),
            format(nrow(tr), big.mark = ","),
            format(sum(tr$expose), big.mark = ",")))

# MI-pooled counterparts, for the side-by-side the appendix needs
mi_path <- file.path(ITT_RESULTS_DIR, sprintf("rolling_landmark_cause_%gy.csv", CAP))
mi <- if (file.exists(mi_path)) read.csv(mi_path) else NULL

rows <- list()
cat("\n  cause         n_events        CC aHR (95% CI)          MI aHR (95% CI)      ratio\n")
for (nm in names(CAUSES)) {
  fit <- fit_rolling(tr, "overall", CAP, cause = CAUSES[[nm]])
  if (is.null(fit)) { cat(sprintf("  %-13s model failed\n", nm)); next }
  b  <- unname(coef(fit)["expose"])
  se <- unname(sqrt(diag(vcov(fit))["expose"]))
  hr <- exp(b); lo <- exp(b - 1.96 * se); hi <- exp(b + 1.96 * se)

  m <- if (!is.null(mi))
    mi[mi$cause == nm & mi$model == "overall" & mi$cap == CAP, ] else NULL
  mi_hr <- if (!is.null(m) && nrow(m)) m$HR[1] else NA_real_

  rows[[length(rows) + 1]] <- data.frame(
    cause = nm, model = "overall", cap = CAP, analysis = "complete_case",
    N_individuals = nrow(d), N_events = fit$nevent, N_rows = fit$n,
    HR = hr, CI_L = lo, CI_H = hi,
    P_Value = 2 * pnorm(-abs(b / se)),
    HR_mi = mi_hr, ratio_cc_mi = hr / mi_hr,
    stringsAsFactors = FALSE)

  cat(sprintf("  %-13s %8s   %5.2f (%.2f-%.2f)   %s   %s\n",
              nm, format(fit$nevent, big.mark = ","), hr, lo, hi,
              if (is.na(mi_hr)) "        --        "
              else sprintf("%5.2f (%.2f-%.2f)", mi_hr, m$CI_L[1], m$CI_H[1]),
              if (is.na(mi_hr)) "" else sprintf("%.3f", hr / mi_hr)))
}

out <- bind_rows(rows)
dest <- file.path(ITT_RESULTS_DIR, sprintf("rolling_complete_case_%gy.csv", CAP))
write.csv(out, dest, row.names = FALSE)
cat(sprintf("\n[46c] wrote %s\n", dest))
cat("\n  Reading this: ratio_cc_mi near 1 means dropping incomplete records does\n")
cat("  not move the estimate, which is what the MAR assumption predicts. A ratio\n")
cat("  far from 1 is evidence against it, not a reason to prefer either fit --\n")
cat("  the complete-case set is smaller AND selected on recording practice.\n")
