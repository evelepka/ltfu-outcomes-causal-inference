# 49. Reproduce the prior-literature designs in our own cohort
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# Published estimates of the LTFU-mortality effect range from 1.68 to 8.9. That
# range is not one quantity measured with varying precision: the studies differ in
# where follow-up starts and in what the comparator is, so they estimate different
# things (see docs/rationale-ccw-vs-landmark.md and the prior-literature .docx).
#
# The only interpretable way to show that design drives the number is to hold the
# cohort, the covariates and the horizon fixed and change ONLY the design. This
# script does that, so the comparison in the Discussion rests on our own data
# rather than on cross-study arithmetic.
#
# It answers Reviewer 2 comment 2 (situate our effect sizes against prior work).
#
# DESIGNS
#   A  origin = treatment initiation, comparator = CURED only
#      the completer / healthy-survivor comparator (Najera-Ortiz 2012)
#   B  origin = treatment initiation, comparator = ALL non-LTFU including deaths
#      during treatment (Cunha 2017, Garcia-Garcia 2002)
#   C  SYMMETRIC end-of-treatment landmark: both arms must be alive at day 180,
#      followed 2 y from there, comparator = CURED
#      (Pablos-Mendez 1996; Kim 2025 uses the same shape at 12 months)
#   D  ours, for reference: rolling landmark, origin = each patient's declaration
#      date, comparator = still in care at that origin. Read from
#      rolling_landmark.csv rather than refit.
#
# In A-C the exposure is "ever LTFU", classified retrospectively, exactly as in the
# studies being reproduced. That is the point: it is what those designs do.
#
# NOTE ON INDEX-EPISODE DERIVATION (invariant 9): the cured flag comes from the
# index episode's programmatic outcome. That is sound here and does NOT hit the
# invariant -- the index outcome is recorded for every patient, so there is no
# differential missingness. This file is on the reviewed allowlist in
# tools/check_conventions.py for that reason.
#
# Output: ITT_Analysis/results/literature_design_reproduction.csv
# ==============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

CAP        <- 2          # years of follow-up from each design's own origin
EOT_DAY    <- 180        # end of the intended course, for design C
CURED_CODE <- "Cura"

# ---------------------------------------------------------------------------
# index-episode programmatic outcome (see the note above)
# ---------------------------------------------------------------------------
build_index_outcome <- function() {
  raw <- read.csv(file.path(DATA_DIR, "Final_table_cleaned.csv"),
                  stringsAsFactors = FALSE)[, c("sinan_clean", "case_type",
                                                "case_outcome", "end_date")]
  raw$end_date <- as.Date(raw$end_date, format = "%B %d, %Y")
  TRANSFER <- c("Transf Outro Municipio", "Transf Outro Estado/Pais")
  novo <- raw[trimws(tools::toTitleCase(tolower(raw$case_type))) == "Novo", ]
  novo <- novo[!is.na(novo$case_outcome) & nzchar(trimws(novo$case_outcome)) &
                 novo$case_outcome != "Mud Diag" &
                 !novo$case_outcome %in% TRANSFER, ]
  novo <- novo[order(novo$end_date), ]
  first <- novo[!duplicated(novo$sinan_clean), c("sinan_clean", "case_outcome")]
  first$index_outcome <- trimws(first$case_outcome)
  first[, c("sinan_clean", "index_outcome")]
}

fit_one <- function(d, label) {
  cv <- drop_constant(d, COVARS)
  f <- tryCatch(
    coxph(as.formula(paste("Surv(t_out, e_out) ~ expose +",
                           paste(cv, collapse = " + "))),
          data = d, ties = "efron"),
    error = function(e) NULL)
  if (is.null(f)) return(NULL)
  list(est = unname(coef(f)["expose"]),
       se  = sqrt(diag(vcov(f))["expose"]),
       n   = nrow(d), ev = sum(d$e_out))
}

# ---------------------------------------------------------------------------
# the three reproduced designs, applied to one imputed dataset
# ---------------------------------------------------------------------------
build_design <- function(d, which) {
  d$cured <- !is.na(d$index_outcome) & d$index_outcome == CURED_CODE
  if (which == "A") {                      # completer comparator, origin day 0
    s <- d[d$is_ltfu | d$cured, ]
    s$expose <- as.integer(s$is_ltfu)
    s$t_out  <- pmin(s$fu_y, CAP)
    s$e_out  <- ifelse(s$fu_y > CAP, 0, s$ev)
  } else if (which == "B") {                # everyone, origin day 0
    s <- d
    s$expose <- as.integer(s$is_ltfu)
    s$t_out  <- pmin(s$fu_y, CAP)
    s$e_out  <- ifelse(s$fu_y > CAP, 0, s$ev)
  } else if (which == "C") {                # symmetric EOT landmark at day 180
    lm_y <- EOT_DAY / 365.25
    s <- d[(d$is_ltfu | d$cured) & d$fu_y > lm_y, ]   # BOTH arms alive at 180 d
    s$expose <- as.integer(s$is_ltfu)
    tt       <- s$fu_y - lm_y
    s$t_out  <- pmin(tt, CAP)
    s$e_out  <- ifelse(tt > CAP, 0, s$ev)
  } else stop("unknown design")
  s[s$t_out > 0, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
outcome_lk <- build_outcome_lookup(verbose = FALSE)
idx_lk     <- build_index_outcome()
imp_files  <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$",
                              full.names = TRUE))
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]
cat(sprintf("[49] literature-design reproduction | %d imputation(s) | cap %.0f y\n",
            length(imp_files), CAP))

preps <- lapply(imp_files, function(p) {
  d <- prepare_rolling(p, outcome_lookup = outcome_lk)
  d <- merge(d, idx_lk, by = "sinan_clean", all.x = TRUE)
  d
})
cat(sprintf("  cohort %s | LTFU %s | cured %s\n",
            format(nrow(preps[[1]]), big.mark = ","),
            format(sum(preps[[1]]$is_ltfu), big.mark = ","),
            format(sum(!is.na(preps[[1]]$index_outcome) &
                         preps[[1]]$index_outcome == CURED_CODE),
                   big.mark = ",")))

LABELS <- c(
  A = "Origin = treatment initiation; comparator = cured only (completer)",
  B = "Origin = treatment initiation; comparator = all non-LTFU incl. on-treatment deaths",
  C = "Symmetric end-of-treatment landmark (both arms alive at day 180); comparator = cured")

rows <- list()
for (w in names(LABELS)) {
  fits <- lapply(preps, function(d) fit_one(build_design(d, w), w))
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) { cat(sprintf("  %s -- not estimable\n", w)); next }
  pl <- pool_loghr(vapply(fits, `[[`, numeric(1), "est"),
                   vapply(fits, `[[`, numeric(1), "se"))
  rows[[length(rows) + 1]] <- data.frame(
    design = w, description = LABELS[[w]], cap = CAP,
    aHR = pl$hr, CI_L = pl$lo, CI_H = pl$hi, P_Value = pl$p, N_imp = pl$M,
    n_rows = round(mean(vapply(fits, `[[`, numeric(1), "n"))),
    n_events = round(mean(vapply(fits, `[[`, numeric(1), "ev"))))
  cat(sprintf("  %s  aHR %5.2f (%.2f-%.2f)  n=%s  events=%s\n", w, pl$hr, pl$lo, pl$hi,
              format(rows[[length(rows)]]$n_rows, big.mark = ","),
              format(rows[[length(rows)]]$n_events, big.mark = ",")))
}

# design D: our primary, read from the rolling output rather than refit
rl <- file.path(ITT_RESULTS_DIR, "rolling_landmark.csv")
if (file.exists(rl)) {
  r <- read.csv(rl)
  r <- r[r$comparator == "in_care" & r$model == "late", ]
  if (nrow(r) == 1) {
    rows[[length(rows) + 1]] <- data.frame(
      design = "D", cap = CAP,
      description = "OURS: origin = patient's LTFU declaration date; comparator = still in care",
      aHR = r$HR, CI_L = r$CI_L, CI_H = r$CI_H, P_Value = r$P_Value,
      N_imp = r$N_imp, n_rows = NA_real_, n_events = NA_real_)
    cat(sprintf("  D  aHR %5.2f (%.2f-%.2f)   <- ours, from rolling_landmark.csv\n",
                r$HR, r$CI_L, r$CI_H))
  }
} else cat("  D -- rolling_landmark.csv absent; run script 42 first\n")

out <- bind_rows(rows)
write.csv(out, file.path(ITT_RESULTS_DIR, "literature_design_reproduction.csv"),
          row.names = FALSE)
cat(sprintf("\n[49] wrote %d rows to literature_design_reproduction.csv\n", nrow(out)))

a <- out$aHR[out$design == "A"]; b <- out$aHR[out$design == "B"]
if (length(a) && length(b))
  cat(sprintf("\nSpread from comparator choice alone (A vs B): %.2f-fold\n", a / b))
