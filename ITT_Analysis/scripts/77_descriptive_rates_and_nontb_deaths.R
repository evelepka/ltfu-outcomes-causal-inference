# ==============================================================================
# 77. Two descriptive requests from review, both previously answered in the
#     response letter without the underlying analysis existing.
#
#   (a) Reviewer 3, comment 3.13: mortality RATES (deaths per 1,000 person-years)
#       alongside the cumulative incidences for the post-LTFU trajectories.
#   (b) Reviewer 4, comment 4.2: a description of deaths NOT attributed to
#       tuberculosis, by age, sex, HIV status, hospitalisation at diagnosis,
#       comorbidity, and underlying cause.
#
# Both are computed in the LTFU cohort with time measured from the LTFU
# declaration date, which is the clock the trajectory figures use.
#
# The cause split uses the paper's TWO-CLASS partition: tuberculosis, and NOT
# tuberculosis as the residual (all deaths minus tuberculosis). That is the
# owner's convention and it is why the two classes sum to the total. TB_ANY_LINE
# is refused, per _rolling.R.
#
# Usage:  Rscript ITT_Analysis/scripts/77_descriptive_rates_and_nontb_deaths.R
# Output: ITT_Analysis/results/ltfu_mortality_rates.csv
#         ITT_Analysis/results/nontb_death_description.csv
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

if (nzchar(Sys.getenv("TB_ANY_LINE")))
  stop("TB_ANY_LINE must be OFF for a manuscript number (see _rolling.R)")

coh <- read.csv(COHORT_CSV, stringsAsFactors = FALSE)
l <- coh[coh$itt_group == "Loss to follow-up", ]
l$time_d  <- suppressWarnings(as.numeric(l$time_d))
l$event_d <- suppressWarnings(as.numeric(l$event_d))
l <- l[!is.na(l$time_d) & !is.na(l$event_d) & l$time_d >= 0, ]
cat(sprintf("[77] LTFU cohort with usable follow-up: %d\n", nrow(l)))

# ---------------------------------------------------------------------------
# (a) rates per 1,000 person-years, overall and by the strata the figure uses
# ---------------------------------------------------------------------------
rate_rows <- function(df, stratum, level) {
  py <- sum(df$time_d, na.rm = TRUE)          # time_d is in YEARS from LTFU
  dd <- sum(df$event_d == 1, na.rm = TRUE)
  # exact Poisson interval on the count, scaled to the person-time
  lo <- if (dd > 0) qchisq(0.025, 2 * dd) / 2 else 0
  hi <- qchisq(0.975, 2 * dd + 2) / 2
  data.frame(stratum = stratum, level = level, n = nrow(df), deaths = dd,
             person_years = py,
             rate_per_1000_py = 1000 * dd / py,
             rate_lo = 1000 * lo / py, rate_hi = 1000 * hi / py,
             stringsAsFactors = FALSE)
}

out_a <- rate_rows(l, "overall", "all patients lost to follow-up")
STRATA <- c(hiv_aids = "HIV status", hosp_admission = "Hospitalised at diagnosis",
            homelessness = "Homelessness", age_group = "Age group",
            alcohol = "Alcohol use")
for (v in names(STRATA)) {
  if (!v %in% names(l)) next
  for (lv in sort(unique(na.omit(l[[v]])))) {
    if (!nzchar(as.character(lv))) next
    out_a <- rbind(out_a, rate_rows(l[!is.na(l[[v]]) & l[[v]] == lv, ],
                                    STRATA[[v]], as.character(lv)))
  }
}
write.csv(out_a, file.path(ITT_RESULTS_DIR, "ltfu_mortality_rates.csv"),
          row.names = FALSE)
cat("\n  (a) mortality rates per 1,000 person-years, from LTFU:\n")
print(as.data.frame(out_a |> mutate(across(where(is.numeric), ~round(., 1)))),
      row.names = FALSE)

# ---------------------------------------------------------------------------
# (b) description of deaths NOT attributed to tuberculosis
# ---------------------------------------------------------------------------
lookup <- build_cause_lookup(verbose = FALSE)
d <- l[l$event_d == 1, , drop = FALSE]
d <- merge(d, lookup[, c("sinan_clean", "tb_hybrid")], by = "sinan_clean",
           all.x = TRUE)
d$tb_hybrid <- !is.na(d$tb_hybrid) & d$tb_hybrid
d$cause_class <- ifelse(d$tb_hybrid, "tuberculosis", "not tuberculosis")
cat(sprintf("\n  (b) deaths among patients lost to follow-up: %d total = %d tuberculosis + %d not tuberculosis\n",
            nrow(d), sum(d$tb_hybrid), sum(!d$tb_hybrid)))

VARS <- c(age_group = "Age group", sex = "Sex", hiv_aids = "HIV status",
          hosp_admission = "Hospitalised at diagnosis", diabetes = "Diabetes",
          alcohol = "Alcohol use", drug_use = "Drug use",
          tobacco_use = "Tobacco use", mental_health = "Mental health condition",
          other_immuno_condition = "Other immunosuppressive condition",
          homelessness = "Homelessness", incarcerated = "Incarcerated",
          clinical_clean = "Clinical form",
          resistance_clean = "Drug-resistance status")

rows_b <- list()
for (v in names(VARS)) {
  if (!v %in% names(d)) next
  for (lv in sort(unique(na.omit(d[[v]])))) {
    if (!nzchar(as.character(lv))) next
    s <- d[!is.na(d[[v]]) & d[[v]] == lv, ]
    rows_b[[length(rows_b) + 1]] <- data.frame(
      variable = VARS[[v]], level = as.character(lv),
      n_not_tb = sum(!s$tb_hybrid),
      pct_of_not_tb = 100 * sum(!s$tb_hybrid) / sum(!d$tb_hybrid),
      n_tb = sum(s$tb_hybrid),
      pct_of_tb = 100 * sum(s$tb_hybrid) / sum(d$tb_hybrid),
      stringsAsFactors = FALSE)
  }
}
out_b <- bind_rows(rows_b)
write.csv(out_b, file.path(ITT_RESULTS_DIR, "nontb_death_description.csv"),
          row.names = FALSE)
cat("\n  columns are the % OF EACH CAUSE CLASS, so each variable sums to ~100 within a column\n")
print(as.data.frame(out_b |> mutate(across(where(is.numeric), ~round(., 1)))),
      row.names = FALSE)
cat(sprintf("\n  wrote %s\n  wrote %s\n",
            file.path(ITT_RESULTS_DIR, "ltfu_mortality_rates.csv"),
            file.path(ITT_RESULTS_DIR, "nontb_death_description.csv")))
