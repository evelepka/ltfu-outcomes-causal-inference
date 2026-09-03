# 49c. Death detection and cause source, by cause class -- for Reviewer 4
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# Reviewer 4 comment 1 asks how many deaths were detected through TBweb and how
# many through SIM, separately for tuberculosis and non-tuberculosis deaths, and
# how accurate each source is.
#
# 49 already answers this, but on two axes the reviewer did not ask for and one
# of which the paper no longer has:
#
#   * it restricts to the LATE window (`time_raw > 0.5`), a split removed from
#     the paper, so the table describes a population that is no longer reported;
#   * it splits by exposure arm. That split answers a different question -- is
#     ascertainment differential, Reviewer 1 comment 7 -- and is already in the
#     Results. Reviewer 4 asks for neither.
#
# It also counts STACKED TRIAL ROWS, not people. A comparator is sampled into
# many trials, so 49's table reports 14,719 comparator deaths against 12,336 that
# exist in the cohort. That is the right denominator for a differential-
# ascertainment argument, made on the analysis population, and the wrong one for
# a descriptive count of where deaths came from.
#
# This is the descriptive table: one row per person, every death in the cohort,
# no window and no arm.
#
# TWO QUESTIONS, KEPT APART
# ---------------------------------------------------------------------------
# The reviewer's sentence runs them together, but they are different:
#   detection   -- was there a death record at all, and in which system
#   cause source -- where the CAUSE came from, given a death was detected
# A death can be detected in both systems while its cause comes from only one.
#
# Usage:  Rscript 49c_death_source_by_cause.R
# Output: ITT_Analysis/results/death_source_by_cause.csv
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

cat("[49c] death detection and cause source, by cause class\n")

lookup <- build_cause_lookup(verbose = TRUE)

# Detection flags, mirroring 49's add_external_cause(). Kept as a copy rather
# than sourcing 49, which would run the whole decomposition.
raw <- read.csv(file.path(DATA_DIR, "Final_table_cleaned.csv"), stringsAsFactors = FALSE)
raw$dod <- as.Date(raw$dod, format = "%B %d, %Y")

dr <- raw[!is.na(raw$dod) & !is.na(raw$cause_of_death_code) &
            nzchar(raw$cause_of_death_code), ]
dr <- dr[order(dr$dod), ]
dr <- dr[!duplicated(dr$sinan_clean, fromLast = TRUE), c("sinan_clean")]
has_icd_ids <- unique(dr)

sim_ids <- unique(raw$sinan_clean[!is.na(raw$dod)])
oc      <- trimws(ifelse(is.na(raw$case_outcome), "", raw$case_outcome))
tbw_ids <- unique(raw$sinan_clean[oc %in% c("Obito TB", "Obito NTB")])

# Cohort deaths, one row per person
coh <- read.csv(COHORT_CSV, stringsAsFactors = FALSE)
d <- coh[coh$event_d == 1, c("sinan_clean", "itt_group"), drop = FALSE]
cat(sprintf("  deaths in the cohort: %s (%s LTFU, %s not)\n",
            format(nrow(d), big.mark = ","),
            format(sum(d$itt_group == "Loss to follow-up"), big.mark = ","),
            format(sum(d$itt_group != "Loss to follow-up"), big.mark = ",")))

d <- merge(d, lookup, by = "sinan_clean", all.x = TRUE)
for (v in c("tb_hybrid", "nontb_hybrid")) d[[v]] <- ifelse(is.na(d[[v]]), FALSE, d[[v]])
d$det_sim   <- d$sinan_clean %in% sim_ids
d$det_tbweb <- d$sinan_clean %in% tbw_ids
d$has_icd   <- d$sinan_clean %in% has_icd_ids

d$cause <- ifelse(d$tb_hybrid, "tuberculosis",
           ifelse(d$nontb_hybrid, "not tuberculosis", "unclassified"))
d$detection <- ifelse(d$det_sim & d$det_tbweb, "both",
              ifelse(d$det_sim, "SIM only",
              ifelse(d$det_tbweb, "TBweb only", "neither")))
d$cause_src <- ifelse(!d$tb_hybrid & !d$nontb_hybrid, "no cause assigned",
              ifelse(d$has_icd, "SIM ICD-10 code", "TBweb outcome"))

pct <- function(x, tot) sprintf("%.1f", 100 * x / tot)
LEV <- c("tuberculosis", "not tuberculosis", "unclassified")

t1 <- d |> count(cause, detection, name = "n") |>
  group_by(cause) |> mutate(pct_of_cause = round(100 * n / sum(n), 1)) |> ungroup()
t2 <- d |> count(cause, cause_src, name = "n") |>
  group_by(cause) |> mutate(pct_of_cause = round(100 * n / sum(n), 1)) |> ungroup()

cat("\n  HOW THE DEATH WAS DETECTED\n")
print(as.data.frame(t1 |> tidyr::pivot_wider(id_cols = cause, names_from = detection,
                                             values_from = n, values_fill = 0) |>
                      mutate(cause = factor(cause, LEV)) |> arrange(cause)),
      row.names = FALSE)
cat("\n  WHERE THE CAUSE CAME FROM\n")
print(as.data.frame(t2 |> tidyr::pivot_wider(id_cols = cause, names_from = cause_src,
                                             values_from = n, values_fill = 0) |>
                      mutate(cause = factor(cause, LEV)) |> arrange(cause)),
      row.names = FALSE)

out <- bind_rows(mutate(t1, table = "detection") |> rename(level = detection),
                 mutate(t2, table = "cause_source") |> rename(level = cause_src))
dest <- file.path(ITT_RESULTS_DIR, "death_source_by_cause.csv")
write.csv(out, dest, row.names = FALSE)

cat(sprintf("\n  totals: %s deaths | detected in SIM %s | in TBweb %s | in both %s\n",
            format(nrow(d), big.mark = ","),
            format(sum(d$det_sim), big.mark = ","),
            format(sum(d$det_tbweb), big.mark = ","),
            format(sum(d$det_sim & d$det_tbweb), big.mark = ",")))
cat(sprintf("[49c] wrote %s\n", dest))
