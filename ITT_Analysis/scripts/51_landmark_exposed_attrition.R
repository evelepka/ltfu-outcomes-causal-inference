# 51. Why the rolling landmark's exposed arm is smaller than the cohort's LTFU
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# The CCW target trial keeps all 20,830 individuals lost to follow-up, because
# both clones of everyone start at treatment initiation. The rolling landmark
# indexes exposure on the DAY of disengagement, so an exposed individual needs a
# day that exists on the 1-180 axis of therapy received. Three groups have no
# such day. Until 2026-09-03 the appendix never said so, and the two designs
# reported different denominators with no explanation.
#
# Counts come from the landmark's own construction in _rolling.R, so they cannot
# drift from the design they describe.
#
# Output: ITT_Analysis/results/landmark_exposed_attrition.csv
# ==============================================================================
.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R")); source(file.path(.here(), "_rolling.R"))

d <- prepare_rolling(file.path(ITT_DATA_DIR, "itt_cohort.csv"),
                     outcome_lookup = build_outcome_lookup())
l  <- d[d$is_ltfu, ]
pa <- l$primary_aband
rest <- l[!pa, ]
lo <- !is.na(rest$dis_d) & rest$dis_d < DAY_MIN
hi <- !is.na(rest$dis_d) & rest$dis_d > DAY_MAX
keep <- !is.na(rest$dis_d) & !lo & !hi

out <- data.frame(
  step = c("ltfu_in_cohort", "primary_abandonment",
           "closure_under_31d", "disengaged_after_day180", "exposed_in_grid"),
  n = c(nrow(l), sum(pa), sum(lo), sum(hi), sum(keep)))
stopifnot(out$n[1] == sum(out$n[2:5]))          # the steps must partition

write.csv(out, file.path(ITT_RESULTS_DIR, "landmark_exposed_attrition.csv"),
          row.names = FALSE)
print(out, row.names = FALSE)
cat(sprintf("\n[51] wrote landmark_exposed_attrition.csv\n"))
