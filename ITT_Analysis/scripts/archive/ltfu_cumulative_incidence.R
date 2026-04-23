# Cumulative Incidence Analysis — LTFU Cohort
# ============================================================
# Retreatment: Competing Risks (Fine-Gray) — death is competing event
# Mortality:   Direct Cumulative Incidence (1 - Kaplan-Meier)
# Milestones:  1 year, 5 years, 10 years, end of study (max follow-up)
# ============================================================

library(survival)
library(cmprsk)
library(dplyr)

df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)
ltfu <- df[df$itt_group == "Loss to follow-up", ]
cat("N LTFU:", nrow(ltfu), "\n")

milestones_yr <- c(1, 5, 10)

# ─────────────────────────────────────────────────────────────
# Helper: extract CIF at milestone with CI
extract_cif <- function(cif_obj, t, cause = 1) {
  idx <- max(which(cif_obj$time <= t), na.rm = TRUE)
  if (length(idx) == 0 || is.na(idx)) return(c(NA, NA, NA))
  est <- cif_obj$est[idx, cause]
  lo  <- cif_obj$lower[idx, cause]
  hi  <- cif_obj$upper[idx, cause]
  round(c(est, lo, hi) * 100, 2)
}

# Helper: extract 1-KM at milestone with CI
extract_km <- function(km_obj, t) {
  s <- summary(km_obj, times = t, extend = TRUE)
  est <- (1 - s$surv) * 100
  lo  <- (1 - s$upper) * 100
  hi  <- (1 - s$lower) * 100
  round(c(est, lo, hi), 2)
}

# ─────────────────────────────────────────────────────────────
# 1. RETREATMENT — Competing Risks CIF
# event_rn: 0=censored, 1=retreatment (event of interest), 2=death (competitor)
cat("Running CIF for retreatment (competing risks)...\n")
cif_retr <- cuminc(
  ftime   = ltfu$time_rn,
  fstatus = ltfu$event_rn,
  cencode = 0
)

# Extract at a given time from cmprsk cuminc object (cause "X 1" = retreatment)
get_cif_at <- function(cif, cause_name, t) {
  obj <- cif[[cause_name]]
  idx <- max(which(obj$time <= t))
  if (length(idx) == 0) return(c(est = NA, lo = NA, hi = NA))
  # cmprsk does not provide CI by default; compute Greenwood-based CI
  # Use est ± 1.96 * sqrt(var)
  est <- obj$est[idx]
  se  <- sqrt(obj$var[idx])
  c(est = round(est * 100, 2),
    lo  = round(max(0, est - 1.96 * se) * 100, 2),
    hi  = round(min(1, est + 1.96 * se) * 100, 2))
}

max_t <- max(ltfu$time_rn, na.rm = TRUE)
times_retr <- c(1, 5, 10, max_t)
labels_retr <- c("1 year", "5 years", "10 years", paste0("End of study (~", round(max_t, 1), " yr)"))

rows_r <- lapply(times_retr, function(t) get_cif_at(cif_retr, "1 1", t))
results_retr <- data.frame(
  Milestone = labels_retr,
  CIF_pct   = sapply(rows_r, `[[`, "est"),
  CI_lo_pct = sapply(rows_r, `[[`, "lo"),
  CI_hi_pct = sapply(rows_r, `[[`, "hi")
)
results_retr$Result <- sprintf("%.2f%% (95%% CI: %.2f%%\u2013%.2f%%)",
                                results_retr$CIF_pct, results_retr$CI_lo_pct, results_retr$CI_hi_pct)
cat("Retreatment CIF results:\n")
print(results_retr[, c("Milestone", "Result")])


# ─────────────────────────────────────────────────────────────
# 2. MORTALITY — Direct Cumulative Incidence (1 - KM)
# event_d: 0=censored/alive, 1=dead
cat("\nRunning 1-KM for mortality...\n")
km_dead <- survfit(Surv(time_d, event_d) ~ 1, data = ltfu)

max_t_d <- max(ltfu$time_d, na.rm = TRUE)
results_mort <- data.frame(
  Milestone = c("1 year", "5 years", "10 years", paste0("End of study (~", round(max_t_d, 1), " yr)")),
  time_yr   = c(milestones_yr, max_t_d)
)

vals_d <- t(sapply(results_mort$time_yr, extract_km, km_obj = km_dead))
results_mort$CumInc_pct <- vals_d[, 1]
results_mort$CI_lo_pct  <- vals_d[, 2]
results_mort$CI_hi_pct  <- vals_d[, 3]
results_mort$Result      <- sprintf("%.2f%% (95%% CI: %.2f%%–%.2f%%)",
                                     results_mort$CumInc_pct,
                                     results_mort$CI_lo_pct,
                                     results_mort$CI_hi_pct)

cat("Mortality cumulative incidence results:\n")
print(results_mort[, c("Milestone", "Result")])

# ─────────────────────────────────────────────────────────────
# 3. SAVE CSV for Python doc generation
write.csv(results_retr[, c("Milestone", "CIF_pct", "CI_lo_pct", "CI_hi_pct", "Result")],
          "ITT_Analysis/results/ci_retreatment.csv", row.names = FALSE)
write.csv(results_mort[, c("Milestone", "CumInc_pct", "CI_lo_pct", "CI_hi_pct", "Result")],
          "ITT_Analysis/results/ci_mortality.csv", row.names = FALSE)

cat("\nCSVs saved.\n")
