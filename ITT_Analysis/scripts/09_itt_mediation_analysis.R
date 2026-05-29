# 09. ITT Mediation Analysis (Time-Dependent Cox for Retreatment)
# ==============================================================================
# Purpose: Evaluate if Retreatment mediates the effect of Abandonment on Mortality.
# Cohort: Landmark ITT (survivors of index treatment).
# Method: Time-Dependent Cox Proportional Hazards Model.
# ==============================================================================

library(dplyr)
library(survival)

cat("Loading ITT cohort for Mediation Analysis...\n")
df <- read.csv("Abandonment Paper/ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# Apply Landmark Filter: Only patients who survived past `end_date`
df_landmark <- df %>% filter(time_d > 0)
cat("Landmark N (survived to end_date):", nrow(df_landmark), "\n")

# Clean Covariates
clean_data <- function(data) {
  na_vals <- c("Missing", "Ignorado", "Unknown", "", "nan")
  data %>% 
  mutate(across(c(race, education, supervised_therapy, alcohol, drug_use_stat, diabetes, hosp_admission, hiv_aids), ~ ifelse(. %in% na_vals, NA, .))) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "≥65")),
    hiv_aids = factor(hiv_aids, levels = c("Negative", "Positive")),
    race_clean = factor(race, levels = c("White", "Black or Mixed", "Other")),
    edu_clean = factor(education, levels = c("≥ 12 years", "8 - 11 years", "≤ 7 years", "None")),
    clinical_clean = factor(clinical_classif, levels = c("Pulmonary", "Extrapulmonary", "Pulmonary and Extrapulmonary or disseminated")),
    dot_status = factor(supervised_therapy, levels = c("No", "Yes")),
    itt_group = factor(itt_group, levels = c("Non-LTFU", "Loss to follow-up")),
    diabetes = factor(diabetes, levels = c("No", "Yes")),
    alcohol = factor(alcohol, levels = c("No", "Yes")),
    drug_use = factor(drug_use_stat, levels = c("No", "Yes")),
    incarcerated = factor(incarceration, levels = c("No", "Yes")),
    homelessness = factor(homelessness, levels = c("No", "Yes")),
    hosp_admission = factor(hosp_admission, levels = c("No", "Yes")),
    # Ensure time variables are properly bounded (time_rn can't be > time_d if they died)
    time_rn = ifelse(event_rn == 1 & time_rn > time_d, time_d, time_rn),
    event_rn = ifelse(event_rn == 1 & time_rn > time_d, 0, event_rn)
  )
}

df_c <- clean_data(df_landmark) %>% na.omit()

cat("N after cleaning/NA removal:", nrow(df_c), "\n")

# Model 1: Total Effect (Standard Cox)
cat("\n--- MODEL 1: TOTAL EFFECT ---\n")
cox_total <- coxph(
    Surv(time_d, event_d) ~ itt_group + age_group + sex + race_clean + edu_clean + 
                            hiv_aids + diabetes + alcohol + drug_use +
                            incarcerated + homelessness + hosp_admission + 
                            clinical_clean + dot_status,
    data = df_c
)
summary(cox_total)

# Prepare Time-Dependent Dataset using tmerge
cat("\nBuilding Time-Varying Dataset...\n")
# Start with base survival
df_tv <- tmerge(
    data1 = df_c,
    data2 = df_c,
    id = sinan_clean,
    death = event(time_d, event_d)
)

# Add time-dependent covariate for retreatment (value becomes 1 AFTER time_rn)
df_tv <- tmerge(
    data1 = df_tv,
    data2 = df_c,
    id = sinan_clean,
    retreatment_status = tdc(time_rn)
)

# Some safety checks on the new dataset
df_tv <- df_tv %>% filter(tstart < tstop)

# Model 2: Direct Effect (Time-Dependent Cox)
cat("\n--- MODEL 2: DIRECT EFFECT (Adj for Time-varying Retreatment) ---\n")
cox_direct <- coxph(
    Surv(tstart, tstop, death) ~ itt_group + retreatment_status + 
                                age_group + sex + race_clean + edu_clean + 
                                hiv_aids + diabetes + alcohol + drug_use +
                                incarcerated + homelessness + hosp_admission + 
                                clinical_clean + dot_status,
    data = df_tv
)
summary(cox_direct)

# Summarize Mediation
coef_tot <- coef(cox_total)["itt_groupLoss to follow-up"]
coef_dir <- coef(cox_direct)["itt_groupLoss to follow-up"]
hr_tot <- exp(coef_tot)
hr_dir <- exp(coef_dir)

prop_mediated <- (coef_tot - coef_dir) / coef_tot * 100

cat("\n================ MEDIATION SUMMARY ================\n")
cat(sprintf("Hazard Ratio (Total Effect): %.2f\n", hr_tot))
cat(sprintf("Hazard Ratio (Direct Effect): %.2f\n", hr_dir))
cat(sprintf("Hazard Ratio of Retreatment itself: %.2f\n", exp(coef(cox_direct)["retreatment_status"])))
cat(sprintf("Proportion of effect mediated by Retreatment: %.1f%%\n", prop_mediated))

# Save summary to a text file for reporting
sink("Abandonment Paper/ITT_Analysis/results/mediation_summary.txt")
cat("Retreatment Mediation Analysis (Time-Dependent Cox)\n")
cat("==================================================\n\n")
cat(sprintf("Total Effect HR (Abandonment -> Death): %.2f\n", hr_tot))
cat(sprintf("Direct Effect HR (Abandonment -> Death, adjusting for Retreatment): %.2f\n", hr_dir))
cat(sprintf("Effect of Retreatment (Retreatment -> Death): %.2f\n", exp(coef(cox_direct)["retreatment_status"])))
cat(sprintf("Proportion Mediated: %.1f%%\n", prop_mediated))
sink()

cat("Mediation script finished successfully.\n")
