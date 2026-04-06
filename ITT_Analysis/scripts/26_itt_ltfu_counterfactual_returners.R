# 26. Counterfactual Mortality for Returners
# ==============================================================================
# Purpose: Estimate what the mortality of "Returners" would have been if they 
#          had the clinical course of patients who never abandoned. 
#          This isolates the severity/relapse penalty.
# ==============================================================================

library(dplyr)
library(survival)

cat("\n--- 1. Loading ITT Cohort ---\n")
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

clean_data <- function(data) {
  na_vals <- c("Missing", "Ignorado", "Unknown", "", "nan")
  data %>% 
  mutate(across(c(race_clean, edu_clean, dot_status, alcohol, drug_use, diabetes, hosp_admission, hiv_aids, clinical_clean, incarcerated, homelessness), ~ as.factor(.))) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "≥65")),
    time_followup = pmin(10, pmax(0.001, time_d_tx))  # Standardize up to 10 years
  ) %>% tidyr::drop_na(age_group, sex, race_clean, edu_clean, hiv_aids, diabetes, alcohol, drug_use, incarcerated, homelessness, hosp_admission, clinical_clean, dot_status)
}

df_c <- clean_data(df)

# Split into Non-LTFU and Returners
df_non_ltfu <- df_c %>% filter(itt_group == "Non-LTFU")
df_returners <- df_c %>% filter(itt_group == "Loss to follow-up" & event_rn == 1)

cat("Reference Cohort (Non-LTFU):", nrow(df_non_ltfu), "\n")
cat("Target Cohort (Returners):", nrow(df_returners), "\n")


cat("\n--- 2. Fitting Baseline Cox Model on Non-LTFU ---\n")
# Model learns how baseline predictors map to death for people who finish treatment normally.
cox_non_ltfu <- coxph(
  Surv(time_followup, event_d) ~ age_group + sex + race_clean +
    edu_clean + hiv_aids + diabetes + alcohol + drug_use +
    incarcerated + homelessness + hosp_admission + clinical_clean + dot_status,
  data = df_non_ltfu,
  x = TRUE
)


cat("\n--- 3. Predicting Counterfactuals for Returners ---\n")
# If the Returners had followed the "Non-LTFU" trajectory, what would their survival be?
surv_predictions <- survfit(cox_non_ltfu, newdata = df_returners)

# Extract expected survival probability at 10 years
# Handling if exact time is not 10, extract max available
t_idx <- which.min(abs(surv_predictions$time - 10))
expected_surv <- surv_predictions$surv[t_idx, ]
expected_mort_mean <- (1 - mean(expected_surv)) * 100


cat("\n--- 4. Calculating Observed Reality for Returners ---\n")
# Calculate their actual observed KM mortality
km_returners <- survfit(Surv(time_followup, event_d) ~ 1, data = df_returners)
obs_t_idx <- which.min(abs(km_returners$time - 10))
observed_mort <- (1 - km_returners$surv[obs_t_idx]) * 100

cat("\n=======================================================\n")
cat("RESULTS: THE SYMPTOMATIC RELAPSE PENALTY\n")
cat("=======================================================\n")
cat(sprintf("Observed 10-Year Mortality for Returners:       %.2f%%\n", observed_mort))
cat(sprintf("Counterfactual Expected (had they not abandoned): %.2f%%\n", expected_mort_mean))
cat(sprintf("Absolute Mortality Penalty associated with Relapse: +%.2f%%\n", observed_mort - expected_mort_mean))
cat("=======================================================\n")
