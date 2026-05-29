# 11. ITT Retreatment Timing Analysis
# ==============================================================================
# Purpose: Evaluate the effect of the "month of retreatment" on mortality.
# Cohort: Restrict to patients who had a retreatment event.
# Method: Cox Proportional Hazards Model from the time of retreatment.
# ==============================================================================

library(dplyr)
library(survival)
library(ggplot2)

cat("Loading ITT cohort for Retreatment Timing Analysis...\n")
df <- read.csv("Abandonment Paper/ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# Filter: Only patients who retreated and survived up to retreatment date
df_rn <- df %>% 
    filter(event_rn == 1) %>%
    mutate(post_rn_time = time_d - time_rn) %>%
    filter(post_rn_time > 0) # Must survive after returning for retreatment to be followed

cat("N of retreaters with follow-up:", nrow(df_rn), "\n")

# Categorize time to retreatment
df_rn <- df_rn %>% mutate(
    time_to_rn_cat = factor(case_when(
        time_rn <= 0.5 ~ "1. Early (0-6 months)",
        time_rn <= 1.0 ~ "2. Intermediate (6-12 months)",
        time_rn <= 3.0 ~ "3. Late (1-3 years)",
        TRUE ~ "4. Very Late (>3 years)"
    ), levels = c("4. Very Late (>3 years)", "3. Late (1-3 years)", "2. Intermediate (6-12 months)", "1. Early (0-6 months)")),
    # Ensure covariates are factors
    sex = factor(sex, levels = c("Female", "Male")),
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "65+")),
    hiv_aids = factor(hiv_aids, levels = c("Negative", "Positive")),
    race = factor(race, levels = c("White", "Black or Mixed", "Other")),
    education = factor(education, levels = c("≥ 12 years", "8 - 11 years", "≤ 7 years", "None")),
    clinical_classif = factor(clinical_classif, levels = c("Pulmonary", "Extrapulmonary", "Pulmonary and Extrapulmonary or disseminated")),
    supervised_therapy = factor(supervised_therapy, levels = c("No", "Yes")),
    itt_group = factor(itt_group, levels = c("Non-LTFU", "Loss to follow-up")),
    diabetes = factor(diabetes, levels = c("No", "Yes")),
    alcohol = factor(alcohol, levels = c("No", "Yes")),
    drug_use_stat = factor(drug_use_stat, levels = c("No", "Yes")),
    incarceration = factor(incarceration, levels = c("No", "Yes")),
    homelessness = factor(homelessness, levels = c("No", "Yes")),
    hosp_admission = factor(hosp_admission, levels = c("No", "Yes"))
) %>% na.omit()

cat("N after covariate cleaning/NA removal:", nrow(df_rn), "\n")

# Distribution of Retreatment Times
cat("\nDistribution of Retreatment Timing:\n")
print(table(df_rn$time_to_rn_cat))

# Multivariate Cox Model for Post-Retreatment Mortality
cat("\n--- MODEL: POST-RETREATMENT SURVIVAL ---\n")
cox_model <- coxph(
    Surv(post_rn_time, event_d) ~ time_to_rn_cat + itt_group + 
                                age_group + sex + race + education + 
                                hiv_aids + diabetes + alcohol + drug_use_stat +
                                incarceration + homelessness + hosp_admission + 
                                clinical_classif + supervised_therapy,
    data = df_rn
)
summary(cox_model)

# Save summary to text file
sink("Abandonment Paper/ITT_Analysis/results/retreatment_timing_summary.txt")
cat("Effect of Retreatment Timing on Post-Retreatment Mortality\n")
cat("========================================================\n\n")
cat("Reference Group for Timing: Very Late (>3 years)\n\n")
print(summary(cox_model)$coef[, c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")])
sink()

cat("\nRetreatment Timing script finished successfully.\n")
