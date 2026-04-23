library(dplyr)
library(survival)
library(mice)
library(broom)

df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# Clean variables into factors
na_vals <- c("Missing", "Ignorado", "Unknown", "", "nan")
df <- df %>%
  mutate(across(c(race_clean, edu_clean, dot_status, alcohol, drug_use, diabetes, hosp_admission, hiv_aids, mental_health, tobacco_use, diagnosis_setting, lab_confirmed_stat), ~ ifelse(. %in% na_vals, NA, .)))

# Standardize Factors
df <- df %>%
  mutate(
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "≥65")),
    sex = as.factor(sex),
    race_clean = relevel(as.factor(race_clean), ref = "White"),
    edu_clean = relevel(as.factor(edu_clean), ref = "≥ 12 years"),
    hiv_aids = as.factor(hiv_aids),
    diabetes = as.factor(diabetes),
    alcohol = as.factor(alcohol),
    drug_use = as.factor(drug_use),
    mental_health = as.factor(mental_health),
    tobacco_use = as.factor(tobacco_use),
    other_immuno_condition = as.factor(other_immuno_condition),
    incarcerated = as.factor(incarcerated),
    homelessness = as.factor(homelessness),
    hosp_admission = as.factor(hosp_admission),
    clinical_clean = relevel(as.factor(clinical_clean), ref = "Pulmonary"),
    dot_status = as.factor(dot_status),
    diagnosis_setting = as.factor(diagnosis_setting),
    bac1_clean = as.factor(bac1_clean),
    sputum_culture_clean = as.factor(sputum_culture_clean),
    resistance_clean = as.factor(resistance_clean)
  )

covariates <- c("age_group", "sex", "race_clean", "edu_clean", "dot_status",
                "alcohol", "drug_use", "diabetes", "hosp_admission", "hiv_aids", "clinical_clean",
                "incarcerated", "homelessness", "mental_health", "tobacco_use", "other_immuno_condition",
                "diagnosis_setting", "bac1_clean", "sputum_culture_clean", "resistance_clean")

# 1. Setup MI Dataset
cat(sprintf("Setting up MICE targeting %d variables over %d records...\n", length(covariates), nrow(df)))

# Impute only the core covariates to save massive time.
df_mi <- df[, c("sinan_clean", covariates)]

# Use m=5, maxit=5
set.seed(42)
# fast PMM
imp <- mice(df_mi, m = 5, maxit = 5, method = "pmm", printFlag = FALSE)

# 2. Iterate Target Trials
res_list <- list()

# Reconstruct list of completed datasets to avoid repeated `complete()` logic inside the inner loop
imp_data <- list()
for (i in 1:5) {
  imp_data[[i]] <- complete(imp, i) %>%
    # Join back non-imputed outcomes and timing
    left_join(df[, c("sinan_clean", "itt_group", "best_start", "end_date", "time_d_tx", "event_d")], by="sinan_clean")
}

for (m in 1:6) {
  cat(sprintf("Running Trial %d...\n", m))
  start_yrs <- ((m - 1) * 30) / 365.25
  end_yrs   <- (m * 30) / 365.25
  
  cox_uv_fits <- list()
  for (i in 1:5) {
    df_i <- imp_data[[i]] %>% 
      filter(time_d_tx > start_yrs) %>%
      mutate(
        tx_duration_yrs = as.numeric(difftime(as.Date(end_date), as.Date(best_start), units="days")) / 365.25,
        eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
      ) %>%
      filter(eligible == 1) %>%
      mutate(
        expose = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs >= start_yrs & tx_duration_yrs < end_yrs, 1, 0),
        time_followup = time_d_tx - start_yrs
      )
    
    f_base <- as.formula(paste("Surv(time_followup, event_d) ~ expose +", paste(covariates, collapse = " + ")))
    cox_uv_fits[[i]] <- coxph(f_base, data = df_i)
  }
  
  # Pool utilizing Rubin's rules
  p_cox <- summary(pool(cox_uv_fits), exponentiate = TRUE, conf.int = TRUE)
  
  expose_row <- p_cox[p_cox$term == "expose", ]
  res_list[[m]] <- data.frame(
    Trial_Month = paste0("Month_", m),
    HR = expose_row$estimate,
    CI_Lower = expose_row$`2.5 %`,
    CI_Upper = expose_row$`97.5 %`,
    P_Value = expose_row$p.value
  )
}

final_res <- do.call(rbind, res_list)
print(final_res)

write.csv(final_res, "ITT_Analysis/results/target_trial_mi_6mo_array_hr.csv", row.names = FALSE)
cat("Successfully exported MI target trial arrays to ITT_Analysis/results/target_trial_mi_6mo_array_hr.csv\n")
