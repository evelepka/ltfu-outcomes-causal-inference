# 30. Sequential Target Trial Emulation (Rolling Cohorts) - 6 Month Array
# ==============================================================================
# Purpose: Compares early to late abandoners sequentially against patients who 
#          stayed on treatment for the exact same duration. Eliminates immortal 
#          time bias without landmarking. 
# ==============================================================================

library(dplyr)
library(survival)
library(broom)

cat("\n--- 1. Loading Cohort & Cleaning Covariates ---\n")
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

clean_data <- function(data) {
  data %>% 
  mutate(across(c(race_clean, edu_clean, dot_status, alcohol, drug_use, diabetes, hosp_admission, hiv_aids, clinical_clean, incarcerated, homelessness), ~ as.factor(.))) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
  ) %>% tidyr::drop_na(age_group, sex, race_clean, edu_clean, hiv_aids, diabetes, alcohol, drug_use, incarcerated, homelessness, hosp_admission, clinical_clean, dot_status)
}

df_c <- clean_data(df) %>%
  mutate(
    date_start = as.Date(best_start),
    date_end = as.Date(end_date),
    tx_duration_yrs = as.numeric(date_end - date_start) / 365.25
  )

cat("\n--- 2. Building 6-Month Sequential Array ---\n")

trial_list <- list()

for (m in 1:6) {
  start_days <- (m - 1) * 30
  end_days   <- m * 30
  
  start_yrs <- start_days / 365.25
  end_yrs   <- end_days / 365.25
  
  # Eligibility: Alive and still on treatment at 'start' 
  df_trial <- df_c %>% filter(time_d_tx > start_yrs) %>%
    mutate(
      eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
    ) %>% filter(eligible == 1) %>%
    mutate(
      # Exposed: Abandoned between 'start' and 'end'
      expose = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs >= start_yrs & tx_duration_yrs < end_yrs, 1, 0),
      trial_month = paste0("Month_", m),
      time_followup = time_d_tx - start_yrs,
      event_d = ifelse(time_followup > 2.0, 0, event_d),
      time_followup = ifelse(time_followup > 2.0, 2.0, time_followup)
    )
  
  trial_list[[m]] <- df_trial
  cat(sprintf("Trial %d (Days %d-%d) | N = %d\n", m, start_days, end_days, nrow(df_trial)))
}

cat("\n--- 3. Fitting Independent Target Trials ---\n")
results <- list()

for (m in 1:6) {
  df_trial <- trial_list[[m]]
  # Calculate exact target trial HR independently to guarantee exact 95% CIs
  fit <- coxph(
    Surv(time_followup, event_d) ~ expose + age_group + sex + race_clean +
      edu_clean + hiv_aids + diabetes + alcohol + drug_use +
      incarcerated + homelessness + hosp_admission + clinical_clean + dot_status,
    data = df_trial
  )
  
  res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "expose")
  
  results[[m]] <- data.frame(
    Trial_Month = paste0("Month_", m),
    HR = res$estimate,
    CI_Lower = res$conf.low,
    CI_Upper = res$conf.high,
    P_Value = res$p.value
  )
}

final_results <- bind_rows(results)
print(final_results)

write.csv(final_results, "ITT_Analysis/results/target_trial_6mo_array_hr.csv", row.names=FALSE)
cat("Results saved to ITT_Analysis/results/target_trial_6mo_array_hr.csv\n")
