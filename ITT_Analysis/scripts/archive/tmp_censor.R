# Temporary array to test non-proportional hazards by truncating follow up
library(dplyr)
library(survival)
library(broom)

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

trial_list <- list()
for (m in 1:6) {
  start_days <- (m - 1) * 30
  end_days   <- m * 30
  start_yrs <- start_days / 365.25
  end_yrs   <- end_days / 365.25
  
  df_trial <- df_c %>% filter(time_d > start_yrs) %>%
    mutate(
      eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
    ) %>% filter(eligible == 1) %>%
    mutate(
      expose = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs >= start_yrs & tx_duration_yrs < end_yrs, 1, 0),
      trial_month = paste0("Month_", m),
      time_followup = time_d - start_yrs
    )
  trial_list[[m]] <- df_trial
}

df_pooled <- bind_rows(trial_list)
df_pooled$trial_month <- factor(df_pooled$trial_month, levels = paste0("Month_", 1:6))

# Test windows
censor_windows <- c(1, 2, 3, 5, 20) # 20 is functionally infinite for this dataset
results <- list()

for (w in censor_windows) {
  df_w <- df_pooled %>%
    mutate(
      event_w = ifelse(time_followup > w, 0, event_d),
      time_w = ifelse(time_followup > w, w, time_followup)
    )
    
  fit <- coxph(
    Surv(time_w, event_w) ~ expose * trial_month + age_group + sex + race_clean +
      edu_clean + hiv_aids + diabetes + alcohol + drug_use +
      incarcerated + homelessness + hosp_admission + clinical_clean + dot_status,
    data = df_w, cluster = sinan_clean
  )
  
  res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(grepl("expose", term) & !grepl(":", term))
  results[[as.character(w)]] <- data.frame(FollowUp_Years = w, HR = res$estimate[1], P = res$p.value[1])
}

print(bind_rows(results))
