library(dplyr)
library(survival)
library(ggplot2)
library(lubridate)
library(tidyr)
library(splines)

df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

clean_data <- function(data) {
  data %>% 
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "≥65")),
    hiv_aids = factor(hiv_aids, levels = c("Negative", "Positive")),
    race_clean = factor(race_clean, levels = c("White", "Black or Mixed", "Other")),
    edu_clean = factor(edu_clean, levels = c("≥ 12 years", "8 - 11 years", "≤ 7 years", "None")),
    clinical_clean = factor(clinical_clean, levels = c("Pulmonary", "Extrapulmonary", "Pulmonary and Extrapulmonary or disseminated")),
    dot_status = factor(dot_status, levels = c("No", "Yes")),
    diabetes = factor(diabetes, levels = c("No", "Yes")),
    alcohol = factor(alcohol, levels = c("No", "Yes")),
    drug_use = factor(drug_use, levels = c("No", "Yes")),
    incarcerated = factor(incarcerated, levels = c("No", "Yes")),
    homelessness = factor(homelessness, levels = c("No", "Yes")),
    hosp_admission = factor(hosp_admission, levels = c("No", "Yes"))
  )
}

df_c <- clean_data(df) %>% tidyr::drop_na(age_group, sex, race_clean, edu_clean, hiv_aids, diabetes, alcohol, drug_use, incarcerated, homelessness, hosp_admission, clinical_clean, dot_status)

df_c$date_start <- as.Date(df_c$best_start)
df_c$date_end   <- as.Date(df_c$end_date)
df_c$tx_duration_yrs <- as.numeric(difftime(df_c$date_end, df_c$date_start, units="days")) / 365.25

# Setup Month 2 Trial
start_yrs <- 30 / 365.25
end_yrs   <- 60 / 365.25

df_trial <- df_c %>% filter(time_d_tx > start_yrs) %>%
  mutate(eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)) %>%
  filter(eligible == 1) %>%
  mutate(
    expose = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs >= start_yrs & tx_duration_yrs < end_yrs, 1, 0),
    time_followup = time_d_tx - start_yrs
  )

# prepare long for G formula (cap at 120 months / 10 yrs)
df_trial <- df_trial %>%
  mutate(
    months_to_event = pmin(120, pmax(1, ceiling(time_followup * 12)))
  ) %>%
  uncount(months_to_event, .id = "month") %>%
  group_by(sinan_clean) %>%
  mutate(death = ifelse(month == max(month) & event_d == 1, 1, 0)) %>%
  ungroup()

fit_death <- glm(
  death ~ expose * ns(month, df = 3) + age_group + sex + race_clean +
    edu_clean + hiv_aids + diabetes + alcohol + drug_use +
    incarcerated + homelessness + hosp_admission + clinical_clean + dot_status,
  data = df_trial, family = binomial(link="cloglog")
)

simulate_survival <- function(model, baseline_data, target_expose) {
  months <- 1:120
  n <- nrow(baseline_data)
  curr_surv <- rep(1, n)
  pop_surv <- rep(1, 120)
  cf_data <- baseline_data %>% mutate(expose = target_expose)
  for (m in months) {
    cf_data$month <- m
    p_death <- predict(model, newdata = cf_data, type = "response")
    curr_surv <- curr_surv * (1 - p_death)
    pop_surv[m] <- mean(curr_surv)
  }
  return(pop_surv)
}

base_snapshot <- df_trial %>% group_by(sinan_clean) %>% slice(1) %>% ungroup()
cat("Simulating Maintained Care...\n")
curve_non <- simulate_survival(fit_death, base_snapshot, 0)
cat("Simulating Abandonment...\n")
curve_aba <- simulate_survival(fit_death, base_snapshot, 1)

plot_df <- data.frame(
  Year = (1:120) / 12,
  Maintained_Care = (1 - curve_non) * 100,
  Abandoned_Treatment = (1 - curve_aba) * 100
)

# Export
write.csv(plot_df, "ITT_Analysis/results/g_formula_target_trial_m2.csv", row.names=FALSE)
cat(sprintf("Year 5 Maintained: %.2f%% | Abandoned %.2f%%\n", plot_df$Maintained_Care[60], plot_df$Abandoned_Treatment[60]))
cat(sprintf("Year 10 Maintained: %.2f%% | Abandoned %.2f%%\n", plot_df$Maintained_Care[120], plot_df$Abandoned_Treatment[120]))
