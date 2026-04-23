# 28. Abandonment Gap Duration Analysis
# ==============================================================================
# Purpose: Evaluate how the duration of the "treatment gap" (abandonment up to retreatment)
#          predicts post-return mortality.
# ==============================================================================

library(dplyr)
library(survival)
library(broom)

cat("\n--- 1. Loading LTFU Cohort ---\n")
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

clean_data <- function(data) {
  data %>% 
  mutate(across(c(race_clean, edu_clean, dot_status, alcohol, drug_use, diabetes, hosp_admission, hiv_aids, clinical_clean, incarcerated, homelessness), ~ as.factor(.))) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
  ) %>% tidyr::drop_na(age_group, sex, race_clean, edu_clean, hiv_aids, diabetes, alcohol, drug_use, incarcerated, homelessness, hosp_admission, clinical_clean, dot_status)
}

# Subset strictly to LTFU who retreat
df_ltfu_rn <- clean_data(df) %>% filter(itt_group == "Loss to follow-up" & event_rn == 1)

cat("\n--- 2. Calculating Gap Duration ---\n")
# gap_duration = date_return - end_date
# return date = best_start + time_rn * 365.25
df_ltfu_rn <- df_ltfu_rn %>%
  mutate(
    date_start = as.Date(best_start),
    date_end = as.Date(end_date),
    date_return = date_start + (time_rn * 365.25),
    gap_days = as.numeric(date_return - date_end),
    gap_months = gap_days / (365.25 / 12)
  ) %>%
  filter(gap_months > 0) # ensure logical positive gap

df_ltfu_rn <- df_ltfu_rn %>%
  mutate(
    gap_cat = case_when(
      gap_months < 2 ~ "<2 months",
      gap_months >= 2 & gap_months <= 6 ~ "2-6 months",
      gap_months > 6 ~ ">6 months"
    )
  )

df_ltfu_rn$gap_cat <- factor(df_ltfu_rn$gap_cat, levels = c("<2 months", "2-6 months", ">6 months"))

cat("Cohort Size after Gap logic:", nrow(df_ltfu_rn), "\n")
print(table(df_ltfu_rn$gap_cat))

cat("\n--- 3. Fitting Prospective Survival From Re-Entry ---\n")
# Time Zero represents the moment of Retreatment.
# Survival time is (time_d - time_rn)
df_post_rn <- df_ltfu_rn %>%
  mutate(
    time_post_rn = pmax(0.001, time_d - time_rn),
    # If they are censored before retreatment technically event_d is 0 and time_d <= time_rn, handle that:
    event_post = ifelse(event_d == 1 & time_d >= time_rn, 1, 0)
  ) %>% filter(time_d >= time_rn)

fit <- coxph(
  Surv(time_post_rn, event_post) ~ gap_cat + age_group + sex + race_clean +
    edu_clean + hiv_aids + diabetes + alcohol + drug_use +
    incarcerated + homelessness + hosp_admission + clinical_clean + dot_status,
  data = df_post_rn
)

res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(grepl("gap_cat", term))
print(res %>% select(term, estimate, conf.low, conf.high, p.value))

write.csv(res, "ITT_Analysis/results/abandonment_duration_hr.csv", row.names=FALSE)
cat("Results saved to ITT_Analysis/results/abandonment_duration_hr.csv\n")
