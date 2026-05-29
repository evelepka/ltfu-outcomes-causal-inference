# 27. Abandonment Timing 180-Day Landmark
# ==============================================================================
# Purpose: Evaluate long-term mortality risk based on WHEN patients abandoned
#          treatment, compared to those who perfectly completed.
#          Uses 180-day landmarking to eliminate early bias.
# ==============================================================================

library(dplyr)
library(survival)
library(broom)
library(ggplot2)

cat("\n--- 1. Loading ITT Cohort ---\n")
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

clean_data <- function(data) {
  data %>% 
  mutate(across(c(race_clean, edu_clean, dot_status, alcohol, drug_use, diabetes, hosp_admission, hiv_aids, clinical_clean, incarcerated, homelessness), ~ as.factor(.))) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
  ) %>% tidyr::drop_na(age_group, sex, race_clean, edu_clean, hiv_aids, diabetes, alcohol, drug_use, incarcerated, homelessness, hosp_admission, clinical_clean, dot_status)
}

df_c <- clean_data(df)

# Construct Exposure Variable
df_c <- df_c %>%
  mutate(
    exposure_grp = case_when(
      itt_group == "Non-LTFU" ~ "Non-LTFU (Completed)",
      itt_group == "Loss to follow-up" & tx_month_grp == "< 2 months" ~ "Abandoned <2 Mo",
      itt_group == "Loss to follow-up" & tx_month_grp == "2 to <4 months" ~ "Abandoned 2-4 Mo",
      itt_group == "Loss to follow-up" & tx_month_grp == "≥ 4 months" ~ "Abandoned >4 Mo",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(exposure_grp))

df_c$exposure_grp <- factor(df_c$exposure_grp, levels = c("Non-LTFU (Completed)", "Abandoned <2 Mo", "Abandoned 2-4 Mo", "Abandoned >4 Mo"))

cat("\n--- 2. 180-Day Landmark Subsetting ---\n")
lm_days <- 180
lm_years <- lm_days / 365.25

df_lm <- df_c %>%
  filter(time_d > lm_years) %>%
  mutate(
    time_d_lm = time_d - lm_years
  )

cat("Original cohort pre-LM:", nrow(df_c), "\n")
cat("Cohort surviving to Day 180:", nrow(df_lm), "\n")


cat("\n--- 3. Fitting Landmark Cox Model For Timing ---\n")
fit <- coxph(
  Surv(time_d_lm, event_d) ~ exposure_grp + age_group + sex + race_clean +
    edu_clean + hiv_aids + diabetes + alcohol + drug_use +
    incarcerated + homelessness + hosp_admission + clinical_clean + dot_status,
  data = df_lm
)

res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(grepl("exposure_grp", term))
print(res %>% select(term, estimate, conf.low, conf.high, p.value))

write.csv(res, "ITT_Analysis/results/abandonment_timing_hr.csv", row.names=FALSE)
cat("Results saved to ITT_Analysis/results/abandonment_timing_hr.csv\n")
