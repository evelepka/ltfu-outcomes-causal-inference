# 29. Post-Dropout Survival Analysis (Moment of Abandonment)
# ==============================================================================
# Purpose: Resolves immortal time bias by resetting Time Zero to the exact
#          day of abandonment. This cleanly isolates the biological survival
#          curves of Early vs Late abandoners from the moment they quit.
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

# Subset to strictly LTFU
df_ltfu <- clean_data(df) %>% filter(itt_group == "Loss to follow-up")

cat("\n--- 2. Resetting Time Zero to Moment of Abandonment ---\n")
df_ltfu_post <- df_ltfu %>%
  mutate(
    date_start = as.Date(best_start),
    date_end = as.Date(end_date),
    time_abandon_yrs = as.numeric(date_end - date_start) / 365.25,
    # Post-dropout survival time
    time_post_dropout = pmax(0.001, time_d - time_abandon_yrs),
    # Ensure they died AFTER abandonment 
    event_post = ifelse(event_d == 1 & time_d >= time_abandon_yrs, 1, 0)
  ) %>%
  filter(time_d >= time_abandon_yrs) 

df_ltfu_post$tx_month_grp <- factor(df_ltfu_post$tx_month_grp, levels = c("< 2 months", "2 to <4 months", "≥ 4 months"))

cat("Total LTFU cohort evaluated:", nrow(df_ltfu_post), "\n")
print(table(df_ltfu_post$tx_month_grp))

cat("\n--- 3. Fitting Post-Dropout Cox Model ---\n")
fit <- coxph(
  Surv(time_post_dropout, event_post) ~ tx_month_grp + age_group + sex + race_clean +
    edu_clean + hiv_aids + diabetes + alcohol + drug_use +
    incarcerated + homelessness + hosp_admission + clinical_clean + dot_status,
  data = df_ltfu_post
)

res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(grepl("tx_month_grp", term))
print(res %>% select(term, estimate, conf.low, conf.high, p.value))

write.csv(res, "ITT_Analysis/results/early_vs_late_post_dropout_hr.csv", row.names=FALSE)
cat("Results saved to ITT_Analysis/results/early_vs_late_post_dropout_hr.csv\n")
