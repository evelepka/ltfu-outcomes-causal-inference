library(dplyr)
library(survival)
library(broom)
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE) %>%
  mutate(across(c(race_clean, edu_clean, dot_status, alcohol, drug_use, diabetes, hosp_admission, hiv_aids, clinical_clean, incarcerated, homelessness), ~ as.factor(.))) %>%
  mutate(sex = factor(sex, levels = c("Female", "Male")), age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "≥65"))) %>% 
  tidyr::drop_na(age_group, sex, race_clean, edu_clean, hiv_aids, diabetes, alcohol, drug_use, incarcerated, homelessness, hosp_admission, clinical_clean, dot_status) %>%
  filter(itt_group %in% c("Non-LTFU", "Loss to follow-up")) %>%
  mutate(is_ltfu = ifelse(itt_group == "Loss to follow-up", 1, 0)) %>%
  filter(time_d > 180 / 365.25) %>% mutate(time_d_lm = time_d - 180 / 365.25)

fit_int <- coxph(Surv(time_d_lm, event_d) ~ is_ltfu * age_group + sex + race_clean + edu_clean + hiv_aids + diabetes + alcohol + drug_use + incarcerated + homelessness + hosp_admission + clinical_clean + dot_status, data = df)
print(tidy(fit_int) %>% filter(grepl("is_ltfu", term)))

for (lvl in levels(df$age_group)) {
  df_sub <- df %>% filter(age_group == lvl)
  fit <- coxph(Surv(time_d_lm, event_d) ~ is_ltfu + sex + race_clean + edu_clean + hiv_aids + diabetes + alcohol + drug_use + incarcerated + homelessness + hosp_admission + clinical_clean + dot_status, data = df_sub)
  cat(lvl, ": HR =", exp(coef(fit)["is_ltfu"]), "\n")
}
