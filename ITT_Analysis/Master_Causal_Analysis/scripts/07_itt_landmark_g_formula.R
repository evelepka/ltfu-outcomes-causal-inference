# 07. ITT Landmark G-Formula Analysis (12-Year Mortality Post-Treatment)
# ==============================================================================
# Purpose: Estimate 12-year causal risk of Death AFTER treatment completion/abandonment.
# Exposure: Abandonment vs. Non-Abandonment (Control).
# Baseline: Treatment End Date (end_date).
# Follow-up: 144 months (12 years).
# Note: Removes immortal time bias by excluding patients who died during treatment.
# ==============================================================================

library(dplyr)
library(survival)
library(ggplot2)
library(lubridate)
library(tidyr)
library(splines)

# 1. Load Data
cat("Loading ITT cohort for Landmark Analysis...\n")
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# Harmonized Recoding
clean_data <- function(data) {
  na_vals <- c("Missing", "Ignorado", "Unknown", "", "nan")
  data %>% 
  mutate(across(c(race, education, supervised_therapy, alcohol, drug_use_stat, diabetes, hosp_admission, hiv_aids), ~ ifelse(. %in% na_vals, NA, .))) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "≥65")),
    hiv_aids = factor(hiv_aids, levels = c("Negative", "Positive")),
    race_clean = factor(race, levels = c("White", "Black or Mixed", "Other")),
    edu_clean = factor(education, levels = c("≥ 12 years", "8 - 11 years", "≤ 7 years", "None")),
    clinical_clean = factor(clinical_classif, levels = c("Pulmonary", "Extrapulmonary", "Pulmonary and Extrapulmonary or disseminated")),
    dot_status = factor(supervised_therapy, levels = c("No", "Yes")),
    itt_group = factor(itt_group, levels = c("Non-LTFU", "Loss to follow-up")),
    diabetes = factor(diabetes, levels = c("No", "Yes")),
    alcohol = factor(alcohol, levels = c("No", "Yes")),
    drug_use = factor(drug_use_stat, levels = c("No", "Yes")),
    incarcerated = factor(incarceration, levels = c("No", "Yes")),
    homelessness = factor(homelessness, levels = c("No", "Yes")),
    hosp_admission = factor(hosp_admission, levels = c("No", "Yes"))
  )
}

# Apply Landmark Filter: Only patients who survived past `end_date`
df_landmark <- df %>% filter(time_d > 0)
cat("Original N:", nrow(df), "\n")
cat("Landmark N (survived to end_date):", nrow(df_landmark), "\n")

df_c <- clean_data(df_landmark) %>% na.omit()
cat("N after cleaning/NA removal:", nrow(df_c), "\n")

# Convert years to months, cap at 144 (12 years)
prepare_long <- function(data, time_col, event_col) {
  data %>%
    mutate(
      months_to_event = pmin(144, pmax(1, round(.data[[time_col]] * 12))),
      event_status = as.numeric(as.character(.data[[event_col]]))
    ) %>%
    uncount(months_to_event, .id = "month") %>%
    group_by(sinan_clean) %>%
    mutate(
      death = ifelse(month == max(month) & event_status == 1, 1, 0)
    ) %>%
    ungroup()
}

# 2. Mortality Analysis
cat("Running Landmark G-formula for 12-year Mortality...\n")
# USE time_d (time from end_date) instead of time_d_tx (time from tx_start)
df_long <- prepare_long(df_c, "time_d", "event_d")

fit_death <- glm(
  death ~ itt_group * ns(month, df = 3) + age_group + sex + race_clean +
    edu_clean + hiv_aids + diabetes + alcohol + drug_use +
    incarcerated + homelessness + hosp_admission + clinical_clean + dot_status,
  data = df_long, family = binomial
)

# 3. Simulate Counterfactual Curves
simulate_survival <- function(model, baseline_data, target_group) {
  months <- 1:144
  n <- nrow(baseline_data)
  curr_surv <- rep(1, n)
  pop_surv <- rep(1, 144)
  
  cf_data <- baseline_data %>% mutate(itt_group = target_group)
  
  for (m in months) {
    cf_data$month <- m
    p_death <- predict(model, newdata = cf_data, type = "response")
    curr_surv <- curr_surv * (1 - p_death)
    pop_surv[m] <- mean(curr_surv)
    if (m %% 24 == 0) cat("Month", m, "complete...\n")
  }
  return(pop_surv)
}

cat("Simulating Non-LTFU...\n")
curve_non <- simulate_survival(fit_death, df_c, "Non-LTFU")
cat("Simulating Loss to follow-up...\n")
curve_aba <- simulate_survival(fit_death, df_c, "Loss to follow-up")

# 4. Results at Milestones
milestones <- c(12, 60, 120, 144)
labels <- c("1 Year", "5 Years", "10 Years", "12 Years")

res_g <- data.frame(
  Time = labels,
  Month = milestones,
  CI_Non_LTFU = 1 - curve_non[milestones],
  CI_LTFU = 1 - curve_aba[milestones]
) %>%
  mutate(
    Risk_Diff = CI_LTFU - CI_Non_LTFU,
    Risk_Ratio = CI_LTFU / CI_Non_LTFU
  )

print(res_g)

# Save
write.csv(res_g, "ITT_Analysis/results/g_formula_mortality_landmark_results.csv", row.names = FALSE)

# 5. Plot (X-axis in Years)
png("ITT_Analysis/results/g_formula_mortality_landmark.png", width = 2400, height = 1800, res = 300)
years <- (0:144) / 12
curve_non_aug <- c(1, curve_non)
curve_aba_aug <- c(1, curve_aba)

plot(years, curve_non_aug, type = "l", col = "blue", lwd = 2, ylim = c(0.85, 1),
     xlab = "Years from Treatment Outcome (Landmark T=0)", ylab = "Counterfactual Survival",
     main = "12-Year Mortality Post-Treatment (Landmark G-Formula)",
     xaxt = "n")
lines(years, curve_aba_aug, col = "red", lwd = 2)
axis(1, at = seq(0, 12, by = 2))
abline(v = milestones/12, lty = 3, col = "gray70")
legend("bottomleft", legend = c("Non-LTFU (Control)", "Loss to follow-up"), col = c("blue", "red"), lwd = 2, bty = "n")
dev.off()

cat("ITT 12-Year Landmark Mortality G-Formula complete.\n")
