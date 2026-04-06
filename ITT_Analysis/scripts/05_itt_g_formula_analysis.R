# 05. ITT G-Formula Analysis (10-Year Mortality Only)
# ==============================================================================
# Purpose: Estimate 10-year causal risk of Death.
# Exposure: Abandonment vs. Non-Abandonment (Control).
# Baseline: Treatment Initiation (tx_start).
# Follow-up: 144 months (12 years).
# ==============================================================================

library(dplyr)
library(survival)
library(ggplot2)
library(lubridate)
library(tidyr)
library(splines)

# 1. Load Data
cat("Loading ITT cohort...\n")
df <- read.csv("Abandonment Paper/ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# Harmonized Recoding
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
    itt_group = factor(itt_group, levels = c("Non-LTFU", "Loss to follow-up")),
    diabetes = factor(diabetes, levels = c("No", "Yes")),
    alcohol = factor(alcohol, levels = c("No", "Yes")),
    drug_use = factor(drug_use, levels = c("No", "Yes")),
    incarcerated = factor(incarcerated, levels = c("No", "Yes")),
    homelessness = factor(homelessness, levels = c("No", "Yes")),
    hosp_admission = factor(hosp_admission, levels = c("No", "Yes"))
  )
}

df_c <- clean_data(df) %>% na.omit()
cat("N after cleaning/NA removal:", nrow(df_c), "\n")

# Convert years to months, cap at 120 (10 years)
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
cat("Running G-formula for 10-year Mortality...\n")
df_long <- prepare_long(df_c, "time_d_tx", "event_d")

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
write.csv(res_g, "Abandonment Paper/ITT_Analysis/results/g_formula_mortality_12y_results.csv", row.names = FALSE)

# 5. Plot (Standardized ggplot2)
cat("Generating Standardized Plot...\n")
curve_non_aug <- c(1, curve_non)
curve_aba_aug <- c(1, curve_aba)

plot_df <- data.frame(
  Year = (0:144) / 12,
  Non_LTFU = (1 - curve_non_aug) * 100,
  LTFU = (1 - curve_aba_aug) * 100
) %>%
pivot_longer(cols = c(Non_LTFU, LTFU), names_to = "group", values_to = "mortality")

png("Abandonment Paper/ITT_Analysis/results/g_formula_mortality_12y.png", width = 2400, height = 1800, res = 300)
ggplot(plot_df, aes(x = Year, y = mortality, color = group)) +
  geom_line(linewidth = 1.2) +
  theme_classic() +
  scale_color_manual(values = c("LTFU" = "red", "Non_LTFU" = "blue"), 
                     labels = c("Loss to follow-up", "Control (Non-LTFU)")) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, 1)) +
  labs(title = "Strategy 1: 12-Year Mortality from Treatment Initiation (Time 0)",
       subtitle = "Cumulative incidence estimates via G-formula",
       x = "Years from Treatment Initiation",
       y = "Cumulative Incidence of Death (%)",
       color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 12))
dev.off()

cat("ITT 12-Year Mortality G-Formula complete.\n")
