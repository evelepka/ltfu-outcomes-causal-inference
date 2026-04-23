# 14. ITT Landmark 180-Day Mortality Comparison
# ==============================================================================
# Purpose: Compare 12-year mortality between LTFU and Control (Non-LTFU)
# Method: 180-day Landmark to mitigate immortal time bias.
# Adjusted for harmonized covariates.
# ==============================================================================

library(dplyr)
library(survival)
library(ggplot2)
library(splines)
library(tidyr)

# 1. Load Data
cat("Loading ITT cohort for Landmark Analysis...\n")
df <- read.csv("Abandonment Paper/ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# 2. Apply 180-day Landmark (0.493 years)
LANDMARK_DAYS <- 180
LANDMARK_YEARS <- LANDMARK_DAYS / 365.25

cat("Applying 180-day Landmark filter...\n")
# Keep only those who survived/were followed at least 180 days from tx_start
df_land <- df %>% filter(time_d_tx > LANDMARK_YEARS)

cat("Original N:", nrow(df), "\n")
cat("Landmark N (survived to day 180):", nrow(df_land), "\n")

# 3. Clean and Recode
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

clean_for_analysis <- function(data) {
  data %>% 
  mutate(
    # Landmark time (reset T=0)
    time_land = time_d_tx - LANDMARK_YEARS
  ) %>%
  clean_data() %>%
  filter(!is.na(itt_group)) %>%
  na.omit()
}

df_c <- clean_for_analysis(df_land)
cat("N after cleaning:", nrow(df_c), "\n")

# 4. Standard Cox Model (Adjusted HR)
cat("\n--- Running Adjusted Cox Landmark Model ---\n")
cox_fit <- coxph(
  Surv(time_land, event_d) ~ itt_group + age_group + sex + race_clean +
    edu_clean + hiv_aids + diabetes + alcohol + drug_use +
    incarcerated + homelessness + hosp_admission + clinical_clean + dot_status,
  data = df_c
)
print(summary(cox_fit))

# 5. G-Formula Analysis (Absolute Risk)
cat("\n--- Running Landmark G-formula (Absolute Risks) ---\n")

# Prepare long data for G-formula
prepare_long_g <- function(data) {
  data %>%
    mutate(
      months_to_event = pmin(144, pmax(1, round(time_land * 12))),
      evt = as.numeric(event_d)
    ) %>%
    uncount(months_to_event, .id = "month") %>%
    group_by(sinan_clean) %>%
    mutate(
      y = ifelse(month == max(month) & evt == 1, 1, 0)
    ) %>%
    ungroup()
}

df_long <- prepare_long_g(df_c)

# Fit pooled logistic model
fit_g <- glm(
  y ~ itt_group * ns(month, df = 3) + age_group + sex + race_clean +
    edu_clean + hiv_aids + diabetes + alcohol + drug_use +
    incarcerated + homelessness + hosp_admission + clinical_clean + dot_status,
  data = df_long, family = binomial
)

# Simulation helper
simulate_pop_risk <- function(model, baseline_data, intervention) {
  months <- 1:138 # 144 - 6 = 138 months of follow-up post-landmark
  n <- nrow(baseline_data)
  surv <- rep(1, n)
  risks <- rep(0, length(months))
  
  sim_data <- baseline_data %>% mutate(itt_group = intervention)
  
  for (m in months) {
    sim_data$month <- m
    p_death <- predict(model, newdata = sim_data, type = "response")
    surv <- surv * (1 - p_death)
    risks[m] <- 1 - mean(surv)
  }
  return(risks)
}

cat("Simulating counterfactuals...\n")
risk_non <- simulate_pop_risk(fit_g, df_c, "Non-LTFU")
risk_ltfu <- simulate_pop_risk(fit_g, df_c, "Loss to follow-up")

# Results at milestones (1, 5, 10 years post-landmark)
ml <- c(6, 54, 114) # +6m landmark = 1yr, 5yr, 10yr milestones from start
res_g <- data.frame(
  Year_Post_Start = c(1, 5, 10),
  Month_Post_Landmark = ml,
  Risk_Non_LTFU = risk_non[ml],
  Risk_LTFU = risk_ltfu[ml]
) %>%
mutate(RD = Risk_LTFU - Risk_Non_LTFU, RR = Risk_LTFU / Risk_Non_LTFU)

print(res_g)

# 6. Save results and Plot
write.csv(res_g, "Abandonment Paper/ITT_Analysis/results/g_formula_mortality_landmark_results.csv", row.names = FALSE)

cat("Generating Plot...\n")
# Create a data frame for plotting (full curves)
# Extract cumulative incidence for both groups
plot_df <- data.frame(
  Month = 1:138,
  Non_LTFU = risk_non,
  LTFU = risk_ltfu
) %>%
pivot_longer(cols = c(Non_LTFU, LTFU), names_to = "group", values_to = "risk")

png("Abandonment Paper/ITT_Analysis/results/landmark_180d_mortalidade_curves.png", width = 2400, height = 1800, res = 300)
ggplot(plot_df, aes(x = (Month + 6)/12, y = risk * 100, color = group)) +
  geom_line(linewidth = 1.2) +
  theme_classic() +
  scale_color_manual(values = c("LTFU" = "red", "Non_LTFU" = "blue"), 
                     labels = c("Loss to follow-up", "Control (Non-LTFU)")) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, 1)) +
  labs(title = "Strategy 2: 12-Year Mortality Post-Landmark (180 days)",
       subtitle = "Estimated via G-formula from the ITT Cohort (N=172,463)",
       x = "Years from Treatment Start",
       y = "Cumulative Incidence of Death (%)",
       color = "Group") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 12))
dev.off()

cat("\nDone. Results saved to ITT_Analysis/results/landmark_180d_mortalidade_curves.png\n")
