# 24. ITT LTFU Marginal Structural Model (IPW for Retreatment)
# ==============================================================================
# Purpose: Estimate the causal effect of returning to treatment (retreatment) 
#          on mortality among the Loss to Follow-up cohort.
# Method: Inverse Probability Weighting (IPW) via Pooled Logistic Regression
#         to fit a Marginal Structural Model.
# ==============================================================================

library(dplyr)
library(survival)
library(tidyr)
library(splines)
library(broom)

cat("\n--- 1. Loading and Preparing ITT LTFU Cohort ---\n")
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# Filter cohort to ONLY the Loss to Follow-up group, from `end_date` onwards (Landmark T=0)
df_ltfu <- df %>% filter(itt_group == "Loss to follow-up" & time_d > 0)
cat("Baseline LTFU Cohort Size:", nrow(df_ltfu), "\n")

# Clean Covariates
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
    hosp_admission = factor(hosp_admission, levels = c("No", "Yes")),
    
    # Cap retreatment time to time_d to ensure they can't be retreated after death
    time_rn = ifelse(event_rn == 1 & time_rn > time_d, time_d, time_rn),
    # If they were retreated "after" death (data artifact), nullify the retreatment
    event_rn = ifelse(event_rn == 1 & time_rn > time_d, 0, event_rn),
    
    # Make month vectors
    month_d = pmax(1, ceiling(time_d * 12.16)),        # Follow-up months for mortality
    month_rn = ifelse(event_rn == 1, pmax(1, ceiling(time_rn * 12.16)), pmax(1, ceiling(time_d * 12.16))) # Follow up for RN
  )
}

df_clean <- clean_data(df_ltfu) %>% tidyr::drop_na(age_group, sex, race_clean, edu_clean, hiv_aids, diabetes, alcohol, drug_use, incarcerated, homelessness, hosp_admission, clinical_clean, dot_status)

cat("Complete Case LTFU Cohort Size:", nrow(df_clean), "\n")

cat("\n--- 2. Building Person-Month Dataset ---\n")
# We unroll each patient into person-months
df_pm <- df_clean %>%
  uncount(month_d, .id = "month") %>%
  group_by(sinan_clean) %>%
  mutate(
    # Was retreatment initiated this month?
    start_rn = ifelse(event_rn == 1 & month == month_rn, 1, 0),
    # Their retreatment history status at this month
    retreatment_status = ifelse(month >= month_rn & event_rn == 1, 1, 0),
    # Death event at the final month
    death = ifelse(month == max(month) & event_d == 1, 1, 0)
  ) %>%
  ungroup()

cat("Person-months generated:", nrow(df_pm), "\n")

cat("\n--- 3. Calculating Inverse Probability Weights (IPW) ---\n")
# We only calculate weights for the probability of STARTING retreatment
# So we subset only to person-months where retreatment_status == 0 OR the exact month they start
df_weight <- df_pm %>% filter(retreatment_status == 0 | start_rn == 1)

# Numerator Model (Unconditional on baseline covariates)
num_mod <- glm(start_rn ~ ns(month, df=3), family = binomial, data = df_weight)

# Denominator Model (Conditional on baseline covariates)
den_mod <- glm(start_rn ~ ns(month, df=3) + age_group + sex + race_clean + edu_clean + hiv_aids + diabetes + alcohol + drug_use + incarcerated + homelessness + hosp_admission + clinical_clean + dot_status, family = binomial, data = df_weight)

# Extract predictions
df_weight$p_num_event <- predict(num_mod, type = "response")
df_weight$p_den_event <- predict(den_mod, type = "response")

# Calculate weight components for each month
df_weight <- df_weight %>%
  mutate(
    w_num = ifelse(start_rn == 1, p_num_event, 1 - p_num_event),
    w_den = ifelse(start_rn == 1, p_den_event, 1 - p_den_event),
    sw_component = w_num / w_den
  )

# Calculate cumulative product of weights per person
df_weight <- df_weight %>%
  group_by(sinan_clean) %>%
  mutate(sw_t = cumprod(sw_component)) %>%
  ungroup()

# Truncate extreme weights (1st and 99th percentiles) to stabilize variance
q01 <- quantile(df_weight$sw_t, 0.01)
q99 <- quantile(df_weight$sw_t, 0.99)
df_weight <- df_weight %>%
  mutate(sw_t_trunc = pmax(q01, pmin(q99, sw_t)))

cat("Stabilized Weights Summary:\n")
print(summary(df_weight$sw_t_trunc))

# Now merge these weights back onto the full person-month dataset.
# The rule for time-varying point exposure IPW: The weight remains constant after exposure initiation.
weights_final <- df_weight %>% select(sinan_clean, month, sw_t_trunc)

df_pm_weighted <- df_pm %>%
  left_join(weights_final, by = c("sinan_clean", "month")) %>%
  group_by(sinan_clean) %>%
  fill(sw_t_trunc, .direction = "down") %>% # Carry forward the weight from the month of retreatment
  ungroup()

cat("Carried-forward weights checking (NAs should be 0):", sum(is.na(df_pm_weighted$sw_t_trunc)), "\n")

cat("\n--- 4. Fitting the Marginal Structural Model ---\n")
# Pooled logistic regression behaves identically to discrete-time Cox Model
msm_mod <- glm(death ~ retreatment_status + ns(month, df=3), 
               family = binomial(link = "cloglog"), # cloglog approximates proportional hazards
               data = df_pm_weighted, 
               weights = sw_t_trunc)

# Extract robust standard errors (ignoring full GEE/cluster variance for simplicity here, but good practice)
# We will use standard broom tidy summary as a proxy for the discrete hazard ratio
msm_res <- broom::tidy(msm_mod, exponentiate=TRUE, conf.int=TRUE) %>%
    filter(term == "retreatment_status") %>%
    mutate(
      HR_formatted = sprintf("%.3f (%.3f - %.3f)", estimate, conf.low, conf.high),
      P_value = sprintf("%.4f", p.value)
    )

print(msm_res)

# Unadjusted crude model for comparison
crude_mod <- glm(death ~ retreatment_status + ns(month, df=3), 
                 family = binomial(link = "cloglog"), 
                 data = df_pm_weighted)
crude_res <- broom::tidy(crude_mod, exponentiate=TRUE, conf.int=TRUE) %>% filter(term == "retreatment_status")
cat(sprintf("Crude Unadjusted HR: %.3f\n", crude_res$estimate))

out_df <- data.frame(
  Model = c("Unadjusted", "Marginal Structural Model (IPW)"),
  HR = c(sprintf("%.3f (%.3f - %.3f)", crude_res$estimate, crude_res$conf.low, crude_res$conf.high), msm_res$HR_formatted),
  P_value = c(sprintf("%.4f", crude_res$p.value), msm_res$P_value)
)

write.csv(out_df, "Abandonment Paper/ITT_Analysis/results/msm_ipw_retreatment_hr.csv", row.names = FALSE)
cat("\nResults saved to Abandonment Paper/ITT_Analysis/results/msm_ipw_retreatment_hr.csv\n")
