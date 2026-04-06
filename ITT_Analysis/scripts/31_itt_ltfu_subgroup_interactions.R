# 31. Subgroup Effect Modification Analysis
# ==============================================================================
# Purpose: Identify if the penalty of abandoning treatment is structurally 
#          worse for specific vulnerability subgroups, evaluated on the strict 
#          180-Day Landmark Cohort to eliminate early survival bias.
# ==============================================================================

library(dplyr)
library(survival)
library(broom)

cat("\n--- 1. Loading and Cleaning Landmark Cohort ---\n")
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
  filter(itt_group %in% c("Non-LTFU", "Loss to follow-up")) %>%
  mutate(
    is_ltfu = ifelse(itt_group == "Loss to follow-up", 1, 0)
  )

lm_years <- 180 / 365.25
df_lm <- df_c %>%
  filter(time_d > lm_years) %>%
  mutate(
    time_d_lm = time_d - lm_years
  )

cat("180-Day Landmark Cohort size:", nrow(df_lm), "\n")


cat("\n--- 2. Defining Subgroup Iteration Logic ---\n")

subgroups <- c("hiv_aids", "homelessness", "incarcerated", "drug_use", "alcohol", "diabetes", "sex")
covariates <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids", "diabetes", "alcohol", "drug_use", "incarcerated", "homelessness", "hosp_admission", "clinical_clean", "dot_status")

results_list <- list()

for (sg in subgroups) {
  cat("\nProcessing Subgroup:", sg, "...\n")
  
  # 1. Interaction P-Value (Pooled Model)
  # Remove the subgroup from baseline covariates so we don't adjust for it twice
  adj_covs <- covariates[covariates != sg]
  form_int <- as.formula(paste("Surv(time_d_lm, event_d) ~ is_ltfu *", sg, "+", paste(adj_covs, collapse = " + ")))
  
  fit_int <- coxph(form_int, data = df_lm)
  res_int <- tidy(fit_int)
  # Find the interaction term row
  int_term <- res_int %>% filter(grepl(paste0("is_ltfu:", sg), term))
  
  # Because an m-level factor has m-1 interactions, we take the minimum p-value as a proxy 
  # or do a formal wald test, but extracting all interaction rows works for transparency
  min_p_int <- min(int_term$p.value)
  
  # 2. Extract Strata HRs (Stratified Models)
  levels_sg <- levels(df_lm[[sg]])
  
  for (lvl in levels_sg) {
    df_sub <- df_lm %>% filter(!!sym(sg) == lvl)
    if (nrow(df_sub) < 100) next # Skip practically empty strata
    
    # Fit independent fully adjusted Cox model inside this stratum
    form_strata <- as.formula(paste("Surv(time_d_lm, event_d) ~ is_ltfu +", paste(adj_covs, collapse = " + ")))
    fit_strata <- coxph(form_strata, data = df_sub)
    
    res_strata <- tidy(fit_strata, exponentiate = TRUE, conf.int = TRUE) %>% 
      filter(term == "is_ltfu")
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Subgroup = sg,
      Stratum = as.character(lvl),
      N_Stratum = nrow(df_sub),
      HR = res_strata$estimate,
      CI_Low = res_strata$conf.low,
      CI_High = res_strata$conf.high,
      P_Value = res_strata$p.value,
      Interaction_P_Value = ifelse(lvl == levels_sg[1], sprintf("(Ref) p=%.4f", min_p_int), sprintf("%.4f", min_p_int))
    )
  }
}

cat("\n--- 3. Finalizing Results ---\n")
final_results <- bind_rows(results_list)
print(final_results)

write.csv(final_results, "ITT_Analysis/results/ltfu_subgroup_effect_modification.csv", row.names = FALSE)
cat("\nResults saved to ITT_Analysis/results/ltfu_subgroup_effect_modification.csv\n")
