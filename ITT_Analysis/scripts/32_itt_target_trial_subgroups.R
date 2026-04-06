# 32. Sequential Target Trial Emulation - Subgroup Interactions
# ==============================================================================
# Purpose: Re-runs the subgroup effect modification mapping, but entirely 
#          inside the pooled 6-Month Sequential Target Trial Emulation framework 
#          to fully eliminate immortal time bias from the subgroup findings.
# ==============================================================================

library(dplyr)
library(survival)
library(broom)

cat("\n--- 1. Loading Cohort & Cleaning Covariates ---\n")
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
  mutate(
    date_start = as.Date(best_start),
    date_end = as.Date(end_date),
    tx_duration_yrs = as.numeric(date_end - date_start) / 365.25
  )

cat("\n--- 2. Building 6-Month Sequential Array ---\n")
trial_list <- list()
for (m in 1:6) {
  start_days <- (m - 1) * 30
  end_days   <- m * 30
  start_yrs <- start_days / 365.25
  end_yrs   <- end_days / 365.25
  
  df_trial <- df_c %>% filter(time_d > start_yrs) %>%
    mutate(
      eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
    ) %>% filter(eligible == 1) %>%
    mutate(
      expose = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs >= start_yrs & tx_duration_yrs < end_yrs, 1, 0),
      trial_month = paste0("Month_", m),
      time_followup = time_d - start_yrs
    )
  trial_list[[m]] <- df_trial
}

df_pooled <- bind_rows(trial_list)
df_pooled$trial_month <- factor(df_pooled$trial_month, levels = paste0("Month_", 1:6))
covariates <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids", "diabetes", 
                "alcohol", "drug_use", "incarcerated", "homelessness", "hosp_admission", 
                "clinical_clean", "dot_status", "trial_month")

cat("\n--- 3. Testing Interactions on Target Trials ---\n")
subgroups <- c("age_group", "sex", "hiv_aids", "homelessness")
results_interactions <- list()
stratified_hrs <- list()

for (sg in subgroups) {
  cat(sprintf("\n[Evaluating Target Trial Subgroup: %s]\n", sg))
  
  # Remove the subgroup from the base covariates to prevent duplication
  base_covs_adj <- setdiff(covariates, sg)
  rhs_interaction <- paste("expose *", sg, "+", paste(base_covs_adj, collapse = " + "))
  f_int <- as.formula(paste("Surv(time_followup, event_d) ~", rhs_interaction))
  
  cat("Fitting pooled interaction model...\n")
  fit_int <- coxph(f_int, data = df_pooled, cluster = sinan_clean)
  
  res_int <- tidy(fit_int, exponentiate = FALSE, conf.int = FALSE)
  # Grab the interaction p-values
  int_terms <- res_int %>% filter(grepl(paste0("expose:", sg), term))
  
  p_val_str <- "Not calculated properly"
  if (nrow(int_terms) > 0) {
      if (nrow(int_terms) == 1) {
          p_val_str <- sprintf("p = %.4f", int_terms$p.value[1])
      } else {
          min_p <- min(int_terms$p.value, na.rm=TRUE)
          p_val_str <- sprintf("min p = %.4f", min_p)
      }
  }
  results_interactions[[sg]] <- p_val_str
  
  # Stratified Model Extraction
  cat("Fitting stratified models...\n")
  lvls <- unique(as.character(df_pooled[[sg]]))
  lvls <- lvls[!is.na(lvls)]
  
  for (lvl in lvls) {
    df_sub <- df_pooled[df_pooled[[sg]] == lvl, ]
    if (nrow(df_sub) == 0) next
    
    rhs_strat <- paste("expose +", paste(base_covs_adj, collapse = " + "))
    f_strat <- as.formula(paste("Surv(time_followup, event_d) ~", rhs_strat))
    
    fit_strat <- NULL
    try({
       fit_strat <- coxph(f_strat, data = df_sub, cluster = sinan_clean)
    }, silent = TRUE)
    
    if (!is.null(fit_strat)) {
        res_strat <- tidy(fit_strat, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "expose")
        if(nrow(res_strat)>0){
            stratified_hrs[[paste(sg, lvl, sep="_")]] <- data.frame(
                Subgroup = sg,
                Level = lvl,
                HR = res_strat$estimate[1],
                CI_L = res_strat$conf.low[1],
                CI_H = res_strat$conf.high[1],
                Interaction_P = p_val_str
            )
        }
    }
  }
}

cat("\n--- 4. Saving Target Trial Subgroup Effect Mod Results ---\n")
final_df <- bind_rows(stratified_hrs)
print(final_df)
write.csv(final_df, "ITT_Analysis/results/target_trial_subgroup_interactions.csv", row.names=FALSE)
cat("\nAnalysis complete! Results saved successfully.\n")
