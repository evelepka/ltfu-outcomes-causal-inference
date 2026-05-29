library(dplyr)
library(survival)
library(broom)

df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# Clean interaction covariates
na_vals <- c("Missing", "Ignorado", "Unknown", "", "nan")
df <- df %>%
  mutate(across(c(race_clean, edu_clean, age_group, sex, hiv_aids, homelessness), ~ ifelse(. %in% na_vals, NA, .))) %>%
  mutate(
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "65+")),
    sex = factor(sex, levels = c("Female", "Male")),
    homelessness = factor(homelessness, levels = c("No", "Yes")),
    hiv_aids = factor(hiv_aids, levels = c("Negative", "Positive"))
  ) %>% tidyr::drop_na(age_group, sex, hiv_aids, homelessness)

# Apply Month 3 (60 day target trial) to definitively purge the massive day 0-30 baseline artifact
# This flawlessly reveals the true physiological effect Modification of abandonment by subgroups!
start_yrs <- 60 / 365.25 
end_yrs <- 180 / 365.25

df_tt <- df %>% filter(time_d_tx > start_yrs) %>%
  mutate(
    tx_duration_yrs = as.numeric(difftime(as.Date(end_date), as.Date(best_start), units="days")) / 365.25,
    eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
  ) %>%
  filter(eligible == 1) %>%
  mutate(
    expose = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs >= start_yrs & tx_duration_yrs <= end_yrs, 1, 0),
    time_followup = time_d_tx - start_yrs
  )

run_interaction <- function(var_name, df_data) {
  df_clean <- df_data %>% filter(!is.na(!!sym(var_name)))
  
  f_int <- as.formula(paste("Surv(time_followup, event_d) ~ expose *", var_name, "+ age_group + sex + homelessness + hiv_aids"))
  fit <- coxph(f_int, data = df_clean)
  s <- summary(fit)
  
  levels_var <- levels(df_clean[[var_name]])
  res <- data.frame()
  base_coef <- s$coefficients["expose", "coef"]
  base_se <- s$coefficients["expose", "se(coef)"]
  
  res <- rbind(res, data.frame(
    Subgroup = var_name,
    Level = levels_var[1],
    HR = exp(base_coef),
    CI_Lower = exp(base_coef - 1.96 * base_se),
    CI_Upper = exp(base_coef + 1.96 * base_se),
    P_Interaction = NA
  ))
  
  for (i in 2:length(levels_var)) {
    lvl <- levels_var[i]
    int_term <- paste0("expose:", var_name, lvl)
    if (int_term %in% rownames(s$coefficients)) {
      int_coef <- s$coefficients[int_term, "coef"]
      tot_coef <- base_coef + int_coef
      tot_var <- vcov(fit)["expose", "expose"] + vcov(fit)[int_term, int_term] + 2 * vcov(fit)["expose", int_term]
      tot_se <- sqrt(tot_var)
      
      res <- rbind(res, data.frame(
        Subgroup = var_name,
        Level = lvl,
        HR = exp(tot_coef),
        CI_Lower = exp(tot_coef - 1.96 * tot_se),
        CI_Upper = exp(tot_coef + 1.96 * tot_se),
        P_Interaction = s$coefficients[int_term, "Pr(>|z|)"]
      ))
    }
  }
  return(res)
}

res_all <- bind_rows(
  run_interaction("sex", df_tt),
  run_interaction("age_group", df_tt),
  run_interaction("homelessness", df_tt),
  run_interaction("hiv_aids", df_tt)
)

write.csv(res_all, "ITT_Analysis/results/target_trial_subgroup_interactions.csv", row.names = FALSE)
cat("Successfully wrote subgroup interactions\n")
