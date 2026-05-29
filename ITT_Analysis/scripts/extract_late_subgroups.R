library(dplyr)
library(survival)

df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)
na_vals <- c("Missing", "Ignorado", "Unknown", "", "nan")
df <- df %>%
  mutate(across(c(race_clean, edu_clean, age_group, sex, hiv_aids, homelessness), ~ ifelse(. %in% na_vals, NA, .))) %>%
  mutate(
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "65+")),
    sex = factor(sex, levels = c("Female", "Male")),
    homelessness = factor(homelessness, levels = c("No", "Yes")),
    hiv_aids = factor(hiv_aids, levels = c("Negative", "Positive"))
  ) %>% tidyr::drop_na(age_group, sex, hiv_aids, homelessness)

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
  ) %>% 
  mutate(
    event_cap = ifelse(time_followup > 5, 0, event_d),
    time_cap = pmin(time_followup, 5)
  ) %>% filter(time_cap > 0)

df_split <- survSplit(Surv(time_cap, event_cap) ~ expose + age_group + sex + homelessness + hiv_aids, data = df_tt, cut = c(0.5), episode = "Epoch")

get_late <- function(var_name) {
  df_late <- df_split %>% filter(Epoch == 2)
  form <- as.formula(paste("Surv(tstart, time_cap, event_cap) ~ expose *", var_name, "+ age_group + sex + homelessness + hiv_aids"))
  fit <- coxph(form, data = df_late)
  s <- summary(fit)
  
  levels_var <- levels(df_late[[var_name]])
  base_coef <- s$coefficients["expose", "coef"]
  base_se <- s$coefficients["expose", "se(coef)"]
  
  res <- data.frame(
    Subgroup = var_name,
    Level = levels_var[1],
    HR = exp(base_coef),
    CI_Lower = exp(base_coef - 1.96 * base_se),
    CI_Upper = exp(base_coef + 1.96 * base_se)
  )
  
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
        CI_Upper = exp(tot_coef + 1.96 * tot_se)
      ))
    }
  }
  return(res)
}

df_sub <- bind_rows(
  get_late("sex"),
  get_late("age_group"),
  get_late("homelessness"),
  get_late("hiv_aids")
)
write.csv(df_sub, "ITT_Analysis/results/target_trial_subgroup_late.csv", row.names=FALSE)
