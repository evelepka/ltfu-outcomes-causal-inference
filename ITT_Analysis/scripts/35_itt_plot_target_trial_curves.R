# 35. Target Trial Adjusted Survival Curves
# ==============================================================================
# Computes the covariate-adjusted cumulative mortality (1 - survival) directly
# from the exact Target Trial Cox proportional hazard models.
# Plots continuous risk curves for Abandonment vs Maintained Care 
# across all 6 sequential trial months.
# ==============================================================================
library(dplyr)
library(survival)
library(ggplot2)
library(gridExtra)
library(stringr)

cat("\n--- 1. Loading Cohort & Covariates ---\n")
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

clean_data <- function(data) {
  data %>% 
  mutate(across(c(race_clean, edu_clean, dot_status, alcohol, drug_use, diabetes, hosp_admission, hiv_aids, clinical_clean, incarcerated, homelessness), ~ as.factor(.))) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "65+"))
  ) %>% tidyr::drop_na(age_group, sex, race_clean, edu_clean, hiv_aids, diabetes, alcohol, drug_use, incarcerated, homelessness, hosp_admission, clinical_clean, dot_status)
}

df_c <- clean_data(df) %>%
  mutate(
    date_start = as.Date(best_start),
    date_end = as.Date(end_date),
    tx_duration_yrs = as.numeric(date_end - date_start) / 365.25
  )

covariates <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids", "diabetes", 
                "alcohol", "drug_use", "incarcerated", "homelessness", "hosp_admission", 
                "clinical_clean", "dot_status")

# Helper function to find the mode
get_mode <- function(v) {
  t <- table(v)
  names(t)[which.max(t)]
}

cat("\n--- 2. Fitting Sequence of Target Trials ---\n")
plot_data <- list()

for (m in 1:6) {
  cat(paste("Processing Target Trial Month", m, "...\n"))
  start_days <- (m - 1) * 30
  end_days   <- m * 30
  start_yrs <- start_days / 365.25
  end_yrs   <- end_days / 365.25
  
  df_trial <- df_c %>% filter(time_d_tx > start_yrs) %>%
    mutate(
      eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
    ) %>% filter(eligible == 1) %>%
    mutate(
      expose = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs >= start_yrs & tx_duration_yrs < end_yrs, 1, 0),
      time_followup = time_d_tx - start_yrs,
      # 2-year administrative censoring
      event_d = ifelse(time_followup > 2.0, 0, event_d),
      time_followup = ifelse(time_followup > 2.0, 2.0, time_followup)
    )
  
  # Fit Model
  f_base <- as.formula(paste("Surv(time_followup, event_d) ~ expose +", paste(covariates, collapse = " + ")))
  fit <- coxph(f_base, data = df_trial)
  
  # Create a reference archetype patient (using the exact factor levels present in the model)
  ref_patient <- df_trial[1:2, c("expose", covariates)]
  ref_patient$expose <- c(0, 1)
  
  # Set to mode of the cohort to represent the "typical" patient baseline hazard
  for (cov in covariates) {
    mod_val <- get_mode(df_trial[[cov]])
    ref_patient[[cov]] <- factor(mod_val, levels = levels(df_trial[[cov]]))
  }
  
  # Predict Adjusted Survival
  sfit <- survfit(fit, newdata = ref_patient)
  
  # sfit$surv is a matrix where column 1 is expose=0, column 2 is expose=1
  df_p0 <- data.frame(
    Trial_Month = paste("Trial", m, "(Abandon in Month", m, ")"),
    Time = sfit$time,
    CumMort = 1 - sfit$surv[, 1],
    Cohort = "Maintained Care"
  )
  df_p1 <- data.frame(
    Trial_Month = paste("Trial", m, "(Abandon in Month", m, ")"),
    Time = sfit$time,
    CumMort = 1 - sfit$surv[, 2],
    Cohort = "Abandoned Treatment"
  )
  
  # Always start at 0
  df_start <- data.frame(
      Trial_Month = paste("Trial", m, "(Abandon in Month", m, ")"),
      Time = 0,
      CumMort = 0,
      Cohort = c("Maintained Care", "Abandoned Treatment")
  )
  
  plot_data[[m]] <- bind_rows(df_start, df_p0, df_p1)
}

df_final <- bind_rows(plot_data)

cat("\n--- 3. Plotting Subplots ---\n")
p1 <- ggplot(df_final, aes(x = Time, y = CumMort, color = Cohort, linetype = Cohort)) +
  geom_step(linewidth = 1.2) +
  facet_wrap(~ Trial_Month, ncol = 3) +
  scale_color_manual(values = c("Maintained Care" = "#2c3e50", "Abandoned Treatment" = "#e74c3c")) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0, 0.10, 0.02)) +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  labs(
    title = "Covariate-Adjusted Cumulative Mortality from Target Trial Emulations",
    subtitle = "Standardized 2-year survival curves adjusting for baseline clinical and demographic severity",
    x = "Years Since Target Trial Enrollment (Start of Month)",
    y = "Adjusted Cumulative Mortality",
    color = "Exposure Group", linetype = "Exposure Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm"),
    panel.grid.minor = element_blank()
  )

ts <- as.numeric(Sys.time())
fname <- sprintf("figure_target_trial_curves_%.0f.png", ts)
artifact_dir <- "/Users/jasonandrews/.gemini/antigravity/brain/c053ef30-5842-41b7-b342-bf735650d865"

ggsave(file.path("ITT_Analysis/results", fname), plot=p1, width=14, height=8, dpi=300, bg="white")
ggsave(file.path(artifact_dir, fname), plot=p1, width=14, height=8, dpi=300, bg="white")

cat("\n--- 4. Updating Walkthrough Cache ---\n")
wt_path <- file.path(artifact_dir, "walkthrough.md")
if(file.exists(wt_path)){
    writeLines(paste0(paste(readLines(wt_path, warn=F), collapse="\n"), 
               "\n![Alternative Figure: Target Trial Curves](", file.path(artifact_dir, fname), ")\n"), wt_path)
}
cat("Success! Generated", fname, "\n")
