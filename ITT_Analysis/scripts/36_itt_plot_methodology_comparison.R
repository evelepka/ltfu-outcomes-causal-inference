# 36. Methodology Contrast: Naive vs Target Trial Survival
# ==============================================================================
# Demonstrates the methodological paradox:
# Panel A shows the crude observational mortality from Day 0, where immortal 
# time bias makes abandonment appear artificially protective early on.
# Panel B shows the adjusted Target Trial causal emulation (Month 2), removing 
# immortal time bias and confounding to reveal the true massive hazard.
# ==============================================================================
library(dplyr)
library(survival)
library(ggplot2)
library(gridExtra)

cat("\n--- 1. Loading Cohort ---\n")
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

get_mode <- function(v) {
  t <- table(v)
  names(t)[which.max(t)]
}

cat("\n--- 2. Panel A: Naive Crude Observational Survival (Immortal Time Bias) ---\n")
# Cap follow-up at 2 years to match TT scale
df_naive <- df %>%
  mutate(
    event_d = ifelse(time_d > 2.0, 0, event_d),
    time_d = ifelse(time_d > 2.0, 2.0, time_d)
  )
fit_naive <- survfit(Surv(time_d, event_d) ~ itt_group, data = df_naive)

df_p_naive <- data.frame(
  Time = fit_naive$time,
  CumMort = 1 - fit_naive$surv,
  strata = rep(names(fit_naive$strata), fit_naive$strata),
  Panel = "A: Naive Observational (Crude from Tx Start)"
)
df_p_naive$strata <- gsub("itt_group=", "", df_p_naive$strata)
# Re-label to match exposing language
df_p_naive$Cohort <- ifelse(df_p_naive$strata == "Loss to follow-up", "Abandoned Treatment", "Maintained Care")

cat("\n--- 3. Panel B: Causal Target Trial (Month 2) ---\n")
# We use Trial 2 (Abandonment in Month 2) as a strong illustrative example
m <- 2
start_yrs <- ((m - 1) * 30) / 365.25
end_yrs   <- (m * 30) / 365.25

df_trial <- df_c %>% filter(time_d > start_yrs) %>%
  mutate(
    eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
  ) %>% filter(eligible == 1) %>%
  mutate(
    expose = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs >= start_yrs & tx_duration_yrs < end_yrs, 1, 0),
    time_followup = time_d - start_yrs,
    event_d = ifelse(time_followup > 2.0, 0, event_d),
    time_followup = ifelse(time_followup > 2.0, 2.0, time_followup)
  )

f_base <- as.formula(paste("Surv(time_followup, event_d) ~ expose +", paste(covariates, collapse = " + ")))
fit_tt <- coxph(f_base, data = df_trial)

ref_patient <- df_trial[1:2, c("expose", covariates)]
ref_patient$expose <- c(0, 1)
for (cov in covariates) {
  ref_patient[[cov]] <- factor(get_mode(df_trial[[cov]]), levels = levels(df_trial[[cov]]))
}

sfit <- survfit(fit_tt, newdata = ref_patient)
df_p0 <- data.frame(Time = sfit$time, CumMort = 1 - sfit$surv[, 1], Cohort = "Maintained Care", Panel = "B: Target Trial Emulation (Adjusted, Month 2)")
df_p1 <- data.frame(Time = sfit$time, CumMort = 1 - sfit$surv[, 2], Cohort = "Abandoned Treatment", Panel = "B: Target Trial Emulation (Adjusted, Month 2)")
df_start <- data.frame(Time = 0, CumMort = 0, Cohort = c("Maintained Care", "Abandoned Treatment"), Panel = "B: Target Trial Emulation (Adjusted, Month 2)")

df_p_tt <- bind_rows(df_start, df_p0, df_p1)

cat("\n--- 4. Assembling and Plotting ---\n")
# Start Naive at 0
df_start_naive <- data.frame(
    Time = 0, CumMort = 0,
    strata = "", Panel = "A: Naive Observational (Crude from Tx Start)",
    Cohort = c("Maintained Care", "Abandoned Treatment")
)
df_plot <- bind_rows(df_p_naive %>% select(-strata), df_start_naive %>% select(-strata), df_p_tt)

p_combined <- ggplot(df_plot, aes(x = Time, y = CumMort, color = Cohort, linetype = Cohort)) +
  geom_step(linewidth = 1.3) +
  facet_wrap(~ Panel, scales="free_x") +
  scale_color_manual(values = c("Maintained Care" = "#2c3e50", "Abandoned Treatment" = "#e74c3c")) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0, 0.15, 0.05), limits=c(0, 0.12)) +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  labs(
    title = "Resolving the Paradox of Treatment Abandonment",
    subtitle = "Panel A shows raw immortal time bias where early death artificially prevents abandonment. Panel B shows the true causal penalty.",
    x = "Years Since Initiation (Panel A)  /  Years Since Month 2 Trial Start (Panel B)",
    y = "Cumulative Mortality",
    color = "Exposure Cohort", linetype = "Exposure Cohort"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    strip.text = element_text(face = "bold", size = 13, hjust = 0),
    strip.background = element_rect(fill = "#ecf0f1", color = NA),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm")
  )

ts <- as.numeric(Sys.time())
fname <- sprintf("figure_methodology_contrast_%.0f.png", ts)
artifact_dir <- "/Users/jasonandrews/.gemini/antigravity/brain/c053ef30-5842-41b7-b342-bf735650d865"

ggsave(file.path("ITT_Analysis/results", fname), plot=p_combined, width=12, height=6.5, dpi=300, bg="white")
ggsave(file.path(artifact_dir, fname), plot=p_combined, width=12, height=6.5, dpi=300, bg="white")

cat("\n--- 5. Updating Walkthrough Cache ---\n")
wt_path <- file.path(artifact_dir, "walkthrough.md")
if(file.exists(wt_path)){
    writeLines(paste0(paste(readLines(wt_path, warn=F), collapse="\n"), 
               "\n![Methodology Contrast](", file.path(artifact_dir, fname), ")\n"), wt_path)
}
cat("Success! Generated", fname, "\n")
