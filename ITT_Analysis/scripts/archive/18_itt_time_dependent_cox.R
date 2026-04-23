# ============================================================================
# TIME-DEPENDENT COX PROPORTIONAL HAZARDS MODEL (PRIMARY CAUSAL ESTIMATOR)
# ============================================================================
# This script runs the definitive causal model to isolate the effect of LTFU 
# on mortality. It completely eliminates Immortal Time Bias without incurring
# Survivor Bias (unlike the Landmark approach), by modeling LTFU as a 
# time-varying covariate (tstart, tstop, event).
# ============================================================================

library(dplyr)
library(survival)
library(broom)

# --- 1. Load Data ---
data_file <- "Abandonment Paper/ITT_Analysis/data/itt_cohort.csv"
if (!file.exists(data_file)) stop("Data not found: ", data_file)
df <- read.csv(data_file)

str(df) # basic check

# --- 2. Cohort Cleaning & Variable Harmonization ---
cat("Cleaning baseline variables...\n")
df <- df %>%
  mutate(
    # Time variables (convert to numeric years for precision if not already)
    time_followup = time_d_tx, 
    death_event = event_d,
    
    # Calculate exact time to LTFU for those who abandoned
    time_to_ltfu = if_else(itt_group == "Loss to follow-up", 
                           as.numeric(difftime(as.Date(end_date), as.Date(best_start), units="days")) / 365.25, 
                           NA_real_)
  )

# Fix edge case: someone abandons exactly on the day they die (or are censored)
# tmerge needs the exposure (LTFU) to occur strictly before the censoring/death event
df <- df %>%
  mutate(
    time_to_ltfu = if_else(!is.na(time_to_ltfu) & time_to_ltfu >= time_followup, 
                           time_followup - 0.001, 
                           time_to_ltfu),
    time_to_ltfu = if_else(!is.na(time_to_ltfu) & time_to_ltfu <= 0, 0.001, time_to_ltfu),
    time_followup = if_else(time_followup <= 0, 0.002, time_followup) # Ensure positive followup
  )

# Harmonize 13 Baseline Confounders as factors with strict reference levels
df$age_group <- factor(df$age_group, levels = c("25 - 44", "15 - 24", "45 - 64", "65+"))
df$sex <- factor(df$sex, levels = c("Male", "Female"))
df$race_clean <- factor(df$race_clean, levels = c("White", "Black or Mixed", "Other"))
df$edu_clean <- factor(df$edu_clean, levels = c("≥ 12 years", "None", "≤ 7 years", "8 - 11 years"))
df$hiv_aids <- factor(df$hiv_aids, levels = c("Negative", "Positive"))
df$dot_status <- factor(df$dot_status, levels = c("Yes", "No"))
df$alcohol <- factor(df$alcohol, levels = c("No", "Yes"))
df$drug_use <- factor(df$drug_use, levels = c("No", "Yes"))
df$incarcerated <- factor(df$incarcerated, levels = c("No", "Yes"))
df$homelessness <- factor(df$homelessness, levels = c("No", "Yes"))
df$hosp_admission <- factor(df$hosp_admission, levels = c("No", "Yes"))
df$clinical_clean <- factor(df$clinical_clean, levels = c("Pulmonary", "Extrapulmonary", "Pulmonary and Extrapulmonary or disseminated"))
df$diabetes <- factor(df$diabetes, levels = c("No", "Yes"))

# Drop missing values in the confounders to ensure clean model matrix
df_clean <- df %>% tidyr::drop_na(age_group, sex, race_clean, edu_clean, hiv_aids, 
                                  dot_status, alcohol, drug_use, incarcerated, 
                                  homelessness, hosp_admission, clinical_clean, diabetes)
cat("Cohort size after complete-case covariate matching for Cox:", nrow(df_clean), "\n")

# --- 3. Splitting Time (Time-Dependent Data Structure) ---
cat("Applying tmerge to split time at exact day of abandonment...\n")
# Everyone starts with ltfu_status = 0 (Non-LTFU)
# At time_to_ltfu, the status flips to 1 for those who abandon.
df_tv <- tmerge(
  data1 = df_clean, 
  data2 = df_clean, 
  id = sinan_clean, 
  death = event(time_followup, death_event),
  ltfu_status = tdc(time_to_ltfu)
)

# Convert ltfu_status to factor for modeling readability
df_tv$ltfu_status_factor <- factor(df_tv$ltfu_status, levels=c(0, 1), labels=c("In Treatment/Cured", "Abandoned (LTFU)"))

cat("Total records in time-varying dataset:", nrow(df_tv), "\n")
cat("Checking time-varying structure:\n")
print(head(df_tv %>% select(sinan_clean, time_to_ltfu, tstart, tstop, death, ltfu_status_factor), 10))


# --- 4. Main Time-Dependent Cox Model ---
cat("\nRunning Time-Dependent Cox Proportional Hazards Model...\n")
cox_td <- coxph(
  Surv(tstart, tstop, death) ~ ltfu_status_factor + age_group + sex + race_clean + 
                               edu_clean + hiv_aids + dot_status + alcohol + drug_use + 
                               incarcerated + homelessness + hosp_admission + 
                               clinical_clean + diabetes, 
  data = df_tv
)

summary_td <- summary(cox_td)
print(summary_td)

# --- 5. Formatting Output to CSV & Word ---
# Extract exactly what we want: the HR and 95% CI for LTFU status
td_res <- broom::tidy(cox_td, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "ltfu_status_factorAbandoned (LTFU)") %>%
  mutate(
    HR_formatted = sprintf("%.2f (%.2f - %.2f)", estimate, conf.low, conf.high),
    P_value = sprintf("%.4f", p.value)
  ) %>%
  select(term, estimate, conf.low, conf.high, HR_formatted, P_value)

cat("\n==============================================\n")
cat("FINAL CAUSAL ESTIMAND: EFFECT OF LTFU ON DEATH\n")
cat("HR:", td_res$HR_formatted, " | P-value:", td_res$P_value, "\n")
cat("==============================================\n")

outfile_csv <- "Abandonment Paper/ITT_Analysis/results/time_dependent_cox_hr.csv"
write.csv(td_res, outfile_csv, row.names=FALSE)

# Generate a minimal Word Document to report the "Gold Standard" causal finding
library(officer)
library(magrittr)

doc <- read_docx() %>%
  body_add_par("Primary Causal Estimator: Time-Dependent Cox Proportional Hazards", style = "heading 1") %>%
  body_add_par("This analysis treats LTFU (Treatment Abandonment) as a time-varying exposure. Every patient enters the cohort at Day 0 (Time of treatment initiation) as Unexposed. For patients who abandon, their person-time strictly before the abandonment date is credited to the Unexposed state, and their person-time after abandonment is credited to the Exposed state.", style = "Normal") %>%
  body_add_par("Advantage over Landmark: Eliminates immortal time bias robustly without conditioning on initial 180-day survival, thus preserving the entire cohort and averting survivor bias.", style = "Normal") %>%
  body_add_par(" ", style = "Normal") %>%
  body_add_par(paste0("Adjusted Hazard Ratio for LTFU: ", td_res$HR_formatted), style = "heading 2") %>%
  body_add_par(paste0("P-value: ", td_res$P_value), style = "Normal") %>%
  body_add_par(" ", style = "Normal") %>%
  body_add_par("Adjusted for 13 baseline demographic, social, behavioral, and clinical confounders exactly analogous to Table 2.", style = "Normal")

out_doc <- "Abandonment Paper/ITT_Analysis/results/Time_Dependent_Cox_Causal_Report.docx"
print(doc, target = out_doc)
cat("Word execution summary saved to:", out_doc, "\n")
