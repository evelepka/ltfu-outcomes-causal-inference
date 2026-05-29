# 13. ITT LTFU Descriptive Outcomes
# ==============================================================================
# Purpose: Recalculate survival and competing risk statistics specifically for 
#          the Loss to Follow-up group based on the strict ITT cohort.
# ==============================================================================

library(dplyr)
library(survival)

cat("Loading ITT cohort...\n")
df <- read.csv("Abandonment Paper/ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# Filter cohort to ONLY the Loss to Follow-up group, from `end_date` onwards (Landmark T=0)
df_ltfu <- df %>% filter(itt_group == "Loss to follow-up" & time_d > 0)
df_ltfu$event_rn_factor <- factor(df_ltfu$event_rn, levels=c(0, 1, 2), labels=c("Censored", "Retreatment", "Death"))

cat("\n=======================================================\n")
cat("Outcomes following loss to follow-up (Strict ITT Update)\n")
cat("=======================================================\n")
cat("Total LTFU Cohort Size:", nrow(df_ltfu), "\n")

# 1. Total Person-Years (PY) and median follow-up (Mortality)
total_py <- sum(df_ltfu$time_d)
median_fu <- median(df_ltfu$time_d)
iqr_fu_lower <- quantile(df_ltfu$time_d, 0.25)
iqr_fu_upper <- quantile(df_ltfu$time_d, 0.75)

cat(sprintf("\n--- FOLLOW-UP ---\n"))
cat(sprintf("Total Person-Years (Mortality): %.1f PY\n", total_py))
cat(sprintf("Median follow-up: %.2f years (IQR %.2f-%.2f)\n", median_fu, iqr_fu_lower, iqr_fu_upper))

# 2. Retreatment Statistics
num_retreated <- sum(df_ltfu$event_rn == 1)
pct_retreated <- (num_retreated / nrow(df_ltfu)) * 100

total_py_rt <- sum(df_ltfu$time_rn)
incidence_rt <- (num_retreated / total_py_rt) * 100000

median_fu_rt <- median(df_ltfu$time_rn)
iqr_fu_rt_l <- quantile(df_ltfu$time_rn, 0.25)
iqr_fu_rt_u <- quantile(df_ltfu$time_rn, 0.75)

df_retreated <- df_ltfu %>% filter(event_rn == 1)
median_time_rt <- median(df_retreated$time_rn)
iqr_time_rt_l <- quantile(df_retreated$time_rn, 0.25)
iqr_time_rt_u <- quantile(df_retreated$time_rn, 0.75)

cat(sprintf("\n--- RETREATMENT ---\n"))
cat(sprintf("Reinitiated treatment: %d (%.1f%%)\n", num_retreated, pct_retreated))
cat(sprintf("Retreatment Incidence Rate: %.1f per 100,000 PY\n", incidence_rt))
cat(sprintf("Total Person-time at risk for RT: %.1f PY\n", total_py_rt))
cat(sprintf("Median follow-up for RT: %.2f years (IQR %.2f-%.2f)\n", median_fu_rt, iqr_fu_rt_l, iqr_fu_rt_u))
cat(sprintf("Median time to retreatment (among retreated): %.2f years (IQR %.2f-%.2f)\n", median_time_rt, iqr_time_rt_l, iqr_time_rt_u))

# 3. Mortality Statistics
num_deaths <- sum(df_ltfu$event_d == 1)
pct_deaths <- (num_deaths / nrow(df_ltfu)) * 100
incidence_death <- (num_deaths / total_py) * 100000

df_died <- df_ltfu %>% filter(event_d == 1)
median_time_death <- median(df_died$time_d)
iqr_time_death_l <- quantile(df_died$time_d, 0.25)
iqr_time_death_u <- quantile(df_died$time_d, 0.75)

cat(sprintf("\n--- MORTALITY ---\n"))
cat(sprintf("All-cause deaths: %d (%.1f%%)\n", num_deaths, pct_deaths))
cat(sprintf("Mortality Incidence Rate: %.1f per 100,000 PY\n", incidence_death))
cat(sprintf("Median time to death (among dead): %.2f years (IQR %.2f-%.2f)\n", median_time_death, iqr_time_death_l, iqr_time_death_u))

# 4. Standard Cumulative Incidence (Kaplan-Meier & Aalen-Johansen)
km_mort <- survfit(Surv(time_d, event_d) ~ 1, data = df_ltfu)
km_rt <- survfit(Surv(time_rn, event_rn_factor) ~ 1, data = df_ltfu)

q_mort <- summary(km_mort, times=c(1, 11.5))
q_rt <- summary(km_rt, times=c(1, 2, 10, 11.5)) 

cat(sprintf("\n--- CUMULATIVE INCIDENCE ---\n"))
cat(sprintf("Mortality at 1 year: %.1f%%\n", (1-q_mort$surv[1])*100))
if(length(q_mort$surv) >= 2) cat(sprintf("Mortality at 11.5 years: %.1f%%\n", (1-q_mort$surv[2])*100))

# For Competing Risks, summary() returns a matrix in pstate. Column for "Retreatment" is 2 if levels are Censored(1), Retreatment(2), Death(3)
# Actually in pstate column names are the factor labels
rt_col <- which(colnames(q_rt$pstate) == "Retreatment")
cat(sprintf("Retreatment at 1 year: %.1f%%\n", q_rt$pstate[1, rt_col]*100))
cat(sprintf("Retreatment at 2 years: %.1f%%\n", q_rt$pstate[2, rt_col]*100))
cat(sprintf("Retreatment at 10 years: %.1f%%\n", q_rt$pstate[3, rt_col]*100))
cat(sprintf("Retreatment at 11.5 years: %.1f%%\n", q_rt$pstate[4, rt_col]*100))


# 5. Stratifications
stratify <- function(var_name, is_mortality=TRUE) {
    if(is_mortality) {
        f <- as.formula(paste("Surv(time_d, event_d) ~", var_name))
        cat(sprintf("\n--- %s Stratification (11.5y Mortality) ---\n", var_name))
        fit <- survfit(f, data = df_ltfu)
        res <- summary(fit, times=11.5)
        
        if(length(res$strata) > 0) {
            for(i in 1:length(res$strata)) {
                cat(sprintf("%s: %.1f%%\n", as.character(res$strata[i]), (1-res$surv[i])*100))
            }
        }
    } else {
        f <- as.formula(paste("Surv(time_rn, event_rn_factor) ~", var_name))
        cat(sprintf("\n--- %s Stratification (10y Retreatment) ---\n", var_name)) # Using 10y for RT context
        fit <- survfit(f, data = df_ltfu)
        res <- summary(fit, times=10)
        
        if(length(res$strata) > 0) {
             rt_idx <- which(colnames(res$pstate) == "Retreatment")
             for(i in 1:length(res$strata)) {
                cat(sprintf("%s: %.1f%%\n", as.character(res$strata[i]), res$pstate[i, rt_idx]*100))
             }
        }
    }
}

# HIV status
stratify("hiv_aids", is_mortality=TRUE)
stratify("hiv_aids", is_mortality=FALSE)

# Homelessness
stratify("homelessness", is_mortality=TRUE)
stratify("homelessness", is_mortality=FALSE)

# Supervised therapy
stratify("supervised_therapy", is_mortality=TRUE)
stratify("supervised_therapy", is_mortality=FALSE)

# Timing of Treatment Abandonment (Months from tx_start to end_date)
df_ltfu$t_start_date <- as.Date(df_ltfu$tx_start)
df_ltfu$t_end_date <- as.Date(df_ltfu$end_date)
df_ltfu$abandon_timing_months <- as.numeric(df_ltfu$t_end_date - df_ltfu$t_start_date) / 30.44

df_ltfu$abandon_timing_cat <- factor(case_when(
    df_ltfu$abandon_timing_months < 2 ~ "1 to <2 months",
    df_ltfu$abandon_timing_months >= 4 ~ ">=4 months",
    TRUE ~ "2-3 months"
), levels = c("1 to <2 months", "2-3 months", ">=4 months"))

stratify("abandon_timing_cat", is_mortality=TRUE)
stratify("abandon_timing_cat", is_mortality=FALSE)

cat("\nDone!\n")
