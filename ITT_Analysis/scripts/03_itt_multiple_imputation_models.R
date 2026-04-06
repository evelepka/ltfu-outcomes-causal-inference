library(mice)
library(dplyr)
library(survival)
library(cmprsk)
library(broom)

# --- Configuration ---
base_path <- "/Users/evelynlepkadelima/Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper"
setwd(base_path)

# Load the newly generated ITT cohort
df_base <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# 1. Clean Variables for MI 
na_vals <- c("Missing", "Ignorado", "Unknown", "", "nan")
df <- df_base %>%
  filter(itt_group == "Loss to follow-up") %>% # Restrict to Loss to follow-up ONLY
  mutate(across(c(race_clean, edu_clean, dot_status, alcohol, drug_use, diabetes, hosp_admission, hiv_aids), ~ ifelse(. %in% na_vals, NA, .))) %>%
  mutate(
    # Ensure times are numeric and non-zero
    time_rn = as.numeric(ifelse(time_rn <= 0, 0.001, time_rn)),
    time_d = as.numeric(ifelse(time_d <= 0, 0.001, time_d)),

    # Combined status for Fine-Gray (Retreatment=1, Death=2, Censored=0)
    fg_status = case_when(
      event_rn == 1 ~ 1,
      event_d == 1 & event_rn == 0 ~ 2,
      TRUE ~ 0
    ),
    fg_time = as.numeric(ifelse(event_rn == 1, time_rn, time_d)),
    
    # Treatment month grouping is already in the CSV
    tx_month_grp = factor(tx_month_grp, levels = c("≥ 4 months", "< 2 months", "2 to <4 months")),
    
    # Clean age_group
    age_group = cut(age_tb, breaks=c(14, 24, 44, 64, 150), labels=c("15-24", "25-44", "45-64", "≥65")),
    
    # Use pre-cleaned factors explicitly
    sex = as.factor(sex),
    race_clean = as.factor(race_clean),
    edu_clean = as.factor(edu_clean),
    clinical_clean = as.factor(clinical_clean),
    
    hiv_aids = as.factor(hiv_aids),
    diabetes = as.factor(diabetes),
    alcohol = as.factor(alcohol),
    drug_use = as.factor(drug_use),
    incarcerated = as.factor(incarcerated),
    homelessness = as.factor(homelessness),
    hosp_admission = as.factor(hosp_admission),
    dot_status = as.factor(dot_status),
    
    # Ensure event types are factors for MI 
    fg_status = as.factor(fg_status),
    event_d = as.factor(event_d)
  )

# Set Reference Levels
df$age_group <- relevel(df$age_group, ref = "15-24")
df$race_clean <- relevel(df$race_clean, ref = "White")
df$edu_clean <- relevel(df$edu_clean, ref = "≥ 12 years")
df$clinical_clean <- relevel(df$clinical_clean, ref = "Pulmonary")
df$tx_month_grp <- relevel(df$tx_month_grp, ref = "≥ 4 months")

# 2. Setup MI Dataset
vars_mi <- c(
  "age_group", "sex", "race_clean", "edu_clean", 
  "hiv_aids", "diabetes", "alcohol", "drug_use", 
  "incarcerated", "homelessness", "hosp_admission", 
  "clinical_clean", "dot_status", "tx_month_grp",
  "fg_status", "fg_time", "event_d", "time_d"
)

df_mi <- df[, vars_mi]

# 3. Run Imputation
set.seed(42)
cat("Starting Imputation (m=5, maxit=5) for N =", nrow(df_mi), "...\n")
imp <- mice(df_mi, m = 5, maxit = 5, method = "pmm", printFlag = TRUE)
cat("Imputation finished. Fitting models...\n")

# 4. Multivariable Modeling 
cat("Fitting models...\n")

# Predictor list
all_vars <- c("age_group", "sex", "race_clean", "edu_clean", 
              "hiv_aids", "diabetes", "alcohol", "drug_use", 
              "incarcerated", "homelessness", "hosp_admission", 
              "clinical_clean", "dot_status", "tx_month_grp")

# --- Mortality (Cox) ---
# Adjusted
fit_cox_mv <- with(imp, coxph(Surv(time_d, as.numeric(as.character(event_d))) ~ age_group + sex + race_clean +
  edu_clean + hiv_aids + diabetes + alcohol + drug_use +
  incarcerated + homelessness + hosp_admission + clinical_clean + dot_status + tx_month_grp))
pool_cox_mv <- summary(pool(fit_cox_mv), exponentiate = TRUE, conf.int = TRUE)

# --- Retreatment (Fine-Gray) ---
fg_fits_mv <- list()
for (i in 1:5) {
  df_i <- complete(imp, i)
  fg_data <- finegray(Surv(fg_time, fg_status) ~ ., data = df_i, etype = "1")
  
  fit_fg_mv <- coxph(
    Surv(fgstart, fgstop, fgstatus) ~ age_group + sex + race_clean +
      edu_clean + hiv_aids + diabetes + alcohol + drug_use +
      incarcerated + homelessness + hosp_admission + clinical_clean + dot_status + tx_month_grp,
    data = fg_data, weights = fgwt
  )
  fg_fits_mv[[i]] <- fit_fg_mv
}
pool_fg_mv <- summary(pool(as.mira(fg_fits_mv)), exponentiate = TRUE, conf.int = TRUE)



# --- Univariate Modeling (for ALL variables) ---
cat("Running univariate models for all vars...\n")
uv_results <- data.frame()

for (v in all_vars) {
  # Mortality Univariate
  cox_uv_fits <- list()
  for (i in 1:5) {
    df_i <- complete(imp, i)
    form_cox <- as.formula(paste("Surv(time_d, as.numeric(as.character(event_d))) ~", v))
    cox_uv_fits[[i]] <- coxph(form_cox, data = df_i)
  }
  p_cox_uv <- summary(pool(as.mira(cox_uv_fits)), exponentiate = TRUE, conf.int = TRUE)
  p_cox_uv$model <- "Cox_Death_Unadjusted"
  
  # Retreatment Univariate
  fg_uv_fits <- list()
  for (i in 1:5) {
    df_i <- complete(imp, i)
    fg_data_i <- finegray(Surv(fg_time, fg_status) ~ ., data = df_i, etype = "1")
    form_fg <- as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", v))
    fit_fg_uv <- coxph(form_fg, data = fg_data_i, weights = fgwt)
    fg_uv_fits[[i]] <- fit_fg_uv
  }
  p_fg_uv <- summary(pool(as.mira(fg_uv_fits)), exponentiate = TRUE, conf.int = TRUE)
  p_fg_uv$model <- "FG_Retr_Unadjusted"
  
  uv_results <- bind_rows(uv_results, p_cox_uv, p_fg_uv)
}

# --- Complete Case Analysis (Sensitivity) ---
cat("Running Complete Case models...\n")
df_cc <- df %>%
  select(all_of(vars_mi)) %>%
  na.omit()

cat("Complete Case N =", nrow(df_cc), "\n")

# CC Mortality (Cox)
fit_cox_cc <- coxph(Surv(time_d, as.numeric(as.character(event_d))) ~ age_group + sex + race_clean +
  edu_clean + hiv_aids + diabetes + alcohol + drug_use +
  incarcerated + homelessness + hosp_admission + clinical_clean + dot_status + tx_month_grp,
  data = df_cc)
res_cox_cc <- summary(fit_cox_cc)

# CC Retreatment (Fine-Gray)
# Ensure fg_status is a factor in df_cc too
df_cc$fg_status <- factor(df_cc$fg_status)
fg_data_cc <- finegray(Surv(fg_time, fg_status) ~ ., data = df_cc, etype = "1")
fit_fg_cc <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age_group + sex + race_clean +
  edu_clean + hiv_aids + diabetes + alcohol + drug_use +
  incarcerated + homelessness + hosp_admission + clinical_clean + dot_status + tx_month_grp,
  data = fg_data_cc, weights = fgwt)
res_fg_cc <- summary(fit_fg_cc)

# 5. Format and Save Results
format_res <- function(summary_obj, model_name) {
  # Handle both pooled (mice) and standard (coxph) summaries
  if ("pool.summary" %in% class(summary_obj) || is.data.frame(summary_obj)) {
    res_df <- as.data.frame(summary_obj)
  } else {
    # It is a standard coxph summary object
    res_df <- as.data.frame(summary_obj$conf.int)
    res_df$p.value <- summary_obj$coefficients[, "Pr(>|z|)"]
    res_df$term <- rownames(res_df)
    colnames(res_df)[1] <- "estimate"
    colnames(res_df)[3] <- "conf.low"
    colnames(res_df)[4] <- "conf.high"
  }
  
  res_df %>%
    mutate(
      model = model_name,
      estimate_str = paste0(
        sprintf("%.2f", estimate), " (",
        sprintf("%.2f", conf.low), "-",
        sprintf("%.2f", conf.high), ")"
      ),
      p_value = case_when(
        p.value < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", p.value)
      )
    ) %>%
    select(model, term, estimate = estimate_str, p_value)
}

# Special formatter for standard summary objects (non-pooled)
format_std <- function(sum_obj, model_name) {
  df_res <- as.data.frame(sum_obj$conf.int)
  # Standard summary$conf.int columns: [exp(coef), exp(-coef), lower .95, upper .95]
  df_res$p <- sum_obj$coefficients[, "Pr(>|z|)"]
  df_res$term <- rownames(df_res)
  
  df_res %>%
    mutate(
      model = model_name,
      estimate_str = paste0(
        sprintf("%.2f", .data[[colnames(df_res)[1]]]), " (",
        sprintf("%.2f", .data[[colnames(df_res)[3]]]), "-",
        sprintf("%.2f", .data[[colnames(df_res)[4]]]), ")"
      ),
      p_val = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
    ) %>%
    select(model, term, estimate = estimate_str, p_value = p_val)
}

final_results <- bind_rows(
  format_res(pool_cox_mv, "Cox_Death_Adjusted_MI"),
  format_res(pool_fg_mv, "FG_Retr_Adjusted_MI"),
  format_std(res_cox_cc, "Cox_Death_Adjusted_CC"),
  format_std(res_fg_cc, "FG_Retr_Adjusted_CC"),
  format_res(uv_results %>% filter(model == "Cox_Death_Unadjusted"), "Cox_Death_Unadjusted"),
  format_res(uv_results %>% filter(model == "FG_Retr_Unadjusted"), "FG_Retr_Unadjusted")
)

write.csv(final_results, "ITT_Analysis/results/multivariable_results_mi_cc.csv", row.names = FALSE)
cat("All pooled MI and CC results saved to ITT_Analysis/results/\n")
