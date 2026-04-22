# 03_itt_multiple_imputation_models.R
# -----------------------------------
# Fits Cox (mortality) and Fine-Gray (retreatment with death as competing risk)
# models on the LTFU subgroup, pooled across m=5 multiply-imputed datasets,
# plus complete-case sensitivity analyses.
#
# Imputation is NOT done here. Run 03a_itt_mi_miceforest.py first to produce
# ITT_Analysis/data/mi/imp_01.csv ... imp_05.csv (miceforest; ~10-30x faster
# than the previous R mice-based pipeline). This script reads those files,
# fits models per imputation, and pools via Rubin's rules through mice::pool().

suppressPackageStartupMessages({
  library(mice)
  library(dplyr)
  library(survival)
  library(broom)
})

# -- Paths -------------------------------------------------------------------
# Locate this script's directory robustly (works under Rscript, source(), and
# interactive eval).
.here <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  frames <- sys.frames()
  for (f in rev(frames)) {
    of <- f$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  getwd()
}
source(file.path(.here(), "_paths.R"))

if (!dir.exists(ITT_MI_DIR)) {
  stop(sprintf(
    "No imputed datasets found at %s. Run 03a_itt_mi_miceforest.py first.",
    ITT_MI_DIR
  ))
}
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
if (length(imp_files) == 0) {
  stop(sprintf("No imp_*.csv files in %s", ITT_MI_DIR))
}
M <- length(imp_files)
cat(sprintf("[mi] Found %d imputed datasets in %s\n", M, ITT_MI_DIR))

# -- Load each imputed dataset and harmonize types --------------------------
# Reference levels match the previous R script.
prepare_imp <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)

  # Times: ensure numeric and non-zero
  d$time_rn <- as.numeric(ifelse(d$time_rn <= 0, 0.001, d$time_rn))
  d$time_d  <- as.numeric(ifelse(d$time_d  <= 0, 0.001, d$time_d))
  d$fg_time <- as.numeric(d$fg_time)

  d$tx_month_grp <- factor(d$tx_month_grp,
                           levels = c("≥ 4 months", "< 2 months", "2 to <4 months"))
  d$tx_month_grp <- relevel(d$tx_month_grp, ref = "≥ 4 months")

  d$age_group    <- factor(d$age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
  d$age_group    <- relevel(d$age_group, ref = "15-24")
  d$race_clean   <- relevel(factor(d$race_clean),   ref = "White")
  d$edu_clean    <- relevel(factor(d$edu_clean),    ref = "≥ 12 years")
  d$clinical_clean <- relevel(factor(d$clinical_clean), ref = "Pulmonary")

  for (v in c("sex", "hiv_aids", "diabetes", "alcohol", "drug_use",
              "incarcerated", "homelessness", "hosp_admission", "dot_status",
              "fg_status", "event_d")) {
    d[[v]] <- as.factor(d[[v]])
  }
  d
}

imp_list_full <- lapply(imp_files, prepare_imp)
cat(sprintf("[mi] Full cohort N = %d (each imputation)\n", nrow(imp_list_full[[1]])))

# This analysis is the within-LTFU risk-factor model. Filter each imputed
# dataset to LTFU for modeling.
imp_list <- lapply(imp_list_full, function(d) {
  d[d$itt_group == "Loss to follow-up", , drop = FALSE]
})
N_LTFU <- nrow(imp_list[[1]])
cat(sprintf("[mi] LTFU subgroup N = %d (each imputation)\n", N_LTFU))

# -- Covariate lists --------------------------------------------------------
all_vars <- c(
  "age_group", "sex", "race_clean", "edu_clean",
  "hiv_aids",  "diabetes", "alcohol",   "drug_use",
  "incarcerated", "homelessness", "hosp_admission",
  "clinical_clean", "dot_status", "tx_month_grp"
)
mv_formula_rhs <- paste(all_vars, collapse = " + ")

# -- Helper: fit a list of models across imputations, then pool -------------
fit_and_pool <- function(fit_fn) {
  fits <- lapply(imp_list, fit_fn)
  pooled <- summary(pool(as.mira(fits)), exponentiate = TRUE, conf.int = TRUE)
  as.data.frame(pooled)
}

# -- 1. Adjusted Cox (mortality) --------------------------------------------
cat("\n[fit] Adjusted Cox (mortality)...\n")
pool_cox_mv <- fit_and_pool(function(d) {
  coxph(as.formula(paste("Surv(time_d, as.numeric(as.character(event_d))) ~", mv_formula_rhs)),
        data = d)
})

# -- 2. Adjusted Fine-Gray (retreatment with death competing) ---------------
cat("[fit] Adjusted Fine-Gray (retreatment)...\n")
pool_fg_mv <- fit_and_pool(function(d) {
  fg_data <- finegray(Surv(fg_time, fg_status) ~ ., data = d, etype = "1")
  coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", mv_formula_rhs)),
        data = fg_data, weights = fgwt)
})

# -- 3. Univariate Cox + Fine-Gray for each variable ------------------------
cat("[fit] Univariate Cox + Fine-Gray for 14 covariates...\n")
uv_results <- list()
for (v in all_vars) {
  cox_uv <- fit_and_pool(function(d) {
    coxph(as.formula(paste("Surv(time_d, as.numeric(as.character(event_d))) ~", v)), data = d)
  })
  cox_uv$model <- "Cox_Death_Unadjusted"

  fg_uv <- fit_and_pool(function(d) {
    fg_data <- finegray(Surv(fg_time, fg_status) ~ ., data = d, etype = "1")
    coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", v)),
          data = fg_data, weights = fgwt)
  })
  fg_uv$model <- "FG_Retr_Unadjusted"

  uv_results[[v]] <- bind_rows(cox_uv, fg_uv)
}
uv_df <- bind_rows(uv_results)

# -- 4. Complete-case sensitivity -------------------------------------------
cat("[fit] Complete-case Cox + Fine-Gray...\n")
# Drop rows with NA in any MV covariate from the first imputed dataset we can
# just as well derive CC from any one of them since the non-missingness
# pattern is preserved. Use the raw itt_cohort.csv to be rigorous.
raw <- read.csv(COHORT_CSV, stringsAsFactors = FALSE) |>
  dplyr::filter(itt_group == "Loss to follow-up")

na_vals <- c("Missing", "Ignorado", "Unknown", "", "nan")
for (v in c("race_clean", "edu_clean", "dot_status", "alcohol", "drug_use",
            "diabetes", "hosp_admission", "hiv_aids")) {
  raw[[v]] <- ifelse(raw[[v]] %in% na_vals, NA, raw[[v]])
}
raw$time_rn <- as.numeric(ifelse(raw$time_rn <= 0, 0.001, raw$time_rn))
raw$time_d  <- as.numeric(ifelse(raw$time_d  <= 0, 0.001, raw$time_d))
raw$fg_status <- with(raw, ifelse(event_rn == 1, 1,
                           ifelse(event_d == 1 & event_rn == 0, 2, 0)))
raw$fg_time   <- with(raw, ifelse(event_rn == 1, time_rn, time_d))

raw$age_group <- cut(raw$age_tb, breaks = c(14, 24, 44, 64, 150),
                     labels = c("15-24", "25-44", "45-64", "≥65"))
raw$age_group <- relevel(factor(raw$age_group), ref = "15-24")
raw$tx_month_grp <- factor(raw$tx_month_grp,
                           levels = c("≥ 4 months", "< 2 months", "2 to <4 months"))
raw$tx_month_grp <- relevel(raw$tx_month_grp, ref = "≥ 4 months")
raw$race_clean   <- relevel(factor(raw$race_clean),   ref = "White")
raw$edu_clean    <- relevel(factor(raw$edu_clean),    ref = "≥ 12 years")
raw$clinical_clean <- relevel(factor(raw$clinical_clean), ref = "Pulmonary")
for (v in c("sex", "hiv_aids", "diabetes", "alcohol", "drug_use",
            "incarcerated", "homelessness", "hosp_admission", "dot_status",
            "fg_status", "event_d")) {
  raw[[v]] <- as.factor(raw[[v]])
}

df_cc <- raw |>
  dplyr::select(all_of(c(all_vars, "fg_status", "fg_time", "event_d", "time_d"))) |>
  tidyr::drop_na()
cat(sprintf("[cc] Complete-case N = %d\n", nrow(df_cc)))

fit_cox_cc <- coxph(as.formula(paste(
  "Surv(time_d, as.numeric(as.character(event_d))) ~", mv_formula_rhs)),
  data = df_cc)
fg_cc_data <- finegray(Surv(fg_time, fg_status) ~ ., data = df_cc, etype = "1")
fit_fg_cc  <- coxph(as.formula(paste(
  "Surv(fgstart, fgstop, fgstatus) ~", mv_formula_rhs)),
  data = fg_cc_data, weights = fgwt)

# -- 5. Format results and save --------------------------------------------
format_pooled <- function(d, model_name) {
  d <- as.data.frame(d)
  d$model <- model_name
  d$estimate_str <- sprintf("%.2f (%.2f-%.2f)", d$estimate, d$conf.low, d$conf.high)
  d$p_value <- ifelse(d$p.value < 0.001, "<0.001", sprintf("%.3f", d$p.value))
  d[, c("model", "term", "estimate_str", "p_value")] |>
    dplyr::rename(estimate = estimate_str)
}

format_cc <- function(fit, model_name) {
  s <- summary(fit)
  d <- as.data.frame(s$conf.int)
  colnames(d) <- c("estimate", "neg_estimate", "conf.low", "conf.high")
  d$p.value <- s$coefficients[, "Pr(>|z|)"]
  d$term <- rownames(d)
  d$model <- model_name
  d$estimate_str <- sprintf("%.2f (%.2f-%.2f)", d$estimate, d$conf.low, d$conf.high)
  d$p_value <- ifelse(d$p.value < 0.001, "<0.001", sprintf("%.3f", d$p.value))
  d[, c("model", "term", "estimate_str", "p_value")] |>
    dplyr::rename(estimate = estimate_str)
}

final <- bind_rows(
  format_pooled(pool_cox_mv, "Cox_Death_Adjusted_MI"),
  format_pooled(pool_fg_mv,  "FG_Retr_Adjusted_MI"),
  format_cc(fit_cox_cc,      "Cox_Death_Adjusted_CC"),
  format_cc(fit_fg_cc,       "FG_Retr_Adjusted_CC"),
  format_pooled(uv_df |> dplyr::filter(model == "Cox_Death_Unadjusted"),
                "Cox_Death_Unadjusted"),
  format_pooled(uv_df |> dplyr::filter(model == "FG_Retr_Unadjusted"),
                "FG_Retr_Unadjusted")
)

out_csv <- file.path(ITT_RESULTS_DIR, "multivariable_results_mi_cc.csv")
dir.create(ITT_RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
write.csv(final, out_csv, row.names = FALSE)
cat(sprintf("\n[done] Pooled MI + complete-case results -> %s\n", out_csv))
