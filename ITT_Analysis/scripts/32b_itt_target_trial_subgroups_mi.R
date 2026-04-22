# 32b. Sequential Target Trial Subgroup HRs — Multiply-Imputed
# ==============================================================================
# MI-pooled version of 32_itt_target_trial_subgroups.R. Reads the five imputed
# datasets produced by 03a_itt_mi_miceforest.py (ITT_Analysis/data/mi/imp_*.csv)
# and, for each subgroup level, fits the 6-month sequential target trial Cox
# model on each imputation, then pools expose HRs via Rubin's rules.
# Output: ITT_Analysis/results/target_trial_subgroup_interactions_mi.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(mice)
  library(broom)
})

# -- Paths -------------------------------------------------------------------
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

# -- Load imputed datasets ---------------------------------------------------
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
if (length(imp_files) == 0) {
  stop(sprintf("No imputed datasets in %s — run 03a_itt_mi_miceforest.py first.",
               ITT_MI_DIR))
}
M <- length(imp_files)
cat(sprintf("[32b] Found %d imputed datasets\n", M))

prepare_imp <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date_start <- as.Date(d$best_start)
  d$date_end   <- as.Date(d$end_date)
  d$tx_duration_yrs <- as.numeric(d$date_end - d$date_start) / 365.25

  d$age_group <- factor(d$age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
  for (v in c("sex", "race_clean", "edu_clean", "hiv_aids", "diabetes",
              "alcohol", "drug_use", "incarcerated", "homelessness",
              "hosp_admission", "clinical_clean", "dot_status")) {
    d[[v]] <- as.factor(d[[v]])
  }
  d
}

imp_list <- lapply(imp_files, prepare_imp)
cat(sprintf("[32b] Each imputed cohort: %d rows\n", nrow(imp_list[[1]])))

COVARIATES <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids",
                "diabetes", "alcohol", "drug_use", "incarcerated",
                "homelessness", "hosp_admission", "clinical_clean",
                "dot_status", "trial_month")

# -- Build the 6-month pooled target-trial array for one imputed df ----------
build_pooled <- function(d) {
  trial_list <- vector("list", 6)
  for (m in 1:6) {
    start_yrs <- (m - 1) * 30 / 365.25
    end_yrs   <-  m      * 30 / 365.25
    tr <- d |>
      dplyr::filter(time_d > start_yrs) |>
      dplyr::mutate(eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)) |>
      dplyr::filter(eligible == 1) |>
      dplyr::mutate(
        expose = ifelse(itt_group == "Loss to follow-up" &
                          tx_duration_yrs >= start_yrs &
                          tx_duration_yrs <  end_yrs, 1, 0),
        trial_month = paste0("Month_", m),
        time_followup = time_d - start_yrs,
        event_d_num = as.numeric(as.character(event_d)),
        # 2-year horizon
        event_d_num = ifelse(time_followup > 2.0, 0, event_d_num),
        time_followup = ifelse(time_followup > 2.0, 2.0, time_followup)
      )
    trial_list[[m]] <- tr
  }
  pooled <- bind_rows(trial_list)
  pooled$trial_month <- factor(pooled$trial_month, levels = paste0("Month_", 1:6))
  pooled
}

# -- Subgroup-stratified MI pooling ------------------------------------------
SUBGROUPS <- c("age_group", "sex", "hiv_aids", "homelessness")
all_results <- list()

for (sg in SUBGROUPS) {
  cat(sprintf("\n[32b] Subgroup: %s\n", sg))
  base_covs <- setdiff(COVARIATES, sg)

  # Determine the levels present (use the first imputation)
  lvls <- sort(unique(as.character(imp_list[[1]][[sg]])))
  lvls <- lvls[!is.na(lvls) & lvls != ""]

  for (lvl in lvls) {
    cat(sprintf("  level: %s\n", lvl))
    fits <- vector("list", M)
    for (i in seq_len(M)) {
      pooled_i <- build_pooled(imp_list[[i]])
      sub_i <- pooled_i[pooled_i[[sg]] == lvl & !is.na(pooled_i[[sg]]), ]
      if (nrow(sub_i) < 100 || sum(sub_i$event_d_num, na.rm = TRUE) < 5) {
        fits[[i]] <- NULL
        next
      }
      rhs <- paste("expose +", paste(base_covs, collapse = " + "))
      f <- as.formula(paste("Surv(time_followup, event_d_num) ~", rhs))
      fit <- tryCatch(
        coxph(f, data = sub_i, cluster = sinan_clean),
        error = function(e) NULL
      )
      fits[[i]] <- fit
    }
    fits <- Filter(Negate(is.null), fits)
    if (length(fits) == 0) {
      cat("    (no fittable model in any imputation)\n")
      next
    }
    pooled_fit <- tryCatch(
      summary(pool(as.mira(fits)), exponentiate = TRUE, conf.int = TRUE),
      error = function(e) { cat("    pool error: ", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(pooled_fit)) next
    expose_row <- pooled_fit[pooled_fit$term == "expose", ]
    if (nrow(expose_row) == 0) next

    all_results[[length(all_results) + 1]] <- data.frame(
      Subgroup = sg,
      Level    = lvl,
      HR       = expose_row$estimate,
      CI_L     = expose_row$conf.low,
      CI_H     = expose_row$conf.high,
      P_Value  = expose_row$p.value,
      N_Imp    = length(fits)
    )
  }
}

final_df <- bind_rows(all_results)
out_path <- file.path(ITT_RESULTS_DIR, "target_trial_subgroup_interactions_mi.csv")
write.csv(final_df, out_path, row.names = FALSE)
cat(sprintf("\n[32b] Wrote %d rows to %s\n", nrow(final_df), out_path))
print(final_df)
