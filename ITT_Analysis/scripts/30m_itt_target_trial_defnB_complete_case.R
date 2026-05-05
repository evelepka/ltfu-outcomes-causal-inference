# 30m. Complete-case sensitivity: target trial under defn-B + grace
# ==============================================================================
# Mirrors 30h (all-cause) and 30i (cause-specific) but operates on individuals
# with COMPLETE covariate data (no imputation). Output CSVs are parallel to the
# MI-pooled primary outputs and consumable by the same downstream summarisers.
#
# Inputs:
#   COHORT_CSV                            (single non-imputed cohort)
#   DATA_DIR/Final_table_cleaned.csv      (for cause-of-death attribution lookup,
#                                          same logic as 30i)
#
# Outputs:
#   target_trial_defnB_cc_early_late_array.csv
#       Trial_Month × model (early|late) × cap (0.5, 2) × N_cc
#   target_trial_defnB_cc_cause_specific.csv
#       Trial_Month × cause × cap (2) × N_cc
#
# Notes:
#   - "Complete case" = drops any individual missing any of the 13 target-trial
#     covariates. Time-varying exposure construction is identical to 30h/30i.
#   - No Rubin's pooling: single fit per (Trial_Month, model, cap[, cause]).
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

.here <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  frames <- sys.frames()
  for (f in rev(frames)) { of <- f$ofile; if (!is.null(of)) return(dirname(normalizePath(of))) }
  getwd()
}
source(file.path(.here(), "_paths.R"))

GRACE_DAYS <- 30
GRACE_YRS  <- GRACE_DAYS / 365.25
SHIFT_YRS  <- 30 / 365.25

COVARS <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids",
            "diabetes", "alcohol", "drug_use", "incarcerated",
            "homelessness", "hosp_admission", "clinical_clean", "dot_status")
RHS <- paste("expose +", paste(COVARS, collapse = " + "))

# ---------------------------------------------------------------------------
# Cause-of-death attribution lookup (mirrors 30i)
# ---------------------------------------------------------------------------
cat("[30m] Building cause attribution lookup ...\n")
raw_path <- file.path(DATA_DIR, "Final_table_cleaned.csv")
raw <- read.csv(raw_path, stringsAsFactors = FALSE)
parse_dt <- function(x) as.Date(x, format = "%B %d, %Y")
raw$end_date <- parse_dt(raw$end_date)
raw$dod      <- parse_dt(raw$dod)

TRANSFER <- c("Transf Outro Municipio", "Transf Outro Estado/Pais")
novo <- raw[trimws(tools::toTitleCase(tolower(raw$case_type))) == "Novo", ]
novo <- novo[!is.na(novo$case_outcome)
             & nzchar(trimws(novo$case_outcome))
             & novo$case_outcome != "Mud Diag"
             & !novo$case_outcome %in% TRANSFER, ]
novo <- novo[order(novo$end_date), ]
first <- novo[!duplicated(novo$sinan_clean), c("sinan_clean", "case_outcome")]

deathrec <- raw[!is.na(raw$dod) & !is.na(raw$cause_of_death_code) &
                  nzchar(raw$cause_of_death_code), ]
deathrec <- deathrec[order(deathrec$dod), ]
deathrec <- deathrec[!duplicated(deathrec$sinan_clean, fromLast = TRUE),
                     c("sinan_clean", "cause_of_death_code")]
deathrec$cause_of_death_code <- toupper(trimws(deathrec$cause_of_death_code))

attr_lookup <- merge(first, deathrec, by = "sinan_clean", all = TRUE)

classify_cod <- function(cod, case_outcome) {
  cod_known <- !is.na(cod) & nzchar(cod)
  tb_strict <- cod_known & (
    grepl("^A1[5-9]", cod) | grepl("^B90", cod) | grepl("^B200", cod)
  )
  resp <- cod_known & grepl("^J[0-9]", cod)
  hiv_other <- cod_known & grepl("^B2[0-4]", cod) & !grepl("^B200", cod)
  cls <- rep("unknown", length(cod))
  cls[cod_known & !tb_strict & !resp & !hiv_other] <- "non_tb"
  cls[hiv_other] <- "hiv_other"
  cls[resp]      <- "respiratory"
  cls[tb_strict] <- "tb_strict"
  miss <- cls == "unknown"
  cls[miss & !is.na(case_outcome) & case_outcome == "Obito TB"]  <- "tb_via_tbweb"
  cls[miss & !is.na(case_outcome) & case_outcome == "Obito NTB"] <- "ntb_via_tbweb"
  cls
}
attr_lookup$cod_class <- classify_cod(attr_lookup$cause_of_death_code,
                                      attr_lookup$case_outcome)
attr_lookup$tb_hybrid     <- attr_lookup$cod_class %in% c("tb_strict", "tb_via_tbweb")
attr_lookup$nontb_hybrid  <- attr_lookup$cod_class %in% c("non_tb", "ntb_via_tbweb")
attr_lookup$tb_simonly    <- attr_lookup$cod_class == "tb_strict"
attr_lookup$nontb_simonly <- attr_lookup$cod_class == "non_tb"

# ---------------------------------------------------------------------------
# Cohort (non-imputed) — apply defn-B shift; restrict to complete cases
# ---------------------------------------------------------------------------
cat("[30m] Loading cohort and computing complete-case subset ...\n")
d <- read.csv(COHORT_CSV, stringsAsFactors = FALSE)
N_total <- nrow(d)

# Treat empty strings and "Unknown" / "Ignorado" as missing for CC filter
recode_missing <- function(x) {
  if (!is.character(x) && !is.factor(x)) return(x)
  v <- as.character(x)
  v[v %in% c("", "Unknown", "Ignorado", "NA", "Sem Informacao")] <- NA
  v
}
for (v in COVARS) d[[v]] <- recode_missing(d[[v]])

cc_mask <- complete.cases(d[, COVARS])
N_cc <- sum(cc_mask)
cat(sprintf("  N total cohort       : %d\n", N_total))
cat(sprintf("  N complete-case      : %d (%.1f%%)\n", N_cc, 100 * N_cc / N_total))
cat(sprintf("  N dropped (missing)  : %d (%.1f%%)\n",
            N_total - N_cc, 100 * (N_total - N_cc) / N_total))

# Per-covariate missingness summary
miss_summary <- sapply(COVARS, function(v) sum(is.na(d[[v]])))
cat("\n  Per-covariate missing counts:\n")
for (v in COVARS) cat(sprintf("    %-18s %6d (%.1f%%)\n",
                              v, miss_summary[v], 100 * miss_summary[v] / N_total))

d <- d[cc_mask, , drop = FALSE]
d$age_group <- factor(d$age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
for (v in c("sex", "race_clean", "edu_clean", "hiv_aids", "diabetes",
            "alcohol", "drug_use", "incarcerated", "homelessness",
            "hosp_admission", "clinical_clean", "dot_status")) {
  d[[v]] <- as.factor(d[[v]])
}

d$date_start <- as.Date(d$best_start)
d$date_end   <- as.Date(d$end_date)
d$tx_duration_yrs      <- as.numeric(d$date_end - d$date_start) / 365.25
d$true_tx_duration_yrs <- pmax(d$tx_duration_yrs - SHIFT_YRS, 1/365.25)

# Merge cause attribution
d <- merge(d, attr_lookup[, c("sinan_clean", "cod_class",
                              "tb_hybrid", "nontb_hybrid",
                              "tb_simonly", "nontb_simonly")],
           by = "sinan_clean", all.x = TRUE)

# ---------------------------------------------------------------------------
# Trial constructors (mirrors 30h / 30i)
# ---------------------------------------------------------------------------
build_trial <- function(d, m) {
  start_yrs   <- (m - 1) * 30 / 365.25
  end_yrs     <-  m      * 30 / 365.25
  origin_yrs  <- start_yrs + GRACE_YRS
  d |>
    dplyr::filter(time_d_tx > origin_yrs) |>
    dplyr::mutate(
      eligible = ifelse(itt_group == "Non-LTFU" |
                          true_tx_duration_yrs >= start_yrs, 1, 0)
    ) |>
    dplyr::filter(eligible == 1) |>
    dplyr::mutate(
      expose = ifelse(itt_group == "Loss to follow-up" &
                        true_tx_duration_yrs >= start_yrs &
                        true_tx_duration_yrs <  end_yrs, 1, 0),
      event_d_num = as.numeric(as.character(event_d)),
      time_raw    = time_d_tx - origin_yrs
    )
}

prep_window <- function(tr, model, cap) {
  e <- tr$event_d_num
  t <- tr$time_raw
  if (model == "early") {
    time_out  <- pmin(t, cap)
    event_out <- ifelse(t <= cap & e == 1, 1, 0)
  } else if (model == "late") {
    time_out  <- pmin(t, cap)
    event_out <- ifelse(t > 0.5 & t <= cap & e == 1, 1, 0)
  } else stop("model must be early or late")
  tr$time_out  <- time_out
  tr$event_out <- event_out
  tr[tr$time_out > 0, , drop = FALSE]
}

prep_cause <- function(tr, cause, cap = 2) {
  e <- tr$event_d_num
  t <- tr$time_raw
  if (cause == "all") {
    cause_event <- e
  } else {
    ind <- tr[[cause]]
    ind[is.na(ind)] <- FALSE
    cause_event <- ifelse(e == 1 & ind, 1, 0)
  }
  time_out  <- pmin(t, cap)
  event_out <- ifelse(t > 0.5 & t <= cap & cause_event == 1, 1, 0)
  tr$time_out  <- time_out
  tr$event_out <- event_out
  tr[tr$time_out > 0, , drop = FALSE]
}

fit_one <- function(sub_df, .label = "") {
  if (nrow(sub_df) < 50) {
    cat(sprintf("    [%s] skip: nrow=%d < 50\n", .label, nrow(sub_df))); return(NULL) }
  if (sum(sub_df$event_out) < 5) {
    cat(sprintf("    [%s] skip: events=%d < 5\n", .label, sum(sub_df$event_out))); return(NULL) }
  if (length(unique(sub_df$expose)) < 2) {
    cat(sprintf("    [%s] skip: only one expose level\n", .label)); return(NULL) }
  f <- tryCatch(
    coxph(as.formula(paste("Surv(time_out, event_out) ~", RHS)), data = sub_df),
    error = function(e) {
      cat(sprintf("    [%s] coxph error: %s\n", .label, conditionMessage(e))); NULL })
  if (is.null(f)) return(NULL)
  if (!"expose" %in% names(coef(f))) {
    cat(sprintf("    [%s] expose dropped from fit\n", .label)); return(NULL) }
  # Only require the expose coefficient (and its variance) to be finite —
  # rare covariate levels can produce non-finite nuisance coefs without
  # corrupting the expose estimate.
  exp_b  <- coef(f)[["expose"]]
  exp_v  <- vcov(f)["expose", "expose"]
  if (!is.finite(exp_b) || !is.finite(exp_v) || exp_v <= 0) {
    cat(sprintf("    [%s] non-finite expose coef/var\n", .label)); return(NULL) }
  exp_se <- sqrt(exp_v)
  z      <- exp_b / exp_se
  p      <- 2 * pnorm(-abs(z))
  data.frame(
    HR    = exp(exp_b),
    CI_L  = exp(exp_b - 1.96 * exp_se),
    CI_H  = exp(exp_b + 1.96 * exp_se),
    P_Value  = p,
    N        = nrow(sub_df),
    N_events = sum(sub_df$event_out)
  )
}

# ---------------------------------------------------------------------------
# All-cause: early (cap=0.5) + late (cap=2) at trial months 1..6
# ---------------------------------------------------------------------------
cat("\n[30m] All-cause (CC) target trial: 6 trial months × {early(0.5), late(2)} ...\n")
all_rows <- list()
for (m in 1:6) {
  tr <- build_trial(d, m)
  for (model in c("early", "late")) {
    cap <- if (model == "early") 0.5 else 2
    sub <- prep_window(tr, model, cap)
    res <- fit_one(sub, sprintf("m=%d %s/%s", m, "all", model))
    if (!is.null(res)) {
      res$Trial_Month <- paste0("Month_", m)
      res$model <- model; res$cap <- cap
      all_rows[[length(all_rows) + 1]] <- res
    }
    cat(sprintf("  m=%d  %-5s  HR=%s\n", m, model,
                if (is.null(res)) "—" else sprintf("%.2f (%.2f–%.2f)", res$HR, res$CI_L, res$CI_H)))
  }
}
all_df <- bind_rows(all_rows)
all_df <- all_df[, c("Trial_Month", "model", "cap", "HR", "CI_L", "CI_H",
                     "P_Value", "N", "N_events")]
all_path <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_cc_early_late_array.csv")
write.csv(all_df, all_path, row.names = FALSE)
cat(sprintf("\n[30m] Wrote %s\n", all_path))

# ---------------------------------------------------------------------------
# Cause-specific: late window (cap=2) for tb_hybrid, nontb_hybrid,
# tb_simonly, nontb_simonly + all (for parity)
# ---------------------------------------------------------------------------
cat("\n[30m] Cause-specific (CC) late mortality (cap=2) ...\n")
CAUSES <- c("all", "tb_hybrid", "nontb_hybrid", "tb_simonly", "nontb_simonly")
cs_rows <- list()
for (m in 1:6) {
  tr <- build_trial(d, m)
  for (cause in CAUSES) {
    sub <- prep_cause(tr, cause, cap = 2)
    res <- fit_one(sub, sprintf("m=%d cause=%s", m, cause))
    if (!is.null(res)) {
      res$Trial_Month <- paste0("Month_", m)
      res$cause <- cause; res$cap <- 2
      cs_rows[[length(cs_rows) + 1]] <- res
    }
    cat(sprintf("  m=%d  %-15s  HR=%s\n", m, cause,
                if (is.null(res)) "—" else sprintf("%.2f (%.2f–%.2f)", res$HR, res$CI_L, res$CI_H)))
  }
}
cs_df <- bind_rows(cs_rows)
cs_df <- cs_df[, c("Trial_Month", "cause", "cap", "HR", "CI_L", "CI_H",
                   "P_Value", "N", "N_events")]
cs_path <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_cc_cause_specific.csv")
write.csv(cs_df, cs_path, row.names = FALSE)
cat(sprintf("\n[30m] Wrote %s\n", cs_path))

# ---------------------------------------------------------------------------
# Side-by-side comparison vs MI primary
# ---------------------------------------------------------------------------
cat("\n[30m] Comparison: CC vs MI primary (late mortality, cap=2)\n")
mi_path <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_mi_early_late_array.csv")
if (file.exists(mi_path)) {
  mi <- read.csv(mi_path, stringsAsFactors = FALSE)
  mi_late <- mi[mi$model == "late" & mi$cap == 2, c("Trial_Month", "HR", "CI_L", "CI_H")]
  cc_late <- all_df[all_df$model == "late" & all_df$cap == 2,
                    c("Trial_Month", "HR", "CI_L", "CI_H", "N")]
  cmp <- merge(mi_late, cc_late, by = "Trial_Month", suffixes = c("_MI", "_CC"))
  cmp <- cmp[order(cmp$Trial_Month), ]
  cat("\n  Trial_Month   MI aHR (95% CI)        CC aHR (95% CI)        N_cc\n")
  for (i in seq_len(nrow(cmp))) {
    cat(sprintf("  %-12s  %5.2f (%.2f–%.2f)    %5.2f (%.2f–%.2f)    %d\n",
                cmp$Trial_Month[i], cmp$HR_MI[i], cmp$CI_L_MI[i], cmp$CI_H_MI[i],
                cmp$HR_CC[i], cmp$CI_L_CC[i], cmp$CI_H_CC[i], cmp$N[i]))
  }
}

cat("\n[30m] done.\n")
