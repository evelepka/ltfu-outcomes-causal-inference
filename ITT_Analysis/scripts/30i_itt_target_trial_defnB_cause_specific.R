# 30i. Cause-specific target-trial Cox under defn-B + grace eligibility
# ==============================================================================
# Mirrors 30e but uses defn-B (LTFU exposure shifted back 30 days, clamped
# to >=1 day for the ~4% violators).
#
# Output: target_trial_defnB_cause_specific.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(mice)
  library(broom)
})

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

GRACE_DAYS <- 30
GRACE_YRS  <- GRACE_DAYS / 365.25
SHIFT_YRS  <- 30 / 365.25

cat("[30i] Building cause attribution lookup from raw data...\n")
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

deathrec <- raw[!is.na(raw$dod) & !is.na(raw$cause_of_death_code) & nzchar(raw$cause_of_death_code), ]
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

attr_lookup$tb_hybrid    <- attr_lookup$cod_class %in% c("tb_strict", "tb_via_tbweb")
attr_lookup$nontb_hybrid <- attr_lookup$cod_class %in% c("non_tb", "ntb_via_tbweb")
attr_lookup$tb_broad     <- attr_lookup$cod_class %in% c("tb_strict", "tb_via_tbweb",
                                                          "respiratory", "hiv_other")
attr_lookup$tb_simonly    <- attr_lookup$cod_class == "tb_strict"
attr_lookup$nontb_simonly <- attr_lookup$cod_class == "non_tb"

cat(sprintf("  Attribution lookup built for %d individuals\n", nrow(attr_lookup)))

imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
stopifnot(length(imp_files) > 0)
M <- length(imp_files)

prepare_imp <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date_start <- as.Date(d$best_start)
  d$date_end   <- as.Date(d$end_date)
  d$tx_duration_yrs <- as.numeric(d$date_end - d$date_start) / 365.25
  d$true_tx_duration_yrs <- pmax(d$tx_duration_yrs - SHIFT_YRS, 1/365.25)
  d$age_group <- factor(d$age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
  for (v in c("sex", "race_clean", "edu_clean", "hiv_aids", "diabetes",
              "alcohol", "drug_use", "incarcerated", "homelessness",
              "hosp_admission", "clinical_clean", "dot_status")) {
    d[[v]] <- as.factor(d[[v]])
  }
  d <- merge(d, attr_lookup[, c("sinan_clean", "cod_class",
                                "tb_hybrid", "nontb_hybrid",
                                "tb_broad", "tb_simonly", "nontb_simonly")],
             by = "sinan_clean", all.x = TRUE)
  d
}
imp_list <- lapply(imp_files, prepare_imp)

COVARS <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids",
            "diabetes", "alcohol", "drug_use", "incarcerated",
            "homelessness", "hosp_admission", "clinical_clean", "dot_status")
RHS <- paste("expose +", paste(COVARS, collapse = " + "))

build_trial_defnB <- function(d, m) {
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

prep_cause <- function(tr, cause, cap = 2) {
  t <- tr$time_raw
  e <- tr$event_d_num
  if (cause == "all") {
    cause_event <- e
  } else {
    indicator <- tr[[cause]]
    indicator[is.na(indicator)] <- FALSE
    cause_event <- ifelse(e == 1 & indicator, 1, 0)
  }
  time_out  <- pmin(t, cap)
  event_out <- ifelse(t > 0.5 & t <= cap & cause_event == 1, 1, 0)
  tr$time_out  <- time_out
  tr$event_out <- event_out
  tr[tr$time_out > 0, , drop = FALSE]
}

fit_pooled <- function(m, cause, cap = 2) {
  fits <- lapply(imp_list, function(d) {
    tr <- build_trial_defnB(d, m)
    if (nrow(tr) < 50) return(NULL)
    tr <- prep_cause(tr, cause, cap)
    if (sum(tr$event_out) < 5) return(NULL)
    tryCatch(
      coxph(as.formula(paste("Surv(time_out, event_out) ~", RHS)), data = tr),
      error = function(e) NULL
    )
  })
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0) return(NULL)
  pooled <- tryCatch(
    summary(pool(as.mira(fits)), exponentiate = TRUE, conf.int = TRUE),
    error = function(e) NULL)
  if (is.null(pooled)) return(NULL)
  ex <- pooled[pooled$term == "expose", ]
  if (nrow(ex) == 0) return(NULL)
  data.frame(
    Trial_Month = paste0("Month_", m),
    cause       = cause, cap = cap,
    HR          = ex$estimate, CI_L = ex$conf.low, CI_H = ex$conf.high,
    P_Value     = ex$p.value, N_imp = length(fits)
  )
}

CAUSES <- c("all", "tb_hybrid", "nontb_hybrid", "tb_broad",
            "tb_simonly", "nontb_simonly")

rows <- list()
for (m in 1:6) {
  for (cause in CAUSES) {
    cat(sprintf("  trial=%d  cause=%-15s\n", m, cause))
    r <- fit_pooled(m, cause, cap = 2)
    if (!is.null(r)) rows[[length(rows) + 1]] <- r
  }
}

final_df <- bind_rows(rows)
out_path <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_cause_specific.csv")
write.csv(final_df, out_path, row.names = FALSE)
cat(sprintf("\n[30i] Wrote %d rows to %s\n", nrow(final_df), out_path))

# Summary
make_wide <- function(df, c_tb, c_ntb) {
  tb  <- df[df$cause == c_tb, ]
  ntb <- df[df$cause == c_ntb, ]
  data.frame(
    Trial_Month = tb$Trial_Month,
    HR_TB    = sprintf("%.2f (%.2f-%.2f)", tb$HR, tb$CI_L, tb$CI_H),
    HR_NonTB = sprintf("%.2f (%.2f-%.2f)", ntb$HR, ntb$CI_L, ntb$CI_H),
    P_TB     = signif(tb$P_Value, 3),
    P_NonTB  = signif(ntb$P_Value, 3)
  )
}
cat("\n[30i] Hybrid (SIM + TBweb) — late mortality, cap = 2 yr (defn B):\n")
print(make_wide(final_df, "tb_hybrid", "nontb_hybrid"))
cat("\n[30i] SIM-only — late mortality (defn B):\n")
print(make_wide(final_df, "tb_simonly", "nontb_simonly"))
