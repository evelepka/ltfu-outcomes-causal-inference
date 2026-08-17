# 30e. Cause-specific target-trial Cox under 30-day grace eligibility
# ==============================================================================
# Tests whether the LTFU effect on mortality is mediated through TB
# (the causal pathway) or is equally seen in non-TB deaths (suggesting
# residual confounding by social/clinical factors).
#
# Two complementary outcome definitions:
#
#   (A) Hybrid TB attribution (primary): combines SIM ICD-10 with
#       TBweb 'Obito TB'/'Obito NTB' case_outcome attribution. Maximizes
#       data use but cause-of-death attribution method differs by arm
#       (Non-LTFU has TBweb attribution at case closure; LTFU does not).
#
#   (B) SIM-only attribution (sensitivity): restricts to deaths with a
#       SIM ICD-10 code. Uniform attribution method across arms; loses
#       data, mostly Obito-on-treatment deaths that already fail the
#       grace eligibility filter anyway.
#
# Definitions used here:
#   TB death = SIM A15-A19 OR B90 OR B200 (HIV w/ mycobacteria)
#              OR (in hybrid) TBweb case_outcome == 'Obito TB'
#   Non-TB death = SIM ICD non-TB-non-respiratory-non-HIV
#              OR (in hybrid) TBweb case_outcome == 'Obito NTB'
#   Unknown = death without SIM code and without TBweb Obito attribution;
#             treated as censored at death time (negative-control direction)
#
# Cause-specific Cox: for the TB outcome, non-TB deaths are censored at
# their actual event time (event = 0); for the non-TB outcome, TB deaths
# are likewise censored. Unknown deaths are censored in BOTH.
#
# Grace eligibility: alive at month m + 30 days; time-at-risk origin
# shifted to month m + 30 days (mirrors 30d).
#
# Output:
#   target_trial_grace_cause_specific.csv
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

# ----------------------------------------------------------------------------
# Build cause attribution lookup from raw data
# ----------------------------------------------------------------------------
cat("[30e] Building cause attribution lookup from raw data...\n")
raw_path <- file.path(DATA_DIR, "Final_table_cleaned.csv")
raw <- read.csv(raw_path, stringsAsFactors = FALSE)
# Raw file uses "August 28, 2015"-style dates
parse_dt <- function(x) as.Date(x, format = "%B %d, %Y")
raw$end_date <- parse_dt(raw$end_date)
raw$dod      <- parse_dt(raw$dod)

# Recover index case_outcome (first valid Novo per individual)
TRANSFER <- c("Transf Outro Municipio", "Transf Outro Estado/Pais")
novo <- raw[trimws(tools::toTitleCase(tolower(raw$case_type))) == "Novo", ]
novo <- novo[!is.na(novo$case_outcome)
             & nzchar(trimws(novo$case_outcome))
             & novo$case_outcome != "Mud Diag"
             & !novo$case_outcome %in% TRANSFER, ]
novo <- novo[order(novo$end_date), ]
first <- novo[!duplicated(novo$sinan_clean), c("sinan_clean", "case_outcome")]
# --- FIX (2026-08-16): take the Obito outcome from ANY episode ---------------
# The block above finds the INDEX (first Novo) episode outcome. An LTFU
# patient's index episode closes as `Abandono`, which can never be
# `Obito TB`/`Obito NTB` -- so index-only lookup discarded TBweb cause of
# death for 1,058 of 1,668 LTFU deaths (63.4%), all recorded on their
# retreatment episode (`Retr Aband`, `Recidiva`). Verified same-death:
# |Obito episode end_date - cohort death_date| median 0 d, 99.9% within 30 d.
# Effect of the defect: the LTFU arm lost ~45% of its deaths to
# cause-specific censoring vs ~5% in the non-LTFU arm, deflating every
# TB-specific hazard ratio.
# SIM precedence is unchanged (classify_cod only consults case_outcome when
# the SIM code is absent), so this only ADDS attribution where there was none.
obito <- raw[!is.na(raw$case_outcome)
             & trimws(raw$case_outcome) %in% c("Obito TB", "Obito NTB"), ]
obito <- obito[order(obito$end_date), ]
obito <- obito[!duplicated(obito$sinan_clean, fromLast = TRUE),
               c("sinan_clean", "case_outcome")]
names(obito)[2] <- "obito_outcome"
first <- merge(first, obito, by = "sinan_clean", all = TRUE)
first$case_outcome <- ifelse(!is.na(first$obito_outcome),
                             first$obito_outcome, first$case_outcome)
first$obito_outcome <- NULL
cat(sprintf("[cause-fix] Obito outcome recovered from any episode for %d individuals\n",
            nrow(obito)))


# Last record with both dod and cause_of_death_code per individual
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
  # Final classes
  cls <- rep("unknown", length(cod))
  cls[cod_known & !tb_strict & !resp & !hiv_other] <- "non_tb"
  cls[hiv_other] <- "hiv_other"
  cls[resp]      <- "respiratory"
  cls[tb_strict] <- "tb_strict"
  # TBweb attribution where SIM is missing
  miss <- cls == "unknown"
  cls[miss & !is.na(case_outcome) & case_outcome == "Obito TB"]  <- "tb_via_tbweb"
  cls[miss & !is.na(case_outcome) & case_outcome == "Obito NTB"] <- "ntb_via_tbweb"
  cls
}
attr_lookup$cod_class <- classify_cod(attr_lookup$cause_of_death_code,
                                      attr_lookup$case_outcome)

# Hybrid attribution
attr_lookup$tb_hybrid    <- attr_lookup$cod_class %in% c("tb_strict", "tb_via_tbweb")
attr_lookup$nontb_hybrid <- attr_lookup$cod_class %in% c("non_tb", "ntb_via_tbweb")
attr_lookup$unk_hybrid   <- attr_lookup$cod_class == "unknown"
# We treat respiratory + hiv_other as TB-broad in a third variant
attr_lookup$tb_broad     <- attr_lookup$cod_class %in% c("tb_strict", "tb_via_tbweb",
                                                          "respiratory", "hiv_other")
# SIM-only attribution
attr_lookup$tb_simonly    <- attr_lookup$cod_class == "tb_strict"
attr_lookup$nontb_simonly <- attr_lookup$cod_class == "non_tb"

cat(sprintf("  Attribution lookup built for %d individuals\n", nrow(attr_lookup)))

# ----------------------------------------------------------------------------
# Load imputed datasets and merge in attribution
# ----------------------------------------------------------------------------
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
stopifnot(length(imp_files) > 0)
M <- length(imp_files)
cat(sprintf("[30e] %d imputed datasets; grace = %d d\n", M, GRACE_DAYS))

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
  d <- merge(d, attr_lookup[, c("sinan_clean", "cod_class",
                                "tb_hybrid", "nontb_hybrid", "unk_hybrid",
                                "tb_broad", "tb_simonly", "nontb_simonly")],
             by = "sinan_clean", all.x = TRUE)
  d
}
imp_list <- lapply(imp_files, prepare_imp)

COVARS <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids",
            "diabetes", "alcohol", "drug_use", "incarcerated",
            "homelessness", "hosp_admission", "clinical_clean", "dot_status")
RHS <- paste("expose +", paste(COVARS, collapse = " + "))

# Build trial with grace eligibility (mirrors 30d)
build_trial_grace <- function(d, m) {
  start_yrs   <- (m - 1) * 30 / 365.25
  end_yrs     <-  m      * 30 / 365.25
  origin_yrs  <- start_yrs + GRACE_YRS
  d |>
    dplyr::filter(time_d_tx > origin_yrs) |>
    dplyr::mutate(
      eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
    ) |>
    dplyr::filter(eligible == 1) |>
    dplyr::mutate(
      expose = ifelse(itt_group == "Loss to follow-up" &
                        tx_duration_yrs >= start_yrs &
                        tx_duration_yrs <  end_yrs, 1, 0),
      event_d_num = as.numeric(as.character(event_d)),
      time_raw    = time_d_tx - origin_yrs
    )
}

# Cause-specific outcome configurator.
# `cause` is one of: "tb_hybrid", "nontb_hybrid", "tb_broad",
#                    "tb_simonly", "nontb_simonly", "all"
# For cause-specific Cox, only deaths from the target cause count as events;
# deaths from other causes are censored at their death time.
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
  # late window: events in (0.5, cap]
  time_out  <- pmin(t, cap)
  event_out <- ifelse(t > 0.5 & t <= cap & cause_event == 1, 1, 0)
  tr$time_out  <- time_out
  tr$event_out <- event_out
  tr[tr$time_out > 0, , drop = FALSE]
}

fit_pooled <- function(m, cause, cap = 2) {
  fits <- lapply(imp_list, function(d) {
    tr <- build_trial_grace(d, m)
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
    cause       = cause,
    cap         = cap,
    HR          = ex$estimate,
    CI_L        = ex$conf.low,
    CI_H        = ex$conf.high,
    P_Value     = ex$p.value,
    N_imp       = length(fits)
  )
}

# Cause definitions to estimate, for the late window (cap=2 yr)
CAUSES <- c(
  "all",            # all-cause (sanity: should match 30d late cap=2)
  "tb_hybrid",      # primary: SIM TB OR TBweb Obito TB
  "nontb_hybrid",   # primary negative control: SIM non-TB OR TBweb Obito NTB
  "tb_broad",       # broad: TB-strict + respiratory + HIV-related
  "tb_simonly",     # sensitivity: SIM TB only
  "nontb_simonly"   # sensitivity: SIM non-TB only
)

rows <- list()
for (m in 1:6) {
  for (cause in CAUSES) {
    cat(sprintf("  trial=%d  cause=%-15s  cap=2.0\n", m, cause))
    r <- fit_pooled(m, cause, cap = 2)
    if (!is.null(r)) rows[[length(rows) + 1]] <- r
  }
}

final_df <- bind_rows(rows)
out_path <- file.path(ITT_RESULTS_DIR, "target_trial_grace_cause_specific_fixedattr.csv")
write.csv(final_df, out_path, row.names = FALSE)
cat(sprintf("\n[30e] Wrote %d rows to %s\n", nrow(final_df), out_path))

# ----------------------------------------------------------------------------
# Print compact summary: TB vs Non-TB side-by-side
# ----------------------------------------------------------------------------
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

cat("\n[30e] Hybrid (SIM + TBweb attribution) — late mortality, cap = 2 yr:\n")
print(make_wide(final_df, "tb_hybrid", "nontb_hybrid"))

cat("\n[30e] SIM-only (sensitivity) — late mortality, cap = 2 yr:\n")
print(make_wide(final_df, "tb_simonly", "nontb_simonly"))

# Counts of events for context
cat("\n[30e] Event counts in trial m=2 (representative), late window:\n")
d1 <- imp_list[[1]]
tr <- build_trial_grace(d1, 2)
cat(sprintf("  trial 2 N (after grace eligibility): %d\n", nrow(tr)))
late_mask <- tr$time_raw > 0.5 & tr$time_raw <= 2
cat(sprintf("  late-window deaths (any cause): %d\n", sum(tr$event_d_num[late_mask] == 1)))
for (cause in CAUSES[-1]) {
  ind <- tr[[cause]]; ind[is.na(ind)] <- FALSE
  n_ev <- sum(late_mask & tr$event_d_num == 1 & ind)
  cat(sprintf("  late-window %s deaths: %d\n", cause, n_ev))
}
