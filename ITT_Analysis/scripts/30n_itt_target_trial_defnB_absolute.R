# 30n. Absolute-scale estimands for the sequential landmark cohort (defn-B + grace)
# ==============================================================================
# Companion to 30h (hazard ratios) and 30i (cause-specific HRs).  Reports, by
# disengagement month m (1-6), the g-formula STANDARDISED 24-month risk under
# LTFU vs continued treatment and the RISK DIFFERENCE (RD), on the same trial
# cohorts used for the HR analysis.  Requested for the epidemiology audience
# because period-specific HRs carry a built-in selection component and are not
# clean causal contrasts (Hernan 2010; Stensrud & Hernan 2020).
#
# Measures (24-month horizon from the post-grace landmark):
#   - all-cause standardised mortality risk (LTFU, continued) + RD
#   - cause-specific cumulative incidence (TB-attributable, non-TB) + RD,
#     via competing-risks g-formula (Aalen-Johansen-type integration of the
#     cause-specific hazards over the standardised survival curve)
#
# Standardisation: Breslow baseline cumulative hazard + linear predictors,
# averaging individual predicted risks over the trial-cohort covariate
# distribution under expose=1 and expose=0 (g-formula).
# Pooled across m=5 imputations (mean of point estimates).
# CIs: nonparametric bootstrap, rotating across imputations to propagate both
# sampling and imputation uncertainty.
#
# Output: target_trial_defnB_absolute_rd.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

.here <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
source(file.path(.here(), "_paths.R"))

GRACE_YRS <- 30 / 365.25
SHIFT_YRS <- 30 / 365.25
CAP       <- 2          # 24-month horizon
B_BOOT    <- 300        # bootstrap reps for CIs
set.seed(42)

# ---- cause attribution lookup (identical logic to 30i) ----------------------
cat("[30n] Building cause attribution lookup from raw data...\n")
raw <- read.csv(file.path(DATA_DIR, "Final_table_cleaned.csv"), stringsAsFactors = FALSE)
parse_dt <- function(x) as.Date(x, format = "%B %d, %Y")
raw$end_date <- parse_dt(raw$end_date); raw$dod <- parse_dt(raw$dod)
TRANSFER <- c("Transf Outro Municipio", "Transf Outro Estado/Pais")
novo <- raw[trimws(tools::toTitleCase(tolower(raw$case_type))) == "Novo", ]
novo <- novo[!is.na(novo$case_outcome) & nzchar(trimws(novo$case_outcome)) &
               novo$case_outcome != "Mud Diag" & !novo$case_outcome %in% TRANSFER, ]
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
  tb_strict <- cod_known & (grepl("^A1[5-9]", cod) | grepl("^B90", cod) | grepl("^B200", cod))
  resp <- cod_known & grepl("^J[0-9]", cod)
  hiv_other <- cod_known & grepl("^B2[0-4]", cod) & !grepl("^B200", cod)
  cls <- rep("unknown", length(cod))
  cls[cod_known & !tb_strict & !resp & !hiv_other] <- "non_tb"
  cls[hiv_other] <- "hiv_other"; cls[resp] <- "respiratory"; cls[tb_strict] <- "tb_strict"
  miss <- cls == "unknown"
  cls[miss & !is.na(case_outcome) & case_outcome == "Obito TB"]  <- "tb_via_tbweb"
  cls[miss & !is.na(case_outcome) & case_outcome == "Obito NTB"] <- "ntb_via_tbweb"
  cls
}
attr_lookup$cod_class    <- classify_cod(attr_lookup$cause_of_death_code, attr_lookup$case_outcome)
attr_lookup$tb_hybrid    <- attr_lookup$cod_class %in% c("tb_strict", "tb_via_tbweb")
attr_lookup$nontb_hybrid <- attr_lookup$cod_class %in% c("non_tb", "ntb_via_tbweb")

# ---- imputations ------------------------------------------------------------
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
M <- length(imp_files)
COVARS <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids",
            "diabetes", "alcohol", "drug_use", "incarcerated",
            "homelessness", "hosp_admission", "clinical_clean", "dot_status")
RHS <- paste("expose +", paste(COVARS, collapse = " + "))

prepare_imp <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date_start <- as.Date(d$best_start); d$date_end <- as.Date(d$end_date)
  d$tx_duration_yrs <- as.numeric(d$date_end - d$date_start) / 365.25
  d$true_tx_duration_yrs <- pmax(d$tx_duration_yrs - SHIFT_YRS, 1/365.25)
  d$age_group <- factor(d$age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
  for (v in COVARS[-1]) d[[v]] <- as.factor(d[[v]])
  merge(d, attr_lookup[, c("sinan_clean", "tb_hybrid", "nontb_hybrid")],
        by = "sinan_clean", all.x = TRUE)
}
imp_list <- lapply(imp_files, prepare_imp)

build_trial <- function(d, m) {
  start_yrs <- (m - 1) * 30 / 365.25; end_yrs <- m * 30 / 365.25
  origin_yrs <- start_yrs + GRACE_YRS
  d |>
    dplyr::filter(time_d_tx > origin_yrs) |>
    dplyr::mutate(eligible = ifelse(itt_group == "Non-LTFU" | true_tx_duration_yrs >= start_yrs, 1, 0)) |>
    dplyr::filter(eligible == 1) |>
    dplyr::mutate(
      expose      = ifelse(itt_group == "Loss to follow-up" &
                             true_tx_duration_yrs >= start_yrs &
                             true_tx_duration_yrs <  end_yrs, 1, 0),
      event_d_num = as.numeric(as.character(event_d)),
      time_raw    = pmin(time_d_tx - origin_yrs, CAP),
      tb_ev       = ifelse(event_d_num == 1 & !is.na(tb_hybrid)    & tb_hybrid,    1, 0),
      ntb_ev      = ifelse(event_d_num == 1 & !is.na(nontb_hybrid) & nontb_hybrid, 1, 0)
    ) |>
    dplyr::filter(time_raw > 0)
}

# g-formula standardised risk at CAP from a cause-specific Cox fit.
# Returns per-individual cumulative hazard contributions on the cohort.
# We compute all-cause survival from the all-cause hazard, and cause-specific
# CIF by integrating cause hazard increments against all-cause survival.
std_risks <- function(tr) {
  # design matrices for expose=0 and expose=1 (covariates held at observed)
  mm0 <- tr; mm1 <- tr; mm0$expose <- 0; mm1$expose <- 1
  fit_all <- coxph(as.formula(paste("Surv(time_raw, event_d_num) ~", RHS)), data = tr, x = FALSE)
  fit_tb  <- coxph(as.formula(paste("Surv(time_raw, tb_ev) ~",  RHS)), data = tr)
  fit_ntb <- coxph(as.formula(paste("Surv(time_raw, ntb_ev) ~", RHS)), data = tr)

  # baseline cumulative hazard (centered) on a monthly grid to CAP
  grid <- seq(0, CAP, by = 1/12)
  bh <- function(fit) {
    bz <- basehaz(fit, centered = TRUE)
    approx(bz$time, bz$hazard, xout = grid, method = "constant",
           yleft = 0, rule = 2)$y
  }
  H_all <- bh(fit_all); H_tb <- bh(fit_tb); H_ntb <- bh(fit_ntb)
  dH_all <- diff(H_all); dH_tb <- diff(H_tb); dH_ntb <- diff(H_ntb)

  lp <- function(fit, newd) as.numeric(predict(fit, newdata = newd, type = "lp"))
  out <- list()
  for (a in c(0, 1)) {
    nd <- if (a == 0) mm0 else mm1
    r_all <- exp(lp(fit_all, nd)); r_tb <- exp(lp(fit_tb, nd)); r_ntb <- exp(lp(fit_ntb, nd))
    # all-cause survival over grid: S_i(t) = exp(-H_all(t)*r_all_i)
    # cause CIF_i(CAP) = sum_t S_i(t_-) * (dH_cause(t)*r_cause_i)
    n <- nrow(nd)
    S_prev <- rep(1, n)
    cif_tb <- rep(0, n); cif_ntb <- rep(0, n)
    Hcum <- rep(0, n)
    for (k in seq_along(dH_all)) {
      S_prev <- exp(-Hcum * r_all)
      cif_tb  <- cif_tb  + S_prev * (dH_tb[k]  * r_tb)
      cif_ntb <- cif_ntb + S_prev * (dH_ntb[k] * r_ntb)
      Hcum <- Hcum + dH_all[k] * r_all
    }
    risk_all <- 1 - exp(-H_all[length(H_all)] * r_all)
    out[[as.character(a)]] <- c(all = mean(risk_all), tb = mean(cif_tb), ntb = mean(cif_ntb))
  }
  c(
    risk_all_cont = out[["0"]]["all"], risk_all_ltfu = out[["1"]]["all"],
    rd_all  = out[["1"]]["all"] - out[["0"]]["all"],
    cif_tb_cont = out[["0"]]["tb"], cif_tb_ltfu = out[["1"]]["tb"],
    rd_tb  = out[["1"]]["tb"] - out[["0"]]["tb"],
    cif_ntb_cont = out[["0"]]["ntb"], cif_ntb_ltfu = out[["1"]]["ntb"],
    rd_ntb = out[["1"]]["ntb"] - out[["0"]]["ntb"]
  )
}

rows <- list()
for (m in 1:6) {
  cat(sprintf("[30n] month %d: point estimate across %d imputations...\n", m, M))
  pt <- rowMeans(sapply(imp_list, function(d) std_risks(build_trial(d, m))))
  # bootstrap CI (rotate imputations)
  cat(sprintf("        bootstrap B=%d...\n", B_BOOT))
  boot <- matrix(NA_real_, nrow = B_BOOT, ncol = length(pt))
  for (b in 1:B_BOOT) {
    d <- imp_list[[(b %% M) + 1]]
    tr <- build_trial(d, m)
    idx <- sample.int(nrow(tr), nrow(tr), replace = TRUE)
    boot[b, ] <- tryCatch(std_risks(tr[idx, , drop = FALSE]), error = function(e) rep(NA, length(pt)))
  }
  ci <- apply(boot, 2, function(x) quantile(x, c(.025, .975), na.rm = TRUE))
  df <- data.frame(Trial_Month = paste0("Month_", m), metric = names(pt),
                   est = as.numeric(pt), lo = ci[1, ], hi = ci[2, ])
  rows[[m]] <- df
}
final <- bind_rows(rows)
out_path <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_absolute_rd.csv")
write.csv(final, out_path, row.names = FALSE)
cat(sprintf("\n[30n] Wrote %d rows to %s\n", nrow(final), out_path))
# headline print
hl <- final[final$metric %in% c("rd_all.all", "rd_tb.tb", "rd_ntb.ntb"), ]
hl$pct <- sprintf("%.1f%% (%.1f to %.1f)", 100*hl$est, 100*hl$lo, 100*hl$hi)
print(hl[, c("Trial_Month", "metric", "pct")])
