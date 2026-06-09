# 30o. Absolute-scale risk differences (all-cause) for the sequential landmark
# cohort analysis — companion to the hazard ratios in 30h (by month) and 32d
# (subgroups). Reports g-formula STANDARDISED 24-month all-cause mortality risk
# under LTFU vs continued treatment and the RISK DIFFERENCE (RD).
#
#   - By disengagement month m=1..6  (defn-B construction, matches Fig 4B)
#       -> target_trial_defnB_absolute_rd.csv
#   - By pre-specified subgroup level (32d stacked construction, matches Fig 4C)
#       -> target_trial_subgroup_absolute_rd.csv
#
# Standardisation: Breslow centred baseline cumulative hazard at 24 months +
# linear predictors, averaging individual predicted risks under expose=1 and
# expose=0 over the trial-cohort covariate distribution. Pooled across m=5
# imputations (mean of point estimates); CIs by nonparametric bootstrap
# rotating across imputations.
# ==============================================================================

suppressPackageStartupMessages({ library(dplyr); library(survival) })
.here <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
source(file.path(.here(), "_paths.R"))

GRACE_YRS <- 30/365.25; SHIFT_YRS <- 30/365.25; CAP <- 2
LATE_LO <- 0.5
MET <- c("risk_cont","risk_ltfu","rd","late_cont","late_ltfu","rd_late")
B_BY_MONTH <- 300; B_SUBGROUP <- 200
set.seed(42)

imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
M <- length(imp_files)
COVARS <- c("age_group","sex","race_clean","edu_clean","hiv_aids","diabetes",
            "alcohol","drug_use","incarcerated","homelessness","hosp_admission",
            "clinical_clean","dot_status")

prep <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$tx_duration_yrs <- as.numeric(as.Date(d$end_date) - as.Date(d$best_start))/365.25
  d$true_tx_duration_yrs <- pmax(d$tx_duration_yrs - SHIFT_YRS, 1/365.25)
  d$age_group <- factor(d$age_group, levels = c("15-24","25-44","45-64","≥65"))
  for (v in COVARS[-1]) d[[v]] <- as.factor(d[[v]])
  d
}
imp_list <- lapply(imp_files, prep)

# g-formula standardised risks + RD (all-cause), PERIOD-SPECIFIC hazards.
# Fits an early-window model (events 0-6mo) and a late-window model (events
# 6-24mo, early deaths censored, at-risk from origin) so the late exposure
# effect matches the late-window aHR (Fig 5C). Then, by g-formula:
#   S_early(6mo; a) = exp(-H_early(0.5) e^{lp_early,a})
#   late_inc(a)     = S_early(6mo; a) * (1 - exp(-(H_late(2)-H_late(0.5)) e^{lp_late,a}))
#   overall_risk(a) = (1 - S_early(6mo; a)) + late_inc(a)
# Marginal (not conditional on 6-mo survival) -> no landmark/collider selection.
.bh <- function(fit, tcut) { bz <- basehaz(fit, centered = TRUE)
  if (any(bz$time <= tcut)) max(bz$hazard[bz$time <= tcut]) else 0 }
std_rd <- function(tr, covars) {
  tr <- tr[tr$time_out > 0, , drop = FALSE]
  if (nrow(tr) < 50 || sum(tr$event_out) < 5) return(setNames(rep(NA_real_, 6), MET))
  rhs <- paste("expose +", paste(covars, collapse = " + "))
  tt <- tr$tr_; ee <- tr$ed
  tr$te <- pmin(tt, LATE_LO); tr$ee <- ifelse(tt <= LATE_LO & ee == 1, 1, 0)   # early events
  tr$tl <- pmin(tt, CAP);     tr$el <- ifelse(tt > LATE_LO & tt <= CAP & ee == 1, 1, 0)  # late events
  fit_e <- tryCatch(coxph(as.formula(paste("Surv(te,ee) ~", rhs)), data = tr), error = function(e) NULL)
  fit_l <- tryCatch(coxph(as.formula(paste("Surv(tl,el) ~", rhs)), data = tr), error = function(e) NULL)
  if (is.null(fit_e) || is.null(fit_l)) return(setNames(rep(NA_real_, 6), MET))
  He <- .bh(fit_e, LATE_LO)
  Hl <- .bh(fit_l, CAP) - .bh(fit_l, LATE_LO)
  nd0 <- tr; nd0$expose <- 0; nd1 <- tr; nd1$expose <- 1
  lpe0 <- predict(fit_e, newdata = nd0, type = "lp"); lpe1 <- predict(fit_e, newdata = nd1, type = "lp")
  lpl0 <- predict(fit_l, newdata = nd0, type = "lp"); lpl1 <- predict(fit_l, newdata = nd1, type = "lp")
  Se0 <- exp(-He * exp(lpe0)); Se1 <- exp(-He * exp(lpe1))
  late0 <- mean(Se0 * (1 - exp(-Hl * exp(lpl0))))
  late1 <- mean(Se1 * (1 - exp(-Hl * exp(lpl1))))
  r0_2 <- mean((1 - Se0) + Se0 * (1 - exp(-Hl * exp(lpl0))))
  r1_2 <- mean((1 - Se1) + Se1 * (1 - exp(-Hl * exp(lpl1))))
  # overall 0-24mo RD AND late-window (6-24mo) marginal RD (period-consistent
  # with the late HR; marginal, not conditional on 6-month survival)
  setNames(c(r0_2, r1_2, r1_2 - r0_2, late0, late1, late1 - late0), MET)
}

boot_ci <- function(builder, covars, B) {
  est <- rowMeans(sapply(imp_list, function(d) std_rd(builder(d), covars)), na.rm = TRUE)
  bb <- matrix(NA_real_, B, 6)
  for (b in 1:B) {
    d <- imp_list[[(b %% M) + 1]]; tr <- builder(d)
    idx <- sample.int(nrow(tr), nrow(tr), replace = TRUE)
    bb[b, ] <- std_rd(tr[idx, , drop = FALSE], covars)
  }
  ci <- apply(bb, 2, quantile, c(.025, .975), na.rm = TRUE)
  data.frame(metric = MET, est = est, lo = ci[1, ], hi = ci[2, ])
}

# ---- (1) By-month RD (defn-B; overall 24-mo) ----
month_builder <- function(m) function(d) {
  s <- (m-1)*30/365.25; e <- m*30/365.25; o <- s + GRACE_YRS
  d |>
    dplyr::filter(time_d_tx > o) |>
    dplyr::mutate(elig = ifelse(itt_group=="Non-LTFU" | true_tx_duration_yrs >= s, 1, 0)) |>
    dplyr::filter(elig == 1) |>
    dplyr::mutate(
      expose = ifelse(itt_group=="Loss to follow-up" & true_tx_duration_yrs>=s & true_tx_duration_yrs<e, 1, 0),
      ed = as.numeric(as.character(event_d)), tr_ = time_d_tx - o,
      time_out = pmin(tr_, CAP), event_out = ifelse(tr_ > CAP, 0, ed))
}
rows <- list()
for (m in 1:6) {
  cat(sprintf("[30o] by-month m=%d\n", m))
  df <- boot_ci(month_builder(m), COVARS, B_BY_MONTH)
  df$Trial_Month <- paste0("Month_", m); rows[[m]] <- df
}
bym <- bind_rows(rows)
write.csv(bym, file.path(ITT_RESULTS_DIR, "target_trial_defnB_absolute_rd.csv"), row.names = FALSE)
cat("[30o] wrote by-month RD\n")
cat("  overall 0-24mo RD:\n"); print(bym[bym$metric=="rd", c("Trial_Month","est","lo","hi")])
cat("  late 6-24mo RD:\n");    print(bym[bym$metric=="rd_late", c("Trial_Month","est","lo","hi")])

# ---- (2) Subgroup RD (32d stacked pooled construction; overall 24-mo) ----
build_pooled <- function(d) {
  tl <- lapply(1:6, function(m) {
    s <- (m-1)*30/365.25; e <- m*30/365.25; o <- s + GRACE_YRS
    d |>
      dplyr::filter(time_d_tx > o) |>
      dplyr::mutate(elig = ifelse(itt_group=="Non-LTFU" | tx_duration_yrs >= s, 1, 0)) |>
      dplyr::filter(elig == 1) |>
      dplyr::mutate(
        expose = ifelse(itt_group=="Loss to follow-up" & tx_duration_yrs>=s & tx_duration_yrs<e, 1, 0),
        trial_month = paste0("Month_", m),
        ed = as.numeric(as.character(event_d)), tr_ = time_d_tx - o,
        time_out = pmin(tr_, CAP), event_out = ifelse(tr_ > CAP, 0, ed))
  })
  out <- bind_rows(tl); out$trial_month <- factor(out$trial_month, levels=paste0("Month_",1:6)); out
}
pooled_list <- lapply(imp_list, build_pooled)
SUBGROUPS <- c("age_group","sex","hiv_aids","homelessness")
srows <- list()
for (sg in SUBGROUPS) {
  base_covs <- c(setdiff(COVARS, sg), "trial_month")
  lvls <- sort(unique(as.character(pooled_list[[1]][[sg]])))
  lvls <- lvls[!is.na(lvls) & lvls != ""]
  for (lvl in lvls) {
    cat(sprintf("[30o] subgroup %s = %s\n", sg, lvl))
    bld <- function(d) { p <- build_pooled(d); p[!is.na(p[[sg]]) & p[[sg]]==lvl, ] }
    # use precomputed pooled for point est speed
    est <- rowMeans(sapply(pooled_list, function(p) std_rd(p[!is.na(p[[sg]]) & p[[sg]]==lvl, ], base_covs)), na.rm=TRUE)
    bb <- matrix(NA_real_, B_SUBGROUP, 6)
    for (b in 1:B_SUBGROUP) {
      p <- pooled_list[[(b %% M)+1]]; sub <- p[!is.na(p[[sg]]) & p[[sg]]==lvl, ]
      idx <- sample.int(nrow(sub), nrow(sub), replace=TRUE)
      bb[b,] <- std_rd(sub[idx,,drop=FALSE], base_covs)
    }
    ci <- apply(bb,2,quantile,c(.025,.975),na.rm=TRUE)
    srows[[length(srows)+1]] <- data.frame(Subgroup=sg, Level=lvl,
      metric=MET, est=est, lo=ci[1,], hi=ci[2,])
  }
}
sub <- bind_rows(srows)
write.csv(sub, file.path(ITT_RESULTS_DIR, "target_trial_subgroup_absolute_rd.csv"), row.names = FALSE)
cat("[30o] wrote subgroup RD\n")
cat("  late 6-24mo RD by subgroup:\n")
print(sub[sub$metric=="rd_late", c("Subgroup","Level","est","lo","hi")])
