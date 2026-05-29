# 57c. Bayesian Cox PH for LTFU vs on-treatment at m=3 landmark
# ==============================================================
# Frequentist Cox (script 57b) didn't converge due to the highly imbalanced
# cohort (138,626 on-tx vs 3,689 LTFU) combined with many one-hot covariates.
# Bayesian Cox via brms (HMC/NUTS) with weakly informative priors should be
# more robust under sparse / quasi-separated data.
#
# Two models per spec (ITT, IPCW):
#   crude     : ltfu only
#   adjusted  : ltfu + key baseline covariates (parsimonious set)
#
# IPCW model uses subject-level weights via `| weights(w)`.
#
# Outputs:
#   ITT_Analysis/results/bayesian_cox_ltfu_vs_ontx_m3.csv
#   ITT_Analysis/results/bayesian_cox_ltfu_vs_ontx_m3.rds

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(readr)
})

BASE <- "/Users/jasonandrews/Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/.shortcut-targets-by-id/18HafZqxrHeLVzpA6oYW-1cro9ygNfwVd/Abandonment Paper"
ITT_PATH  <- file.path(BASE, "ITT_Analysis/data/analytic_itt_m3.csv")
IPCW_PATH <- file.path(BASE, "ITT_Analysis/data/analytic_ipcw_m3.csv")
OUT_CSV   <- file.path(BASE, "ITT_Analysis/results/bayesian_cox_ltfu_vs_ontx_m3.csv")
OUT_RDS   <- file.path(BASE, "ITT_Analysis/results/bayesian_cox_ltfu_vs_ontx_m3.rds")

# parsimonious covariate set (chosen for clinical relevance + identifiability)
COVARS <- c("age_tb","sex","hiv_aids","hosp_admission","clinical_clean",
            "homelessness","alcohol","drug_use","dot_status")

prep <- function(d) {
  d <- d |>
    filter(t_d > 0) |>
    mutate(
      cens = 1 - ev_d,                  # brms cens(): 1 = right-censored
      ltfu = as.integer(ltfu),
      across(all_of(setdiff(COVARS, "age_tb")), ~ factor(.x))
    )
  d <- d[complete.cases(d[, c("t_d","ev_d","ltfu", COVARS, "w")]), ]
  d
}

itt_full  <- prep(read_csv(ITT_PATH,  show_col_types = FALSE))
ipcw_full <- prep(read_csv(IPCW_PATH, show_col_types = FALSE))

# Subsample the on-treatment arm 5:1 to keep brms cox tractable.
# LTFU coefficient is unbiased by subsampling on-tx (loses precision only).
set.seed(20260528)
subsample_ontx <- function(d, ratio = 5) {
  ltfu <- d[d$ltfu==1, ]
  ontx <- d[d$ltfu==0, ]
  n_take <- min(nrow(ontx), nrow(ltfu) * ratio)
  ontx_sub <- ontx[sample.int(nrow(ontx), n_take), ]
  rbind(ltfu, ontx_sub)
}
itt  <- subsample_ontx(itt_full)
ipcw <- subsample_ontx(ipcw_full)
cat(sprintf("ITT  N = %s (events=%s; LTFU events=%s, on-tx events=%s)\n",
            format(nrow(itt),  big.mark=","), sum(itt$ev_d),
            sum(itt$ev_d[itt$ltfu==1]), sum(itt$ev_d[itt$ltfu==0])))
cat(sprintf("IPCW N = %s (events=%s; LTFU events=%s, on-tx events=%s)\n",
            format(nrow(ipcw), big.mark=","), sum(ipcw$ev_d),
            sum(ipcw$ev_d[ipcw$ltfu==1]), sum(ipcw$ev_d[ipcw$ltfu==0])))

# weakly informative priors for log-HRs
priors <- prior(normal(0, 1), class = "b")

# Two models per spec: crude and adjusted
formula_crude    <- bf(t_d | cens(cens) ~ ltfu)
formula_adjusted <- bf(t_d | cens(cens) ~ ltfu + age_tb + sex + hiv_aids +
                       hosp_admission + clinical_clean + homelessness +
                       alcohol + drug_use + dot_status)

# IPCW versions add weights
formula_crude_w    <- bf(t_d | cens(cens) + weights(w) ~ ltfu)
formula_adjusted_w <- bf(t_d | cens(cens) + weights(w) ~ ltfu + age_tb + sex + hiv_aids +
                          hosp_admission + clinical_clean + homelessness +
                          alcohol + drug_use + dot_status)

fit_one <- function(d, fm, name) {
  cat(sprintf("\n[fit] %s ...\n", name))
  t0 <- Sys.time()
  mod <- brm(formula = fm,
             data = d,
             family = brmsfamily("cox"),
             prior = priors,
             chains = 2, iter = 1000, warmup = 500,
             cores = 2, seed = 20260528, refresh = 200,
             control = list(adapt_delta = 0.9))
  cat(sprintf("  done in %.1f min\n", as.numeric(Sys.time()-t0, units="mins")))
  mod
}

models <- list(
  itt_crude     = fit_one(itt,  formula_crude,    "ITT crude"),
  itt_adj       = fit_one(itt,  formula_adjusted, "ITT adjusted"),
  ipcw_crude    = fit_one(ipcw, formula_crude_w,  "IPCW crude (weighted)"),
  ipcw_adj      = fit_one(ipcw, formula_adjusted_w,"IPCW adjusted (weighted)")
)

summarise_ltfu <- function(mod, name) {
  s <- as_draws_df(mod, variable = "b_ltfu")
  hr <- exp(s$b_ltfu)
  data.frame(
    model         = name,
    posterior_median_HR = round(median(hr), 3),
    ci_lo_HR      = round(quantile(hr, 0.025), 3),
    ci_hi_HR      = round(quantile(hr, 0.975), 3),
    p_HR_gt_1     = round(mean(hr > 1), 3)
  )
}

results <- do.call(rbind, list(
  summarise_ltfu(models$itt_crude,   "ITT crude"),
  summarise_ltfu(models$itt_adj,     "ITT adjusted"),
  summarise_ltfu(models$ipcw_crude,  "IPCW crude"),
  summarise_ltfu(models$ipcw_adj,    "IPCW adjusted")
))
print(results)
write.csv(results, OUT_CSV, row.names = FALSE)
saveRDS(models, OUT_RDS)
cat(sprintf("\nwrote %s\n", OUT_CSV))
