# 52b. Figure 3 (defn-B primary) — Causal mortality panels
# ==============================================================================
# Mirrors 52 (and inherits its polish) but uses defn-B primary inputs.
#
# Panel A change (2026-05-05): the time-varying HR(t) is now produced as
# per-month-bin Cox fits per imputation, each with the same 13-covariate
# adjustment set used by the target trial (30h). The earlier attempt to use
# a single spline tt(expose) Cox hit R's 24 GB memory limit on the 193K-row
# counting-process frame; per-bin Cox fits avoid the tt() expansion entirely.
# Per-imputation log-HR per bin is pooled across M=5 imputations with Rubin's
# rules. This isolates the time-handling artefact (immortal-time bias) from
# baseline confounding — naïve unadjusted Cox conflates the two.
#
# The unadjusted per-bin curve (the prior implementation) is retained as a
# supplementary figure (Figure_3_supp_unadj_vs_adj_HRt.{png,pdf}) so the contrast
# is visible.
#
# Inputs:
#   ITT_MI_DIR/imp_*.csv                              (5 imputed cohorts)
#   COHORT_CSV                                        (for unadj supp curve)
#   ITT_RESULTS_DIR/target_trial_defnB_mi_early_late_array.csv  (panel B, from 30h)
#   ITT_RESULTS_DIR/target_trial_defnB_subgroups_mi.csv         (panel C, from 32e)
#
# Outputs:
#   Figure_3_causal_mortality_defnB.{png,pdf}          (primary, adj spline pA)
#   Figure_3a_HR_over_time_defnB.csv                   (unadj per-bin, kept)
#   Figure_3a_HR_over_time_defnB_adjusted.csv          (MI-pooled spline curve)
#   Figure_3_supp_unadj_vs_adj_HRt.{png,pdf}           (supp comparison)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(splines)
  library(patchwork)
  library(scales)
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

OUT_PNG       <- file.path(ITT_RESULTS_DIR, "Figure_3_causal_mortality_defnB.png")
OUT_PDF       <- file.path(ITT_RESULTS_DIR, "Figure_3_causal_mortality_defnB.pdf")
SUPP_PNG      <- file.path(ITT_RESULTS_DIR, "Figure_3_supp_unadj_vs_adj_HRt.png")
SUPP_PDF      <- file.path(ITT_RESULTS_DIR, "Figure_3_supp_unadj_vs_adj_HRt.pdf")
ADJ_CSV       <- file.path(ITT_RESULTS_DIR, "Figure_3a_HR_over_time_defnB_adjusted.csv")
UNADJ_CSV     <- file.path(ITT_RESULTS_DIR, "Figure_3a_HR_over_time_defnB.csv")

SHIFT_DAYS <- 30
SHIFT_YRS  <- SHIFT_DAYS / 365.25
HORIZON    <- 2.0

COVARS <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids",
            "diabetes", "alcohol", "drug_use", "incarcerated",
            "homelessness", "hosp_admission", "clinical_clean", "dot_status")

# Optional: limit number of imputations for quick sanity testing.
# Set to NULL to use all.
M_LIMIT <- NULL

# Bin width (months) for the per-bin time-varying HR. 1 month is the original
# resolution; 2 months is more robust early on.
BIN_WIDTH_MO <- 1
MIN_EVENTS_PER_BIN <- 5

# ---------------------------------------------------------------------------
# Panel A (PRIMARY): adjusted per-bin HR(t) with MI pooling
# ---------------------------------------------------------------------------
cat("[fig3-defnB] Panel A: adjusted per-bin HR(t) with MI pooling ...\n")

prepare_imp_data <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$age_group <- factor(d$age_group,
                        levels = c("15-24", "25-44", "45-64", "≥65"))
  for (v in c("sex", "race_clean", "edu_clean", "hiv_aids", "diabetes",
              "alcohol", "drug_use", "incarcerated", "homelessness",
              "hosp_admission", "clinical_clean", "dot_status")) {
    d[[v]] <- as.factor(d[[v]])
  }
  d
}

# Build defn-B counting-process frame (mirrors the original 52b lines 43-70
# but operates on an imputed cohort).
build_counting_process <- function(d) {
  d$patient_id <- seq_len(nrow(d))
  d$tx_dur <- as.numeric(difftime(as.Date(d$end_date),
                                  as.Date(d$best_start), units = "days")) / 365.25
  d$tx_dur_true <- pmax(d$tx_dur - SHIFT_YRS, 1/365.25)

  d$event_d   <- ifelse(d$time_d_tx > HORIZON, 0, d$event_d)
  d$time_d_tx <- ifelse(d$time_d_tx > HORIZON, HORIZON, d$time_d_tx)

  is_ltfu <- d$itt_group == "Loss to follow-up"

  d_p1 <- d
  d_p1$tstart <- 0
  d_p1$tstop  <- ifelse(is_ltfu & d_p1$tx_dur_true < d_p1$time_d_tx,
                        d_p1$tx_dur_true, d_p1$time_d_tx)
  d_p1$event  <- ifelse(is_ltfu & d_p1$tx_dur_true < d_p1$time_d_tx,
                        0, d_p1$event_d)
  d_p1$expose <- 0

  d_p2 <- d[is_ltfu & d$tx_dur_true < d$time_d_tx, ]
  d_p2$tstart <- d_p2$tx_dur_true
  d_p2$tstop  <- d_p2$time_d_tx
  d_p2$event  <- d_p2$event_d
  d_p2$expose <- 1

  bind_rows(d_p1, d_p2) |>
    dplyr::filter(round(tstop - tstart, 4) > 0)
}

# Fit one Cox in a single bin and return (logHR, SE) for the expose coefficient.
fit_bin_adjusted <- function(d_bin) {
  if (sum(d_bin$event) < MIN_EVENTS_PER_BIN ||
      length(unique(d_bin$expose)) < 2) return(NULL)
  rhs <- paste("expose +", paste(COVARS, collapse = " + "))
  f <- tryCatch(coxph(as.formula(paste("Surv(tstart, tstop, event) ~", rhs)),
                      data = d_bin, id = patient_id),
                error = function(e) NULL)
  if (is.null(f) || !"expose" %in% names(coef(f))) return(NULL)
  if (any(!is.finite(coef(f))) || any(!is.finite(diag(vcov(f))))) return(NULL)
  data.frame(logHR = coef(f)[["expose"]],
             SE    = sqrt(vcov(f)["expose", "expose"]),
             n_events = sum(d_bin$event))
}

# Per-imputation: build counting-process frame, split into monthly bins, fit one
# adjusted Cox per bin. Returns a long-form (imp, bin, logHR, SE) frame.
fit_per_bin_one_imp <- function(path, imp_idx) {
  d <- prepare_imp_data(path)
  d_split <- build_counting_process(d)
  BREAKS <- seq(0, HORIZON, by = (BIN_WIDTH_MO / 12))
  d_pw <- survSplit(Surv(tstart, tstop, event) ~ .,
                    data = d_split,
                    cut = BREAKS[-c(1, length(BREAKS))],
                    episode = "month_bin")
  d_pw$month_bin <- as.integer(d_pw$month_bin)
  bins <- sort(unique(d_pw$month_bin))
  out <- list()
  for (b in bins) {
    res <- fit_bin_adjusted(d_pw[d_pw$month_bin == b, ])
    if (!is.null(res)) {
      res$bin <- b; res$imp <- imp_idx
      out[[length(out) + 1]] <- res
    }
  }
  bind_rows(out)
}

# --- MI loop ---
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$",
                             full.names = TRUE))
stopifnot(length(imp_files) > 0)
if (!is.null(M_LIMIT)) imp_files <- imp_files[seq_len(min(M_LIMIT, length(imp_files)))]
M <- length(imp_files)
cat(sprintf("[fig3-defnB] Per-bin adjusted Cox on %d imputed cohort(s); bin = %d mo\n",
            M, BIN_WIDTH_MO))

curves_long <- list()
for (i in seq_along(imp_files)) {
  t0 <- Sys.time()
  res <- tryCatch(fit_per_bin_one_imp(imp_files[i], i),
                  error = function(e) {
                    cat(sprintf("  imp %d: failed: %s\n", i, conditionMessage(e)))
                    NULL
                  })
  if (is.null(res) || nrow(res) == 0) next
  curves_long[[length(curves_long) + 1]] <- res
  cat(sprintf("  imp %d/%d done in %.1f s; %d bins fit\n",
              i, M, as.numeric(difftime(Sys.time(), t0, units = "secs")),
              nrow(res)))
}
stopifnot(length(curves_long) > 0)
curves_long <- bind_rows(curves_long)

# --- Pool with Rubin's rules per bin ---
M_used <- length(unique(curves_long$imp))
pooled <- curves_long |>
  dplyr::group_by(bin) |>
  dplyr::summarise(
    Qbar = mean(logHR),
    Ubar = mean(SE^2),
    B    = if (dplyr::n() > 1) var(logHR) else 0,
    n_events_max = max(n_events),
    n_imp = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Tvar  = Ubar + (1 + 1/M_used) * B,
    HR    = exp(Qbar),
    CI_L  = exp(Qbar - 1.96 * sqrt(Tvar)),
    CI_H  = exp(Qbar + 1.96 * sqrt(Tvar)),
    t = (bin - 0.5) * (BIN_WIDTH_MO / 12)         # bin midpoint in years
  ) |>
  dplyr::filter(is.finite(HR), is.finite(CI_H), CI_H < 50)

pooled$month_mid <- pooled$t * 12
write.csv(pooled[, c("bin", "t", "month_mid", "HR", "CI_L", "CI_H",
                     "n_events_max", "n_imp")],
          ADJ_CSV, row.names = FALSE)
cat(sprintf("[fig3-defnB] Wrote %s (%d bins pooled)\n", ADJ_CSV, nrow(pooled)))

pA <- ggplot(pooled, aes(x = month_mid, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_L, ymax = CI_H),
                width = 0.10, linewidth = 0.5,
                color = "#c0392b", alpha = 0.5) +
  geom_smooth(method = "loess", span = 0.55, se = FALSE,
              color = "#c0392b", linewidth = 1.4) +
  geom_point(color = "#c0392b", size = 2.4) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                minor_breaks = c(0.33, 0.75, 1.5, 3, 6),
                limits = c(0.15, 10)) +
  scale_x_continuous(breaks = seq(0, HORIZON * 12, 3),
                     limits = c(0, HORIZON * 12)) +
  labs(title = "A. Adjusted time-varying hazard ratio",
       x = "Time since treatment start (months)",
       y = "Adjusted HR") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
        panel.grid.minor.y = element_line(color = "grey94", linewidth = 0.25))

# ---------------------------------------------------------------------------
# Unadjusted per-bin HR(t) — retained for the supplementary comparison
# ---------------------------------------------------------------------------
cat("[fig3-defnB] Unadjusted per-bin HR(t) (for supplementary) ...\n")
df <- read.csv(COHORT_CSV, stringsAsFactors = FALSE)
df$patient_id  <- seq_len(nrow(df))
df$tx_dur      <- as.numeric(difftime(as.Date(df$end_date),
                                      as.Date(df$best_start), units = "days")) / 365.25
df$tx_dur_true <- pmax(df$tx_dur - SHIFT_YRS, 1/365.25)
df$event_d     <- ifelse(df$time_d_tx > HORIZON, 0, df$event_d)
df$time_d_tx   <- ifelse(df$time_d_tx > HORIZON, HORIZON, df$time_d_tx)

df_p1 <- df
df_p1$tstart <- 0
df_p1$tstop  <- ifelse(df_p1$itt_group == "Loss to follow-up" &
                         df_p1$tx_dur_true < df_p1$time_d_tx,
                       df_p1$tx_dur_true, df_p1$time_d_tx)
df_p1$event <- ifelse(df_p1$itt_group == "Loss to follow-up" &
                        df_p1$tx_dur_true < df_p1$time_d_tx,
                      0, df_p1$event_d)
df_p1$expose <- 0

df_p2 <- df[df$itt_group == "Loss to follow-up" & df$tx_dur_true < df$time_d_tx, ]
df_p2$tstart <- df_p2$tx_dur_true
df_p2$tstop  <- df_p2$time_d_tx
df_p2$event  <- df_p2$event_d
df_p2$expose <- 1

df_split <- bind_rows(df_p1, df_p2) |>
  dplyr::filter(round(tstop - tstart, 4) > 0)

BREAKS <- seq(0, HORIZON, by = 1/12)
df_piecewise <- survSplit(Surv(tstart, tstop, event) ~ .,
                          data = df_split, cut = BREAKS[-c(1, length(BREAKS))],
                          episode = "month_bin")
df_piecewise$month_bin <- factor(df_piecewise$month_bin,
                                 levels = seq_len(length(BREAKS) - 1))

bucket_hr_unadj <- function(d, bin) {
  sub <- d[d$month_bin == bin, ]
  if (sum(sub$event) < 3 || length(unique(sub$expose)) < 2) return(NULL)
  f <- tryCatch(coxph(Surv(tstart, tstop, event) ~ expose,
                      data = sub, id = patient_id),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  s <- summary(f)
  data.frame(
    month_mid = (as.numeric(as.character(bin)) - 0.5) / 12,
    HR   = s$coefficients[1, "exp(coef)"],
    CI_L = s$conf.int[1, "lower .95"],
    CI_H = s$conf.int[1, "upper .95"],
    n_events = sum(sub$event)
  )
}
bins   <- levels(df_piecewise$month_bin)
hr_tbl <- do.call(rbind, lapply(bins, function(b) bucket_hr_unadj(df_piecewise, b)))
hr_tbl <- hr_tbl[is.finite(hr_tbl$HR) & is.finite(hr_tbl$CI_H) & hr_tbl$CI_H < 50, ]
write.csv(hr_tbl, UNADJ_CSV, row.names = FALSE)
cat(sprintf("[fig3-defnB] Unadjusted curve estimated at %d bins; wrote %s\n",
            nrow(hr_tbl), UNADJ_CSV))

# ---------------------------------------------------------------------------
# Panel B: Target-trial HR by month — Early vs Late (defn-B + grace)
# ---------------------------------------------------------------------------
cat("[fig3-defnB] Panel B: early-vs-late HR array (defn-B + grace) ...\n")
tt_csv <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_mi_early_late_array.csv")
if (!file.exists(tt_csv)) stop(sprintf("Missing %s (run 30h)", tt_csv))
dfB_raw <- read.csv(tt_csv, stringsAsFactors = FALSE)
dfB <- dfB_raw |>
  dplyr::filter((model == "early" & cap == 0.5) | (model == "late" & cap == 2)) |>
  dplyr::mutate(
    Month = as.numeric(gsub("Month_", "", Trial_Month)),
    Window = dplyr::recode(model,
                           "early" = "Early (0–6 months)",
                           "late"  = "Late (6–24 months)")
  )
dfB$Window <- factor(dfB$Window, levels = c("Early (0–6 months)", "Late (6–24 months)"))

palette_EL <- c("Early (0–6 months)" = "#3498db",
                "Late (6–24 months)" = "#e74c3c")

pB <- ggplot(dfB, aes(x = Month, y = HR, color = Window, group = Window)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_L, ymax = CI_H),
                width = 0.14, linewidth = 0.7) +
  geom_line(linewidth = 0.9, alpha = 0.8) +
  geom_point(size = 3.5) +
  scale_color_manual(values = palette_EL) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4),
                minor_breaks = c(0.33, 0.75, 1.5, 3)) +
  scale_x_continuous(breaks = 1:6, labels = as.character(1:6)) +
  labs(title = "B. Adjusted mortality HR by month of LTFU",
       x = "Month of LTFU",
       y = "AHR",
       color = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
        panel.grid.minor.y = element_line(color = "grey94", linewidth = 0.25))

# ---------------------------------------------------------------------------
# Panel C: Subgroup forest — LATE mortality (defn-B + grace, cap = 2)
# ---------------------------------------------------------------------------
cat("[fig3-defnB] Panel C: subgroup forest (defn-B + grace, 24mo cap) ...\n")
sub_csv <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_subgroups_mi.csv")
if (!file.exists(sub_csv)) stop(sprintf("Missing %s (run 32e)", sub_csv))
dfC_all <- read.csv(sub_csv, stringsAsFactors = FALSE)
dfC <- dfC_all |> dplyr::filter(model == "late", cap == 2)

make_rowlabel <- function(subgroup, level) {
  ifelse(subgroup == "age_group",    paste0("Age: ", level, " years"),
  ifelse(subgroup == "sex",          level,
  ifelse(subgroup == "hiv_aids",     ifelse(level == "Negative", "HIV-negative", "HIV-positive"),
  ifelse(subgroup == "homelessness", ifelse(level == "Yes", "Experiencing homelessness", "Housed"),
  paste(subgroup, level, sep = " — ")))))
}
dfC$rowlabel <- make_rowlabel(dfC$Subgroup, dfC$Level)

dfC$Subgroup_clean <- factor(
  dplyr::recode(dfC$Subgroup,
    "age_group"    = "Age",
    "sex"          = "Sex",
    "hiv_aids"     = "HIV status",
    "homelessness" = "Homelessness"),
  levels = c("Age", "Sex", "HIV status", "Homelessness"))
dfC <- dfC |> dplyr::arrange(Subgroup_clean, Level) |>
  dplyr::mutate(rowlabel = factor(rowlabel, levels = rev(unique(rowlabel))))

dfC$hr_text <- sprintf("%.2f (%.2f–%.2f)", dfC$HR, dfC$CI_L, dfC$CI_H)

y_top <- length(levels(dfC$rowlabel)) + 0.8

pC_text <- ggplot(dfC, aes(y = rowlabel)) +
  geom_text(aes(x = 0, label = rowlabel, color = Subgroup_clean),
            hjust = 0, size = 3.8, fontface = "plain") +
  geom_text(aes(x = 1, label = hr_text),
            hjust = 1, size = 3.7, color = "grey25") +
  annotate("text", x = 0, y = y_top, label = "Characteristic",
           hjust = 0, fontface = "bold", size = 3.9, color = "grey20") +
  annotate("text", x = 1, y = y_top, label = "AHR (95% CI)",
           hjust = 1, fontface = "bold", size = 3.9, color = "grey20") +
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  labs(title = "C. Mortality AHR by subgroup") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13,
                                   margin = margin(b = 8)),
        plot.margin = margin(5, 2, 20, 8))

pC_forest <- ggplot(dfC, aes(x = HR, y = rowlabel, color = Subgroup_clean)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5,
             color = "grey50") +
  geom_errorbarh(aes(xmin = CI_L, xmax = CI_H),
                 height = 0.25, linewidth = 0.8) +
  geom_point(size = 3.4) +
  scale_x_log10(breaks = c(0.5, 1, 2, 3, 4),
                minor_breaks = c(0.75, 1.5, 2.5, 3.5),
                limits = c(0.7, 5)) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.5))) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "AHR", y = NULL, color = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        legend.margin = margin(t = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey88", linewidth = 0.35),
        panel.grid.minor.x = element_line(color = "grey94", linewidth = 0.25),
        plot.margin = margin(30, 8, 5, 0))

pC <- pC_text + pC_forest + plot_layout(widths = c(1.4, 1))

# ---------------------------------------------------------------------------
# Compose primary figure: A on top, B and C side-by-side beneath
# ---------------------------------------------------------------------------
fig3 <- pA / (pB | pC) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

ggsave(OUT_PNG, fig3, width = 14, height = 10, dpi = 300, bg = "white")
ggsave(OUT_PDF, fig3, width = 14, height = 10, bg = "white")
cat(sprintf("[fig3-defnB] Wrote %s\n", OUT_PNG))
cat(sprintf("[fig3-defnB] Wrote %s\n", OUT_PDF))

# ---------------------------------------------------------------------------
# Supplementary figure: unadjusted vs adjusted HR(t), side-by-side
# ---------------------------------------------------------------------------
cat("[fig3-defnB] Supplementary: unadjusted vs adjusted HR(t) ...\n")

pUn <- ggplot(hr_tbl, aes(x = month_mid, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_L, ymax = CI_H),
                width = 0.02, linewidth = 0.4, color = "grey50", alpha = 0.6) +
  geom_smooth(method = "loess", span = 0.55, se = TRUE,
              color = "#7f8c8d", fill = "#bdc3c7",
              linewidth = 1.2, alpha = 0.25) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                minor_breaks = c(0.33, 0.75, 1.5, 3, 6),
                limits = c(0.15, 10)) +
  scale_x_continuous(breaks = seq(0, HORIZON, 0.25),
                     labels = function(x) as.character(round(x * 12)),
                     limits = c(0, HORIZON)) +
  labs(title = "Unadjusted (per-month Cox)",
       x = "Time since treatment start (months)", y = "Crude HR") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35))

pAd <- ggplot(pooled, aes(x = month_mid, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_ribbon(aes(ymin = CI_L, ymax = CI_H),
              fill = "#e74c3c", alpha = 0.20) +
  geom_line(color = "#c0392b", linewidth = 1.4) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                minor_breaks = c(0.33, 0.75, 1.5, 3, 6),
                limits = c(0.15, 10)) +
  scale_x_continuous(breaks = seq(0, HORIZON * 12, 3),
                     limits = c(0, HORIZON * 12)) +
  labs(title = "Adjusted (MI-pooled spline tt(expose))",
       x = "Time since treatment start (months)", y = "Adjusted HR") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35))

supp_fig <- (pUn | pAd) +
  plot_annotation(
    title = "Time-varying HR for LTFU vs. on-treatment, anchored at treatment initiation",
    subtitle = "Both panels use the same defn-B counting-process specification; the right panel additionally adjusts for the target-trial covariate set and pools across imputations.",
    theme = theme(plot.background = element_rect(fill = "white", color = NA),
                  plot.title = element_text(face = "bold"))
  )

ggsave(SUPP_PNG, supp_fig, width = 12, height = 5, dpi = 300, bg = "white")
ggsave(SUPP_PDF, supp_fig, width = 12, height = 5, bg = "white")
cat(sprintf("[fig3-defnB] Wrote %s\n", SUPP_PNG))
cat(sprintf("[fig3-defnB] Wrote %s\n", SUPP_PDF))
