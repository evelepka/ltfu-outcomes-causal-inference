# 52. Figure 3 — Causal mortality panels
# ==============================================================================
# Three-panel figure:
#   A. Crude time-varying mortality (Nelson-Aalen with abandonment as a
#      time-varying exposure). Shows the immortal-time-bias pattern where
#      abandoners appear falsely protected early and catastrophic late.
#   B. Sequential target-trial Cox HR by month of abandonment (MI-pooled),
#      showing late abandonment still carries a large causal hazard.
#   C. Forest plot of subgroup-stratified target-trial HRs (MI-pooled),
#      highlighting heterogeneity in the mortality penalty.
# Inputs:
#   - itt_cohort.csv (for panel A)
#   - target_trial_mi_6mo_array_hr.csv        (panel B, from script 30b)
#   - target_trial_subgroup_interactions_mi.csv (panel C, from script 32b)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
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

OUT_PNG <- file.path(ITT_RESULTS_DIR, "Figure_3_causal_mortality.png")
OUT_PDF <- file.path(ITT_RESULTS_DIR, "Figure_3_causal_mortality.pdf")

# ---------------------------------------------------------------------------
# Panel A: Crude time-varying mortality (counting-process Nelson-Aalen)
# ---------------------------------------------------------------------------
cat("[fig3] Panel A: crude time-varying cumulative hazard ...\n")
df <- read.csv(COHORT_CSV, stringsAsFactors = FALSE)
df$patient_id <- seq_len(nrow(df))
df$tx_dur <- as.numeric(difftime(as.Date(df$end_date),
                                 as.Date(df$best_start), units = "days")) / 365.25

HORIZON <- 2.0
df$event_d  <- ifelse(df$time_d_tx > HORIZON, 0, df$event_d)
df$time_d_tx <- ifelse(df$time_d_tx > HORIZON, HORIZON, df$time_d_tx)

# Everyone enters "Maintained Care" at time 0, stays until abandonment or end
df_p1 <- df
df_p1$tstart <- 0
df_p1$tstop  <- ifelse(df_p1$itt_group == "Loss to follow-up" &
                         df_p1$tx_dur < df_p1$time_d_tx,
                       df_p1$tx_dur, df_p1$time_d_tx)
df_p1$event <- ifelse(df_p1$itt_group == "Loss to follow-up" &
                        df_p1$tx_dur < df_p1$time_d_tx,
                      0, df_p1$event_d)
df_p1$group <- "Maintained care"

# LTFU individuals spawn a second row starting at abandonment
df_p2 <- df[df$itt_group == "Loss to follow-up" & df$tx_dur < df$time_d_tx, ]
df_p2$tstart <- df_p2$tx_dur
df_p2$tstop  <- df_p2$time_d_tx
df_p2$event  <- df_p2$event_d
df_p2$group  <- "Abandoned treatment"

df_split <- bind_rows(df_p1, df_p2) |>
  dplyr::filter(round(tstop - tstart, 4) > 0)

fit <- survfit(Surv(tstart, tstop, event) ~ group,
               data = df_split, id = patient_id)

df_plot_A <- data.frame(
  Time = fit$time,
  CumHaz = fit$cumhaz,
  strata = rep(names(fit$strata), fit$strata)
)
df_plot_A$strata <- gsub("group=", "", df_plot_A$strata)
df_plot_A <- bind_rows(
  data.frame(Time = 0, CumHaz = 0, strata = c("Maintained care", "Abandoned treatment")),
  df_plot_A
)

palette_A <- c("Maintained care" = "#2c3e50", "Abandoned treatment" = "#e74c3c")
pA <- ggplot(df_plot_A, aes(x = Time, y = CumHaz,
                            color = strata, linetype = strata)) +
  geom_step(linewidth = 1.15) +
  scale_color_manual(values = palette_A) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 0.14),
                     breaks = seq(0, 0.14, 0.02)) +
  scale_x_continuous(breaks = seq(0, HORIZON, 0.5)) +
  labs(title = "A. Crude time-varying mortality",
       subtitle = "Cumulative hazard; everyone starts in 'Maintained care' and moves to 'Abandoned treatment' on the day they drop out.\nAbandoners' curve only starts after abandonment — classic immortal-time-bias pattern.",
       x = "Years since treatment start",
       y = "Cumulative hazard of mortality",
       color = NULL, linetype = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = c(0.28, 0.88),
        legend.key.width = unit(1.4, "cm"),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9))

# ---------------------------------------------------------------------------
# Panel B: Target-trial HR by month of abandonment (MI-pooled)
# ---------------------------------------------------------------------------
cat("[fig3] Panel B: target trial HR array ...\n")
tt_csv <- file.path(ITT_RESULTS_DIR, "target_trial_mi_6mo_array_hr.csv")
if (!file.exists(tt_csv)) stop(sprintf("Missing %s (run 30b_itt_target_trial_mi_generator.R)", tt_csv))
dfB <- read.csv(tt_csv, stringsAsFactors = FALSE)
dfB$Month <- as.numeric(gsub("Month_", "", dfB$Trial_Month))
dfB$sig <- ifelse(dfB$CI_Lower > 1 | dfB$CI_Upper < 1, "significant", "not sig.")

pB <- ggplot(dfB, aes(x = Month, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper),
                width = 0.12, linewidth = 0.7, color = "#34495e") +
  geom_line(color = "#34495e", linewidth = 0.7, alpha = 0.6) +
  geom_point(aes(color = sig), size = 3.5) +
  scale_color_manual(values = c("significant" = "#e74c3c",
                                "not sig." = "#95a5a6")) +
  scale_y_continuous(breaks = seq(0.5, 2.5, 0.5)) +
  scale_x_continuous(breaks = 1:6, labels = paste("Mo", 1:6)) +
  labs(title = "B. Mortality HR by month of abandonment",
       subtitle = "Sequential target-trial emulation, 2-year horizon, MI-pooled (m=5)",
       x = "Month of loss to follow-up",
       y = "Hazard ratio (abandon vs. stay)",
       color = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9))

# ---------------------------------------------------------------------------
# Panel C: Subgroup forest (MI-pooled, from 32b)
# ---------------------------------------------------------------------------
cat("[fig3] Panel C: subgroup forest ...\n")
sub_csv <- file.path(ITT_RESULTS_DIR, "target_trial_subgroup_interactions_mi.csv")
if (!file.exists(sub_csv)) stop(sprintf("Missing %s (run 32b_itt_target_trial_subgroups_mi.R)", sub_csv))
dfC <- read.csv(sub_csv, stringsAsFactors = FALSE)

# Pretty labels
dfC <- dfC |>
  dplyr::mutate(
    Subgroup_clean = dplyr::recode(Subgroup,
      "age_group"    = "Age group",
      "sex"          = "Sex",
      "hiv_aids"     = "HIV status",
      "homelessness" = "Homelessness"
    ),
    Level_clean = dplyr::recode(Level,
      "Negative" = "HIV-negative", "Positive" = "HIV-positive",
      "No" = "Not homeless", "Yes" = "Homeless"
    ),
    rowlabel = sprintf("%s — %s", Subgroup_clean, Level_clean)
  )

# Order rows by subgroup category then level
dfC$Subgroup_clean <- factor(dfC$Subgroup_clean,
                             levels = c("Age group", "Sex", "HIV status", "Homelessness"))
dfC <- dfC |> dplyr::arrange(Subgroup_clean, Level) |>
  dplyr::mutate(rowlabel = factor(rowlabel, levels = rev(unique(rowlabel))))

pC <- ggplot(dfC, aes(x = HR, y = rowlabel, color = Subgroup_clean)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = CI_L, xmax = CI_H),
                 height = 0.15, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)", HR, CI_L, CI_H)),
            hjust = -0.15, size = 3.2, color = "black") +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4),
                limits = c(0.4, 8)) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "C. Mortality HR by subgroup (target trial)",
       subtitle = "Sequential target-trial emulation stratified by baseline subgroup; MI-pooled",
       x = "Hazard ratio (log scale)", y = NULL, color = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9),
        axis.text.y = element_text(size = 10))

# ---------------------------------------------------------------------------
# Compose: A on top, B and C side-by-side beneath
# ---------------------------------------------------------------------------
fig3 <- pA / (pB | pC) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

ggsave(OUT_PNG, fig3, width = 14, height = 10, dpi = 300, bg = "white")
ggsave(OUT_PDF, fig3, width = 14, height = 10, bg = "white")
cat(sprintf("[fig3] Wrote %s\n", OUT_PNG))
cat(sprintf("[fig3] Wrote %s\n", OUT_PDF))
