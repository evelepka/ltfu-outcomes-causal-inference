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
# Panel A: Time-varying hazard RATIO HR(t) via piecewise Cox
# ---------------------------------------------------------------------------
# Counting-process format: everyone enters as "not abandoned" at t=0; LTFU
# individuals transition to "abandoned" at their tx_duration. We then split
# follow-up at monthly boundaries and fit Cox with expose × month
# interaction, extracting HR per interval. The pattern should show HR < 1
# very early (immortal-time artefact — abandoners haven't crossed over yet
# at short times) then HR > 1 as abandoners accumulate risk.
# ---------------------------------------------------------------------------
cat("[fig3] Panel A: HR(t) piecewise Cox ...\n")
df <- read.csv(COHORT_CSV, stringsAsFactors = FALSE)
df$patient_id <- seq_len(nrow(df))
df$tx_dur <- as.numeric(difftime(as.Date(df$end_date),
                                 as.Date(df$best_start), units = "days")) / 365.25

HORIZON <- 2.0
df$event_d  <- ifelse(df$time_d_tx > HORIZON, 0, df$event_d)
df$time_d_tx <- ifelse(df$time_d_tx > HORIZON, HORIZON, df$time_d_tx)

# Everyone enters "not abandoned" at t=0; stays until abandonment or end.
df_p1 <- df
df_p1$tstart <- 0
df_p1$tstop  <- ifelse(df_p1$itt_group == "Loss to follow-up" &
                         df_p1$tx_dur < df_p1$time_d_tx,
                       df_p1$tx_dur, df_p1$time_d_tx)
df_p1$event <- ifelse(df_p1$itt_group == "Loss to follow-up" &
                        df_p1$tx_dur < df_p1$time_d_tx,
                      0, df_p1$event_d)
df_p1$expose <- 0

# LTFU individuals spawn a second row starting at abandonment (expose=1)
df_p2 <- df[df$itt_group == "Loss to follow-up" & df$tx_dur < df$time_d_tx, ]
df_p2$tstart <- df_p2$tx_dur
df_p2$tstop  <- df_p2$time_d_tx
df_p2$event  <- df_p2$event_d
df_p2$expose <- 1

df_split <- bind_rows(df_p1, df_p2) |>
  dplyr::filter(round(tstop - tstart, 4) > 0)

# Further split each interval at monthly boundaries (24 buckets over 2 yrs).
BREAKS <- seq(0, HORIZON, by = 1/12)
df_piecewise <- survSplit(Surv(tstart, tstop, event) ~ .,
                          data = df_split, cut = BREAKS[-c(1, length(BREAKS))],
                          episode = "month_bin")
# month_bin = 1 corresponds to (0, 1/12]; label by midpoint in months
df_piecewise$month_bin <- factor(df_piecewise$month_bin,
                                 levels = seq_len(length(BREAKS) - 1))

fit_pw <- coxph(Surv(tstart, tstop, event) ~ expose:strata(month_bin),
                data = df_piecewise, id = patient_id, control = coxph.control(timefix = FALSE))

# Fit per bucket separately so we get clean per-bucket HR + CI
bucket_hr <- function(d, bin) {
  sub <- d[d$month_bin == bin, ]
  if (sum(sub$event) < 3 || length(unique(sub$expose)) < 2) return(NULL)
  f <- tryCatch(coxph(Surv(tstart, tstop, event) ~ expose,
                      data = sub, id = patient_id),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  s <- summary(f)
  data.frame(
    month_mid = (as.numeric(as.character(bin)) - 0.5) / 12,
    HR = s$coefficients[1, "exp(coef)"],
    CI_L = s$conf.int[1, "lower .95"],
    CI_H = s$conf.int[1, "upper .95"],
    n_events = sum(sub$event)
  )
}

bins <- levels(df_piecewise$month_bin)
hr_tbl <- do.call(rbind, lapply(bins, function(b) bucket_hr(df_piecewise, b)))
hr_tbl <- hr_tbl[is.finite(hr_tbl$HR) & is.finite(hr_tbl$CI_H) & hr_tbl$CI_H < 50, ]
cat(sprintf("[fig3] Piecewise HR(t) estimated at %d time buckets\n", nrow(hr_tbl)))

write.csv(hr_tbl,
          file.path(ITT_RESULTS_DIR, "Figure_3a_HR_over_time.csv"),
          row.names = FALSE)

pA <- ggplot(hr_tbl, aes(x = month_mid, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  # LOESS smooth on the log10 scale (the scale_y_log10 transforms data first,
  # so geom_smooth fits on log-HR). Points are kept faint for transparency.
  geom_point(color = "#c0392b", size = 1.6, alpha = 0.35) +
  geom_smooth(method = "loess", span = 0.55, se = TRUE,
              color = "#c0392b", fill = "#e74c3c",
              linewidth = 1.4, alpha = 0.22) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                minor_breaks = c(0.33, 0.75, 1.5, 3, 6),
                limits = c(0.15, 10)) +
  scale_x_continuous(breaks = seq(0, HORIZON, 0.25),
                     labels = function(x) paste0(round(x*12), "mo"),
                     limits = c(0, HORIZON)) +
  labs(title = "A. Time-varying hazard ratio (crude, counting-process)",
       subtitle = "Piecewise monthly HR for abandoners vs. those still on treatment at that instant; LOESS smooth with 95% band.\nEarly months: apparent HR < 1 (immortal-time artefact). Later months: HR grows, reflecting accumulated abandonment mortality.",
       x = "Time since treatment start",
       y = "HR (log scale)") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11),
        panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
        panel.grid.minor.y = element_line(color = "grey94", linewidth = 0.25))

# ---------------------------------------------------------------------------
# Panel B: Target-trial HR by month of abandonment — Early vs Late (MI-pooled)
# ---------------------------------------------------------------------------
cat("[fig3] Panel B: early-vs-late HR array ...\n")
tt_csv <- file.path(ITT_RESULTS_DIR, "target_trial_mi_early_late_array.csv")
if (!file.exists(tt_csv)) stop(sprintf("Missing %s (run 30c_itt_target_trial_early_late.R)", tt_csv))
dfB_raw <- read.csv(tt_csv, stringsAsFactors = FALSE)
# Keep the Early (cap=0.5) and Late (cap=5) rows
dfB <- dfB_raw |>
  dplyr::filter((model == "early" & cap == 0.5) | (model == "late" & cap == 5)) |>
  dplyr::mutate(
    Month = as.numeric(gsub("Month_", "", Trial_Month)),
    Window = dplyr::recode(model,
                           "early" = "Early (0–6 months)",
                           "late"  = "Late (6–60 months)")
  )
dfB$Window <- factor(dfB$Window, levels = c("Early (0–6 months)", "Late (6–60 months)"))

palette_EL <- c("Early (0–6 months)" = "#3498db",
                "Late (6–60 months)" = "#e74c3c")

pB <- ggplot(dfB, aes(x = Month, y = HR, color = Window, group = Window)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_L, ymax = CI_H),
                width = 0.14, linewidth = 0.7) +
  geom_line(linewidth = 0.9, alpha = 0.8) +
  geom_point(size = 3.5) +
  scale_color_manual(values = palette_EL) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4),
                minor_breaks = c(0.33, 0.75, 1.5, 3)) +
  scale_x_continuous(breaks = 1:6, labels = paste("Mo", 1:6)) +
  labs(title = "B. Mortality HR by month of LTFU",
       subtitle = "Sequential target-trial emulation; MI-pooled; Early = deaths in first 6 months from treatment start; Late = deaths 6-60 months",
       x = "Month of loss to follow-up",
       y = "Hazard ratio (log scale)",
       color = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11),
        panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
        panel.grid.minor.y = element_line(color = "grey94", linewidth = 0.25))

# ---------------------------------------------------------------------------
# Panel C: Subgroup forest — LATE mortality (MI-pooled, from 32b, 5-yr cap)
# ---------------------------------------------------------------------------
cat("[fig3] Panel C: subgroup forest (late, 5y cap) ...\n")
sub_csv <- file.path(ITT_RESULTS_DIR, "target_trial_subgroup_interactions_mi.csv")
if (!file.exists(sub_csv)) stop(sprintf("Missing %s (run 32b_itt_target_trial_subgroups_mi.R)", sub_csv))
dfC_all <- read.csv(sub_csv, stringsAsFactors = FALSE)
dfC <- dfC_all |> dplyr::filter(model == "late", cap == 5)

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

# Two-panel forest: left panel = row labels + HR text (sans-serif);
# right panel = point + CI bars. patchwork shares the y-axis ordering
# via the common factor in dfC$rowlabel.
dfC$hr_text <- sprintf("%.2f (%.2f–%.2f)", dfC$HR, dfC$CI_L, dfC$CI_H)

# Add a horizontal "group header" row between subgroup categories
# (we'll approximate via colored row labels; ggplot doesn't do true headers).
pC_text <- ggplot(dfC, aes(y = rowlabel)) +
  # Row label on the far left, colored by subgroup
  geom_text(aes(x = 0, label = rowlabel, color = Subgroup_clean),
            hjust = 0, size = 3.6, fontface = "plain") +
  # HR text on the right of the text panel, dark grey
  geom_text(aes(x = 1, label = hr_text),
            hjust = 1, size = 3.5, color = "grey25") +
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  labs(title = "C. Late-mortality HR by subgroup",
       subtitle = "Sequential target-trial emulation; deaths 6–60 months from treatment start; MI-pooled") +
  theme_void(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13,
                                   margin = margin(b = 2)),
        plot.subtitle = element_text(size = 9, color = "grey40",
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
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Hazard ratio (log scale)", y = NULL, color = NULL) +
  theme_classic(base_size = 13) +
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
