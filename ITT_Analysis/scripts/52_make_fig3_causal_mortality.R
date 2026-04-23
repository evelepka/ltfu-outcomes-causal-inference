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
  geom_smooth(method = "loess", span = 0.55, se = TRUE,
              color = "#c0392b", fill = "#e74c3c",
              linewidth = 1.4, alpha = 0.22) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                minor_breaks = c(0.33, 0.75, 1.5, 3, 6),
                limits = c(0.15, 10)) +
  scale_x_continuous(breaks = seq(0, HORIZON, 0.25),
                     labels = function(x) as.character(round(x * 12)),
                     limits = c(0, HORIZON)) +
  labs(title = "A. Time-varying hazard ratio",
       x = "Time since treatment start (months)",
       y = "HR") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13),
        panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
        panel.grid.minor.y = element_line(color = "grey94", linewidth = 0.25))

# ---------------------------------------------------------------------------
# Panel B: Target-trial HR by month of abandonment — Early vs Late (MI-pooled)
# ---------------------------------------------------------------------------
cat("[fig3] Panel B: early-vs-late HR array (24mo cap) ...\n")
tt_csv <- file.path(ITT_RESULTS_DIR, "target_trial_mi_early_late_array.csv")
if (!file.exists(tt_csv)) stop(sprintf("Missing %s (run 30c_itt_target_trial_early_late.R)", tt_csv))
dfB_raw <- read.csv(tt_csv, stringsAsFactors = FALSE)
# Early (cap=0.5) vs Late (cap=2 = 24 months)
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
# Panel C: Subgroup forest — LATE mortality (MI-pooled, from 32b, 5-yr cap)
# ---------------------------------------------------------------------------
cat("[fig3] Panel C: subgroup forest (late, 24mo cap) ...\n")
sub_csv <- file.path(ITT_RESULTS_DIR, "target_trial_subgroup_interactions_mi.csv")
if (!file.exists(sub_csv)) stop(sprintf("Missing %s (run 32b_itt_target_trial_subgroups_mi.R)", sub_csv))
dfC_all <- read.csv(sub_csv, stringsAsFactors = FALSE)
dfC <- dfC_all |> dplyr::filter(model == "late", cap == 2)

# Row labels: subgroup+level → single human-readable label
make_rowlabel <- function(subgroup, level) {
  ifelse(subgroup == "age_group",    paste0("Age: ", level, " years"),
  ifelse(subgroup == "sex",          level,
  ifelse(subgroup == "hiv_aids",     ifelse(level == "Negative", "HIV-negative", "HIV-positive"),
  ifelse(subgroup == "homelessness", ifelse(level == "Yes", "Experiencing homelessness", "Housed"),
  paste(subgroup, level, sep = " — ")))))
}
dfC$rowlabel <- make_rowlabel(dfC$Subgroup, dfC$Level)

# Ordering: subgroup category, then within-category level order
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

# Header row — use the top y slot for column headings
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
