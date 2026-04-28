# 55. Appendix figure — defn-A vs defn-B side-by-side comparison
# ==============================================================================
# Two-panel figure for the appendix:
#   A. Late-mortality aHR (cap = 2 yr) by trial month, defn-A grace vs
#      defn-B grace lines overlaid.
#   B. Cause-specific TB aHR by trial month, defn-A vs defn-B
#      (hybrid attribution, late mortality cap=2 yr).
#
# Output: Figure_S_defn_comparison.png/pdf
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(scales)
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

OUT_PNG <- file.path(ITT_RESULTS_DIR, "Figure_S_defn_comparison.png")
OUT_PDF <- file.path(ITT_RESULTS_DIR, "Figure_S_defn_comparison.pdf")

# ---- Panel A: late mortality
defB <- read.csv(file.path(ITT_RESULTS_DIR, "target_trial_defnB_mi_early_late_array.csv"),
                 stringsAsFactors = FALSE)
defA <- read.csv(file.path(ITT_RESULTS_DIR, "target_trial_grace_mi_early_late_array.csv"),
                 stringsAsFactors = FALSE)
defB$Month <- as.numeric(gsub("Month_", "", defB$Trial_Month))
defA$Month <- as.numeric(gsub("Month_", "", defA$Trial_Month))
defB$Defn <- "Defn B (last visit; primary)"
defA$Defn <- "Defn A (recorded end_date)"

dfA <- bind_rows(
  defB |> filter(model == "late", cap == 2),
  defA |> filter(model == "late", cap == 2)
)
dfA$Defn <- factor(dfA$Defn,
                   levels = c("Defn B (last visit; primary)",
                              "Defn A (recorded end_date)"))

palette_defn <- c("Defn B (last visit; primary)" = "#2c3e50",
                  "Defn A (recorded end_date)"   = "#e67e22")

pA <- ggplot(dfA, aes(x = Month, y = HR, color = Defn, group = Defn)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_L, ymax = CI_H),
                width = 0.14, linewidth = 0.7,
                position = position_dodge(width = 0.25)) +
  geom_line(linewidth = 0.9, alpha = 0.8,
            position = position_dodge(width = 0.25)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.25)) +
  scale_color_manual(values = palette_defn) +
  scale_y_log10(breaks = c(0.5, 1, 2, 3, 5),
                minor_breaks = c(0.75, 1.5, 2.5, 4),
                limits = c(0.6, 5)) +
  scale_x_continuous(breaks = 1:6, labels = paste("Mo", 1:6)) +
  labs(title = "A. Late-mortality aHR by trial month — definition comparison",
       subtitle = "Sequential target-trial emulation, MI-pooled, grace-period eligibility; late = 6–24 mo from grace-shifted origin",
       x = "Trial month",
       y = "Hazard ratio (log scale)",
       color = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9),
        panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
        panel.grid.minor.y = element_line(color = "grey94", linewidth = 0.25))

# ---- Panel B: cause-specific TB-cause aHR comparison
csB <- read.csv(file.path(ITT_RESULTS_DIR, "target_trial_defnB_cause_specific.csv"),
                stringsAsFactors = FALSE)
csA <- read.csv(file.path(ITT_RESULTS_DIR, "target_trial_grace_cause_specific.csv"),
                stringsAsFactors = FALSE)
csB$Month <- as.numeric(gsub("Month_", "", csB$Trial_Month))
csA$Month <- as.numeric(gsub("Month_", "", csA$Trial_Month))
csB$Defn <- "Defn B (last visit; primary)"
csA$Defn <- "Defn A (recorded end_date)"

dfB_panel <- bind_rows(
  csB |> filter(cause == "tb_hybrid", cap == 2),
  csA |> filter(cause == "tb_hybrid", cap == 2)
)
dfB_panel$Defn <- factor(dfB_panel$Defn,
                          levels = c("Defn B (last visit; primary)",
                                     "Defn A (recorded end_date)"))

pB <- ggplot(dfB_panel, aes(x = Month, y = HR, color = Defn, group = Defn)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_L, ymax = CI_H),
                width = 0.14, linewidth = 0.7,
                position = position_dodge(width = 0.25)) +
  geom_line(linewidth = 0.9, alpha = 0.8,
            position = position_dodge(width = 0.25)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.25)) +
  scale_color_manual(values = palette_defn) +
  scale_y_log10(breaks = c(0.5, 1, 2, 3, 5, 10),
                minor_breaks = c(0.75, 1.5, 2.5, 4, 7),
                limits = c(0.5, 12)) +
  scale_x_continuous(breaks = 1:6, labels = paste("Mo", 1:6)) +
  labs(title = "B. Cause-specific TB-cause aHR by trial month — definition comparison",
       subtitle = "Hybrid attribution (SIM ICD-10 + TBweb Obito TB); late mortality 6–24 mo",
       x = "Trial month",
       y = "TB-cause hazard ratio (log scale)",
       color = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9),
        panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
        panel.grid.minor.y = element_line(color = "grey94", linewidth = 0.25))

fig <- pA / pB +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Appendix Figure A1. Definition-A vs definition-B comparison of late-mortality aHRs",
    subtitle = "The trial-month axis is interpreted differently between the two definitions: Defn-B Month m corresponds approximately to Defn-A Month m+1 (i.e., a 30-day shift). Hence both lines convey the same underlying biological signal, just differently labeled.",
    theme = theme(plot.background = element_rect(fill = "white", color = NA),
                  plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 10, color = "grey25",
                                               margin = margin(b = 6)))
  )

ggsave(OUT_PNG, fig, width = 11, height = 9, dpi = 300, bg = "white")
ggsave(OUT_PDF, fig, width = 11, height = 9, bg = "white")
cat(sprintf("[fig-appendix] Wrote %s\n", OUT_PNG))
cat(sprintf("[fig-appendix] Wrote %s\n", OUT_PDF))
