# 54b. Appendix Figure A2 — cause-specific (defn-A grace, sensitivity)
# ==============================================================================
# Mirrors 54 but uses defn-A grace-period cause-specific data.
# Output: Figure_S_defnA_cause_specific.png/pdf
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

OUT_PNG <- file.path(ITT_RESULTS_DIR, "Figure_S_defnA_cause_specific.png")
OUT_PDF <- file.path(ITT_RESULTS_DIR, "Figure_S_defnA_cause_specific.pdf")

cs <- read.csv(file.path(ITT_RESULTS_DIR, "target_trial_grace_cause_specific.csv"),
               stringsAsFactors = FALSE)
cs$Month <- as.numeric(gsub("Month_", "", cs$Trial_Month))
cs <- cs |> filter(cap == 2)

mk_panel <- function(df_panel, title, subtitle, c_tb, c_ntb) {
  d_tb  <- df_panel |> filter(cause == c_tb)  |> mutate(Cause = "TB-cause death")
  d_ntb <- df_panel |> filter(cause == c_ntb) |> mutate(Cause = "Non-TB cause death")
  d <- bind_rows(d_tb, d_ntb)
  d$Cause <- factor(d$Cause, levels = c("TB-cause death", "Non-TB cause death"))
  ggplot(d, aes(x = Month, y = HR, color = Cause, group = Cause)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(aes(ymin = CI_L, ymax = CI_H),
                  width = 0.14, linewidth = 0.7,
                  position = position_dodge(width = 0.25)) +
    geom_line(linewidth = 0.9, alpha = 0.7,
              position = position_dodge(width = 0.25)) +
    geom_point(size = 3.5, position = position_dodge(width = 0.25)) +
    scale_color_manual(values = c("TB-cause death" = "#c0392b",
                                   "Non-TB cause death" = "#2980b9")) +
    scale_y_log10(breaks = c(0.5, 1, 2, 3, 5, 10),
                  minor_breaks = c(0.75, 1.5, 2.5, 4, 7),
                  limits = c(0.5, 12)) +
    scale_x_continuous(breaks = 1:6, labels = paste("Mo", 1:6)) +
    labs(title = title, subtitle = subtitle,
         x = "Trial month (defn-A: recorded end_date)",
         y = "Hazard ratio (log scale)",
         color = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9),
          panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
          panel.grid.minor.y = element_line(color = "grey94", linewidth = 0.25))
}

pA <- mk_panel(cs,
               "A. Hybrid attribution (defn-A grace; sensitivity)",
               "TB-cause aHR (red) vs non-TB aHR (blue); late mortality 6–24 mo",
               "tb_hybrid", "nontb_hybrid")

pB <- mk_panel(cs,
               "B. SIM-only attribution (defn-A grace; sensitivity)",
               "Restricted to deaths with SIM ICD-10 codes",
               "tb_simonly", "nontb_simonly")

fig <- pA / pB +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Appendix Figure A2. Cause-specific mortality under definition A (sensitivity)",
    subtitle = "Same analysis as main Figure 4 but using the recorded end_date as the disengagement date (defn A) instead of end_date − 30 d (defn B). Trial-month labels shifted by one month relative to the primary; magnitudes broadly comparable.",
    theme = theme(plot.background = element_rect(fill = "white", color = NA),
                  plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9.5, color = "grey25",
                                               margin = margin(b = 8)))
  )

ggsave(OUT_PNG, fig, width = 11, height = 9, dpi = 300, bg = "white")
ggsave(OUT_PDF, fig, width = 11, height = 9, bg = "white")
cat(sprintf("[fig-appendix-A2] Wrote %s\n", OUT_PNG))
