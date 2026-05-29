# 54. Figure 4 — Cause-specific mortality (TB vs non-TB) target-trial
# ==============================================================================
# Two-panel main-manuscript figure visualising the central causal-interpretation
# evidence: LTFU exposure raises TB-cause mortality strongly, while the
# negative-control non-TB cause hazards remain near unity (or only modestly
# elevated). Inputs:
#   - target_trial_defnB_cause_specific.csv  (defn-B + grace; primary)
#   - target_trial_grace_cause_specific.csv  (defn-A + grace; sensitivity overlay)
#
# Panel A. Hybrid attribution (SIM ICD-10 + TBweb Obito TB/NTB) — primary.
# Panel B. SIM-only attribution (uniform method across arms) — sensitivity.
#
# In each panel: TB-cause aHR (red) vs non-TB aHR (blue) by trial month,
# both with 95% CIs. Reference line at 1.
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

OUT_PNG <- file.path(ITT_RESULTS_DIR, "Figure_4_cause_specific.png")
OUT_PDF <- file.path(ITT_RESULTS_DIR, "Figure_4_cause_specific.pdf")
# Appendix sensitivity (SIM-only) goes out separately so the main figure is panel-A only
OUT_APP_PNG <- file.path(ITT_RESULTS_DIR, "Figure_S_cause_specific_simonly.png")
OUT_APP_PDF <- file.path(ITT_RESULTS_DIR, "Figure_S_cause_specific_simonly.pdf")

cs_path <- file.path(ITT_RESULTS_DIR, "target_trial_defnB_cause_specific.csv")
if (!file.exists(cs_path)) stop(sprintf("Missing %s (run 30i)", cs_path))
cs <- read.csv(cs_path, stringsAsFactors = FALSE)
cs$Month <- as.numeric(gsub("Month_", "", cs$Trial_Month))
cs <- cs |> dplyr::filter(cap == 2)

mk_panel <- function(df_panel, title, subtitle, c_tb, c_ntb) {
  d_tb  <- df_panel |> dplyr::filter(cause == c_tb)  |> dplyr::mutate(Cause = "TB-cause death")
  d_ntb <- df_panel |> dplyr::filter(cause == c_ntb) |> dplyr::mutate(Cause = "Non-TB cause death")
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
         x = "Month of disengagement (defn-B)",
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
               "Hybrid attribution (SIM ICD-10 + TBweb Obito) primary",
               "TB-cause aHR (red) vs non-TB aHR (blue); late mortality 6 to 24 months from disengagement",
               "tb_hybrid", "nontb_hybrid")

pB <- mk_panel(cs,
               "SIM-only attribution (uniform across arms) sensitivity",
               "Restricted to deaths with SIM ICD-10 codes; same trials, same covariates",
               "tb_simonly", "nontb_simonly")

# Main manuscript Figure 5: Panel A (hybrid attribution) only
fig_main <- pA +
  plot_annotation(
    title = "Figure 5. Cause-specific mortality after disengagement",
    subtitle = "Sequential target-trial emulation, MI-pooled, with defn-B + grace eligibility. The TB-cause hazard ratio (red) rising sharply across trial months while the non-TB negative-control hazard (blue) remains near unity supports a causal interpretation: the LTFU effect is mediated through interrupted TB therapy rather than purely confounded by social/clinical predictors of disengagement.",
    theme = theme(plot.background = element_rect(fill = "white", color = NA),
                  plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9.5, color = "grey25",
                                               margin = margin(b = 8)))
  )

ggsave(OUT_PNG, fig_main, width = 11, height = 6, dpi = 300, bg = "white")
ggsave(OUT_PDF, fig_main, width = 11, height = 6, bg = "white")
cat(sprintf("[fig5-main] Wrote %s\n", OUT_PNG))
cat(sprintf("[fig5-main] Wrote %s\n", OUT_PDF))

# Appendix figure: Panel B (SIM-only attribution sensitivity)
fig_app <- pB +
  plot_annotation(
    title = "Figure S. SIM-only attribution sensitivity for cause-specific mortality",
    subtitle = "Same target-trial emulation as Figure 5, restricted to deaths with available SIM ICD-10 codes",
    theme = theme(plot.background = element_rect(fill = "white", color = NA),
                  plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9.5, color = "grey25",
                                               margin = margin(b = 8)))
  )

ggsave(OUT_APP_PNG, fig_app, width = 11, height = 6, dpi = 300, bg = "white")
ggsave(OUT_APP_PDF, fig_app, width = 11, height = 6, bg = "white")
cat(sprintf("[fig5-app]  Wrote %s\n", OUT_APP_PNG))
cat(sprintf("[fig5-app]  Wrote %s\n", OUT_APP_PDF))
