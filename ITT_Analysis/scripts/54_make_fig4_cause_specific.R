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
    scale_x_continuous(breaks = 1:6, labels = as.character(1:6)) +
    labs(title = title, subtitle = subtitle,
         x = "Month of LTFU",
         y = "AHR (log scale)",
         color = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9),
          panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
          panel.grid.minor.y = element_line(color = "grey94", linewidth = 0.25))
}

pA <- mk_panel(cs,
               "A. Adjusted cause-specific mortality — hybrid attribution (primary)",
               "TB-cause AHR (red) vs non-TB AHR (blue); late mortality 6–24 mo from disengagement; SIM ICD-10 + TBweb Obito",
               "tb_hybrid", "nontb_hybrid")

pB <- mk_panel(cs,
               "B. Adjusted cause-specific mortality — SIM-only attribution (sensitivity)",
               "Restricted to deaths with SIM ICD-10 codes (uniform across arms); same trials, same covariates",
               "tb_simonly", "nontb_simonly")

fig4 <- pA / pB +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Figure 4. Adjusted cause-specific mortality after LTFU",
    subtitle = "Sequential target-trial emulation, MI-pooled; defn-B + grace eligibility. The TB-cause AHR (red) rising sharply across trial months while the non-TB negative-control AHR (blue) remains near unity supports a causal interpretation: the LTFU effect is mediated through interrupted TB therapy rather than purely confounded by social/clinical predictors of disengagement.",
    theme = theme(plot.background = element_rect(fill = "white", color = NA),
                  plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 9.5, color = "grey25",
                                               margin = margin(b = 8)))
  )

ggsave(OUT_PNG, fig4, width = 11, height = 9, dpi = 300, bg = "white")
ggsave(OUT_PDF, fig4, width = 11, height = 9, bg = "white")
cat(sprintf("[fig4] Wrote %s\n", OUT_PNG))
cat(sprintf("[fig4] Wrote %s\n", OUT_PDF))
