# 44b. Figure: timing curve stratified by DOT status.
# The average effect is similar in both groups; the SHAPE is resolvable only under
# DOT, where the inferred disengagement day is close to the last observed dose.
# For self-administered patients monthly dispensing blurs the timing variable by
# up to a month, attenuating the curve toward flatness.
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
td <- read.csv(file.path(ITT_RESULTS_DIR, "rolling_landmark_timing_by_dot.csv")) |>
  mutate(dot = factor(ifelse(dot == "Yes", "Directly observed therapy",
                             "Self-administered (no DOT)"),
                      levels = c("Directly observed therapy",
                                 "Self-administered (no DOT)")))
p <- ggplot(td, aes(day, HR, colour = dot, fill = dot)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  geom_ribbon(aes(ymin = CI_L, ymax = CI_H), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("#1B7837", "#762A83"), name = NULL) +
  scale_fill_manual(values = c("#1B7837", "#762A83"), name = NULL) +
  scale_y_log10(breaks = c(0.5, 1, 1.5, 2, 3, 4, 6)) +
  scale_x_continuous(breaks = seq(0, 180, 30),
                     sec.axis = sec_axis(~ . / 30, name = "Month of therapy",
                                         breaks = 0:6)) +
  labs(x = "Day of last observed contact (inferred disengagement)",
       y = "Adjusted hazard ratio, late (6-24 mo) mortality (log scale)",
       subtitle = paste0("Average effect is similar in both groups; only the DOT ",
                         "curve resolves the SHAPE.\nUnder self-administration, ",
                         "monthly dispensing blurs the timing variable.")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 10, colour = "grey25"))
for (ext in c("png", "pdf")) {
  f <- file.path(ITT_RESULTS_DIR, paste0("Figure_rolling_timing_by_dot.", ext))
  ggsave(f, p, width = 9, height = 5.8, dpi = 300, bg = "white")
  cat(sprintf("[44b] wrote %s\n", basename(f)))
}
