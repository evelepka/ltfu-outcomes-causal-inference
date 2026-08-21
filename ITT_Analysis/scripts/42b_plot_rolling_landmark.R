# 42b. Figure: continuous timing curve from the rolling landmark,
#      with the monthly landmark estimates overlaid for comparison.
# ==============================================================================
# The monthly design can only place six points, at the midpoints of its
# disengagement windows. The rolling design estimates a continuous function of
# disengagement day, which is what the "is there a safe point to disengage?"
# question actually asks.
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  getwd()
}
source(file.path(.here(), "_paths.R"))

tim <- read.csv(file.path(ITT_RESULTS_DIR, "rolling_landmark_timing.csv")) |>
  filter(comparator == "in_care")

# monthly landmark, late window, cap = 2 -> plotted at window midpoints
mon <- read.csv(file.path(ITT_RESULTS_DIR, "target_trial_defnB_mi_early_late_array.csv")) |>
  filter(model == "late", cap == 2) |>
  mutate(m = as.integer(sub("Month_", "", Trial_Month)),
         day = (m - 0.5) * 30) |>
  filter(m <= 6)

p <- ggplot(tim, aes(day, HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  geom_ribbon(aes(ymin = CI_L, ymax = CI_H), fill = "#2C7FB8", alpha = 0.20) +
  geom_line(colour = "#2C7FB8", linewidth = 1.1) +
  geom_pointrange(data = mon, aes(day, HR, ymin = CI_L, ymax = CI_H),
                  colour = "#D95F02", size = 0.45, linewidth = 0.7) +
  scale_y_log10(breaks = c(0.5, 1, 1.5, 2, 3, 4, 5)) +
  scale_x_continuous(breaks = seq(0, 180, 30),
                     sec.axis = sec_axis(~ . / 30, name = "Month of therapy",
                                         breaks = 0:6)) +
  labs(x = "Day of last observed contact (inferred disengagement)",
       y = "Adjusted hazard ratio for late (6–24 mo) mortality (log scale)",
       subtitle = paste0("Blue: rolling landmark, continuous in day of disengagement ",
                         "(origin = declaration date)\n",
                         "Orange: monthly landmark, plotted at window midpoints")) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 10, colour = "grey25"))

for (ext in c("png", "pdf")) {
  f <- file.path(ITT_RESULTS_DIR, paste0("Figure_rolling_landmark_timing.", ext))
  ggsave(f, p, width = 9, height = 5.6, dpi = 300, bg = "white")
  cat(sprintf("[42b] wrote %s\n", f))
}
