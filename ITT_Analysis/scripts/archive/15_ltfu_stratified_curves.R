# 15_ltfu_stratified_curves.R
# Stratified mortality analysis (HIV+ vs HIV-) within the LTFU group only

library(dplyr)
library(survival)
library(ggplot2)
library(tidyr)
library(scales)

# 1. Load Data
cohort <- read.csv("Abandonment Paper/ITT_Analysis/data/itt_cohort.csv")

# Filter for LTFU only (itt_group == "Loss to follow-up")
ltfu_data <- cohort %>% filter(itt_group == "Loss to follow-up")

# 2. Prepare for Survival Analysis
# 'event_d' is mortality (1=death, 0=censored)
# 'time_d_tx' is time in years from treatment start
fit <- survfit(Surv(time_d_tx, event_d) ~ hiv_aids, data = ltfu_data)

# 3. Extract Values at 1, 5, 10, and 12 Years
summary_years <- summary(fit, times = c(1, 5, 10, 12))
results_table <- data.frame(
  strata = summary_years$strata,
  year = summary_years$time,
  mortality_pct = (1 - summary_years$surv) * 100,
  lower_ci = (1 - summary_years$upper) * 100,
  upper_ci = (1 - summary_years$lower) * 100
)

print("--- Cumulative Mortality Incidence (%) in LTFU Group by HIV ---")
print(results_table)

# 4. Plotting
# For better visualization, we'll use ggplot2
tidy_surv <- function(s_fit) {
  d <- data.frame(time = s_fit$time, surv = s_fit$surv, lower = s_fit$lower, upper = s_fit$upper, strata = rep(names(s_fit$strata), s_fit$strata))
  
  split_d <- split(d, d$strata)
  d_ext <- do.call(rbind, lapply(split_d, function(df_stratum) {
    last_row <- tail(df_stratum, 1)
    last_row$time <- 12.0
    rbind(
      data.frame(time=0, surv=1, lower=1, upper=1, strata=unique(df_stratum$strata)),
      df_stratum,
      last_row
    )
  }))
  rownames(d_ext) <- NULL
  d_ext
}

plot_df <- tidy_surv(fit)
plot_df$mortality <- 1 - plot_df$surv
plot_df$hiv_status <- gsub("hiv_aids=", "", plot_df$strata)

p <- ggplot(plot_df, aes(x = time, y = mortality, color = hiv_status)) +
  geom_step(linewidth = 1.2) +
  theme_classic(base_size = 14) +
  labs(
    title = "Cumulative incidence of mortality by HIV status",
    subtitle = "Kaplan-Meier estimates",
    x = "Years since loss to follow-up",
    y = "Cumulative incidence"
  ) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0, 0.40)) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 12)) +
  scale_color_manual(name = NULL, 
                     values = c("Negative" = "#377eb8", "Positive" = "#e41a1c", "Unknown" = "#999999"),
                     labels = c("Negative" = "HIV-negative", "Positive" = "HIV-positive", "Unknown" = "Unknown")) +
  theme(plot.title = element_text(face="bold"), legend.position = "bottom")

# Save Result
out_dir <- "Abandonment Paper/ITT_Analysis/results/manuscript figures/stratified_mortality/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(paste0(out_dir, "Mortality_KM_12y_by_HIV.png"), p, width = 8, height = 6, dpi = 300)
write.csv(results_table, paste0(out_dir, "ltfu_mortality_hiv_values.csv"), row.names = FALSE)

print("Script completed. Plot and values saved.")
