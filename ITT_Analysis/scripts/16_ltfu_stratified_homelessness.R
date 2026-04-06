# 16_ltfu_stratified_homelessness.R
# Stratified mortality analysis (Homelessness Yes vs No) within the LTFU group only

library(dplyr)
library(survival)
library(ggplot2)
library(tidyr)
library(scales)

# 1. Load Data
cohort <- read.csv("Abandonment Paper/ITT_Analysis/data/itt_cohort.csv")

# Filter for LTFU only
ltfu_data <- cohort %>% filter(itt_group == "Loss to follow-up")

# 2. Survival Analysis - Homelessness
fit_home <- survfit(Surv(time_d_tx, event_d) ~ homelessness, data = ltfu_data)

# Log-rank test
sdiff <- survdiff(Surv(time_d_tx, event_d) ~ homelessness, data = ltfu_data)
pval <- 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
pval_text <- if(pval < 0.001) "p < 0.001" else paste0("p = ", round(pval, 3))

# 3. Extract Values at 1, 5, 10, and 11.9 Years
get_vals <- function(f, t) {
  s <- summary(f, times = t)
  data.frame(strata = s$strata, year = s$time, mortality_pct = (1 - s$surv) * 100)
}
results_table <- get_vals(fit_home, c(1, 5, 10, 11.9))

print("--- Cumulative Mortality Incidence (%) in LTFU Group by Homelessness ---")
print(results_table)

# 4. Plotting
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

plot_df <- tidy_surv(fit_home)
plot_df$mortality <- 1 - plot_df$surv
plot_df$homelessness_status <- gsub("homelessness=", "", plot_df$strata)

p <- ggplot(plot_df, aes(x = time, y = mortality, color = homelessness_status)) +
  geom_step(linewidth = 1.2) +
  theme_classic(base_size = 14) +
  labs(
    title = "Cumulative incidence of mortality by homelessness status",
    subtitle = "Kaplan-Meier estimates",
    x = "Years since loss to follow-up",
    y = "Cumulative incidence"
  ) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0, 0.40)) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 12)) +
  scale_color_manual(name = NULL, 
                     values = c("No" = "#377eb8", "Yes" = "#e41a1c"),
                     labels = c("No" = "Not homeless", "Yes" = "Homeless")) +
  theme(plot.title = element_text(face="bold"), legend.position = "bottom")

# Save Results
out_dir <- "Abandonment Paper/ITT_Analysis/results/manuscript figures/stratified_mortality/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(paste0(out_dir, "Mortality_KM_12y_by_Homelessness.png"), p, width = 8, height = 6, dpi = 300)
write.csv(results_table, paste0(out_dir, "ltfu_mortality_homelessness_values.csv"), row.names = FALSE)

cat(paste0("PVAL:", pval_text))
print("Script completed. Plot and values saved.")
