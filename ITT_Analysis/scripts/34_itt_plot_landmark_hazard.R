# 34. Generate Landmark Cumulative Hazard Plot
# ==============================================================================
# Shows the strictly diverging mortality hazards of Early Abandonment vs 
# those remaining in care, landmarked at Month 2 to strictly eliminate 
# immortal time bias from the baseline unadjusted curves.
# ==============================================================================
library(dplyr)
library(ggplot2)
library(survival)

cat("\n--- 1. Loading and Landmarking Cohort ---\n")
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors=FALSE)

# Landmark at 2 Months (60 days) to eliminate early immortal time bias
lm_days <- 60
lm_yrs <- lm_days/365.25

df$date_start <- as.Date(df$best_start)
df$date_end <- as.Date(df$end_date)
df$tx_duration_yrs <- as.numeric(difftime(df$date_end, df$date_start, units="days")) / 365.25

# Only include patients who survived to the landmark
df_lm <- df %>% filter(time_d > lm_yrs)

# Define exposure exactly at landmark: Did you abandon before Month 2?
df_lm <- df_lm %>%
  mutate(
    expose = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs <= lm_yrs, 
                   "Early Abandonment (Months 1-2)", 
                   "Persistent Care (Active or Cured)"),
    time_followup = time_d - lm_yrs,
    # Cap follow-up at 5 years
    event_d = ifelse(time_followup > 5, 0, event_d),
    time_followup = ifelse(time_followup > 5, 5, time_followup)
  )

cat("\n--- 2. Fitting Nelson-Aalen Target Curve ---\n")
fit <- survfit(Surv(time_followup, event_d) ~ expose, data=df_lm)

# Extract exactly at events for step plotting
df_p <- data.frame(
  time = fit$time,
  cumhaz = fit$cumhaz,
  strata = rep(names(fit$strata), fit$strata)
)
df_p$strata <- gsub("expose=", "", df_p$strata)

cat("\n--- 3. Plotting ---\n")
p1 <- ggplot(df_p, aes(x = time, y = cumhaz, color = strata, linetype = strata)) +
  geom_step(linewidth = 1.5) +
  scale_color_manual(values = c("Early Abandonment (Months 1-2)" = "#e74c3c", 
                               "Persistent Care (Active or Cured)" = "#2c3e50")) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0, 0.15, 0.05), limits=c(0, 0.15)) +
  scale_x_continuous(breaks = 0:5) +
  labs(
    title = "Diverging Mortality: Landmark Cumulative Hazard",
    subtitle = "Nelson-Aalen estimate from 2-Month Survival Landmark",
    x = "Years Since Landmark (Month 2)",
    y = "Cumulative Hazard of Mortality",
    color = "Cohort Trajectory", linetype = "Cohort Trajectory"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face="bold", size=18),
    legend.position = c(0.25, 0.85),
    legend.background = element_rect(fill="white", color="white"),
    legend.key.width = unit(2, "cm")
  )

ts <- as.numeric(Sys.time())
fname <- sprintf("figure_landmark_hazard_%.0f.png", ts)
artifact_dir <- "/Users/jasonandrews/.gemini/antigravity/brain/c053ef30-5842-41b7-b342-bf735650d865"

ggsave(file.path("ITT_Analysis/results", fname), plot=p1, width=10, height=7, dpi=300, bg="white")
ggsave(file.path(artifact_dir, fname), plot=p1, width=10, height=7, dpi=300, bg="white")

cat("\n--- 4. Updating Walkthrough Cache ---\n")
wt_path <- file.path(artifact_dir, "walkthrough.md")
if(file.exists(wt_path)){
    writeLines(paste0(paste(readLines(wt_path, warn=F), collapse="\n"), 
               "\n![Alternative Figure: Cumulative Hazard](", file.path(artifact_dir, fname), ")\n"), wt_path)
}
cat("Success! Generated", fname, "\n")
