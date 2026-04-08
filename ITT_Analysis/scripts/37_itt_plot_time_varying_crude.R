# 37. Time-Varying Crude Observational Dynamics
# ==============================================================================
# Calculates crude mortality by treating "Abandonment" as a time-varying 
# exposure (counting process format tstart -> tstop).
# Everyone contributes to "Maintained Care" until the exact day they drop out,
# meaning Maintained Care fully captures early severe/HIV deaths, verifying
# why abandoners retrospectively appear paradoxically "lower risk" overall.
# ==============================================================================
library(dplyr)
library(survival)
library(ggplot2)

cat("\n--- 1. Loading Cohort ---\n")
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# Give everyone a unique ID so the counting process doesn't get confused
df$patient_id <- 1:nrow(df)
df$date_start <- as.Date(df$best_start)
df$date_end <- as.Date(df$end_date)
# TX Duration = Time until abandonment (or time until treatment completed)
df$tx_dur <- as.numeric(difftime(df$date_end, df$date_start, units="days")) / 365.25

# Standardize Follow-up Horizon (2 Years)
df$event_d <- ifelse(df$time_d > 2.0, 0, df$event_d)
df$time_d <- ifelse(df$time_d > 2.0, 2.0, df$time_d)

cat("\n--- 2. Building Counting Process Data (tstart, tstop) ---\n")
# Everyone starts in Maintained Care!
df_p1 <- df
df_p1$tstart <- 0
# They remain in Maintained Care until the end of follow-up, OR until the day they abandon
df_p1$tstop <- ifelse(df_p1$itt_group == "Loss to follow-up" & df_p1$tx_dur < df_p1$time_d, 
                      df_p1$tx_dur, df_p1$time_d)
# If they drop out, they DO NOT die in Maintained Care (event=0 for this interval)
df_p1$event <- ifelse(df_p1$itt_group == "Loss to follow-up" & df_p1$tx_dur < df_p1$time_d, 
                      0, df_p1$event_d)
df_p1$group <- "Maintained Care (On Treatment or Completed)"


# For those who abandon, they spawn a second observation row starting at abandonment!
df_p2 <- df[df$itt_group == "Loss to follow-up" & df$tx_dur < df$time_d, ]
df_p2$tstart <- df_p2$tx_dur
# They remain in the Abandoned state until the end of follow-up
df_p2$tstop <- df_p2$time_d
df_p2$event <- df_p2$event_d
df_p2$group <- "Abandoned Treatment"

# Bind and throw away any logical noise (0 time intervals)
df_split <- bind_rows(df_p1, df_p2) %>% filter(tstop > tstart)

cat("\n--- 3. Fitting Time-Varying Nelson Aalen Estimator ---\n")
# By using tstart, tstop we correctly dynamically assign person-time
fit <- survfit(Surv(tstart, tstop, event) ~ group, data = df_split, id = patient_id)

df_plot <- data.frame(
  Time = fit$time,
  CumMort = fit$cumhaz, # Cumulative Hazard closely approximates Crude Cumulative Incidence
  strata = rep(names(fit$strata), fit$strata)
)
df_plot$strata <- gsub("group=", "", df_plot$strata)
df_start <- data.frame(Time = 0, CumMort = 0, strata = c("Maintained Care (On Treatment or Completed)", "Abandoned Treatment"))
df_plot <- bind_rows(df_plot, df_start)

cat("\n--- 4. Visualizing Dynamic Trajectories ---\n")
p1 <- ggplot(df_plot, aes(x = Time, y = CumMort, color = strata, linetype = strata)) +
  geom_step(linewidth = 1.3) +
  scale_color_manual(values = c("Maintained Care (On Treatment or Completed)" = "#2c3e50", "Abandoned Treatment" = "#e74c3c")) +
  scale_y_continuous(labels = scales::percent, breaks=seq(0, 0.15, 0.02), limits=c(0, 0.14)) +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  labs(
    title = "Time-Varying Crude Mortality (Person-Time Accounting)",
    subtitle = "Everyone starts in 'Maintained Care'. Only once an individual abandons do they enter the 'Abandoned' track.\nMaintained Care fully absorbs immediate early deaths.",
    x = "Years Since Tuberculosis Diagnosis",
    y = "Cumulative Hazard of Mortality (Crude)",
    color = "Current Patient State", linetype = "Current Patient State"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = c(0.3, 0.85),
    legend.title = element_blank(),
    legend.key.width = unit(2, "cm")
  )

ts <- as.numeric(Sys.time())
fname <- sprintf("figure_time_varying_crude_%.0f.png", ts)
artifact_dir <- "/Users/jasonandrews/.gemini/antigravity/brain/c053ef30-5842-41b7-b342-bf735650d865"

ggsave(file.path("ITT_Analysis/results", fname), plot=p1, width=11, height=7, dpi=300, bg="white")
ggsave(file.path(artifact_dir, fname), plot=p1, width=11, height=7, dpi=300, bg="white")

cat("\n--- 5. Updating Walkthrough Cache ---\n")
wt_path <- file.path(artifact_dir, "walkthrough.md")
if(file.exists(wt_path)){
    writeLines(paste0(paste(readLines(wt_path, warn=F), collapse="\n"), 
               "\n![Time-Varying Crude Plot](", file.path(artifact_dir, fname), ")\n"), wt_path)
}
cat("Success! Generated", fname, "\n")
