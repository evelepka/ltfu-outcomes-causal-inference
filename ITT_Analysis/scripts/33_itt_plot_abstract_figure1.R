# 33. Generate Abstract Figure 1 (Target Trials & Competing Risks & Rainclouds)
# ==============================================================================
# Generates a FOUR-panel publication figure.
# Panel A: 3 Stacked Raincloud plots (Abandonment, Retreatment, Death timing)
# Panel B: Cumulative Risk of Mortality (5 Year)
# Panel C: Target Trial Emulation HR mapping across Month 1-6
# Panel D: Subgroup Forest Plot
# ==============================================================================

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(stringr)
library(survival)

cat("\n--- 1. Loading Data ---\n")
df_cohort <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)
df_timing <- read.csv("ITT_Analysis/results/target_trial_6mo_array_hr.csv", stringsAsFactors = FALSE)
df_sub <- read.csv("ITT_Analysis/results/target_trial_subgroup_interactions.csv", stringsAsFactors = FALSE)

# Clean Timing Data
df_timing$Month_Num <- as.numeric(gsub("Month_", "", df_timing$Trial_Month))

# Clean Subgroup Data
df_sub <- df_sub %>%
  mutate(
    Subgroup = recode(Subgroup,
                      "age_group" = "Age Group",
                      "sex" = "Sex",
                      "hiv_aids" = "HIV/AIDS",
                      "homelessness" = "Homelessness"),
    Level = factor(Level, levels = c("15-24", "25-44", "45-64", "≥65", "Female", "Male", "Negative", "Positive", "No", "Yes")),
    Plot_Group = paste0(Subgroup, "\n(", Interaction_P, ")")
  )
df_sub$Subgroup <- factor(df_sub$Subgroup, levels = c("Age Group", "Homelessness", "Sex", "HIV/AIDS"))


cat("\n--- 2. Generating Panel A: Raincloud Stack ---\n")
# Filter to LTFU cohort
ltfu <- df_cohort %>% filter(itt_group == "Loss to follow-up")

ltfu$abandon_months <- as.numeric(difftime(as.Date(ltfu$end_date), as.Date(ltfu$best_start), units="days")) / 30.4
df_ab <- ltfu %>% filter(abandon_months <= 6 & abandon_months >= 0)
df_rn <- ltfu %>% filter(event_rn == 1 & time_rn <= 12)
df_d <- ltfu %>% filter(event_d == 1 & time_d <= 12)

build_raincloud <- function(df, x_var, x_label, title_text, col_fill, x_limits, x_breaks) {
  max_dens <- max(density(df[[x_var]], na.rm=TRUE)$y)
  df$jitter_y <- runif(nrow(df), -max_dens * 0.25, -max_dens * 0.05)
  
  ggplot(df, aes_string(x = x_var)) +
    geom_density(fill = col_fill, color = NA, alpha = 0.6) +
    geom_point(aes(y = jitter_y), color = col_fill, alpha = 0.15, size = 0.6, position="jitter") +
    geom_hline(yintercept=0, color="black", linewidth=0.5) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks) +
    labs(title = title_text, x = x_label, y = "Density") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )
}

pRain1 <- build_raincloud(df_ab, "abandon_months", "Months to Abandonment", "Timing of Abandonment", "#e74c3c", c(0, 6), seq(0, 6, 1))
pRain2 <- build_raincloud(df_rn, "time_rn", "Years to Retreatment", "Timing of Return to Care", "#f1c40f", c(0, 12), seq(0, 12, 2))
pRain3 <- build_raincloud(df_d, "time_d", "Years to Death", "Timing of Mortality", "#2c3e50", c(0, 12), seq(0, 12, 2))

pA <- arrangeGrob(
  textGrob("A: Event Distributions", gp=gpar(fontsize=16, fontface="bold"), hjust=0, x=0.05),
  pRain1, pRain2, pRain3, 
  ncol = 1, heights = c(0.15, 1, 1, 1)
)


cat("\n--- 3. Generating Panel B: Cumulative Risk (KM 5yr) ---\n")
fit_all <- survfit(Surv(time_d, event_d) ~ 1, data=ltfu)
fit_hiv <- survfit(Surv(time_d, event_d) ~ 1, data=ltfu[!is.na(ltfu$hiv_aids) & ltfu$hiv_aids=="Positive", ])
fit_hom <- survfit(Surv(time_d, event_d) ~ 1, data=ltfu[!is.na(ltfu$homelessness) & ltfu$homelessness=="Yes", ])

extract_km <- function(fit, label) {
  df <- data.frame(
    time = c(0, fit$time),
    cum_mort = c(0, 1 - fit$surv),
    group = label
  ) 
  # Bound at 5 years
  df <- df %>% filter(time <= 5)
  return(df)
}

df_km <- bind_rows(
  extract_km(fit_all, "Overall LTFU Cohort"),
  extract_km(fit_hiv, "People living with HIV"),
  extract_km(fit_hom, "People experiencing homelessness")
)
df_km$group <- factor(df_km$group, levels=c("Overall LTFU Cohort", "People experiencing homelessness", "People living with HIV"))

pB <- ggplot(df_km, aes(x = time, y = cum_mort, color = group, linetype = group)) +
  geom_step(linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent, limits = c(0, max(df_km$cum_mort) * 1.1)) +
  scale_x_continuous(breaks = 0:5, limits = c(0, 5)) +
  scale_color_manual(values=c("Overall LTFU Cohort" = "black", "People living with HIV" = "#e74c3c", "People experiencing homelessness" = "#3498db")) +
  scale_linetype_manual(values=c("solid", "dashed", "dotdash")) +
  labs(
    title = "B: Cumulative Risk of Mortality",
    subtitle = "Survival post-abandonment (up to 5 years)",
    x = "Years from Abandonment",
    y = "Cumulative Mortality",
    color = "", linetype = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = c(0.4, 0.85),
    legend.background = element_rect(fill="white", color=NA, alpha=0.8),
    panel.grid.minor = element_blank()
  )


cat("\n--- 4. Generating Panel C: Target Trial Timing ---\n")
pC <- ggplot(df_timing, aes(x = Month_Num, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, color = "#2c3e50", linewidth=0.8) +
  geom_line(color = "#34495e", linewidth = 1, alpha = 0.5) +
  geom_point(size = 4, color = "#e74c3c") +
  scale_y_continuous(limits = c(0.8, 5.0), breaks = seq(1, 5, 1)) +
  scale_x_continuous(breaks = 1:6, labels = paste("Month", 1:6)) +
  labs(
    title = "C: Penalty by Month of Abandonment",
    subtitle = "Sequential Target Trial Hazard Ratios",
    x = "Exact Month of Treatment Abandonment",
    y = "Adjusted Hazard Ratio (aHR)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


cat("\n--- 5. Generating Panel D: Subgroups (Forest Plot) ---\n")
pD <- ggplot(df_sub, aes(x = HR, y = Level, color = Subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_errorbar(aes(xmin = CI_L, xmax = CI_H), width = 0.2, linewidth=0.8) +
  geom_point(size = 4) +
  facet_grid(Plot_Group ~ ., scales = "free_y", space = "free_y", switch="y") +
  scale_x_continuous(limits = c(0.8, 5.0), breaks = seq(1, 5, 1)) +
  scale_color_manual(values = c("Age Group" = "#3498db", "Homelessness" = "#e67e22", "Sex" = "#9b59b6", "HIV/AIDS" = "#1abc9c")) +
  labs(
    title = "D: Competing Risks Modification",
    subtitle = "Target Trial Interactions",
    x = "Adjusted Hazard Ratio (aHR)",
    y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", hjust=1),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.minor = element_blank()
  )


cat("\n--- 6. Assembling and Saving ---\n")
# Combine into a 2x2 grid
final_plot <- arrangeGrob(pA, pB, pC, pD, ncol=2, nrow=2, heights = c(1, 1), widths = c(1, 1.1))

# Dynamic timestamped filename
ts <- as.numeric(Sys.time())
fname <- sprintf("figure1_target_trials_%.0f.png", ts)
artifact_dir <- "/Users/jasonandrews/.gemini/antigravity/brain/c053ef30-5842-41b7-b342-bf735650d865"
out_path_repo <- file.path("ITT_Analysis/results", fname)
out_path_artifact <- file.path(artifact_dir, fname)

ggsave(out_path_repo, plot = final_plot, width = 16, height = 12, dpi = 300, bg="white")
ggsave(out_path_artifact, plot = final_plot, width = 16, height = 12, dpi = 300, bg="white")

cat("\n--- 7. Updating Walkthrough Cache ---\n")
walkthrough_path <- file.path(artifact_dir, "walkthrough.md")
if (file.exists(walkthrough_path)) {
  wt_text <- readLines(walkthrough_path, warn=FALSE)
  wt_blob <- paste(wt_text, collapse = "\n")
  
  if (grepl("!\\[Figure 1\\]", wt_blob)) {
      new_blob <- str_replace(wt_blob, "!\\[Figure 1\\]\\(/Users/jasonandrews/[^)]+\\.png\\)", 
                                       paste0("![Figure 1](", out_path_artifact, ")"))
  } else {
      new_blob <- paste0(wt_blob, "\n\n![Figure 1](", out_path_artifact, ")\n")
  }
  
  writeLines(new_blob, walkthrough_path)
}

cat("Success! Generated", fname, "\n")
