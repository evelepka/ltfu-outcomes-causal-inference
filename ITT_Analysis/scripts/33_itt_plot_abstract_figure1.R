# 33. Generate Abstract Figure 1 (Target Trials & Competing Risks & Rainclouds)
# ==============================================================================
# Generates a three-panel publication figure.
# Panel A: 3 Stacked Raincloud plots (Abandonment, Retreatment, Death timing)
# Panel B: Target Trial Emulation HR mapping across Month 1-6
# Panel C: Subgroup Forest Plot
# ==============================================================================

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(stringr)

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

# 1. Timing of Abandonment (Months 0-6)
ltfu$abandon_months <- as.numeric(difftime(as.Date(ltfu$end_date), as.Date(ltfu$best_start), units="days")) / 30.4
df_ab <- ltfu %>% filter(abandon_months <= 6 & abandon_months >= 0)

# 2. Timing of Retreatment (Years 0-2)
df_rn <- ltfu %>% filter(event_rn == 1 & time_rn <= 2)

# 3. Timing of Death (Years 0-12)
df_d <- ltfu %>% filter(event_d == 1 & time_d <= 12)

# Helper function to generate pseudo-raincloud (density + bottom jitter)
build_raincloud <- function(df, x_var, x_label, title_text, col_fill, x_limits, x_breaks) {
  # Calculate density height to scale the jitter appropriately
  max_dens <- max(density(df[[x_var]], na.rm=TRUE)$y)
  # Create a jitter strip scaled to 25% of the max density, directly underneath
  df$jitter_y <- runif(nrow(df), -max_dens * 0.25, -max_dens * 0.05)
  
  ggplot(df, aes_string(x = x_var)) +
    geom_density(fill = col_fill, color = NA, alpha = 0.6) +
    geom_point(aes(y = jitter_y), color = col_fill, alpha = 0.15, size = 0.6, position="jitter") +
    geom_hline(yintercept=0, color="black", size=0.5) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks) +
    labs(title = title_text, x = x_label, y = "Density") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )
}

pRain1 <- build_raincloud(df_ab, "abandon_months", "Months to Abandonment", "Timing of Abandonment", "#e74c3c", c(0, 6), seq(0, 6, 1))
pRain2 <- build_raincloud(df_rn, "time_rn", "Years to Retreatment", "Timing of Return to Care", "#f1c40f", c(0, 2), seq(0, 2, 0.5))
pRain3 <- build_raincloud(df_d, "time_d", "Years to Death", "Timing of Mortality", "#2c3e50", c(0, 12), seq(0, 12, 2))

pA <- arrangeGrob(
  textGrob("A: Event Distributions", gp=gpar(fontsize=16, fontface="bold"), hjust=0, x=0.05),
  pRain1, pRain2, pRain3, 
  ncol = 1, heights = c(0.15, 1, 1, 1)
)


cat("\n--- 3. Generating Panel B: Target Trial Timing ---\n")
pB <- ggplot(df_timing, aes(x = Month_Num, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, color = "#2c3e50", size=0.8) +
  geom_line(color = "#34495e", size = 1, alpha = 0.5) +
  geom_point(size = 4, color = "#e74c3c") +
  scale_y_continuous(limits = c(0.8, 5.0), breaks = seq(1, 5, 1)) +
  scale_x_continuous(breaks = 1:6, labels = paste("Month", 1:6)) +
  labs(
    title = "B: Penalty by Month of Abandonment",
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


cat("\n--- 4. Generating Panel C: Subgroups (Forest Plot) ---\n")
pC <- ggplot(df_sub, aes(x = HR, y = Level, color = Subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
  geom_errorbarh(aes(xmin = CI_L, xmax = CI_H), height = 0.2, size=0.8) +
  geom_point(size = 4) +
  facet_grid(Plot_Group ~ ., scales = "free_y", space = "free_y", switch="y") +
  scale_x_continuous(limits = c(0.8, 5.0), breaks = seq(1, 5, 1)) +
  scale_color_manual(values = c("Age Group" = "#3498db", "Homelessness" = "#e67e22", "Sex" = "#9b59b6", "HIV/AIDS" = "#1abc9c")) +
  labs(
    title = "C: Competing Risks Modification",
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


cat("\n--- 5. Assembling and Saving ---\n")
# Combine the three panels side-by-side
final_plot <- arrangeGrob(pA, pB, pC, widths = c(1, 1, 1.2), ncol=3)

# Dynamic timestamped filename
ts <- as.numeric(Sys.time())
fname <- sprintf("figure1_target_trials_%.0f.png", ts)
artifact_dir <- "/Users/jasonandrews/.gemini/antigravity/brain/c053ef30-5842-41b7-b342-bf735650d865"
out_path_repo <- file.path("ITT_Analysis/results", fname)
out_path_artifact <- file.path(artifact_dir, fname)

ggsave(out_path_repo, plot = final_plot, width = 18, height = 7, dpi = 300, bg="white")
ggsave(out_path_artifact, plot = final_plot, width = 18, height = 7, dpi = 300, bg="white")

cat("\n--- 6. Updating Walkthrough Cache ---\n")
walkthrough_path <- file.path(artifact_dir, "walkthrough.md")
if (file.exists(walkthrough_path)) {
  wt_text <- readLines(walkthrough_path, warn=FALSE)
  wt_blob <- paste(wt_text, collapse = "\n")
  
  # Regex cache busting
  if (grepl("!\\[Figure 1\\]", wt_blob)) {
      new_blob <- str_replace(wt_blob, "!\\[Figure 1\\]\\(/Users/jasonandrews/[^)]+\\.png\\)", 
                                       paste0("![Figure 1](", out_path_artifact, ")"))
  } else {
      # Append to the end if not found
      new_blob <- paste0(wt_blob, "\n\n![Figure 1](", out_path_artifact, ")\n")
  }
  
  writeLines(new_blob, walkthrough_path)
  cat("Successfully injected", fname, "into walkthrough.md\n")
}

cat("Success! Generated", fname, "\n")
