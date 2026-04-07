# 33. Generate Abstract Figure 1 (Target Trials & Competing Risks)
# ==============================================================================
# Generates a two-panel publication figure matching the Target Trial outputs.
# Panel A: Target Trial Emulation HR mapping across Month 1-6
# Panel B: Subgroup Forest Plot
# ==============================================================================

library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)

cat("\n--- 1. Loading Data ---\n")
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
    # Group label for interaction text
    Plot_Group = paste0(Subgroup, "\n(", Interaction_P, ")")
  )

# Reverse order so Age is at the top of the plot
df_sub$Subgroup <- factor(df_sub$Subgroup, levels = c("Age Group", "Homelessness", "Sex", "HIV/AIDS"))

cat("\n--- 2. Generating Panel A: Timing ---\n")
pA <- ggplot(df_timing, aes(x = Month_Num, y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, color = "#2c3e50", size=0.8) +
  geom_line(color = "#34495e", size = 1, alpha = 0.5) +
  geom_point(size = 4, color = "#e74c3c") +
  scale_y_continuous(limits = c(0.8, 5.0), breaks = seq(1, 5, 1)) +
  scale_x_continuous(breaks = 1:6, labels = paste("Month", 1:6)) +
  labs(
    title = "A: Relative Penalty by Month of Abandonment",
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

cat("\n--- 3. Generating Panel B: Subgroups (Forest Plot) ---\n")
pB <- ggplot(df_sub, aes(x = HR, y = Level, color = Subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
  geom_errorbarh(aes(xmin = CI_L, xmax = CI_H), height = 0.2, size=0.8) +
  geom_point(size = 4) +
  facet_grid(Plot_Group ~ ., scales = "free_y", space = "free_y", switch="y") +
  scale_x_continuous(limits = c(0.8, 5.0), breaks = seq(1, 5, 1)) +
  scale_color_manual(values = c("Age Group" = "#3498db", "Homelessness" = "#e67e22", "Sex" = "#9b59b6", "HIV/AIDS" = "#1abc9c")) +
  labs(
    title = "B: Competing Risks Modification",
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

cat("\n--- 4. Assembling and Saving ---\n")
final_plot <- arrangeGrob(pA, pB, widths = c(1, 1.2), ncol=2)

# Dynamic timestamped filename
ts <- as.numeric(Sys.time())
fname <- sprintf("figure1_target_trials_%.0f.png", ts)
artifact_dir <- "/Users/jasonandrews/.gemini/antigravity/brain/c053ef30-5842-41b7-b342-bf735650d865"
out_path_repo <- file.path("ITT_Analysis/results", fname)
out_path_artifact <- file.path(artifact_dir, fname)

ggsave(out_path_repo, final_plot, width = 14, height = 7, dpi = 300, bg="white")
ggsave(out_path_artifact, final_plot, width = 14, height = 7, dpi = 300, bg="white")

cat("\n--- 5. Updating Walkthrough Cache ---\n")
walkthrough_path <- file.path(artifact_dir, "walkthrough.md")
if (file.exists(walkthrough_path)) {
  wt_text <- readLines(walkthrough_path, warn=FALSE)
  wt_blob <- paste(wt_text, collapse = "\n")
  
  # Regex cache busting
  if (grepl("!\\[Figure 1\\]", wt_blob)) {
      new_blob <- str_replace(wt_blob, "!\\[Figure 1\\]\\(/Users/jasonandrews/[^)]+\\.png\\)", 
                                       paste0("![Figure 1](", out_path_artifact, ")"))
  } else {
      # Append to the end
      new_blob <- paste0(wt_blob, "\n\n![Figure 1](", out_path_artifact, ")\n")
  }
  
  writeLines(new_blob, walkthrough_path)
  cat("Successfully injected", fname, "into walkthrough.md\n")
}

cat("Success! Generated", fname, "\n")
