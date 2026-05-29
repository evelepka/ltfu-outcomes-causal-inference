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
library(patchwork)
library(cowplot)

cat("\n--- 1. Loading Data ---\n")
df_cohort <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)
df_timing <- read.csv("ITT_Analysis/results/target_trial_piecewise_hr.csv", stringsAsFactors = FALSE)
df_sub <- read.csv("ITT_Analysis/results/target_trial_subgroup_late.csv", stringsAsFactors = FALSE)

# Clean Timing Data
# Map Epoch to cleaner labels for plotting
df_timing <- df_timing %>%
  mutate(
    Month_Num = as.numeric(gsub("Month ", "", Trial)),
    Epoch = factor(Epoch, levels = c("Early (Months 0-6)", "Late (Months 6-60)"))
  )

# Clean Subgroup Data
df_sub <- df_sub %>%
  mutate(
    Subgroup = recode(Subgroup,
                      "age_group" = "Age Group",
                      "sex" = "Sex",
                      "hiv_aids" = "HIV/AIDS",
                      "homelessness" = "Homelessness"),
    Level = factor(Level, levels = c("15-24", "25-44", "45-64", "65+", "Female", "Male", "Negative", "Positive", "No", "Yes")),
    Plot_Group = Subgroup
  )
df_sub$Subgroup <- factor(df_sub$Subgroup, levels = c("Age Group", "Homelessness", "Sex", "HIV/AIDS"))


cat("\n--- 2. Generating Panel A: Raincloud Stack ---\n")
# Filter to LTFU cohort
ltfu <- df_cohort %>% filter(itt_group == "Loss to follow-up")

ltfu$abandon_months <- as.numeric(difftime(as.Date(ltfu$end_date), as.Date(ltfu$best_start), units="days")) / 30.4
df_ab <- ltfu %>% filter(abandon_months <= 6 & abandon_months >= 0)
df_rn <- ltfu %>% filter(event_rn == 1 & time_rn <= 12)
df_d <- ltfu %>% filter(event_d == 1 & time_d <= 12)

build_raincloud <- function(df, x_var, x_label, col_fill, x_limits, x_breaks) {
  max_dens <- max(density(df[[x_var]], na.rm=TRUE)$y)
  df$jitter_y <- runif(nrow(df), -max_dens * 0.25, -max_dens * 0.05)
  
  ggplot(df, aes_string(x = x_var)) +
    geom_density(fill = col_fill, color = NA, alpha = 0.6) +
    geom_point(aes(y = jitter_y), color = col_fill, alpha = 0.15, size = 0.6, position="jitter") +
    geom_hline(yintercept=0, color="black", linewidth=0.5) +
    scale_x_continuous(limits = x_limits, breaks = x_breaks) +
    labs(x = x_label, y = "Density") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(t=10, r=10, b=10, l=10)
    )
}

pRain1 <- build_raincloud(df_ab, "abandon_months", "Months to LTFU", "#e74c3c", c(0, 6), seq(0, 6, 1))
pRain2 <- build_raincloud(df_rn, "time_rn", "Years to Retreatment", "#f1c40f", c(0, 12), seq(0, 12, 2))
pRain3 <- build_raincloud(df_d, "time_d", "Years to Death", "#2c3e50", c(0, 12), seq(0, 12, 2))

pA_body <- cowplot::plot_grid(pRain1, pRain2, pRain3, ncol = 1, align = "v", axis = "l", rel_heights=c(1, 1, 1))


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
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  labs(
    x = "Years from LTFU",
    y = "Cumulative Mortality",
    color = "", linetype = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.25, 0.85),
    legend.background = element_rect(fill="white", color=NA),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=10, r=10, b=10, l=10)
  )


cat("\n--- 4. Generating Panel C: Target Trial Timing ---\n")
pC <- ggplot(df_timing, aes(x = Month_Num, y = HR, color = Epoch, group = Epoch)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, linewidth=0.8) +
  geom_line(linewidth = 1, alpha = 0.5) +
  geom_point(size = 4) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  coord_cartesian(ylim = c(0.4, 5.0)) +
  scale_x_continuous(breaks = 1:6, labels = paste("Month", 1:6)) +
  scale_color_manual(values = c("Early (Months 0-6)" = "#3498db", "Late (Months 6-60)" = "#e74c3c")) +
  labs(
    x = "Month of LTFU",
    y = "Adjusted Hazard Ratio (aHR)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t=10, r=10, b=10, l=10)
  )


cat("\n--- 5. Generating Panel D: Subgroups (Forest Plot) ---\n")
pD <- ggplot(df_sub, aes(x = HR, y = Level, color = Subgroup)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper), width = 0.2, linewidth=0.8) +
  geom_point(size = 4) +
  facet_grid(Plot_Group ~ ., scales = "free_y", space = "free_y", switch="y") +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8)) +
  coord_cartesian(xlim = c(0.75, 10.0)) +
  scale_color_manual(values = c("Age Group" = "#3498db", "Homelessness" = "#e67e22", "Sex" = "#9b59b6", "HIV/AIDS" = "#1abc9c")) +
  labs(
    x = "Adjusted Hazard Ratio (aHR)",
    y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", hjust=1),
    panel.spacing = unit(0.5, "lines"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=10, r=10, b=10, l=10)
  )


cat("\n--- 6. Assembling and Saving ---\n")
caption_text <- "Figure 1. (A) Timing of LTFU, return to care and mortality, among those LTFU; (B) Cumulative risk of mortality following LTFU overall, among people experiencing homelessness and among people living with HIV; (C) Early (0-6 months) and late (>6 months) mortality hazard by month of LTFU; (D) Late mortality hazard risk by subgroup."

# Absolute positioning of textGrobs with top anchoring to prevent descender cut-off
tA <- textGrob("A: Event Distributions\nTiming of LTFU", x = unit(0.05, "npc"), y = unit(0.9, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 16))
tB <- textGrob("B: Cumulative Risk of Mortality\nMortality post-LTFU (up to 5 years)", x = unit(0.05, "npc"), y = unit(0.9, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 16))
tC <- textGrob("C: Mortality Hazard by Month of LTFU\nSequential Target Trial Hazard Ratios", x = unit(0.05, "npc"), y = unit(0.9, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 16))
tD <- textGrob("D: Mortality Risk by Subgroup\nTarget Trial Interactions", x = unit(0.05, "npc"), y = unit(0.9, "npc"), just = c("left", "top"), gp = gpar(fontface = "bold", fontsize = 16))

# Assemble blocks with native plot widths (no cowplot margin padding synchronization)
blkA <- arrangeGrob(tA, pA_body, ncol=1, heights=c(0.12, 0.88))
blkC <- arrangeGrob(tC, pC, ncol=1, heights=c(0.12, 0.88))
blkB <- arrangeGrob(tB, pB, ncol=1, heights=c(0.12, 0.88))
blkD <- arrangeGrob(tD, pD, ncol=1, heights=c(0.12, 0.88))

col1 <- arrangeGrob(blkA, blkC, ncol=1, heights=c(1.3, 1))
col2 <- arrangeGrob(blkB, blkD, ncol=1, heights=c(1.3, 1))
matrix_grid <- arrangeGrob(col1, col2, ncol=2, widths=c(1, 1.1))

caption_grob <- textGrob(str_wrap(caption_text, width = 140), x=unit(0.05, "npc"), y=unit(0.9, "npc"), just=c("left", "top"), gp=gpar(fontsize=14))
final_plot <- arrangeGrob(matrix_grid, caption_grob, ncol=1, heights=c(1, 0.08))

# Dynamic timestamped filename
ts <- as.numeric(Sys.time())
fname <- sprintf("figure1_target_trials_%.0f.png", ts)
fname_pdf <- sprintf("figure1_target_trials_%.0f.pdf", ts)
artifact_dir <- "/tmp"
out_path_repo <- file.path("ITT_Analysis/results", fname)
out_path_artifact <- file.path(artifact_dir, fname)
out_path_pdf <- file.path("ITT_Analysis/results", fname_pdf)
out_path_web <- file.path("ITT_Analysis/results", "Figure_1_Final_Web.png")

ggsave(out_path_repo, plot = final_plot, width = 16, height = 12, dpi = 600, bg="white")
ggsave(out_path_pdf, plot = final_plot, width = 16, height = 12, bg="white")
ggsave(out_path_web, plot = final_plot, width = 16, height = 12, dpi = 60, bg="white")

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
