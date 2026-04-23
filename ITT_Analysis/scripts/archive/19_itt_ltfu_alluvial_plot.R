library(dplyr)
library(ggplot2)
library(tidyr)
# Try loading ggalluvial, install if missing
if(!require(ggalluvial)) {
    install.packages("ggalluvial", repos="http://cran.us.r-project.org")
    library(ggalluvial)
}

cat("Loading ITT cohort data...\n")
df <- read.csv("Abandonment Paper/ITT_Analysis/data/harmonized/itt_cohort.csv", stringsAsFactors = FALSE)

# Filter for the Abandonment group (LTFU)
df_ltfu <- df %>% filter(itt_group == "LTFU")
cat("Total individuals in LTFU group:", nrow(df_ltfu), "\n")

# Calculate medians for annotations
median_time_to_retreatment <- median(df_ltfu$time_r[df_ltfu$event_r == 1], na.rm = TRUE)
median_time_to_death_no_retreatment <- median(df_ltfu$time_d[df_ltfu$event_d == 1 & df_ltfu$event_r == 0], na.rm = TRUE)
median_time_to_death_with_retreatment <- median(df_ltfu$time_d[df_ltfu$event_d == 1 & df_ltfu$event_r == 1], na.rm = TRUE)

cat("Median time to retreatment (months):", median_time_to_retreatment, "\n")
cat("Median time to death (no retreatment, months):", median_time_to_death_no_retreatment, "\n")
cat("Median time to death (with retreatment, months):", median_time_to_death_with_retreatment, "\n")

# Categorize flows
df_flow <- df_ltfu %>%
  mutate(
    Retreatment_Status = ifelse(event_r == 1, "Retreated", "No Retreatment"),
    Final_Status = ifelse(event_d == 1, "Died", "Survived")
  ) %>%
  group_by(Retreatment_Status, Final_Status) %>%
  summarise(freq = n(), .groups = "drop")

print(df_flow)

# Create the alluvial plot
p <- ggplot(df_flow,
       aes(y = freq, axis1 = Retreatment_Status, axis2 = Final_Status)) +
  geom_alluvium(aes(fill = Final_Status), width = 1/12) +
  geom_stratum(width = 1/4, fill = "grey90", color = "grey20") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_fill_manual(values = c("Died" = "#E31A1C", "Survived" = "#1F78B4")) +
  scale_x_discrete(limits = c("Retreatment Status", "Final Survival Status"), expand = c(.05, .05)) +
  theme_minimal() +
  labs(title = "Clinical Flow Post-Abandonment",
       subtitle = sprintf("Following %d individuals who abandoned their first TB treatment\nMedian time to Retreatment: %.1f months\nMedian time to Death (overall): %.1f months", 
                          nrow(df_ltfu), median_time_to_retreatment, median(df_ltfu$time_d[df_ltfu$event_d == 1], na.rm=TRUE)),
       y = "Number of Individuals",
       fill = "Final Status") +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12))

output_path <- "Abandonment Paper/ITT_Analysis/results/19_ltfu_post_abandonment_flow.png"
ggsave(output_path, p, width = 8, height = 6, dpi = 300, bg="white")
cat("Plot saved to", output_path, "\n")

