# 19_Create_Stratified_Retreatment.R
# Stratified Cumulative Incidence of Retreatment (Accounting for Competing Risk of Death)
# Groups: HIV, Homelessness, Timing of LTFU

library(cmprsk)
library(ggplot2)
library(dplyr)
library(scales)

# 1. Load Data
cohort <- read.csv("Abandonment Paper/ITT_Analysis/data/itt_cohort.csv")
ltfu <- cohort %>% filter(itt_group == "Loss to follow-up")

# Directories
out_dir <- "Abandonment Paper/ITT_Analysis/results/manuscript figures/stratified_retreatment/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Helper function to extract CIF at exact time t for cause 1 (Retreatment)
get_cif_at <- function(cif, group_name, t) {
  obj <- cif[[paste0(group_name, " 1")]]
  if (is.null(obj)) return(NA)
  idx <- max(which(obj$time <= t), na.rm=TRUE)
  if (length(idx) == 0 || is.infinite(idx)) return(NA)
  return(obj$est[idx] * 100)
}

# Helper to build plot df for cause 1
get_plot_df <- function(cif, groups) {
  max_t <- max(ltfu$time_rn, na.rm=TRUE)
  do.call(rbind, lapply(groups, function(g) {
    obj <- cif[[paste0(g, " 1")]]
    if (is.null(obj)) return(NULL)
    d <- data.frame(time = obj$time, est = obj$est, group = as.character(g))
    last_est <- tail(obj$est, 1)
    # Add time 0 and max_t
    rbind(
      data.frame(time = 0, est = 0, group = as.character(g)),
      d,
      data.frame(time = max_t, est = last_est, group = as.character(g))
    )
  }))
}

# ------------------------------------------------------------------------------
# A) HIV Status
# ------------------------------------------------------------------------------
cat("Processing HIV...\n")
ltfu_hiv <- ltfu %>% filter(!is.na(hiv_aids) & hiv_aids != "")
cif_hiv <- cuminc(ftime = ltfu_hiv$time_rn, fstatus = ltfu_hiv$event_rn, group = ltfu_hiv$hiv_aids, cencode = 0)
pval_hiv <- cif_hiv$Tests["1", "pv"] # Gray's test p-value for cause 1
pval_text_hiv <- if(pval_hiv < 0.001) "p < 0.001" else paste0("p = ", round(pval_hiv, 3))

groups_hiv <- unique(ltfu_hiv$hiv_aids)
df_hiv <- get_plot_df(cif_hiv, groups_hiv)
df_hiv$group <- factor(df_hiv$group, levels = c("Negative", "Positive", "Unknown"))

p_hiv <- ggplot(df_hiv, aes(x = time, y = est, color = group)) +
  geom_step(linewidth = 1.2) +
  theme_classic() +
  annotate("text", x = 0.5, y = 0.75, label = paste0("Gray's test ", pval_text_hiv), size = 4.5, hjust = 0, fontface = "italic") +
  labs(
    title = "Cumulative incidence of retreatment by HIV status",
    subtitle = "Aalen-Johansen estimates accounting for competing risk of death",
    x = "Years since loss to follow-up",
    y = "Cumulative incidence"
  ) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0, 0.80)) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 12)) +
  scale_color_manual(name = NULL, 
                     values = c("Negative" = "#377eb8", "Positive" = "#e41a1c", "Unknown" = "#999999"),
                     labels = c("Negative" = "HIV-negative", "Positive" = "HIV-positive", "Unknown" = "Unknown")) +
  theme(plot.title = element_text(face="bold"), legend.position = "bottom")

ggsave(paste0(out_dir, "Retreatment_CIF_12y_by_HIV.png"), plot = p_hiv, width = 8, height = 6, dpi = 300)

# Table
res_hiv <- data.frame(
  Group = groups_hiv,
  Year_1 = sapply(groups_hiv, function(g) get_cif_at(cif_hiv, g, 1)),
  Year_2 = sapply(groups_hiv, function(g) get_cif_at(cif_hiv, g, 2)),
  Year_5 = sapply(groups_hiv, function(g) get_cif_at(cif_hiv, g, 5)),
  Year_10 = sapply(groups_hiv, function(g) get_cif_at(cif_hiv, g, 10)),
  Year_11.9 = sapply(groups_hiv, function(g) get_cif_at(cif_hiv, g, 11.9))
)
write.csv(res_hiv, paste0(out_dir, "retreatment_cif_hiv_values.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# B) Homelessness
# ------------------------------------------------------------------------------
cat("Processing Homelessness...\n")
ltfu_home <- ltfu %>% filter(!is.na(homelessness) & homelessness != "")
cif_home <- cuminc(ftime = ltfu_home$time_rn, fstatus = ltfu_home$event_rn, group = ltfu_home$homelessness, cencode = 0)
pval_home <- cif_home$Tests["1", "pv"]
pval_text_home <- if(pval_home < 0.001) "p < 0.001" else paste0("p = ", round(pval_home, 3))

groups_home <- unique(ltfu_home$homelessness)
df_home <- get_plot_df(cif_home, groups_home)
df_home$group <- factor(df_home$group, levels = c("No", "Yes"))

p_home <- ggplot(df_home, aes(x = time, y = est, color = group)) +
  geom_step(linewidth = 1.2) +
  theme_classic() +
  annotate("text", x = 0.5, y = 0.75, label = paste0("Gray's test ", pval_text_home), size = 4.5, hjust = 0, fontface = "italic") +
  labs(
    title = "Cumulative incidence of retreatment by homelessness status",
    subtitle = "Aalen-Johansen estimates accounting for competing risk of death",
    x = "Years since loss to follow-up",
    y = "Cumulative incidence"
  ) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0, 0.80)) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 12)) +
  scale_color_manual(name = NULL, 
                     values = c("No" = "#377eb8", "Yes" = "#e41a1c"),
                     labels = c("No" = "Not homeless", "Yes" = "Homeless")) +
  theme(plot.title = element_text(face="bold"), legend.position = "bottom")

ggsave(paste0(out_dir, "Retreatment_CIF_12y_by_Homelessness.png"), plot = p_home, width = 8, height = 6, dpi = 300)

# Table
res_home <- data.frame(
  Group = groups_home,
  Year_1 = sapply(groups_home, function(g) get_cif_at(cif_home, g, 1)),
  Year_2 = sapply(groups_home, function(g) get_cif_at(cif_home, g, 2)),
  Year_5 = sapply(groups_home, function(g) get_cif_at(cif_home, g, 5)),
  Year_10 = sapply(groups_home, function(g) get_cif_at(cif_home, g, 10)),
  Year_11.9 = sapply(groups_home, function(g) get_cif_at(cif_home, g, 11.9))
)
write.csv(res_home, paste0(out_dir, "retreatment_cif_homelessness_values.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# C) Timing of LTFU
# ------------------------------------------------------------------------------
cat("Processing Timing of LTFU...\n")
ltfu_tim <- ltfu %>% filter(!is.na(tx_month_grp) & tx_month_grp != "")
ltfu_tim$tx_month_grp <- factor(ltfu_tim$tx_month_grp, levels = c("< 1 month", "1 to <2 months", "2 to <4 months", "≥ 4 months"))
cif_tim <- cuminc(ftime = ltfu_tim$time_rn, fstatus = ltfu_tim$event_rn, group = ltfu_tim$tx_month_grp, cencode = 0)
pval_tim <- cif_tim$Tests["1", "pv"]
pval_text_tim <- if(pval_tim < 0.001) "p < 0.001" else paste0("p = ", round(pval_tim, 3))

groups_tim <- unique(as.character(ltfu_tim$tx_month_grp))
df_tim <- get_plot_df(cif_tim, groups_tim)
df_tim$group <- factor(df_tim$group, levels = levels(ltfu_tim$tx_month_grp))

p_tim <- ggplot(df_tim, aes(x = time, y = est, color = group)) +
  geom_step(linewidth = 1.2) +
  theme_classic() +
  annotate("text", x = 0.5, y = 0.75, label = paste0("Gray's test ", pval_text_tim), size = 4.5, hjust = 0, fontface = "italic") +
  labs(
    title = "Cumulative incidence of retreatment by timing of loss to follow-up",
    subtitle = "Aalen-Johansen estimates accounting for competing risk of death",
    x = "Years since loss to follow-up",
    y = "Cumulative incidence"
  ) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0, 0.80)) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 12)) +
  scale_color_manual(name = NULL, 
                     values = c("< 1 month" = "#e41a1c", "1 to <2 months" = "#ff7f00", "2 to <4 months" = "#4daf4a", "≥ 4 months" = "#377eb8"),
                     labels = c("< 1 month" = "< 1 month", "1 to <2 months" = "1-2 months", "2 to <4 months" = "2-4 months", "≥ 4 months" = "> 4 months")) +
  theme(plot.title = element_text(face="bold"), legend.position = "bottom")

ggsave(paste0(out_dir, "Retreatment_CIF_12y_by_TimingLTFU.png"), plot = p_tim, width = 8, height = 6, dpi = 300)

# Table
res_tim <- data.frame(
  Group = groups_tim,
  Year_1 = sapply(groups_tim, function(g) get_cif_at(cif_tim, g, 1)),
  Year_2 = sapply(groups_tim, function(g) get_cif_at(cif_tim, g, 2)),
  Year_5 = sapply(groups_tim, function(g) get_cif_at(cif_tim, g, 5)),
  Year_10 = sapply(groups_tim, function(g) get_cif_at(cif_tim, g, 10)),
  Year_11.9 = sapply(groups_tim, function(g) get_cif_at(cif_tim, g, 11.9))
)
write.csv(res_tim, paste0(out_dir, "retreatment_cif_timing_values.csv"), row.names = FALSE)

cat("Successfully generated 3 stratified plots and tables for Retreatment.\n")
