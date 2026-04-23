# ==============================================================================
# MAIN MANUSCRIPT FIGURES: CUMULATIVE INCIDENCE (MORTALITY & RETREATMENT)
# Generates high-resolution KM and CIF plots with 95% Confidence Intervals
# for the LTFU Cohort.
# ==============================================================================

library(survival)
library(cmprsk)
library(ggplot2)
library(dplyr)
library(scales)

# --- 1. Load Data ---
cat("Loading data...\n")
df <- read.csv("Abandonment Paper/ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)
# Filter strictly to LTFU cohort as requested
ltfu <- df %>% filter(itt_group == "Loss to follow-up")
cat("LTFU Cohort Size =", nrow(ltfu), "\n")

# Prepare output directory
out_dir <- "Abandonment Paper/ITT_Analysis/results/manuscript figures/cumulative_incidence/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 2. FIGURE 1: MORTALITY (KAPLAN-MEIER WITH 95% CI) ---
cat("Generating Mortality KM curve...\n")
km_fit <- survfit(Surv(time_d, event_d) ~ 1, data = ltfu)

# Extract data for ggplot
km_data <- data.frame(
  time = km_fit$time,
  surv = km_fit$surv,
  lower = km_fit$lower,
  upper = km_fit$upper
)
# Convert to cumulative incidence (1 - survival)
km_data <- km_data %>%
  mutate(
    cuminc = 1 - surv,
    ci_lower = 1 - upper, # Flipped because 1 - max_surv = min_cuminc
    ci_upper = 1 - lower
  )
# Add time 0
km_data <- rbind(data.frame(time=0, surv=1, lower=1, upper=1, cuminc=0, ci_lower=0, ci_upper=0), km_data)

p_mort <- ggplot(km_data, aes(x = time)) +
  geom_step(aes(y = cuminc), color = "#b2182b", linewidth = 1.2) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0, 0.40)) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 12)) +
  labs(
    title = "Cumulative incidence of mortality",
    subtitle = "Kaplan-Meier estimates",
    x = "Years since loss to follow-up",
    y = "Cumulative incidence"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face="bold"))

ggsave(paste0(out_dir, "Cumulative_Incidence_Mortality_12y.png"), plot = p_mort, width = 8, height = 6, dpi = 300)

# --- 3. RETREATMENT (COMPETING RISKS CIF WITH 95% CI) ---
cat("Generating Retreatment CIF curves...\n")
# fstatus: 0=censored, 1=retreatment, 2=death
cif_fit <- cuminc(ftime = ltfu$time_rn, fstatus = ltfu$event_rn, cencode = 0)

# Extracting retreatment ("1 1") and death ("1 2")
retr_obj <- cif_fit[["1 1"]]
death_obj <- cif_fit[["1 2"]]

cif_data <- rbind(
  data.frame(time = retr_obj$time, est = retr_obj$est, Event = "Retreatment"),
  data.frame(time = death_obj$time, est = death_obj$est, Event = "Death")
)

# Add time 0
cif_data <- rbind(
  data.frame(time=0, est=0, Event="Retreatment"),
  data.frame(time=0, est=0, Event="Death"),
  cif_data
)

cif_data$Event <- factor(cif_data$Event, levels=c("Retreatment", "Death"))

# FIGURE 2: RETREATMENT CIF 12 YEARS
p_retr_12y <- ggplot(cif_data, aes(x = time, color = Event)) +
  geom_step(aes(y = est), linewidth = 1.2) +
  scale_color_manual(name = NULL, values = c("Retreatment" = "#2166ac", "Death" = "#b2182b")) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0, 0.60)) +
  scale_x_continuous(breaks = seq(0, 12, by = 2), limits = c(0, 12)) +
  labs(
    title = "Cumulative incidence of retreatment",
    subtitle = "Aalen-Johansen estimates accounting for competing risk of death",
    x = "Years since loss to follow-up",
    y = "Cumulative incidence"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face="bold"), legend.position = "bottom")

ggsave(paste0(out_dir, "Cumulative_Incidence_Retreatment_12y.png"), plot = p_retr_12y, width = 8, height = 6, dpi = 300)

# FIGURE 3: RETREATMENT CIF 24 MONTHS (ZOOMED)
cif_data_months <- cif_data %>% mutate(time_months = time * 12)

p_retr_2y <- ggplot(cif_data_months, aes(x = time_months, color = Event)) +
  geom_step(aes(y = est), linewidth = 1.2) +
  scale_color_manual(name = NULL, values = c("Retreatment" = "#2166ac", "Death" = "#b2182b")) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0, 0.50)) +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0, 24)) +
  labs(
    title = "Early cumulative incidence of retreatment",
    subtitle = "Aalen-Johansen estimates accounting for competing risk of death",
    x = "Months since loss to follow-up",
    y = "Cumulative incidence"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face="bold"), legend.position = "bottom")

ggsave(paste0(out_dir, "Cumulative_Incidence_Retreatment_24m.png"), plot = p_retr_2y, width = 8, height = 6, dpi = 300)

cat("Successfully generated and saved all three publication-ready figures.\n")
