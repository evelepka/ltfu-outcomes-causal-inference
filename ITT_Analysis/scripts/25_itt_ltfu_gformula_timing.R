# 25. ITT LTFU G-Formula for Retreatment Timing
# ==============================================================================
# Purpose: Estimate how the timing of return to treatment impacts 10-year mortality.
# Method: Parametric G-computation on the LTFU cohort.
# ==============================================================================

library(dplyr)
library(survival)
library(tidyr)
library(splines)
library(ggplot2)

cat("\n--- 1. Loading and Preparing ITT LTFU Cohort ---\n")
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

df_ltfu <- df %>% filter(itt_group == "Loss to follow-up" & time_d > 0)

clean_data <- function(data) {
  data %>% 
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    age_group = factor(age_group, levels = c("15-24", "25-44", "45-64", "≥65")),
    hiv_aids = factor(hiv_aids, levels = c("Negative", "Positive")),
    race_clean = factor(race_clean, levels = c("White", "Black or Mixed", "Other")),
    edu_clean = factor(edu_clean, levels = c("≥ 12 years", "8 - 11 years", "≤ 7 years", "None")),
    clinical_clean = factor(clinical_clean, levels = c("Pulmonary", "Extrapulmonary", "Pulmonary and Extrapulmonary or disseminated")),
    dot_status = factor(dot_status, levels = c("No", "Yes")),
    diabetes = factor(diabetes, levels = c("No", "Yes")),
    alcohol = factor(alcohol, levels = c("No", "Yes")),
    drug_use = factor(drug_use, levels = c("No", "Yes")),
    incarcerated = factor(incarcerated, levels = c("No", "Yes")),
    homelessness = factor(homelessness, levels = c("No", "Yes")),
    hosp_admission = factor(hosp_admission, levels = c("No", "Yes")),
    
    time_rn = ifelse(event_rn == 1 & time_rn > time_d, time_d, time_rn),
    event_rn = ifelse(event_rn == 1 & time_rn > time_d, 0, event_rn),
    
    # Cap at 10 years (120 months) + small buffer
    month_d = pmin(120, pmax(1, ceiling(time_d * 12))),
    month_rn = ifelse(event_rn == 1, pmin(120, pmax(1, ceiling(time_rn * 12))), pmin(120, pmax(1, ceiling(time_d * 12))))
  )
}

df_clean <- clean_data(df_ltfu) %>% tidyr::drop_na(age_group, sex, race_clean, edu_clean, hiv_aids, diabetes, alcohol, drug_use, incarcerated, homelessness, hosp_admission, clinical_clean, dot_status)

cat("Complete Case LTFU Cohort Size:", nrow(df_clean), "\n")

cat("\n--- 2. Building Person-Month Dataset & Fitting Outcome Model ---\n")
df_pm <- df_clean %>%
  uncount(month_d, .id = "month") %>%
  group_by(sinan_clean) %>%
  mutate(
    retreatment_status = ifelse(month >= month_rn & event_rn == 1, 1, 0),
    death = ifelse(month == max(month) & event_d == 1, 1, 0)
  ) %>%
  ungroup()

# Outcome model (pooled cloglog approximates hazard)
fit_death <- glm(
  death ~ retreatment_status * ns(month, df=3) + age_group + sex + race_clean +
    edu_clean + hiv_aids + diabetes + alcohol + drug_use +
    incarcerated + homelessness + hosp_admission + clinical_clean + dot_status,
  family = binomial(link = "cloglog"), data = df_pm
)

cat("Model fitted.\n")

cat("\n--- 3. Simulating Counterfactual Curves ---\n")
simulate_survival <- function(model, baseline_data, intervention_type) {
  # intervention_type: "natural", "never", "month_3", "month_6"
  months <- 1:120
  n <- nrow(baseline_data)
  curr_surv <- rep(1, n)
  pop_mort <- rep(0, 120)
  
  cf_data <- baseline_data
  
  for (m in months) {
    cf_data$month <- m
    if (intervention_type == "never") {
      cf_data$retreatment_status <- 0
    } else if (intervention_type == "month_3") {
      cf_data$retreatment_status <- ifelse(m >= 3, 1, 0)
    } else if (intervention_type == "month_6") {
      cf_data$retreatment_status <- ifelse(m >= 6, 1, 0)
    } else {
      # Natural course: observe actual data
      cf_data$retreatment_status <- ifelse(m >= cf_data$month_rn & cf_data$event_rn == 1, 1, 0)
    }
    
    p_death <- predict(model, newdata = cf_data, type = "response")
    curr_surv <- curr_surv * (1 - p_death)
    pop_mort[m] <- 1 - mean(curr_surv)
  }
  return(pop_mort)
}

cat("Simulating 'Never Return' Scenario...\n")
curve_never <- simulate_survival(fit_death, df_clean, "never")
cat("Simulating 'Return at Month 3' Scenario...\n")
curve_m3 <- simulate_survival(fit_death, df_clean, "month_3")
cat("Simulating 'Return at Month 6' Scenario...\n")
curve_m6 <- simulate_survival(fit_death, df_clean, "month_6")
cat("Simulating 'Natural Course' Scenario...\n")
curve_nc <- simulate_survival(fit_death, df_clean, "natural")

cat("\n--- 4. Saving & Plotting Results ---\n")
plot_df <- data.frame(
  Year = (1:120) / 12,
  Never_Return = curve_never * 100,
  Return_Month_3 = curve_m3 * 100,
  Return_Month_6 = curve_m6 * 100,
  Natural_Course = curve_nc * 100
)

# Export csv
write.csv(plot_df, "ITT_Analysis/results/g_formula_retreatment_timing.csv", row.names=FALSE)

plot_long <- plot_df %>% pivot_longer(cols = -Year, names_to = "Scenario", values_to = "Mortality")

plot_long$Scenario <- factor(plot_long$Scenario, 
                             levels = c("Never_Return", "Natural_Course", "Return_Month_6", "Return_Month_3"),
                             labels = c("Never Return", "Natural Course (Observed)", "Mandatory Return at 6 Months", "Mandatory Return at 3 Months"))

ts <- as.numeric(Sys.time())
img_filename <- sprintf("gformula_timing_%s.png", round(ts))
img_path <- file.path("ITT_Analysis/results", img_filename)

png(img_path, width=2400, height=1800, res=300)
p <- ggplot(plot_long, aes(x = Year, y = Mortality, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  theme_classic() +
  scale_color_manual(values = c("Never Return" = "#d73027", 
                                "Natural Course (Observed)" = "#969696", 
                                "Mandatory Return at 6 Months" = "#4575b4", 
                                "Mandatory Return at 3 Months" = "#313695")) +
  scale_y_continuous(limits = c(0, max(plot_df$Never_Return)+1), breaks = seq(0, 15, 2)) +
  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(title = "Causal Impact of Retreatment Timing on Mortality",
       subtitle = "10-Year Cumulative Incidence via Parametric G-Formula (LTFU Cohort)",
       x = "Years since Treatment Abandonment",
       y = "Cumulative Incidence of Death (%)") +
  theme(legend.position = "bottom", legend.title=element_blank(),
        plot.title = element_text(face="bold", size=14),
        axis.title = element_text(size=12))
print(p)
dev.off()

cat(sprintf("\nPlot saved to %s\n", img_path))

# Dynamic Walkthrough Update Subroutine
abs_img_path <- normalizePath(img_path)
walkthrough_path <- "/Users/jasonandrews/.gemini/antigravity/brain/c053ef30-5842-41b7-b342-bf735650d865/walkthrough.md"

if (file.exists(walkthrough_path)) {
    wt_text <- readLines(walkthrough_path)
    # Search for an existing plot link replacing it if matched
    img_pattern <- "!\\[.*\\]\\(.*/gformula_timing_.*\\.png\\)"
    new_img_link <- sprintf("![G-Formula Timing Plot](%s)", abs_img_path)
    if (any(grepl(img_pattern, wt_text))) {
        wt_text <- gsub(img_pattern, new_img_link, wt_text)
    } else {
        # append if not found
        wt_text <- c(wt_text, "", new_img_link)
    }
    writeLines(wt_text, walkthrough_path)
    cat("Automatically updated walkthrough.md with the new cache-busted image.\n")
} else {
    cat("Warning: walkthrough.md not found, could not update image link.\n")
}
