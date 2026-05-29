library(dplyr)
library(survival)
library(muhaz)
library(survRM2)
library(ggplot2)
library(patchwork)

df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)

# Base restrictions and cleanup
df <- df %>% filter(time_d_tx >= 0, itt_group %in% c("Loss to follow-up", "Non-LTFU")) %>%
  mutate(
    event_cap = ifelse(time_d_tx > 5, 0, event_d),
    time_cap = pmin(time_d_tx, 5),
    is_abandon = ifelse(itt_group == "Loss to follow-up", 1, 0)
  )

# ---------------------------------------------------------
# 1. Smoothed Hazard Function using muhaz
# ---------------------------------------------------------
cat("Calculating Smoothed Hazards...\n")
df_ltfu <- filter(df, is_abandon == 1)
df_mnt  <- filter(df, is_abandon == 0)

smooth_ltfu <- muhaz(df_ltfu$time_cap, df_ltfu$event_cap, max.time = 5)
smooth_mnt  <- muhaz(df_mnt$time_cap, df_mnt$event_cap, max.time = 5)

hz_df <- data.frame(
  time = c(smooth_ltfu$est.grid, smooth_mnt$est.grid),
  hazard = c(smooth_ltfu$haz.est, smooth_mnt$haz.est),
  group = c(rep("Abandoned", length(smooth_ltfu$est.grid)), rep("Maintained", length(smooth_mnt$est.grid)))
)

p_hazard <- ggplot(hz_df, aes(x = time, y = hazard, color = group)) +
  geom_line(linewidth = 1.2) +
  theme_classic() +
  scale_color_manual(values = c("Abandoned" = "#d95f02", "Maintained" = "#1b9e77")) +
  labs(title = "Panel A: Instantaneous Hazard Rate",
       subtitle = "Illustrates the visual crossover 'delay' mechanism of abandonment",
       x = "Years Since Treatment Start", y = "Smoothed Hazard") +
  theme(legend.position = "bottom", legend.title = element_blank())

# ---------------------------------------------------------
# 2. Dynamic Point-in-time RMST
# ---------------------------------------------------------
cat("Calculating Dynamic RMST Trajectories...\n")
taus <- seq(0.5, 5, by = 0.5)
rmst_diffs <- numeric(length(taus))
rmst_se <- numeric(length(taus))

# Maintained is Arm=1, Abandoned is Arm=0
df$arm <- ifelse(df$is_abandon == 0, 1, 0)

for (i in seq_along(taus)) {
  m <- rmst2(time = df$time_cap, status = df$event_cap, arm = df$arm, tau = taus[i])
  diff_res <- m$unadjusted.result[1,] # line 1 is RMST(1)-RMST(0)
  rmst_diffs[i] <- diff_res[1]
  
  # CI tracking
  rmst_se[i] <- (diff_res[3] - diff_res[2]) / (2 * 1.96)
}

df_rmst <- data.frame(tau = taus, diff = rmst_diffs, se = rmst_se) %>%
  mutate(
    ci_l = diff - 1.96 * se,
    ci_u = diff + 1.96 * se
  )

p_rmst <- ggplot(df_rmst, aes(x = tau, y = diff * 365.25)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = ci_l * 365.25, ymax = ci_u * 365.25), fill = "blue", alpha = 0.2) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_point(size = 3, color = "blue", fill = "white", shape = 21) +
  theme_classic() +
  labs(title = "Panel B: Days of Life Gained over Time",
       subtitle = "Cumulative life-days gained by maintaining treatment across the follow-up window",
       x = "Follow-up Window Length (Years)", 
       y = "Days of Life Gained (RMST Diff)")

# ---------------------------------------------------------
# 3. Piecewise Cox Proportional Hazards
# ---------------------------------------------------------
cat("Running Piecewise Cox Regression...\n")
# We will split at 6 months (0.5 years) and evaluate the crude HR in [0, 0.5) vs [0.5, 5.0]
df_split <- survSplit(Surv(time_cap, event_cap) ~ ., data = df, cut = c(0.5), episode = "tgroup")

# tgroup: 1 = [0, 0.5), 2 = [0.5, 5]
fit_pw <- coxph(Surv(tstart, time_cap, event_cap) ~ is_abandon:strata(tgroup), data = df_split)
s_pw <- summary(fit_pw)

pw_res <- data.frame(
  Epoch = c("Months 0-6", "Months 6-60"),
  HR = exp(s_pw$coefficients[, "coef"]),
  CI_Lower = exp(s_pw$coefficients[, "coef"] - 1.96 * s_pw$coefficients[, "se(coef)"]),
  CI_Upper = exp(s_pw$coefficients[, "coef"] + 1.96 * s_pw$coefficients[, "se(coef)"]),
  P_Value = s_pw$coefficients[, "Pr(>|z|)"]
)
print(pw_res)
write.csv(pw_res, "ITT_Analysis/results/piecewise_cox_results.csv", row.names=FALSE)

cat("Saving Visualization...\n")
layout <- p_hazard | p_rmst
out_name <- "ITT_Analysis/results/figure_delayed_impact.png"
ggsave(out_name, plot = layout, width = 12, height = 5, dpi = 300)
cat("Successfully completed delayed impact sequence!\n")
