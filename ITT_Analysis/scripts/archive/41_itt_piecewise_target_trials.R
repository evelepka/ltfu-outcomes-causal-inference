library(dplyr)
library(survival)
library(ggplot2)
library(patchwork)

ts <- as.character(as.integer(Sys.time()))

df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)
df <- df %>% filter(time_d_tx >= 0)

max_yrs <- 5.0

# PANEL A: Instantaneous Hazard Ratio over time (Piecewise Spline)
df_crude_split <- df %>% filter(itt_group %in% c("Loss to follow-up", "Non-LTFU")) %>%
  mutate(event_cap = ifelse(time_d_tx > max_yrs, 0, event_d),
         time_cap = pmin(time_d_tx, max_yrs),
         is_abandon = ifelse(itt_group == "Loss to follow-up", 1, 0)) %>%
  filter(time_cap > 0)

cuts <- seq(0.25, 4.75, by=0.25)
df_split <- survSplit(Surv(time_cap, event_cap) ~ is_abandon, data = df_crude_split, cut = cuts, episode = "Epoch")
fit_a <- coxph(Surv(tstart, time_cap, event_cap) ~ is_abandon:strata(Epoch), data = df_split)
s_a <- summary(fit_a)

df_hr <- data.frame(
  time = c(cuts * 12, 60), # Map exactly to Months natively!
  hr = exp(s_a$coefficients[, "coef"]),
  ci_l = exp(s_a$coefficients[, "coef"] - 1.96 * s_a$coefficients[, "se(coef)"]),
  ci_u = exp(s_a$coefficients[, "coef"] + 1.96 * s_a$coefficients[, "se(coef)"])
)

pA <- ggplot(df_hr, aes(x = time, y = hr)) +
  geom_hline(yintercept = 1, linetype = "solid", color = "gray20", linewidth=0.8) +
  geom_ribbon(aes(ymin = ci_l, ymax = ci_u), fill = "#8e44ad", alpha = 0.2) +
  geom_line(color = "#8e44ad", linewidth = 1.2) +
  geom_point(color = "#8e44ad", size=1.5) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = "gray85", linetype = "dashed")) +
  scale_y_log10(limits = c(0.01, 20), breaks = c(0.05, 0.1, 0.5, 1, 2, 5, 10)) +
  scale_x_continuous(breaks = seq(0, 60, by=12), minor_breaks=seq(0, 60, by=3)) +
  labs(title = "Panel A: Dynamic Hazard Ratio over time (3-Month Granular Epochs)",
       subtitle = "Calculated via continuous Piecewise Cox mapping across crude timeline",
       x = "Months Since Treatment Start", y = "Hazard Ratio (log scale)")


# PANELS B & C: Piecewise Epochs for Target Trials
split_yrs <- 0.5 
res_list <- list()
for (m in 1:6) {
  start_yrs <- ((m - 1) * 30) / 365.25
  end_yrs   <- (m * 30) / 365.25
  
  df_tt <- df %>% filter(time_d_tx > start_yrs) %>%
    mutate(
      tx_duration_yrs = as.numeric(difftime(as.Date(end_date), as.Date(best_start), units="days")) / 365.25,
      eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
    ) %>%
    filter(eligible == 1) %>%
    mutate(
      is_abandon = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs >= start_yrs & tx_duration_yrs < end_yrs, 1, 0),
      time_followup = time_d_tx - start_yrs
    ) %>%
    mutate(
      event_cap = ifelse(time_followup > max_yrs, 0, event_d),
      time_cap = pmin(time_followup, max_yrs)
    ) %>% filter(time_cap > 0)
    
  df_split_tt <- survSplit(Surv(time_cap, event_cap) ~ is_abandon, data = df_tt, cut = c(split_yrs), episode = "Epoch")
  fit <- coxph(Surv(tstart, time_cap, event_cap) ~ is_abandon:strata(Epoch), data = df_split_tt)
  s <- summary(fit)
  
  res_list[[length(res_list) + 1]] <- data.frame(
    Trial = paste("Month", m), Epoch = "Early post-LTFU",
    HR = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=1", "coef"]),
    CI_Lower = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=1", "coef"] - 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=1", "se(coef)"]),
    CI_Upper = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=1", "coef"] + 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=1", "se(coef)"])
  )
  res_list[[length(res_list) + 1]] <- data.frame(
    Trial = paste("Month", m), Epoch = "Late post-LTFU",
    HR = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=2", "coef"]),
    CI_Lower = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=2", "coef"] - 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=2", "se(coef)"]),
    CI_Upper = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=2", "coef"] + 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=2", "se(coef)"])
  )
}

df_res <- do.call(rbind, res_list)
df_res$Trial <- factor(df_res$Trial, levels = rev(paste("Month", 1:6)))

p_early <- ggplot(df_res %>% filter(Epoch == "Early post-LTFU"), aes(x = HR, y = Trial)) +
  geom_vline(xintercept = 1, linetype = "solid", color = "gray20", linewidth=0.8) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper), width = 0.2, color = "#2c3e50") +
  geom_point(size = 4, color = "#1b9e77") +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = "gray85", linetype = "dashed")) +
  scale_x_log10(limits=c(0.1, 5.0), breaks=c(0.1, 0.5, 1, 2, 5)) +
  labs(title = "Panel B: Early post-LTFU (months 0-6)",
       x = "Hazard Ratio (log scale)", y = "")

p_late <- ggplot(df_res %>% filter(Epoch == "Late post-LTFU"), aes(x = HR, y = Trial)) +
  geom_vline(xintercept = 1, linetype = "solid", color = "gray20", linewidth=0.8) +
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper), width = 0.2, color = "#2c3e50") +
  geom_point(size = 4, color = "#d95f02") +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = "gray85", linetype = "dashed"),
        axis.text.y = element_blank()) +
  scale_x_log10(limits=c(0.5, 10), breaks=c(0.5, 1, 2, 5, 10)) +
  labs(title = "Panel C: Late post-LTFU (months 6-60)",
       x = "Hazard Ratio (log scale)", y = "")

layout <- pA / (p_early | p_late)
out_name <- sprintf("ITT_Analysis/results/figure3_piecewise_%s.png", ts)
ggsave(out_name, plot = layout, width = 10, height = 8, dpi = 300)
cat(out_name)
