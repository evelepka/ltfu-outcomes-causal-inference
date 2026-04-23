library(dplyr)
library(survival)
library(ggplot2)
library(patchwork)
library(broom)

ts <- as.character(as.integer(Sys.time()))

df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)
df <- df %>% filter(time_d_tx >= 0)

max_yrs <- 2.0

# Crude
df_crude <- df %>% mutate(
  event_cap = ifelse(time_d_tx > max_yrs, 0, event_d),
  time_cap = pmin(time_d_tx, max_yrs)
)
fit_crude <- survfit(Surv(time_cap, event_cap) ~ itt_group, data = df_crude)
td_crude <- tidy(fit_crude)
td_crude$strata <- gsub("itt_group=", "", td_crude$strata)

p_crude <- ggplot(td_crude, aes(x = time, y = 1 - estimate, color = strata)) +
  geom_step(linewidth = 1) +
  theme_classic() +
  scale_color_manual(values = c("Loss to follow-up" = "#d95f02", "Non-LTFU" = "#1b9e77")) +
  labs(title = "Panel A: Crude 2-Year Mortality", subtitle = "Anchored at Tx Day 0",
       x = "Years Since Treatment Start", y = "Cumulative Mortality") +
  theme(legend.position = "bottom", legend.title=element_blank(), plot.title = element_text(face="bold", size=10)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.15)) +
  scale_x_continuous(limits = c(0, max_yrs))

# Target Trial Function
generate_tt_plot <- function(data, m) {
  start_yrs <- ((m - 1) * 30) / 365.25
  end_yrs   <- (m * 30) / 365.25
  
  df_tt <- data %>% filter(time_d_tx > start_yrs) %>%
    mutate(
      tx_duration_yrs = as.numeric(difftime(as.Date(end_date), as.Date(best_start), units="days")) / 365.25,
      eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
    ) %>%
    filter(eligible == 1) %>%
    mutate(
      expose_label = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs >= start_yrs & tx_duration_yrs < end_yrs, "Abandoned", "Maintained"),
      time_followup = time_d_tx - start_yrs
    ) %>%
    mutate(
      event_cap = ifelse(time_followup > max_yrs, 0, event_d),
      time_cap = pmin(time_followup, max_yrs),
      is_abandon = ifelse(expose_label == "Abandoned", 1, 0)
    ) %>% filter(time_cap > 0)
    
  df_split <- survSplit(Surv(time_cap, event_cap) ~ is_abandon, data = df_tt, cut = c(0.5), episode = "Epoch")
  fit_cox <- coxph(Surv(tstart, time_cap, event_cap) ~ is_abandon:strata(Epoch), data = df_split)
  s <- summary(fit_cox)
  
  e1_hr <- exp(s$coefficients["is_abandon:strata(Epoch)Epoch=1", "coef"])
  e1_l  <- exp(s$coefficients["is_abandon:strata(Epoch)Epoch=1", "coef"] - 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=1", "se(coef)"])
  e1_u  <- exp(s$coefficients["is_abandon:strata(Epoch)Epoch=1", "coef"] + 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=1", "se(coef)"])
  
  e2_hr <- exp(s$coefficients["is_abandon:strata(Epoch)Epoch=2", "coef"])
  e2_l  <- exp(s$coefficients["is_abandon:strata(Epoch)Epoch=2", "coef"] - 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=2", "se(coef)"])
  e2_u  <- exp(s$coefficients["is_abandon:strata(Epoch)Epoch=2", "coef"] + 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=2", "se(coef)"])
  
  anno_text_early <- sprintf("Early aHR: %.2f\n(%.2f-%.2f)", e1_hr, e1_l, e1_u)
  anno_text_late  <- sprintf("Late aHR: %.2f\n(%.2f-%.2f)", e2_hr, e2_l, e2_u)

  fit <- survfit(Surv(time_cap, event_cap) ~ expose_label, data = df_tt)
  td <- tidy(fit)
  td$strata <- gsub("expose_label=", "", td$strata)
  
  p <- ggplot(td, aes(x = time, y = 1 - estimate, color = strata)) +
    geom_vline(xintercept = 0.5, linetype="dashed", color="gray50") +
    geom_step(linewidth = 1) +
    annotate("text", x = 0.45, y = 0.14, label = anno_text_early, hjust = 1, vjust = 1, size = 2.5) +
    annotate("text", x = 0.55, y = 0.14, label = anno_text_late, hjust = 0, vjust = 1, size = 2.5) +
    theme_classic() +
    scale_color_manual(values = c("Abandoned" = "#d95f02", "Maintained" = "#1b9e77")) +
    labs(title = sprintf("Trial %d: Abandonment M%d", m, m),
         subtitle = sprintf("Baseline T=0 shifts to Day %d", (m-1)*30),
         x = "Years since T=0", y = "Cum. Mort.") +
    theme(legend.position = "none", plot.title = element_text(face="bold", size=9), axis.title.y = element_blank()) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 0.15)) +
    scale_x_continuous(limits = c(0, max_yrs))
  return(p)
}

p2 <- generate_tt_plot(df, 2)
p3 <- generate_tt_plot(df, 3)
p4 <- generate_tt_plot(df, 4)
p5 <- generate_tt_plot(df, 5)

layout <- (p_crude / (p2 | p3) / (p4 | p5)) + plot_layout(heights = c(1.5, 1, 1))
out_name <- sprintf("ITT_Analysis/results/figure2_crude_panels_2yr_%s.png", ts)
ggsave(out_name, plot = layout, width = 10, height = 10, dpi = 300)
cat(out_name)
