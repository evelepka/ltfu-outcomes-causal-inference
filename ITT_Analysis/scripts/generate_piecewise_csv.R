df <- read.csv("ITT_Analysis/data/itt_cohort.csv")
library(survival)

res_list <- list()
for (m in 1:6) {
  start_yrs <- ((m - 1) * 30) / 365.25
  end_yrs   <- (m * 30) / 365.25
  
  df_tt <- df[df$time_d_tx > start_yrs, ]
  df_tt$tx_duration_yrs <- as.numeric(difftime(as.Date(df_tt$end_date), as.Date(df_tt$best_start), units="days")) / 365.25
  df_tt$eligible <- ifelse(df_tt$itt_group == "Non-LTFU" | df_tt$tx_duration_yrs >= start_yrs, 1, 0)
  df_tt <- df_tt[df_tt$eligible == 1, ]
  df_tt$is_abandon <- ifelse(df_tt$itt_group == "Loss to follow-up" & df_tt$tx_duration_yrs >= start_yrs & df_tt$tx_duration_yrs < end_yrs, 1, 0)
  df_tt$time_followup <- df_tt$time_d_tx - start_yrs
  df_tt$event_cap <- ifelse(df_tt$time_followup > 5, 0, df_tt$event_d)
  df_tt$time_cap <- pmin(df_tt$time_followup, 5)
  df_tt <- df_tt[df_tt$time_cap > 0, ]
  
  df_split <- survSplit(Surv(time_cap, event_cap) ~ is_abandon, data = df_tt, cut = c(0.5), episode = "Epoch")
  fit <- coxph(Surv(tstart, time_cap, event_cap) ~ is_abandon:strata(Epoch), data = df_split)
  s <- summary(fit)
  
  # Epoch 1
  res_list[[length(res_list) + 1]] <- data.frame(
    Trial = paste("Month", m), Epoch = "Early (Months 0-6)",
    HR = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=1", "coef"]),
    CI_Lower = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=1", "coef"] - 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=1", "se(coef)"]),
    CI_Upper = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=1", "coef"] + 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=1", "se(coef)"])
  )
  # Epoch 2
  res_list[[length(res_list) + 1]] <- data.frame(
    Trial = paste("Month", m), Epoch = "Late (Months 6-60)",
    HR = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=2", "coef"]),
    CI_Lower = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=2", "coef"] - 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=2", "se(coef)"]),
    CI_Upper = exp(s$coefficients["is_abandon:strata(Epoch)Epoch=2", "coef"] + 1.96 * s$coefficients["is_abandon:strata(Epoch)Epoch=2", "se(coef)"])
  )
}
df_res <- do.call(rbind, res_list)
write.csv(df_res, "ITT_Analysis/results/target_trial_piecewise_hr.csv", row.names=FALSE)
