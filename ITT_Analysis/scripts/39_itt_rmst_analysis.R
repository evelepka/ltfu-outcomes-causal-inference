library(dplyr)
library(survival)
library(survRM2)

cat("Loading ITT cohort for RMST analysis...\n")
df <- read.csv("ITT_Analysis/data/itt_cohort.csv", stringsAsFactors = FALSE)
df <- df %>% filter(time_d_tx >= 0)

max_yrs <- 5

# Overall Crude RMST
df_crude <- df %>% mutate(
  event_cap = ifelse(time_d_tx > max_yrs, 0, event_d),
  time_cap = pmin(time_d_tx, max_yrs),
  arm = ifelse(itt_group == "Non-LTFU", 1, 0) # 1 = Maintained Care, 0 = Abandoned
)

# RMST requires time in days or years, we use years.
rmst_obj <- rmst2(time = df_crude$time_cap, 
                  status = df_crude$event_cap, 
                  arm = df_crude$arm, 
                  tau = max_yrs)
cat("\n============================================\n")
cat("CRUDE 5-YEAR RMST (Anchored at Tx Start Day 0)\n")
cat("============================================\n")
print(rmst_obj$unadjusted.result)

# Target Trial Arrays Function
generate_tt_rmst <- function(data, m, tau_yrs=5) {
  start_yrs <- ((m - 1) * 30) / 365.25
  end_yrs   <- (m * 30) / 365.25
  
  df_tt <- data %>% filter(time_d_tx > start_yrs) %>%
    mutate(
      tx_duration_yrs = as.numeric(difftime(as.Date(end_date), as.Date(best_start), units="days")) / 365.25,
      eligible = ifelse(itt_group == "Non-LTFU" | tx_duration_yrs >= start_yrs, 1, 0)
    ) %>%
    filter(eligible == 1) %>%
    mutate(
      arm = ifelse(itt_group == "Loss to follow-up" & tx_duration_yrs >= start_yrs & tx_duration_yrs < end_yrs, 0, 1),
      time_followup = time_d_tx - start_yrs
    ) %>%
    mutate(
      event_cap = ifelse(time_followup > tau_yrs, 0, event_d),
      time_cap = pmin(time_followup, tau_yrs)
    )
    
  rmst_tt <- rmst2(time = df_tt$time_cap, 
                   status = df_tt$event_cap, 
                   arm = df_tt$arm, 
                   tau = tau_yrs)
  
  diff_res <- rmst_tt$unadjusted.result[1,] # RMST(arm=1) - RMST(arm=0)
  
  return(data.frame(
    Trial_Month = m,
    Diff_Years = diff_res[1],
    Diff_Lower = diff_res[2],
    Diff_Upper = diff_res[3],
    P_Value = diff_res[4]
  ))
}

cat("\n============================================\n")
cat("SEQUENTIAL TARGET TRIAL 5-YEAR RMST \n")
cat("============================================\n")
res_list <- list()
for (m in 1:6) {
  res_list[[m]] <- generate_tt_rmst(df, m, 5)
}
final_res <- do.call(rbind, res_list)
print(final_res)

write.csv(final_res, "ITT_Analysis/results/target_trial_rmst.csv", row.names=FALSE)
cat("Exported RMST object successfully.\n")
