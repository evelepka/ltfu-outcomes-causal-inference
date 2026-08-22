library(survival)
df <- read.csv("ITT_Analysis/data/itt_cohort.csv")
df_ltfu <- df[df$itt_group == "Loss to follow-up", ]
sf <- survfit(Surv(time_d, event_d) ~ 1, data = df_ltfu)
print(max(df_ltfu$time_d, na.rm=T))
print(summary(sf, times = c(1, 5, 10, 11)))
