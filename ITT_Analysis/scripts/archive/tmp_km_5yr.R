library(survival)
df <- read.csv("ITT_Analysis/data/itt_cohort.csv")
df_ltfu <- df[df$itt_group == "Loss to follow-up", ]

# 5 year overall
sf <- survfit(Surv(time_d, event_d) ~ 1, data = df_ltfu)
s.all <- summary(sf, times = 5)
cat("Overall 5-year mortality:", (1 - s.all$surv) * 100, "%\n")

# 5 year HIV
sf_hiv <- survfit(Surv(time_d, event_d) ~ hiv_aids, data = df_ltfu)
s.hiv <- summary(sf_hiv, times = 5)
cat("HIV- 5-year:", (1 - s.hiv$surv[1]) * 100, "%, HIV+ 5-year:", (1 - s.hiv$surv[2]) * 100, "%\n")

# 5 year homelessness
sf_homeless <- survfit(Surv(time_d, event_d) ~ homelessness, data = df_ltfu)
s.homeless <- summary(sf_homeless, times = 5)
cat("Not homeless 5-year:", (1 - s.homeless$surv[1]) * 100, "%, Homeless 5-year:", (1 - s.homeless$surv[2]) * 100, "%\n")

# 5 year age
sf_age <- survfit(Surv(time_d, event_d) ~ age_group, data = df_ltfu)
s.age <- summary(sf_age, times = 5)
cat("Age 5-year:", (1 - s.age$surv) * 100, "%\n")
