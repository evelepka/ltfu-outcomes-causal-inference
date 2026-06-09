# 72. Figure 4A: crude time-varying HR(t) for LTFU vs on-treatment, in 2-month
# intervals, from a SINGLE piecewise Cox with an exposure-by-interval interaction
# (counting-process/start-stop, defn-B disengagement). One joint model estimates
# every interval (no per-bin convergence dropouts).
# Output: Figure_4a_HR_over_time_2mo.csv  (consumed by 71_make_fig4_causal_combined.py)
suppressPackageStartupMessages({ library(dplyr); library(survival) })
.here <- function(){a<-grep("^--file=",commandArgs(FALSE),value=TRUE)
  if(length(a)) return(dirname(normalizePath(sub("^--file=","",a[1])))); getwd()}
source(file.path(.here(),"_paths.R"))
SHIFT <- 30/365.25; HORIZON <- 2.0; INTERVAL_MO <- 2

df <- read.csv(COHORT_CSV, stringsAsFactors=FALSE); df$patient_id <- seq_len(nrow(df))
df$tx_dur <- as.numeric(difftime(as.Date(df$end_date), as.Date(df$best_start), units="days"))/365.25
df$tx_dur_true <- pmax(df$tx_dur - SHIFT, 1/365.25)
df$event_d  <- ifelse(df$time_d_tx > HORIZON, 0, df$event_d)
df$time_d_tx <- ifelse(df$time_d_tx > HORIZON, HORIZON, df$time_d_tx)
# counting-process exposure transition at the (defn-B) disengagement time
p1 <- df; p1$tstart <- 0
p1$tstop <- ifelse(p1$itt_group=="Loss to follow-up" & p1$tx_dur_true < p1$time_d_tx, p1$tx_dur_true, p1$time_d_tx)
p1$event <- ifelse(p1$itt_group=="Loss to follow-up" & p1$tx_dur_true < p1$time_d_tx, 0, p1$event_d)
p1$expose <- 0
p2 <- df[df$itt_group=="Loss to follow-up" & df$tx_dur_true < df$time_d_tx, ]
p2$tstart <- p2$tx_dur_true; p2$tstop <- p2$time_d_tx; p2$event <- p2$event_d; p2$expose <- 1
ds <- bind_rows(p1, p2) |> dplyr::filter(round(tstop - tstart, 4) > 0)

w <- INTERVAL_MO/12; BR <- seq(0, HORIZON, by=w)
dp <- survSplit(Surv(tstart,tstop,event) ~ ., data=ds, cut=BR[-c(1,length(BR))], episode="bin")
dp <- dp[(dp$tstop - dp$tstart) > 1e-7, ]
dp$bin <- factor(dp$bin, levels=seq_len(length(BR)-1))
fit <- coxph(Surv(tstart,tstop,event) ~ expose:bin, data=dp, ties="breslow")
ci <- summary(fit)$conf.int
rows <- list()
for (i in seq_len(length(BR)-1)) {
  nm <- paste0("expose:bin", i)
  if (nm %in% rownames(ci) && is.finite(ci[nm,"exp(coef)"]))
    rows[[i]] <- data.frame(month_mid=(i-0.5)*w, HR=ci[nm,"exp(coef)"],
                            CI_L=ci[nm,"lower .95"], CI_H=ci[nm,"upper .95"])
}
out <- bind_rows(rows); out <- out[is.finite(out$CI_H), ]
f <- file.path(ITT_RESULTS_DIR, "Figure_4a_HR_over_time_2mo.csv")
write.csv(out, f, row.names=FALSE)
cat(sprintf("[72] wrote %s (%d intervals)\n", f, nrow(out))); print(out)
