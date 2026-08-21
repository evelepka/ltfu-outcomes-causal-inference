# 50. Rolling landmark: subgroup absolute risk differences under BOTH windows
# ==============================================================================
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# The standardized subgroup risk differences currently in the manuscript are
# computed over the full 0-24 month window, and in the two highest-mortality
# strata they are not credible: age >=65 gives -4.59 pp and homeless +1.31,
# against CCW's +0.30 and -0.37.
#
# Option B under discussion (2026-08-19) is to report the 6-24 month window as
# the main analysis, as the June submission did. June reported late-window RDs
# for PLHIV, homeless, housed and the youngest group -- all plausible -- but
# never reported one for >=65, so nobody knows whether that stratum behaves in
# the late window or stays negative. That single number decides whether option B
# is complete or whether >=65 has to be omitted from the absolute scale under
# any design.
#
# This script computes both windows side by side so the comparison is one table.
#
# CAVEAT ON THE LATE-WINDOW ESTIMAND -- read before quoting these numbers.
# `prep_outcome(model = "late")` counts events in (0.5, cap] with NO left
# truncation, so it does not condition on surviving six months; everyone is in
# the denominator from time zero. But deaths before six months are censored at
# their own event time, so the quantity is the risk of late death in a world
# where early death does not remove people. Reviewer 3 objected to exactly this
# in their comment on the standardized risk difference: "is it the effect of
# LTFU on death at 6-24 months, had everyone who died during TB treatment
# lived? This would not be a realistic estimand." Whatever these numbers show,
# the estimand needs stating and justifying, not just reporting.
#
# Usage:  Rscript 50_itt_rolling_subgroup_windows.R
#         N_IMP=1 Rscript 50_itt_rolling_subgroup_windows.R    # fast check
#
# Output: ITT_Analysis/results/rolling_subgroup_windows.csv
# ==============================================================================
suppressPackageStartupMessages({ library(splines); library(dplyr) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

HORIZON <- 2
SUBGRPS <- c("age_group", "hiv_aids", "homelessness")
WINDOWS <- c("overall", "late")     # 0-24 mo, and 6-24 mo

# ---------------------------------------------------------------------------
# Standardize to the exposed patients within one subgroup level and one window.
# Same machinery as script 45, with the window as an argument and the subgroup
# variable dropped from the covariate set (both arms come from within the level).
# ---------------------------------------------------------------------------

# Inside a small subgroup, and especially in the late window, some covariate
# levels run out of events entirely; their coefficients then go to infinity and
# coxph warns. Measured on imp_01: `race_clean == "Other"` has ZERO late-window
# events among both PLHIV and 15-24s, and `geo4 == "IntermediarioAdjacente"`
# has one among 15-24s.
#
# Fix is the standard one and is applied per model, not globally: any level with
# fewer than MIN_EVENTS_PER_LEVEL events is merged into the level with the most
# events, which then serves as the reference. Nothing is dropped from the data —
# only the granularity of the adjustment is reduced, and only where the data
# cannot support it. Which variables were collapsed is recorded per row.
MIN_EVENTS_PER_LEVEL <- 5

collapse_sparse_levels <- function(sub, vars) {
  collapsed <- character(0)
  for (v in vars) {
    if (!is.factor(sub[[v]]) && !is.character(sub[[v]])) next
    f  <- as.character(sub[[v]])
    ev <- tapply(sub$event_out, f, sum); ev <- ev[!is.na(ev)]
    if (length(ev) < 2 || all(ev >= MIN_EVENTS_PER_LEVEL)) next
    keep <- names(ev)[ev >= MIN_EVENTS_PER_LEVEL]
    ref  <- names(ev)[which.max(ev)]
    if (!length(keep)) next
    sub[[v]] <- factor(ifelse(f %in% keep, f, ref), levels = union(ref, keep))
    sub[[v]] <- relevel(sub[[v]], ref = ref)
    collapsed <- c(collapsed, v)
  }
  list(d = sub, collapsed = paste(collapsed, collapse = ";"))
}

subgroup_rd <- function(tr, var, level, window, horizon = HORIZON) {
  sub <- if (identical(var, "overall")) tr
         else tr[!is.na(tr[[var]]) & tr[[var]] == level, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  sub <- prep_outcome(sub, window, horizon)
  if (sum(sub$event_out) < 20) return(NULL)

  g <- collapse_sparse_levels(sub, setdiff(COVARS, var)); sub <- g$d

  cv <- drop_constant(sub, setdiff(COVARS, var))
  rhs <- paste(c("expose", cv, "ns(trial_day, df = 3)"), collapse = " + ")
  warned <- character(0)
  fit <- withCallingHandlers(
    tryCatch(coxph(as.formula(paste("Surv(time_out, event_out) ~", rhs)),
                   data = sub, cluster = pid, ties = "efron", x = FALSE, model = TRUE),
             error = function(e) NULL),
    warning = function(w) { warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (is.null(fit)) return(NULL)

  bh <- basehaz(fit, centered = FALSE)
  i  <- findInterval(horizon, bh$time)
  H0 <- if (i == 0) 0 else bh$hazard[max(i, 1)]

  ex <- sub[sub$expose == 1, , drop = FALSE]
  if (!nrow(ex)) return(NULL)
  d1 <- ex; d1$expose <- 1
  d0 <- ex; d0$expose <- 0
  r1 <- mean(1 - exp(-H0 * exp(predict(fit, newdata = d1, type = "lp", reference = "zero"))))
  r0 <- mean(1 - exp(-H0 * exp(predict(fit, newdata = d0, type = "lp", reference = "zero"))))

  data.frame(subgroup = var, level = as.character(level), window = window,
             risk1 = 100 * r1, risk0 = 100 * r0, rd = 100 * (r1 - r0),
             hr = unname(exp(coef(fit)["expose"])),
             n_exposed = nrow(ex), n_events = sum(sub$event_out),
             collapsed = g$collapsed,
             nonconverged = if (any(grepl("may be infinite", warned))) "yes" else "",
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]
cat(sprintf("[50] subgroup RDs, both windows | %d imputation(s)\n", length(imp_files)))

outcome_lk <- build_outcome_lookup()
prepped <- lapply(imp_files, prepare_rolling, outcome_lookup = outcome_lk)
stacks  <- lapply(prepped, build_rolling, comparator = "in_care", carry = character(0))

# The cells to estimate: an overall row first, as the reference against which
# the subgroup rows and the published 0-24 figure (1.375 pp) can be read.
CELLS <- c(list(list(var = "overall", level = "all")),
           unlist(lapply(SUBGRPS, function(v)
             lapply(sort(unique(na.omit(as.character(stacks[[1]][[v]])))),
                    function(l) list(var = v, level = l))), recursive = FALSE))

rows <- list()
for (cell in CELLS) {
  var <- cell$var
  for (lv in cell$level) {
    for (w in WINDOWS) {
      per <- Filter(Negate(is.null), lapply(stacks, subgroup_rd, var = var, level = lv, window = w))
      if (!length(per)) { cat(sprintf("  %-14s %-10s %-8s SKIP\n", var, lv, w)); next }
      b <- bind_rows(per)
      r <- data.frame(subgroup = var, level = lv, window = w,
                      risk1 = mean(b$risk1), risk0 = mean(b$risk0), rd = mean(b$rd),
                      hr = exp(mean(log(b$hr))), n_exposed = round(mean(b$n_exposed)),
                      n_events = round(mean(b$n_events)), M = nrow(b),
                      collapsed = paste(unique(b$collapsed[nzchar(b$collapsed)]), collapse = ";"),
                      nonconverged = if (any(nzchar(b$nonconverged))) "yes" else "",
                      stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- r
      cat(sprintf("  %-14s %-10s %-8s RD %+7.3f  risk1 %6.2f%%  risk0 %6.2f%%  HR %5.2f%s%s\n",
                  var, lv, w, r$rd, r$risk1, r$risk0, r$hr,
                  if (nzchar(r$collapsed)) sprintf("  [collapsed: %s]", r$collapsed) else "",
                  if (nzchar(r$nonconverged)) "  <== NON-CONVERGED" else ""))
    }
  }
}

out <- bind_rows(rows)

# ---------------------------------------------------------------------------
# Cluster bootstrap, same construction as script 45: resample PATIENTS in the
# prepped data and REBUILD the stack, with a distinct seed per replicate and a
# randomly drawn imputation. Restricted to one window (default `late`, the one
# under discussion for reporting) because every extra cell multiplies runtime.
# ---------------------------------------------------------------------------
B <- as.integer(Sys.getenv("B", unset = "0"))
BOOT_WINDOW <- Sys.getenv("BOOT_WINDOW", unset = "late")
if (B > 0) {
  # The five stacked datasets are ~384k rows each and are only needed for the
  # point estimates above; each bootstrap replicate builds its own. Holding them
  # through the loop exhausted RAM and pushed the machine into swap: a first
  # attempt ran at 26% CPU efficiency and stalled after 50 replicates. Free them
  # before starting, and collect between replicates.
  rm(stacks); invisible(gc())
  cat(sprintf("\n  cluster bootstrap on the '%s' window, B=%d\n", BOOT_WINDOW, B))
  set.seed(5050)
  imp_pick <- sample.int(length(prepped), B, replace = TRUE)
  bs <- vector("list", B); t0 <- Sys.time()
  for (b in seq_len(B)) {
    d  <- prepped[[imp_pick[b]]]
    dd <- d[sample.int(nrow(d), nrow(d), replace = TRUE), , drop = FALSE]
    dd$pid <- seq_len(nrow(dd))
    tr <- build_rolling(dd, comparator = "in_care", carry = character(0),
                        seed = SEED + 1000L * b)
    bs[[b]] <- bind_rows(lapply(CELLS, function(cell)
      subgroup_rd(tr, cell$var, cell$level, BOOT_WINDOW)))
    rm(tr, dd); if (b %% 10 == 0) invisible(gc())
    if (b %% 25 == 0) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("    rep %d/%d  %.1fs/rep  ETA %.1f min\n", b, B, el/b, el/b*(B-b)/60))
    }
  }
  bb <- bind_rows(Filter(Negate(is.null), bs))
  ci <- bb |> group_by(subgroup, level) |>
    summarise(rd_lo = quantile(rd, .025, na.rm = TRUE),
              rd_hi = quantile(rd, .975, na.rm = TRUE),
              rd_boot_mean = mean(rd, na.rm = TRUE), B_ok = n(), .groups = "drop") |>
    mutate(window = BOOT_WINDOW)
  out <- left_join(out, ci, by = c("subgroup", "level", "window"))
  cat("\n  bootstrap diagnostic (bias = boot mean - point):\n")
  for (i in which(out$window == BOOT_WINDOW)) {
    bias <- out$rd_boot_mean[i] - out$rd[i]
    cat(sprintf("    %-14s %-10s RD %+6.3f (%+.3f to %+.3f)  bias %+.3f  %s\n",
                out$subgroup[i], out$level[i], out$rd[i], out$rd_lo[i], out$rd_hi[i], bias,
                if (isTRUE(abs(bias) < 0.1 * abs(out$rd[i]))) "acceptable" else "CHECK"))
  }
}

write.csv(out, file.path(ITT_RESULTS_DIR, "rolling_subgroup_windows.csv"), row.names = FALSE)
cat(sprintf("\n[50] wrote %s\n", file.path(ITT_RESULTS_DIR, "rolling_subgroup_windows.csv")))

cat("\n  side by side (RD, percentage points):\n")
w <- out |> select(subgroup, level, window, rd) |>
  tidyr::pivot_wider(names_from = window, values_from = rd)
print(as.data.frame(w), row.names = FALSE, digits = 3)
cat("\n  NOTE: the late-window estimand censors early deaths -- see the header,\n")
cat("  and Reviewer 3's objection to it, before quoting any of these.\n")
