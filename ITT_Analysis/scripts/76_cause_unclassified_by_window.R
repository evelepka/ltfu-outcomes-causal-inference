# ==============================================================================
# 76. What share of deaths could not be assigned to a cause class, by arm and
#     by follow-up window?
#
# WHY THIS EXISTS
# The Discussion limitation used to read "Cause of death could not be classified
# for 11.4% of late-window deaths among patients lost to follow-up and 16.1%
# among those continuing in care." Those shares were computed on the LATE window
# (0.5-2 y from the trial origin), which the paper has removed. Quoting them
# beside five-year results would be wrong, so they are recomputed here.
#
# "Could not be classified" uses the SAME definition as script 49's ascertainment
# table: !tb_hybrid & !nontb_hybrid. Note this is the narrow sense. Under the
# two-class partition the paper reports, the residual class absorbs these deaths,
# so they are not missing from the decomposition -- they are deaths whose cause is
# not specific enough to place in the tuberculosis class.
#
# The script FIRST reproduces the old late-window figures. If it cannot, its
# five-year figures should not be trusted either, so both are always printed.
#
# TB_ANY_LINE is left OFF, per the instruction in _rolling.R: "Sensitivity
# analysis only. Leave OFF for anything the manuscript reports."
#
# Usage:  Rscript ITT_Analysis/scripts/76_cause_unclassified_by_window.R
# Output: ITT_Analysis/results/rolling_cause_unclassified.csv
# ==============================================================================
suppressPackageStartupMessages({ library(dplyr) })

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R"))
source(file.path(.here(), "_rolling.R"))

if (nzchar(Sys.getenv("TB_ANY_LINE")))
  stop("TB_ANY_LINE must be OFF for a manuscript number (see _rolling.R)")

lookup <- build_cause_lookup(verbose = FALSE)

imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$",
                             full.names = TRUE))
stopifnot(length(imp_files) > 0)

# Imputation 1 only, matching script 49, which computes its ascertainment table
# from stacks[[1]]. Cause is not imputed, so pooling would not change these
# shares materially; using the same input makes the reproduction check exact.
prepped <- prepare_rolling(imp_files[1], cause_lookup = lookup)
s1 <- build_rolling(prepped, comparator = "in_care",
                    carry = c("tb_hybrid", "nontb_hybrid"))

WINDOWS <- list(
  list(label = "late (0.5-2 y, as previously reported)", lo = 0.5,  hi = 2),
  list(label = "whole 5 y from origin (current paper)",  lo = 0.0,  hi = 5),
  list(label = "whole 2 y from origin",                  lo = 0.0,  hi = 2)
)

rows <- list()
for (w in WINDOWS) {
  sub <- s1[s1$event_d_num == 1 & s1$time_raw > w$lo & s1$time_raw <= w$hi, ,
            drop = FALSE]
  sub$unclassified <- !sub$tb_hybrid & !sub$nontb_hybrid
  agg <- sub |>
    group_by(expose) |>
    summarise(n_deaths = n(),
              n_unclassified = sum(unclassified),
              pct_unclassified = 100 * mean(unclassified),
              n_tb = sum(tb_hybrid),
              n_nontb = sum(nontb_hybrid),
              .groups = "drop")
  agg$arm    <- ifelse(agg$expose == 1, "LTFU", "in care")
  agg$window <- w$label
  agg$lo <- w$lo; agg$hi <- w$hi
  rows[[length(rows) + 1]] <- agg

  cat(sprintf("\n  %s\n", w$label))
  for (i in seq_len(nrow(agg))) {
    cat(sprintf("    %-8s  n=%6d   unclassified %5d (%.1f%%)   tb %5d   nontb %5d\n",
                agg$arm[i], agg$n_deaths[i], agg$n_unclassified[i],
                agg$pct_unclassified[i], agg$n_tb[i], agg$n_nontb[i]))
  }
}

out <- bind_rows(rows)[, c("window", "lo", "hi", "arm", "expose", "n_deaths",
                           "n_unclassified", "pct_unclassified", "n_tb",
                           "n_nontb")]
write.csv(out, file.path(ITT_RESULTS_DIR, "rolling_cause_unclassified.csv"),
          row.names = FALSE)
cat(sprintf("\n  wrote %s\n",
            file.path(ITT_RESULTS_DIR, "rolling_cause_unclassified.csv")))

lw <- out[out$lo == 0.5, ]
cat(sprintf("\n  REPRODUCTION CHECK vs the old manuscript sentence (11.4%% LTFU, 16.1%% in care):\n"))
cat(sprintf("    LTFU    %.1f%%   in care %.1f%%\n",
            lw$pct_unclassified[lw$arm == "LTFU"],
            lw$pct_unclassified[lw$arm == "in care"]))
