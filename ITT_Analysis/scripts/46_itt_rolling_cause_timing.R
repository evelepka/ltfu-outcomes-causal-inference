# 46. Rolling landmark: cause-specific TIMING curves
# ==============================================================================
# Scripts 42/44 give the all-cause timing curve; 43 gives cause-specific effects
# averaged over disengagement day. This crosses the two: how does the TB-specific
# effect vary with WHEN the patient disengaged, versus the non-TB effect?
#
# WHY IT MATTERS.  The paper's causal claim is that the excess mortality runs
# through untreated tuberculosis. If that is right, the TB curve should carry the
# timing structure -- disengaging with more of the course left should be worse --
# while non-TB mortality, which reflects the social and clinical vulnerability
# that leads people to disengage rather than the missing therapy itself, should be
# comparatively FLAT in disengagement day. A TB curve with structure sitting on a
# flat non-TB curve is much stronger evidence than two elevated averages.
#
# Restricted to DOT patients as the primary read, because script 44 showed the
# timing variable is precisely measured only under DOT: the inferred day is the
# last observed contact, and under self-administration monthly dispensing blurs it
# by up to a month, attenuating any shape. The all-patient version is reported too.
#
# Output: ITT_Analysis/results/rolling_landmark_cause_timing.csv
# ==============================================================================
suppressPackageStartupMessages({ library(splines) })
.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1])))); getwd()
}
source(file.path(.here(), "_paths.R")); source(file.path(.here(), "_rolling.R"))

CAUSES <- c("tb_hybrid", "nontb_hybrid")
GRID   <- seq(DAY_MIN, DAY_MAX, by = 5)

outcome_lk <- build_outcome_lookup()
imp_files <- sort(list.files(ITT_MI_DIR, pattern = "^imp_\\d+\\.csv$", full.names = TRUE))
n_imp <- as.integer(Sys.getenv("N_IMP", unset = length(imp_files)))
imp_files <- imp_files[seq_len(min(n_imp, length(imp_files)))]
cat(sprintf("[46] cause-specific timing | %d imputation(s)\n", length(imp_files)))
lookup <- build_cause_lookup()

stacks <- lapply(imp_files, function(p) {
  d <- prepare_rolling(p, cause_lookup = lookup, outcome_lookup = outcome_lk,
                       extra_factors = "dot_status")
  build_rolling(d, comparator = "in_care", carry = c(CAUSES, "dot_status"))
})

curve_for <- function(subs, cause, covars) {
  tf <- Filter(Negate(is.null),
               lapply(subs, fit_rolling, model = "late", cap = HORIZON_Y,
                      cause = cause, covars = covars, timing = TRUE))
  if (!length(tf)) return(NULL)
  bs <- ns(subs[[1]]$trial_day, df = 3)
  per <- lapply(tf, function(f) {
    b <- coef(f); V <- vcov(f); j <- grep("^expose", names(b))
    X <- cbind(1, predict(bs, GRID))
    list(est = as.vector(X %*% b[j]),
         var = rowSums((X %*% V[j, j, drop = FALSE]) * X))
  })
  E <- do.call(cbind, lapply(per, `[[`, "est"))
  Vr <- do.call(cbind, lapply(per, `[[`, "var"))
  out <- lapply(seq_along(GRID), function(r) {
    q <- pool_loghr(E[r, ], sqrt(Vr[r, ]))
    if (is.null(q)) return(NULL)
    data.frame(day = GRID[r], HR = q$hr, CI_L = q$lo, CI_H = q$hi, N_imp = q$M)
  })
  bind_rows(out)
}

rows <- list()
for (pop in c("DOT", "all")) {
  subs <- if (pop == "DOT")
    lapply(stacks, function(st) st[!is.na(st$dot_status) &
                                    as.character(st$dot_status) == "Yes", ]) else stacks
  covars <- if (pop == "DOT") setdiff(COVARS, "dot_status") else COVARS
  cat(sprintf("\n--- %s (exposed=%s) ---\n", pop,
              format(sum(subs[[1]]$expose), big.mark = ",")))
  for (cs in CAUSES) {
    cv <- curve_for(subs, cs, covars)
    if (is.null(cv)) { cat(sprintf("  %-14s not estimable\n", cs)); next }
    rows[[length(rows) + 1]] <- cv |> mutate(population = pop, cause = cs)
    rng <- max(cv$HR) / min(cv$HR)
    cat(sprintf("  %-14s day1 %.2f  day76 %.2f  day151 %.2f   peak/trough %.2f\n",
                cs, cv$HR[cv$day == 1], cv$HR[cv$day == 76], cv$HR[cv$day == 151], rng))
  }
}
res <- bind_rows(rows)
write.csv(res, file.path(ITT_RESULTS_DIR, "rolling_landmark_cause_timing.csv"),
          row.names = FALSE)

# the discriminating comparison: does TB carry more timing structure than non-TB?
for (pop in c("DOT", "all")) {
  a <- res[res$population == pop & res$cause == "tb_hybrid", ]
  b <- res[res$population == pop & res$cause == "nontb_hybrid", ]
  if (!nrow(a) || !nrow(b)) next
  # Report the ratio, do NOT convert it into a verdict. Peak/trough is a crude
  # summary dominated by the curve endpoints, where the intervals are widest, and
  # any threshold on it would be arbitrary. On this data the TB and non-TB curves
  # have SIMILAR shapes and differ mainly in LEVEL, so the timing structure is not
  # demonstrably TB-specific. See ADR-0004.
  cat(sprintf("\n%s: peak/trough  TB %.2f  vs  non-TB %.2f  (ratio %.2f)\n", pop,
              max(a$HR)/min(a$HR), max(b$HR)/min(b$HR),
              (max(a$HR)/min(a$HR)) / (max(b$HR)/min(b$HR))))
}
cat(sprintf("\n[46] wrote %d rows\n", nrow(res)))
