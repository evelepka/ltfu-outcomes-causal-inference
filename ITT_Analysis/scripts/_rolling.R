# Shared helpers for the rolling (continuous-time) landmark analysis.
# ==============================================================================
# Sourced by 42_itt_rolling_landmark.R (all-cause + timing) and
# 43_itt_rolling_cause_subgroup.R (cause-specific + subgroups).
#
# WHY A SHARED FILE: the cause-attribution defect fixed on 2026-08-16 existed in
# TWO places at once, because 30e and 30i are near-duplicates and only one was
# ever read. Anything used by more than one rolling script belongs here, so a fix
# lands once.
# ==============================================================================

suppressPackageStartupMessages({ library(dplyr); library(survival) })

GRACE_D   <- 30      # days of absence required to meet the LTFU definition
DAY_MIN   <- 1
DAY_MAX   <- 180     # end of the intended 6-month course
K_COMP    <- 20      # comparators sampled per exposed patient per trial
HORIZON_Y <- 2
SEED      <- 2026

GEO_REF <- "Urbano"  # reference level for geo4; custody folds in here

# `geo4` is a MODEL covariate as of 2026-08-18 (owner decision), addressing
# Reviewer 1 comment 11 on spatial/administrative clustering. Adjustment
# attenuates the primary estimate by a few percent; the magnitude is recorded in
# docs/, not here, because this file is mirrored to a PUBLIC repository.
COVARS <- c("age_group", "sex", "race_clean", "edu_clean", "hiv_aids",
            "diabetes", "alcohol", "drug_use", "incarcerated",
            "homelessness", "hosp_admission", "clinical_clean", "dot_status",
            "geo4")

# ---------------------------------------------------------------------------
# Cause-of-death attribution
# ---------------------------------------------------------------------------
# Mirrors classify_cod() in 30e/30i INCLUDING the 2026-08-16 fix: the TBweb
# Obito outcome is taken from ANY episode, not just the index one. An LTFU
# patient's index episode closes as `Abandono`, so an index-only lookup
# discarded the cause for 1,058 of 1,668 LTFU deaths (63.4%) recorded on
# retreatment episodes.
#
# SIM ICD-10 takes precedence; TBweb only fills gaps. Respiratory (^J) and
# non-B200 HIV (^B20-B24) deaths belong to NEITHER hybrid class and are
# therefore censored in both cause-specific analyses, as is unknown cause.
build_cause_lookup <- function(verbose = TRUE) {
  raw_path <- file.path(DATA_DIR, "Final_table_cleaned.csv")
  raw <- read.csv(raw_path, stringsAsFactors = FALSE)
  parse_dt <- function(x) as.Date(x, format = "%B %d, %Y")
  raw$end_date <- parse_dt(raw$end_date)
  raw$dod      <- parse_dt(raw$dod)

  TRANSFER <- c("Transf Outro Municipio", "Transf Outro Estado/Pais")
  novo <- raw[trimws(tools::toTitleCase(tolower(raw$case_type))) == "Novo", ]
  novo <- novo[!is.na(novo$case_outcome) & nzchar(trimws(novo$case_outcome)) &
                 novo$case_outcome != "Mud Diag" &
                 !novo$case_outcome %in% TRANSFER, ]
  novo  <- novo[order(novo$end_date), ]
  first <- novo[!duplicated(novo$sinan_clean), c("sinan_clean", "case_outcome")]

  # --- the 2026-08-16 fix: Obito outcome from ANY episode ------------------
  obito <- raw[!is.na(raw$case_outcome) &
                 trimws(raw$case_outcome) %in% c("Obito TB", "Obito NTB"), ]
  obito <- obito[order(obito$end_date), ]
  obito <- obito[!duplicated(obito$sinan_clean, fromLast = TRUE),
                 c("sinan_clean", "case_outcome")]
  names(obito)[2] <- "obito_outcome"
  first <- merge(first, obito, by = "sinan_clean", all = TRUE)
  first$case_outcome <- ifelse(!is.na(first$obito_outcome),
                               first$obito_outcome, first$case_outcome)
  first$obito_outcome <- NULL
  if (verbose) cat(sprintf("  [cause] Obito recovered from any episode: %d\n",
                           nrow(obito)))

  dr <- raw[!is.na(raw$dod) & !is.na(raw$cause_of_death_code) &
              nzchar(raw$cause_of_death_code), ]
  dr <- dr[order(dr$dod), ]
  dr <- dr[!duplicated(dr$sinan_clean, fromLast = TRUE),
           c("sinan_clean", "cause_of_death_code")]
  dr$cause_of_death_code <- toupper(trimws(dr$cause_of_death_code))

  a <- merge(first, dr, by = "sinan_clean", all = TRUE)
  cod <- ifelse(is.na(a$cause_of_death_code), "", a$cause_of_death_code)
  known     <- nzchar(cod)
  tb_strict <- known & grepl("^(A1[5-9]|B90|B200)", cod)
  resp      <- known & grepl("^J[0-9]", cod)
  hiv_other <- known & grepl("^B2[0-4]", cod) & !grepl("^B200", cod)

  cls <- rep("unknown", nrow(a))
  cls[known & !tb_strict & !resp & !hiv_other] <- "non_tb"
  cls[hiv_other] <- "hiv_other"; cls[resp] <- "respiratory"
  cls[tb_strict] <- "tb_strict"
  oc <- ifelse(is.na(a$case_outcome), "", a$case_outcome)
  cls[cls == "unknown" & oc == "Obito TB"]  <- "tb_via_tbweb"
  cls[cls == "unknown" & oc == "Obito NTB"] <- "ntb_via_tbweb"

  data.frame(sinan_clean = a$sinan_clean,
             tb_hybrid    = cls %in% c("tb_strict", "tb_via_tbweb"),
             nontb_hybrid = cls %in% c("non_tb", "ntb_via_tbweb"),
             tb_broad     = cls %in% c("tb_strict", "tb_via_tbweb",
                                       "respiratory", "hiv_other"),
             tb_simonly    = cls == "tb_strict",
             nontb_simonly = cls == "non_tb",
             stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Index-episode outcome lookup (light: 4 columns only)
# ---------------------------------------------------------------------------
# `Abandono Primario` (never established on treatment) must be identified from
# the TBweb OUTCOME CODE, not inferred from dates. Only 377 LTFU patients are
# Abandono Primario; a date-based rule (`back-shifted day < 1`) captures 1,530 and
# wrongly discards 1,169 genuine `Abandono` patients who DID start treatment and
# stopped almost immediately. Those belong in the day-1 bin.
build_outcome_lookup <- function(verbose = TRUE) {
  raw <- read.csv(file.path(DATA_DIR, "Final_table_cleaned.csv"),
                  stringsAsFactors = FALSE)[, c("sinan_clean", "case_type",
                                                "case_outcome", "end_date")]
  raw$end_date <- as.Date(raw$end_date, format = "%B %d, %Y")
  TRANSFER <- c("Transf Outro Municipio", "Transf Outro Estado/Pais")
  novo <- raw[trimws(tools::toTitleCase(tolower(raw$case_type))) == "Novo", ]
  novo <- novo[!is.na(novo$case_outcome) & nzchar(trimws(novo$case_outcome)) &
                 novo$case_outcome != "Mud Diag" &
                 !novo$case_outcome %in% TRANSFER, ]
  novo <- novo[order(novo$end_date), ]
  first <- novo[!duplicated(novo$sinan_clean), c("sinan_clean", "case_outcome")]
  first$primary_aband <- trimws(first$case_outcome) == "Abandono Primario"
  if (verbose) cat(sprintf("  [outcome] Abandono Primario: %d individuals\n",
                           sum(first$primary_aband)))
  first[, c("sinan_clean", "primary_aband")]
}

# ---------------------------------------------------------------------------
# Geography (Reviewer 1 comment 11)
# ---------------------------------------------------------------------------
# Municipality -> IBGE 2017 rural/urban typology, via the crosswalk built by
# 47_build_ibge_typology.py. That script keys on the RAW TBweb `tx_city` string
# for every distinct value present in the data, so this function only ever does
# an exact match -- accent/case normalisation is locale-dependent in R and is
# deliberately NOT repeated here.
#
# `geo4` folds custody into the reference level: `geo_class == "Prison"` is 100%
# collinear with the `incarcerated` covariate (all 20,669 custody patients are
# incarcerated; 99.97% of incarcerated are custody). Giving it its own level
# triggers a non-convergence warning and contaminates the exposure estimate.
#
# The city is taken from the INDEX episode. That is correct here and does NOT hit
# invariant 8: `tx_city` is recorded on every episode including the index, so
# there is no differential missingness (0 missing in either arm). Using the modal
# or most-recent city WOULD hit it -- 9.86% of LTFU patients have >1 city versus
# 1.35% of non-LTFU, because LTFU patients are the ones with retreatment records.
.rolling_cache <- new.env(parent = emptyenv())

build_geo_lookup <- function(verbose = TRUE) {
  if (!is.null(.rolling_cache$geo)) return(.rolling_cache$geo)
  xw_path <- file.path(PROJECT_ROOT, "ITT_Analysis", "external",
                       "municipality_typology_sp.csv")
  if (!file.exists(xw_path))
    stop("[geo] missing crosswalk: run 47_build_ibge_typology.py first")
  xw <- read.csv(xw_path, stringsAsFactors = FALSE)[, c("municipality", "geo4")]
  names(xw)[1] <- "tx_city"   # crosswalk column is `municipality`

  raw <- read.csv(file.path(DATA_DIR, "Final_table_cleaned.csv"),
                  stringsAsFactors = FALSE)[, c("sinan_clean", "case_type",
                                                "case_outcome", "end_date",
                                                "tx_city")]
  raw$end_date <- as.Date(raw$end_date, format = "%B %d, %Y")
  TRANSFER <- c("Transf Outro Municipio", "Transf Outro Estado/Pais")
  novo <- raw[trimws(tools::toTitleCase(tolower(raw$case_type))) == "Novo", ]
  novo <- novo[!is.na(novo$case_outcome) & nzchar(trimws(novo$case_outcome)) &
                 novo$case_outcome != "Mud Diag" &
                 !novo$case_outcome %in% TRANSFER, ]
  novo <- novo[order(novo$end_date), ]
  first <- novo[!duplicated(novo$sinan_clean), c("sinan_clean", "tx_city")]
  # fall back to any episode when the index episode has no city recorded
  anyc <- raw[order(raw$end_date), ]
  anyc <- anyc[!duplicated(anyc$sinan_clean), c("sinan_clean", "tx_city")]
  names(anyc)[2] <- "tx_city_any"
  first <- merge(first, anyc, by = "sinan_clean", all = TRUE)
  first$tx_city <- ifelse(is.na(first$tx_city), first$tx_city_any, first$tx_city)

  out <- merge(first[, c("sinan_clean", "tx_city")], xw, by = "tx_city",
               all.x = TRUE)
  cov <- mean(!is.na(out$geo4))
  if (cov < 0.999)
    stop(sprintf(paste0("[geo] crosswalk covers only %.2f%% of patients; ",
                        "re-run 47_build_ibge_typology.py"), 100 * cov))
  if (verbose) cat(sprintf("  [geo] %d patients mapped (%.2f%%)\n",
                           sum(!is.na(out$geo4)), 100 * cov))
  out <- out[, c("sinan_clean", "geo4")]
  .rolling_cache$geo <- out
  out
}

# ---------------------------------------------------------------------------
# Load one imputed dataset; derive day-level clocks
# ---------------------------------------------------------------------------
prepare_rolling <- function(path, cause_lookup = NULL, extra_factors = NULL,
                            outcome_lookup = NULL) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date_start <- as.Date(d$best_start)
  d$date_end   <- as.Date(d$end_date)
  d$is_ltfu    <- d$itt_group == "Loss to follow-up"

  d$tx_end_d <- as.numeric(d$date_end - d$date_start)
  # defnB: recorded closure is the case-CLOSURE date, ~30 d after last contact.
  # Do NOT floor at 1. `pmax(tx_end_d - 30, 1)` piled 1,530 patients (7.3% of
  # LTFU) onto day 1 -- these are `Abandono Primario`, notified but never
  # meaningfully established on treatment, and they swamped the genuine early
  # days (494 at days 2-7, 473 at 8-14). They have no valid "day of
  # disengagement", so they are excluded from the day-indexed curve and analysed
  # separately by build_primary_abandonment().
  d$dis_d_raw <- ifelse(d$is_ltfu, d$tx_end_d - GRACE_D, NA_real_)
  # `Abandono Primario` comes from the OUTCOME CODE, not from the date. Patients
  # who DID start and stopped at once have a back-shifted day <= 0 only because
  # the closure sits 30 d after their single visit; they are floored to day 1.
  if (!is.null(outcome_lookup)) {
    d <- merge(d, outcome_lookup, by = "sinan_clean", all.x = TRUE)
    d$primary_aband <- d$is_ltfu & !is.na(d$primary_aband) & d$primary_aband
  } else {
    d$primary_aband <- rep(FALSE, nrow(d))
  }
  d$dis_d <- ifelse(d$is_ltfu & !d$primary_aband,
                    pmax(d$dis_d_raw, DAY_MIN), d$dis_d_raw)
  # when does this patient stop being "in care"?
  d$care_end_d <- ifelse(d$is_ltfu, pmax(d$dis_d_raw, 0), d$tx_end_d)
  d$fu_y <- d$time_d_tx
  d$ev   <- as.numeric(as.character(d$event_d))

  d$age_group <- factor(d$age_group, levels = c("15-24", "25-44", "45-64", "≥65"))
  # geo4 is excluded here: it is merged in below and its reference level is set
  # explicitly, so it must not be blindly as.factor()'d before it exists
  for (v in setdiff(COVARS, c("age_group", "geo4"))) d[[v]] <- as.factor(d[[v]])
  for (v in extra_factors) if (!is.null(d[[v]])) d[[v]] <- as.factor(d[[v]])

  # calendar period, for the subgroup analysis
  d$period <- factor(ifelse(as.integer(format(d$date_start, "%Y")) <= 2019,
                            "2013-2019", "2020-2023"))

  if (!is.null(cause_lookup)) {
    d <- merge(d, cause_lookup, by = "sinan_clean", all.x = TRUE)
    for (v in c("tb_hybrid", "nontb_hybrid", "tb_broad",
                "tb_simonly", "nontb_simonly")) {
      d[[v]] <- ifelse(is.na(d[[v]]), FALSE, d[[v]]) & (d$ev == 1)
    }
  }
  # geography (Reviewer 1 comment 11) -- always merged, so no calling script can
  # forget it and silently fit without a covariate that is now in COVARS
  d <- merge(d, build_geo_lookup(verbose = FALSE), by = "sinan_clean", all.x = TRUE)
  d$geo4 <- factor(ifelse(is.na(d$geo4), GEO_REF, as.character(d$geo4)))
  d$geo4 <- relevel(d$geo4, ref = GEO_REF)

  d$pid <- seq_len(nrow(d))
  d
}

# ---------------------------------------------------------------------------
# Stack the rolling trials
# ---------------------------------------------------------------------------
# trial   = a disengagement DAY d
# origin  = d + 30 days, i.e. the patient's own LTFU declaration date, so
#           follow-up NEVER begins before the definition is met
# exposed = disengaged on day d, alive at the origin
# compare = alive AND still in care at the origin (matched on time since
#           treatment initiation). Comparing against people who already
#           completed therapy is the immortal-time trap this design avoids.
build_rolling <- function(d, comparator = c("in_care", "alive_only"),
                          carry = character(0), seed = SEED,
                          defn = c("defnB", "defnA")) {
  comparator <- match.arg(comparator); defn <- match.arg(defn)
  # defnB: disengagement = closure - 30, grace 30  -> origin = closure
  # defnA: disengagement = closure,      grace 0   -> origin = closure
  # Identical origins, so the ESTIMATE is convention-invariant; only the
  # x-axis value (how much therapy the patient received) differs by 30 days.
  if (defn == "defnA") {
    d$dis_d <- ifelse(d$is_ltfu & !d$primary_aband,
                      pmax(d$tx_end_d, DAY_MIN + GRACE_D), d$dis_d)
    grace <- 0; dmin <- DAY_MIN + GRACE_D; dmax <- DAY_MAX + GRACE_D
  } else {
    grace <- GRACE_D; dmin <- DAY_MIN; dmax <- DAY_MAX
  }
  set.seed(seed)
  # primary abandonment is excluded here and handled by
  # build_primary_abandonment(); see the note in prepare_rolling()
  days <- sort(unique(floor(d$dis_d[d$is_ltfu & !d$primary_aband &
                                      d$dis_d >= dmin & d$dis_d <= dmax])))
  n <- nrow(d); out <- vector("list", length(days))
  for (i in seq_along(days)) {
    dd <- days[i]; origin_y <- (dd + grace) / 365.25
    alive <- d$fu_y > origin_y
    exp_l <- d$is_ltfu & !d$primary_aband & floor(d$dis_d) == dd & alive
    exp_l[is.na(exp_l)] <- FALSE
    if (!any(exp_l)) next
    pool_l <- if (comparator == "in_care") {
      alive & (d$care_end_d / 365.25) > origin_y & !exp_l
    } else alive & !exp_l
    pool_l[is.na(pool_l)] <- FALSE
    pool <- which(pool_l); exp_i <- which(exp_l)
    if (length(pool) < 2) next
    n_take <- min(length(pool), K_COMP * length(exp_i))
    cmp_i  <- if (length(pool) <= n_take) pool else sample(pool, n_take)
    idx <- c(exp_i, cmp_i)
    out[[i]] <- data.frame(
      pid = d$pid[idx], trial_day = dd,
      expose = c(rep(1L, length(exp_i)), rep(0L, length(cmp_i))),
      time_raw = d$fu_y[idx] - origin_y,
      event_d_num = d$ev[idx],
      d[idx, c(COVARS, carry), drop = FALSE],
      stringsAsFactors = FALSE)
  }
  bind_rows(out)
}

# ---------------------------------------------------------------------------
# Primary abandonment ("Abandono Primario"): notified but never meaningfully
# established on treatment, so there is no meaningful "day of disengagement".
# These trials are indexed by the DECLARATION day (their closure date, 1-30 d),
# which is also the origin, so arms remain aligned on time since notification.
# Interpretation differs from the timed curve and it must be reported separately.
# ---------------------------------------------------------------------------
build_primary_abandonment <- function(d, carry = character(0), seed = SEED) {
  set.seed(seed)
  days <- sort(unique(d$tx_end_d[d$primary_aband & d$tx_end_d >= 1 &
                                   d$tx_end_d <= 2 * GRACE_D]))
  n <- nrow(d); out <- vector("list", length(days))
  for (i in seq_along(days)) {
    dd <- days[i]; origin_y <- dd / 365.25
    alive <- d$fu_y > origin_y
    exp_l <- d$primary_aband & d$tx_end_d == dd & alive
    exp_l[is.na(exp_l)] <- FALSE
    if (!any(exp_l)) next
    pool_l <- alive & (d$care_end_d / 365.25) > origin_y & !exp_l
    pool_l[is.na(pool_l)] <- FALSE
    pool <- which(pool_l); exp_i <- which(exp_l)
    if (length(pool) < 2) next
    n_take <- min(length(pool), K_COMP * length(exp_i))
    cmp_i <- if (length(pool) <= n_take) pool else sample(pool, n_take)
    idx <- c(exp_i, cmp_i)
    out[[i]] <- data.frame(
      pid = d$pid[idx], trial_day = dd,
      expose = c(rep(1L, length(exp_i)), rep(0L, length(cmp_i))),
      time_raw = d$fu_y[idx] - origin_y,
      event_d_num = d$ev[idx],
      d[idx, c(COVARS, carry), drop = FALSE],
      stringsAsFactors = FALSE)
  }
  bind_rows(out)
}

# ---------------------------------------------------------------------------
# Outcome windows -- copied verbatim from 30h so estimates stay comparable.
# "late" counts events in (0.5, cap] WITHOUT left-truncating, so it does not
# condition on surviving to 6 months.
#
# cause: NULL for all-cause, else the name of a logical column. A death from a
# competing cause is censored at its own event time (event_out = 0), matching
# the cause-specific Cox in 30e/30i.
# ---------------------------------------------------------------------------
prep_outcome <- function(tr, model, cap, cause = NULL) {
  t <- tr$time_raw; e <- tr$event_d_num
  if (!is.null(cause)) e <- ifelse(tr[[cause]], 1, 0)
  if (model == "overall") {
    tr$time_out <- pmin(t, cap); tr$event_out <- ifelse(t > cap, 0, e)
  } else if (model == "early") {
    tr$time_out <- pmin(t, 0.5); tr$event_out <- ifelse(t <= 0.5 & e == 1, 1, 0)
  } else if (model == "late") {
    tr$time_out <- pmin(t, cap); tr$event_out <- ifelse(t > 0.5 & t <= cap & e == 1, 1, 0)
  } else stop("unknown model")
  tr[tr$time_out > 0, , drop = FALSE]
}

drop_constant <- function(tr, vars) {
  vars[vapply(vars, function(v)
    length(unique(tr[[v]][!is.na(tr[[v]])])) > 1, logical(1))]
}

fit_rolling <- function(tr, model, cap, cause = NULL, covars = COVARS,
                        timing = FALSE, min_events = 10) {
  tr <- prep_outcome(tr, model, cap, cause)
  if (sum(tr$event_out) < min_events) return(NULL)
  rhs <- paste(c("expose", drop_constant(tr, covars)), collapse = " + ")
  if (timing) rhs <- paste(rhs, "+ expose:splines::ns(trial_day, df = 3)")
  f <- as.formula(paste("Surv(time_out, event_out) ~", rhs, "+ strata(trial_day)"))
  tryCatch(coxph(f, data = tr, cluster = pid, ties = "efron"),
           error = function(e) NULL)
}

# Rubin pooling on the log-HR scale
pool_loghr <- function(est, se) {
  ok <- is.finite(est) & is.finite(se); est <- est[ok]; se <- se[ok]
  M <- length(est); if (!M) return(NULL)
  qbar <- mean(est); ubar <- mean(se^2)
  bvar <- if (M > 1) stats::var(est) else 0
  s <- sqrt(ubar + (1 + 1 / M) * bvar)
  list(hr = exp(qbar), lo = exp(qbar - 1.96 * s), hi = exp(qbar + 1.96 * s),
       p = 2 * stats::pnorm(-abs(qbar / s)), M = M)
}

pooled_expose <- function(fits) {
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) return(NULL)
  est <- vapply(fits, function(f) unname(coef(f)["expose"]), numeric(1))
  se  <- vapply(fits, function(f) sqrt(diag(vcov(f))["expose"]), numeric(1))
  pool_loghr(est, se)
}
