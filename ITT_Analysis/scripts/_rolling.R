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
# Default 2, overridable. The manuscript reports five years throughout, so the
# cause-specific hazard ratios in the Results have to be derived at five as
# well -- until 2026-08-21 they were 2-year values sitting beside 5-year
# absolute risks in the same paragraph. Left at 2 by default so nothing else
# in the rolling family shifts silently.
HORIZON_Y <- as.numeric(Sys.getenv("HORIZON_Y", unset = "2"))
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
# MUNICIPALITY ADJUSTMENT (Reviewer 1, comment 1, R2 round)
# ---------------------------------------------------------------------------
# Patients in the same municipality are not independent: programme quality
# varies and municipality predicts disengagement strongly (ICC 8.4%, OR 3.17
# between the 95th and 5th percentile municipality, dAIC 1106; patient-level
# glmer, lme4 1.1.37). The CCW carries a municipality intercept in its
# censoring-weight model; this is the same adjustment for the landmark, so that
# a reviewer comparing the two designs does not find them adjusted differently
# -- the reason COVARS is shared between them in the first place.
#
# It is a RIDGE term, not a COVARS entry. As a plain factor the 640 indicators
# would be unpenalised fixed effects -- 218 municipalities have no LTFU event
# at all -- and every subgroup routine would treat it as an ordinary covariate.
# Ridge with theta = 1/sigma^2 is a Gaussian random intercept of variance
# sigma^2, and unlike coxme it leaves basehaz() and predict(newdata) working,
# which the g-formula standardisation in script 45 depends on.
#
# sigma^2 = 0.1151 is the ML estimate from a PATIENT-LEVEL model of the outcome
# this model predicts: glmer(died within 5 y ~ 14 covariates + (1|muni)) on the
# 120,366 patients with five years of follow-up or an earlier death. It is not
# taken from the stacked data, where the same patient recurs in up to 20 trials
# and a municipality effect is confounded with that repetition (coxme on the
# stack returns 0.63 and 1.76 in two windows of the same data).
MUNI_ADJUST <- identical(Sys.getenv("MUNI_ADJUST", unset = "0"), "1")
MUNI_SIGMA2 <- as.numeric(Sys.getenv("MUNI_SIGMA2", unset = "0.1151"))

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

  # OPT-IN, default OFF: TB_ANY_LINE=1 counts tuberculosis recorded ANYWHERE on
  # the death certificate, not only as the underlying cause. The analysis
  # extract carries CAUSABAS alone -- the Methods say so and call the estimate
  # conservative for it -- but the full certificate is in
  # `Banco de dados/LINKAGE SIM (1).xlsx`; 47b_build_tb_anyline_flag.py distils
  # it. Adds ~23% more tuberculosis deaths, sitting mostly under HIV, COPD and
  # lung cancer, i.e. exactly the deaths held out of both classes as ambiguous.
  # SECONDARY ANALYSIS. Superseded 2026-08-31: the owner promoted this to the
  # main text as a labelled secondary analysis (Methods 67, Results 138,
  # Discussion 151, appendix section 3.4), because Discussion 151 already
  # asserted that the 61% tuberculosis share is a lower bound BECAUSE
  # contributing mentions are not counted -- and this measures it: 26% more
  # deaths classified, share 61% -> 63%, risk ratio 2.37 -> 2.26 because the
  # extra deaths fall in both arms.
  #
  # It is still NOT the primary rule. Leave OFF by default: every headline
  # number, and every row of rolling_cause_cif.csv, is the underlying cause.
  # When running it, always set OUT_SUFFIX so the primary outputs survive.
  if (nzchar(Sys.getenv("TB_ANY_LINE"))) {
    fl <- file.path(DATA_DIR, "tb_any_line_flag.csv")
    if (!file.exists(fl)) stop("TB_ANY_LINE=1 but ", fl, " is missing; run 47b first")
    anyl <- read.csv(fl, colClasses = "character")
    anyl <- anyl[anyl$tb_anyline == "1", "sinan_clean"]
    add <- as.character(a$sinan_clean) %in% trimws(anyl)
    cat(sprintf("  [cause] TB_ANY_LINE=1: %d deaths reclassified as tuberculosis\n",
                sum(add & !tb_strict)))
    tb_strict <- tb_strict | add
  }
  # NEGATIVE CONTROL OUTCOME, added 2026-08-21. External causes (ICD-10 chapter
  # XX, V01-Y98: assault, self-harm, falls, poisoning, transport). Chosen because
  # it is strongly patterned by the confounders that drive disengagement --
  # homelessness 5.8x, alcohol 3.6x, drug use 2.4x on crude rates -- while
  # untreated tuberculosis has no plausible pathway to it. "All non-TB deaths"
  # could not serve this role, because loss to follow-up genuinely does cause
  # non-TB death.
  #
  # ADDITIVE ON PURPOSE: these deaths STAY inside nontb_hybrid. Removing them
  # would silently move every existing non-TB estimate and break the two-class
  # partition Figure 5 depends on. The negative control is a separate analysis,
  # not a third bar.
  #
  # Analyse it on the CAUSE-SPECIFIC HAZARD, never the CIF: a real TB excess
  # depletes the pool available to die of external causes, so a CIF would be
  # pushed down mechanically and would manufacture a reassuring null.
  external  <- known & grepl("^[VWXY]", cod)
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
             external      = external,
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

# City of the INDEX episode. Shared by build_geo_lookup() and
# build_muni_lookup() so the typology and the municipality term can never be
# derived from different city assignments for the same patient.
.index_city <- function() {
  if (!is.null(.rolling_cache$city)) return(.rolling_cache$city)
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
  first$tx_city <- ifelse(is.na(first$tx_city), first$tx_city_any,
                          first$tx_city)
  out <- first[, c("sinan_clean", "tx_city")]
  .rolling_cache$city <- out
  out
}

# Municipality as an adjustment term; see the MUNI_ADJUST note above.
build_muni_lookup <- function(verbose = TRUE) {
  out <- .index_city()
  if (verbose) cat(sprintf("  [muni] %d patients, %d municipalities\n",
                           nrow(out), length(unique(out$tx_city))))
  out
}

build_geo_lookup <- function(verbose = TRUE) {
  if (!is.null(.rolling_cache$geo)) return(.rolling_cache$geo)
  xw_path <- file.path(PROJECT_ROOT, "ITT_Analysis", "external",
                       "municipality_typology_sp.csv")
  if (!file.exists(xw_path))
    stop("[geo] missing crosswalk: run 47_build_ibge_typology.py first")
  xw <- read.csv(xw_path, stringsAsFactors = FALSE)[, c("municipality", "geo4")]
  names(xw)[1] <- "tx_city"   # crosswalk column is `municipality`

  out <- merge(.index_city(), xw, by = "tx_city", all.x = TRUE)
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
  # NO FLOOR (owner decision 2026-08-24). dis_d is the inferred disengagement day
  # and the timing analysis is indexed on months of therapy received, so a day of
  # <= 0 has no position on that axis. Leaving it negative lets the day-grid filter
  # below (dis_d >= dmin) drop these patients, and keeps origin_y anchored on the
  # real declaration date rather than on a floored one. Flooring to DAY_MIN both
  # asserted a timing the dates do not support and pushed the origin past the
  # actual declaration, reintroducing a few days of immortal time. Do not restore
  # the pmax(): the CCW keeps these patients on purpose, this design cannot.
  d$dis_d <- d$dis_d_raw
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
                "tb_simonly", "nontb_simonly", "external")) {
      d[[v]] <- ifelse(is.na(d[[v]]), FALSE, d[[v]]) & (d$ev == 1)
    }
  }
  # geography (Reviewer 1 comment 11) -- always merged, so no calling script can
  # forget it and silently fit without a covariate that is now in COVARS
  d <- merge(d, build_geo_lookup(verbose = FALSE), by = "sinan_clean", all.x = TRUE)
  d$geo4 <- factor(ifelse(is.na(d$geo4), GEO_REF, as.character(d$geo4)))
  d$geo4 <- relevel(d$geo4, ref = GEO_REF)

  # municipality: always attached, like geo4, so no calling script can forget
  # it. Used only when MUNI_ADJUST is set; see the note at the top.
  d <- merge(d, build_muni_lookup(verbose = FALSE), by = "sinan_clean",
             all.x = TRUE)
  d$muni <- factor(ifelse(is.na(d$tx_city), "__NA__", as.character(d$tx_city)))
  # Pooled municipality, derived here and only here, from the full cohort.
  # Attaching it alongside `muni` means every downstream data frame carries it,
  # including the ones scripts build themselves (script 49 does) rather than
  # through build_rolling().
  d <- with_muni_p(d)

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
    # same rule as defnB above: no floor, so the day grid excludes them
    d$dis_d <- ifelse(d$is_ltfu & !d$primary_aband, d$tx_end_d, d$dis_d)
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
      d[idx, c(COVARS, "muni", "muni_p",
               setdiff(carry, c("muni", "muni_p"))), drop = FALSE],
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
      d[idx, c(COVARS, "muni", "muni_p",
               setdiff(carry, c("muni", "muni_p"))), drop = FALSE],
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

# The municipality random-intercept term, or "" when MUNI_ADJUST is off or the
# data carry fewer than two municipalities (a subgroup stratum can collapse to
# one).
#
# frailty(), NOT ridge(). `ridge()` begins with cbind(...), which turns a factor
# into its INTEGER CODES: ridge(muni) fits one numeric column holding the
# alphabetical rank of the city, not 640 indicators. The symptom is an aHR
# identical to five decimals across a 40-fold range of theta, because a single
# column is barely penalised either way. frailty() expands the factor properly;
# its effective degrees of freedom then move with theta (88 at sigma^2 = 0.05,
# 377 at 2.0), which is the check that the term is real.
#
# theta is FIXED at 1/sigma^2 rather than estimated. Estimating it on the
# stacked data is not meaningful: the same patient recurs in up to 20 trials, so
# a municipality effect is confounded with that repetition.
# Municipalities below MUNI_MIN patients are pooled into one "__small__"
# level before the term is built.
#
# WHY. Sao Paulo has 640 municipalities but ~7.6 EFFECTIVE clusters: 429 of
# them hold fewer than 50 patients between them. Each still adds a parameter
# the penalised Newton step has to place, and on the thinner causes that is
# what stops the fit converging -- observed non-convergence at 600 inner
# iterations for non-TB, where the fits that do converge need 9 to 62. Pooling
# the small ones removes parameters that carry almost no information. The
# script already applies the same idea to covariates via
# collapse_sparse_levels(); this is that rule, applied to municipality.
MUNI_MIN <- as.integer(Sys.getenv("MUNI_MIN", unset = "2000"))

muni_pooled <- function(x) {
  x <- as.character(x); x[is.na(x)] <- "__NA__"
  n <- table(x)
  small <- names(n)[n < MUNI_MIN]
  x[x %in% small] <- "__small__"
  factor(x)
}

muni_term <- function(tr) {
  if (!MUNI_ADJUST) return("")
  if (is.null(tr$muni) || length(unique(tr$muni[!is.na(tr$muni)])) < 2) return("")
  # sparse = FALSE is required, not a preference. With the default sparse
  # representation coxph accepts no `init` at all -- every length is refused --
  # and the thin causes then have to start from zero, which is exactly what
  # script 46d documents as divergent (600 inner iterations without settling).
  # Dense, the term takes an ordinary `init` and converges in 4 steps. Dense is
  # only affordable because MUNI_MIN pools the small municipalities: 640 levels
  # dense costs 21 minutes a fit, 30 levels costs 31 seconds.
  sprintf(" + frailty(muni_p, distribution = 'gaussian', theta = %.10f, sparse = FALSE)",
          MUNI_SIGMA2)
}

# Attach the pooled factor. Callers run this on the data they are about to fit,
# so the pooling reflects that data set -- a bootstrap replicate pools by its
# own counts, as collapse_sparse_levels() does for covariates.
with_muni_p <- function(tr) {
  # Built whenever `muni` is present, NOT only when MUNI_ADJUST is on.
  # build_rolling() carries muni_p unconditionally, so gating its creation on
  # the flag made every landmark script fail with "undefined columns selected"
  # in the default configuration -- the one anyone cloning this repository
  # runs. The column is inert when the flag is off: nothing reads it.
  if (!is.null(tr$muni)) tr$muni_p <- muni_pooled(tr$muni)
  tr
}

# CERTAINTY-UNIT CLUSTER BOOTSTRAP (R2 round, Reviewer 1 comment 1)
# The landmark's intervals came from resampling patients, which treats patients
# in one municipality as independent. These resample by municipality with Sao
# Paulo city -- 33.5% of the cohort, ~7.6 effective clusters in all -- as a
# certainty unit: it enters every replicate with its patients resampled, and
# the other municipalities are resampled as whole blocks. Same design as
# CCW_analysis/ccw_random_intercept.py --boot-unit certainty, so both designs'
# intervals reflect the same clustering. Uses the REAL municipality (`muni`),
# never the pooled `muni_p`: pooling is for the model, not the sampling unit.
#
#   blk <- certainty_blocks(prepped[[k]])   # once per imputation
#   idx <- certainty_draw(blk)              # per replicate; d[idx, ] as usual
certainty_blocks <- function(d) {
  if (is.null(d$muni)) stop("certainty_blocks: data carry no `muni` column")
  m <- as.character(d$muni); m[is.na(m)] <- "__NA__"
  rows <- split(seq_len(nrow(d)), m)
  cap <- names(rows)[which.max(lengths(rows))]
  list(capital = rows[[cap]], others = unname(rows[names(rows) != cap]),
       capital_name = cap)
}

certainty_draw <- function(blk) {
  k <- length(blk$others)
  c(blk$capital[sample.int(length(blk$capital), length(blk$capital), replace = TRUE)],
    unlist(blk$others[sample.int(k, k, replace = TRUE)], use.names = FALSE))
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
  # No pooling here either: muni_p arrives from prepare_rolling(). See fit_cox.
  mterm <- muni_term(tr)
  rhs <- paste0(rhs, mterm)
  f <- as.formula(paste("Surv(time_out, event_out) ~", rhs, "+ strata(trial_day)"))
  # coxph cannot carry a penalised term and a robust cluster at once: both
  # frailty(muni) + cluster = pid and frailty(muni) + cluster(pid) fail with
  # "length of 'dimnames' [2] not equal to array extent". The cluster does not
  # affect the point estimate, only the variance, so the municipality fit keeps
  # the estimate and its standard error is NOT usable. Callers that quote an
  # interval must take it from boot_expose(), which resamples patients as whole
  # blocks and so handles the same repetition that cluster = pid approximates.
  if (nzchar(mterm)) {
    tryCatch(coxph(f, data = tr, ties = "efron"), error = function(e) NULL)
  } else {
    tryCatch(coxph(f, data = tr, cluster = pid, ties = "efron"),
             error = function(e) NULL)
  }
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

# One place that knows coxph cannot hold a penalised term and a robust cluster
# at the same time. Scripts that build their own formula call this instead of
# coxph() directly, so the rule lives here and not in five copies.
#
# When the municipality term is present the fit carries no cluster, so its
# standard error is NOT usable; the caller's interval must come from a
# bootstrap that resamples patients as blocks.
# The cluster is added as a FORMULA term, not as the `cluster =` argument.
# coxph resolves that argument in the environment of the FORMULA, not of the
# caller, so `cluster = data$pid` inside a function silently fails to find
# `data` and the fit dies. A `+ cluster(pid)` term is resolved from the model
# frame instead, and gives the identical coefficient and standard error.
#
# Errors are reported, not swallowed. An earlier version returned NULL on
# error, so a broken fit surfaced only as an unrelated failure further
# downstream; the cause took a separate debugging run to find.
fit_cox <- function(f, data, mterm, ...) {
  # Deliberately does NOT pool here. Pooling in two places was a real bug: the
  # stack was pooled on the whole cohort and then re-pooled inside the fit on
  # whatever subset it received, so a by-month analysis dropped cities that
  # cleared the threshold overall, and predict() on the un-re-pooled rows then
  # died with "factor ... has new levels BARUERI, BAURU, ...". The pooled
  # factor is a covariate definition: it is derived ONCE, in prepare_rolling(),
  # from the full cohort, and never recomputed.
  if (!nzchar(mterm)) f <- stats::update(f, . ~ . + cluster(pid))
  # The error is written straight to stderr, NOT raised as a warning: an
  # earlier version raised one and then muffled every warning three lines
  # later, so the message it existed to produce never appeared and a failing
  # fit looked like a silent NULL.
  withCallingHandlers(
    tryCatch(coxph(f, data = data, ties = "efron", ...),
             error = function(e) {
               cat(sprintf("      [fit_cox] ERRO: %s\n", conditionMessage(e)),
                   file = stderr())
               NULL
             }),
    warning = function(w) invokeRestart("muffleWarning"))
}

# ---------------------------------------------------------------------------
# PATIENT BLOCK BOOTSTRAP for the exposure log-HR
# ---------------------------------------------------------------------------
# WHY. A patient enters up to K_COMP + 1 trials, so stacked rows are not
# independent. Without a municipality term `cluster = pid` handles that
# analytically; with one, coxph refuses the combination (see fit_rolling), and
# an unclustered model standard error would be too narrow. Resampling PATIENTS
# as whole blocks addresses the same dependence without relying on the model's
# variance at all, and it is the convention the CCW side of this project
# already uses.
#
# Imputations are CYCLED, b %% M, as in ccw_v3.bootstrap(): one pass of B
# replicates then propagates sampling and imputation uncertainty together.
#
# Returns the same shape as pool_loghr() so callers can swap one for the other.
boot_expose <- function(stacks, model, cap, cause = NULL, covars = COVARS,
                        B = 300, seed = SEED, verbose = FALSE) {
  M <- length(stacks)
  blocks <- lapply(stacks, function(tr) split(seq_len(nrow(tr)), tr$pid))
  point <- pooled_expose(lapply(stacks, fit_rolling, model = model, cap = cap,
                                cause = cause, covars = covars))
  if (is.null(point)) return(NULL)
  set.seed(seed)
  est <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    k  <- (b - 1) %% M + 1
    bl <- blocks[[k]]
    idx <- unlist(bl[sample.int(length(bl), length(bl), replace = TRUE)],
                  use.names = FALSE)
    fit <- fit_rolling(stacks[[k]][idx, , drop = FALSE], model = model,
                       cap = cap, cause = cause, covars = covars)
    if (!is.null(fit)) est[b] <- unname(coef(fit)["expose"])
    if (verbose && b %% 25 == 0)
      cat(sprintf("      boot %d/%d\n", b, B), file = stderr())
  }
  ok <- is.finite(est)
  if (sum(ok) < 20) return(NULL)
  q <- stats::quantile(est[ok], c(0.025, 0.975), names = FALSE)
  list(hr = point$hr, lo = exp(q[1]), hi = exp(q[2]),
       p = NA_real_, M = M, B_ok = sum(ok), source = "patient block bootstrap")
}

# Which interval a caller should quote. With a municipality term the model's
# standard error carries no cluster and is too narrow (see fit_rolling), so the
# interval must come from the patient block bootstrap. Without one, the
# cluster-robust Rubin pooling is kept, so existing numbers are untouched.
# Validated on 2026-10-17: with MUNI_ADJUST off the bootstrap reproduces the
# robust interval to within Monte-Carlo error (width ratio 0.94 overall, 1.03
# late, B = 200), so the change of method is not what moves any estimate.
expose_ci <- function(stacks, model, cap, cause = NULL, covars = COVARS,
                      B = 300, seed = SEED, verbose = FALSE) {
  if (MUNI_ADJUST)
    boot_expose(stacks, model = model, cap = cap, cause = cause,
                covars = covars, B = B, seed = seed, verbose = verbose)
  else
    pooled_expose(lapply(stacks, fit_rolling, model = model, cap = cap,
                         cause = cause, covars = covars))
}

pooled_expose <- function(fits) {
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) return(NULL)
  est <- vapply(fits, function(f) unname(coef(f)["expose"]), numeric(1))
  se  <- vapply(fits, function(f) sqrt(diag(vcov(f))["expose"]), numeric(1))
  pool_loghr(est, se)
}
