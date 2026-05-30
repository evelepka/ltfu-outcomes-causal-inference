# 51. Figure 2 — Cumulative incidence of retreatment up to 24 months,
#                stratified by HIV, month of abandonment, homelessness,
#                and hospitalization at index.
# ==============================================================================
# Four-panel Aalen-Johansen cumulative incidence figure. Death is the competing
# risk, modeled via cmprsk::cuminc. Analysis is restricted to the LTFU subgroup
# (N=21,619) with time origin at abandonment and horizon = 2.0 years.
# ==============================================================================

suppressPackageStartupMessages({
  library(cmprsk)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
})

.here <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  frames <- sys.frames()
  for (f in rev(frames)) {
    of <- f$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  getwd()
}
source(file.path(.here(), "_paths.R"))

HORIZON <- 2.0  # 24 months
OUT_PNG <- file.path(ITT_RESULTS_DIR, "Figure_2_stratified_retreatment_24mo.png")
OUT_PDF <- file.path(ITT_RESULTS_DIR, "Figure_2_stratified_retreatment_24mo.pdf")

cohort <- read.csv(COHORT_CSV, stringsAsFactors = FALSE)
ltfu <- cohort |> dplyr::filter(itt_group == "Loss to follow-up")
cat(sprintf("[fig2] LTFU N = %d\n", nrow(ltfu)))

# --- Helpers ---------------------------------------------------------------
cif_plot_df <- function(cif_obj, group_levels) {
  do.call(rbind, lapply(group_levels, function(g) {
    obj <- cif_obj[[paste0(g, " 1")]]  # cause 1 = retreatment
    if (is.null(obj)) return(NULL)
    d <- data.frame(time = obj$time, est = obj$est,
                    group = as.character(g))
    d <- d[d$time <= HORIZON, , drop = FALSE]
    rbind(
      data.frame(time = 0, est = 0, group = as.character(g)),
      d,
      data.frame(time = HORIZON, est = tail(d$est, 1), group = as.character(g))
    )
  }))
}

make_panel <- function(df_plot, palette, legend_title, subtitle, pval_text,
                       group_order = NULL) {
  if (!is.null(group_order)) {
    df_plot$group <- factor(df_plot$group, levels = group_order)
  } else {
    df_plot$group <- factor(df_plot$group)
  }
  ggplot(df_plot, aes(x = time, y = est, color = group)) +
    geom_step(linewidth = 1.15) +
    annotate("text", x = 0.05, y = 0.58, label = pval_text,
             size = 4.5, hjust = 0, fontface = "italic") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, 0.60)) +
    scale_x_continuous(breaks = seq(0, HORIZON, by = 0.5),
                       labels = seq(0, 24, by = 6),
                       limits = c(0, HORIZON)) +
    scale_color_manual(name = legend_title, values = palette) +
    labs(subtitle = subtitle,
         x = "Months since loss to follow-up",
         y = "Cumulative incidence of retreatment") +
    theme_classic(base_size = 14) +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(face = "bold", size = 14))
}

gray_p <- function(cif_obj) {
  p <- cif_obj$Tests["1", "pv"]
  if (p < 0.001) "Gray's test p < 0.001" else sprintf("Gray's test p = %.3f", p)
}

# --- A) HIV ----------------------------------------------------------------
d_hiv <- ltfu |> dplyr::filter(!is.na(hiv_aids) & hiv_aids != "" & hiv_aids != "Unknown")
cif_hiv <- cuminc(d_hiv$time_rn, d_hiv$event_rn, group = d_hiv$hiv_aids, cencode = 0)
pA <- make_panel(
  cif_plot_df(cif_hiv, unique(d_hiv$hiv_aids)),
  palette = c("Negative" = "#2c3e50", "Positive" = "#e74c3c"),
  legend_title = "HIV status", subtitle = "A. HIV status",
  pval_text = gray_p(cif_hiv),
  group_order = c("Negative", "Positive")
)

# --- B) Month of abandonment ------------------------------------------------
d_tim <- ltfu |>
  dplyr::filter(!is.na(tx_month_grp) & tx_month_grp != "") |>
  dplyr::mutate(tx_month_grp = factor(tx_month_grp,
                                       levels = c("< 2 months", "2 to <4 months", "≥ 4 months")))
cif_tim <- cuminc(d_tim$time_rn, d_tim$event_rn,
                  group = d_tim$tx_month_grp, cencode = 0)
pB <- make_panel(
  cif_plot_df(cif_tim, levels(d_tim$tx_month_grp)),
  palette = c("< 2 months" = "#e41a1c",
              "2 to <4 months" = "#ff7f00",
              "≥ 4 months" = "#377eb8"),
  legend_title = "Month of LTFU", subtitle = "B. Month of LTFU",
  pval_text = gray_p(cif_tim),
  group_order = c("< 2 months", "2 to <4 months", "≥ 4 months")
)

# --- C) Homelessness --------------------------------------------------------
d_home <- ltfu |> dplyr::filter(!is.na(homelessness) & homelessness != "")
cif_home <- cuminc(d_home$time_rn, d_home$event_rn,
                   group = d_home$homelessness, cencode = 0)
pC <- make_panel(
  cif_plot_df(cif_home, unique(d_home$homelessness)),
  palette = c("No" = "#2c3e50", "Yes" = "#e74c3c"),
  legend_title = "Homelessness",
  subtitle = "C. Homelessness",
  pval_text = gray_p(cif_home),
  group_order = c("No", "Yes")
)

# --- D) Hospitalization at index -------------------------------------------
d_hosp <- ltfu |> dplyr::filter(!is.na(hosp_admission) & hosp_admission != "")
cif_hosp <- cuminc(d_hosp$time_rn, d_hosp$event_rn,
                   group = d_hosp$hosp_admission, cencode = 0)
pD <- make_panel(
  cif_plot_df(cif_hosp, unique(d_hosp$hosp_admission)),
  palette = c("No" = "#2c3e50", "Yes" = "#e74c3c"),
  legend_title = "Hospitalization at index",
  subtitle = "D. Hospitalization at index",
  pval_text = gray_p(cif_hosp),
  group_order = c("No", "Yes")
)

# --- Compose ---------------------------------------------------------------
fig2 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

ggsave(OUT_PNG, fig2, width = 12, height = 10, dpi = 300, bg = "white")
ggsave(OUT_PDF, fig2, width = 12, height = 10, bg = "white")
cat(sprintf("[fig2] Wrote %s\n", OUT_PNG))
cat(sprintf("[fig2] Wrote %s\n", OUT_PDF))

# --- Side outputs: per-panel CIF values at selected horizons ---------------
horizons <- c(0.5, 1, 1.5, 2)
get_cif_at <- function(cif_obj, g, t) {
  obj <- cif_obj[[paste0(g, " 1")]]
  if (is.null(obj)) return(NA)
  idx <- max(which(obj$time <= t), na.rm = TRUE)
  if (length(idx) == 0 || is.infinite(idx)) return(NA)
  obj$est[idx] * 100
}

tab <- function(cif_obj, groups, label) {
  d <- data.frame(
    panel = label,
    group = groups,
    month_6  = sapply(groups, function(g) get_cif_at(cif_obj, g, horizons[1])),
    month_12 = sapply(groups, function(g) get_cif_at(cif_obj, g, horizons[2])),
    month_18 = sapply(groups, function(g) get_cif_at(cif_obj, g, horizons[3])),
    month_24 = sapply(groups, function(g) get_cif_at(cif_obj, g, horizons[4]))
  )
  d
}

summary_tbl <- bind_rows(
  tab(cif_hiv,  sort(unique(d_hiv$hiv_aids)),        "HIV"),
  tab(cif_tim,  levels(d_tim$tx_month_grp),          "Month_abandonment"),
  tab(cif_home, sort(unique(d_home$homelessness)),   "Homelessness"),
  tab(cif_hosp, sort(unique(d_hosp$hosp_admission)), "Hospitalization")
)
write.csv(summary_tbl,
          file.path(ITT_RESULTS_DIR, "Figure_2_cif_values_24mo.csv"),
          row.names = FALSE)
cat("[fig2] Wrote per-panel CIF values at 6/12/18/24 months\n")
