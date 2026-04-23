# 53. Figure 4 — Within-LTFU multivariable forest (Mortality + Retreatment)
# ==============================================================================
# Combined table/forest rendering of Cox_Death_Adjusted_MI (left) and
# FG_Retr_Adjusted_MI (right) from multivariable_results_mi_cc.csv.
# Rows are grouped into three sections:
#   - Demographic   (age, race, education)
#   - Psychosocial  (homelessness, incarceration)
#   - Biomedical    (sex, HIV, diabetes, alcohol, drug use, hospital admission,
#                    clinical form, DOT, LTFU duration)
# Layout (5 patchwork columns sharing y-axis via identical facet_grid):
#   labels | AHR (95% CI) text | mortality forest | SHR (95% CI) text | retreatment forest
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
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

OUT_PNG <- file.path(ITT_RESULTS_DIR, "Figure_4_within_ltfu_forest.png")
OUT_PDF <- file.path(ITT_RESULTS_DIR, "Figure_4_within_ltfu_forest.pdf")

mv <- read.csv(file.path(ITT_RESULTS_DIR, "multivariable_results_mi_cc.csv"),
               stringsAsFactors = FALSE)

# Parse "1.23 (1.01-1.40)" → HR, CI_L, CI_H
parse_est <- function(s) {
  m <- regmatches(s, regexec("^([0-9.]+) \\(([0-9.]+)-([0-9.]+)\\)$", s))[[1]]
  if (length(m) == 4) {
    return(c(as.numeric(m[2]), as.numeric(m[3]), as.numeric(m[4])))
  }
  c(NA_real_, NA_real_, NA_real_)
}
parsed <- t(sapply(mv$estimate, parse_est))
mv$HR   <- parsed[, 1]
mv$CI_L <- parsed[, 2]
mv$CI_H <- parsed[, 3]

# Term → (Group, Label, display order) mapping.
# Order within each group reflects the display order top → bottom.
TERM_MAP <- list(
  list(term = "age_group25-44",                                                group = "Demographic",  label = "Age: 25-44 years"),
  list(term = "age_group45-64",                                                group = "Demographic",  label = "Age: 45-64 years"),
  list(term = "age_group≥65",                                             group = "Demographic",  label = "Age: ≥65 years"),
  list(term = "race_cleanBlack or Mixed",                                      group = "Demographic",  label = "Race: Black or Mixed"),
  list(term = "race_cleanOther",                                               group = "Demographic",  label = "Race: Other"),
  list(term = "edu_clean≤ 7 years",                                       group = "Demographic",  label = "Education: ≤7 years"),
  list(term = "edu_clean8 - 11 years",                                         group = "Demographic",  label = "Education: 8-11 years"),
  list(term = "edu_cleanNone",                                                 group = "Demographic",  label = "No formal education"),
  list(term = "homelessnessYes",                                               group = "Psychosocial", label = "Experiencing homelessness"),
  list(term = "incarceratedYes",                                               group = "Psychosocial", label = "Incarcerated"),
  list(term = "sexMale",                                                       group = "Biomedical",   label = "Male"),
  list(term = "hiv_aidsPositive",                                              group = "Biomedical",   label = "HIV-positive"),
  list(term = "diabetesYes",                                                   group = "Biomedical",   label = "Diabetes"),
  list(term = "alcoholYes",                                                    group = "Biomedical",   label = "Alcohol use"),
  list(term = "drug_useYes",                                                   group = "Biomedical",   label = "Drug use"),
  list(term = "hosp_admissionYes",                                             group = "Biomedical",   label = "Hospital admission"),
  list(term = "clinical_cleanExtrapulmonary",                                  group = "Biomedical",   label = "Extrapulmonary TB"),
  list(term = "clinical_cleanPulmonary and Extrapulmonary or disseminated",    group = "Biomedical",   label = "Mixed/disseminated TB"),
  list(term = "dot_statusYes",                                                 group = "Biomedical",   label = "Directly-observed therapy"),
  list(term = "tx_month_grp< 2 months",                                        group = "Biomedical",   label = "LTFU in <2 months"),
  list(term = "tx_month_grp2 to <4 months",                                    group = "Biomedical",   label = "LTFU in 2-4 months")
)
meta <- do.call(rbind, lapply(TERM_MAP, as.data.frame, stringsAsFactors = FALSE))
meta$order <- seq_len(nrow(meta))

# Join meta onto both models, keep display order
make_df <- function(model_name) {
  d <- mv |>
    dplyr::filter(model == model_name) |>
    dplyr::inner_join(meta, by = "term") |>
    dplyr::arrange(order)
  d$Group <- factor(d$group, levels = c("Demographic", "Psychosocial", "Biomedical"))
  # Factor levels top-to-bottom (reverse = bottom appears first in discrete axis)
  d$Label <- factor(d$label, levels = rev(unique(d$label)))
  d$hr_text <- sprintf("%.2f (%.2f–%.2f)", d$HR, d$CI_L, d$CI_H)
  d
}

cox <- make_df("Cox_Death_Adjusted_MI")
fg  <- make_df("FG_Retr_Adjusted_MI")

cat(sprintf("[fig4] %d Cox rows, %d FG rows\n", nrow(cox), nrow(fg)))

# Common facet + theme machinery -----------------------------------------------
common_strip <- theme(
  strip.placement       = "outside",
  strip.background      = element_blank(),
  panel.spacing.y       = unit(8, "pt")
)

# 1. Labels panel (leftmost) — strip on left shows group names
p_labels <- ggplot(cox, aes(y = Label)) +
  geom_text(aes(x = 0, label = Label),
            hjust = 0, size = 3.8, color = "grey20") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  facet_grid(rows = vars(Group), scales = "free_y", space = "free_y",
             switch = "y") +
  labs(title = "Characteristic") +
  theme_void(base_size = 11) +
  common_strip +
  theme(
    plot.title         = element_text(face = "bold", size = 12, hjust = 0,
                                      margin = margin(b = 10)),
    strip.text.y.left  = element_text(angle = 0, face = "bold", hjust = 0,
                                      size = 10.5, color = "#2c3e50",
                                      margin = margin(r = 10)),
    plot.margin        = margin(5, 0, 5, 5)
  )

# 2. Cox AHR text panel
p_cox_text <- ggplot(cox, aes(y = Label)) +
  geom_text(aes(x = 0, label = hr_text),
            hjust = 0.5, size = 3.4, color = "grey25") +
  scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  facet_grid(rows = vars(Group), scales = "free_y", space = "free_y") +
  labs(title = "AHR (95% CI)") +
  theme_void(base_size = 11) +
  common_strip +
  theme(
    plot.title   = element_text(face = "bold", size = 12, hjust = 0.5,
                                margin = margin(b = 10)),
    strip.text   = element_blank(),
    plot.margin  = margin(5, 0, 5, 0)
  )

# 3. Cox forest panel
p_cox_forest <- ggplot(cox, aes(x = HR, y = Label)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             linewidth = 0.5, color = "grey50") +
  geom_errorbarh(aes(xmin = CI_L, xmax = CI_H),
                 height = 0.28, linewidth = 0.7, color = "#2c3e50") +
  geom_point(size = 2.8, color = "#2c3e50") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4),
                minor_breaks = c(0.33, 0.75, 1.5, 3),
                limits = c(0.2, 6)) +
  facet_grid(rows = vars(Group), scales = "free_y", space = "free_y") +
  labs(title = "A. Adjusted mortality HR",
       x = "AHR", y = NULL) +
  theme_classic(base_size = 11) +
  common_strip +
  theme(
    plot.title          = element_text(face = "bold", size = 12, hjust = 0,
                                       margin = margin(b = 10)),
    axis.text.y         = element_blank(),
    axis.ticks.y        = element_blank(),
    axis.line.y         = element_blank(),
    strip.text          = element_blank(),
    panel.grid.major.x  = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.x  = element_line(color = "grey96", linewidth = 0.2),
    plot.margin         = margin(5, 6, 5, 0)
  )

# 4. FG SHR text panel
p_fg_text <- ggplot(fg, aes(y = Label)) +
  geom_text(aes(x = 0, label = hr_text),
            hjust = 0.5, size = 3.4, color = "grey25") +
  scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  facet_grid(rows = vars(Group), scales = "free_y", space = "free_y") +
  labs(title = "SHR (95% CI)") +
  theme_void(base_size = 11) +
  common_strip +
  theme(
    plot.title   = element_text(face = "bold", size = 12, hjust = 0.5,
                                margin = margin(b = 10)),
    strip.text   = element_blank(),
    plot.margin  = margin(5, 0, 5, 0)
  )

# 5. FG forest panel
p_fg_forest <- ggplot(fg, aes(x = HR, y = Label)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             linewidth = 0.5, color = "grey50") +
  geom_errorbarh(aes(xmin = CI_L, xmax = CI_H),
                 height = 0.28, linewidth = 0.7, color = "#b03a2e") +
  geom_point(size = 2.8, color = "#b03a2e") +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 2, 3),
                minor_breaks = c(0.75, 1.25, 1.75, 2.5),
                limits = c(0.5, 3.5)) +
  facet_grid(rows = vars(Group), scales = "free_y", space = "free_y") +
  labs(title = "B. Adjusted retreatment SHR",
       x = "SHR", y = NULL) +
  theme_classic(base_size = 11) +
  common_strip +
  theme(
    plot.title          = element_text(face = "bold", size = 12, hjust = 0,
                                       margin = margin(b = 10)),
    axis.text.y         = element_blank(),
    axis.ticks.y        = element_blank(),
    axis.line.y         = element_blank(),
    strip.text          = element_blank(),
    panel.grid.major.x  = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.x  = element_line(color = "grey96", linewidth = 0.2),
    plot.margin         = margin(5, 8, 5, 0)
  )

fig4 <- p_labels + p_cox_text + p_cox_forest + p_fg_text + p_fg_forest +
  plot_layout(widths = c(1.7, 0.75, 1.4, 0.75, 1.4)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  )

ggsave(OUT_PNG, fig4, width = 15, height = 9, dpi = 300, bg = "white")
ggsave(OUT_PDF, fig4, width = 15, height = 9, bg = "white")
cat(sprintf("[fig4] Wrote %s\n", OUT_PNG))
cat(sprintf("[fig4] Wrote %s\n", OUT_PDF))
