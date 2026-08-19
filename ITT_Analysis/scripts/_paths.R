# Shared project path resolution for R scripts in this repo.
#
# Source this file at the top of any script that needs to read from
# Data/ or ITT_Analysis/data/. Never hardcode a user-specific path.
#
# Resolution order:
#   1. TB_ABANDONMENT_ROOT environment variable (set for non-standard setups)
#   2. Known Google Drive mounts for current collaborators
#
# A candidate only counts if it actually LOOKS like the project root, i.e. it
# contains ITT_Analysis/data. A bare dir.exists() check was not enough: on a
# machine where none of the candidates were mounted, the old script-relative
# fallback resolved PROJECT_ROOT to the user's HOME directory and every
# downstream path silently became ~/ITT_Analysis/... — reads failed with
# confusing errors and writes would have landed in the home directory.
# There is deliberately NO fallback now. Not finding the root is a hard error.

.is_project_root <- function(p) {
  nzchar(p) && dir.exists(file.path(p, "ITT_Analysis", "data"))
}

find_project_root <- function() {
  env_root <- Sys.getenv("TB_ABANDONMENT_ROOT", unset = "")
  if (nzchar(env_root)) {
    if (!.is_project_root(env_root)) {
      stop(sprintf(
        "TB_ABANDONMENT_ROOT is set to '%s' but that is not a project root\n(no ITT_Analysis/data inside it).",
        env_root
      ), call. = FALSE)
    }
    return(normalizePath(env_root))
  }

  home <- path.expand("~")
  candidates <- c(
    file.path(home, "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/TB SP 2026/LTFU Paper"),
    file.path(home, "Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper"),
    file.path(home, "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper")
  )
  for (cand in candidates) {
    if (.is_project_root(cand)) return(cand)
  }

  stop(paste0(
    "Could not locate the project root (the Google Drive folder holding Data/ and ITT_Analysis/).\n",
    "Tried:\n  - ", paste(candidates, collapse = "\n  - "), "\n",
    "Mount Google Drive, or set TB_ABANDONMENT_ROOT to the folder containing ITT_Analysis/data."
  ), call. = FALSE)
}

PROJECT_ROOT <- find_project_root()
DATA_DIR <- file.path(PROJECT_ROOT, "Data")
ITT_DATA_DIR <- file.path(PROJECT_ROOT, "ITT_Analysis", "data")
ITT_RESULTS_DIR <- file.path(PROJECT_ROOT, "ITT_Analysis", "results")
ITT_MI_DIR <- file.path(ITT_DATA_DIR, "mi")

COHORT_CSV <- file.path(ITT_DATA_DIR, "itt_cohort.csv")

message(sprintf("[paths] PROJECT_ROOT = %s", PROJECT_ROOT))
