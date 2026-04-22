# Shared project path resolution for R scripts in this repo.
#
# Source this file at the top of any script that needs to read from
# Data/ or ITT_Analysis/data/. Never hardcode a user-specific path.
#
# Resolution order:
#   1. TB_ABANDONMENT_ROOT environment variable (set for non-standard setups)
#   2. Known Google Drive mounts for current collaborators
#   3. Script-relative fallback (../../ from ITT_Analysis/scripts/)

find_project_root <- function() {
  env_root <- Sys.getenv("TB_ABANDONMENT_ROOT", unset = "")
  if (nzchar(env_root) && dir.exists(env_root)) {
    return(normalizePath(env_root))
  }
  home <- path.expand("~")
  candidates <- c(
    file.path(home, "Library/CloudStorage/GoogleDrive-jasonandr@gmail.com/My Drive/Abandonment Paper"),
    file.path(home, "Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper")
  )
  for (c in candidates) {
    if (dir.exists(c)) return(c)
  }
  # Last resort: walk up from the working directory
  normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
}

PROJECT_ROOT <- find_project_root()
DATA_DIR <- file.path(PROJECT_ROOT, "Data")
ITT_DATA_DIR <- file.path(PROJECT_ROOT, "ITT_Analysis", "data")
ITT_RESULTS_DIR <- file.path(PROJECT_ROOT, "ITT_Analysis", "results")
ITT_MI_DIR <- file.path(ITT_DATA_DIR, "mi")

COHORT_CSV <- file.path(ITT_DATA_DIR, "itt_cohort.csv")

message(sprintf("[paths] PROJECT_ROOT = %s", PROJECT_ROOT))
