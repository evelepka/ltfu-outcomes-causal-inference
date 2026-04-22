# ==============================================================================
# MASTER PATH CONFIGURATION (R)
# ------------------------------------------------------------------------------
# Change only the base_path if the project folder moves.
# ==============================================================================

base_path <- "/Users/evelynlepkadelima/Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper"

# Helpful directory variables
data_dir <- file.path(base_path, "Data")
fig_dir <- file.path(base_path, "figures")
doc_dir <- file.path(base_path, "Drafts")

# Load necessary libraries often used
if (!require("dplyr")) install.packages("dplyr")
if (!require("lubridate")) install.packages("lubridate")

cat("R Project Path configured to:", base_path, "\n")
