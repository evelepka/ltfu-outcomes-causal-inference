import os
from pathlib import Path

# ==============================================================================
# MASTER PATH CONFIGURATION (PYTHON)
# ------------------------------------------------------------------------------
# Change only the BASE_PATH if the project folder moves.
# ==============================================================================

BASE_PATH = Path("/Users/evelynlepkadelima/Library/CloudStorage/GoogleDrive-evelynlepka@gmail.com/My Drive/Abandonment Outcomes/Abandonment Paper")

DATA_DIR = BASE_PATH / "Data"
CODE_DIR = BASE_PATH / "code"
FIG_DIR  = BASE_PATH / "figures"
DOC_DIR  = BASE_PATH / "Drafts"
TABLES_DIR = BASE_PATH / "Tables and results"

# Ensure directories exist
for d in [DATA_DIR, FIG_DIR, DOC_DIR, TABLES_DIR]:
    d.mkdir(parents=True, exist_ok=True)

print(f"Project Base Path configured to: {BASE_PATH}")
