import subprocess
import os
import sys
from pathlib import Path

# Load configuration for paths
try:
    from config import BASE_PATH, CODE_DIR
except ImportError:
    print("Error: config.py not found. Please ensure you are running this from the code/ directory.")
    sys.exit(1)

def run_script(script_path, interpreter="python3"):
    """Helper to run a script and check execution."""
    print(f"\n>>> Running {interpreter} {script_path.name}...")
    try:
        if interpreter == "Rscript":
            # For R, we don't use -u but we want to see output
            result = subprocess.run([interpreter, str(script_path)], 
                                  capture_output=False, text=True, check=True)
        else:
            result = subprocess.run([interpreter, "-u", str(script_path)], 
                                  capture_output=False, text=True, check=True)
        print(f"DONE: {script_path.name}")
    except subprocess.CalledProcessError as e:
        print(f"FAILED: {script_path.name} with exit code {e.returncode}")
        # Optionally sys.exit(1) if you want the whole pipeline to stop on error
        pass

def main():
    print("="*60)
    print("MASTER ANALYSIS PIPELINE: Tuberculosis Abandonment Outcomes")
    print("="*60)
    print(f"Base Directory: {BASE_PATH}")

    # Change directory to CODE_DIR for imports to work
    os.chdir(CODE_DIR)

    # --- PHASE 1: Data Cleaning ---
    run_script(CODE_DIR / "00_clean_sinan.py")

    # --- PHASE 2: Eligibility & Cohort Selection ---
    run_script(CODE_DIR / "06_eligibility_cohort.py")

    # --- PHASE 3: Recoding & Table 1 Baseline ---
    run_script(CODE_DIR / "07_analysis_table1_survival.py")
    run_script(CODE_DIR / "generate_table1_marginal.py")

    # --- PHASE 4: Statistical Modeling (R) ---
    run_script(CODE_DIR / "impute_and_analyze.R", "Rscript")
    run_script(CODE_DIR / "analyze_mortality_macro_causes.R", "Rscript")

    # --- PHASE 5: Tables & Figures ---
    run_script(CODE_DIR / "make_tables_mortality.py")
    run_script(CODE_DIR / "make_tables_retreatment.py")
    run_script(CODE_DIR / "make_table_mortality_causes.py")
    run_script(CODE_DIR / "plot_strategic_panel.py")

    # --- PHASE 6: Manuscript Update ---
    # Choice of version (v3 is currently the most advanced)
    run_script(CODE_DIR / "integral_manuscript_update_v3.py")

    print("\n" + "="*60)
    print("FULL PIPELINE EXECUTION COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()
