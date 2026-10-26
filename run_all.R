# ==============================================================================
# run_all.R
# Run the full pipeline for one or both crops.
#
# Usage (from the repository root):
#   Rscript run_all.R                 # both crops, all steps
#   Rscript run_all.R maize           # one crop, all steps
#   Rscript run_all.R cassava 5 8     # one crop, steps 5 to 8
# ==============================================================================

args  <- commandArgs(trailingOnly = TRUE)
crops <- if (length(args) >= 1) args[1] else c("maize", "cassava")
from  <- if (length(args) >= 2) as.integer(args[2]) else 1
to    <- if (length(args) >= 3) as.integer(args[3]) else 8

steps <- c(
  "01_extract_env_at_occurrences.R",
  "02_fit_ensemble_model.R",
  "03_evaluate_model.R",
  "04_project_future.R",
  "05_baseline_suitability_summary.R",
  "06_future_suitability_area.R",
  "07_change_detection.R",
  "08_plot_future_maps.R"
)

for (crop in crops) {
  Sys.setenv(CROP = crop)
  for (step in steps[from:to]) {
    message("\n>>> [", crop, "] ", step)
    # Each script runs in a fresh environment so steps stay independent
    source(file.path("R", step), local = new.env())
  }
}
