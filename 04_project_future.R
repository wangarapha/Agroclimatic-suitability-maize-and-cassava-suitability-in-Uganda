# ==============================================================================
# 04_project_future.R
# Project the fitted ensemble model onto each future climate scenario.
#
# Each sub-folder of data/climate/future/ holds the bioclimatic rasters for one
# scenario/period (e.g. "BCC-CSM2-MR_ssp126_2021-2040" or "SSP126_2021-2040").
# Rasters are processed exactly like the current climate (align, aggregate,
# project to UTM 36N).
#
# Inputs : biomod2/<crop>_biomod_ensemble.rds; data/climate/future/<scenario>/
# Outputs: projections/SSPxxx_YYYY-YYYY.tif (EMwmean suitability, 0-1)
#
# Author : Rapheal Wangalwa
# ==============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(terra)
  library(biomod2)
})

cat("Crop:", CROP_LABEL, "\n")

bm_ensemble <- readRDS(FILES$ensemble_out)

scenario_folders <- list.dirs(PATHS$env_future, full.names = TRUE, recursive = FALSE)
if (length(scenario_folders) == 0) stop("No scenario folders found in ", PATHS$env_future)
print(basename(scenario_folders))

# Output names drop everything except SSP and period, so two folders for the
# same SSP/period (e.g. two GCMs) would overwrite each other
out_labels <- vapply(basename(scenario_folders), scenario_label, character(1))
if (anyDuplicated(out_labels)) {
  stop("Several scenario folders map to the same SSP/period: ",
       paste(unique(out_labels[duplicated(out_labels)]), collapse = ", "))
}

for (scen_path in scenario_folders) {
  scen_name <- basename(scen_path)
  cat("\n==============================\n")
  cat("Processing scenario:", scen_name, "\n")
  cat("==============================\n")

  # 1. Environmental stack for this scenario
  env_utm <- prepare_env_utm(scen_path)

  # 2. Ensemble forecasting
  future_em <- BIOMOD_EnsembleForecasting(
    bm.em         = bm_ensemble,
    proj.name     = scen_name,
    new.env       = env_utm,
    models.chosen = "all",
    metric.binary = "ROC",
    metric.filter = "ROC",
    do.stack      = TRUE,
    output.format = ".img",
    on_0_1000     = FALSE
  )

  # 3. Save the EMwmean layer
  preds <- get_predictions(future_em)
  layer_name <- grep("EMwmean", names(preds), value = TRUE)[1]
  if (is.na(layer_name)) {
    warning("No EMwmean layer found for ", scen_name, "; skipping save.")
    next
  }

  # Files are named "SSPxxx_YYYY-YYYY.tif" so later scripts can parse
  # scenario and period from the file name
  out_file <- file.path(OUT$projections, paste0(scenario_label(scen_name), ".tif"))
  writeRaster(preds[[layer_name]], filename = out_file, overwrite = TRUE)
  cat("  - Saved EMwmean raster to:", out_file, "\n")
}

cat("\nAll scenarios processed.\n")
