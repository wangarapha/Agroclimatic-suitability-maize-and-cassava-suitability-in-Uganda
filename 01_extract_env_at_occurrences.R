# ==============================================================================
# 01_extract_env_at_occurrences.R
# Summarise current bioclimatic conditions at crop occurrence points.
#
# Inputs : occurrence CSV (longitude, latitude); current climate rasters;
#          optional occurrence shapefile (e.g. GBIF points)
# Outputs: tables/<crop>_env_variable_stats_at_occurrences.csv
#          tables/<crop>_points_with_bioclim_values.csv   (if shapefile given)
#
# Author : Rapheal Wangalwa
# ==============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "utils.R"))

suppressPackageStartupMessages(library(terra))

cat("Crop:", CROP_LABEL, "\n")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
occurrence <- read.csv(PATHS$occurrence_csv)
occurrence <- occurrence[, c("longitude", "latitude")]

env_layers <- load_env_stack(PATHS$env_current)
print(env_layers)

# ------------------------------------------------------------------------------
# Summary statistics of each variable at occurrence points
# ------------------------------------------------------------------------------
# Assumes longitude/latitude are in the same CRS as the climate rasters
occ_pts <- vect(occurrence, geom = c("longitude", "latitude"), crs = crs(env_layers))

env_at_occ <- terra::extract(env_layers, occ_pts)
env_at_occ <- as.data.frame(env_at_occ)[, -1, drop = FALSE]  # drop ID column

env_stats <- data.frame(
  Variable = colnames(env_at_occ),
  Mean = sapply(env_at_occ, mean, na.rm = TRUE),
  SD   = sapply(env_at_occ, sd,   na.rm = TRUE),
  Min  = sapply(env_at_occ, min,  na.rm = TRUE),
  Max  = sapply(env_at_occ, max,  na.rm = TRUE)
)
print(env_stats)
write_table(env_stats, out_table("env_variable_stats_at_occurrences.csv"))

# ------------------------------------------------------------------------------
# Optional: attach bioclim values to the attributes of a point shapefile
# ------------------------------------------------------------------------------
if (!is.na(PATHS$occurrence_shp) && file.exists(PATHS$occurrence_shp)) {
  occ_shp <- project(vect(PATHS$occurrence_shp), crs(env_layers))
  bio_vals <- terra::extract(env_layers, occ_shp)  # first column is ID

  occ_env_df <- cbind(as.data.frame(occ_shp), bio_vals[, -1])
  write_table(occ_env_df, out_table("points_with_bioclim_values.csv"))
} else {
  cat("No occurrence shapefile configured for", CROP_LABEL, "- skipping.\n")
}
