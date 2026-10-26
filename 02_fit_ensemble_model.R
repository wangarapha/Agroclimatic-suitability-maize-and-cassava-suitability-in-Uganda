# ==============================================================================
# 02_fit_ensemble_model.R
# Fit the biomod2 ensemble suitability model for one crop under current
# climate and project it to the current climate.
#
# Workflow
#   1. Load occurrences and current bioclimatic rasters (aggregated, UTM 36N)
#   2. Sample pseudo-absences with the surface range envelope (SRE) strategy
#   3. Drop pseudo-absences within 5 km of any presence
#   4. Fit XGBOOST, RF and MAXNET with k-fold cross-validation
#   5. Build EMmean / EMwmean ensembles (models with validation AUC >= 0.6)
#   6. Project the ensemble to current climate and save the EMwmean raster
#
# Inputs : occurrence CSV; data/climate/current/*.tif
# Outputs: biomod2/ (model files), biomod2/<crop>_biomod_*.rds,
#          baseline/<crop>_ensemble_EMwmean_current.tif,
#          figures/<crop>_presence_absence_points.png
#
# Author : Rapheal Wangalwa
# ==============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(biomod2)
})

cat("Crop:", CROP_LABEL, "\n")

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
occurrence <- read.csv(PATHS$occurrence_csv)
str(occurrence)
occurrence <- occurrence[, c("longitude", "latitude")]

# Current climate: aligned, aggregated (factor AGG_FACT) and projected to UTM
env_utm <- prepare_env_utm(PATHS$env_current)
print(env_utm)

occ_sf  <- st_as_sf(occurrence, coords = c("longitude", "latitude"), crs = 4326)
occ_utm <- st_transform(occ_sf, UTM_CRS)

# ------------------------------------------------------------------------------
# 2. Pseudo-absence sampling (SRE)
# ------------------------------------------------------------------------------
pa_args <- list(
  resp.var       = rep(1, nrow(occurrence)),
  expl.var       = env_utm,
  resp.xy        = st_coordinates(occ_utm),
  resp.name      = CROP_LABEL,
  dir.name       = OUT$biomod,
  filter.raster  = TRUE,
  PA.nb.rep      = PARAMS$pa_nb_rep,
  PA.nb.absences = nrow(occurrence),
  PA.strategy    = PARAMS$pa_strategy,
  PA.sre.quant   = PARAMS$pa_sre_quant
)
if (!is.null(PARAMS$pa_seed)) pa_args$seed.val <- PARAMS$pa_seed

pa_data <- do.call(BIOMOD_FormatingData, pa_args)

# Presences are stored first, followed by the pseudo-absences (NA response).
# NOTE: this assumes no presences were removed by filter.raster = TRUE;
# see docs/CODE_REVIEW_NOTES.md.
n_pres     <- nrow(occurrence)
n_total    <- length(pa_data@data.species)
pres_idx   <- seq_len(n_pres)
absence_idx <- (n_pres + 1):n_total
cat("Presence indices: 1 to", n_pres, "\n")
cat("Pseudo-absence indices:", n_pres + 1, "to", n_total, "\n")

pa_data@data.species[pres_idx]    <- 1
pa_data@data.species[absence_idx] <- 0

all_points <- data.frame(
  x        = pa_data@coord[, 1],
  y        = pa_data@coord[, 2],
  presence = ifelse(pa_data@data.species == 1, "Presence", "Absence")
)

# ------------------------------------------------------------------------------
# 3. Remove pseudo-absences within PARAMS$pa_min_dist_m of any presence
# ------------------------------------------------------------------------------
points_sf   <- st_as_sf(all_points, coords = c("x", "y"), crs = 32636)
presence_sf <- filter(points_sf, presence == "Presence")
absence_sf  <- filter(points_sf, presence == "Absence")

too_close <- st_is_within_distance(absence_sf, presence_sf, dist = PARAMS$pa_min_dist_m)
absence_filtered <- absence_sf[lengths(too_close) == 0, ]

points_filtered <- bind_rows(presence_sf, absence_filtered)
coords <- st_coordinates(points_filtered)

df_clean <- points_filtered %>%
  st_drop_geometry() %>%
  mutate(x = coords[, 1], y = coords[, 2])

cat("Presences:", sum(df_clean$presence == "Presence"),
    "| Pseudo-absences kept:", sum(df_clean$presence == "Absence"), "\n")

p_points <- ggplot(df_clean, aes(x = x, y = y, color = presence)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Presence" = "red", "Absence" = "blue")) +
  labs(title = "Presence and Absence Points",
       subtitle = paste(CROP_LABEL, "- points used for model fitting"),
       x = "Easting", y = "Northing", color = "Point Type") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_points)
save_plot(p_points, out_figure("presence_absence_points"), width = 8, height = 8)

# biomod2 expects a numeric response (1 = presence, 0 = absence)
df_clean$presence <- ifelse(df_clean$presence == "Presence", 1, 0)

# Re-format as presence/absence data (no further pseudo-absence sampling)
bm_data <- BIOMOD_FormatingData(
  resp.var  = df_clean$presence,
  expl.var  = env_utm,
  resp.xy   = df_clean[, c("x", "y")],
  resp.name = CROP_LABEL,
  dir.name  = OUT$biomod
)
bm_data

# ------------------------------------------------------------------------------
# 4. Single models
# ------------------------------------------------------------------------------
bm_models <- BIOMOD_Modeling(
  bm.format    = bm_data,
  modeling.id  = paste0(CROP_LABEL, "_Model"),
  models       = PARAMS$models,
  CV.strategy  = PARAMS$cv_strategy,
  CV.k         = PARAMS$cv_k,
  CV.balance   = PARAMS$cv_balance,
  metric.eval  = PARAMS$metric_eval,
  var.import   = PARAMS$var_import,
  scale.models = PARAMS$scale_models,
  prevalence   = PARAMS$prevalence,
  seed.val     = PARAMS$seed_models
)

# ------------------------------------------------------------------------------
# 5. Ensemble models
# ------------------------------------------------------------------------------
bm_ensemble <- BIOMOD_EnsembleModeling(
  bm.mod                = bm_models,
  models.chosen         = "all",
  em.by                 = "all",                    # ensemble over all runs
  em.algo               = PARAMS$em_algo,
  metric.select         = PARAMS$em_metric_select,  # weight/select by AUC (ROC)
  metric.select.thresh  = PARAMS$em_select_thresh,
  metric.select.dataset = PARAMS$em_select_dataset,
  metric.eval           = PARAMS$metric_eval,
  var.import            = PARAMS$em_var_import,
  EMci.alpha            = PARAMS$em_ci_alpha,
  EMwmean.decay         = PARAMS$em_wmean_decay,
  nb.cpu                = 1,
  seed.val              = PARAMS$em_seed,
  do.progress           = TRUE
)

# Keep the fitted objects for scripts 03 and 04
saveRDS(bm_models,   FILES$model_out)
saveRDS(bm_ensemble, FILES$ensemble_out)

# ------------------------------------------------------------------------------
# 6. Project the ensemble to current climate
# ------------------------------------------------------------------------------
current_proj <- BIOMOD_EnsembleForecasting(
  bm.em         = bm_ensemble,
  proj.name     = "Current",
  new.env       = env_utm,
  models.chosen = "all",
  metric.binary = "ROC",
  metric.filter = "ROC",
  do.stack      = TRUE,
  on_0_1000     = FALSE,
  output.format = ".img"
)

ensemble_preds <- get_predictions(current_proj)
print(names(ensemble_preds))

writeRaster(ensemble_preds[[EMWMEAN_LAYER]], FILES$baseline_tif, overwrite = TRUE)
cat("Saved current EMwmean suitability to:", FILES$baseline_tif, "\n")
