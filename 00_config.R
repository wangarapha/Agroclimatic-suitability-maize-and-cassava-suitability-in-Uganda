# ==============================================================================
# 00_config.R
# Project configuration: crop selection, input/output paths and all analysis
# parameters used in the paper.
#
# Every pipeline script sources this file first. Edit paths here (not inside
# the individual scripts). Parameter values reproduce the published analysis;
# change them only if you intend to deviate from the paper.
# ==============================================================================

suppressPackageStartupMessages(library(here))

# ------------------------------------------------------------------------------
# Crop selection
# ------------------------------------------------------------------------------
# Set the crop with the CROP environment variable, e.g. from a shell:
#   CROP=cassava Rscript R/02_fit_ensemble_model.R
# or from R before sourcing a script:
#   Sys.setenv(CROP = "cassava")
CROP <- tolower(Sys.getenv("CROP", unset = "maize"))

CROPS <- list(
  maize = list(
    label          = "Maize",
    occurrence_csv = "maize_occurrence.csv",     # columns: longitude, latitude
    occurrence_shp = "maize_points_gbif.shp"     # optional, used in script 01
  ),
  cassava = list(
    label          = "Cassava",
    occurrence_csv = "cassava_occurrence.csv",
    occurrence_shp = NA_character_
  )
)

if (!CROP %in% names(CROPS)) {
  stop("Unknown CROP '", CROP, "'. Use one of: ",
       paste(names(CROPS), collapse = ", "))
}
CROP_LABEL <- CROPS[[CROP]]$label

# ------------------------------------------------------------------------------
# Input paths (see data/README.md for the expected layout)
# ------------------------------------------------------------------------------
DATA_DIR <- here("data")

PATHS <- list(
  occurrence_csv = file.path(DATA_DIR, "occurrence", CROPS[[CROP]]$occurrence_csv),
  occurrence_shp = if (is.na(CROPS[[CROP]]$occurrence_shp)) NA_character_ else
    file.path(DATA_DIR, "occurrence", CROPS[[CROP]]$occurrence_shp),
  env_current    = file.path(DATA_DIR, "climate", "current"),  # bioclim .tif files
  env_future     = file.path(DATA_DIR, "climate", "future"),   # one sub-folder per scenario
  uganda_shp     = file.path(DATA_DIR, "boundaries", "Uganda.shp"),
  subregions_shp = file.path(DATA_DIR, "boundaries", "Ug_Subregions.shp")
)

# ------------------------------------------------------------------------------
# Output paths (everything is written under outputs/<crop>/)
# ------------------------------------------------------------------------------
OUT_ROOT <- here("outputs", CROP)

OUT <- list(
  biomod      = file.path(OUT_ROOT, "biomod2"),      # biomod2 model files
  baseline    = file.path(OUT_ROOT, "baseline"),     # current-climate suitability
  projections = file.path(OUT_ROOT, "projections"),  # future suitability (.tif only)
  change      = file.path(OUT_ROOT, "change"),       # transition rasters
  tables      = file.path(OUT_ROOT, "tables"),
  figures     = file.path(OUT_ROOT, "figures")
)
invisible(lapply(OUT, dir.create, recursive = TRUE, showWarnings = FALSE))

# Key intermediate files shared between scripts
FILES <- list(
  baseline_tif = file.path(OUT$baseline, paste0(CROP, "_ensemble_EMwmean_current.tif")),
  model_out    = file.path(OUT$biomod, paste0(CROP, "_biomod_models.rds")),
  ensemble_out = file.path(OUT$biomod, paste0(CROP, "_biomod_ensemble.rds"))
)

# Small helpers to build output file names: out_table("x.csv") -> tables/<crop>_x.csv
out_table  <- function(name) file.path(OUT$tables,  paste0(CROP, "_", name))
out_figure <- function(name) file.path(OUT$figures, paste0(CROP, "_", name))

# ------------------------------------------------------------------------------
# Spatial settings
# ------------------------------------------------------------------------------
UTM_CRS   <- "EPSG:32636"  # WGS 84 / UTM zone 36N (Uganda)
WGS84_CRS <- "EPSG:4326"
AGG_FACT  <- 2             # aggregation factor applied to climate rasters before modelling

# ------------------------------------------------------------------------------
# Modelling parameters (biomod2), as used in the paper
# ------------------------------------------------------------------------------
PARAMS <- list(
  # Pseudo-absence sampling
  pa_nb_rep      = 2,
  pa_strategy    = "sre",
  pa_sre_quant   = 0.025,
  pa_min_dist_m  = 5000,   # pseudo-absences closer than this to a presence are dropped
  pa_seed        = NULL,   # NULL = not seeded, as in the original analysis (see docs/CODE_REVIEW_NOTES.md)

  # Single models
  models         = c("XGBOOST", "RF", "MAXNET"),
  cv_strategy    = "kfold",
  cv_k           = 2,      # 2-fold cross-validation
  cv_balance     = "presences",
  metric_eval    = c("TSS", "ROC"),
  var_import     = 2,      # permutations for variable importance
  scale_models   = TRUE,
  prevalence     = 0.5,
  seed_models    = 123,

  # Ensemble
  em_algo           = c("EMmean", "EMwmean"),
  em_metric_select  = "ROC",
  em_select_thresh  = 0.6, # keep single models with validation AUC >= 0.6
  em_select_dataset = "validation",
  em_var_import     = 3,
  em_ci_alpha       = 0.05,
  em_wmean_decay    = "proportional",
  em_seed           = 42
)

# Name of the weighted-mean ensemble layer produced by biomod2
EMWMEAN_LAYER <- paste0(CROP_LABEL, "_EMwmeanByROC_mergedData_mergedRun_mergedAlgo")

# ------------------------------------------------------------------------------
# Suitability classes
# ------------------------------------------------------------------------------
# Five suitability classes on the 0-1 ensemble probability scale.
# Used by 06_future_suitability_area.R (terra::classify with a from/to/becomes matrix)
SUIT_CLASS_MATRIX <- matrix(c(
  0,   0.1999999, 1,  # Not suitable        (< 0.2)
  0.2, 0.3999999, 2,  # Marginally suitable (0.2 - 0.399)
  0.4, 0.5999999, 3,  # Moderately suitable (0.4 - 0.599)
  0.6, 0.7999999, 4,  # Highly suitable     (0.6 - 0.799)
  0.8, 1.0,       5   # Very highly suitable (0.8 - 1.0)
), ncol = 3, byrow = TRUE)

# Used by 05_baseline_suitability_summary.R (breaks -> right-closed intervals)
SUIT_BREAKS <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)

SUIT_CLASS_NAMES <- c("Not Suitable", "Marginally Suitable", "Moderately Suitable",
                      "Highly Suitable", "Very Highly Suitable")

# Threshold separating "optimal" (>= 0.6) from "not optimal" suitability,
# used in the change-detection analysis
OPTIMAL_THRESHOLD <- 0.6

# ------------------------------------------------------------------------------
# Future scenarios
# ------------------------------------------------------------------------------
SSP_LEVELS    <- c("SSP126", "SSP245", "SSP370", "SSP585")
PERIOD_LEVELS <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
SSP_COLOURS   <- c(SSP126 = "green3", SSP245 = "blue", SSP370 = "orange", SSP585 = "red")
