# ==============================================================================
# install_packages.R
# Install the R packages used by the pipeline.
#
# For exact reproducibility we recommend renv:
#   install.packages("renv"); renv::init()   # then renv::snapshot()
# ==============================================================================

pkgs <- c(
  # Spatial data
  "terra", "sf",
  # Modelling (biomod2 and the algorithms it calls)
  "biomod2", "xgboost", "randomForest", "maxnet", "pROC",
  # Data handling
  "dplyr", "tidyr", "here",
  # Plotting
  "ggplot2", "patchwork", "scales", "treemapify", "ggalluvial"
)

missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing)

# biomod2 >= 4.2 is required (API used: CV.strategy, em.algo, get_predictions(), ...)
stopifnot(packageVersion("biomod2") >= "4.2")
