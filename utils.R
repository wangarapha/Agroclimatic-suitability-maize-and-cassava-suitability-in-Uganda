# ==============================================================================
# utils.R
# Helper functions shared by the pipeline scripts.
# ==============================================================================

# ------------------------------------------------------------------------------
# Environmental rasters
# ------------------------------------------------------------------------------

#' Align a raster to a reference raster (CRS, extent and resolution).
#' Reprojects if the CRS differs, then crops and resamples (bilinear).
align_raster_to_ref <- function(r, ref) {
  if (!terra::same.crs(r, ref)) {
    r <- terra::project(r, terra::crs(ref))
    cat("  - Reprojected raster to match CRS of the reference raster\n")
  }
  if (is.null(terra::intersect(terra::ext(r), terra::ext(ref)))) {
    stop("Raster extents do not overlap with the reference raster.")
  }
  r <- terra::crop(r, terra::ext(ref))
  terra::resample(r, ref, method = "bilinear")
}

#' Read all .tif files in a folder and stack them on the grid of the first file.
load_env_stack <- function(folder) {
  env_files <- list.files(folder, pattern = "\\.tif$", full.names = TRUE)
  if (length(env_files) == 0) stop("No .tif files found in ", folder)

  ref <- terra::rast(env_files[1])
  cat("Aligning", length(env_files), "rasters in:", folder, "\n")
  aligned <- lapply(env_files, function(f) align_raster_to_ref(terra::rast(f), ref))
  terra::rast(aligned)
}

#' Prepare a modelling stack: align, aggregate and project to UTM.
#' The same function is used for current and future climate so that both are
#' processed identically.
prepare_env_utm <- function(folder, agg_fact = AGG_FACT, utm_crs = UTM_CRS) {
  env <- load_env_stack(folder)
  if (!is.null(agg_fact) && agg_fact > 1) {
    env <- terra::aggregate(env, fact = agg_fact)
  }
  terra::project(env, utm_crs)
}

# ------------------------------------------------------------------------------
# Clipping / reprojection
# ------------------------------------------------------------------------------

#' Crop and mask a raster to the Uganda boundary, then project it to WGS84.
clip_to_uganda_wgs84 <- function(r, ug, crs_out = WGS84_CRS) {
  ug_proj <- if (!terra::same.crs(r, ug)) terra::project(ug, r) else ug
  r <- terra::mask(terra::crop(r, ug_proj), ug_proj)
  terra::project(r, crs_out)
}

# ------------------------------------------------------------------------------
# Scenario names
# ------------------------------------------------------------------------------

#' Build a standard "SSPxxx_YYYY-YYYY" label from a folder or file name such as
#' "BCC-CSM2-MR_ssp126_2021-2040" or "SSP126_2021-2040".
#' Returns the input unchanged if no SSP/period can be found.
scenario_label <- function(x) {
  ssp    <- toupper(regmatches(x, regexpr("[Ss][Ss][Pp][0-9]{3}", x)))
  period <- regmatches(x, regexpr("[0-9]{4}-[0-9]{4}", x))
  if (length(ssp) == 0 || length(period) == 0) return(x)
  paste0(ssp, "_", period)
}

#' Split "SSPxxx_YYYY-YYYY" labels into ordered Scenario and Period factors.
add_scenario_period <- function(df, col = "Raster") {
  labels <- vapply(df[[col]], scenario_label, character(1))
  df$Scenario <- factor(sub("_.*$", "", labels), levels = SSP_LEVELS, ordered = TRUE)
  df$Period   <- factor(sub("^[^_]*_", "", labels), levels = PERIOD_LEVELS, ordered = TRUE)
  df
}

# ------------------------------------------------------------------------------
# Map helpers
# ------------------------------------------------------------------------------

lon_lab <- function(x) paste0(round(x, 1), "°E")

lat_lab <- function(y) {
  lab <- round(y, 1)
  paste0(abs(lab), ifelse(lab > 0, "°N", ifelse(lab < 0, "°S", "°")))
}

#' Axis breaks and padded limits (5 %) for faceted maps of a long data frame
#' with x/y columns.
map_extent <- function(df, pad = 0.05) {
  x_range <- range(df$x)
  y_range <- range(df$y)
  list(
    x_breaks = pretty(x_range, n = 3),
    y_breaks = pretty(y_range, n = 4),
    x_limits = x_range + c(-1, 1) * pad * diff(x_range),
    y_limits = y_range + c(-1, 1) * pad * diff(y_range)
  )
}

#' Uganda boundary as an sf object in WGS84.
uganda_sf_wgs84 <- function(ug, fallback_crs) {
  ug_sf <- sf::st_as_sf(ug)
  if (is.na(sf::st_crs(ug_sf))) sf::st_crs(ug_sf) <- sf::st_crs(fallback_crs)
  sf::st_transform(ug_sf, 4326)
}

# ------------------------------------------------------------------------------
# Saving
# ------------------------------------------------------------------------------

#' Save a ggplot as PNG (300 dpi) and PDF with the same base name.
#' `pdf_height` lets the PDF use a different height, as in some paper figures.
save_plot <- function(plot, path_no_ext, width, height, dpi = 300,
                      pdf_height = height, rds = FALSE) {
  ggplot2::ggsave(paste0(path_no_ext, ".png"), plot,
                  width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(path_no_ext, ".pdf"), plot,
                  width = width, height = pdf_height)
  if (rds) saveRDS(plot, paste0(path_no_ext, ".rds"))
  invisible(path_no_ext)
}

write_table <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE)
  cat("Saved:", path, "\n")
  invisible(path)
}
