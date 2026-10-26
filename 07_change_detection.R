# ==============================================================================
# 07_change_detection.R
# Change in optimal suitability (ensemble probability >= 0.6) between the
# baseline and each future scenario.
#
# Transition classes
#   0 = Persistently unsuitable / persistent non-prime zone  (0 -> 0)
#   1 = Loss of optimal suitability / loss of prime zone      (1 -> 0)
#   2 = Persistent optimal suitability / persistent prime zone (1 -> 1)
#   3 = Gain of optimal suitability / gain of prime zone      (0 -> 1)
#
# Inputs : baseline/<crop>_ensemble_EMwmean_current.tif (script 02)
#          projections/SSPxxx_YYYY-YYYY.tif             (script 04)
#          Uganda boundary and sub-region shapefiles
# Outputs: change/<crop>_{baseline,future}_binary_0_6.tif,
#          change/<crop>_suitability_transitions_0_6.tif,
#          figures/<crop>_suitability_transitions,
#          tables/<crop>_transitions_area_{national,by_subregion}.csv
#
# Author : Rapheal Wangalwa
# ==============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

cat("Crop:", CROP_LABEL, "\n")

# ------------------------------------------------------------------------------
# 1. Load and clip rasters to Uganda (WGS84)
# ------------------------------------------------------------------------------
future_files <- list.files(OUT$projections, pattern = "\\.tif$", full.names = TRUE)
cat("Future maps found:", length(future_files), "\n")

uganda <- vect(PATHS$uganda_shp)

baseline_wgs <- clip_to_uganda_wgs84(rast(FILES$baseline_tif), uganda)
names(baseline_wgs) <- "baseline"

future_stack_wgs <- rast(lapply(future_files, function(f) clip_to_uganda_wgs84(rast(f), uganda)))
names(future_stack_wgs) <- tools::file_path_sans_ext(basename(future_files))

# ------------------------------------------------------------------------------
# 2. Binary maps: 0 = not optimal (< threshold), 1 = optimal (>= threshold)
# ------------------------------------------------------------------------------
binary_fun <- function(x) ifelse(is.na(x), NA, ifelse(x < OPTIMAL_THRESHOLD, 0, 1))

baseline_bin <- app(baseline_wgs, binary_fun)
names(baseline_bin) <- "baseline_bin"

future_bin_stack <- rast(lapply(seq_len(nlyr(future_stack_wgs)),
                                function(i) app(future_stack_wgs[[i]], binary_fun)))
names(future_bin_stack) <- names(future_stack_wgs)

# ------------------------------------------------------------------------------
# 3. Transition classes
# ------------------------------------------------------------------------------
# combo = baseline * 2 + future  ->  0:(0->0)  1:(0->1)  2:(1->0)  3:(1->1)
# recoded to the transition classes listed in the header
transition_rcl <- matrix(c(
  0, 0,  # 0 -> 0  persistently unsuitable
  1, 3,  # 0 -> 1  gain
  2, 1,  # 1 -> 0  loss
  3, 2   # 1 -> 1  persistent optimal
), ncol = 2, byrow = TRUE)

change_class_stack <- rast(lapply(seq_len(nlyr(future_bin_stack)), function(i) {
  combo <- app(c(baseline_bin, future_bin_stack[[i]]), fun = function(v) {
    ifelse(is.na(v[1]) | is.na(v[2]), NA, v[1] * 2 + v[2])
  })
  classify(combo, rcl = transition_rcl, include.lowest = TRUE)
}))
names(change_class_stack) <- names(future_bin_stack)

writeRaster(baseline_bin,
            file.path(OUT$change, paste0(CROP, "_baseline_binary_0_6.tif")), overwrite = TRUE)
writeRaster(future_bin_stack,
            file.path(OUT$change, paste0(CROP, "_future_binary_0_6.tif")), overwrite = TRUE)
writeRaster(change_class_stack,
            file.path(OUT$change, paste0(CROP, "_suitability_transitions_0_6.tif")), overwrite = TRUE)

# ------------------------------------------------------------------------------
# 4. Faceted transition map
# ------------------------------------------------------------------------------
map_class_labels <- c("Persistently unsuitable", "Loss of optimal suitability",
                      "Persistent optimal suitability", "Gain of optimal suitability")
map_class_colours <- setNames(c("#001D3D", "red", "lightyellow", "blue"), map_class_labels)

uganda_sf <- uganda_sf_wgs84(uganda, baseline_wgs)

change_long <- as.data.frame(change_class_stack, xy = TRUE, na.rm = TRUE) |>
  pivot_longer(cols = -c(x, y), names_to = "model", values_to = "change_code") |>
  mutate(change_cat = factor(change_code, levels = 0:3, labels = map_class_labels))

ext_info <- map_extent(change_long)

p_change <- ggplot(change_long, aes(x = x, y = y, fill = change_cat)) +
  geom_raster() +
  geom_sf(data = uganda_sf, inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 1) +
  facet_wrap(~ model, ncol = 4) +
  scale_fill_manual(name = "Suitability transitions", values = map_class_colours) +
  scale_x_continuous(breaks = ext_info$x_breaks, labels = lon_lab) +
  scale_y_continuous(breaks = ext_info$y_breaks, labels = lat_lab) +
  coord_sf(xlim = ext_info$x_limits, ylim = ext_info$y_limits, expand = TRUE, crs = 4326) +
  labs(x = "Longitude", y = "Latitude",
       title = paste("Transitions in", CROP_LABEL, "optimal suitability (baseline vs future)")) +
  theme_bw() +
  theme(
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", colour = "black"),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(size = 9),
    axis.text.x      = element_text(color = "black", size = 8, face = "bold"),
    axis.text.y      = element_text(color = "black", size = 8, face = "bold"),
    axis.title       = element_text(face = "bold", size = 12)
  )
print(p_change)
save_plot(p_change, out_figure("suitability_transitions"), width = 11, height = 12)

# ------------------------------------------------------------------------------
# 5. Area (km2) and percentage of each transition class - national
# ------------------------------------------------------------------------------
# Transitions are projected back to UTM (nearest neighbour) to get areas in km2
change_class_utm <- project(change_class_stack, UTM_CRS, method = "near")
res_xy        <- res(change_class_utm)
cell_area_km2 <- abs(res_xy[1] * res_xy[2]) / 1e6
cat("Cell area (km2):", cell_area_km2, "\n")

class_codes  <- 0:3
table_class_labels <- c("Persistent non-prime zone", "Loss of prime zone",
                        "Persistent prime zone", "Gain of prime zone")

area_national <- bind_rows(lapply(seq_len(nlyr(change_class_utm)), function(i) {
  r <- change_class_utm[[i]]
  bind_rows(lapply(class_codes, function(k) {
    n_cells <- as.numeric(global(r == k, "sum", na.rm = TRUE)[1, 1])
    data.frame(model = names(change_class_utm)[i], class_code = k,
               count = n_cells, area_km2 = n_cells * cell_area_km2)
  }))
})) |>
  mutate(class_label = factor(class_code, levels = class_codes, labels = table_class_labels)) |>
  mutate(area_km2 = round(area_km2, 2)) |>
  group_by(model) |>
  mutate(percentage = round(area_km2 / sum(area_km2, na.rm = TRUE) * 100, 2)) |>
  ungroup() |>
  select(model, class_code, class_label, count, area_km2, percentage)

print(area_national)
write_table(area_national, out_table("transitions_area_national.csv"))

# ------------------------------------------------------------------------------
# 6. Area and percentage of each transition class - by sub-region
# ------------------------------------------------------------------------------
subregions_utm <- project(vect(PATHS$subregions_shp), UTM_CRS)
subregions_utm$sub_id <- seq_len(nrow(subregions_utm))
attr_df <- as.data.frame(subregions_utm)[, c("sub_id", "Sub_Region")]

area_sub <- bind_rows(lapply(seq_len(nlyr(change_class_utm)), function(i) {
  r <- change_class_utm[[i]]
  bind_rows(lapply(class_codes, function(k) {
    # Number of cells of class k in each sub-region.
    # NOTE: touches = TRUE counts boundary cells in every polygon they touch.
    e <- terra::extract(r == k, subregions_utm, fun = sum, na.rm = TRUE, touches = TRUE)
    colnames(e)[1:2] <- c("sub_id", "count")
    e$count[is.na(e$count)] <- 0

    df <- merge(e, attr_df, by = "sub_id")
    df$area_km2   <- df$count * cell_area_km2
    df$class_code <- k
    df$model      <- names(change_class_utm)[i]
    df
  }))
})) |>
  mutate(class_label = factor(class_code, levels = class_codes, labels = table_class_labels)) |>
  group_by(model, Sub_Region) |>
  mutate(total_area_km2 = sum(area_km2, na.rm = TRUE),
         percentage = ifelse(total_area_km2 > 0, area_km2 / total_area_km2 * 100, NA_real_)) |>
  ungroup() |>
  mutate(area_km2 = round(area_km2, 2), percentage = round(percentage, 2)) |>
  select(model, Sub_Region, class_code, class_label, count, area_km2, percentage)

print(area_sub)
write_table(area_sub, out_table("transitions_area_by_subregion.csv"))
