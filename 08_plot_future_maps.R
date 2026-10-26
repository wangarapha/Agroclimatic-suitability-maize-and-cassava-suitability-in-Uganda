# ==============================================================================
# 08_plot_future_maps.R
# Faceted map of continuous suitability (0-1) for every future scenario
# (4 SSPs x 4 periods), clipped to Uganda.
#
# Inputs : projections/SSPxxx_YYYY-YYYY.tif (script 04); Uganda boundary
# Outputs: figures/<crop>_future_suitability_maps.{png,pdf}
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
  library(grid)  # unit()
})

cat("Crop:", CROP_LABEL, "\n")

# ------------------------------------------------------------------------------
# Load, clip to Uganda and reproject to WGS84
# ------------------------------------------------------------------------------
r_files <- list.files(OUT$projections, pattern = "\\.tif$", full.names = TRUE)
r_stack <- rast(r_files)
names(r_stack) <- tools::file_path_sans_ext(basename(r_files))

uganda <- vect(PATHS$uganda_shp)

cat("Clipping and reprojecting rasters...\n")
r_stack_wgs <- rast(lapply(seq_len(nlyr(r_stack)),
                           function(i) clip_to_uganda_wgs84(r_stack[[i]], uganda)))
names(r_stack_wgs) <- names(r_stack)

uganda_sf <- uganda_sf_wgs84(uganda, r_stack)

df_long <- as.data.frame(r_stack_wgs, xy = TRUE, na.rm = TRUE) |>
  pivot_longer(cols = -c(x, y), names_to = "model", values_to = "value")

ext_info <- map_extent(df_long)

# ------------------------------------------------------------------------------
# Plot: 4 x 4 panels with a shared binned legend
# ------------------------------------------------------------------------------
p <- ggplot(df_long, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = uganda_sf, inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 1) +
  facet_wrap(~ model, ncol = 4) +
  scale_fill_viridis_b(
    name   = "Suitability",
    option = "H",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),  # 10 colour bands of width 0.1
    guide  = guide_colorsteps(
      direction       = "horizontal",
      barwidth        = unit(20, "lines"),
      barheight       = unit(1, "lines"),
      frame.colour    = "black",
      frame.linewidth = 1,
      ticks.colour    = NA,
      title.position  = "bottom",
      title.hjust     = 0.5
    )
  ) +
  scale_x_continuous(breaks = ext_info$x_breaks, labels = lon_lab) +
  scale_y_continuous(breaks = ext_info$y_breaks, labels = lat_lab) +
  coord_sf(xlim = ext_info$x_limits, ylim = ext_info$y_limits, expand = TRUE, crs = 4326) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", colour = "black"),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.text      = element_text(size = 10),
    axis.text.x      = element_text(color = "black", size = 8, face = "bold"),
    axis.text.y      = element_text(color = "black", size = 8, face = "bold"),
    axis.title       = element_text(face = "bold", size = 12)
  )
print(p)
save_plot(p, out_figure("future_suitability_maps"), width = 11, height = 12)
