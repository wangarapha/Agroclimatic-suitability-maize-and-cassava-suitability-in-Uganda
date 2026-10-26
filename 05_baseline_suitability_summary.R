# ==============================================================================
# 05_baseline_suitability_summary.R
# Area of each suitability class under current (baseline) climate, at country,
# region and sub-region level, with summary figures.
#
# Classes (ensemble probability): <=0.2 Not, 0.2-0.4 Marginally,
# 0.4-0.6 Moderately, 0.6-0.8 Highly, 0.8-1.0 Very highly suitable.
#
# Inputs : baseline/<crop>_ensemble_EMwmean_current.tif (from script 02)
#          Uganda sub-region shapefile (fields: Sub_Region, Ug_Region)
# Outputs: tables/<crop>_suitability_*.csv
#          figures/<crop>_baseline_{donut,treemap,classified_map,alluvial,lollipop}
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
  library(scales)
  library(treemapify)
  library(ggalluvial)
})

cat("Crop:", CROP_LABEL, "\n")

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------
suitability_raster <- rast(FILES$baseline_tif)
uganda_shape <- st_read(PATHS$subregions_shp, quiet = TRUE)
uganda_shape <- st_transform(uganda_shape, crs(suitability_raster))

# Pixel area in km2 (raster is in UTM metres)
res_m          <- res(suitability_raster)
pixel_area_km2 <- res_m[1] * res_m[2] / 1e6
cat(sprintf("Pixel resolution: %.0f x %.0f m (%.6f km2)\n", res_m[1], res_m[2], pixel_area_km2))

# ------------------------------------------------------------------------------
# 2. Classify into five suitability classes and clip to Uganda
# ------------------------------------------------------------------------------
raster_class_labels <- c("Very Low", "Low", "Moderate", "High", "Very High")

classified_raster <- classify(
  suitability_raster,
  rcl = matrix(c(SUIT_BREAKS[-length(SUIT_BREAKS)], SUIT_BREAKS[-1], 1:5), ncol = 3),
  include.lowest = TRUE
)
levels(classified_raster) <- data.frame(ID = 1:5, Class = raster_class_labels)

uganda_vect   <- vect(uganda_shape)
masked_raster <- mask(crop(classified_raster, uganda_vect), uganda_vect)

# ------------------------------------------------------------------------------
# 3. Pixel counts per class and sub-region
# ------------------------------------------------------------------------------
# One row per sub-region (same order as the polygons), one count column per class.
# The counts are bound to the polygon attributes by hand because
# extract(..., bind = TRUE) fails with a summary function in some terra versions.
class_counts <- terra::extract(masked_raster, uganda_vect,
                               fun = function(x) table(factor(x, levels = 1:5)))
extracted_df <- cbind(as.data.frame(uganda_vect),
                      as.data.frame(class_counts)[, -1, drop = FALSE])

# The 5 count columns (classes 1-5, in order) are named after the raster
# category column, e.g. "Class", "Class.1", ... "Class.4"
count_cols <- grep("^Class(\\.|$)", names(extracted_df), value = TRUE)
stopifnot(length(count_cols) == 5)
names(extracted_df)[match(count_cols, names(extracted_df))] <- paste0("Count_", 1:5)

for (k in 1:5) {
  extracted_df[[paste0("Area_km2_Class", k)]] <- extracted_df[[paste0("Count_", k)]] * pixel_area_km2
}

extracted_df <- extracted_df %>%
  mutate(
    Total_Count    = rowSums(across(starts_with("Count_"))),
    Total_Area_km2 = Total_Count * pixel_area_km2
  )
for (k in 1:5) {
  extracted_df[[paste0("Prop_Class", k)]] <- extracted_df[[paste0("Area_km2_Class", k)]] / extracted_df$Total_Area_km2
}

# ------------------------------------------------------------------------------
# 4. Country, region and sub-region summaries
# ------------------------------------------------------------------------------
country_totals <- extracted_df %>%
  summarise(across(starts_with("Area_km2_Class"), \(x) sum(x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Class", values_to = "Area_km2") %>%
  mutate(
    Country_Total_km2 = sum(Area_km2),
    Country_Percent   = 100 * Area_km2 / Country_Total_km2
  )
print(country_totals)

region_totals <- extracted_df %>%
  group_by(Ug_Region) %>%
  summarise(across(starts_with("Area_km2_Class"), \(x) sum(x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(starts_with("Area_km2_Class"), names_to = "Class", values_to = "Area_km2") %>%
  group_by(Ug_Region) %>%
  mutate(
    Region_Total_km2 = sum(Area_km2),
    Region_Percent   = 100 * Area_km2 / Region_Total_km2
  ) %>%
  ungroup()
print(region_totals)

subregion_long <- extracted_df %>%
  as_tibble() %>%
  select(Sub_Region, Ug_Region, starts_with("Area_km2_Class")) %>%
  pivot_longer(starts_with("Area_km2_Class"), names_to = "Class_ID", values_to = "Area_km2") %>%
  mutate(
    Class_ID   = as.integer(gsub("Area_km2_Class", "", Class_ID)),
    Class_Name = factor(SUIT_CLASS_NAMES[Class_ID], levels = SUIT_CLASS_NAMES)
  ) %>%
  select(-Class_ID) %>%
  group_by(Sub_Region) %>%
  mutate(
    Subregion_Total_km2 = sum(Area_km2, na.rm = TRUE),
    Subregion_Percent   = 100 * Area_km2 / Subregion_Total_km2
  ) %>%
  ungroup()
print(subregion_long)

# All three aggregation levels in one table
sub_sum <- subregion_long %>%
  group_by(Level = "Subregion", Group = Sub_Region, Class_Name) %>%
  summarise(Area_km2 = sum(Area_km2, na.rm = TRUE), .groups = "drop")
reg_sum <- subregion_long %>%
  group_by(Level = "Region", Group = Ug_Region, Class_Name) %>%
  summarise(Area_km2 = sum(Area_km2, na.rm = TRUE), .groups = "drop")
cty_sum <- subregion_long %>%
  group_by(Class_Name) %>%
  summarise(Area_km2 = sum(Area_km2, na.rm = TRUE), .groups = "drop") %>%
  mutate(Level = "Country", Group = "Uganda") %>%
  select(Level, Group, Class_Name, Area_km2)
all_sum <- bind_rows(sub_sum, reg_sum, cty_sum)

write_table(country_totals, out_table("suitability_country_totals.csv"))
write_table(region_totals,  out_table("suitability_region_totals.csv"))
write_table(subregion_long, out_table("suitability_subregion_long.csv"))
write_table(all_sum,        out_table("suitability_all_levels_summary.csv"))
write_table(subregion_long, out_table("suitability_treemap_data.csv"))

# ------------------------------------------------------------------------------
# 5. Figures
# ------------------------------------------------------------------------------
class_colours <- c("#001D3D", "lightgreen", "yellow", "#FF7F00", "firebrick4")

# 5a. Donut chart: share of national area per class
country_plot_df <- subregion_long %>%
  group_by(Class_Name) %>%
  summarise(Area_km2 = sum(Area_km2, na.rm = TRUE), .groups = "drop") %>%
  mutate(Percent = 100 * Area_km2 / sum(Area_km2),
         Label   = paste0(sprintf("%.1f", Percent), "%"))

donut_plot <- ggplot(country_plot_df, aes(x = 2, y = Area_km2, fill = Class_Name)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  xlim(1.2, 2.5) +  # creates the hole
  scale_fill_manual(values = class_colours) +
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5),
            size = 6, fontface = "bold") +
  labs(title = paste(CROP_LABEL, "Suitability Area in Uganda (%)"),
       fill = "Suitability class") +
  theme_void() +
  theme(legend.text     = element_text(size = 10),
        legend.position = "none",
        legend.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.title      = element_text(hjust = 0.5, face = "bold", size = 14))
print(donut_plot)
save_plot(donut_plot, out_figure("baseline_suitability_donut"), width = 6, height = 6)

# 5b. Treemap: sub-region composition by class
treemap_plot <- ggplot(subregion_long, aes(
  area     = Area_km2,
  fill     = Class_Name,
  subgroup = Sub_Region,
  label    = paste0(Sub_Region, "\n", sprintf("%.1f%%", Subregion_Percent))
)) +
  geom_treemap(color = "white") +
  geom_treemap_subgroup_border(color = "white") +
  geom_treemap_subgroup_text(place = "topleft", grow = TRUE, alpha = 0.4) +
  geom_treemap_text(place = "centre", grow = TRUE, reflow = TRUE, min.size = 2.5) +
  labs(title = "Subregion Suitability Area Composition", fill = "Suitability class") +
  scale_fill_manual(values = class_colours) +
  theme_void() +
  theme(legend.text     = element_text(size = 10),
        legend.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        legend.position = "bottom",
        plot.title      = element_text(hjust = 0.5, face = "bold", size = 14))
print(treemap_plot)
save_plot(treemap_plot, out_figure("baseline_suitability_treemap"), width = 8, height = 10)

# 5c. Classified suitability map with sub-region labels
classified_ug_wgs84 <- project(masked_raster, WGS84_CRS)
classified_df <- as.data.frame(classified_ug_wgs84, xy = TRUE, na.rm = TRUE)
colnames(classified_df) <- c("lon", "lat", "Class")

# Label positions computed in UTM (safe point-on-surface), then back to WGS84
uganda_labels <- uganda_shape %>%
  st_transform(32636) %>%
  mutate(geometry = st_point_on_surface(geometry)) %>%
  st_transform(4326)

suit_map <- ggplot() +
  geom_raster(data = classified_df,
              aes(x = lon, y = lat, fill = factor(Class, levels = raster_class_labels))) +
  geom_sf(data = uganda_shape, fill = NA, color = "black", linewidth = 0.8) +
  geom_sf_label(data = uganda_labels, aes(label = Sub_Region),
                size = 7, fontface = "bold", fill = "white", alpha = 0.8,
                color = "black", linewidth = 0.2,
                label.padding = unit(0.05, "lines")) +
  scale_fill_manual(values = class_colours) +
  coord_sf(crs = 4326, expand = TRUE) +
  labs(title = paste(CROP_LABEL, "Suitability Map (Uganda)"),
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(axis.title      = element_blank(),
        axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.border    = element_blank(),
        panel.grid      = element_blank(),
        legend.position = "none")
print(suit_map)
save_plot(suit_map, out_figure("baseline_classified_suitability_map"), width = 8, height = 10)
# A SpatRaster cannot be stored with saveRDS(), so keep it as a GeoTIFF
writeRaster(classified_ug_wgs84,
            file.path(OUT$baseline, paste0(CROP, "_classified_suitability_wgs84.tif")),
            overwrite = TRUE)

# 5d. Alluvial plot: area flow from sub-regions to classes
alluv_df <- select(subregion_long, Sub_Region, Class_Name, Area_km2)

alluvial_plot <- ggplot(alluv_df, aes(axis1 = Sub_Region, axis2 = Class_Name, y = Area_km2)) +
  geom_alluvium(aes(fill = Class_Name), alpha = 1) +
  geom_stratum(width = 0.3, alpha = 0.3, color = "black", linewidth = 0.5) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3, fontface = "bold") +
  scale_fill_viridis_d(option = "H", drop = FALSE) +
  scale_x_discrete(limits = c("Subregion", "Suitability class")) +
  labs(y = "Area (km²)", fill = "Suitability class",
       title = "Area Flow from Subregions to Suitability Classes") +
  theme_minimal() +
  theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        axis.text.x     = element_text(face = "bold", color = "black"),
        axis.title.y    = element_text(face = "bold", color = "black"))
print(alluvial_plot)
save_plot(alluvial_plot, out_figure("baseline_suitability_alluvial"), width = 10, height = 12)
saveRDS(alluv_df, out_table("suitability_alluvial_data.rds"))

# 5e. Lollipop plot: sub-region contributions within each class
lollipop_colours <- c(
  "Not Suitable"         = "red",
  "Marginally Suitable"  = "maroon4",
  "Moderately Suitable"  = "darkorange",
  "Highly Suitable"      = "blue",
  "Very Highly Suitable" = "blue4"
)

lollipop_plot <- ggplot(subregion_long,
                        aes(x = reorder(Sub_Region, Area_km2), y = Area_km2, color = Class_Name)) +
  geom_segment(aes(xend = Sub_Region, yend = 0), linewidth = 1.1, alpha = 0.8) +
  geom_point(size = 3.5) +
  geom_text(aes(label = comma(round(Area_km2))), hjust = -0.5, size = 3, show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ Class_Name, scales = "free_y", ncol = 2) +
  scale_color_manual(values = lollipop_colours, drop = FALSE) +
  labs(x = "Subregion", y = "Area (km²)", color = "Suitability class",
       title = "Subregional Area Contribution Within Suitability Classes") +
  theme_minimal() +
  theme(plot.title       = element_text(hjust = 0.5, face = "bold", color = "black"),
        strip.text       = element_text(face = "bold", color = "black", size = 10),
        legend.position  = "right",
        axis.text.y      = element_text(color = "black"),
        panel.grid.major = element_line(color = "darkgrey", linewidth = 0.05),
        panel.grid.minor = element_blank())
print(lollipop_plot)
save_plot(lollipop_plot, out_figure("baseline_suitability_lollipop"), width = 10, height = 8)
