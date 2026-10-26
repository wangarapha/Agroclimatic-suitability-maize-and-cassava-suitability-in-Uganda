# ==============================================================================
# 06_future_suitability_area.R
# Reclassify the future suitability maps into five classes and compute the
# area of each class for Uganda and for each sub-region.
#
# "Suitable" in the summary tables and figures = Highly + Very highly suitable
# (ensemble probability >= 0.6).
#
# Inputs : projections/SSPxxx_YYYY-YYYY.tif (from script 04)
#          Uganda boundary and sub-region shapefiles
# Outputs: tables/<crop>_future_*.csv
#          figures/<crop>_future_suitable_area_line, _subregion_heatmap
#
# Author : Rapheal Wangalwa
# ==============================================================================

source(here::here("R", "00_config.R"))
source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

cat("Crop:", CROP_LABEL, "\n")

tif_files <- list.files(OUT$projections, pattern = "\\.tif$", full.names = TRUE)
if (length(tif_files) == 0) stop("No projections found in ", OUT$projections)

uganda_shp    <- vect(PATHS$uganda_shp)
uganda_subreg <- vect(PATHS$subregions_shp)

# ==============================================================================
# Part A. National level
# ==============================================================================
results <- list()

for (file_path in tif_files) {
  cat("Processing:", basename(file_path), "\n")
  r <- rast(file_path)

  uganda_proj <- project(uganda_shp, r)
  r_clipped   <- mask(crop(r, uganda_proj, snap = "out"), uganda_proj)
  if (all(is.na(values(r_clipped)))) next

  cell_area_km2 <- prod(res(r_clipped)) / 1e6
  r_class  <- classify(r_clipped, SUIT_CLASS_MATRIX, include.lowest = TRUE)
  freq_tab <- as.data.frame(freq(r_class))

  r_name <- tools::file_path_sans_ext(basename(file_path))
  result_df <- data.frame(
    Raster      = r_name,
    Class_ID    = 1:5,
    Class_Name  = SUIT_CLASS_NAMES,
    Pixel_Count = 0
  )
  # Fill in counts for the classes present (values outside 1-5 are ignored)
  for (i in seq_len(nrow(freq_tab))) {
    idx <- which(result_df$Class_ID == freq_tab$value[i])
    if (length(idx) > 0) result_df$Pixel_Count[idx] <- freq_tab$count[i]
  }

  result_df$Area_km2   <- result_df$Pixel_Count * cell_area_km2
  result_df$Percentage <- result_df$Area_km2 / sum(result_df$Area_km2) * 100
  results[[r_name]] <- result_df
}

national_df <- bind_rows(results, .id = "Raster_Original") %>%
  group_by(Raster) %>%
  mutate(
    Total_Pixels        = sum(Pixel_Count),
    Total_Area_km2      = sum(Area_km2),
    Percentage_of_Total = (Area_km2 / Total_Area_km2) * 100
  ) %>%
  ungroup()
print(national_df)
write_table(national_df, out_table("future_suitability_national.csv"))

# Highly + very highly suitable area per scenario/period
national_suitable <- national_df %>%
  filter(Class_Name %in% c("Highly Suitable", "Very Highly Suitable")) %>%
  group_by(Raster) %>%
  summarise(Suitable_Area_km2   = sum(Area_km2),
            Suitable_Percentage = sum(Percentage))
write_table(national_suitable, out_table("future_suitable_area_national.csv"))

plot_data <- national_suitable %>%
  add_scenario_period("Raster") %>%
  select(Raster, Scenario, Period, Suitable_Area_km2, Suitable_Percentage) %>%
  arrange(Scenario, Period)
print(plot_data)
write_table(plot_data, out_table("future_suitable_area_by_scenario_period.csv"))

# Line plot. NOTE: geom_smooth() draws a LOESS-smoothed curve through the four
# periods, not the raw values.
line_plot <- ggplot(plot_data, aes(x = Period, y = Suitable_Percentage,
                                   color = Scenario, group = Scenario)) +
  geom_smooth(se = FALSE, linewidth = 1.5, span = 0.6) +
  scale_color_manual(values = SSP_COLOURS) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(title = CROP_LABEL, x = "Time Period", y = "Suitable Area (%)",
       color = "Climate Scenario") +
  theme_bw() +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.text.x        = element_text(color = "black", face = "bold"),
    axis.text.y        = element_text(color = "black", face = "bold"),
    legend.position    = "right",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "black", linetype = "dotted", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "black", linetype = "dotted", linewidth = 0.5),
    axis.title.x       = element_text(face = "bold", size = 12),
    axis.title.y       = element_text(face = "bold", size = 12),
    panel.border       = element_rect(color = "black", linewidth = 1.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.8),
    plot.background    = element_rect(fill = "white", color = "black", linewidth = 1)
  ) +
  ylim(min(plot_data$Suitable_Percentage) - 2,
       max(plot_data$Suitable_Percentage) + 2)
print(line_plot)
save_plot(line_plot, out_figure("future_suitable_area_line"), width = 10, height = 6, rds = TRUE)

# ==============================================================================
# Part B. Sub-region level
# ==============================================================================
uganda_subreg$ID <- seq_len(nrow(uganda_subreg))
results <- list()

for (file_path in tif_files) {
  cat("Processing:", basename(file_path), "\n")
  r <- rast(file_path)

  subreg_proj <- project(uganda_subreg, r)
  r_clipped   <- mask(crop(r, subreg_proj, snap = "out"), subreg_proj)
  if (all(is.na(values(r_clipped)))) next

  cell_area_km2 <- prod(res(r_clipped)) / 1e6
  r_class <- classify(r_clipped, SUIT_CLASS_MATRIX, include.lowest = TRUE)

  # Class value of every cell, tagged with its sub-region
  vals <- terra::extract(r_class, subreg_proj)
  colnames(vals)[2] <- "Class_ID"
  vals$Sub_Region <- subreg_proj$Sub_Region[vals$ID]
  vals <- filter(vals, !is.na(Class_ID))

  region_class_df <- vals %>%
    group_by(Sub_Region, Class_ID) %>%
    summarise(Pixel_Count = n(), .groups = "drop")

  # Make sure every sub-region x class combination appears (0 if absent)
  r_name <- tools::file_path_sans_ext(basename(file_path))
  result_df <- expand.grid(Sub_Region = unique(subreg_proj$Sub_Region), Class_ID = 1:5) %>%
    left_join(region_class_df, by = c("Sub_Region", "Class_ID")) %>%
    mutate(
      Pixel_Count = if_else(is.na(Pixel_Count), 0L, Pixel_Count),
      Class_Name  = SUIT_CLASS_NAMES[Class_ID],
      Raster      = r_name,
      Area_km2    = Pixel_Count * cell_area_km2
    ) %>%
    group_by(Raster, Sub_Region) %>%
    mutate(
      Total_Area_Subregion    = sum(Area_km2),
      Percentage_of_Subregion = if_else(Total_Area_Subregion > 0,
                                        Area_km2 / Total_Area_Subregion * 100, NA_real_)
    ) %>%
    ungroup()

  results[[r_name]] <- result_df
}

subregion_df <- bind_rows(results) %>% mutate(Class_Order = Class_ID)
print(table(subregion_df$Class_Name))
write_table(subregion_df, out_table("future_suitability_5class_by_subregion.csv"))

# Wide table: "area km2 (percent %)" per class, one row per raster x sub-region
subregion_wide <- subregion_df %>%
  mutate(Area_Pct_Combined = paste0(round(Area_km2, 2), " km² (",
                                    round(Percentage_of_Subregion, 1), "%)")) %>%
  select(Raster, Sub_Region, Class_Name, Area_Pct_Combined) %>%
  pivot_wider(id_cols = c(Raster, Sub_Region), names_from = Class_Name,
              values_from = Area_Pct_Combined, values_fill = "0 km² (0%)")
write_table(subregion_wide, out_table("future_suitability_5class_combined_format.csv"))

# Highly + very highly suitable area per sub-region
subregion_suitable <- subregion_df %>%
  filter(Class_Name %in% c("Highly Suitable", "Very Highly Suitable")) %>%
  group_by(Raster, Sub_Region) %>%
  summarise(
    Suitable_Pixel_Count             = sum(Pixel_Count),
    Suitable_Area_km2                = sum(Area_km2),
    Suitable_Percentage_of_Subregion = sum(Percentage_of_Subregion),
    .groups = "drop"
  ) %>%
  arrange(Raster, Sub_Region)
print(subregion_suitable)
write_table(subregion_suitable, out_table("future_suitable_area_by_subregion.csv"))

plot_data_sub <- subregion_suitable %>%
  add_scenario_period("Raster") %>%
  select(Raster, Scenario, Period, Sub_Region, Suitable_Area_km2,
         Suitable_Percentage = Suitable_Percentage_of_Subregion) %>%
  mutate(Sub_Region = factor(Sub_Region, levels = sort(unique(Sub_Region)))) %>%
  arrange(Sub_Region, Scenario, Period)
write_table(plot_data_sub, out_table("future_suitable_area_by_subregion_scenario_period.csv"))

# Heatmap of suitable area (%) by scenario and period, one panel per sub-region
# NOTE: values above 95 % fall outside the colour scale limits
plot_heatmap <- ggplot(plot_data_sub, aes(x = Period, y = Scenario, fill = Suitable_Percentage)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(Suitable_Percentage, 1)), size = 3) +
  scale_fill_viridis_b(
    name   = "Suitable Area (%)",
    option = "H",
    guide  = guide_colorsteps(
      barheight       = unit(16, "lines"),
      barwidth        = unit(0.8, "lines"),
      frame.colour    = "black",
      frame.linewidth = 1,
      ticks.colour    = NA
    ),
    limits = c(0, 95),
    breaks = seq(0, 95, by = 5)
  ) +
  facet_wrap(~ Sub_Region, ncol = 3) +
  labs(x = "Time Period", y = "Scenario") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 30, hjust = 0.5, color = "black", size = 8, face = "bold"),
    axis.text.y      = element_text(color = "black", size = 8, face = "bold"),
    axis.ticks.x     = element_line(color = "black"),
    axis.ticks.y     = element_line(color = "black"),
    axis.title.x     = element_text(face = "bold", size = 12),
    axis.title.y     = element_text(face = "bold", size = 12),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    legend.text      = element_text(size = 10),
    legend.title     = element_text(face = "bold", hjust = 0.5)
  )
print(plot_heatmap)
save_plot(plot_heatmap, out_figure("future_suitable_area_subregion_heatmap"),
          width = 11, height = 12, rds = TRUE)
