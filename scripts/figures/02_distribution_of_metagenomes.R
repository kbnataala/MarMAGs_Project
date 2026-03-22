################################################################################
# Script: 02_distribution_of_metagenomes.R
#
# Purpose:
#   Documentation template for generating summary plots describing the
#   environmental and geographic distribution of metagenomes / MAG sources.
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution. It is not presented as a fully reproducible
#   pipeline component.
#
# Expected inputs:
#   - data/processed/libs_data_edited.csv
#   - data/processed/MarMAGs_Dataframe.csv
#
# Expected outputs:
#   - results/figures/MarMAGs_Figure2.pdf (pie_charts_2c)
#   - results/figures/MarMAGs_Figure2.pdf (ocean_sea_barplot_2b)
#   - results/figures/MarMAGs_Figure2.pdf (IHO_sea_barplot_2d)
#   - results/figures/07_Marine_Region.svg
#   - results/figures/05_territory.svg
#   - results/figures/MarMAGs_Figure2.pdf (world_map_2a)
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(RColorBrewer)
  library(gridExtra)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

base_dir <- "."
input_dir <- file.path(base_dir, "data", "processed")
output_dir <- file.path(base_dir, "results", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

libs_file <- file.path(input_dir, "libs_data_edited.csv")
marmags_file <- file.path(input_dir, "MarMAGs_Dataframe.csv")

# ------------------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------------------

if (!file.exists(libs_file)) {
  stop("Missing input file: ", libs_file)
}

if (!file.exists(marmags_file)) {
  stop("Missing input file: ", marmags_file)
}

libs <- read.csv(libs_file, header = TRUE, stringsAsFactors = FALSE)
marmags_dataframe <- read.csv(marmags_file, header = TRUE, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

create_pie_chart <- function(data, category_column, title) {
  freq_table <- data %>%
    count(.data[[category_column]]) %>%
    mutate(
      percentage = n / sum(n) * 100,
      category_clean = ifelse(
        is.na(.data[[category_column]]) | .data[[category_column]] == "",
        "NA",
        .data[[category_column]]
      ),
      labels = paste0(category_clean, " (", round(percentage, 2), "%)")
    ) %>%
    arrange(desc(percentage))

  n_cols <- max(3, min(8, nrow(freq_table)))
  palette_colors <- brewer.pal(n = n_cols, "Set2")
  assigned_colors <- setNames(palette_colors[seq_len(min(nrow(freq_table), length(palette_colors)))],
                              freq_table$labels[seq_len(min(nrow(freq_table), length(palette_colors)))])

  ggplot(freq_table, aes(x = "", y = percentage, fill = labels)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = assigned_colors) +
    labs(title = title, fill = category_column) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
}

create_bar_plot <- function(data, column_name, title = "", top_n = 10) {
  freq_table <- data %>%
    mutate(
      value_clean = ifelse(
        is.na(.data[[column_name]]) | .data[[column_name]] == "",
        "NA",
        .data[[column_name]]
      )
    ) %>%
    count(value_clean, sort = TRUE) %>%
    mutate(
      percentage = (n / sum(n)) * 100,
      value_clean = fct_reorder(value_clean, n)
    ) %>%
    slice_max(n, n = top_n)

  ggplot(freq_table, aes(x = value_clean, y = n, fill = value_clean)) +
    geom_col(color = "black", width = 0.7) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")),
              hjust = -0.2, size = 4) +
    coord_flip() +
    scale_fill_brewer(palette = "Set3") +
    expand_limits(y = max(freq_table$n) * 1.1) +
    labs(title = title, x = NULL, y = "Frequency") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      legend.position = "none"
    )
}

save_svg_plot <- function(filename, width, height, plot_expr) {
  svg(filename, width = width, height = height)
  print(plot_expr)
  dev.off()
}

# ------------------------------------------------------------------------------
# Pie charts
# ------------------------------------------------------------------------------

salinity_pie <- create_pie_chart(libs, "Salinity.Categories", "Salinity Categories")
temperature_pie <- create_pie_chart(libs, "Temperature.Categories", "Temperature Categories")
depth_pie <- create_pie_chart(libs, "Depth.Categories", "Depth Categories")

svg(file.path(output_dir, "02_pie_charts.svg"), width = 16, height = 6)
grid.arrange(temperature_pie, salinity_pie, depth_pie, ncol = 3)
dev.off()

# ------------------------------------------------------------------------------
# Bar plots
# ------------------------------------------------------------------------------

ocean_sea_plot <- create_bar_plot(libs, "Ocean_and_Sea", top_n = 10)
save_svg_plot(file.path(output_dir, "03_ocean_sea.svg"), 14, 8, ocean_sea_plot)

iho_sea_plot <- create_bar_plot(libs, "IHO_Sea", top_n = 10)
save_svg_plot(file.path(output_dir, "06_IHO_sea.svg"), 12, 8, iho_sea_plot)

marine_region_plot <- create_bar_plot(libs, "Marine_Region", top_n = 12)
save_svg_plot(file.path(output_dir, "07_Marine_Region.svg"), 18, 8, marine_region_plot)

territory_plot <- create_bar_plot(libs, "Territory", top_n = 12)
save_svg_plot(file.path(output_dir, "05_territory.svg"), 18, 8, territory_plot)

# ------------------------------------------------------------------------------
# World map of genome sources
# ------------------------------------------------------------------------------

required_cols <- c("Genome_source", "sample_latitude", "sample_longitude")
missing_cols <- setdiff(required_cols, colnames(marmags_dataframe))
if (length(missing_cols) > 0) {
  stop("Missing required columns in MarMAGs_Dataframe.csv: ",
       paste(missing_cols, collapse = ", "))
}

fig_map <- marmags_dataframe[, required_cols]
fig_map <- fig_map[!is.na(fig_map$sample_latitude) & !is.na(fig_map$sample_longitude), ]

world_map <- map_data("world")

sample_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "gray85",
    color = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data = fig_map,
    aes(
      x = sample_longitude,
      y = sample_latitude,
      color = factor(Genome_source, levels = c("MarMAGs", "GEM_catalog", "OceanDNA"))
    ),
    size = 2.5,
    alpha = 0.6
  ) +
  scale_color_manual(
    values = c("MarMAGs" = "#E33539", "GEM_catalog" = "#1f77b4", "OceanDNA" = "#2ca02c"),
    drop = FALSE
  ) +
  coord_cartesian(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", color = "Genome source") +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "bottom"
  )

save_svg_plot(file.path(output_dir, "04_world_map.svg"), 12, 8, sample_map)
