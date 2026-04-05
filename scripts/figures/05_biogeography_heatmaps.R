################################################################################
# Script: 05_biogeography_heatmaps.R
#
# Purpose:
#   Documentation template for generating heatmaps of phylum prevalence across
#   environmental and geographic gradients.
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution. It is not presented as a fully reproducible
#   pipeline component.
#
# Expected input:
#   - path/to/input/table
#
# Expected outputs:
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S2)
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S3)
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S4)
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S4)
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S6)
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S14a)
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S14b)
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S14c)
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S15)
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S16)
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

base_dir <- "."
input_file <- file.path(base_dir, "results", "tables", "merged_data.csv")
output_dir <- file.path(base_dir, "results", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

merged_data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

required_cols <- c(
  "GTDB_Tk_Phylum",
  "Depth_Category",
  "Salinity_Category",
  "temp_Category",
  "ocean_sea_name",
  "IHO_Sea"
)

missing_cols <- setdiff(required_cols, colnames(merged_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns in merged_data.csv: ",
       paste(missing_cols, collapse = ", "))
}

# ------------------------------------------------------------------------------
# Optional detection of CBB-capable genomes
# ------------------------------------------------------------------------------

# Update these candidate column names to match your real table if needed
cbb_candidate_cols <- c(
  "CBB",
  "CBB_cycle",
  "Carbon_fixation_CBB",
  "Calvin_Benson_Bassham",
  "CBB_presence"
)

cbb_col <- cbb_candidate_cols[cbb_candidate_cols %in% colnames(merged_data)][1]

if (is.na(cbb_col) || length(cbb_col) == 0) {
  message("No explicit CBB column found. CBB-specific heatmaps will be skipped.")
  has_cbb <- FALSE
} else {
  has_cbb <- TRUE
}

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

compute_prevalence_table <- function(data, gradient_col, top_n = NULL, remove_unknown = TRUE) {

  required <- c("GTDB_Tk_Phylum", gradient_col)
  missing <- setdiff(required, colnames(data))
  if (length(missing) > 0) {
    stop("Missing columns for prevalence table: ", paste(missing, collapse = ", "))
  }

  df <- data %>%
    filter(!is.na(GTDB_Tk_Phylum), !is.na(.data[[gradient_col]]))

  if (remove_unknown) {
    df <- df %>%
      filter(
        GTDB_Tk_Phylum != "Unknown",
        .data[[gradient_col]] != "Unknown"
      )
  }

  # Count per phylum × category
  counts <- df %>%
    count(.data[[gradient_col]], GTDB_Tk_Phylum, name = "n")

  # Total genomes per gradient category
  totals <- df %>%
    count(.data[[gradient_col]], name = "total")

  prevalence <- counts %>%
    left_join(totals, by = gradient_col) %>%
    mutate(prevalence = n / total)

  # Keep top phyla if requested
  if (!is.null(top_n)) {
    top_phyla <- prevalence %>%
      group_by(GTDB_Tk_Phylum) %>%
      summarise(total_n = sum(n), .groups = "drop") %>%
      slice_max(total_n, n = top_n) %>%
      pull(GTDB_Tk_Phylum)

    prevalence <- prevalence %>%
      filter(GTDB_Tk_Phylum %in% top_phyla)
  }

  prevalence
}

plot_prevalence_heatmap <- function(prevalence_df, gradient_col, title = NULL) {

  prevalence_df <- prevalence_df %>%
    group_by(GTDB_Tk_Phylum) %>%
    mutate(total_n = sum(n)) %>%
    ungroup() %>%
    mutate(
      GTDB_Tk_Phylum = fct_reorder(GTDB_Tk_Phylum, total_n),
      gradient_value = .data[[gradient_col]]
    )

  ggplot(prevalence_df, aes(x = gradient_value, y = GTDB_Tk_Phylum, fill = prevalence)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "white",
      high = "#2166AC",
      labels = percent_format(accuracy = 1)
    ) +
    labs(
      x = gradient_col,
      y = "Phylum",
      fill = "Prevalence",
      title = title
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

save_heatmap <- function(plot_obj, filename, width = 10, height = 7) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = plot_obj,
    width = width,
    height = height
  )
}

# ------------------------------------------------------------------------------
# All MAGs heatmaps
# ------------------------------------------------------------------------------

temp_prev <- compute_prevalence_table(merged_data, "temp_Category", top_n = 20)
temp_plot <- plot_prevalence_heatmap(temp_prev, "temp_Category", "Phylum prevalence across temperature categories")
save_heatmap(temp_plot, "05_heatmap_temperature_all.svg", width = 10, height = 8)

sal_prev <- compute_prevalence_table(merged_data, "Salinity_Category", top_n = 20)
sal_plot <- plot_prevalence_heatmap(sal_prev, "Salinity_Category", "Phylum prevalence across salinity categories")
save_heatmap(sal_plot, "05_heatmap_salinity_all.svg", width = 10, height = 8)

depth_prev <- compute_prevalence_table(merged_data, "Depth_Category", top_n = 20)
depth_plot <- plot_prevalence_heatmap(depth_prev, "Depth_Category", "Phylum prevalence across depth categories")
save_heatmap(depth_plot, "05_heatmap_depth_all.svg", width = 12, height = 8)

ocean_prev <- compute_prevalence_table(merged_data, "ocean_sea_name", top_n = 15)
ocean_plot <- plot_prevalence_heatmap(ocean_prev, "ocean_sea_name", "Phylum prevalence across ocean/sea categories")
save_heatmap(ocean_plot, "05_heatmap_ocean_all.svg", width = 12, height = 8)

sea_prev <- compute_prevalence_table(merged_data, "IHO_Sea", top_n = 15)
sea_plot <- plot_prevalence_heatmap(sea_prev, "IHO_Sea", "Phylum prevalence across IHO sea categories")
save_heatmap(sea_plot, "05_heatmap_iho_sea_all.svg", width = 14, height = 8)

# ------------------------------------------------------------------------------
# Optional CBB-only heatmaps
# ------------------------------------------------------------------------------

if (has_cbb) {
  cbb_data <- merged_data %>%
    filter(!is.na(.data[[cbb_col]])) %>%
    filter(.data[[cbb_col]] == 1 | .data[[cbb_col]] == TRUE | .data[[cbb_col]] == "Yes")

  if (nrow(cbb_data) > 0) {
    temp_prev_cbb <- compute_prevalence_table(cbb_data, "temp_Category", top_n = 20)
    temp_plot_cbb <- plot_prevalence_heatmap(temp_prev_cbb, "temp_Category", "CBB-capable phyla across temperature categories")
    save_heatmap(temp_plot_cbb, "05_heatmap_temperature_cbb.svg", width = 10, height = 8)

    sal_prev_cbb <- compute_prevalence_table(cbb_data, "Salinity_Category", top_n = 20)
    sal_plot_cbb <- plot_prevalence_heatmap(sal_prev_cbb, "Salinity_Category", "CBB-capable phyla across salinity categories")
    save_heatmap(sal_plot_cbb, "05_heatmap_salinity_cbb.svg", width = 10, height = 8)

    depth_prev_cbb <- compute_prevalence_table(cbb_data, "Depth_Category", top_n = 20)
    depth_plot_cbb <- plot_prevalence_heatmap(depth_prev_cbb, "Depth_Category", "CBB-capable phyla across depth categories")
    save_heatmap(depth_plot_cbb, "05_heatmap_depth_cbb.svg", width = 12, height = 8)
  } else {
    message("CBB column found, but no CBB-positive rows detected.")
  }
}
