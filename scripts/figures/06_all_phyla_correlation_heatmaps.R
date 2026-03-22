################################################################################
# Script: 06_all_phyla_correlation_heatmaps.R
#
# Purpose:
#   Documentation template for generating correlation heatmaps among phyla
#   across environmental gradients using prevalence profiles.
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution. It is not presented as a fully reproducible
#   pipeline component.
#
# Expected input:
#   - results/tables/merged_data.csv
#
# Expected outputs:
#   - results/figures/06_depth_phyla_correlation_heatmap.png
#   - results/figures/06_salinity_phyla_correlation_heatmap.png
#   - results/figures/06_temperature_phyla_correlation_heatmap.png
#   - results/figures/06_all_phyla_correlation_combined.png
#   - results/figures/06_all_phyla_correlation_combined.svg
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

base_dir <- "."
input_file <- file.path(base_dir, "results", "tables", "merged_data.csv")
output_dir <- file.path(base_dir, "results", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

merged_data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

required_cols <- c(
  "GTDB_Tk_Phylum",
  "Depth_Category",
  "Salinity_Category",
  "temp_Category"
)

missing_cols <- setdiff(required_cols, colnames(merged_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns in merged_data.csv: ",
       paste(missing_cols, collapse = ", "))
}

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

build_presence_matrix <- function(data, gradient_col, phylum_col = "GTDB_Tk_Phylum") {
  df <- data %>%
    filter(
      !is.na(.data[[gradient_col]]),
      !is.na(.data[[phylum_col]]),
      .data[[gradient_col]] != "Unknown",
      .data[[phylum_col]] != "Unknown"
    ) %>%
    count(.data[[gradient_col]], .data[[phylum_col]], name = "n") %>%
    pivot_wider(
      names_from = .data[[phylum_col]],
      values_from = n,
      values_fill = 0
    )

  as.data.frame(df)
}

top_variable_phyla <- function(mat, top_n = 50) {
  vars <- apply(mat, 2, var, na.rm = TRUE)
  vars <- sort(vars, decreasing = TRUE)
  names(vars)[seq_len(min(top_n, length(vars)))]
}

make_corr_heatmap <- function(data, gradient_col, output_prefix, top_n = 50) {

  wide_df <- build_presence_matrix(data, gradient_col)
  if (ncol(wide_df) < 3) {
    stop("Not enough phyla columns to compute correlation for ", gradient_col)
  }

  gradient_values <- wide_df[[gradient_col]]
  mat <- as.matrix(wide_df[, -1, drop = FALSE])

  # Keep top variable phyla only
  keep_phyla <- top_variable_phyla(mat, top_n = top_n)
  mat_sub <- mat[, keep_phyla, drop = FALSE]

  # Compute correlation across phyla profiles
  cor_mat <- cor(mat_sub, method = "spearman", use = "pairwise.complete.obs")

  # Diverging palette
  hm_cols <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

  hm <- pheatmap(
    cor_mat,
    color = hm_cols,
    breaks = seq(-1, 1, length.out = 101),
    border_color = NA,
    fontsize = 9,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45,
    main = paste("Phylum correlation across", gradient_col),
    silent = TRUE
  )

  png(
    filename = file.path(output_dir, paste0(output_prefix, ".png")),
    width = 2200,
    height = 1800,
    res = 220
  )
  grid.newpage()
  grid.draw(hm$gtable)
  dev.off()

  return(hm$gtable)
}

# ------------------------------------------------------------------------------
# Individual heatmaps
# ------------------------------------------------------------------------------

depth_gt <- make_corr_heatmap(
  data = merged_data,
  gradient_col = "Depth_Category",
  output_prefix = "06_depth_phyla_correlation_heatmap",
  top_n = 50
)

salinity_gt <- make_corr_heatmap(
  data = merged_data,
  gradient_col = "Salinity_Category",
  output_prefix = "06_salinity_phyla_correlation_heatmap",
  top_n = 50
)

temperature_gt <- make_corr_heatmap(
  data = merged_data,
  gradient_col = "temp_Category",
  output_prefix = "06_temperature_phyla_correlation_heatmap",
  top_n = 50
)

# ------------------------------------------------------------------------------
# Combined figure
# ------------------------------------------------------------------------------

combined <- arrangeGrob(
  grobs = list(depth_gt, salinity_gt, temperature_gt),
  ncol = 3
)

png(
  filename = file.path(output_dir, "06_all_phyla_correlation_combined.png"),
  width = 3600,
  height = 1600,
  res = 220
)
grid.newpage()
grid.draw(combined)
dev.off()

svg(
  filename = file.path(output_dir, "06_all_phyla_correlation_combined.svg"),
  width = 18,
  height = 8
)
grid.newpage()
grid.draw(combined)
dev.off()
