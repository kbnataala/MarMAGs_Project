################################################################################
# Script: 04_mag_quality_metrics.R
#
# Purpose:
#   Documentation template for visualizing MAG quality metrics across genome
#   sources (MarMAGs, GEM_catalog, OceanDNA).
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution. It is not presented as a fully reproducible
#   pipeline component.
#
# Expected input:
#   - results/tables/MarMAGs_Dataframe.csv
#
# Expected outputs:
#   - results/figures/04_quality_completeness.svg
#   - results/figures/04_quality_contamination.svg
#   - results/figures/04_quality_strain_heterogeneity.svg
#   - results/figures/04_quality_score.svg
#   - results/figures/04_quality_combined.png
#   - results/figures/04_quality_combined.svg
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(ggpubr)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

base_dir <- "."
input_file <- file.path(base_dir, "results", "tables", "MarMAGs_Dataframe.csv")
output_dir <- file.path(base_dir, "results", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

marmags_dataframe <- read.csv(
  file = input_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

required_cols <- c(
  "Genome_source",
  "CheckM_Completeness",
  "CheckM_Contamination",
  "CheckM_Strain_heterogeneity",
  "CheckM_Quality_score"
)

missing_cols <- setdiff(required_cols, colnames(marmags_dataframe))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input file: ",
       paste(missing_cols, collapse = ", "))
}

fig_data <- marmags_dataframe %>%
  dplyr::select(
    Genome_source,
    CheckM_Completeness,
    CheckM_Contamination,
    CheckM_Strain_heterogeneity,
    CheckM_Quality_score
  )

# Standardize factor order
fig_data$Genome_source <- factor(
  fig_data$Genome_source,
  levels = c("MarMAGs", "GEM_catalog", "OceanDNA")
)

# ------------------------------------------------------------------------------
# Helper function for violin + boxplot panels
# ------------------------------------------------------------------------------

make_quality_plot <- function(data, y_var, y_label, output_name, y_positions) {

  stat_pvalue_data <- data.frame(
    group1 = c("OceanDNA", "GEM_catalog", "OceanDNA"),
    group2 = c("GEM_catalog", "MarMAGs", "MarMAGs"),
    p.adj = c("***", "***", "***"),
    y.position = y_positions
  )

  p <- ggplot(
    data,
    aes(x = Genome_source, y = .data[[y_var]])
  ) +
    geom_violin(aes(color = Genome_source), fill = "white", trim = FALSE, linewidth = 1) +
    geom_boxplot(aes(color = Genome_source), fill = "white", width = 0.15, linewidth = 0.8) +
    scale_color_manual(values = c("MarMAGs" = "#1f77b4", "GEM_catalog" = "#d62728", "OceanDNA" = "#2ca02c")) +
    stat_pvalue_manual(
      data = stat_pvalue_data,
      label = "p.adj",
      size = 8,
      bracket.size = 0.8
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20),
      axis.line = element_line(linewidth = 0.8)
    ) +
    labs(
      x = "Genome source",
      y = y_label
    )

  ggsave(
    filename = file.path(output_dir, paste0(output_name, ".svg")),
    plot = p,
    width = 8,
    height = 8
  )

  return(p)
}

# ------------------------------------------------------------------------------
# Individual panels
# ------------------------------------------------------------------------------

completeness_plot <- make_quality_plot(
  data = fig_data,
  y_var = "CheckM_Completeness",
  y_label = "CheckM Completeness (%)",
  output_name = "04_quality_completeness",
  y_positions = c(108, 112, 116)
)

contamination_plot <- make_quality_plot(
  data = fig_data,
  y_var = "CheckM_Contamination",
  y_label = "CheckM Contamination (%)",
  output_name = "04_quality_contamination",
  y_positions = c(11, 12, 13)
)

strain_plot <- make_quality_plot(
  data = fig_data,
  y_var = "CheckM_Strain_heterogeneity",
  y_label = "CheckM Strain Heterogeneity",
  output_name = "04_quality_strain_heterogeneity",
  y_positions = c(120, 125, 130)
)

quality_score_plot <- make_quality_plot(
  data = fig_data,
  y_var = "CheckM_Quality_score",
  y_label = "CheckM Quality Score",
  output_name = "04_quality_score",
  y_positions = c(105, 107, 109)
)

# ------------------------------------------------------------------------------
# Combined figure
# ------------------------------------------------------------------------------

combined_plot <- plot_grid(
  quality_score_plot,
  completeness_plot,
  contamination_plot,
  strain_plot,
  label_size = 24,
  ncol = 4,
  labels = c("(a)", "(b)", "(c)", "(d)")
)

ggsave(
  filename = file.path(output_dir, "04_quality_combined.png"),
  plot = combined_plot,
  width = 24,
  height = 8,
  dpi = 600
)

ggsave(
  filename = file.path(output_dir, "04_quality_combined.svg"),
  plot = combined_plot,
  width = 24,
  height = 8
)
