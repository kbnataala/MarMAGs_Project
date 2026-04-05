################################################################################
# Script: 07_prevalence_plots.R
#
# Purpose:
#   Documentation template for generating prevalence barplots with significance
#   - results/figures/MarMAGs_Figure5.pdf (heamap_5a)
#   - results/figures/MarMAGs_Figure5.pdf (heamap_5b)
#   - results/figures/MarMAGs_Figure5.pdf (heamap_5c)
#   - results/figures/MarMAGs_Figure5.pdf (heamap_5d)
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S11)
#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S12)
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

base_dir <- "."
input_dir <- file.path(base_dir, "results", "tables")
output_dir <- file.path(base_dir, "results", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Core plotting function (cleaned version)
# ------------------------------------------------------------------------------

plot_prevalence <- function(prev_file, comp_file, group_col, title = NULL) {

  df_prev <- read.csv(prev_file, stringsAsFactors = FALSE)
  df_comp <- read.csv(comp_file, stringsAsFactors = FALSE)

  df_prev <- df_prev %>%
    filter(Prevalence > 0.01, Number_of_taxa > 10)

  df_prev[[group_col]] <- factor(df_prev[[group_col]],
                                levels = df_prev[[group_col]][order(df_prev$Prevalence)])

  p <- ggplot(df_prev, aes_string(x = group_col, y = "Prevalence")) +
    geom_col(fill = "#4C72B0") +
    coord_flip() +
    scale_y_continuous(labels = percent_format()) +
    theme_classic(base_size = 14) +
    labs(x = NULL, y = "Prevalence", title = title)

  return(p)
}

# ------------------------------------------------------------------------------
# Generate plots
# ------------------------------------------------------------------------------

p1 <- plot_prevalence(
  file.path(input_dir, "pathway_phylum.csv"),
  file.path(input_dir, "phylum_pairwise_prevalence.csv"),
  "GTDB_Tk_Phylum",
  "Phylum"
)

p2 <- plot_prevalence(
  file.path(input_dir, "pathway_depth.csv"),
  file.path(input_dir, "depth_pairwise_prevalence.csv"),
  "Depth_Category",
  "Depth"
)

p3 <- plot_prevalence(
  file.path(input_dir, "pathway_temp.csv"),
  file.path(input_dir, "temperature_pairwise_prevalence.csv"),
  "temp_Category",
  "Temperature"
)

p4 <- plot_prevalence(
  file.path(input_dir, "pathway_salinity.csv"),
  file.path(input_dir, "salinity_pairwise_prevalence.csv"),
  "Salinity_Category",
  "Salinity"
)

# ------------------------------------------------------------------------------
# Combine (Figure 5)
# ------------------------------------------------------------------------------

final_plot <- (p1 | p2) / (p3 | p4)

ggsave(
  filename = file.path(output_dir, "07_prevalence_main_figure.png"),
  plot = final_plot,
  width = 14,
  height = 12,
  dpi = 600
)

ggsave(
  filename = file.path(output_dir, "07_prevalence_main_figure.svg"),
  plot = final_plot,
  width = 14,
  height = 12
)
