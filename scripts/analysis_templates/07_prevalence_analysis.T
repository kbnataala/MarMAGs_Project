################################################################################
# Script: 07_prevalence_analysis.R
#
# Purpose:
#   Documentation template for calculating pathway prevalence and performing
#   pairwise proportion tests across taxonomic and environmental groupings.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

base_dir <- "."
input_file <- file.path(base_dir, "results", "tables", "merged_data.csv")
output_dir <- file.path(base_dir, "results", "tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

data_table <- read.csv(input_file, stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Prevalence calculation
# ------------------------------------------------------------------------------

calculate_prevalence <- function(data, group_col, pathway_col = "CBB") {
  data %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      Number_of_taxa = n(),
      Pathway_presence = sum(.data[[pathway_col]], na.rm = TRUE),
      Prevalence = Pathway_presence / Number_of_taxa,
      Pathway = pathway_col,
      .groups = "drop"
    )
}

# ------------------------------------------------------------------------------
# Pairwise comparison (generalized)
# ------------------------------------------------------------------------------

run_pairwise_tests <- function(prev_data, group_col, output_name) {

  groups <- unique(prev_data[[group_col]])
  pairs <- combn(groups, 2, simplify = FALSE)

  results <- lapply(pairs, function(pair) {

    g1 <- pair[1]
    g2 <- pair[2]

    row1 <- prev_data %>% filter(.data[[group_col]] == g1)
    row2 <- prev_data %>% filter(.data[[group_col]] == g2)

    if (nrow(row1) == 0 || nrow(row2) == 0) return(NULL)

    x <- c(row1$Pathway_presence, row2$Pathway_presence)
    n <- c(row1$Number_of_taxa, row2$Number_of_taxa)

    test <- prop.test(x, n, correct = FALSE)

    data.frame(
      Pathway = unique(prev_data$Pathway),
      Comparison = paste(g1, "vs", g2),
      Prevalence_1 = x[1]/n[1],
      Prevalence_2 = x[2]/n[2],
      p_value = test$p.value,
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    )
  })

  results_df <- bind_rows(results)

  write.csv(
    results_df,
    file.path(output_dir, paste0(output_name, "_pairwise_prevalence.csv")),
    row.names = FALSE
  )

  return(results_df)
}

# ------------------------------------------------------------------------------
# Run analyses (examples)
# ------------------------------------------------------------------------------

phylum_prev <- calculate_prevalence(data_table, "GTDB_Tk_Phylum")
write.csv(phylum_prev, file.path(output_dir, "pathway_phylum.csv"), row.names = FALSE)

run_pairwise_tests(phylum_prev, "GTDB_Tk_Phylum", "phylum")

depth_prev <- calculate_prevalence(data_table, "Depth_Category")
write.csv(depth_prev, file.path(output_dir, "pathway_depth.csv"), row.names = FALSE)

run_pairwise_tests(depth_prev, "Depth_Category", "depth")

sal_prev <- calculate_prevalence(data_table, "Salinity_Category")
write.csv(sal_prev, file.path(output_dir, "pathway_salinity.csv"), row.names = FALSE)

run_pairwise_tests(sal_prev, "Salinity_Category", "salinity")

temp_prev <- calculate_prevalence(data_table, "temp_Category")
write.csv(temp_prev, file.path(output_dir, "pathway_temp.csv"), row.names = FALSE)

run_pairwise_tests(temp_prev, "temp_Category", "temperature")
