################################################################################
# Script: 05_multivariate_biogeography.R
#
# Purpose:
#   Documentation template for testing how environmental and geographic
#   gradients are associated with phylum-level presence/absence patterns using
#   multivariate binomial GLMs (mvabund::manyglm).
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution. It is not presented as a fully reproducible
#   pipeline component.
#
# Expected inputs:
#   - results/tables/merged_data.csv
#   - results/tables/bin_libs.csv
#
# Expected outputs:
#   - results/tables/depth_univariate_results.csv
#   - results/tables/salinity_univariate_results.csv
#   - results/tables/temp_univariate_results.csv
#   - results/tables/ocean_univariate_results.csv
#   - results/tables/sea_univariate_results.csv
#   - results/tables/*_anova_summary.xlsx
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(mvabund)
  library(writexl)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

base_dir <- "."
input_dir <- file.path(base_dir, "results", "tables")
output_dir <- input_dir

merged_file <- file.path(input_dir, "merged_data.csv")
bin_libs_file <- file.path(input_dir, "bin_libs.csv")

# ------------------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------------------

if (!file.exists(merged_file)) {
  stop("Missing input file: ", merged_file)
}

if (!file.exists(bin_libs_file)) {
  stop("Missing input file: ", bin_libs_file)
}

merged_data <- read.csv(merged_file, header = TRUE, stringsAsFactors = FALSE)
bin_libs <- read.csv(bin_libs_file, header = TRUE, stringsAsFactors = FALSE)

required_merged_cols <- c("Genome", "GTDB_Tk_Phylum")
required_binlib_cols <- c("bin", "library")

missing_merged <- setdiff(required_merged_cols, colnames(merged_data))
missing_binlib <- setdiff(required_binlib_cols, colnames(bin_libs))

if (length(missing_merged) > 0) {
  stop("Missing required columns in merged_data.csv: ",
       paste(missing_merged, collapse = ", "))
}

if (length(missing_binlib) > 0) {
  stop("Missing required columns in bin_libs.csv: ",
       paste(missing_binlib, collapse = ", "))
}

# Merge once, use everywhere
merged_data <- merge(
  merged_data,
  bin_libs,
  by.x = "Genome",
  by.y = "bin",
  all.x = TRUE
)

# ------------------------------------------------------------------------------
# Helper function
# ------------------------------------------------------------------------------

run_manyglm_by_gradient <- function(data,
                                    gradient_col,
                                    output_prefix,
                                    remove_unknown = FALSE) {

  required_cols <- c("library", "GTDB_Tk_Phylum", gradient_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns for gradient '", gradient_col, "': ",
         paste(missing_cols, collapse = ", "))
  }

  subset_df <- data %>%
    select(all_of(required_cols)) %>%
    filter(
      !is.na(library),
      !is.na(GTDB_Tk_Phylum),
      !is.na(.data[[gradient_col]])
    )

  if (remove_unknown) {
    subset_df <- subset_df %>%
      filter(
        library != "Unknown",
        GTDB_Tk_Phylum != "Unknown",
        .data[[gradient_col]] != "Unknown"
      )
  }

  # Build sample-by-phyla presence matrix
  presence_matrix <- subset_df %>%
    distinct(library, GTDB_Tk_Phylum) %>%
    mutate(presence = 1) %>%
    pivot_wider(
      names_from = GTDB_Tk_Phylum,
      values_from = presence,
      values_fill = 0
    )

  # Add metadata
  metadata <- subset_df %>%
    distinct(library, .data[[gradient_col]])

  full_data <- presence_matrix %>%
    left_join(metadata, by = "library") %>%
    column_to_rownames("library")

  comm_data <- full_data %>% select(-all_of(gradient_col))
  gradient_group <- full_data[[gradient_col]]

  # mvabund analysis
  comm_mv <- mvabund(comm_data)
  mv_model <- manyglm(comm_mv ~ gradient_group, family = "binomial")

  result <- anova.manyglm(
    mv_model,
    p.uni = "adjusted",
    resamp = "montecarlo",
    show.warning = TRUE
  )

  print(result)

  # Extract univariate results
  test_stat <- t(result$uni.test["gradient_group", ])
  p_val <- t(result$uni.p["gradient_group", ])

  univariate_df <- data.frame(
    Taxon = colnames(result$uni.test),
    Deviance = as.numeric(test_stat),
    P_value = as.numeric(p_val)
  ) %>%
    arrange(P_value)

  # Save outputs
  write.csv(
    univariate_df,
    file.path(output_dir, paste0(output_prefix, "_univariate_results.csv")),
    row.names = FALSE
  )

  write_xlsx(
    list(
      Multivariate = as.data.frame(result$table),
      Univariate = univariate_df
    ),
    file.path(output_dir, paste0(output_prefix, "_anova_summary.xlsx"))
  )

  invisible(list(
    model = mv_model,
    result = result,
    univariate = univariate_df
  ))
}

# ------------------------------------------------------------------------------
# Run analyses
# ------------------------------------------------------------------------------

depth_results <- run_manyglm_by_gradient(
  data = merged_data,
  gradient_col = "Depth_Category",
  output_prefix = "depth",
  remove_unknown = FALSE
)

salinity_results <- run_manyglm_by_gradient(
  data = merged_data,
  gradient_col = "Salinity_Category",
  output_prefix = "salinity",
  remove_unknown = FALSE
)

temp_results <- run_manyglm_by_gradient(
  data = merged_data,
  gradient_col = "temp_Category",
  output_prefix = "temp",
  remove_unknown = FALSE
)

ocean_results <- run_manyglm_by_gradient(
  data = merged_data,
  gradient_col = "ocean_sea_name",
  output_prefix = "ocean",
  remove_unknown = TRUE
)

sea_results <- run_manyglm_by_gradient(
  data = merged_data,
  gradient_col = "IHO_Sea",
  output_prefix = "sea",
  remove_unknown = TRUE
)
