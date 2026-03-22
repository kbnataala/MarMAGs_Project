################################################################################
# Script: 11_nri_nti_stats_by_gradient.R
#
# Purpose:
#   Documentation template for computing paired and independent t-test summaries
#   for NRI and NTI across environmental or geographic categories.
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution.
#
# Expected input:
#   - results/tables/<category_table>.csv
#
# Expected outputs:
#   - results/tables/<prefix>_independent_t_test_results.csv
#   - results/tables/<prefix>_paired_t_test_results.csv
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

base_dir <- "."
input_dir <- file.path(base_dir, "results", "tables")
output_dir <- input_dir

student_ttest <- function(mean1, sd1, n1, mean2, sd2, n2) {
  pooled_var <- (((n1 - 1) * sd1^2) + ((n2 - 1) * sd2^2)) / (n1 + n2 - 2)
  se <- sqrt(pooled_var * (1 / n1 + 1 / n2))
  t_val <- (mean1 - mean2) / se
  df <- n1 + n2 - 2
  p_val <- 2 * pt(-abs(t_val), df)
  crit_t <- qt(0.975, df)

  c(t_value = t_val, p_value = p_val, df = df, critical_t = crit_t)
}

run_independent_tests <- function(data, category_col, prefix) {
  data <- data %>%
    mutate(metrics = toupper(metrics))

  nri_data <- subset(data, metrics == "NRI")
  nti_data <- subset(data, metrics == "NTI")

  compare_groups <- function(metric_data, metric_name) {
    levels_vec <- unique(metric_data[[category_col]])
    result_list <- list()

    for (i in 1:(length(levels_vec) - 1)) {
      for (j in (i + 1):length(levels_vec)) {
        g1 <- subset(metric_data, .data[[category_col]] == levels_vec[i])
        g2 <- subset(metric_data, .data[[category_col]] == levels_vec[j])

        stats <- student_ttest(
          g1$nri_and_nti, g1$std_dev, g1$number,
          g2$nri_and_nti, g2$std_dev, g2$number
        )

        result_list[[length(result_list) + 1]] <- data.frame(
          Metric = metric_name,
          Group1 = levels_vec[i],
          Group2 = levels_vec[j],
          t_value = round(stats["t_value"], 3),
          p_value = signif(stats["p_value"], 3),
          df = stats["df"],
          critical_t = round(stats["critical_t"], 2)
        )
      }
    }
    bind_rows(result_list)
  }

  final_results <- bind_rows(
    compare_groups(nri_data, "NRI"),
    compare_groups(nti_data, "NTI")
  )

  write.csv(
    final_results,
    file.path(output_dir, paste0(prefix, "_independent_t_test_results.csv")),
    row.names = FALSE
  )

  final_results
}

run_paired_tests <- function(data, category_col, prefix) {
  data <- data %>%
    mutate(metrics = toupper(metrics))

  nri_data <- subset(data, metrics == "NRI")
  nti_data <- subset(data, metrics == "NTI")

  nri_data <- nri_data[order(nri_data[[category_col]]), ]
  nti_data <- nti_data[order(nti_data[[category_col]]), ]

  paired_ttest <- function(nri, nti, sd_nri, sd_nti, n, category_value) {
    diff <- nti - nri
    sd_diff <- sqrt(sd_nri^2 + sd_nti^2)
    se_diff <- sd_diff / sqrt(n)
    t_val <- diff / se_diff
    df <- n - 1
    p_val <- 2 * pt(-abs(t_val), df)
    crit_t <- qt(0.975, df)

    data.frame(
      Category = category_value,
      Mean_NRI = nri,
      Mean_NTI = nti,
      Mean_Difference = diff,
      t_value = round(t_val, 3),
      p_value = signif(p_val, 3),
      df = df,
      critical_t = round(crit_t, 2)
    )
  }

  paired_results <- mapply(
    paired_ttest,
    nri = nri_data$nri_and_nti,
    nti = nti_data$nri_and_nti,
    sd_nri = nri_data$std_dev,
    sd_nti = nti_data$std_dev,
    n = nri_data$number,
    category_value = nri_data[[category_col]],
    SIMPLIFY = FALSE
  )

  paired_results_df <- bind_rows(paired_results)

  write.csv(
    paired_results_df,
    file.path(output_dir, paste0(prefix, "_paired_t_test_results.csv")),
    row.names = FALSE
  )

  paired_results_df
}

run_nri_nti_tests <- function(input_csv, category_col, prefix) {
  input_file <- file.path(input_dir, input_csv)
  if (!file.exists(input_file)) stop("Missing input file: ", input_file)

  data <- read.csv(input_file, stringsAsFactors = FALSE)
  names(data) <- tolower(gsub("\\.", "_", names(data)))

  required_cols <- c(category_col, "metrics", "nri_and_nti", "std_dev", "number")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", input_csv, ": ",
         paste(missing_cols, collapse = ", "))
  }

  data$number <- as.numeric(as.character(data$number))
  data$nri_and_nti <- as.numeric(as.character(data$nri_and_nti))
  data$std_dev <- as.numeric(as.character(data$std_dev))

  list(
    independent = run_independent_tests(data, category_col, prefix),
    paired = run_paired_tests(data, category_col, prefix)
  )
}

# Examples
# run_nri_nti_tests("temperature_categories.csv", "temperature_categories", "temperature")
# run_nri_nti_tests("depth_categories.csv", "depth_categories", "depth")
# run_nri_nti_tests("salinity_categories.csv", "salinity", "salinity")
# run_nri_nti_tests("ocean.csv", "ocean", "ocean")
# run_nri_nti_tests("sea.csv", "sea", "sea")
