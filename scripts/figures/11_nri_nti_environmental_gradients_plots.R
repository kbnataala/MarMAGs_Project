################################################################################
# Script: 11_nri_nti_environmental_gradients_plots.R
#
# Purpose:
#   Documentation template for plotting divergent NRI/NTI barplots with
#   paired-test stars and optional compact letter displays (CLD).
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution.
#
# Expected inputs:
#   - results/tables/<category_table>.csv
#   - results/tables/<prefix>_paired_t_test_results.csv
#   - results/tables/<prefix>_independent_t_test_results.csv
#
# Expected outputs:
#   - results/figures/11_<prefix>_nri_nti.svg
#   - results/figures/11_merged_nri_nti_panels.svg
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(multcompView)
  library(patchwork)
})

base_dir <- "."
input_dir <- file.path(base_dir, "results", "tables")
output_dir <- file.path(base_dir, "results", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

p_to_stars <- function(p) {
  ifelse(is.na(p), "ns",
         ifelse(p <= 0.001, "***",
                ifelse(p <= 0.01, "**",
                       ifelse(p <= 0.05, "*", "ns"))))
}

cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
  pcol <- intersect(c("p_value", "p.value", "p-value"), names(indep_df))[1]
  if (is.na(pcol)) return(setNames(rep("", length(cats)), cats))

  sub <- indep_df %>%
    filter(Metric == metric, Group1 %in% cats, Group2 %in% cats)

  if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))

  pv <- sub[[pcol]]
  names(pv) <- paste(pmin(sub$Group1, sub$Group2), pmax(sub$Group1, sub$Group2), sep = "-")
  cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters

  out <- setNames(rep("", length(cats)), cats)
  out[names(cld)] <- cld
  out
}

plot_divergent_nri_nti <- function(
  data_csv,
  paired_csv = NULL,
  indep_csv = NULL,
  category_col,
  x_label,
  category_levels = NULL,
  output_prefix = NULL
) {
  data_file <- file.path(input_dir, data_csv)
  if (!file.exists(data_file)) stop("Missing input file: ", data_file)

  df <- read.csv(data_file, stringsAsFactors = FALSE)
  names(df) <- tolower(gsub("\\.", "_", names(df)))
  category_col <- tolower(category_col)

  required <- c(category_col, "metrics", "nri_and_nti", "std_dev")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns in ", data_csv, ": ",
         paste(missing, collapse = ", "))
  }

  df$metrics <- toupper(df$metrics)

  if (!is.null(category_levels)) {
    df[[category_col]] <- factor(df[[category_col]], levels = category_levels)
  } else {
    df[[category_col]] <- factor(df[[category_col]])
  }

  paired_ann <- NULL
  if (!is.null(paired_csv)) {
    paired_file <- file.path(input_dir, paired_csv)
    if (file.exists(paired_file)) {
      paired <- read.csv(paired_file, stringsAsFactors = FALSE)
      names(paired) <- tolower(gsub("\\.", "_", names(paired)))
      cat_col2 <- intersect(c("category", category_col), names(paired))[1]
      pcol <- intersect(c("p_value", "p.value", "p-value"), names(paired))[1]

      if (!is.na(cat_col2) && !is.na(pcol)) {
        names(paired)[names(paired) == cat_col2] <- "category"
        paired_ann <- paired %>%
          mutate(stars = p_to_stars(.data[[pcol]]))
      }
    }
  }

  p <- ggplot(df, aes(x = .data[[category_col]], y = nri_and_nti, fill = metrics)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.65) +
    geom_errorbar(
      aes(ymin = nri_and_nti - std_dev, ymax = nri_and_nti + std_dev),
      position = position_dodge(width = 0.6),
      width = 0.2
    ) +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7") +
    scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +
    coord_flip() +
    labs(x = x_label, y = "NRI and NTI", fill = "Metrics") +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )

  if (!is.null(paired_ann)) {
    pos_df <- df %>%
      group_by(.data[[category_col]]) %>%
      summarise(
        top = max(nri_and_nti + ifelse(nri_and_nti >= 0, std_dev, 0), na.rm = TRUE),
        bot = min(nri_and_nti - ifelse(nri_and_nti < 0, std_dev, 0), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(y_pos = ifelse(top <= 0, bot - 0.2, top + 0.2))

    names(pos_df)[1] <- "category"
    paired_ann <- left_join(paired_ann, pos_df, by = "category")

    p <- p +
      geom_text(
        data = paired_ann,
        aes(x = category, y = y_pos, label = stars),
        inherit.aes = FALSE,
        fontface = "bold",
        size = 5
      )
  }

  if (!is.null(output_prefix)) {
    ggsave(
      file.path(output_dir, paste0("11_", output_prefix, "_nri_nti.svg")),
      p,
      width = 10,
      height = 10
    )
  }

  p
}

# Examples
plot_temp <- plot_divergent_nri_nti(
  data_csv = "temperature_categories.csv",
  paired_csv = "temperature_paired_t_test_results.csv",
  indep_csv = "temperature_independent_t_test_results.csv",
  category_col = "temperature_categories",
  x_label = "Temperature",
  output_prefix = "temperature"
)

plot_depth <- plot_divergent_nri_nti(
  data_csv = "depth_categories.csv",
  paired_csv = "depth_paired_t_test_results.csv",
  indep_csv = "depth_independent_t_test_results.csv",
  category_col = "depth_categories",
  x_label = "Depth",
  output_prefix = "depth"
)

plot_salinity <- plot_divergent_nri_nti(
  data_csv = "salinity_categories.csv",
  paired_csv = "salinity_paired_t_test_results.csv",
  indep_csv = "salinity_independent_t_test_results.csv",
  category_col = "salinity",
  x_label = "Salinity",
  output_prefix = "salinity"
)

plot_ocean <- plot_divergent_nri_nti(
  data_csv = "ocean.csv",
  paired_csv = "ocean_paired_t_test_results.csv",
  indep_csv = "ocean_independent_t_test_results.csv",
  category_col = "ocean",
  x_label = "Ocean / Sea",
  output_prefix = "ocean"
)

merged_plot <- (plot_temp | plot_depth) / (plot_salinity | plot_ocean)

ggsave(
  file.path(output_dir, "11_merged_nri_nti_panels.svg"),
  merged_plot,
  width = 14,
  height = 14
)
