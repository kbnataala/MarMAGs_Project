################################################################################
# Script: 10_tree_reduction_sensitivity_nri_nti.R
#
# Purpose:
#   Documentation template for evaluating robustness of NRI and NTI to tree-size
#   reduction. The script summarizes:
#     1. Standardized slopes of metric value vs tree size
#     2. Rank-preservation relative to the full tree
#     3. Faceted regression panels by region
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution. It is not presented as a fully reproducible
#   pipeline component.
#
# Expected inputs:
#   - results/tables/NRI_sea.csv
#   - results/tables/NTI_sea.csv
#
# Expected outputs:
#   - results/supplementary_tables/MarMAGs_additional_File_04.docx(Fig. S31)
#   - results/supplementary_tables/MarMAGs_additional_File_04.docx(Fig. S32)
#   - results/supplementary_tables/MarMAGs_additional_File_04.docx(Fig. S33)
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(stringr)
  library(scales)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# Paths and settings
# ------------------------------------------------------------------------------

base_dir <- "."
input_dir <- file.path(base_dir, "results", "tables")
figure_dir <- file.path(base_dir, "results", "figures")
table_dir <- input_dir

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  file.path(input_dir, "NRI_sea.csv"),
  file.path(input_dir, "NTI_sea.csv")
)

size_col <- "Tree_percentage"
size_label <- "Tree size (%)"
delta_std <- 0.20

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing input files: ", paste(missing_files, collapse = ", "))
}

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

parse_metric <- function(path) {
  nm <- basename(path)
  if (str_detect(nm, regex("NTI", ignore_case = TRUE))) "NTI" else "NRI"
}

read_and_pivot <- function(path, size_col) {
  dat <- readr::read_csv(path, show_col_types = FALSE)

  if (!size_col %in% names(dat)) {
    stop("Column '", size_col, "' not found in ", path)
  }

  dat %>%
    rename(TreeSize = !!sym(size_col)) %>%
    pivot_longer(-TreeSize, names_to = "Region", values_to = "Value") %>%
    mutate(
      Metric = parse_metric(path),
      SourceFile = basename(path)
    )
}

theme_clean <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      legend.title = element_text(face = "bold")
    )
}

# ------------------------------------------------------------------------------
# Read and tidy data
# ------------------------------------------------------------------------------

all_long <- purrr::map_df(files, read_and_pivot, size_col = size_col) %>%
  mutate(
    Region = Region %>%
      str_replace_all("[.]+", "_") %>%
      str_replace_all("__+", "_")
  ) %>%
  drop_na(TreeSize, Value)

if (nrow(all_long) == 0) {
  stop("No valid rows found after reshaping input files.")
}

# ------------------------------------------------------------------------------
# Core regression summaries
# ------------------------------------------------------------------------------

reg <- all_long %>%
  group_by(Metric, Region) %>%
  group_modify(~{
    df <- .x %>% drop_na(TreeSize, Value)

    if (nrow(df) < 3) {
      return(tibble(
        n = nrow(df),
        slope_raw = NA_real_,
        slope_se = NA_real_,
        p_value = NA_real_,
        r2 = NA_real_,
        beta_std = NA_real_,
        se_std = NA_real_,
        ci95_lo_std = NA_real_,
        ci95_hi_std = NA_real_,
        TOST_std_equiv = NA
      ))
    }

    m <- lm(Value ~ TreeSize, data = df)
    mt <- tidy(m)
    mg <- glance(m)
    s <- mt %>% filter(term == "TreeSize")

    sd_x <- sd(df$TreeSize, na.rm = TRUE)
    sd_y <- sd(df$Value, na.rm = TRUE)

    beta_std <- if (is.finite(sd_x) && is.finite(sd_y) && sd_y > 0) {
      s$estimate * sd_x / sd_y
    } else {
      NA_real_
    }

    se_std <- if (is.finite(sd_x) && is.finite(sd_y) && sd_y > 0) {
      s$std.error * sd_x / sd_y
    } else {
      NA_real_
    }

    n <- nrow(df)

    if (is.finite(beta_std) && is.finite(se_std) && se_std > 0 && n > 2) {
      dfree <- n - 2
      p1 <- 1 - pt((beta_std + delta_std) / se_std, df = dfree)
      p2 <- 1 - pt((delta_std - beta_std) / se_std, df = dfree)
      equiv <- (p1 < 0.05) & (p2 < 0.05)
    } else {
      equiv <- NA
    }

    tibble(
      n = n,
      slope_raw = s$estimate,
      slope_se = s$std.error,
      p_value = s$p.value,
      r2 = mg$r.squared,
      beta_std = beta_std,
      se_std = se_std,
      ci95_lo_std = beta_std - 1.96 * se_std,
      ci95_hi_std = beta_std + 1.96 * se_std,
      TOST_std_equiv = equiv
    )
  }) %>%
  ungroup() %>%
  group_by(Metric) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# ------------------------------------------------------------------------------
# Rank preservation relative to full tree
# ------------------------------------------------------------------------------

rank_pres <- all_long %>%
  group_by(Metric) %>%
  group_modify(~{
    dat <- .x
    base <- dat %>%
      filter(TreeSize == 100) %>%
      select(Region, Value100 = Value)

    dat %>%
      filter(TreeSize != 100) %>%
      group_by(TreeSize) %>%
      group_modify(~{
        cur <- left_join(.x %>% select(Region, Value), base, by = "Region") %>%
          drop_na()

        tibble(
          rho = if (nrow(cur) >= 3) {
            suppressWarnings(cor(cur$Value, cur$Value100, method = "spearman"))
          } else {
            NA_real_
          }
        )
      }) %>%
      ungroup()
  }) %>%
  ungroup()

# ------------------------------------------------------------------------------
# Figure 1: Forest plot of standardized slopes
# ------------------------------------------------------------------------------

reg_for_plot <- reg %>%
  mutate(Region = forcats::fct_reorder(Region, beta_std, .fun = median, na.rm = TRUE))

p_forest <- ggplot(
  reg_for_plot,
  aes(y = Region, x = beta_std, xmin = ci95_lo_std, xmax = ci95_hi_std, shape = Metric)
) +
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = -delta_std, xmax = delta_std, alpha = 0.08) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  labs(
    x = expression(Standardized~slope~(beta[std])),
    y = NULL,
    title = "Effect of tree size on NRI and NTI"
  ) +
  theme_clean() +
  theme(legend.position = "right")

# ------------------------------------------------------------------------------
# Figure 2: Rank preservation curves
# ------------------------------------------------------------------------------

p_rank <- rank_pres %>%
  ggplot(aes(x = as.numeric(TreeSize), y = rho, color = Metric, group = Metric)) +
  geom_hline(yintercept = 1, linetype = 3) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = size_label,
    y = expression(Spearman~rho),
    title = "Rank preservation relative to full tree"
  ) +
  theme_clean() +
  theme(legend.position = "right")

# ------------------------------------------------------------------------------
# Figure 3: Faceted regression panels
# ------------------------------------------------------------------------------

make_panel_plot <- function(metric_label) {
  df_m <- all_long %>% filter(Metric == metric_label)

  ann_m <- reg %>%
    filter(Metric == metric_label) %>%
    mutate(
      lab = sprintf(
        "β = %.3f\nR² = %.2f\nFDR-p = %.2g\nTOST(±%.2f): %s",
        slope_raw,
        r2,
        p_adj_BH,
        delta_std,
        ifelse(isTRUE(TOST_std_equiv), "equivalent", "not equivalent")
      )
    )

  ann_pos <- df_m %>%
    group_by(Region) %>%
    summarise(
      x_pos = min(TreeSize, na.rm = TRUE) + 0.02 * diff(range(TreeSize, na.rm = TRUE)),
      y_pos = min(Value, na.rm = TRUE) + 0.06 * diff(range(Value, na.rm = TRUE)),
      .groups = "drop"
    )

  ann_m <- left_join(ann_m, ann_pos, by = "Region")

  ggplot(df_m, aes(TreeSize, Value)) +
    geom_point(size = 1.2, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
    facet_wrap(~Region, scales = "free_y") +
    geom_text(
      data = ann_m,
      aes(x = x_pos, y = y_pos, label = lab),
      hjust = 0,
      vjust = 0,
      size = 3.2,
      lineheight = 1
    ) +
    labs(
      x = size_label,
      y = metric_label,
      title = paste(metric_label, "vs tree size")
    ) +
    theme_clean()
}

p_nri <- make_panel_plot("NRI")
p_nti <- make_panel_plot("NTI")
fig_panels <- (p_nri / p_nti) + plot_annotation(tag_levels = "a")

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------

ggsave(file.path(figure_dir, "10_forest_std_slopes.svg"), p_forest, width = 6, height = 6)
ggsave(file.path(figure_dir, "10_rank_preservation.svg"), p_rank, width = 6, height = 4)

n_regions <- all_long %>% distinct(Region) %>% nrow()
rows_est <- ceiling(n_regions / 3)
fig3_height <- max(5, rows_est * 2 * 1.8)

ggsave(file.path(figure_dir, "10_regression_panels.svg"), fig_panels, width = 8, height = fig3_height)

readr::write_csv(reg, file.path(table_dir, "10_core_regression_standardized.csv"))
readr::write_csv(rank_pres, file.path(table_dir, "10_rank_preservation.csv"))

message("Saved tree-reduction sensitivity figures and tables.")
