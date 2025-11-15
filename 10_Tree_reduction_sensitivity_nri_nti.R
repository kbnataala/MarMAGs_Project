#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/benchmarking_nti_nri_analysis")

# NTI/NRI vs Tree Size — Robustness to Reduced Trees (Nature-style)
# Author: <you>
#
# Streamlined Analysis — Three Core Outputs Only (Final)
# 1) Forest plot of STANDARDIZED slopes (NTI vs NRI) with equivalence band (±Δ_std)
# 2) Rank-preservation curves (Spearman rho) vs % retained
# 3) Faceted regression panels (scatter + OLS) for NRI and NTI
#
# All other analyses/figures have been removed.

# -------------------------
# Load packages
# -------------------------
packages <- c("tidyverse", "broom", "stringr", "scales", "patchwork")
for (p in packages) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(tidyverse)
library(broom)
library(stringr)
library(scales)
library(patchwork)

# Optional: fonts for a Nature-like look
if (!requireNamespace("showtext", quietly = TRUE)) install.packages("showtext")
if (!requireNamespace("systemfonts", quietly = TRUE)) install.packages("systemfonts")
library(showtext)
library(systemfonts)
showtext_auto(enable = TRUE)
base_family <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "sans"

# -------------------------
# SETTINGS — edit these to match your files
# -------------------------
files <- c(
  "NRI_sea.csv",
  "NTI_sea.csv"
  # add more here if needed
)
size_col   <- "Tree_percentage"     # column name for % of tree retained
size_label <- "Tree size (%)"

# Equivalence bound in STANDARDIZED slope units (small effect)
DELTA_STD <- 0.20

# Output directory
out_dir <- "figures_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------
# Helpers
# -------------------------
parse_metric <- function(path) {
  nm <- basename(path)
  if (str_detect(nm, regex("NTI", ignore_case = TRUE))) "NTI" else "NRI"
}

read_and_pivot <- function(path, size_col) {
  dat <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot(size_col %in% names(dat))
  dat %>%
    rename(TreeSize = !!sym(size_col)) %>%
    pivot_longer(-TreeSize, names_to = "Region", values_to = "Value") %>%
    mutate(Metric = parse_metric(path), SourceFile = basename(path))
}

# Minimal Nature-style theme
theme_nature <- function(base_size = 8, font_family = base_family) {
  theme_classic(base_size = base_size, base_family = font_family) %+replace%
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.ticks.length = unit(2, "pt"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# -------------------------
# Read & tidy
# -------------------------
all_long <- purrr::map_df(files, read_and_pivot, size_col = size_col) %>%
  mutate(Region = Region %>% str_replace_all("[.]+", "_") %>% str_replace_all("__+", "_")) %>%
  drop_na(TreeSize, Value)
stopifnot(nrow(all_long) > 0)

# -------------------------
# Core computations
# -------------------------
# (A) Per-panel regression: raw and STANDARDIZED slopes, FDR, TOST (standardized)
reg <- all_long %>%
  group_by(Metric, Region) %>%
  group_modify(~{
    df <- .x %>% drop_na(TreeSize, Value)
    m  <- lm(Value ~ TreeSize, data = df)
    mt <- tidy(m)
    mg <- glance(m)
    s  <- mt %>% filter(term == "TreeSize")
    # standardize slope: beta_std = beta_raw * SD(x) / SD(y)
    sd_x <- sd(df$TreeSize, na.rm = TRUE)
    sd_y <- sd(df$Value,    na.rm = TRUE)
    beta_std <- if (is.finite(sd_x) && is.finite(sd_y) && sd_y > 0) s$estimate * sd_x / sd_y else NA_real_
    se_std   <- if (is.finite(sd_x) && is.finite(sd_y) && sd_y > 0) s$std.error * sd_x / sd_y else NA_real_
    n <- nrow(df)
    # TOST in standardized space
    if (is.finite(beta_std) && is.finite(se_std) && se_std > 0 && n > 2) {
      dfree <- n - 2
      p1 <- 1 - pt((beta_std + DELTA_STD)/se_std, df = dfree)
      p2 <- 1 - pt((DELTA_STD - beta_std)/se_std, df = dfree)
      equiv <- (p1 < 0.05) & (p2 < 0.05)
    } else { equiv <- NA }
    tibble(n = n,
           slope_raw = s$estimate, slope_se = s$std.error, p_value = s$p.value,
           r2 = mg$r.squared,
           beta_std = as.numeric(beta_std), se_std = as.numeric(se_std),
           ci95_lo_std = beta_std - 1.96*se_std,
           ci95_hi_std = beta_std + 1.96*se_std,
           TOST_std_equiv = equiv)
  }) %>% ungroup() %>%
  group_by(Metric) %>% mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>% ungroup()

# (B) Rank preservation vs 100% per metric
rank_pres <- all_long %>%
  group_by(Metric) %>%
  group_modify(~{
    dat <- .x
    base <- dat %>% filter(TreeSize == 100) %>% select(Region, Value100 = Value)
    dat %>% filter(TreeSize != 100) %>%
      group_by(TreeSize) %>%
      group_modify(~{
        cur <- left_join(.x %>% select(Region, Value), base, by = "Region") %>% drop_na()
        tibble(rho = if (nrow(cur) >= 3) suppressWarnings(cor(cur$Value, cur$Value100, method = "spearman")) else NA_real_)
      }) %>% ungroup()
  }) %>% ungroup()

# -------------------------
# FIGURES (three only)
# -------------------------
# 1) Forest of standardized slopes (NTI vs NRI), with equivalence band ±DELTA_STD
reg_for_plot <- reg %>% mutate(Region = forcats::fct_reorder(Region, beta_std, .fun = median, na.rm = TRUE))

p_forest <- ggplot(reg_for_plot, aes(y = Region, x = beta_std, xmin = ci95_lo_std, xmax = ci95_hi_std, shape = Metric)) +
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = -DELTA_STD, xmax =  DELTA_STD, alpha = 0.08) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  scale_shape_discrete(name = "Metric") +
  labs(x = expression(Standardized~slope~(beta[std])~per~SD(Tree~size)), y = NULL,
       title = "") +
  theme_nature() +
  theme(legend.position = "right")

# 2) Rank-preservation curves
p_rank <- rank_pres %>%
  ggplot(aes(x = as.numeric(TreeSize), y = rho, color = Metric, group = Metric)) +
  geom_hline(yintercept = 1, linetype = 3) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_color_discrete(name = "Metric") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = size_label, y = expression(Spearman~rho),
       title = "") +
  theme_nature() +
  theme(legend.position = "right")

# 3) Faceted regression panels (scatter + OLS) for NRI and NTI
make_panel_plot <- function(metric_label) {
  df_m <- all_long %>% filter(Metric == metric_label)
  ann_m <- reg %>% filter(Metric == metric_label) %>%
    mutate(lab = sprintf("%s
%s
%s
%s",
                         sprintf("β = %.3f", slope_raw),
                         sprintf("R² = %.2f", r2),
                         sprintf("FDR-p = %.2g", p_adj_BH),
                         ifelse(isTRUE(TOST_std_equiv),
                                sprintf("TOST(±%.2f): equivalent", DELTA_STD),
                                sprintf("TOST(±%.2f): not equivalent", DELTA_STD))))
  # place labels bottom-left to avoid overlap with points/line
  ann_pos <- df_m %>% group_by(Region) %>% summarise(
    x_pos = min(TreeSize, na.rm = TRUE) + 0.02 * diff(range(TreeSize, na.rm = TRUE)),
    y_pos = min(Value,    na.rm = TRUE) + 0.06 * diff(range(Value,    na.rm = TRUE)), .groups = "drop")
  ann_m <- left_join(ann_m, ann_pos, by = "Region")
  
  ggplot(df_m, aes(TreeSize, Value)) +
    geom_point(size = 1.2, alpha = 0.85, stroke = 0) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
    facet_wrap(~ Region, scales = "free_y") +
    geom_text(data = ann_m, aes(x = x_pos, y = y_pos, label = lab),
              hjust = 0, vjust = 0, size = 2.6, lineheight = 1.0) +
    labs(x = size_label, y = metric_label) +
    theme_nature()
}


p_nri <- make_panel_plot("NRI") + ggtitle("NRI vs Tree size")
p_nti <- make_panel_plot("NTI") + ggtitle("NTI vs Tree size")

# combine vertically
fig_panels <- (p_nri / p_nti) + plot_annotation(tag_levels = "a")

# -------------------------
# EXPORTS (only the three figures + minimal tables)
# -------------------------
use_cairo <- isTRUE(capabilities("cairo"))
dev_fun <- if (use_cairo) cairo_pdf else grDevices::pdf

# Figures
ggsave(file.path(out_dir, "FIG1_Forest_StdSlopes.pdf"), p_forest, width = 5.4, height = 5.8, device = dev_fun)
ggsave(file.path(out_dir, "FIG1_Forest_StdSlopes.tiff"), p_forest, width = 5.4, height = 5.8, dpi = 600, compression = "lzw")

ggsave(file.path(out_dir, "FIG2_RankPreservation.pdf"), p_rank, width = 5.0, height = 3.6, device = dev_fun)
ggsave(file.path(out_dir, "FIG2_RankPreservation.tiff"), p_rank, width = 5.0, height = 3.6, dpi = 600, compression = "lzw")

# Height for panels depends on number of facets
n_regions <- all_long %>% distinct(Region) %>% nrow()
rows_est  <- ceiling(n_regions / 3)
fig3_height <- max(4.5, rows_est * 1.8 * 2)  # two rows of panels (NRI, NTI)

ggsave(file.path(out_dir, "FIG3_RegressionPanels.pdf"), fig_panels, width = 8.2, height = fig3_height, device = dev_fun)
ggsave(file.path(out_dir, "FIG3_RegressionPanels.tiff"), fig_panels, width = 7.2, height = fig3_height, dpi = 600, compression = "lzw")

# Tables (minimal)
readr::write_csv(reg,       file.path(out_dir, "table_core_regression_standardized.csv"))
readr::write_csv(rank_pres, file.path(out_dir, "table_rank_preservation.csv"))

message("Saved streamlined core figures and tables to ", normalizePath(out_dir))














###############
# Streamlined Analysis — Three Core Outputs Only (Final)
# 1) Forest plot of STANDARDIZED slopes (NTI vs NRI) with equivalence band (±Δ_std)
# 2) Rank-preservation curves (Spearman rho) vs % retained
# 3) Faceted regression panels (scatter + OLS) for NRI and NTI
#
# All other analyses/figures have been removed.

# -------------------------
# Load packages
# -------------------------
packages <- c("tidyverse", "broom", "stringr", "scales", "patchwork")
for (p in packages) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(tidyverse)
library(broom)
library(stringr)
library(scales)
library(patchwork)

# Optional: fonts for a Nature-like look
if (!requireNamespace("showtext", quietly = TRUE)) install.packages("showtext")
if (!requireNamespace("systemfonts", quietly = TRUE)) install.packages("systemfonts")
library(showtext)
library(systemfonts)
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
library(svglite)
showtext_auto(enable = TRUE)
if (exists("showtext_opts")) showtext_opts(dpi = 600)
base_family <- if ("Arial" %in% systemfonts::system_fonts()$family) "Arial" else "sans"
BASE_SIZE <- 12  # increased overall text size

# -------------------------
# SETTINGS — edit these to match your files
# -------------------------
files <- c(
  "NRI_sea.csv",
  "NTI_sea.csv"
  # add more here if needed
)
size_col   <- "Tree_percentage"     # column name for % of tree retained
size_label <- "Tree size (%)"

# Equivalence bound in STANDARDIZED slope units (small effect)
DELTA_STD <- 0.20

# Output directory
out_dir <- "figures_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------
# Helpers
# -------------------------
parse_metric <- function(path) {
  nm <- basename(path)
  if (str_detect(nm, regex("NTI", ignore_case = TRUE))) "NTI" else "NRI"
}

read_and_pivot <- function(path, size_col) {
  dat <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot(size_col %in% names(dat))
  dat %>%
    rename(TreeSize = !!sym(size_col)) %>%
    pivot_longer(-TreeSize, names_to = "Region", values_to = "Value") %>%
    mutate(Metric = parse_metric(path), SourceFile = basename(path))
}

# Minimal Nature-style theme
theme_nature <- function(base_size = BASE_SIZE, font_family = base_family) {
  theme_classic(base_size = base_size, base_family = font_family) %+replace%
    theme(
      axis.text = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size),
      axis.ticks.length = unit(2, "pt"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.position = "none",
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.margin = margin(4, 4, 4, 4)
    )
}

# -------------------------
# Read & tidy
# -------------------------
all_long <- purrr::map_df(files, read_and_pivot, size_col = size_col) %>%
  mutate(Region = Region %>% str_replace_all("[.]+", "_") %>% str_replace_all("__+", "_")) %>%
  drop_na(TreeSize, Value)
stopifnot(nrow(all_long) > 0)

# -------------------------
# Core computations
# -------------------------
# (A) Per-panel regression: raw and STANDARDIZED slopes, FDR, TOST (standardized)
reg <- all_long %>%
  group_by(Metric, Region) %>%
  group_modify(~{
    df <- .x %>% drop_na(TreeSize, Value)
    m  <- lm(Value ~ TreeSize, data = df)
    mt <- tidy(m)
    mg <- glance(m)
    s  <- mt %>% filter(term == "TreeSize")
    # standardize slope: beta_std = beta_raw * SD(x) / SD(y)
    sd_x <- sd(df$TreeSize, na.rm = TRUE)
    sd_y <- sd(df$Value,    na.rm = TRUE)
    beta_std <- if (is.finite(sd_x) && is.finite(sd_y) && sd_y > 0) s$estimate * sd_x / sd_y else NA_real_
    se_std   <- if (is.finite(sd_x) && is.finite(sd_y) && sd_y > 0) s$std.error * sd_x / sd_y else NA_real_
    n <- nrow(df)
    # TOST in standardized space
    if (is.finite(beta_std) && is.finite(se_std) && se_std > 0 && n > 2) {
      dfree <- n - 2
      p1 <- 1 - pt((beta_std + DELTA_STD)/se_std, df = dfree)
      p2 <- 1 - pt((DELTA_STD - beta_std)/se_std, df = dfree)
      equiv <- (p1 < 0.05) & (p2 < 0.05)
    } else { equiv <- NA }
    tibble(n = n,
           slope_raw = s$estimate, slope_se = s$std.error, p_value = s$p.value,
           r2 = mg$r.squared,
           beta_std = as.numeric(beta_std), se_std = as.numeric(se_std),
           ci95_lo_std = beta_std - 1.96*se_std,
           ci95_hi_std = beta_std + 1.96*se_std,
           TOST_std_equiv = equiv)
  }) %>% ungroup() %>%
  group_by(Metric) %>% mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>% ungroup()

# (B) Rank preservation vs 100% per metric
rank_pres <- all_long %>%
  group_by(Metric) %>%
  group_modify(~{
    dat <- .x
    base <- dat %>% filter(TreeSize == 100) %>% select(Region, Value100 = Value)
    dat %>% filter(TreeSize != 100) %>%
      group_by(TreeSize) %>%
      group_modify(~{
        cur <- left_join(.x %>% select(Region, Value), base, by = "Region") %>% drop_na()
        tibble(rho = if (nrow(cur) >= 3) suppressWarnings(cor(cur$Value, cur$Value100, method = "spearman")) else NA_real_)
      }) %>% ungroup()
  }) %>% ungroup()

# -------------------------
# FIGURES (three only)
# -------------------------
# 1) Forest of standardized slopes (NTI vs NRI), with equivalence band ±DELTA_STD
reg_for_plot <- reg %>% mutate(Region = forcats::fct_reorder(Region, beta_std, .fun = median, na.rm = TRUE))

p_forest <- ggplot(reg_for_plot, aes(y = Region, x = beta_std, xmin = ci95_lo_std, xmax = ci95_hi_std, shape = Metric)) +
  annotate("rect", ymin = -Inf, ymax = Inf, xmin = -DELTA_STD, xmax =  DELTA_STD, alpha = 0.08) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  scale_shape_discrete(name = "Metric") +
  labs(x = expression(Standardized~slope~(beta[std])~per~SD(Tree~size)), y = NULL,
       title = "Effect of tree size on NTI vs NRI (standardized slopes)") +
  theme_nature() +
  theme(text = element_text(size = BASE_SIZE + 1),
        legend.position = "right",
        legend.title = element_text(size = BASE_SIZE + 1),
        legend.text  = element_text(size = BASE_SIZE))

# 2) Rank-preservation curves
p_rank <- rank_pres %>%
  ggplot(aes(x = as.numeric(TreeSize), y = rho, color = Metric, group = Metric)) +
  geom_hline(yintercept = 1, linetype = 3) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
  scale_color_discrete(name = "Metric") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = size_label, y = expression(Spearman~rho),
       title = "Between-region rank preservation vs full tree") +
  theme_nature() +
  theme(text = element_text(size = BASE_SIZE + 1),
        legend.position = "right",
        legend.title = element_text(size = BASE_SIZE + 1),
        legend.text  = element_text(size = BASE_SIZE))

# 3) Faceted regression panels (scatter + OLS) for NRI and NTI
make_panel_plot <- function(metric_label) {
  df_m <- all_long %>% filter(Metric == metric_label)
  ann_m <- reg %>% filter(Metric == metric_label) %>%
    mutate(lab = sprintf("%s
%s
%s
%s",
                         sprintf("β = %.3f", slope_raw),
                         sprintf("R² = %.2f", r2),
                         sprintf("FDR-p = %.2g", p_adj_BH),
                         ifelse(isTRUE(TOST_std_equiv),
                                sprintf("TOST(±%.2f): equivalent", DELTA_STD),
                                sprintf("TOST(±%.2f): not equivalent", DELTA_STD))))
  # place labels bottom-left to avoid overlap with points/line
  ann_pos <- df_m %>% group_by(Region) %>% summarise(
    x_pos = min(TreeSize, na.rm = TRUE) + 0.02 * diff(range(TreeSize, na.rm = TRUE)),
    y_pos = min(Value,    na.rm = TRUE) + 0.06 * diff(range(Value,    na.rm = TRUE)), .groups = "drop")
  ann_m <- left_join(ann_m, ann_pos, by = "Region")
  
  ggplot(df_m, aes(TreeSize, Value)) +
    geom_point(size = 1.2, alpha = 0.85, stroke = 0) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
    facet_wrap(~ Region, scales = "free_y") +
    geom_text(data = ann_m, aes(x = x_pos, y = y_pos, label = lab),
              hjust = 0, vjust = 0, size = 3.8, lineheight = 1.0) +
    labs(x = size_label, y = metric_label) +
    theme_nature()
}


p_nri <- make_panel_plot("NRI") + ggtitle("NRI vs Tree size")
p_nti <- make_panel_plot("NTI") + ggtitle("NTI vs Tree size")

# combine vertically
fig_panels <- (p_nri / p_nti) + plot_annotation(tag_levels = "a")

# -------------------------
# EXPORTS (only the three figures + minimal tables)
# -------------------------
use_cairo <- isTRUE(capabilities("cairo"))
dev_fun <- if (use_cairo) cairo_pdf else grDevices::pdf

# Figures
ggsave(file.path(out_dir, "FIG1_Forest_StdSlopes.pdf"), p_forest, width = 5.4, height = 5.8, device = dev_fun)
ggsave(file.path(out_dir, "FIG1_Forest_StdSlopes.tiff"), p_forest, width = 5.4, height = 5.8, dpi = 600, compression = "lzw")
# SVG
ggsave(file.path(out_dir, "FIG1_Forest_StdSlopes.svg"), p_forest, width = 5.4, height = 5.8, device = svglite::svglite)

ggsave(file.path(out_dir, "FIG2_RankPreservation.pdf"), p_rank, width = 5.0, height = 3.6, device = dev_fun)
ggsave(file.path(out_dir, "FIG2_RankPreservation.tiff"), p_rank, width = 5.0, height = 3.6, dpi = 600, compression = "lzw")
# SVG
ggsave(file.path(out_dir, "FIG2_RankPreservation.svg"), p_rank, width = 5.0, height = 3.6, device = svglite::svglite)

# Height for panels depends on number of facets
n_regions <- all_long %>% distinct(Region) %>% nrow()
rows_est  <- ceiling(n_regions / 3)
fig3_height <- max(4.5, rows_est * 1.8 * 2)  # two rows of panels (NRI, NTI)

ggsave(file.path(out_dir, "FIG3_RegressionPanels.pdf"), fig_panels, width = 7.2, height = fig3_height, device = dev_fun)
ggsave(file.path(out_dir, "FIG3_RegressionPanels.tiff"), fig_panels, width = 7.2, height = fig3_height, dpi = 600, compression = "lzw")
# SVG
ggsave(file.path(out_dir, "FIG3_RegressionPanels.svg"), fig_panels, width = 10.2, height = fig3_height, device = svglite::svglite)

# Tables (minimal)
readr::write_csv(reg,       file.path(out_dir, "table_core_regression_standardized.csv"))
readr::write_csv(rank_pres, file.path(out_dir, "table_rank_preservation.csv"))

message("Saved streamlined core figures and tables to ", normalizePath(out_dir))
