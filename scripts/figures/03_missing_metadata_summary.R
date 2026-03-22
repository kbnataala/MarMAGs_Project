################################################################################
# Script: 03_missing_metadata_summary.R
#
# Purpose:
#   Documentation template for plotting field-wise metadata missingness.
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution.
#
# Expected input:
#   - results/tables/MarMAGs_Additional_File_01.csv
#
# Expected outputs:
#   - results/figures/03_metadata_missingness.png
#   - results/figures/03_metadata_missingness.svg
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

base_dir <- "."
input_file <- file.path(base_dir, "results", "tables", "MarMAGs_Additional_File_01.csv")
output_dir <- file.path(base_dir, "results", "figures")
out_prefix <- file.path(output_dir, "03_metadata_missingness")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

# ------------------------------------------------------------------------------
# Helper: define missing values
# ------------------------------------------------------------------------------

is_missing <- function(x) {
  y <- if (is.character(x)) stringr::str_trim(x) else x
  is.na(y) | (is.character(y) & y == "")
}

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

raw <- suppressMessages(readr::read_csv(input_file, show_col_types = FALSE))

if (ncol(raw) < 3) {
  stop("Input table must contain at least 3 columns.")
}

# Exclude the first two columns from missingness analysis
dat <- raw[, -(1:min(2, ncol(raw))), drop = FALSE]

# Create missingness summary
field_summary <- tibble(
  Field   = names(dat),
  Missing = colSums(as.data.frame(lapply(dat, is_missing)), na.rm = TRUE)
) %>%
  mutate(
    Total = nrow(dat),
    Pct   = Missing / Total,
    Label = sprintf("%d (%.1f%%)", Missing, 100 * Pct)
  ) %>%
  arrange(Pct) %>%
  mutate(Field = factor(Field, levels = Field))

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

p <- ggplot(field_summary, aes(x = Pct, y = Field)) +
  geom_col(fill = "#2C7FB8", width = 0.75) +
  geom_text(
    aes(label = Label),
    hjust = -0.05,
    size = 4,
    fontface = "bold"
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Missing values",
    y = "Metadata field"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(size = 14)
  ) +
  coord_cartesian(xlim = c(0, max(field_summary$Pct) + 0.12))

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

plot_width <- 14
plot_height <- max(8, 0.28 * nrow(field_summary) + 2)

ggsave(
  filename = paste0(out_prefix, ".png"),
  plot = p,
  width = plot_width,
  height = plot_height,
  dpi = 600
)

ggsave(
  filename = paste0(out_prefix, ".svg"),
  plot = p,
  width = plot_width,
  height = plot_height
)

message("Missingness figure saved to: ", output_dir)
