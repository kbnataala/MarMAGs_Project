################################################################################
# Script: 08_cf_bubble_plot.R
#
# Purpose:
#   Documentation template for visualizing carbon fixation gene presence across
#   top marine phyla using a bubble plot.
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution.

#   - results/supplementary_tables/MarMAGs_additional_File_02.docx(Fig. S10)
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

base_dir <- "."
input_dir <- file.path(base_dir, "results", "tables")
output_dir <- file.path(base_dir, "results", "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

taxonomy_file <- file.path(input_dir, "MarMAGs_phylum.csv")

pathway_files <- list(
  CBB   = file.path(input_dir, "CBB_pathway.csv"),
  rTCA  = file.path(input_dir, "rTCA_pathway.csv"),
  `3-HPB` = file.path(input_dir, "3HPB_bicycle.csv"),
  DCHB  = file.path(input_dir, "DCHB_pathway.csv"),
  HPHB  = file.path(input_dir, "HPHB_pathway.csv"),
  bWLJ  = file.path(input_dir, "bWJL_pathway.csv"),
  aWLJ  = file.path(input_dir, "aWJL_pathway.csv")
)

output_prefix <- file.path(output_dir, "08_cf_bubble_plot")

# ------------------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------------------

if (!file.exists(taxonomy_file)) {
  stop("Missing taxonomy file: ", taxonomy_file)
}

missing_files <- pathway_files[!file.exists(pathway_files)]
if (length(missing_files) > 0) {
  stop("Missing pathway files: ", paste(missing_files, collapse = ", "))
}

# ------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------

taxonomy <- read_csv(taxonomy_file, show_col_types = FALSE) %>%
  rename(Genome_ID = 1, Phylum = 2) %>%
  drop_na()

read_pa <- function(path, pathway_name) {
  read_csv(path, show_col_types = FALSE) %>%
    rename(Genome_ID = 1) %>%
    pivot_longer(-Genome_ID, names_to = "Gene", values_to = "Present") %>%
    filter(!is.na(Present), Present > 0) %>%
    mutate(Pathway = pathway_name) %>%
    select(Genome_ID, Gene, Pathway)
}

long_pa <- purrr::imap_dfr(pathway_files, read_pa)

# ------------------------------------------------------------------------------
# Merge + filter
# ------------------------------------------------------------------------------

joined <- long_pa %>%
  left_join(taxonomy, by = "Genome_ID") %>%
  drop_na(Phylum)

# Top 30 phyla
top_phyla <- taxonomy %>%
  count(Phylum, name = "n_genomes") %>%
  arrange(desc(n_genomes)) %>%
  slice_head(n = 30) %>%
  pull(Phylum)

joined <- joined %>% filter(Phylum %in% top_phyla)

# ------------------------------------------------------------------------------
# Aggregation
# ------------------------------------------------------------------------------

agg <- joined %>%
  distinct(Genome_ID, Phylum, Gene, Pathway) %>%
  count(Phylum, Gene, Pathway, name = "n_genomes")

phylum_totals <- taxonomy %>%
  filter(Phylum %in% top_phyla) %>%
  count(Phylum, name = "tot_genomes")

agg <- agg %>%
  left_join(phylum_totals, by = "Phylum") %>%
  mutate(percent = 100 * n_genomes / tot_genomes)

# ------------------------------------------------------------------------------
# Axis construction
# ------------------------------------------------------------------------------

axis_df <- agg %>%
  distinct(Pathway, Gene) %>%
  arrange(Pathway, Gene) %>%
  mutate(GeneAxis = paste(Pathway, Gene, sep = " | "))

agg <- agg %>%
  mutate(
    Phylum = factor(Phylum, levels = top_phyla),
    GeneAxis = factor(paste(Pathway, Gene, sep = " | "),
                      levels = rev(axis_df$GeneAxis))
  )

axis_labels <- setNames(axis_df$Gene, axis_df$GeneAxis)

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

p <- ggplot(agg, aes(x = Phylum, y = GeneAxis)) +
  geom_point(
    aes(size = percent, fill = Pathway),
    shape = 21, color = "grey30", alpha = 0.9
  ) +
  scale_y_discrete(labels = axis_labels) +
  scale_size(range = c(2, 14), name = "% genomes") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Phylum (top 30)",
    y = "Gene",
    fill = "Pathway"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

ggsave(paste0(output_prefix, ".png"), p, width = 16, height = 18, dpi = 300)
ggsave(paste0(output_prefix, ".svg"), p, width = 16, height = 18)

p
