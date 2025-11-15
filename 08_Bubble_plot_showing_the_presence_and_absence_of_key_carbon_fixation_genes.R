#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/09_CF_bubble_plot")

# CF bubble plot across top 30 phyla (robust CSV)
# - top counts closer to bubbles
# - reduced bottom padding
# - saves PNG, PDF, SVG
# ------------------------------------------------------------
library(tidyverse)
library(scales)    # pretty_breaks / comma
library(svglite)   # SVG output
library(Cairo)     # Cairo PDF device

# ===== USER SETTINGS =========================================================
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/09_CF_bubble_plot")

tax_path <- "MarMAGs_phylum.csv"

pa_paths <- c(
  CBB    = "CBB_pathway.csv",
  rTCA   = "rTCA_pathway.csv",
  `3-HPB`= "3HPB_bicycle.csv",
  DCHB   = "DCHB_pathway.csv",
  HPHB   = "HPHB_pathway.csv",
  bWLJ   = "bWJL_pathway.csv",
  aWLJ   = "aWJL_pathway.csv"
)

out_base <- "cf_bubble_plot"

# colors
path_pal <- c(
  CBB  = "#1b9e77", rTCA = "#d95f02", `3-HPB` = "#7570b3",
  DCHB = "#e7298a", HPHB = "#66a61e", bWLJ = "#e6ab02", aWLJ = "#a6761d"
)
# ============================================================================

# ---------- helper: auto-find a file if the provided path is missing ----------
auto_find <- function(path, patterns) {
  if (file.exists(path)) return(path)
  for (pat in patterns) {
    hits <- list.files(pattern = pat, ignore.case = TRUE, recursive = FALSE)
    if (length(hits) > 0) {
      message("Resolved '", path, "' to '", hits[1], "' via pattern \"", pat, "\"")
      return(hits[1])
    }
  }
  stop("Missing file and could not auto-find: ", path)
}

pathway_patterns <- list(
  CBB   = c("^CBB.*\\.csv$"),
  rTCA  = c("^rTCA.*\\.csv$"),
  `3-HPB`= c("^3\\s*[-–—]?\\s*HPB.*\\.csv$", "^03_?3HPB.*\\.csv$", "^3HPB.*\\.csv$"),
  DCHB  = c("^DCHB.*\\.csv$"),
  HPHB  = c("^HPHB.*\\.csv$"),
  bWLJ  = c("^bW[LJ][JL].*\\.csv$", "^bWLJ.*\\.csv$", "^bWJL.*\\.csv$"),
  aWLJ  = c("^aW[LJ][JL].*\\.csv$", "^aWLJ.*\\.csv$", "^aWJL.*\\.csv$")
)

# ---------- readers & cleaners ----------
read_tax <- function(p) {
  tax <- readr::read_csv(p, show_col_types = FALSE)
  nms <- names(tax)
  gcol <- nms[which(tolower(nms) %in% c("genome_id","genome","id","genome name","genome_name"))][1]
  pcol <- nms[which(tolower(nms) == "phylum")][1]
  if (is.na(gcol) || is.na(pcol))
    stop("Taxonomy CSV must contain columns for Genome_ID and Phylum.")
  tax |>
    rename(Genome_ID = all_of(gcol), Phylum = all_of(pcol)) |>
    select(Genome_ID, Phylum) |>
    drop_na()
}

read_pa <- function(path, pathway_label) {
  tb <- readr::read_csv(path, show_col_types = FALSE)
  if (!"Genome_ID" %in% names(tb)) names(tb)[1] <- "Genome_ID"
  tb |>
    pivot_longer(-Genome_ID, names_to = "Gene", values_to = "Present") |>
    mutate(Present = suppressWarnings(as.numeric(Present))) |>
    filter(!is.na(Present), Present > 0) |>
    mutate(Pathway = pathway_label) |>
    select(Genome_ID, Gene, Pathway)
}

canonicalize_pathway <- function(x) {
  x |>
    stringr::str_replace_all("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-") |>
    stringr::str_replace_all("\\s*-\\s*", "-") |>
    stringr::str_replace_all("^3\\s*HPB$", "3-HPB") |>
    stringr::str_trim()
}

# ---------- resolve any renamed/moved files ----------
if (!file.exists(tax_path)) {
  tax_hits <- list.files(pattern = "phylum.*\\.csv$", ignore.case = TRUE)
  if (length(tax_hits) == 0) stop("Cannot find taxonomy CSV (expected '", tax_path, "').")
  message("Resolved taxonomy to: ", tax_hits[1]); tax_path <- tax_hits[1]
}
for (nm in names(pa_paths)) pa_paths[nm] <- auto_find(pa_paths[nm], pathway_patterns[[nm]])

# ---------- load & combine ----------
taxonomy <- read_tax(tax_path)
long_pa  <- purrr::imap_dfr(pa_paths, ~ read_pa(.x, .y)) |>
  mutate(Pathway = canonicalize_pathway(Pathway))

exp_levels <- c("CBB","rTCA","3-HPB","DCHB","HPHB","bWLJ","aWLJ")
long_pa <- long_pa |> mutate(Pathway = factor(Pathway, levels = exp_levels))

joined <- long_pa |>
  left_join(taxonomy, by = "Genome_ID") |>
  drop_na(Phylum)

# Top 30 phyla by # genomes (from taxonomy)
top_phyla <- taxonomy |>
  count(Phylum, name = "n_genomes") |>
  arrange(desc(n_genomes)) |>
  slice_head(n = 30) |>
  pull(Phylum)

joined <- joined |> filter(Phylum %in% top_phyla)

# Aggregate presences
agg <- joined |>
  distinct(Genome_ID, Phylum, Gene, Pathway) |>
  count(Phylum, Gene, Pathway, name = "n_genomes")

# Totals per phylum and % sizing
phylum_totals <- taxonomy |>
  filter(Phylum %in% top_phyla) |>
  count(Phylum, name = "tot_genomes")

agg <- agg |>
  left_join(phylum_totals, by = "Phylum") |>
  mutate(percent = 100 * n_genomes / tot_genomes)

# ---------- UNIQUE y-axis (avoid duplicated gene names across pathways) ----------
axis_df <- agg |>
  distinct(Pathway, Gene) |>
  arrange(Pathway, Gene) |>
  mutate(GeneAxis = paste(Pathway, Gene, sep = " | "))

gene_axis_levels <- axis_df$GeneAxis
axis_labels      <- setNames(axis_df$Gene, axis_df$GeneAxis)

agg <- agg |>
  mutate(
    Phylum   = factor(Phylum, levels = top_phyla),
    GeneAxis = factor(paste(Pathway, Gene, sep = " | "),
                      levels = rev(gene_axis_levels)),
    Pathway  = factor(Pathway, levels = exp_levels)
  )

# ---------- size legend (4 bubbles) ----------
size_breaks <- pretty_breaks(n = 4)(agg$percent) |> unique() |> sort()
if (length(size_breaks) == 0) size_breaks <- c(5, 15, 30, 60)

# ---------- phylum counts (top labels) ----------
count_labels <- phylum_totals |>
  filter(Phylum %in% top_phyla) |>
  mutate(
    Phylum = factor(Phylum, levels = top_phyla),
    label  = paste0("n=", tot_genomes)
  )

# ---- knobs for bigger text (no change to plot dimensions) ----
base_text <- 16      # was 12
label_size <- 4.5    # was 3

# ---------- plot ----------
p <- ggplot(agg, aes(x = Phylum, y = GeneAxis)) +
  geom_point(
    aes(size = percent, fill = Pathway),
    shape = 21, color = "grey20", alpha = 0.9
  ) +
  # counts at top of panel (vertical)
  geom_text(
    data = count_labels,
    aes(x = Phylum, y = Inf, label = label),
    inherit.aes = FALSE,
    angle = 270, hjust = 0.4, vjust = 0.05,
    size = label_size
  ) +
  scale_y_discrete(
    labels = axis_labels,
    expand = expansion(mult = c(0.01, 0.08))
  ) +
  scale_fill_manual(values = path_pal, name = "Carbon fixation pathways", drop = FALSE) +
  scale_size(
    range = c(2, 14),
    breaks = size_breaks,
    name   = "Bubble size\n(% genomes)"
  ) +
  labs(x = "Phylum (top 30 by # genomes)", y = "Gene") +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 5, shape = 21)),
    size = guide_legend(order = 2, override.aes = list(fill = "white", colour = "grey30"))
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = base_text) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = base_text),  # bigger tick labels
    axis.text.y = element_text(size = base_text),
    axis.title  = element_text(face = "bold", size = base_text + 2),
    legend.title = element_text(face = "bold", size = base_text + 1),
    legend.text  = element_text(size = base_text),
    legend.key.size = unit(7, "mm"),
    panel.grid.major = element_line(linetype = "dotted", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    plot.margin = margin(t = 26, r = 22, b = 8, l = 12)  # tiny bump for larger text
  )

# Figure size (unchanged)
w <- max(16, 0.4 * length(levels(agg$Phylum)) + 2)
h <- max(18, 0.35 * length(levels(agg$GeneAxis)) + 2)

# ---------- SAVE ----------
ggsave(paste0(out_base, ".png"), p, width = w, height = h, dpi = 300)
ggsave(paste0(out_base, ".pdf"), p, width = w, height = h, device = cairo_pdf)
ggsave(paste0(out_base, ".svg"), p, width = w, height = h)

p
