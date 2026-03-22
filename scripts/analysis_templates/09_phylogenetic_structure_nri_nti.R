################################################################################
# Script: 09_phylogenetic_structure_nri_nti.R
#
# Purpose:
#   Documentation template for calculating phylogenetic structure metrics
#   Net Relatedness Index (NRI) and Nearest Taxon Index (NTI) from a
#   community matrix and phylogenetic tree.
#
# Repository role:
#   This script is provided as a documentation template and may require
#   adaptation for local execution. It is not presented as a fully reproducible
#   pipeline component.
#
# Expected inputs:
#   - results/tables/territory_matrix.csv
#   - results/tables/gtdbtk_unrooted_bacteria.tree_trimmed_tree_X_8029
#
# Expected output:
#   - results/tables/nri_nti_territory.csv
################################################################################

suppressPackageStartupMessages({
  library(picante)
  library(ape)
})

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

base_dir <- "."
input_dir <- file.path(base_dir, "results", "tables")
output_dir <- input_dir

community_file <- file.path(input_dir, "territory_matrix.csv")
tree_file <- file.path(input_dir, "gtdbtk_unrooted_bacteria.tree_trimmed_tree_X_8029")
output_file <- file.path(output_dir, "nri_nti_territory.csv")

# ------------------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------------------

if (!file.exists(community_file)) {
  stop("Missing community matrix file: ", community_file)
}

if (!file.exists(tree_file)) {
  stop("Missing phylogenetic tree file: ", tree_file)
}

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

community_data <- read.csv(community_file, row.names = 1, check.names = FALSE)
phylo_tree <- read.tree(tree_file)

community_matrix <- as.matrix(community_data)

if (nrow(community_matrix) == 0 || ncol(community_matrix) == 0) {
  stop("Community matrix is empty.")
}

# ------------------------------------------------------------------------------
# Match taxa between matrix and tree
# ------------------------------------------------------------------------------

shared_taxa <- intersect(colnames(community_matrix), phylo_tree$tip.label)

if (length(shared_taxa) == 0) {
  stop("No shared taxa found between community matrix columns and tree tip labels.")
}

community_matrix <- community_matrix[, shared_taxa, drop = FALSE]
phylo_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, shared_taxa))

# ------------------------------------------------------------------------------
# Calculate NRI and NTI
# ------------------------------------------------------------------------------

tree_dist <- cophenetic(phylo_tree)

nri_results <- ses.mpd(
  samp = community_matrix,
  dis = tree_dist,
  null.model = "taxa.labels",
  runs = 999
)

nti_results <- ses.mntd(
  samp = community_matrix,
  dis = tree_dist,
  null.model = "taxa.labels",
  runs = 999
)

# ------------------------------------------------------------------------------
# Combine outputs
# ------------------------------------------------------------------------------

colnames(nri_results) <- paste0("NRI_", colnames(nri_results))
colnames(nti_results) <- paste0("NTI_", colnames(nti_results))

nri_results$Community <- rownames(nri_results)
nti_results$Community <- rownames(nti_results)

combined_results <- merge(
  nri_results,
  nti_results,
  by = "Community",
  all = TRUE
)

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

write.csv(combined_results, output_file, row.names = FALSE)

message("NRI and NTI results written to: ", output_file)
