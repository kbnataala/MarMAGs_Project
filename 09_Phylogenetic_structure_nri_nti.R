# Install and load necessary packages
#install.packages("picante")
library(picante)

# Load the community data (your abundance matrix) and the phylogenetic tree
# Assuming 'community_data.csv' is your abundance matrix file
community_data <- read.csv("territory_matrix.csv", row.names = 1)
phylo_tree <- read.tree("gtdbtk_unrooted_bacteria.tree_trimmed_tree_X_8029")

# Convert the data to a matrix
community_matrix <- as.matrix(community_data)

# Calculate NRI and NTI for each community (each phylum is treated as a 'community')
# NRI is calculated using 'ses.mpd' (Mean Pairwise Distance)
nri_results <- ses.mpd(community_matrix, cophenetic(phylo_tree), null.model="taxa.labels", runs=999)

# NTI is calculated using 'ses.mntd' (Mean Nearest Taxon Distance)
nti_results <- ses.mntd(community_matrix, cophenetic(phylo_tree), null.model="taxa.labels", runs=999)

# View the results
print(nri_results)
print(nti_results)


# Combine the NTI and NRI results into a single data frame
# Assuming `nri_results` and `nti_results` are already loaded as data frames

# Rename columns to avoid duplicate names when combining
colnames(nri_results) <- paste0("NRI_", colnames(nri_results))
colnames(nti_results) <- paste0("NTI_", colnames(nti_results))

# Add a Phylum column to each data frame for merging
nri_results$Phylum <- rownames(nri_results)
nti_results$Phylum <- rownames(nti_results)

# Merge the NRI and NTI results by Phylum
combined_results <- merge(nri_results, nti_results, by = "Phylum", all = TRUE)


# Write the combined data to a CSV file
write.csv(combined_results, "nri_nti_territory.csv", row.names = FALSE)

# Print confirmation
cat("The NRI and NTI results have been saved to 'nri_nti_class.csv'")


