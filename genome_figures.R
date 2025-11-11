#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures")

library(ggplot2)
library(ggmap)
library(sp)
library(maps)
library(cowplot)

#loading the final marmags dataframe
marmags_dataframe <- read.csv(file = "01_MarMAGs_Dataframe_03062024_.csv", header = TRUE, stringsAsFactors = FALSE)

fig1a <- marmags_dataframe[, c("Genome_source", "sample_latitude", "sample_longitude")]


# Remove rows with missing latitude or longitude
fig1a <- fig1a[!is.na(fig1a$sample_latitude) & !is.na(fig1a$sample_longitude), ]

# Define the base world map
mapworld <- borders("world", colour = "gray85", fill = "gray80")

svg("/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/04_world_map.svg", width = 12, height = 8)

# Plot the data
ggplot(fig1a) +
  mapworld +
  ylim(-90, 90) + # Limits for latitude
  geom_point(aes(x = sample_longitude, 
                 y = sample_latitude, 
                 color = factor(Genome_source, levels = c("MarMAGs", "GEM_catalog", "OceanDNA"))), 
             size = 3, shape = 16, alpha = 0.6) +
  scale_color_manual(values = c("#E33539", "#1f77b4", "#2ca02c"),  # Ensure colors match reordered factors
                     breaks = c("MarMAGs", "GEM_catalog", "OceanDNA")) +  # Reorder legend
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(
    text = element_text(size = 8, family = "sans", face = "plain"),
    axis.title = element_text(size = 16, face = "bold"), # Axis titles
    axis.text = element_text(size = 16), # Axis tick labels (latitude, longitude)
    legend.title = element_text(size = 16, face = "bold"), # Legend title
    legend.text = element_text(size = 16), # Legend items
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    legend.position = "bottom"
  )


dev.off()

#####################################################################################

#loading the final marmags dataframe
marmags_dataframe <- read.csv(file = "01_MarMAGs_Dataframe_03062024_.csv", header = TRUE, stringsAsFactors = FALSE)

fig2_violin_plots <- marmags_dataframe[, c("Genome", "Genome_ID", "OTU_cluster", "Genome_source", "CheckM_Completeness",
                               "CheckM_Contamination", "CheckM_Strain_heterogeneity", "CheckM_Quality_score",
                               "CheckM_Quality", "BBTools_n_contigs", "BBTools_ctg_N50", "BBTools_gc_avg" )]

######################################################
#1- CheckM_Completeness
fig2_violin_plots$CheckM_Completeness

#Levene's Test for Homogeneity of Variances
library(car)
levene_test <- leveneTest(CheckM_Completeness ~ Genome_source, data = fig2_violin_plots)
print(levene_test)

# Perform Welch's ANOVA (for unequal variances)
anova_welch <- oneway.test(CheckM_Completeness ~ Genome_source, data = fig2_violin_plots)
print(anova_welch)


# Pairwise Welch t-tests
pairwise_results <- pairwise.t.test(
  x = fig2_violin_plots$CheckM_Completeness,
  g = fig2_violin_plots$Genome_source,
  p.adjust.method = "bonferroni", # Adjust for multiple comparisons
  pool.sd = FALSE # Welch t-test
)
print(pairwise_results)


library(ggplot2)
library(ggpubr)

# Define the data for p-value annotations
stat_pvalue_data <- data.frame(
  group1 = c("OceanDNA", "GEM_catalog", "OceanDNA"),
  group2 = c("GEM_catalog", "MarMAGs", "MarMAGs"),
  p.adj = c("***", "***", "***"),
  y.position = c(108, 112, 116) # Adjust these values based on your data range
)

svg("/Thesis/work_package_3/Data_analysis/01_figures/05_genome_properties/01_completeness.svg", width = 10, height = 10)

# Beautified violin plot
completeness <- ggplot(fig2_violin_plots, aes(x = factor(Genome_source, levels = c("MarMAGs", "GEM_catalog", "OceanDNA")), y = CheckM_Completeness)) +
  geom_violin(aes(color = Genome_source), fill = "white", trim = FALSE, size = 1) +
  geom_boxplot(aes(color = Genome_source), fill = "white", width = 0.15, size = 0.8) +
  scale_color_manual(values = c("#1f77b4", "#d62728", "#2ca02c")) +  # Ensure colors match reordered factors
  stat_pvalue_manual(data = stat_pvalue_data, label = "p.adj", size = 10, label.size = 10, bracket.size = 1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 23),
    axis.title = element_text(size = 23),
    axis.line = element_line(size = 1.0)
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 23),
    axis.title = element_text(size = 23),
    axis.line = element_line(size = 1.0)
  ) +
  labs(x = "Genome Source", y = "CheckM Completeness (%)")

print(completeness)
dev.off()

######################################################
#2- CheckM_Contamination
fig2_violin_plots$CheckM_Contamination

#Levene's Test for Homogeneity of Variances
library(car)
levene_test <- leveneTest(CheckM_Contamination ~ Genome_source, data = fig2_violin_plots)
print(levene_test)

# Perform Welch's ANOVA (for unequal variances)
anova_welch <- oneway.test(CheckM_Contamination ~ Genome_source, data = fig2_violin_plots)
print(anova_welch)


# Pairwise Welch t-tests
pairwise_results <- pairwise.t.test(
  x = fig2_violin_plots$CheckM_Contamination,
  g = fig2_violin_plots$Genome_source,
  p.adjust.method = "bonferroni", # Adjust for multiple comparisons
  pool.sd = FALSE # Welch t-test
)
print(pairwise_results)


library(ggplot2)
library(ggpubr)

# Define the data for p-value annotations
stat_pvalue_data <- data.frame(
  group1 = c("OceanDNA", "GEM_catalog", "OceanDNA"),
  group2 = c("GEM_catalog", "MarMAGs", "MarMAGs"),
  p.adj = c("***", "***", "***"),
  y.position = c(11, 12, 13) # Adjust these values based on your data range
)

svg("/Thesis/work_package_3/Data_analysis/01_figures/05_genome_properties/02_contamination.svg", width = 10, height = 10)

# Beautified violin plot
contamination <- ggplot(fig2_violin_plots, aes(x = factor(Genome_source, levels = c("MarMAGs", "GEM_catalog", "OceanDNA")), y = CheckM_Contamination)) +
  geom_violin(aes(color = Genome_source), fill = "white", trim = FALSE, size = 1) +
  geom_boxplot(aes(color = Genome_source), fill = "white", width = 0.15, size = 0.8) +
  scale_color_manual(values = c("#1f77b4", "#d62728", "#2ca02c")) +  # Ensure colors match reordered factors
  stat_pvalue_manual(data = stat_pvalue_data, label = "p.adj", size = 10, label.size = 10, bracket.size = 1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 23),
    axis.title = element_text(size = 23),
    axis.line = element_line(size = 1.0)
  ) +
  labs(x = "Genome Source", y = "CheckM Contamination (%)")

print(contamination)
dev.off()
######################################################
#3- Strain_heterogeneity
fig2_violin_plots$CheckM_Strain_heterogeneity

#Levene's Test for Homogeneity of Variances
library(car)
levene_test <- leveneTest(CheckM_Strain_heterogeneity ~ Genome_source, data = fig2_violin_plots)
print(levene_test)

# Perform Welch's ANOVA (for unequal variances)
anova_welch <- oneway.test(CheckM_Strain_heterogeneity ~ Genome_source, data = fig2_violin_plots)
print(anova_welch)


# Pairwise Welch t-tests
pairwise_results <- pairwise.t.test(
  x = fig2_violin_plots$CheckM_Strain_heterogeneity,
  g = fig2_violin_plots$Genome_source,
  p.adjust.method = "bonferroni", # Adjust for multiple comparisons
  pool.sd = FALSE # Welch t-test
)
print(pairwise_results)


library(ggplot2)
library(ggpubr)

# Define the data for p-value annotations
stat_pvalue_data <- data.frame(
  group1 = c("OceanDNA", "GEM_catalog", "OceanDNA"),
  group2 = c("GEM_catalog", "MarMAGs", "MarMAGs"),
  p.adj = c("***", "***", "***"),
  y.position = c(120, 125, 130) # Adjust these values based on your data range
)

svg("/Thesis/work_package_3/Data_analysis/01_figures/05_genome_properties/03_strainheterogeneity.svg", width = 10, height = 10)

# Beautified violin plot
Strain_heterogeneity <- ggplot(fig2_violin_plots, aes(x = factor(Genome_source, levels = c("MarMAGs", "GEM_catalog", "OceanDNA")), y = CheckM_Strain_heterogeneity)) +
  geom_violin(aes(color = Genome_source), fill = "white", trim = FALSE, size = 1) +
  geom_boxplot(aes(color = Genome_source), fill = "white", width = 0.15, size = 0.8) +
  scale_color_manual(values = c("#1f77b4", "#d62728", "#2ca02c")) +  # Ensure colors match reordered factors
  stat_pvalue_manual(data = stat_pvalue_data, label = "p.adj", size = 10, label.size = 10, bracket.size = 1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 23),
    axis.title = element_text(size = 23),
    axis.line = element_line(size = 1.0)
  ) +
  labs(x = "Genome Source", y = "CheckM Strain Heterogeneity")

print(Strain_heterogeneity)
dev.off()
######################################################
#4- Quality score
fig2_violin_plots$CheckM_Quality_score

#Levene's Test for Homogeneity of Variances
library(car)
levene_test <- leveneTest(CheckM_Quality_score ~ Genome_source, data = fig2_violin_plots)
print(levene_test)

# Perform Welch's ANOVA (for unequal variances)
anova_welch <- oneway.test(CheckM_Quality_score ~ Genome_source, data = fig2_violin_plots)
print(anova_welch)


# Pairwise Welch t-tests
pairwise_results <- pairwise.t.test(
  x = fig2_violin_plots$CheckM_Quality_score,
  g = fig2_violin_plots$Genome_source,
  p.adjust.method = "bonferroni", # Adjust for multiple comparisons
  pool.sd = FALSE # Welch t-test
)
print(pairwise_results)


library(ggplot2)
library(ggpubr)

# Define the data for p-value annotations
stat_pvalue_data <- data.frame(
  group1 = c("OceanDNA", "GEM_catalog", "OceanDNA"),
  group2 = c("GEM_catalog", "MarMAGs", "MarMAGs"),
  p.adj = c("***", "***", "***"),
  y.position = c(105, 107, 109) # Adjust these values based on your data range
)

svg("/Thesis/work_package_3/Data_analysis/01_figures/05_genome_properties/04_qualityscore.svg", width = 10, height = 10)

# Beautified violin plot
Quality_score <- ggplot(fig2_violin_plots, aes(x = factor(Genome_source, levels = c("MarMAGs", "GEM_catalog", "OceanDNA")), y = CheckM_Quality_score)) +
  geom_violin(aes(color = Genome_source), fill = "white", trim = FALSE, size = 1) +
  geom_boxplot(aes(color = Genome_source), fill = "white", width = 0.15, size = 0.8) +
  scale_color_manual(values = c("#1f77b4", "#d62728", "#2ca02c")) +  # Ensure colors match reordered factors
  stat_pvalue_manual(data = stat_pvalue_data, label = "p.adj", size = 10, label.size = 10, bracket.size = 1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 23),
    axis.title = element_text(size = 23),
    axis.line = element_line(size = 1.0)
  ) +
  labs(x = "Genome Source", y = "CheckM Quality Score")

print(Quality_score)
dev.off()
######################################################
#5- Number of Contigs
fig2_violin_plots$BBTools_n_contigs

#Levene's Test for Homogeneity of Variances
library(car)
levene_test <- leveneTest(BBTools_n_contigs ~ Genome_source, data = fig2_violin_plots)
print(levene_test)

# Perform Welch's ANOVA (for unequal variances)
anova_welch <- oneway.test(BBTools_n_contigs ~ Genome_source, data = fig2_violin_plots)
print(anova_welch)


# Pairwise Welch t-tests
pairwise_results <- pairwise.t.test(
  x = fig2_violin_plots$BBTools_n_contigs,
  g = fig2_violin_plots$Genome_source,
  p.adjust.method = "bonferroni", # Adjust for multiple comparisons
  pool.sd = FALSE # Welch t-test
)
print(pairwise_results)


library(ggplot2)
library(ggpubr)

# Define the data for p-value annotations
stat_pvalue_data <- data.frame(
  group1 = c("OceanDNA", "GEM_catalog", "OceanDNA"),
  group2 = c("GEM_catalog", "MarMAGs", "MarMAGs"),
  p.adj = c("***", "***", "***"),
  y.position = c(4200, 4400, 4600) # Adjust these values based on your data range
)

svg("/Thesis/work_package_3/Data_analysis/01_figures/05_genome_properties/05_numContigs.svg", width = 10, height = 10)

# Beautified violin plot
n_contigs <- ggplot(fig2_violin_plots, aes(x = factor(Genome_source, levels = c("MarMAGs", "GEM_catalog", "OceanDNA")), y = BBTools_n_contigs)) +
  geom_violin(aes(color = Genome_source), fill = "white", trim = FALSE, size = 1) +
  geom_boxplot(aes(color = Genome_source), fill = "white", width = 0.15, size = 0.8) +
  scale_color_manual(values = c("#1f77b4", "#d62728", "#2ca02c")) +  # Ensure colors match reordered factors
  stat_pvalue_manual(data = stat_pvalue_data, label = "p.adj", size = 10, label.size = 10, bracket.size = 1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.line = element_line(size = 1.0)
  ) +
  labs(x = "Genome Source", y = "Number of Contigs")

print(n_contigs)
dev.off()
######################################################
#6- Contigs N50
fig2_violin_plots$BBTools_ctg_N50

#Levene's Test for Homogeneity of Variances
library(car)
levene_test <- leveneTest(BBTools_ctg_N50 ~ Genome_source, data = fig2_violin_plots)
print(levene_test)

# Perform Welch's ANOVA (for unequal variances)
anova_welch <- oneway.test(BBTools_ctg_N50 ~ Genome_source, data = fig2_violin_plots)
print(anova_welch)


# Pairwise Welch t-tests
pairwise_results <- pairwise.t.test(
  x = fig2_violin_plots$BBTools_ctg_N50,
  g = fig2_violin_plots$Genome_source,
  p.adjust.method = "bonferroni", # Adjust for multiple comparisons
  pool.sd = FALSE # Welch t-test
)
print(pairwise_results)


library(ggplot2)
library(ggpubr)

# Define the data for p-value annotations
stat_pvalue_data <- data.frame(
  group1 = c("OceanDNA", "GEM_catalog", "OceanDNA"),
  group2 = c("GEM_catalog", "MarMAGs", "MarMAGs"),
  p.adj = c("***", "***", "***"),
  y.position = c(1320, 1400, 1450) # Adjust these values based on your data range
)

svg("/Thesis/work_package_3/Data_analysis/01_figures/05_genome_properties/06_Contigs_N50.svg", width = 10, height = 10)

# Beautified violin plot
contigs_N50 <- ggplot(fig2_violin_plots, aes(x = factor(Genome_source, levels = c("MarMAGs", "GEM_catalog", "OceanDNA")), y = BBTools_ctg_N50)) +
  geom_violin(aes(color = Genome_source), fill = "white", trim = FALSE, size = 1) +
  geom_boxplot(aes(color = Genome_source), fill = "white", width = 0.15, size = 0.8) +
  scale_color_manual(values = c("#1f77b4", "#d62728", "#2ca02c")) +  # Ensure colors match reordered factors
  stat_pvalue_manual(data = stat_pvalue_data, label = "p.adj", size = 10, label.size = 10, bracket.size = 1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.line = element_line(size = 1.0)
  ) +
  labs(x = "Genome Source", y = "Contigs N50")

print(contigs_N50)
dev.off()
######################################################
#7- Contigs N50
fig2_violin_plots$BBTools_gc_avg

#Levene's Test for Homogeneity of Variances
library(car)
levene_test <- leveneTest(BBTools_gc_avg ~ Genome_source, data = fig2_violin_plots)
print(levene_test)

# Perform Welch's ANOVA (for unequal variances)
anova_welch <- oneway.test(BBTools_gc_avg ~ Genome_source, data = fig2_violin_plots)
print(anova_welch)


# Pairwise Welch t-tests
pairwise_results <- pairwise.t.test(
  x = fig2_violin_plots$BBTools_gc_avg,
  g = fig2_violin_plots$Genome_source,
  p.adjust.method = "bonferroni", # Adjust for multiple comparisons
  pool.sd = FALSE # Welch t-test
)
print(pairwise_results)


library(ggplot2)
library(ggpubr)

# Define the data for p-value annotations
stat_pvalue_data <- data.frame(
  group1 = c("OceanDNA", "OceanDNA"),
  group2 = c("GEM_catalog", "MarMAGs"),
  p.adj = c("***", "***"),
  y.position = c(0.85, 0.9) # Adjust these values based on your data range
)

svg("/Thesis/work_package_3/Data_analysis/01_figures/05_genome_properties/07_GC_average.svg", width = 10, height = 10)

# Beautified violin plot
gc_average <- ggplot(fig2_violin_plots, aes(x = factor(Genome_source, levels = c("MarMAGs", "GEM_catalog", "OceanDNA")), y = BBTools_gc_avg)) +
  geom_violin(aes(color = Genome_source), fill = "white", trim = FALSE, size = 1) +
  geom_boxplot(aes(color = Genome_source), fill = "white", width = 0.15, size = 0.8) +
  scale_color_manual(values = c("#1f77b4", "#d62728", "#2ca02c")) +  # Ensure colors match reordered factors
  stat_pvalue_manual(data = stat_pvalue_data, label = "p.adj", size = 10, label.size = 10, bracket.size = 1) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.line = element_line(size = 1.0)
  ) +
  labs(x = "Genome Source", y = "GC Average")

print(gc_average)
dev.off()
######################################################
library(svglite)


# Arrange plots in the desired layout

# Combine all
combined_plot <- plot_grid(Quality_score, completeness, contamination, Strain_heterogeneity, label_size = 34, ncol = 4, labels = c("(a)", "(b)", "(c)", "(d)"))

ggsave("/Thesis/work_package_3/Data_analysis/01_figures/05_genome_properties/08_combined_figure.png", combined_plot, width = 28, height = 10, dpi = 600)
ggsave("/Thesis/work_package_3/Data_analysis/01_figures/05_genome_properties/08_combined_figure.svg", combined_plot, width = 28, height = 10, dpi = 600)





































# ######################################################
# 
# # Load required libraries
# library(maps)
# library(ggplot2)
# library(dplyr)
# library(raster)
# library(tidyr)
# 
# # Step 1: Load your dataset
# ##loading the final marmags dataframe
# marmags_dataframe <- read.csv(file = "01_MarMAGs_Dataframe_03062024_.csv", header = TRUE, stringsAsFactors = FALSE)
# fig2 <- marmags_dataframe[, c("Genome_source", "sample_latitude", "sample_longitude")]
# 
# # Step 2: Preprocess the data
# # Rename the columns to match your data structure (if necessary)
# colnames(fig2) <- c("Genome_source", "sample_latitude", "sample_longitude")
# 
# # Select only Latitude and Longitude columns and filter out rows with NA
# mag_coordinates <- fig2 %>%
#   filter(!is.na(sample_latitude) & !is.na(sample_longitude)) %>%
#   select(sample_latitude, sample_longitude)
# 
# # Count the number of MAGs per unique coordinate
# mag_abundance <- mag_coordinates %>%
#   group_by(sample_latitude, sample_longitude) %>%
#   summarise(Count = n(), .groups = "drop")
# 
# # Step 1: Define the grid of coordinates
# lon <- seq(-180, 180, 0.2)  # Longitude range with 0.2° resolution
# lat <- seq(-90, 90, 0.2)    # Latitude range with 0.2° resolution
# 
# # Create a grid of coordinates
# coordinate <- expand.grid(longitude = lon, latitude = lat)
# 
# # Step 2: Map MAG abundance onto the grid
# # Assume mag_abundance contains columns: sample_latitude, sample_longitude, and Count
# # Use raster::extract to map the MAG abundance data onto the grid
# coordinate_raster <- coordinate %>%
#   rowwise() %>%
#   mutate(
#     MAG_count = sum(
#       mag_abundance$Count[
#         abs(mag_abundance$sample_latitude - latitude) <= 0.1 &
#           abs(mag_abundance$sample_longitude - longitude) <= 0.1
#       ], na.rm = TRUE
#     )
#   ) %>%
#   ungroup()
# 
# # Remove rows with NA MAG_count
# coordinate_raster <- coordinate_raster %>% filter(!is.na(MAG_count))
# 
# # Step 3: Plot the heatmap
# # Load world map data
# world <- ne_coastline(scale = "medium", returnclass = "sp")
# 
# # Generate the heatmap
# ggplot() +
#   geom_tile(data = coordinate_raster, aes(x = longitude, y = latitude, fill = MAG_count)) +
#   geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = NA, color = "gray50") +
#   scale_fill_gradientn(
#     colours = rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", 
#                     "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")),
#     values = scales::rescale(c(0, 0.5, 0.65, 0.75, 0.85, 1)),
#     name = "MAG Abundance"
#   ) +
#   coord_equal() +
#   theme_minimal() +
#   theme(axis.title = element_blank()) +
#   labs(title = "Global MAG Abundance Heatmap")
# 
# # Save the plot as a PDF
# ggsave("MAG_Abundance_Heatmap.pdf", width = 10, height = 6)



# # Step 3: Generate the heat map
# # Load world map data
# world_map <- map_data("world")
# 
# # Plot the heat map
# ggplot() +
#   geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
#                fill = "gray90", color = "gray50") +
#   geom_tile(data = mag_abundance_binned, aes(x = sample_longitude, y = sample_latitude, fill = Count)) +
#   scale_fill_gradient(low = "yellow", high = "red", name = "MAG Count") +
#   coord_fixed(ratio = 1.3) +
#   theme_minimal() +
#   labs(title = "Global Distribution of MAGs",
#        x = "Longitude",
#        y = "Latitude") +
#   theme(axis.title = element_text(size = 12),
#         axis.text = element_text(size = 10),
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10))
# 
# 
# mag_abundance_binned <- mag_abundance %>%
#   mutate(
#     bin_latitude = round(sample_latitude, 1),  # Round latitude to 1 decimal
#     bin_longitude = round(sample_longitude, 1) # Round longitude to 1 decimal
#   ) %>%
#   group_by(bin_latitude, bin_longitude) %>%
#   summarise(Count = sum(Count), .groups = "drop")
# 
# ggplot() +
#   geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
#                fill = "gray90", color = "gray50") +
#   geom_tile(data = mag_abundance_binned, aes(x = bin_longitude, y = bin_latitude, fill = Count)) +
#   scale_fill_gradient(low = "yellow", high = "red", name = "MAG Count") +
#   coord_fixed(ratio = 1.3) +
#   theme_minimal() +
#   labs(title = "Global Distribution of MAGs (Binned)",
#        x = "Longitude",
#        y = "Latitude") +
#   theme(axis.title = element_text(size = 12),
#         axis.text = element_text(size = 10),
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10))
# 
# 
# 
# 
# summary(mag_abundance)
# 
# ggplot(data = mag_abundance, aes(x = sample_longitude, y = sample_latitude)) +
#   geom_point() +
#   labs(title = "Coordinate Points", x = "Longitude", y = "Latitude") +
#   theme_minimal()
# 

# #violin plots
# 
# dodge <- position_dodge(width = 0.4)
# 
# ggplot(data = fig2_violin_plots, aes(x = Genome_source, y = CheckM_Quality_score, fill = CheckM_Quality)) +
#   geom_violin(position = dodge)+
#   geom_boxplot(width=.1, outlier.colour=NA, position = dodge) 
# 
# #############
# 
# 
# #################
# 
# library(ggplot2)
# theme_set(
#   theme_classic() +
#     theme(legend.position = "top")
# )
# 
# ggplot(fig2_violin_plots, aes(x = Genome_source, y = BBTools_gc_avg)) +
#   geom_violin(aes(fill = Genome_source), trim = FALSE) +
#   geom_boxplot(aes(fill = Genome_source), width = 0.15) +
#   scale_fill_manual(values = c("#E33539","#C18C00","#780062")) +  # Use scale_fill_manual
#   theme(legend.position = "none")
# 
# head(fig2_violin_plots)
# summary(fig2_violin_plots)
# ########################
# 
# library(dplyr)
# 
# # Shapiro-Wilk test for normality for each group in Genome_source
# normality_test <- fig2_violin_plots %>%
#   group_by(Genome_source) %>%
#   summarise(p_value = shapiro.test(BBTools_gc_avg)$p.value)
# 
# print(normality_test)
# 
# 
# 
# ####################
# 
# library(ggpubr)
# 
# # Create the violin plot with statistical significance annotations
# ggplot(fig2_violin_plots, aes(x = Genome_source, y = BBTools_gc_avg)) +
#   geom_violin(aes(fill = Genome_source), trim = FALSE) +
#   geom_boxplot(aes(fill = Genome_source), width = 0.15) +
#   scale_fill_manual(values = c("#E33539", "#C18C00", "#780062")) +
#   stat_compare_means(aes(group = Genome_source), 
#                      method = "anova", # Use "kruskal.test" if assumptions are violated
#                      label = "p.signif", 
#                      label.y = 0.95) +
#   theme(legend.position = "none")
# 
# 
# 
# kruskal_test <- kruskal.test(BBTools_gc_avg ~ Genome_source, data = fig2_violin_plots)
# print(kruskal_test)
# 
# ######################################################
