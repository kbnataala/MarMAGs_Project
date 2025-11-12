#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures")

library(ggplot2)
library(ggmap)
library(sp)
library(maps)
library(cowplot)
library(car)
library(ggpubr)
library(svglite)


#loading the final marmags dataframe
marmags_dataframe <- read.csv(file = "01_MarMAGs_Dataframe_03062024_.csv", header = TRUE, stringsAsFactors = FALSE)

fig2_violin_plots <- marmags_dataframe[, c("Genome", "Genome_ID", "OTU_cluster", "Genome_source", "CheckM_Completeness",
                               "CheckM_Contamination", "CheckM_Strain_heterogeneity", "CheckM_Quality_score",
                               "CheckM_Quality", "BBTools_n_contigs", "BBTools_ctg_N50", "BBTools_gc_avg" )]

######################################################
#1- CheckM_Completeness
fig2_violin_plots$CheckM_Completeness

#Levene's Test for Homogeneity of Variances

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
### Figure S22

# Arrange plots in the desired layout

# Combine all
combined_plot <- plot_grid(Quality_score, completeness, contamination, Strain_heterogeneity, label_size = 34, ncol = 4, labels = c("(a)", "(b)", "(c)", "(d)"))

ggsave("/Thesis/work_package_3/Data_analysis/01_figures/05_genome_properties/08_combined_figure.png", combined_plot, width = 28, height = 10, dpi = 600)
ggsave("/Thesis/work_package_3/Data_analysis/01_figures/05_genome_properties/08_combined_figure.svg", combined_plot, width = 28, height = 10, dpi = 600)




