#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence")

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Read the input data
data_table <- read.csv(file = "merged_data.csv", header = TRUE, stringsAsFactors = FALSE)

# Define a function to calculate the prevalence
calculate_prevalence <- function(data, taxonomy_level, pathway_column = "DCHB") {
  data %>%
    group_by(across(all_of(taxonomy_level)), .add = TRUE) %>%
    summarize(
      Pathway = first(!!sym(pathway_column)), # Use the pathway column name
      Number_of_taxa = n(),
      Pathway_presence = sum(!!sym(pathway_column)),
      Prevalence = Pathway_presence / Number_of_taxa
    ) %>%
    ungroup()
}

# Specify taxonomy levels to calculate the metrics
#taxonomy_levels <- c("GTDB_Tk_Domain", "GTDB_Tk_Phylum", "GTDB_Tk_Class", "GTDB_Tk_Order", "GTDB_Tk_Family", "GTDB_Tk_Genus")

# Specify taxonomy levels to calculate the metrics
#taxonomy_levels <- c("OTU_cluster", "GTDB_Tk_Phylum")
taxonomy_levels <- c("GTDB_Tk_Class")

# Generate the results for each taxonomy level
results <- lapply(taxonomy_levels, function(level) {
  calculate_prevalence(data_table, level)
})

# Combine all results into one data frame
final_output <- bind_rows(results, .id = "Taxa_level")

# Write the output to a file
write.csv(final_output, "DCHB_class_prevalence.csv", row.names = FALSE)

#write.csv(final_output, "CBB_per_OTU_prevalence_output.csv", row.names = FALSE)

# Load dplyr just to be safe
library(dplyr)

# Summarize presence/absence counts of DCHB per phylum
contingency_table <- data_table %>%
  group_by(GTDB_Tk_Phylum, DCHB) %>%
  summarize(Count = n()) %>%
  tidyr::pivot_wider(names_from = DCHB, values_from = Count, values_fill = 0)

# Optional: rename columns for clarity
colnames(contingency_table) <- c("Phylum", "Absent", "Present")

################################################

# Define a function to calculate prevalence per OTU cluster
calculate_prevalence_otu <- function(data, pathway_column = "CBB") {
  data %>%
    group_by(OTU_cluster) %>%
    summarize(
      Number_of_MAGs = n(), # Total number of MAGs in the OTU cluster
      Pathway_presence = sum(!!sym(pathway_column)), # Count of MAGs with the pathway
      Prevalence = Pathway_presence / Number_of_MAGs, # Prevalence
      Pathway = ifelse(Pathway_presence > 0, 1, 0) # Assign 1 if pathway is present in any MAG
    ) %>%
    ungroup()
}

# Calculate prevalence per OTU cluster
otu_prevalence <- calculate_prevalence_otu(data_table)

# Remove the Pathway column if redundant
otu_prevalence <- otu_prevalence %>% select(-Pathway)

# Write the output to a file
write.csv(otu_prevalence, "otu_prevalence_output.csv", row.names = FALSE)

######################################
######################################




######################################
######################################

# 1 # CF Pathways and Phylums

# Read the input data
prev_data <- read.csv(file = "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_phylum.csv", header = TRUE, stringsAsFactors = FALSE)

prev_data <- prev_data[prev_data$Prevalence > 0, ]


svg("/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/01_Prevalence_pathways_phylum.svg",width = 16, height = 28)
# Plot the heatmap with the new color range
prev_plot1 <- ggplot(prev_data, aes(x = GTDB_Tk_Phylum, y = Pathway)) +
  geom_tile(aes(fill = Prevalence), color = "black") +
  # scale_fill_gradientn(colors = rev(c("#000000", "#330002", 
  #                                     "#4D0000", "#6A0010", "#871224", "#A4313A", 
  #                                     "#C24D51", "#E06769", "#FF8282"))) +
  scale_fill_continuous(low = "#93160E", high = "lightgreen") +
  # scale_fill_gradientn(colors = rev(c("#FFEBD5", "#FFD2AB",
  #                                     "#FFB481", "#FF9661", "#FF652D", "#DB4520",
  #                                     "#B72B16", "#93160E", "#7A0809"))) +
    #7A0809
  #93160E
  #B72B16
  #DB4520
  #FF652D
  #FF9661
  #FFB481
  #FFD2AB
  #FFEBD5
  
  #7A0809
  #93160E
  #B72B16
  #DB4520
  #FF652D
  #FF9661
  #FFB481
  #FFD2AB
  #FFEBD5
  
  #  scale_fill_gradientn(colors = rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", 
#                                      "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", 
#                                      "#66C2A5", "#3288BD", "#5E4FA2"))) +
  theme_minimal() +
  coord_flip() +  
  labs(title = "",
       x = "Phylum", y = "Pathway", fill = "Prevalence") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), # Increase title size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Standardize x-axis text
    axis.text.y = element_text(size = 18),  # Standardize y-axis text
    axis.title.x = element_text(size = 20),  # X-axis title size
    axis.title.y = element_text(size = 20),  # Y-axis title size
    legend.key.width = unit(2, "cm"),  # Standard legend box size
    legend.key.height = unit(2, "cm"),
    legend.title = element_text(size = 22),  # Legend title size
    legend.text = element_text(size = 18),  # Legend text size
    plot.tag = element_text(size = 30, face = "bold")  # Increase label size (A-E)
  ) 

print(prev_plot1)

dev.off()

# Setup
library(dplyr)
pathways <- unique(prev_data$Pathway)
GTDB_Tk_Phylum <- unique(prev_data$GTDB_Tk_Phylum)
pairs <- combn(GTDB_Tk_Phylum, 2, simplify = FALSE)

# Create empty results table
results_table <- data.frame(
  Pathway = character(),
  Comparison = character(),
  Prevalence_1 = numeric(),
  Prevalence_2 = numeric(),
  p_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Loop through pathways and salinity pairs
for (pw in pathways) {
  sub_df <- prev_data[prev_data$Pathway == pw, ]
  
  for (pair in pairs) {
    sal1 <- pair[1]
    sal2 <- pair[2]
    
    row1 <- sub_df[sub_df$GTDB_Tk_Phylum == sal1, ]
    row2 <- sub_df[sub_df$GTDB_Tk_Phylum == sal2, ]
    
    if (nrow(row1) == 0 || nrow(row2) == 0) next
    
    x <- c(row1$Pathway_presence, row2$Pathway_presence)
    n <- c(row1$Number_of_taxa, row2$Number_of_taxa)
    
    test <- prop.test(x, n, correct = FALSE)
    
    results_table <- rbind(results_table, data.frame(
      Pathway = pw,
      Comparison = paste(sal1, "vs", sal2),
      Prevalence_1 = round(x[1]/n[1], 6),
      Prevalence_2 = round(x[2]/n[2], 6),
      p_value = round(test$p.value, 5),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    ))
  }
}

# View the table
print(results_table)

# Optional: export to CSV
write.csv(results_table, "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/GTDB_Tk_Phylum_pairwise_prevalence_comparison.csv", row.names = FALSE)

######################################
# 2 # CF Pathways and Class

# Read the input data
prev_data <- read.csv(file = "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_class.csv", header = TRUE, stringsAsFactors = FALSE)

prev_data <- prev_data[prev_data$Prevalence > 0, ]


svg("/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/01_Prevalence_pathways_class.svg",width = 12, height = 32)
# Plot the heatmap with the new color range
ggplot(prev_data, aes(x = GTDB_Tk_Class, y = Pathway)) +
  geom_tile(aes(fill = Prevalence), color = "black") +
  scale_fill_gradientn(colors = rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", 
                                      "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", 
                                      "#66C2A5", "#3288BD", "#5E4FA2"))) +
  theme_minimal() +
  coord_flip() +  
  labs(title = "Prevalence of CF Pathway Across Class",
       x = "Phylum",
       y = "Pathway",
       fill = "Prevalence") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), # Increase title size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Standardize x-axis text
    axis.text.y = element_text(size = 18),  # Standardize y-axis text
    axis.title.x = element_text(size = 20),  # X-axis title size
    axis.title.y = element_text(size = 20),  # Y-axis title size
    legend.key.width = unit(2, "cm"),  # Standard legend box size
    legend.key.height = unit(2, "cm"),
    legend.title = element_text(size = 22),  # Legend title size
    legend.text = element_text(size = 18),  # Legend text size
    plot.tag = element_text(size = 30, face = "bold")  # Increase label size (A-E)
  )

dev.off()

# Setup
library(dplyr)
pathways <- unique(prev_data$Pathway)
class <- unique(prev_data$GTDB_Tk_Class)
pairs <- combn(class, 2, simplify = FALSE)

# Create empty results table
results_table <- data.frame(
  Pathway = character(),
  Comparison = character(),
  Prevalence_1 = numeric(),
  Prevalence_2 = numeric(),
  p_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Loop through pathways and salinity pairs
for (pw in pathways) {
  sub_df <- prev_data[prev_data$Pathway == pw, ]
  
  for (pair in pairs) {
    sal1 <- pair[1]
    sal2 <- pair[2]
    
    row1 <- sub_df[sub_df$GTDB_Tk_Class == sal1, ]
    row2 <- sub_df[sub_df$GTDB_Tk_Class == sal2, ]
    
    if (nrow(row1) == 0 || nrow(row2) == 0) next
    
    x <- c(row1$Pathway_presence, row2$Pathway_presence)
    n <- c(row1$Number_of_taxa, row2$Number_of_taxa)
    
    test <- prop.test(x, n, correct = FALSE)
    
    results_table <- rbind(results_table, data.frame(
      Pathway = pw,
      Comparison = paste(sal1, "vs", sal2),
      Prevalence_1 = round(x[1]/n[1], 6),
      Prevalence_2 = round(x[2]/n[2], 6),
      p_value = round(test$p.value, 5),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    ))
  }
}

# View the table
print(results_table)

# Optional: export to CSV
write.csv(results_table, "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/class_pairwise_prevalence_comparison.csv", row.names = FALSE)

######################################
# 3 # CF Pathways and Depth

# Read the input data
prev_data <- read.csv(file = "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_depth.csv", header = TRUE, stringsAsFactors = FALSE)

prev_data <- prev_data[prev_data$Prevalence > 0, ]


svg("/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/03_Prevalence_pathways_depth.svg",width = 12, height = 8)
# Plot the heatmap with the new color range
prev_plot2 <- ggplot(prev_data, aes(x = Depth_Category, y = Pathway)) +
  geom_tile(aes(fill = Prevalence), color = "black") +
  scale_fill_gradientn(colors = rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", 
                                      "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", 
                                      "#66C2A5", "#3288BD", "#5E4FA2"))) +
  theme_minimal() +
  coord_flip() +  
  labs(title = "",
       x = "Depth",
       y = "Pathway",
       fill = "Prevalence") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), # Increase title size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Standardize x-axis text
    axis.text.y = element_text(size = 18),  # Standardize y-axis text
    axis.title.x = element_text(size = 20),  # X-axis title size
    axis.title.y = element_text(size = 20),  # Y-axis title size
    legend.key.width = unit(2, "cm"),  # Standard legend box size
    legend.key.height = unit(2, "cm"),
    legend.title = element_text(size = 22),  # Legend title size
    legend.text = element_text(size = 18),  # Legend text size
    plot.tag = element_text(size = 30, face = "bold")  # Increase label size (A-E)
  )

print(prev_plot2)

dev.off()

# Setup
library(dplyr)
pathways <- unique(prev_data$Pathway)
depth <- unique(prev_data$Depth_Category)
pairs <- combn(depth, 2, simplify = FALSE)

# Create empty results table
results_table <- data.frame(
  Pathway = character(),
  Comparison = character(),
  Prevalence_1 = numeric(),
  Prevalence_2 = numeric(),
  p_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Loop through pathways and salinity pairs
for (pw in pathways) {
  sub_df <- prev_data[prev_data$Pathway == pw, ]
  
  for (pair in pairs) {
    sal1 <- pair[1]
    sal2 <- pair[2]
    
    row1 <- sub_df[sub_df$Depth_Category == sal1, ]
    row2 <- sub_df[sub_df$Depth_Category == sal2, ]
    
    if (nrow(row1) == 0 || nrow(row2) == 0) next
    
    x <- c(row1$Pathway_presence, row2$Pathway_presence)
    n <- c(row1$Number_of_taxa, row2$Number_of_taxa)
    
    test <- prop.test(x, n, correct = FALSE)
    
    results_table <- rbind(results_table, data.frame(
      Pathway = pw,
      Comparison = paste(sal1, "vs", sal2),
      Prevalence_1 = round(x[1]/n[1], 6),
      Prevalence_2 = round(x[2]/n[2], 6),
      p_value = round(test$p.value, 5),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    ))
  }
}

# View the table
print(results_table)

# Optional: export to CSV
write.csv(results_table, "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/depth_pairwise_prevalence_comparison.csv", row.names = FALSE)


######################################
# 4 # CF Pathways and IHO Sea

# Read the input data
prev_data <- read.csv(file = "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_IHOsea.csv", header = TRUE, stringsAsFactors = FALSE)

prev_data <- prev_data[prev_data$Prevalence > 0, ]

svg("/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/04_Prevalence_pathways_IHOsea.svg",width = 16, height = 28)
# Plot the heatmap with the new color range
ggplot(prev_data, aes(x = IHO_Sea, y = Pathway)) +
  geom_tile(aes(fill = Prevalence), color = "black") +
  scale_fill_gradientn(colors = rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", 
                                      "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", 
                                      "#66C2A5", "#3288BD", "#5E4FA2"))) +
  theme_minimal() +
  coord_flip() +  
  labs(title = "",
       x = "IHO Sea",
       y = "Pathway",
       fill = "Prevalence") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), # Increase title size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Standardize x-axis text
    axis.text.y = element_text(size = 18),  # Standardize y-axis text
    axis.title.x = element_text(size = 20),  # X-axis title size
    axis.title.y = element_text(size = 20),  # Y-axis title size
    legend.key.width = unit(2, "cm"),  # Standard legend box size
    legend.key.height = unit(2, "cm"),
    legend.title = element_text(size = 22),  # Legend title size
    legend.text = element_text(size = 18),  # Legend text size
    plot.tag = element_text(size = 30, face = "bold")  # Increase label size (A-E)
  )

dev.off()

# Setup
library(dplyr)
pathways <- unique(prev_data$Pathway)
IHO_Sea <- unique(prev_data$IHO_Sea)
pairs <- combn(IHO_Sea, 2, simplify = FALSE)

# Create empty results table
results_table <- data.frame(
  Pathway = character(),
  Comparison = character(),
  Prevalence_1 = numeric(),
  Prevalence_2 = numeric(),
  p_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Loop through pathways and salinity pairs
for (pw in pathways) {
  sub_df <- prev_data[prev_data$Pathway == pw, ]
  
  for (pair in pairs) {
    sal1 <- pair[1]
    sal2 <- pair[2]
    
    row1 <- sub_df[sub_df$IHO_Sea == sal1, ]
    row2 <- sub_df[sub_df$IHO_Sea == sal2, ]
    
    if (nrow(row1) == 0 || nrow(row2) == 0) next
    
    x <- c(row1$Pathway_presence, row2$Pathway_presence)
    n <- c(row1$Number_of_taxa, row2$Number_of_taxa)
    
    test <- prop.test(x, n, correct = FALSE)
    
    results_table <- rbind(results_table, data.frame(
      Pathway = pw,
      Comparison = paste(sal1, "vs", sal2),
      Prevalence_1 = round(x[1]/n[1], 6),
      Prevalence_2 = round(x[2]/n[2], 6),
      p_value = round(test$p.value, 5),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    ))
  }
}

# View the table
print(results_table)

# Optional: export to CSV
write.csv(results_table, "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/IHO_Sea_pairwise_prevalence_comparison.csv", row.names = FALSE)

######################################
# 5 # CF Pathways and Marine Regions

# Read the input data
prev_data <- read.csv(file = "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_MarRegion.csv", header = TRUE, stringsAsFactors = FALSE)

prev_data <- prev_data[prev_data$Prevalence > 0, ]


svg("/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/05_Prevalence_pathways_MarRegion.svg",width = 20, height = 36)
# Plot the heatmap with the new color range
ggplot(prev_data, aes(x = MarRegion, y = Pathway)) +
  geom_tile(aes(fill = Prevalence), color = "black") +
  scale_fill_gradientn(colors = rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", 
                                      "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", 
                                      "#66C2A5", "#3288BD", "#5E4FA2"))) +
  theme_minimal() +
  coord_flip() +  
  labs(title = "",
       x = "Marine Regions",
       y = "Pathway",
       fill = "Prevalence") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), # Increase title size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Standardize x-axis text
    axis.text.y = element_text(size = 18),  # Standardize y-axis text
    axis.title.x = element_text(size = 20),  # X-axis title size
    axis.title.y = element_text(size = 20),  # Y-axis title size
    legend.key.width = unit(2, "cm"),  # Standard legend box size
    legend.key.height = unit(2, "cm"),
    legend.title = element_text(size = 22),  # Legend title size
    legend.text = element_text(size = 18),  # Legend text size
    plot.tag = element_text(size = 30, face = "bold")  # Increase label size (A-E)
  )

dev.off()

# Setup
library(dplyr)
pathways <- unique(prev_data$Pathway)
MarRegion <- unique(prev_data$MarRegion)
pairs <- combn(MarRegion, 2, simplify = FALSE)

# Create empty results table
results_table <- data.frame(
  Pathway = character(),
  Comparison = character(),
  Prevalence_1 = numeric(),
  Prevalence_2 = numeric(),
  p_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Loop through pathways and salinity pairs
for (pw in pathways) {
  sub_df <- prev_data[prev_data$Pathway == pw, ]
  
  for (pair in pairs) {
    sal1 <- pair[1]
    sal2 <- pair[2]
    
    row1 <- sub_df[sub_df$MarRegion == sal1, ]
    row2 <- sub_df[sub_df$MarRegion == sal2, ]
    
    if (nrow(row1) == 0 || nrow(row2) == 0) next
    
    x <- c(row1$Pathway_presence, row2$Pathway_presence)
    n <- c(row1$Number_of_taxa, row2$Number_of_taxa)
    
    test <- prop.test(x, n, correct = FALSE)
    
    results_table <- rbind(results_table, data.frame(
      Pathway = pw,
      Comparison = paste(sal1, "vs", sal2),
      Prevalence_1 = round(x[1]/n[1], 6),
      Prevalence_2 = round(x[2]/n[2], 6),
      p_value = round(test$p.value, 5),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    ))
  }
}

# View the table
print(results_table)

# Optional: export to CSV
write.csv(results_table, "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/MarRegion_pairwise_prevalence_comparison.csv", row.names = FALSE)

######################################
# 6 # CF Pathways and Ocean and Sea

# Read the input data
prev_data <- read.csv(file = "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_ocean_sea.csv", header = TRUE, stringsAsFactors = FALSE)

prev_data <- prev_data[prev_data$Prevalence > 0, ]


svg("/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/06_Prevalence_pathways_ocean_sea.svg",width = 14, height = 14)
# Plot the heatmap with the new color range
prev_plot3 <- ggplot(prev_data, aes(x = ocean_sea_name, y = Pathway)) +
  geom_tile(aes(fill = Prevalence), color = "black") +
  scale_fill_gradientn(colors = rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", 
                                      "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", 
                                      "#66C2A5", "#3288BD", "#5E4FA2"))) +
  theme_minimal() +
  coord_flip() +  
  labs(title = "",
       x = "Ocean and Sea",
       y = "Pathway",
       fill = "Prevalence") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), # Increase title size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Standardize x-axis text
    axis.text.y = element_text(size = 18),  # Standardize y-axis text
    axis.title.x = element_text(size = 20),  # X-axis title size
    axis.title.y = element_text(size = 20),  # Y-axis title size
    legend.key.width = unit(2, "cm"),  # Standard legend box size
    legend.key.height = unit(2, "cm"),
    legend.title = element_text(size = 22),  # Legend title size
    legend.text = element_text(size = 18),  # Legend text size
    plot.tag = element_text(size = 30, face = "bold")  # Increase label size (A-E)
  )

print(prev_plot3)

dev.off()

# Setup
library(dplyr)
pathways <- unique(prev_data$Pathway)
ocean_sea_name <- unique(prev_data$ocean_sea_name)
pairs <- combn(ocean_sea_name, 2, simplify = FALSE)

# Create empty results table
results_table <- data.frame(
  Pathway = character(),
  Comparison = character(),
  Prevalence_1 = numeric(),
  Prevalence_2 = numeric(),
  p_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Loop through pathways and salinity pairs
for (pw in pathways) {
  sub_df <- prev_data[prev_data$Pathway == pw, ]
  
  for (pair in pairs) {
    sal1 <- pair[1]
    sal2 <- pair[2]
    
    row1 <- sub_df[sub_df$ocean_sea_name == sal1, ]
    row2 <- sub_df[sub_df$ocean_sea_name == sal2, ]
    
    if (nrow(row1) == 0 || nrow(row2) == 0) next
    
    x <- c(row1$Pathway_presence, row2$Pathway_presence)
    n <- c(row1$Number_of_taxa, row2$Number_of_taxa)
    
    test <- prop.test(x, n, correct = FALSE)
    
    results_table <- rbind(results_table, data.frame(
      Pathway = pw,
      Comparison = paste(sal1, "vs", sal2),
      Prevalence_1 = round(x[1]/n[1], 6),
      Prevalence_2 = round(x[2]/n[2], 6),
      p_value = round(test$p.value, 5),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    ))
  }
}

# View the table
print(results_table)

# Optional: export to CSV
write.csv(results_table, "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/ocean_sea_name_pairwise_prevalence_comparison.csv", row.names = FALSE)

######################################
# CF Pathways and Salinity

# Read the input data
prev_data <- read.csv(file = "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_salinity.csv", header = TRUE, stringsAsFactors = FALSE)

prev_data <- prev_data[prev_data$Prevalence > 0, ]


svg("/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/07_Prevalence_pathways_salinity.svg",width = 12, height = 8)
# Plot the heatmap with the new color range
prev_plot4 <- ggplot(prev_data, aes(x = Salinity_Category, y = Pathway)) +
  geom_tile(aes(fill = Prevalence), color = "black") +
  scale_fill_gradientn(colors = rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", 
                                      "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", 
                                      "#66C2A5", "#3288BD", "#5E4FA2"))) +
  theme_minimal() +
  coord_flip() +  
  labs(title = "",
       x = "Salinity",
       y = "Pathway",
       fill = "Prevalence") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), # Increase title size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Standardize x-axis text
    axis.text.y = element_text(size = 18),  # Standardize y-axis text
    axis.title.x = element_text(size = 20),  # X-axis title size
    axis.title.y = element_text(size = 20),  # Y-axis title size
    legend.key.width = unit(2, "cm"),  # Standard legend box size
    legend.key.height = unit(2, "cm"),
    legend.title = element_text(size = 22),  # Legend title size
    legend.text = element_text(size = 18),  # Legend text size
    plot.tag = element_text(size = 30, face = "bold")  # Increase label size (A-E)
  )

print(prev_plot4)

dev.off()

# Setup
library(dplyr)
pathways <- unique(prev_data$Pathway)
salinities <- unique(prev_data$Salinity_Category)
pairs <- combn(salinities, 2, simplify = FALSE)

# Create empty results table
results_table <- data.frame(
  Pathway = character(),
  Comparison = character(),
  Prevalence_1 = numeric(),
  Prevalence_2 = numeric(),
  p_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Loop through pathways and salinity pairs
for (pw in pathways) {
  sub_df <- prev_data[prev_data$Pathway == pw, ]
  
  for (pair in pairs) {
    sal1 <- pair[1]
    sal2 <- pair[2]
    
    row1 <- sub_df[sub_df$Salinity_Category == sal1, ]
    row2 <- sub_df[sub_df$Salinity_Category == sal2, ]
    
    if (nrow(row1) == 0 || nrow(row2) == 0) next
    
    x <- c(row1$Pathway_presence, row2$Pathway_presence)
    n <- c(row1$Number_of_taxa, row2$Number_of_taxa)
    
    test <- prop.test(x, n, correct = FALSE)
    
    results_table <- rbind(results_table, data.frame(
      Pathway = pw,
      Comparison = paste(sal1, "vs", sal2),
      Prevalence_1 = round(x[1]/n[1], 6),
      Prevalence_2 = round(x[2]/n[2], 6),
      p_value = round(test$p.value, 5),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    ))
  }
}

# View the table
print(results_table)

# Optional: export to CSV
write.csv(results_table, "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/salinity_pairwise_prevalence_comparison.csv", row.names = FALSE)

######################################
# CF Pathways and Temperature

# Read the input data
prev_data <- read.csv(file = "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_temp.csv", header = TRUE, stringsAsFactors = FALSE)

prev_data <- prev_data[prev_data$Prevalence > 0, ]


svg("/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/08_Prevalence_pathways_temperature.svg",width = 12, height = 8)
# Plot the heatmap with the new color range
prev_plot5 <- ggplot(prev_data, aes(x = temp_Category, y = Pathway)) +
  geom_tile(aes(fill = Prevalence), color = "black") +
  scale_fill_gradientn(colors = rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", 
                                      "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", 
                                      "#66C2A5", "#3288BD", "#5E4FA2"))) +
  theme_minimal() +
  coord_flip() +  
  labs(title = "",
       x = "Temperature",
       y = "Pathway",
       fill = "Prevalence") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), # Increase title size
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),  # Standardize x-axis text
    axis.text.y = element_text(size = 18),  # Standardize y-axis text
    axis.title.x = element_text(size = 20),  # X-axis title size
    axis.title.y = element_text(size = 20),  # Y-axis title size
    legend.key.width = unit(2, "cm"),  # Standard legend box size
    legend.key.height = unit(2, "cm"),
    legend.title = element_text(size = 22),  # Legend title size
    legend.text = element_text(size = 18),  # Legend text size
    plot.tag = element_text(size = 30, face = "bold")  # Increase label size (A-E)
  )

print(prev_plot5)

dev.off()

# Setup
library(dplyr)
pathways <- unique(prev_data$Pathway)
temp_Category <- unique(prev_data$temp_Category)
pairs <- combn(temp_Category, 2, simplify = FALSE)

# Create empty results table
results_table <- data.frame(
  Pathway = character(),
  Comparison = character(),
  Prevalence_1 = numeric(),
  Prevalence_2 = numeric(),
  p_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

# Loop through pathways and salinity pairs
for (pw in pathways) {
  sub_df <- prev_data[prev_data$Pathway == pw, ]
  
  for (pair in pairs) {
    sal1 <- pair[1]
    sal2 <- pair[2]
    
    row1 <- sub_df[sub_df$temp_Category == sal1, ]
    row2 <- sub_df[sub_df$temp_Category == sal2, ]
    
    if (nrow(row1) == 0 || nrow(row2) == 0) next
    
    x <- c(row1$Pathway_presence, row2$Pathway_presence)
    n <- c(row1$Number_of_taxa, row2$Number_of_taxa)
    
    test <- prop.test(x, n, correct = FALSE)
    
    results_table <- rbind(results_table, data.frame(
      Pathway = pw,
      Comparison = paste(sal1, "vs", sal2),
      Prevalence_1 = round(x[1]/n[1], 6),
      Prevalence_2 = round(x[2]/n[2], 6),
      p_value = round(test$p.value, 5),
      Significant = ifelse(test$p.value < 0.05, "Yes", "No")
    ))
  }
}

# View the table
print(results_table)

# Optional: export to CSV
write.csv(results_table, "/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/temp_pairwise_prevalence_comparison.csv", row.names = FALSE)

######################################

# Simulating the 4 smaller heatmaps (replace with real plots)
prev_plot1 <- prev_plot1 + labs(tag = "A")
prev_plot2 <- prev_plot2 + labs(tag = "B")
prev_plot3 <- prev_plot3 + labs(tag = "E")
prev_plot4 <- prev_plot4 + labs(tag = "C")
prev_plot5 <- prev_plot5 + labs(tag = "D")


# Arrange the heatmaps: A (big) on the left, B-E (small) on the right
heatmaps_layout <- (prev_plot1 | (prev_plot2 / prev_plot4 / prev_plot5 / prev_plot3)) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Save as SVG
svg("/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/merged_heatmaps.svg",
    width = 20, height = 16)
print(heatmaps_layout)
dev.off()

# Display the layout in RStudio
heatmaps_layout

################
#######################################################
#######################################################
#####
##### New code for Prevalence Plots
#######################################################
#######################################################

# ======= Libraries =======
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggsignif)

# ======= 1) Read & prep prevalence (≥1% AND >10 taxa) =======
df_prev <- read.csv(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_phylum.csv",
  stringsAsFactors = FALSE
)

cbb_df2 <- df_prev %>%
  filter(Pathway == "CBB",
         Prevalence >= 0.01,
         Number_of_taxa > 10) %>%
  arrange(Prevalence) %>%
  mutate(
    Phylum = factor(GTDB_Tk_Phylum, levels = GTDB_Tk_Phylum),
    frac_label = sprintf("%.2f", Prevalence)  # fraction labels on bars
    # If you prefer percent labels on bars (while keeping axis 0–1), use:
    # pct_label = percent(Prevalence, accuracy = 0.1)
  )

# ======= 2) Read & prep pairwise sigs (only among displayed phyla) =======
df_comp <- read.csv(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/GTDB_Tk_Phylum_pairwise_prevalence_comparison.csv",
  stringsAsFactors = FALSE
)

STAT <- df_comp %>%
  filter(Pathway == "CBB", Significant == "Yes") %>%
  separate(Comparison, into = c("Group1","Group2"), sep = " vs ") %>%
  filter(Group1 %in% levels(cbb_df2$Phylum),
         Group2 %in% levels(cbb_df2$Phylum)) %>%
  mutate(
    labels = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Stagger y-positions for significance brackets within [0, 1)
if (nrow(STAT) > 0) {
  y_cap   <- 0.99
  y_start <- min(y_cap - 0.02, max(cbb_df2$Prevalence, na.rm = TRUE) + 0.02)
  STAT$y_pos <- if (nrow(STAT) == 1) y_start else seq(y_start, y_cap, length.out = nrow(STAT))
}

# Reshape for geom_signif(manual=TRUE)
STAT2 <- if (nrow(STAT) > 0) {
  STAT %>%
    transmute(
      xmin = Group1,
      xmax = Group2,
      annotations = labels,
      y_position = y_pos
    )
} else {
  STAT # empty
}

# ======= 3) Plot (y-axis strictly 0–1) =======
p <- ggplot(cbb_df2, aes(x = Phylum, y = Prevalence)) +
  # Set bar fill color
  geom_col(width = 0.6, fill = "#4C72B0") +  # Oxford Blue
  # Bar labels (fractions). Switch to aes(label = pct_label) if you prefer %.
  geom_text(aes(label = frac_label), vjust = -0.5, size = 4, fontface = "bold") +
  # Only add significance if present
  { if (nrow(STAT2) > 0) geom_signif(
    inherit.aes = FALSE,
    data = STAT2,
    aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position),
    manual = TRUE,
    tip_length = 0.02,
    textsize = 6,
    color = "black"   # always black lines & stars
  ) else NULL } +
  scale_y_continuous(
    limits = c(0, 1),
    labels = label_number(accuracy = 0.01),  # show 0.00–1.00 on axis
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_cartesian(clip = "on") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x   = element_text(angle = 90, hjust = 1, size = 12, face = "bold"),
    axis.text.y   = element_text(size = 20),
    axis.title.y  = element_text(size = 20, face = "bold"),
    plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 20, hjust = 0.5),
    panel.border  = element_rect(color = "black", fill = NA),
    axis.line     = element_line(color = "black")
  ) +
  labs(
    x        = NULL,
    y        = "CBB Pathway Prevalence (fraction)",
    title    = "CBB Pathway Prevalence Among Major Phyla",
    subtitle = "Phyla with Prevalence ≥ 0.01 and Number_of_taxa > 10; *p≤0.05, **p≤0.01, ***p≤0.001"
  )

print(p)

# ======= 4) Save publication-ready files =======
# PNG (RGB)
ggsave("CBB_phylum_prevalence_fraction.png", plot = p, width = 8, height = 5, dpi = 300)
# SVG
ggsave("CBB_phylum_prevalence_fraction.svg", plot = p, width = 8, height = 5)
# PDF with Cairo (better color embedding for some printers)
# install.packages("Cairo") if needed
ggsave("CBB_phylum_prevalence_fraction.pdf", plot = p, width = 8, height = 5, device = cairo_pdf)

####################################################################################
####################################################################################
########### Function to plot the figure for multiple tables

# ==============================
# Reusable plotting function
# ==============================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggsignif)
  library(rlang)
})

plot_prevalence_by_taxon <- function(
    prev_csv,
    comp_csv,
    pathway            = "CBB",
    taxon_col          = "GTDB_Tk_Phylum",
    min_prev           = 0.01,
    min_n_taxa         = 10,
    label_mode         = c("fraction","percent"),
    bar_fill           = "#1D3557",
    bracket_color      = "black",
    # ---- New controls for y-axis scaling ----
    auto_ylim          = FALSE,     # FALSE = keep 0–1; TRUE = adapt per-plot
    y_limit            = c(0, 1),   # used when auto_ylim = FALSE
    y_top_manual       = NULL,      # override top (e.g., 0.15). If set, auto_ylim is ignored
    top_pad            = 0.02,      # space above tallest bar before first bracket
    bracket_step       = 0.03,      # vertical step between brackets
    # -----------------------------------------
    outfile_prefix     = NULL,
    width              = 8,
    height             = 8,
    dpi                = 900,
    save_formats       = c("png","svg","pdf")
) {
  label_mode <- match.arg(label_mode)
  
  df_prev <- read.csv(prev_csv, stringsAsFactors = FALSE)
  df_comp <- read.csv(comp_csv, stringsAsFactors = FALSE)
  
  taxon_sym <- rlang::sym(taxon_col)
  
  cbb_df2 <- df_prev %>%
    filter(Pathway == pathway,
           Prevalence >= min_prev,
           Number_of_taxa > min_n_taxa) %>%
    arrange(Prevalence) %>%
    mutate(
      !!taxon_sym := as.character(!!taxon_sym),
      GroupVar = !!taxon_sym,
      GroupVar = factor(GroupVar, levels = unique(GroupVar))
    )
  
  if (nrow(cbb_df2) == 0) {
    stop("No rows left after filtering. Check min_prev / min_n_taxa / pathway / taxon_col.")
  }
  
  cbb_df2 <- cbb_df2 %>%
    mutate(
      bar_label = if (label_mode == "percent")
        percent(Prevalence, accuracy = 0.1)
      else
        sprintf("%.2f", Prevalence)
    )
  
  # --- Sig table
  STAT <- df_comp %>%
    filter(Pathway == pathway, Significant == "Yes") %>%
    separate(Comparison, into = c("Group1","Group2"), sep = " vs ") %>%
    filter(Group1 %in% levels(cbb_df2$GroupVar),
           Group2 %in% levels(cbb_df2$GroupVar)) %>%
    mutate(
      labels = case_when(
        p_value <= 0.001 ~ "***",
        p_value <= 0.01  ~ "**",
        p_value <= 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  
  # --- Decide y-axis top (and bracket positions) ---
  max_prev <- max(cbb_df2$Prevalence, na.rm = TRUE)
  n_comp   <- nrow(STAT)
  
  # Manual y-top wins
  if (!is.null(y_top_manual)) {
    y_top <- y_top_manual
  } else if (isTRUE(auto_ylim)) {
    # adapt axis to data + brackets, rounded up to 0.01
    y_start <- max_prev + top_pad
    y_cap_needed <- if (n_comp <= 1) y_start else y_start + (n_comp - 1)*bracket_step
    y_top_raw <- max(y_cap_needed + top_pad, max_prev * 1.12)  # a bit of headroom
    y_top <- min(1, ceiling(y_top_raw * 100) / 100)            # round up to 2 d.p., cap at 1
    # ensure it's not silly-small
    y_top <- max(y_top, 0.1)
  } else {
    y_top <- y_limit[2]
  }
  
  # Now that we know y_top, compute bracket positions and cap them safely
  if (n_comp > 0) {
    cap <- y_top - 0.01
    y_start <- min(max_prev + top_pad, cap)
    y_end   <- min(cap, y_start + (n_comp - 1)*bracket_step)
    STAT$y_pos <- if (n_comp == 1) y_start else seq(y_start, y_end, length.out = n_comp)
  }
  
  STAT2 <- if (n_comp > 0) {
    STAT %>% transmute(
      xmin = Group1,
      xmax = Group2,
      annotations = labels,
      y_position = y_pos
    )
  } else STAT
  
  # --- Plot
  p <- ggplot(cbb_df2, aes(x = GroupVar, y = Prevalence)) +
    geom_col(width = 0.7, fill = bar_fill) +
    geom_text(aes(label = bar_label), vjust = -0.4, size = 12, fontface = "bold") +
    { if (nrow(STAT2) > 0) geom_signif(
      inherit.aes = FALSE,
      data = STAT2,
      aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position),
      manual = TRUE,
      tip_length = 0.005,
      textsize = 22,
      color = bracket_color
    ) else NULL } +
    scale_y_continuous(
      limits = c(0, y_top),
      labels = if (label_mode == "percent") percent_format(accuracy = 1) else label_number(accuracy = 0.01),
      expand = expansion(mult = c(0, 0.01))
    ) +
    coord_cartesian(clip = "on") +
    theme_classic(base_size = 22) +
    theme(
      axis.text.x   = element_text(angle = 90, hjust = 1, size = 28, face = "bold"),
      axis.text.y   = element_text(size = 28),
      axis.title.y  = element_text(size = 28, face = "bold"),
      plot.title    = element_text(size = 28, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 28, hjust = 0.5),
      panel.border  = element_rect(color = "black", fill = NA),
      axis.line     = element_line(color = "black")
    ) +
    labs(
      x        = NULL,
      y        = if (label_mode == "percent") "Prevalence (%)" else "Prevalence",
      #title    = paste(pathway, "Pathway Prevalence by", taxon_col),
      #subtitle = sprintf("Filters: Prevalence ≥ %s & Number_of_taxa > %s; y-top = %.2f", min_prev, min_n_taxa, y_top)
    )
  
  if (!is.null(outfile_prefix)) {
    if ("png" %in% save_formats) ggsave(paste0(outfile_prefix, ".png"), plot = p, width = width, height = height, dpi = dpi)
    if ("svg" %in% save_formats) ggsave(paste0(outfile_prefix, ".svg"), plot = p, width = width, height = height)
    if ("pdf" %in% save_formats) {
      dev <- tryCatch(grDevices::cairo_pdf, error = function(e) grDevices::pdf)
      ggsave(paste0(outfile_prefix, ".pdf"), plot = p, width = width, height = height, device = dev)
    }
  }
  
  return(p)
}

########################################

cbb_phylum <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_phylum.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/GTDB_Tk_Phylum_pairwise_prevalence_comparison.csv",
  pathway        = "CBB",
  taxon_col      = "GTDB_Tk_Phylum",
  min_prev       = 0.01,
  min_n_taxa     = 10,
  label_mode     = "fraction",   # or "percent"
  bar_fill       = "#4C72B0",    # change this freely for bar color
  bracket_color  = "black",
  #auto_ylim = FALSE,            # keep 0–1
  #y_limit   = c(0, 12),
  auto_ylim = TRUE,
  #y_top_manual = 0.12,
  outfile_prefix = "CBB_phylum_fraction"  # saves .png/.svg/.pdf
)


cbb_depth <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_depth.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/depth_pairwise_prevalence_comparison.csv",
  pathway        = "CBB",
  taxon_col      = "Depth_Category",
  min_prev       = 0.01,
  min_n_taxa     = 10,
  label_mode     = "fraction",   # or "percent"
  bar_fill       = "#4C72B0",    # change this freely for bar color
  bracket_color  = "black",
  auto_ylim = TRUE,
  outfile_prefix = "CBB_phylum_fraction"  # saves .png/.svg/.pdf
)

cbb_temp <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_temp.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/temp_pairwise_prevalence_comparison.csv",
  pathway        = "CBB",
  taxon_col      = "temp_Category",
  min_prev       = 0.01,
  min_n_taxa     = 10,
  label_mode     = "fraction",   # or "percent"
  bar_fill       = "#4C72B0",    # change this freely for bar color
  bracket_color  = "black",
  auto_ylim = TRUE,
  outfile_prefix = "CBB_phylum_fraction"  # saves .png/.svg/.pdf
)


cbb_sal <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_salinity.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/salinity_pairwise_prevalence_comparison.csv",
  pathway        = "CBB",
  taxon_col      = "Salinity_Category",
  min_prev       = 0.01,
  min_n_taxa     = 10,
  label_mode     = "fraction",   # or "percent"
  bar_fill       = "#4C72B0",    # change this freely for bar color
  bracket_color  = "black",
  auto_ylim = TRUE,
  outfile_prefix = "CBB_phylum_fraction"  # saves .png/.svg/.pdf
)


cbb_ocean <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_ocean_sea.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/ocean_sea_name_pairwise_prevalence_comparison.csv",
  pathway        = "CBB",
  taxon_col      = "ocean_sea_name",
  min_prev       = 0.01,
  min_n_taxa     = 10,
  label_mode     = "fraction",   # or "percent"
  bar_fill       = "#4C72B0",    # change this freely for bar color
  bracket_color  = "black",
  outfile_prefix = "CBB_phylum_fraction.png"  # saves .png/.svg/.pdf
)

ggsave("CBB_prevalence_ocean.svg", plot = cbb_ocean, width = 15, height = 17)

cbb_sea <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_IHOsea.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/IHO_Sea_pairwise_prevalence_comparison.csv",
  pathway        = "CBB",
  taxon_col      = "IHO_Sea",
  min_prev       = 0.03,
  min_n_taxa     = 10,
  label_mode     = "fraction",   # or "percent"
  bar_fill       = "#4C72B0",    # change this freely for bar color
  bracket_color  = "black",
  outfile_prefix = "CBB_phylum_fraction"  # saves .png/.svg/.pdf
)

ggsave("CBB_prevalence_sea.svg", plot = cbb_sea, width = 15, height = 18)

##################################################
##################################################

plot_prevalence_by_taxon <- function(
    prev_csv,
    comp_csv,
    pathway            = "CBB",
    taxon_col          = "GTDB_Tk_Phylum",
    min_prev           = 0.01,
    min_n_taxa         = 10,
    label_mode         = c("fraction","percent"),
    bar_fill           = "#1D3557",
    bracket_color      = "black",
    # ---- y-axis & bracket spacing ----
    auto_ylim          = FALSE,
    y_limit            = c(0, 1),
    y_top_manual       = NULL,
    top_pad            = 0.015,   # was 0.02  → tighter
    bracket_step       = 0.025,   # was 0.03  → tighter
    # ----------------------------------
    outfile_prefix     = NULL,
    width              = 8,
    height             = 8,
    dpi                = 900,
    save_formats       = c("png","svg","pdf")
) {
  label_mode <- match.arg(label_mode)
  
  df_prev <- read.csv(prev_csv, stringsAsFactors = FALSE)
  df_comp <- read.csv(comp_csv, stringsAsFactors = FALSE)
  
  taxon_sym <- rlang::sym(taxon_col)
  
  cbb_df2 <- df_prev %>%
    dplyr::filter(Pathway == pathway,
                  Prevalence >= min_prev,
                  Number_of_taxa > min_n_taxa) %>%
    dplyr::arrange(Prevalence) %>%
    dplyr::mutate(
      !!taxon_sym := as.character(!!taxon_sym),
      GroupVar = !!taxon_sym,
      GroupVar = factor(GroupVar, levels = unique(GroupVar))
    )
  
  if (nrow(cbb_df2) == 0) {
    stop("No rows left after filtering. Check min_prev / min_n_taxa / pathway / taxon_col.")
  }
  
  cbb_df2 <- cbb_df2 %>%
    dplyr::mutate(
      bar_label = if (label_mode == "percent")
        scales::percent(Prevalence, accuracy = 0.1)
      else
        sprintf("%.2f", Prevalence)
    )
  
  # --- Sig table
  STAT <- df_comp %>%
    dplyr::filter(Pathway == pathway, Significant == "Yes") %>%
    tidyr::separate(Comparison, into = c("Group1","Group2"), sep = " vs ") %>%
    dplyr::filter(Group1 %in% levels(cbb_df2$GroupVar),
                  Group2 %in% levels(cbb_df2$GroupVar)) %>%
    dplyr::mutate(
      labels = dplyr::case_when(
        p_value <= 0.001 ~ "***",
        p_value <= 0.01  ~ "**",
        p_value <= 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  
  # --- Decide y-axis top (and bracket positions) ---
  max_prev <- max(cbb_df2$Prevalence, na.rm = TRUE)
  n_comp   <- nrow(STAT)
  
  if (!is.null(y_top_manual)) {
    y_top <- y_top_manual
  } else if (isTRUE(auto_ylim)) {
    y_start <- max_prev + top_pad
    y_cap_needed <- if (n_comp <= 1) y_start else y_start + (n_comp - 1)*bracket_step
    y_top_raw <- max(y_cap_needed + top_pad, max_prev * 1.12)
    y_top <- min(1, ceiling(y_top_raw * 100) / 100)
    y_top <- max(y_top, 0.1)
  } else {
    y_top <- y_limit[2]
  }
  
  if (n_comp > 0) {
    cap <- y_top - 0.01
    y_start <- min(max_prev + top_pad, cap)
    y_end   <- min(cap, y_start + (n_comp - 1)*bracket_step)
    STAT$y_pos <- if (n_comp == 1) y_start else seq(y_start, y_end, length.out = n_comp)
  }
  
  STAT2 <- if (n_comp > 0) {
    dplyr::transmute(STAT,
                     xmin = Group1,
                     xmax = Group2,
                     annotations = labels,
                     y_position = y_pos
    )
  } else STAT
  
  # --- Build the signif layer with version-safe args ---------------------------
  have_margin_top <- requireNamespace("ggsignif", quietly = TRUE) &&
    utils::packageVersion("ggsignif") >= "0.6.4"
  
  signif_layer <- if (nrow(STAT2) > 0) {
    if (have_margin_top) {
      ggsignif::geom_signif(
        inherit.aes = FALSE,
        data = STAT2,
        aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position),
        manual      = TRUE,
        tip_length  = 0.002,  # shorter tips
        vjust       = 1.1,    # push stars DOWN, closer to line
        margin_top  = 0.005,  # tighten gap above line
        textsize    = 20,
        linewidth   = 0.7,
        color       = bracket_color
      )
    } else {
      # Fallback for older ggsignif without margin_top
      ggsignif::geom_signif(
        inherit.aes = FALSE,
        data = STAT2,
        aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position),
        manual      = TRUE,
        tip_length  = 0.002,
        vjust       = 1.15,   # a tad more push if margin_top is unavailable
        textsize    = 20,
        linewidth   = 0.7,
        color       = bracket_color
      )
    }
  } else NULL
  
  # --- Plot --------------------------------------------------------------------
  p <- ggplot(cbb_df2, aes(x = GroupVar, y = Prevalence)) +
    geom_col(width = 0.7, fill = bar_fill) +
    geom_text(aes(label = bar_label), vjust = -0.4, size = 12, fontface = "bold") +
    signif_layer +
    scale_y_continuous(
      limits = c(0, y_top),
      labels = if (label_mode == "percent")
        scales::percent_format(accuracy = 1)
      else
        scales::label_number(accuracy = 0.01),
      expand = expansion(mult = c(0, 0.01))
    ) +
    coord_cartesian(clip = "on") +
    theme_classic(base_size = 22) +
    theme(
      axis.text.x   = element_text(angle = 90, hjust = 1, size = 28, face = "bold"),
      axis.text.y   = element_text(size = 28),
      axis.title.y  = element_text(size = 28, face = "bold"),
      plot.title    = element_text(size = 28, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 28, hjust = 0.5),
      panel.border  = element_rect(color = "black", fill = NA),
      axis.line     = element_line(color = "black")
    ) +
    labs(
      x = NULL,
      y = if (label_mode == "percent") "Prevalence (%)" else "Prevalence"
    )
  
  if (!is.null(outfile_prefix)) {
    if ("png" %in% save_formats) ggsave(paste0(outfile_prefix, ".png"), plot = p, width = width, height = height, dpi = dpi)
    if ("svg" %in% save_formats) ggsave(paste0(outfile_prefix, ".svg"), plot = p, width = width, height = height)
    if ("pdf" %in% save_formats) {
      dev <- tryCatch(grDevices::cairo_pdf, error = function(e) grDevices::pdf)
      ggsave(paste0(outfile_prefix, ".pdf"), plot = p, width = width, height = height, device = dev)
    }
  }
  
  return(p)
}

##################################################
##################################################

cbb_class <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_class.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/class_pairwise_prevalence_comparison.csv",
  pathway        = "CBB",
  taxon_col      = "GTDB_Tk_Class",
  min_prev       = 0.025,
  min_n_taxa     = 10,
  label_mode     = "fraction",   # or "percent"
  bar_fill       = "#4C72B0",    # change this freely for bar color
  bracket_color  = "black",
  outfile_prefix = "CBB_prevalence_class"  # saves .png/.svg/.pdf
)

ggsave("CBB_prevalence_class.svg", plot = cbb_class, width = 25, height = 39)

#dev.off()

#install.packages("patchwork")  # if needed
library(patchwork)
library(ggplot2)

# ---- Helper to assemble and save a grid of ggplots ----
compose_grid <- function(
    plots,                 # list of ggplot objects
    ncol = 2,              # number of columns
    nrow = NULL,           # number of rows (optional)
    title = NULL,          # optional overall title
    tag_levels = "A",      # panel tags A, B, C... ("a", "1", "i" also work)
    tag_prefix = "(",      # e.g., "(" -> "(A)"
    tag_suffix = ")",      # e.g., ")" -> "(A)"
    legend = c("keep","collect","none"), # usually "none" for your single-color bars
    file = NULL,           # e.g., "combined.png" (if NULL, just returns the plot)
    width = 174,           # width in mm (174 mm ~ 2-column width)
    height = 120,          # height in mm
    dpi = 300
) {
  legend <- match.arg(legend)
  
  # Build layout
  p <- wrap_plots(plotlist = plots, ncol = ncol, nrow = nrow,
                  guides = if (legend == "collect") "collect" else "keep")
  
  # Add global title + panel tags
  p <- p + plot_annotation(
    title = title,
    tag_levels = tag_levels,
    tag_prefix = tag_prefix,
    tag_suffix = tag_suffix
  )
  
  # Optionally drop all legends
  if (legend == "none") {
    p <- p & theme(legend.position = "none")
  }
  
  # Save if requested
  if (!is.null(file)) {
    ggsave(file, p, width = width, height = height, dpi = dpi, units = "mm")
  }
  
  return(p)
}

# ---- Example usage ----
# Suppose you already have ggplot objects like:
# p_phylum, p_class, p_order, p_family, p_genus, p_species

plots <- list(cbb_phylum, cbb_depth, cbb_temp, cbb_sal)

combo <- compose_grid(
  plots = plots,
  ncol = 2,                 # 3 columns x 2 rows here
  title = "",
  tag_levels = "a",
  legend = "none",          # your bars have no legend
  file = "CBB_panels_patchwork.png",
  width = 180, height = 240 # tweak to your journal size
)

print(combo)

ggsave("CBB_prevalence_merge.png", plot = combo, width = 13, height = 18)

#########################################
# 
# plots <- list(cbb_class)
# 
# class <- compose_grid(
#   plots = plots,
#   ncol = 1,                 # 3 columns x 2 rows here
#   title = "",
#   tag_levels = "",
#   legend = "none",          # your bars have no legend
#   file = "CBB_panels_patchwork.png",
#   width = 180, height = 240 # tweak to your journal size
# )
# 
# print(combo)
# 
# ggsave("CBB_prevalence_class.svg", plot = class, width = 10, height = 16)
# 
# print(cbb_class)

# cbb_phylum <- plot_prevalence_by_taxon(
#   prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_phylum.csv",
#   comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/GTDB_Tk_Phylum_pairwise_prevalence_comparison.csv",
#   pathway        = "CBB",
#   taxon_col      = "GTDB_Tk_Phylum",
#   min_prev       = 0.01,
#   min_n_taxa     = 10,
#   label_mode     = "fraction",   # or "percent"
#   bar_fill       = "#4C72B0",    # change this freely for bar color
#   bracket_color  = "black",
#   outfile_prefix = "CBB_phylum_fraction"  # saves .png/.svg/.pdf
# )


##############################################################
##############################################################

# ==============================
# Reusable plotting function (descending order + aligned brackets)
# ==============================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(scales);  library(rlang); library(grid)
})

plot_prevalence_by_taxon <- function(
    prev_csv,
    comp_csv,
    pathway            = "CBB",
    taxon_col          = "GTDB_Tk_Phylum",
    min_prev           = 0.01,
    min_n_taxa         = 10,
    label_mode         = c("fraction","percent"),
    bar_fill           = "#1D3557",
    bracket_color      = "grey20",
    auto_ylim          = FALSE,
    y_limit            = c(0, 1),
    y_top_manual       = NULL,
    # spacing just above the highest involved bar (per bracket)
    top_pad            = 0.015,
    # vertical step added when two brackets overlap in x-range
    bracket_step       = 0.03,
    outfile_prefix     = NULL,
    width              = 8,
    height             = 8,
    dpi                = 900,
    save_formats       = c("png","svg","pdf"),
    sort_bars          = c("descending","ascending","none"),
    x_edge_pad         = 0.12,
    base_size          = 26,
    axis_text_x_size   = 24,
    axis_text_y_size   = 24,
    axis_title_size    = 26,
    bar_label_size     = 18,
    bar_label_vjust    = -0.45,
    # significance style
    sig_linewidth      = 1.0,
    # these are *maximums*; actual tick length is reduced if near a bar
    sig_tick_max       = 0.02,
    sig_tick_min       = 0.004,
    sig_label_nudge_y  = 0.01,
    sig_label_size     = 8,
    sig_label_fill     = NA,
    sig_label_alpha    = 1,
    sig_label_padding  = unit(0.08, "lines"),
    # split the bracket under the label
    sig_gap_frac       = 0.25,
    sig_gap_max        = 0.6
) {
  label_mode <- match.arg(label_mode)
  sort_bars  <- match.arg(sort_bars)
  
  df_prev <- read.csv(prev_csv, stringsAsFactors = FALSE)
  df_comp <- read.csv(comp_csv, stringsAsFactors = FALSE)
  
  taxon_sym <- rlang::sym(taxon_col)
  
  # filter + (re)order bars
  cbb_df2 <- df_prev %>%
    filter(Pathway == pathway,
           Prevalence >= min_prev,
           Number_of_taxa > min_n_taxa) %>%
    { if (sort_bars == "descending") arrange(., desc(Prevalence))
      else if (sort_bars == "ascending") arrange(., Prevalence) else . } %>%
    mutate(
      !!taxon_sym := as.character(!!taxon_sym),
      GroupVar    = !!taxon_sym,
      GroupVar    = factor(GroupVar, levels = unique(GroupVar))
    )
  
  stopifnot(nrow(cbb_df2) > 0)
  
  cbb_df2 <- cbb_df2 %>%
    mutate(
      bar_label = if (label_mode == "percent")
        percent(Prevalence, accuracy = 0.1)
      else sprintf("%.2f", Prevalence),
      xnum = as.numeric(GroupVar)
    )
  
  lvl_order <- levels(cbb_df2$GroupVar)
  
  # significance table -> numeric x aligned to bars
  STAT <- df_comp %>%
    filter(Pathway == pathway, Significant == "Yes") %>%
    tidyr::separate(Comparison, c("Group1","Group2"), sep = " vs ") %>%
    filter(Group1 %in% lvl_order, Group2 %in% lvl_order) %>%
    transmute(
      xmin = as.numeric(factor(Group1, levels = lvl_order)),
      xmax = as.numeric(factor(Group2, levels = lvl_order)),
      annotations = dplyr::case_when(
        p_value <= 0.001 ~ "***",
        p_value <= 0.01  ~ "**",
        p_value <= 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  
  # --- SMART BRACKET PLACEMENT ---
  max_prev <- max(cbb_df2$Prevalence, na.rm = TRUE)
  n_comp   <- nrow(STAT)
  
  # y-axis top
  if (!is.null(y_top_manual)) {
    y_top <- y_top_manual
  } else if (isTRUE(auto_ylim)) {
    y_top <- max(max_prev * 1.12, max_prev + top_pad + (max(0, n_comp-1)) * bracket_step + 0.02)
    y_top <- min(1, ceiling(y_top * 100) / 100)
    y_top <- max(y_top, 0.1)
  } else {
    y_top <- y_limit[2]
  }
  
  # Bar heights by x index
  bar_heights <- setNames(cbb_df2$Prevalence, cbb_df2$xnum)
  
  if (n_comp > 0) {
    # place shorter brackets first (closer to bars)
    STAT <- STAT %>% arrange((xmax - xmin), xmin)
    
    placed <- data.frame(xmin=numeric(), xmax=numeric(), y=numeric())
    y_vals <- numeric(nrow(STAT))
    local_max_vec <- numeric(nrow(STAT))
    
    for (i in seq_len(nrow(STAT))) {
      rng <- STAT$xmin[i]:STAT$xmax[i]
      local_max <- max(bar_heights[as.character(rng)], na.rm = TRUE)
      local_max_vec[i] <- local_max
      
      # start just above the tallest bar spanned
      y_candidate <- local_max + top_pad
      
      # if overlapping existing brackets at the same/too-close height, step up
      if (nrow(placed) > 0) {
        repeat {
          overlap <- which(!(STAT$xmax[i] < placed$xmin | STAT$xmin[i] > placed$xmax))
          if (length(overlap) == 0) break
          needed <- max(placed$y[overlap]) + bracket_step
          if (y_candidate < needed - 1e-9) y_candidate <- needed else break
        }
      }
      
      # cap to stay within axis
      y_vals[i] <- pmin(y_candidate, y_top - 0.01)
      placed <- rbind(placed, data.frame(xmin=STAT$xmin[i], xmax=STAT$xmax[i], y=y_vals[i]))
    }
    
    STAT$y <- y_vals
    
    # split under label (gap) and compute safe per-bracket tick length
    STAT$xmid <- (STAT$xmin + STAT$xmax)/2
    STAT$gap  <- pmin(sig_gap_max, sig_gap_frac * (STAT$xmax - STAT$xmin))
    STAT$xleft_end    <- pmax(STAT$xmin,  STAT$xmid - STAT$gap/2)
    STAT$xright_start <- pmin(STAT$xmax,  STAT$xmid + STAT$gap/2)
    
    # shorten tick so it never hits the local top bar
    STAT$local_max <- local_max_vec
    STAT$tick_down <- pmax(sig_tick_min,
                           pmin(sig_tick_max, 0.8 * (STAT$y - STAT$local_max - 0.002)))
  }
  
  # dynamic label nudge (keeps labels clear even when brackets are dense)
  sig_label_nudge_y_eff <- max(sig_label_nudge_y, 0.45 * bracket_step)
  
  # --- PLOT ---
  p <- ggplot(cbb_df2, aes(x = xnum, y = Prevalence)) +
    geom_col(width = 0.7, fill = bar_fill) +
    geom_text(aes(label = bar_label),
              vjust = bar_label_vjust, size = bar_label_size, fontface = "bold") +
    { if (n_comp > 0) geom_segment(
      data = STAT, aes(x = xmin, xend = xleft_end, y = y, yend = y),
      linewidth = sig_linewidth, lineend = "round", color = bracket_color
    ) } +
    { if (n_comp > 0) geom_segment(
      data = STAT, aes(x = xright_start, xend = xmax, y = y, yend = y),
      linewidth = sig_linewidth, lineend = "round", color = bracket_color
    ) } +
    { if (n_comp > 0) geom_segment(
      data = STAT, aes(x = xmin, xend = xmin, y = y, yend = y - tick_down),
      linewidth = sig_linewidth, lineend = "round", color = bracket_color
    ) } +
    { if (n_comp > 0) geom_segment(
      data = STAT, aes(x = xmax, xend = xmax, y = y, yend = y - tick_down),
      linewidth = sig_linewidth, lineend = "round", color = bracket_color
    ) } +
    { if (n_comp > 0) geom_label(
      data = STAT,
      aes(x = xmid, y = y + sig_label_nudge_y_eff, label = annotations),
      size = sig_label_size, fontface = "bold",
      fill = sig_label_fill, alpha = sig_label_alpha,
      label.size = 0, label.padding = sig_label_padding
    ) } +
    scale_x_continuous(
      breaks = seq_along(lvl_order),
      labels = lvl_order,
      expand = expansion(add = x_edge_pad)
    ) +
    scale_y_continuous(
      limits = c(0, y_top),
      labels = if (label_mode == "percent") percent_format(accuracy = 1) else label_number(accuracy = 0.01),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = base_size) +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                  size = axis_text_x_size, face = "bold"),
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_size, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line    = element_line(color = "black")
    ) +
    labs(x = NULL, y = if (label_mode == "percent") "Prevalence (%)" else "Prevalence")
  
  # saving
  if (!is.null(outfile_prefix)) {
    if ("png" %in% save_formats) ggsave(paste0(outfile_prefix, ".png"), plot = p, width = width, height = height, dpi = dpi)
    if ("svg" %in% save_formats) ggsave(paste0(outfile_prefix, ".svg"), plot = p, width = width, height = height)
    if ("pdf" %in% save_formats) {
      dev <- tryCatch(grDevices::cairo_pdf, error = function(e) grDevices::pdf)
      ggsave(paste0(outfile_prefix, ".pdf"), plot = p, width = width, height = height, device = dev)
    }
  }
  
  p
}

# #######################
# 
# cbb_phylum <- plot_prevalence_by_taxon(
#   prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_phylum.csv",
#   comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/GTDB_Tk_Phylum_pairwise_prevalence_comparison.csv",
#   pathway   = "CBB",
#   taxon_col = "GTDB_Tk_Phylum",
#   sort_bars = "descending",
#   
#   auto_ylim = TRUE,
#   top_pad = 0.020,          # start the first bracket a bit higher
#   bracket_step = 0.040,     # more vertical spacing between brackets
#   
#   sig_label_nudge_y = 0.016,# lift asterisks off the bracket line
#   sig_gap_frac = 0.35,      # bigger gap under the asterisk
#   sig_gap_max  = 0.9,       # allow wide splits for long brackets
#   sig_tick_max = 0.012,     # shorten ticks overall
#   sig_tick_min = 0.005,     # but keep a visible minimum
#   
#   bar_label_vjust = -0.35,  # keeps the bold “0.23/0.16…” from nudging into brackets
#   bar_label_size  = 16,     # slightly smaller label so it doesn’t dominate
#   x_edge_pad = 0.10         # a touch tighter on the sides
# )
# 
# cbb_depth <- plot_prevalence_by_taxon(
#   prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_depth.csv",
#   comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/depth_pairwise_prevalence_comparison.csv",
#   pathway   = "CBB",
#   taxon_col = "Depth_Category",
#   sort_bars = "descending",
#   
#   auto_ylim = TRUE,
#   top_pad = 0.015,
#   bracket_step = 0.030,
#   
#   sig_label_nudge_y = 0.014,
#   sig_gap_frac = 0.30,
#   sig_gap_max  = 0.7,
#   sig_tick_max = 0.012,
#   sig_tick_min = 0.005,
#   
#   bar_label_vjust = -0.40,
#   bar_label_size  = 16,
#   x_edge_pad = 0.12
# )
# 
# cbb_temp <- plot_prevalence_by_taxon(
#   prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_temp.csv",
#   comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/temp_pairwise_prevalence_comparison.csv",
#   pathway   = "CBB",
#   taxon_col = "temp_Category",
#   sort_bars = "descending",
#   auto_ylim = TRUE, top_pad = 0.015, bracket_step = 0.030,
#   sig_label_nudge_y = 0.014, sig_gap_frac = 0.30, sig_gap_max = 0.7,
#   sig_tick_max = 0.012, sig_tick_min = 0.005,
#   bar_label_vjust = -0.40, bar_label_size = 16, x_edge_pad = 0.12
# )
# 
# cbb_sal <- plot_prevalence_by_taxon(
#   prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_salinity.csv",
#   comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/salinity_pairwise_prevalence_comparison.csv",
#   pathway   = "CBB",
#   taxon_col = "Salinity_Category",
#   sort_bars = "descending",
#   auto_ylim = TRUE, top_pad = 0.015, bracket_step = 0.030,
#   sig_label_nudge_y = 0.014, sig_gap_frac = 0.30, sig_gap_max = 0.7,
#   sig_tick_max = 0.012, sig_tick_min = 0.005,
#   bar_label_vjust = -0.40, bar_label_size = 16, x_edge_pad = 0.12
# )
# 
# #######################
# 
# # Phylum
# cbb_phylum <- plot_prevalence_by_taxon(
#   prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_phylum.csv",
#   comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/GTDB_Tk_Phylum_pairwise_prevalence_comparison.csv",
#   pathway        = "CBB",
#   taxon_col      = "GTDB_Tk_Phylum",
#   min_prev       = 0.01,
#   min_n_taxa     = 10,
#   label_mode     = "fraction",
#   bar_fill       = "#4C72B0",
#   bracket_color  = "black",
#   auto_ylim      = TRUE,
#   outfile_prefix = "CBB_prevalence_phylum"
# )
# 
# # Depth
# cbb_depth <- plot_prevalence_by_taxon(
#   prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_depth.csv",
#   comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/depth_pairwise_prevalence_comparison.csv",
#   pathway        = "CBB",
#   taxon_col      = "Depth_Category",
#   min_prev       = 0.01,
#   min_n_taxa     = 10,
#   label_mode     = "fraction",
#   bar_fill       = "#4C72B0",
#   bracket_color  = "black",
#   auto_ylim      = TRUE,
#   outfile_prefix = "CBB_prevalence_depth"
# )
# 
# # Temperature
# cbb_temp <- plot_prevalence_by_taxon(
#   prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_temp.csv",
#   comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/temp_pairwise_prevalence_comparison.csv",
#   pathway        = "CBB",
#   taxon_col      = "temp_Category",
#   min_prev       = 0.01,
#   min_n_taxa     = 10,
#   label_mode     = "fraction",
#   bar_fill       = "#4C72B0",
#   bracket_color  = "black",
#   auto_ylim      = TRUE,
#   outfile_prefix = "CBB_prevalence_temp"
# )
# 
# # Salinity
# cbb_sal <- plot_prevalence_by_taxon(
#   prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_salinity.csv",
#   comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/salinity_pairwise_prevalence_comparison.csv",
#   pathway        = "CBB",
#   taxon_col      = "Salinity_Category",
#   min_prev       = 0.01,
#   min_n_taxa     = 10,
#   label_mode     = "fraction",
#   bar_fill       = "#4C72B0",
#   bracket_color  = "black",
#   auto_ylim      = TRUE,
#   outfile_prefix = "CBB_prevalence_salinity"
# )
# 
# # Ocean/Sea name
# cbb_ocean <- plot_prevalence_by_taxon(
#   prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_ocean_sea.csv",
#   comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/ocean_sea_name_pairwise_prevalence_comparison.csv",
#   pathway        = "CBB",
#   taxon_col      = "ocean_sea_name",
#   min_prev       = 0.01,
#   min_n_taxa     = 10,
#   label_mode     = "fraction",
#   bar_fill       = "#4C72B0",
#   bracket_color  = "black",
#   auto_ylim      = TRUE,
#   outfile_prefix = "CBB_prevalence_ocean"
# )
# 
# # IHO Sea
# cbb_sea <- plot_prevalence_by_taxon(
#   prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_IHOsea.csv",
#   comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/IHO_Sea_pairwise_prevalence_comparison.csv",
#   pathway        = "CBB",
#   taxon_col      = "IHO_Sea",
#   min_prev       = 0.03,
#   min_n_taxa     = 10,
#   label_mode     = "fraction",
#   bar_fill       = "#4C72B0",
#   bracket_color  = "black",
#   auto_ylim      = TRUE,
#   outfile_prefix = "CBB_prevalence_IHOSea"
# )
# 
# # If you still want extra standalone SVGs with custom sizes:
# ggsave("CBB_prevalence_ocean.svg", plot = cbb_ocean, width = 15, height = 17, units = "in", device = svglite::svglite)
# ggsave("CBB_prevalence_sea.svg",   plot = cbb_sea,   width = 15, height = 18, units = "in", device = svglite::svglite)
# 
# #######################################################
# # Step 1: Calculate prevalence per phylum
# prevalence_per_phylum <- data_table %>%
#   group_by(GTDB_Tk_Phylum) %>%
#   summarize(
#     Total_MAGs = n(),
#     Pathway_presence = sum(CBB), # Replace 'CBB' with the relevant column
#     Prevalence = Pathway_presence / Total_MAGs
#   )
# 
# # Step 2: Count the number of unique gOTUs per phylum
# gotu_counts_per_phylum <- data_table %>%
#   group_by(GTDB_Tk_Phylum) %>%
#   summarize(Number_of_gOTUs = n_distinct(OTU_cluster))
# 
# # Merge prevalence and gOTU counts
# phylum_data <- prevalence_per_phylum %>%
#   inner_join(gotu_counts_per_phylum, by = "GTDB_Tk_Phylum")
# 
# # Step 3: Calculate weighted prevalence for each phylum
# phylum_data <- phylum_data %>%
#   mutate(Weighted_Prevalence = Prevalence * Number_of_gOTUs)
# 
# # Step 4: Calculate total number of gOTUs across all phyla
# total_gOTUs <- sum(phylum_data$Number_of_gOTUs)
# 
# # Step 5: Calculate WAP
# WAP <- sum(phylum_data$Weighted_Prevalence) / total_gOTUs
# 
# # Write output to a file (optional)
# write.csv(phylum_data, "phylum_prevalence_output.csv", row.names = FALSE)
# 
# # Print the WAP
# print(paste("Weighted Average Prevalence (WAP):", WAP))
# 
####################################
####################################

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(scales);  library(rlang)
})

plot_prevalence_by_taxon <- function(
    prev_csv,
    comp_csv,
    pathway            = "CBB",
    taxon_col          = "GTDB_Tk_Phylum",
    min_prev           = 0.01,
    min_n_taxa         = 10,
    label_mode         = c("fraction","percent"),
    label_digits       = 2,      # NEW: decimals in bar labels
    bar_fill           = "#4C72B0",
    bracket_color      = "black",
    # axis & spacing
    auto_ylim          = TRUE,
    y_limit            = c(0, 1),
    y_top_manual       = NULL,
    x_edge_pad         = 0.10,
    # bracket placement (no overlaps with bars or each other)
    top_pad            = 0.024,
    bracket_step       = 0.055,
    bar_clearance      = 0.012,
    # stars
    star_above         = 0.012,
    star_gap_below     = 0.008,
    # ticks
    sig_linewidth      = 1.0,
    sig_tick_min       = 0.004,
    sig_tick_max       = 0.010,
    # text sizes
    base_size          = 27,
    axis_text_x_size   = 34,
    axis_text_y_size   = 34,
    axis_title_size    = 34,
    bar_label_size     = 13,
    bar_label_vjust    = -0.35,
    star_size          = 17.5,
    # saving
    outfile_prefix     = NULL,
    width              = 8, height = 8, dpi = 900,
    save_formats       = c("png","svg","pdf"),
    sort_bars          = c("descending","ascending","none")
){
  label_mode <- match.arg(label_mode)
  sort_bars  <- match.arg(sort_bars)
  
  df_prev <- read.csv(prev_csv, stringsAsFactors = FALSE)
  df_comp <- read.csv(comp_csv, stringsAsFactors = FALSE)
  taxon_sym <- rlang::sym(taxon_col)
  
  # ---- filter + order bars ----
  bars <- df_prev %>%
    filter(Pathway == pathway,
           Prevalence >= min_prev,
           Number_of_taxa > min_n_taxa) %>%
    { if (sort_bars == "descending") arrange(., desc(Prevalence))
      else if (sort_bars == "ascending") arrange(., Prevalence) else . } %>%
    mutate(
      !!taxon_sym := as.character(!!taxon_sym),
      GroupVar    = !!taxon_sym,
      GroupVar    = factor(GroupVar, levels = unique(!!taxon_sym))
    )
  stopifnot(nrow(bars) > 0)
  
  # --- bar labels: now controlled by label_digits ---
  bars <- bars %>%
    mutate(
      bar_label = if (label_mode == "percent")
        percent(Prevalence, accuracy = 10^(-label_digits))
      else
        sprintf(paste0("%.", label_digits, "f"), Prevalence),
      xnum = as.numeric(GroupVar)
    )
  lvl_order <- levels(bars$GroupVar)
  
  # ---- comparisons -> numeric x ----
  STAT <- df_comp %>%
    filter(Pathway == pathway, Significant == "Yes") %>%
    tidyr::separate(Comparison, c("Group1","Group2"), sep = " vs ") %>%
    filter(Group1 %in% lvl_order, Group2 %in% lvl_order) %>%
    transmute(
      xmin = as.numeric(factor(Group1, levels = lvl_order)),
      xmax = as.numeric(factor(Group2, levels = lvl_order)),
      annotations = case_when(
        p_value <= 0.001 ~ "***",
        p_value <= 0.01  ~ "**",
        p_value <= 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  if (nrow(STAT) == 0) {
    message("No significant comparisons; drawing bars only.")
    STAT <- NULL
  }
  
  max_prev <- max(bars$Prevalence, na.rm = TRUE)
  
  # ---- place brackets smartly (above tallest spanned bar; stack if overlap) ----
  if (!is.null(STAT)) {
    STAT <- STAT %>%
      mutate(span = abs(xmax - xmin)) %>%
      arrange(span, xmin)
    
    bar_heights <- setNames(bars$Prevalence, as.character(bars$xnum))
    placed <- data.frame(xmin=integer(0), xmax=integer(0), y=numeric(0))
    STAT$y <- NA_real_
    STAT$local_max <- NA_real_
    
    for (i in seq_len(nrow(STAT))) {
      rng <- seq(STAT$xmin[i], STAT$xmax[i])
      local_max <- max(bar_heights[as.character(rng)], na.rm = TRUE)
      y_cand <- local_max + top_pad
      
      if (nrow(placed) > 0) {
        repeat {
          ov <- with(placed, !(STAT$xmax[i] < xmin | STAT$xmin[i] > xmax))
          if (!any(ov)) break
          need <- max(placed$y[ov]) + bracket_step
          if (y_cand < need - 1e-9) y_cand <- need else break
        }
      }
      STAT$local_max[i] <- local_max
      STAT$y[i] <- y_cand
      placed <- rbind(placed, data.frame(xmin=STAT$xmin[i], xmax=STAT$xmax[i], y=STAT$y[i]))
    }
    
    # ticks with hard clearance from bar tops
    raw_tick <- 0.75*(STAT$y - STAT$local_max - 0.002)
    max_allowed <- STAT$y - (STAT$local_max + bar_clearance)
    STAT$tick_down <- pmax(sig_tick_min, pmin(sig_tick_max, pmin(raw_tick, max_allowed)))
    
    # star positions: above own line, but kept below next higher overlapping line
    STAT$xmid <- (STAT$xmin + STAT$xmax)/2
    STAT$label_y <- NA_real_
    ord <- order(STAT$y)  # bottom-up
    for (i in ord) {
      y_star <- STAT$y[i] + star_above
      blockers <- which(STAT$y > STAT$y[i] &
                          STAT$xmin <= STAT$xmid[i] &
                          STAT$xmax >= STAT$xmid[i])
      if (length(blockers)) y_star <- min(y_star, min(STAT$y[blockers]) - star_gap_below)
      STAT$label_y[i] <- y_star
    }
  }
  
  # ---- y axis top ----
  if (!is.null(y_top_manual)) {
    y_top <- y_top_manual
  } else if (isTRUE(auto_ylim)) {
    top_need <- if (!is.null(STAT)) max(STAT$label_y) + 0.03 else max_prev*1.12
    y_top <- max(top_need, max_prev*1.12)
    y_top <- min(1, ceiling(y_top*100)/100)
    y_top <- max(y_top, 0.1)
  } else {
    y_top <- y_limit[2]
  }
  
  # ---- plot ----
  p <- ggplot(bars, aes(x = as.numeric(GroupVar), y = Prevalence)) +
    geom_col(width = 0.7, fill = bar_fill) +
    geom_text(aes(label = bar_label),
              vjust = bar_label_vjust, size = bar_label_size) +
    { if (!is.null(STAT)) geom_segment(
      data = STAT, aes(x = xmin, xend = xmax, y = y, yend = y),
      linewidth = sig_linewidth, lineend = "round", color = bracket_color
    ) } +
    { if (!is.null(STAT)) geom_segment(
      data = STAT, aes(x = xmin, xend = xmin, y = y, yend = y - tick_down),
      linewidth = sig_linewidth, lineend = "round", color = bracket_color
    ) } +
    { if (!is.null(STAT)) geom_segment(
      data = STAT, aes(x = xmax, xend = xmax, y = y, yend = y - tick_down),
      linewidth = sig_linewidth, lineend = "round", color = bracket_color
    ) } +
    { if (!is.null(STAT)) geom_text(
      data = STAT, aes(x = xmid, y = label_y, label = annotations),
      size = star_size, fontface = "bold"
    ) } +
    scale_x_continuous(
      breaks = seq_along(lvl_order),
      labels = lvl_order,
      expand = expansion(add = x_edge_pad)
    ) +
    scale_y_continuous(
      limits = c(0, y_top),
      labels = if (label_mode == "percent") percent_format(accuracy = 1)
      else label_number(accuracy = 0.01),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = base_size) +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                  size = axis_text_x_size, face = "bold"),
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_size, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line    = element_line(color = "black")
    ) +
    labs(x = NULL, y = if (label_mode == "percent") "Prevalence (%)" else "Prevalence")
  
  # save if requested
  if (!is.null(outfile_prefix)) {
    if ("png" %in% save_formats)
      ggsave(paste0(outfile_prefix, ".png"), plot = p, width = width, height = height, dpi = dpi)
    if ("svg" %in% save_formats)
      ggsave(paste0(outfile_prefix, ".svg"), plot = p, width = width, height = height)
    if ("pdf" %in% save_formats) {
      dev <- tryCatch(grDevices::cairo_pdf, error = function(e) grDevices::pdf)
      ggsave(paste0(outfile_prefix, ".pdf"), plot = p, width = width, height = height, device = dev)
    }
  }
  
  p
}

##################################################################################################

cbb_phylum <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_phylum.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/GTDB_Tk_Phylum_pairwise_prevalence_comparison.csv",
  pathway   = "CBB",
  taxon_col = "GTDB_Tk_Phylum",
  min_prev  = 0.01,
  min_n_taxa = 10,
  sort_bars = "descending",
  auto_ylim = TRUE,
  label_mode = "fraction",
  label_digits = 3,        # << three decimals in the bar labels
  
  # spacing tuned for dense brackets
  top_pad = 0.032,
  bracket_step = 0.040,
  bar_clearance = 0.012,
  
  # fine look/feel
  sig_tick_min = 0.004,
  sig_tick_max = 0.008,
  star_above   = 0.012,
  star_gap_below = 0.006,
  
  # optional polish
  x_edge_pad = 0.10,
  bar_label_size = 14,
  bar_label_vjust = -0.30
)

ggsave("CBB_prevalence_phylum.svg", plot = cbb_phylum, width = 13, height = 19.8)

###############################################################################################################

cbb_depth <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_depth.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/depth_pairwise_prevalence_comparison.csv",
  pathway   = "CBB",
  taxon_col = "Depth_Category",
  min_prev  = 0.01,
  min_n_taxa = 10,
  sort_bars = "descending",
  auto_ylim = TRUE,
  label_mode = "fraction",
  label_digits = 3,        # << three decimals in the bar labels
    
  # spacing tuned for dense brackets
  top_pad = 0.032,
  bracket_step = 0.040,
  bar_clearance = 0.012,
  
  # fine look/feel
  sig_tick_min = 0.004,
  sig_tick_max = 0.008,
  star_above   = 0.012,
  star_gap_below = 0.006,
  
  # optional polish
  x_edge_pad = 0.10,
  bar_label_size = 14,
  bar_label_vjust = -0.30
)

ggsave("CBB_prevalence_depth.svg", plot = cbb_depth, width = 13, height = 19.8)

###############################################################################################################

cbb_temp <- plot_prevalence_by_taxon(
  prev_csv  = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_temp.csv",
  comp_csv  = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/temp_pairwise_prevalence_comparison.csv",
  pathway   = "CBB",
  taxon_col = "temp_Category",
  min_prev  = 0.01,
  min_n_taxa = 10,
  sort_bars = "descending",
  auto_ylim = TRUE,
  label_mode = "fraction",
  label_digits = 3,        # << three decimals in the bar labels
  
  # spacing tuned for dense brackets
  top_pad = 0.032,
  bracket_step = 0.040,
  bar_clearance = 0.012,
  
  # fine look/feel
  sig_tick_min = 0.004,
  sig_tick_max = 0.008,
  star_above   = 0.012,
  star_gap_below = 0.006,
  
  # optional polish
  x_edge_pad = 0.10,
  bar_label_size = 14,
  bar_label_vjust = -0.30
)

ggsave("CBB_prevalence_temp.svg", plot = cbb_temp, width = 13, height = 14.8)

###############################################################################################################

cbb_sal <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_salinity.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/salinity_pairwise_prevalence_comparison.csv",
  pathway   = "CBB",
  taxon_col = "Salinity_Category",
  min_prev  = 0.01,
  min_n_taxa = 10,
  sort_bars = "descending",
  auto_ylim = TRUE,
  label_mode = "fraction",
  label_digits = 3,        # << three decimals in the bar labels
    
  # spacing tuned for dense brackets
  top_pad = 0.032,
  bracket_step = 0.040,
  bar_clearance = 0.012,
  
  # fine look/feel
  sig_tick_min = 0.004,
  sig_tick_max = 0.008,
  star_above   = 0.012,
  star_gap_below = 0.006,
  
  # optional polish
  x_edge_pad = 0.10,
  bar_label_size = 14,
  bar_label_vjust = -0.30
)

ggsave("CBB_prevalence_salinity.svg", plot = cbb_sal, width = 13, height = 14.8)

###############################################################################################################

cbb_ocean <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_ocean_sea.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/ocean_sea_name_pairwise_prevalence_comparison.csv",
  pathway   = "CBB",
  taxon_col = "ocean_sea_name",
  min_prev  = 0.01,
  min_n_taxa = 10,
  sort_bars = "descending",
  auto_ylim = TRUE,
  label_mode = "fraction",
  label_digits = 3,        # << three decimals in the bar labels
    
  # more separation for the big fan-out from the top bar
  top_pad       = 0.032,
  bracket_step  = 0.048,
  bar_clearance = 0.012,
  
  # keep stars distinct
  star_above     = 0.010,
  star_gap_below = 0.006,
  star_size      = 15.8,
  
  # tidy ticks
  sig_tick_min = 0.004,
  sig_tick_max = 0.006,
  
  # minor polish
  x_edge_pad      = 0.10,
  bar_label_size  = 14,
  bar_label_vjust = -0.30
)

ggsave("CBB_prevalence_ocean.svg", plot = cbb_ocean, width = 25, height = 40.8)

# # --- dimensions: keep tall aspect but journal-friendly size ---
# width_in  <- 15.5   # ~two-column width
# height_in <- 29.0  # tall panel; matches your layout
# 
# # PDF (Cairo if available for crisp text)
# pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
# ggsave("CBB_prevalence_ocean.pdf",
#        plot = cbb_ocean, width = width_in, height = height_in, units = "in",
#        device = pdf_dev)
# 
# # SVG (vector; scales cleanly)
# ggsave("CBB_prevalence_ocean.svg",
#        plot = cbb_ocean, width = width_in, height = height_in, units = "in")
# 
# # TIFF — 1200 dpi, LZW compression (via ragg; robust at large sizes)
# if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")
# ggsave("CBB_prevalence_ocean.tif",
#        plot = cbb_ocean,
#        width = width_in, height = height_in, units = "in",
#        device = ragg::agg_tiff,
#        res = 1200,
#        compression = "lzw")

###############################################################################################################

cbb_sea <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_IHOsea.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/IHO_Sea_pairwise_prevalence_comparison.csv",
  pathway   = "CBB",
  taxon_col = "IHO_Sea",
  min_prev  = 0.03,
  min_n_taxa = 10,
  sort_bars = "descending",
  auto_ylim = TRUE,
  label_mode = "fraction",
  label_digits = 3,
  
  # more separation for the big fan-out from the top bar
  top_pad       = 0.028,
  bracket_step  = 0.048,
  bar_clearance = 0.012,
  
  # keep stars distinct
  star_above     = 0.010,
  star_gap_below = 0.006,
  star_size      = 15.8,
  
  # tidy ticks
  sig_tick_min = 0.004,
  sig_tick_max = 0.006,
  
  # minor polish
  x_edge_pad      = 0.10,
  bar_label_size  = 14,
  bar_label_vjust = -0.30
)

ggsave("CBB_prevalence_sea.svg", plot = cbb_sea, width = 25, height = 32.8)

###############################################################################################################

cbb_class <- plot_prevalence_by_taxon(
  prev_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/input_files/pathway_class.csv",
  comp_csv       = "D:/Thesis/work_package_3/Data_analysis/01_figures/01_prevalence/output_figures/class_pairwise_prevalence_comparison.csv",
  pathway   = "CBB",
  taxon_col = "GTDB_Tk_Class",
  min_prev  = 0.03,
  min_n_taxa = 10,
  sort_bars = "descending",
  auto_ylim = TRUE,
  label_mode = "fraction",
  label_digits = 3,

  # more separation for the big fan-out from the top bar
  top_pad       = 0.018,
  bracket_step  = 0.068,
  bar_clearance = 0.022,

  # keep stars distinct
  star_above     = 0.010,
  star_gap_below = 0.006,
  star_size      = 20.8,

  # tidy ticks
  sig_tick_min = 0.004,
  sig_tick_max = 0.006,

  # minor polish
  x_edge_pad      = 0.10,
  bar_label_size  = 14,
  bar_label_vjust = -0.30
)

ggsave("CBB_prevalence_class.svg", plot = cbb_class, width = 30, height = 40)

###############################################################################################################

cbb_class <- plot_prevalence_by_taxon(

  pathway        = "CBB",
  taxon_col      = "",
  min_prev       = 0.025,
  min_n_taxa     = 10,
  label_mode     = "fraction",   # or "percent"
  bar_fill       = "#4C72B0",    # change this freely for bar color
  bracket_color  = "black",
  outfile_prefix = "CBB_prevalence_class"  # saves .png/.svg/.pdf
)

ggsave("CBB_prevalence_class.svg", plot = cbb_class, width = 25, height = 39)



#######################################

# install.packages("patchwork") # if needed
library(patchwork)

compose_grid <- function(
    plots,
    ncol = 2,
    nrow = NULL,
    title = NULL,
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")",
    legend = c("keep","collect","none"),
    file = NULL,           # e.g., "CBB_panels_patchwork.png"
    width_mm = 180,        # typical 2-column width ~ 174–180 mm
    height_mm = 240,
    dpi = 900
) {
  legend <- match.arg(legend)
  
  p <- wrap_plots(plotlist = plots, ncol = ncol, nrow = nrow,
                  guides = if (legend == "collect") "collect" else "keep") +
    plot_annotation(
      title = title,
      tag_levels = tag_levels,
      tag_prefix = tag_prefix,
      tag_suffix = tag_suffix
    )
  
  if (legend == "none") p <- p & theme(legend.position = "none")
  
  if (!is.null(file)) {
    # PNG via ragg
    if (grepl("\\.png$", file, ignore.case = TRUE)) {
      ggsave(file, p,
             width = width_mm, height = height_mm, units = "mm",
             dpi = dpi, device = ragg::agg_png)
    } else if (grepl("\\.svg$", file, ignore.case = TRUE)) {
      ggsave(file, p,
             width = width_mm, height = height_mm, units = "mm",
             device = svglite::svglite)
    } else if (grepl("\\.pdf$", file, ignore.case = TRUE)) {
      dev_fun <- tryCatch(grDevices::cairo_pdf, error = function(e) grDevices::pdf)
      ggsave(file, p,
             width = width_mm, height = height_mm, units = "mm",
             device = dev_fun)
    } else {
      ggsave(file, p, width = width_mm, height = height_mm, units = "mm", dpi = dpi)
    }
  }
  
  p
}

plots <- list(cbb_phylum, cbb_depth, cbb_temp, cbb_sal)

combo <- compose_grid(
  plots = plots,
  ncol = 2,
  title = "",
  tag_levels = "a",
  legend = "none",
  file = "CBB_panels_patchwork.png",
  width_mm = 180, height_mm = 240, dpi = 900
)

ggsave("CBB_prevalence_merge.svg", plot = combo, width = 25, height = 32.8)

# # Optional additional export
# ggsave("CBB_prevalence_merge.png", plot = combo, width = 24, height = 30,
#        units = "in", dpi = 900, device = ragg::agg_png)



