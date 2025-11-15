#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/04_alpha_diversity")

# Load required libraries
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(reshape2)
library(tibble)  

#################################################################################
###### depth
# Load required packages
library(tidyverse)
library(mvabund)

# STEP 1: Subset and clean your data

# Read data
merged_data <- read.csv("merged_data.csv")
bin_libs <- read.csv("bin_libs.csv")

#merging the OTU and metadata 
merged_data <- merge(merged_data, bin_libs, by.x = "Genome", by.y = "bin", all.x = TRUE)


subset_df <- merged_data %>%
  select(library, GTDB_Tk_Phylum, Depth_Category) %>%
  filter(!is.na(library) & !is.na(GTDB_Tk_Phylum) & !is.na(Depth_Category))

# STEP 2: Build sample-by-phyla presence matrix
presence_matrix <- subset_df %>%
  distinct(library, GTDB_Tk_Phylum) %>%       # 1 row per taxon per sample
  mutate(presence = 1) %>%                    # tag presence
  pivot_wider(names_from = GTDB_Tk_Phylum,
              values_from = presence,
              values_fill = 0)                # fill absences with 0

# STEP 3: Add Depth_Category metadata
metadata <- subset_df %>%
  distinct(library, Depth_Category)

# Merge matrix and metadata
full_data <- presence_matrix %>%
  left_join(metadata, by = "library") %>%
  column_to_rownames("library")

# STEP 4: Separate community matrix and metadata
comm_data <- full_data %>% select(-Depth_Category)
depth_group <- full_data$Depth_Category

# STEP 5: Run mvabund analysis (many GLMs with binomial family)
comm_mv <- mvabund(comm_data)
mv_model <- manyglm(comm_mv ~ depth_group, family = "binomial")

# # STEP 6: Statistical test (overall + per-taxon effects)
# depth_anova_result <- anova.manyglm(mv_model, p.uni = "adjusted")

# STEP 6: Statistical test with proper univariate p-values
depth_anova_result_01 <- anova.manyglm(
  mv_model,
  p.uni = "adjusted",          # get adjusted univariate p-values
  resamp = "montecarlo",       # use classic permutation, not PIT-trap
  show.warning = TRUE
)

# Output results
#print(depth_anova_result)
print(depth_anova_result_01)

# Extract and transpose
test_stat <- t(depth_anova_result_01$uni.test["depth_group", ])
p_val <- t(depth_anova_result_01$uni.p["depth_group", ])

# Combine into a tidy table
univariate_df <- data.frame(
  Taxon = colnames(depth_anova_result_01$uni.test),
  Deviance = as.numeric(test_stat),
  P_value = as.numeric(p_val)
)

# Optional: Sort by p-value
univariate_df <- univariate_df[order(univariate_df$P_value), ]

# Save to CSV
write.csv(
  univariate_df,
  "depth_univariate_results.csv",
  row.names = FALSE
)

library(writexl)

write_xlsx(
  list(
    Multivariate = depth_anova_result_01$table,
    Univariate = univariate_df
  ),
  "depth_anova_summary.xlsx"
)

######################
###### salinity ###

# Load required packages
library(tidyverse)
library(mvabund)

# STEP 1: Subset and clean your data

# Read data
merged_data <- read.csv("merged_data.csv")
bin_libs <- read.csv("bin_libs.csv")

#merging the OTU and metadata 
merged_data <- merge(merged_data, bin_libs, by.x = "Genome", by.y = "bin", all.x = TRUE)

subset_df <- merged_data %>%
  select(library, GTDB_Tk_Phylum, Salinity_Category) %>%
  filter(!is.na(library) & !is.na(GTDB_Tk_Phylum) & !is.na(Salinity_Category))

# STEP 2: Build sample-by-phyla presence matrix
presence_matrix <- subset_df %>%
  distinct(library, GTDB_Tk_Phylum) %>%       # 1 row per taxon per sample
  mutate(presence = 1) %>%                    # tag presence
  pivot_wider(names_from = GTDB_Tk_Phylum,
              values_from = presence,
              values_fill = 0)                # fill absences with 0

# STEP 3: Add Depth_Category metadata
metadata <- subset_df %>%
  distinct(library, Salinity_Category)

# Merge matrix and metadata
full_data <- presence_matrix %>%
  left_join(metadata, by = "library") %>%
  column_to_rownames("library")

# STEP 4: Separate community matrix and metadata
comm_data <- full_data %>% select(-Salinity_Category)
salinity_group <- full_data$Salinity_Category

# STEP 5: Run mvabund analysis (many GLMs with binomial family)
comm_mv <- mvabund(comm_data)
mv_model <- manyglm(comm_mv ~ salinity_group, family = "binomial")

# STEP 6: Statistical test (overall + per-taxon effects)
#salinity_anova_result <- anova.manyglm(mv_model, p.uni = "adjusted")

# STEP 6: Statistical test with proper univariate p-values
salinity_anova_result_01 <- anova.manyglm(
  mv_model,
  p.uni = "adjusted",          # get adjusted univariate p-values
  resamp = "montecarlo",       # use classic permutation, not PIT-trap
  show.warning = TRUE
)

# Output results
#print(salinity_anova_result)
print(salinity_anova_result_01)

# Extract and transpose
test_stat <- t(salinity_anova_result_01$uni.test["salinity_group", ])
p_val <- t(salinity_anova_result_01$uni.p["salinity_group", ])

# Combine into a tidy table
univariate_df <- data.frame(
  Taxon = colnames(salinity_anova_result_01$uni.test),
  Deviance = as.numeric(test_stat),
  P_value = as.numeric(p_val)
)

# Optional: Sort by p-value
univariate_df <- univariate_df[order(univariate_df$P_value), ]

# Save to CSV
write.csv(
  univariate_df,
  "salinity_univariate_results.csv",
  row.names = FALSE
)

library(writexl)

write_xlsx(
  list(
    Multivariate = salinity_anova_result_01$table,
    Univariate = univariate_df
  ),
  "salinity_anova_summary.xlsx"
)


######################
###### temperature ###

# Load required packages
library(tidyverse)
library(mvabund)

# STEP 1: Subset and clean your data

# Read data
merged_data <- read.csv("merged_data.csv")
bin_libs <- read.csv("bin_libs.csv")

#merging the OTU and metadata 
merged_data <- merge(merged_data, bin_libs, by.x = "Genome", by.y = "bin", all.x = TRUE)

subset_df <- merged_data %>%
  select(library, GTDB_Tk_Phylum, temp_Category) %>%
  filter(!is.na(library) & !is.na(GTDB_Tk_Phylum) & !is.na(temp_Category))

# STEP 2: Build sample-by-phyla presence matrix
presence_matrix <- subset_df %>%
  distinct(library, GTDB_Tk_Phylum) %>%       # 1 row per taxon per sample
  mutate(presence = 1) %>%                    # tag presence
  pivot_wider(names_from = GTDB_Tk_Phylum,
              values_from = presence,
              values_fill = 0)                # fill absences with 0

# STEP 3: Add Depth_Category metadata
metadata <- subset_df %>%
  distinct(library, temp_Category)

# Merge matrix and metadata
full_data <- presence_matrix %>%
  left_join(metadata, by = "library") %>%
  column_to_rownames("library")

# STEP 4: Separate community matrix and metadata
comm_data <- full_data %>% select(-temp_Category)
temp_group <- full_data$temp_Category

# STEP 5: Run mvabund analysis (many GLMs with binomial family)
comm_mv <- mvabund(comm_data)
mv_model <- manyglm(comm_mv ~ temp_group, family = "binomial")

# # STEP 6: Statistical test (overall + per-taxon effects)
# temp_anova_result <- anova.manyglm(mv_model, p.uni = "adjusted")

# STEP 6: Statistical test with proper univariate p-values
temp_anova_result_01 <- anova.manyglm(
  mv_model,
  p.uni = "adjusted",          # get adjusted univariate p-values
  resamp = "montecarlo",       # use classic permutation, not PIT-trap
  show.warning = TRUE
)

# Output results
#print(temp_anova_result)
print(temp_anova_result_01)


# Extract and transpose
test_stat <- t(temp_anova_result_01$uni.test["temp_group", ])
p_val <- t(temp_anova_result_01$uni.p["temp_group", ])

# Combine into a tidy table
univariate_df <- data.frame(
  Taxon = colnames(temp_anova_result_01$uni.test),
  Deviance = as.numeric(test_stat),
  P_value = as.numeric(p_val)
)

# Optional: Sort by p-value
univariate_df <- univariate_df[order(univariate_df$P_value), ]

# Save to CSV
write.csv(
  univariate_df,
  "temp_univariate_results.csv",
  row.names = FALSE
)

library(writexl)

write_xlsx(
  list(
    Multivariate = temp_anova_result_01$table,
    Univariate = univariate_df
  ),
  "temp_anova_summary.xlsx"
)

########################################
######################
###### oceans#

# Load required packages
library(tidyverse)
library(mvabund)

# STEP 1: Subset and clean your data

# Read data
merged_data <- read.csv("merged_data.csv")
bin_libs <- read.csv("bin_libs.csv")

#merging the OTU and metadata 
merged_data <- merge(merged_data, bin_libs, by.x = "Genome", by.y = "bin", all.x = TRUE)

# subset_df <- merged_data %>%
#   select(library, GTDB_Tk_Phylum, ocean_sea_name) %>%
#   filter(!is.na(library) & !is.na(GTDB_Tk_Phylum) & !is.na(ocean_sea_name))

subset_df <- merged_data %>%
  select(library, GTDB_Tk_Phylum, ocean_sea_name) %>%
  filter(
    !is.na(library) & library != "Unknown",
    !is.na(GTDB_Tk_Phylum) & GTDB_Tk_Phylum != "Unknown",
    !is.na(ocean_sea_name) & ocean_sea_name != "Unknown"
  )



# STEP 2: Build sample-by-phyla presence matrix
presence_matrix <- subset_df %>%
  distinct(library, GTDB_Tk_Phylum) %>%       # 1 row per taxon per sample
  mutate(presence = 1) %>%                    # tag presence
  pivot_wider(names_from = GTDB_Tk_Phylum,
              values_from = presence,
              values_fill = 0)                # fill absences with 0

# STEP 3: Add Depth_Category metadata
metadata <- subset_df %>%
  distinct(library, ocean_sea_name)

# Merge matrix and metadata
full_data <- presence_matrix %>%
  left_join(metadata, by = "library") %>%
  column_to_rownames("library")

# STEP 4: Separate community matrix and metadata
comm_data <- full_data %>% select(-ocean_sea_name)
temp_group <- full_data$ocean_sea_name

# STEP 5: Run mvabund analysis (many GLMs with binomial family)
comm_mv <- mvabund(comm_data)
mv_model <- manyglm(comm_mv ~ temp_group, family = "binomial")

# # STEP 6: Statistical test (overall + per-taxon effects)
# temp_anova_result <- anova.manyglm(mv_model, p.uni = "adjusted")

# STEP 6: Statistical test with proper univariate p-values
ocean_anova_result_01 <- anova.manyglm(
  mv_model,
  p.uni = "adjusted",          # get adjusted univariate p-values
  resamp = "montecarlo",       # use classic permutation, not PIT-trap
  show.warning = TRUE
)

# Output results
#print(temp_anova_result)
print(ocean_anova_result_01)


# Extract and transpose
test_stat <- t(ocean_anova_result_01$uni.test["temp_group", ])
p_val <- t(ocean_anova_result_01$uni.p["temp_group", ])

# Combine into a tidy table
univariate_df <- data.frame(
  Taxon = colnames(ocean_anova_result_01$uni.test),
  Deviance = as.numeric(test_stat),
  P_value = as.numeric(p_val)
)

# Optional: Sort by p-value
univariate_df <- univariate_df[order(univariate_df$P_value), ]

# Save to CSV
write.csv(
  univariate_df,
  "ocean_univariate_results.csv",
  row.names = FALSE
)

library(writexl)

write_xlsx(
  list(
    Multivariate = ocean_anova_result_01$table,
    Univariate = univariate_df
  ),
  "ocean_anova_summary.xlsx"
)


#################################################
######################
###### sea #

# Load required packages
library(tidyverse)
library(mvabund)

# STEP 1: Subset and clean your data

# Read data
merged_data <- read.csv("merged_data.csv")
bin_libs <- read.csv("bin_libs.csv")

#merging the OTU and metadata 
merged_data <- merge(merged_data, bin_libs, by.x = "Genome", by.y = "bin", all.x = TRUE)

# subset_df <- merged_data %>%
#   select(library, GTDB_Tk_Phylum, ocean_sea_name) %>%
#   filter(!is.na(library) & !is.na(GTDB_Tk_Phylum) & !is.na(ocean_sea_name))

subset_df <- merged_data %>%
  select(library, GTDB_Tk_Phylum, IHO_Sea) %>%
  filter(
    !is.na(library) & library != "Unknown",
    !is.na(GTDB_Tk_Phylum) & GTDB_Tk_Phylum != "Unknown",
    !is.na(IHO_Sea) & IHO_Sea != "Unknown"
  )



# STEP 2: Build sample-by-phyla presence matrix
presence_matrix <- subset_df %>%
  distinct(library, GTDB_Tk_Phylum) %>%       # 1 row per taxon per sample
  mutate(presence = 1) %>%                    # tag presence
  pivot_wider(names_from = GTDB_Tk_Phylum,
              values_from = presence,
              values_fill = 0)                # fill absences with 0

# STEP 3: Add Depth_Category metadata
metadata <- subset_df %>%
  distinct(library, IHO_Sea)

# Merge matrix and metadata
full_data <- presence_matrix %>%
  left_join(metadata, by = "library") %>%
  column_to_rownames("library")

# STEP 4: Separate community matrix and metadata
comm_data <- full_data %>% select(-IHO_Sea)
temp_group <- full_data$IHO_Sea

# STEP 5: Run mvabund analysis (many GLMs with binomial family)
comm_mv <- mvabund(comm_data)
mv_model <- manyglm(comm_mv ~ temp_group, family = "binomial")

# # STEP 6: Statistical test (overall + per-taxon effects)
# temp_anova_result <- anova.manyglm(mv_model, p.uni = "adjusted")

# STEP 6: Statistical test with proper univariate p-values
sea_anova_result_01 <- anova.manyglm(
  mv_model,
  p.uni = "adjusted",          # get adjusted univariate p-values
  resamp = "montecarlo",       # use classic permutation, not PIT-trap
  show.warning = TRUE
)

# Output results
#print(temp_anova_result)
print(sea_anova_result_01)


# Extract and transpose
test_stat <- t(sea_anova_result_01$uni.test["temp_group", ])
p_val <- t(sea_result_01$uni.p["temp_group", ])

# Combine into a tidy table
univariate_df <- data.frame(
  Taxon = colnames(sea_anova_result_01$uni.test),
  Deviance = as.numeric(test_stat),
  P_value = as.numeric(p_val)
)

# Optional: Sort by p-value
univariate_df <- univariate_df[order(univariate_df$P_value), ]

# Save to CSV
write.csv(
  univariate_df,
  "sea_univariate_results.csv",
  row.names = FALSE
)

library(writexl)

write_xlsx(
  list(
    Multivariate = sea_anova_result_01$table,
    Univariate = univariate_df
  ),
  "sea_anova_summary.xlsx"
)


#################################################

#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/03_general_biogeography")

# Load required libraries
library(Hmisc)
library(ggplot2)
library(dplyr)


# Read data
merged_data <- read.csv("merged_data.csv")

# Remove unknown salinity categories
merged_data <- merged_data %>% filter(temp_Category != "Unknown")

# Ensure `occurence` is numeric
merged_data$occurence <- as.numeric(merged_data$occurence)

# Summarize data by Salinity_Category and GTDB_Tk_Phylum
summarized_data <- merged_data %>%
  group_by(temp_Category, GTDB_Tk_Phylum) %>%
  summarise(pathway_present = sum(occurence, na.rm = TRUE), .groups = "drop") %>%
  filter(!temp_Category %in% c("", NA))  # Remove empty or NA values

# Print results
print(summarized_data)


# Apply log10 transformation to pathway_present (to compress wide ranges)
summarized_data <- summarized_data %>%
  mutate(pathway_present_log = log10(pathway_present + 1))  # log(1 + value) to avoid log(0)

# Convert to wide format for clustering
library(tidyr)
data_wide <- summarized_data %>%
  pivot_wider(names_from = temp_Category, values_from = pathway_present_log, values_fill = 0)

# Perform hierarchical clustering on phyla
distance_matrix <- dist(data_wide[,-1])  # Compute distance matrix (excluding phylum names)
clustering <- hclust(distance_matrix, method = "ward.D2")  # Hierarchical clustering

# Order phyla based on clustering
ordered_phyla <- unique(data_wide$GTDB_Tk_Phylum[clustering$order])
summarized_data$GTDB_Tk_Phylum <- factor(summarized_data$GTDB_Tk_Phylum, 
                                         levels = unique(ordered_phyla))

library(ggplot2)
library(dplyr)
library(forcats)

# Clean & order
plot_df <- summarized_data %>%
  filter(!is.na(GTDB_Tk_Phylum), GTDB_Tk_Phylum != "") %>%
  mutate(
    # (Optional) enforce your preferred temperature order:
    # temp_Category = fct_relevel(temp_Category, "Cold", "Temperate", "Warm", "Hot"),
    GTDB_Tk_Phylum = fct_drop(as.factor(GTDB_Tk_Phylum)),
    GTDB_Tk_Phylum = fct_reorder(GTDB_Tk_Phylum, pathway_present_log,
                                 .fun = median, .desc = TRUE)
  )

# Plot (Temperature on x, Phylum on y) — legend on the RIGHT
p <- ggplot(plot_df, aes(x = temp_Category, y = GTDB_Tk_Phylum)) +
  geom_tile(
    aes(fill = pathway_present_log),
    color = "white",
    size = 0.3,
    width = 0.98,
    height = 0.98
  ) +
  scale_x_discrete(expand = c(0, 0), drop = TRUE) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(
    option = "magma", direction = -1, na.value = "grey90",
    name = "Pathway presence (log)"
  ) +
  coord_fixed() +
  labs(title = NULL, x = "Temperature", y = "Phylum") +
  theme_minimal(base_size = 22, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 24),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.text = element_text(color = "black"),
    axis.ticks.length = unit(3, "mm"),
    plot.margin = margin(20, 20, 20, 20, "pt"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(face = "bold", size = 22),
    legend.text  = element_text(size = 20),
    legend.key.height = unit(6, "mm"),
    legend.key.width  = unit(10, "mm")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      # Right-side legend: tall & thin bar
      barheight = unit(120, "mm"),
      barwidth  = unit(6, "mm"),
      ticks = TRUE,
      label.theme = element_text(size = 18),
      title.theme = element_text(face = "bold", size = 22)
    )
  )

p

# SVG export — add some width to make space for the right legend
ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/03_general_biogeography/01_temp_category_phylum.svg",
  plot = p, device = "svg",
  width = 420, height = 600, units = "mm", dpi = 1200
)


##########################################

# Read data
merged_data <- read.csv("merged_data.csv")

# Remove unknown salinity categories
merged_data <- merged_data %>% filter(Salinity_Category != "Unknown")

# Ensure `occurence` is numeric
merged_data$occurence <- as.numeric(merged_data$occurence)

# Summarize data by Salinity_Category and GTDB_Tk_Phylum
summarized_data <- merged_data %>%
  group_by(Salinity_Category, GTDB_Tk_Phylum) %>%
  summarise(pathway_present = sum(occurence, na.rm = TRUE), .groups = "drop") %>%
  filter(!Salinity_Category %in% c("", NA))  # Remove empty or NA values

# Print results
print(summarized_data)

# Apply log10 transformation to pathway_present (to compress wide ranges)
summarized_data <- summarized_data %>%
  mutate(pathway_present_log = log10(pathway_present + 1))  # log(1 + value) to avoid log(0)

# Convert to wide format for clustering
library(tidyr)
data_wide <- summarized_data %>%
  pivot_wider(names_from = Salinity_Category, values_from = pathway_present_log, values_fill = 0)

# Perform hierarchical clustering on phyla
distance_matrix <- dist(data_wide[,-1])  # Compute distance matrix (excluding phylum names)
clustering <- hclust(distance_matrix, method = "ward.D2")  # Hierarchical clustering

# Order phyla based on clustering
ordered_phyla <- unique(data_wide$GTDB_Tk_Phylum[clustering$order])
summarized_data$GTDB_Tk_Phylum <- factor(summarized_data$GTDB_Tk_Phylum, 
                                         levels = unique(ordered_phyla))

library(ggplot2)
library(dplyr)
library(forcats)

# Remove NA/empty phyla, order phyla by median signal (for a structured x-axis)
plot_df <- summarized_data %>%
  filter(!is.na(GTDB_Tk_Phylum), GTDB_Tk_Phylum != "") %>%
  mutate(
    # (Optional) lock in your preferred depth order:
    # Depth_Category = fct_relevel(Depth_Category, "Shallow", "Mid", "Deep"),
    GTDB_Tk_Phylum = fct_drop(as.factor(GTDB_Tk_Phylum)),
    GTDB_Tk_Phylum = fct_reorder(GTDB_Tk_Phylum, pathway_present_log,
                                 .fun = median, .desc = TRUE)
  )

# Plot (Salinity on x, Phylum on y) — legend on the RIGHT
p <- ggplot(plot_df, aes(x = Salinity_Category, y = GTDB_Tk_Phylum)) +
  geom_tile(
    aes(fill = pathway_present_log),
    color = "white",
    size = 0.3,
    width = 0.98,
    height = 0.98
  ) +
  scale_x_discrete(expand = c(0, 0), drop = TRUE) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(
    option = "magma", direction = -1, na.value = "grey90",
    name = "Pathway presence (log)"
  ) +
  coord_fixed() +
  labs(title = NULL, x = "Salinity", y = "Phylum") +
  theme_minimal(base_size = 22, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 24),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.text = element_text(color = "black"),
    axis.ticks.length = unit(3, "mm"),
    plot.margin = margin(20, 20, 20, 20, "pt"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(face = "bold", size = 22),
    legend.text  = element_text(size = 20),
    legend.key.height = unit(6, "mm"),
    legend.key.width  = unit(10, "mm")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(120, "mm"),   # tall & thin for right side
      barwidth  = unit(6, "mm"),
      ticks = TRUE,
      label.theme = element_text(size = 18),
      title.theme = element_text(face = "bold", size = 22)
    )
  )

p

# High-resolution export — taller height helps many phyla on the y-axis
ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/03_general_biogeography/02_salinity_category_phylum.tiff",
  plot = p, width = 320, height = 600, units = "mm",
  dpi = 600, compression = "lzw"
)

# SVG export — add some width to make space for the right legend
ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/03_general_biogeography/01_temp_category_phylum.svg",
  plot = p, device = "svg",
  width = 420, height = 600, units = "mm", dpi = 1200
)# SVG export — add some width to make space for the right legend
ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/03_general_biogeography/02_salinity_category_phylum.svg",
  plot = p, device = "svg",
  width = 420, height = 600, units = "mm", dpi = 1200
)

#dev.off()

##########################################

# Read data
merged_data <- read.csv("merged_data.csv")

# Remove unknown salinity categories
merged_data <- merged_data %>% filter(Depth_Category != "Unknown")

# Ensure `occurence` is numeric
merged_data$occurence <- as.numeric(merged_data$occurence)

# Summarize data by Salinity_Category and GTDB_Tk_Phylum
summarized_data <- merged_data %>%
  group_by(Depth_Category, GTDB_Tk_Phylum) %>%
  summarise(pathway_present = sum(occurence, na.rm = TRUE), .groups = "drop") %>%
  filter(!Depth_Category %in% c("", NA))  # Remove empty or NA values

# Print results
print(summarized_data)

# Apply log10 transformation to pathway_present (to compress wide ranges)
summarized_data <- summarized_data %>%
  mutate(pathway_present_log = log10(pathway_present + 1))  # log(1 + value) to avoid log(0)

# Convert to wide format for clustering
library(tidyr)
data_wide <- summarized_data %>%
  pivot_wider(names_from = Salinity_Category, values_from = pathway_present_log, values_fill = 0)

# Perform hierarchical clustering on phyla
distance_matrix <- dist(data_wide[,-1])  # Compute distance matrix (excluding phylum names)
clustering <- hclust(distance_matrix, method = "ward.D2")  # Hierarchical clustering

# Order phyla based on clustering
ordered_phyla <- unique(data_wide$GTDB_Tk_Phylum[clustering$order])
summarized_data$GTDB_Tk_Phylum <- factor(summarized_data$GTDB_Tk_Phylum, 
                                         levels = unique(ordered_phyla))

library(ggplot2)
library(dplyr)
library(forcats)

# Remove NA/empty phyla, order phyla by median signal (for a structured x-axis)
plot_df <- summarized_data %>%
  filter(!is.na(GTDB_Tk_Phylum), GTDB_Tk_Phylum != "") %>%
  mutate(
    # (Optional) lock in your preferred depth order:
    # Depth_Category = fct_relevel(Depth_Category, "Shallow", "Mid", "Deep"),
    GTDB_Tk_Phylum = fct_drop(as.factor(GTDB_Tk_Phylum)),
    GTDB_Tk_Phylum = fct_reorder(GTDB_Tk_Phylum, pathway_present_log,
                                 .fun = median, .desc = TRUE)
  )

# Plot (Depth on x, Phylum on y) — legend on the RIGHT
p <- ggplot(plot_df, aes(x = Depth_Category, y = GTDB_Tk_Phylum)) +
  geom_tile(
    aes(fill = pathway_present_log),
    color = "white",
    size = 0.3,
    width = 0.98,
    height = 0.98
  ) +
  scale_x_discrete(expand = c(0, 0), drop = TRUE) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(
    option = "magma", direction = -1, na.value = "grey90",
    name = "Pathway presence (log)"
  ) +
  coord_fixed() +
  labs(title = NULL, x = "Depth", y = "Phylum") +
  theme_minimal(base_size = 22, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 24),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.text = element_text(color = "black"),
    axis.ticks.length = unit(3, "mm"),
    plot.margin = margin(20, 20, 20, 20, "pt"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(face = "bold", size = 22),
    legend.text  = element_text(size = 20),
    legend.key.height = unit(6, "mm"),
    legend.key.width  = unit(10, "mm")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(120, "mm"),  # tall & thin for right side
      barwidth  = unit(6, "mm"),
      ticks = TRUE,
      label.theme = element_text(size = 18),
      title.theme = element_text(face = "bold", size = 22)
    )
  )

p

# High-resolution export — taller height helps many phyla on the y-axis
ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/03_general_biogeography/03_depth_category_phylum_flipped.tiff",
  plot = p, width = 320, height = 600, units = "mm",
  dpi = 600, compression = "lzw"
)

ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/03_general_biogeography/03_depth_category_phylum_flipped.svg",
  plot = p, device = "svg",
  width = 420, height = 600, units = "mm", dpi = 1200
)


##########################################
## Oceans heatmap + hierarchical clustering
##########################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# ---- knobs you can tweak ------------------------------------------------------
GROUP_COL      <- "ocean_sea_name"  # <-- CHANGE if your column is different (e.g., "Ocean")
TOP_N          <- 50
dist_method    <- "euclidean"       # "euclidean", "manhattan", ...
hclust_method  <- "ward.D2"         # "ward.D2", "complete", "average", ...
cluster_rows   <- TRUE              # cluster phyla (y)
cluster_cols   <- TRUE              # cluster oceans (x)

# ---- read & prep --------------------------------------------------------------
merged_data <- read.csv("merged_data.csv")

merged_data <- merged_data %>%
  mutate(
    occurence = as.numeric(occurence),
    GTDB_Tk_Phylum = trimws(as.character(GTDB_Tk_Phylum)),
    !!GROUP_COL := trimws(as.character(.data[[GROUP_COL]]))
  ) %>%
  filter(!is.na(.data[[GROUP_COL]]),
         .data[[GROUP_COL]] != "",
         .data[[GROUP_COL]] != "Unknown")

# summarize occurrences per Ocean x Phylum
summarized_data <- merged_data %>%
  group_by(.data[[GROUP_COL]], GTDB_Tk_Phylum) %>%
  summarise(pathway_present = sum(occurence, na.rm = TRUE), .groups = "drop") %>%
  mutate(pathway_present_log = log10(pathway_present + 1)) %>%
  # unify the grouping column name for plotting/pivoting
  rename(Group = !!GROUP_COL)

# ---- limit to top N phyla (frequency first, then total) -----------------------
phylum_scores <- summarized_data %>%
  group_by(GTDB_Tk_Phylum) %>%
  summarise(
    freq_across_groups = sum(pathway_present > 0, na.rm = TRUE),
    total_occ          = sum(pathway_present, na.rm = TRUE),
    .groups = "drop"
  )

top_phyla <- phylum_scores %>%
  arrange(desc(freq_across_groups), desc(total_occ)) %>%
  slice_head(n = TOP_N) %>%
  pull(GTDB_Tk_Phylum)

summarized_top <- summarized_data %>%
  filter(GTDB_Tk_Phylum %in% top_phyla)

# ---- wide matrix for clustering (dup-safe) ------------------------------------
data_wide <- summarized_top %>%
  pivot_wider(
    names_from  = Group,
    values_from = pathway_present_log,
    values_fill = 0
  ) %>%
  group_by(GTDB_Tk_Phylum) %>%                                # collapse duplicate phylum rows
  summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)), .groups = "drop")

stopifnot(nrow(data_wide) > 0)

mat <- as.matrix(data_wide[, -1, drop = FALSE])               # rows = phyla, cols = oceans
rownames(mat) <- data_wide$GTDB_Tk_Phylum

# ---- clustering ---------------------------------------------------------------
if (cluster_rows && nrow(mat) > 1) {
  d_row  <- dist(mat, method = dist_method)
  hc_row <- hclust(d_row, method = hclust_method)
  ordered_phyla <- unique(rownames(mat)[hc_row$order])
} else {
  ordered_phyla <- rownames(mat)
}

if (cluster_cols && ncol(mat) > 1) {
  d_col  <- dist(t(mat), method = dist_method)
  hc_col <- hclust(d_col, method = hclust_method)
  ordered_groups <- unique(colnames(mat)[hc_col$order])
} else {
  ordered_groups <- colnames(mat)
}

# ---- long data with clustered orders ------------------------------------------
plot_df <- summarized_top %>%
  filter(GTDB_Tk_Phylum %in% ordered_phyla,
         Group          %in% ordered_groups) %>%
  mutate(
    GTDB_Tk_Phylum = factor(GTDB_Tk_Phylum, levels = ordered_phyla),
    Group          = factor(Group,          levels = ordered_groups)
  )

# ---- plot (Ocean on x, Phylum on y) -------------------------------------------
p <- ggplot(plot_df, aes(x = Group, y = GTDB_Tk_Phylum)) +
  geom_tile(
    aes(fill = pathway_present_log),
    color  = "white",
    size   = 0.3,
    width  = 0.98,
    height = 0.98
  ) +
  scale_x_discrete(expand = c(0, 0), drop = TRUE) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(
    option = "magma", direction = -1, na.value = "grey90",
    name = "Pathway presence (log)"
  ) +
  coord_fixed() +
  labs(title = NULL, x = "Ocean", y = "Phylum") +
  theme_minimal(base_size = 22, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 42),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 33),
    axis.text.y = element_text(size = 33),
    axis.text   = element_text(color = "black"),
    axis.ticks.length = unit(3, "mm"),
    plot.margin = margin(20, 20, 20, 20, "pt"),
    legend.position  = "right",
    legend.direction = "vertical",
    legend.title = element_text(face = "bold", size = 33),
    legend.text  = element_text(size = 33),
    legend.key.height = unit(6, "mm"),
    legend.key.width  = unit(10, "mm")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(120, "mm"),
      barwidth  = unit(6, "mm"),
      ticks = TRUE,
      label.theme = element_text(size = 18),
      title.theme = element_text(face = "bold", size = 22)
    )
  )

p

# ---- export (sizes that won’t OOM) -------------------------------------------
ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/03_general_biogeography/XX_ocean_phylum_top50.tiff",
  plot = p, width = 720, height = 900, units = "mm",
  dpi = 900, compression = "lzw"
)

ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/03_general_biogeography/XX_ocean_phylum_top50.svg",
  plot = p, device = "svg",
  width = 720, height = 900, units = "mm", dpi = 1200
)



############################################
############################################

##########################################

library(dplyr)
library(tidyr)
library(ggplot2)

# --- clustering knobs ----------------------------------------------------------
dist_method   <- "euclidean"   # "euclidean", "manhattan", ...
hclust_method <- "ward.D2"     # "ward.D2", "complete", "average", ...
cluster_rows  <- TRUE
cluster_cols  <- TRUE

# --- Read & prep --------------------------------------------------------------
merged_data <- read.csv("merged_data.csv")

# Keep known seas only; ensure numeric & tidy labels
merged_data <- merged_data %>%
  filter(IHO_Sea != "Unknown", IHO_Sea != "", !is.na(IHO_Sea)) %>%
  mutate(
    occurence       = as.numeric(occurence),
    GTDB_Tk_Phylum  = trimws(as.character(GTDB_Tk_Phylum)),
    IHO_Sea         = trimws(as.character(IHO_Sea))
  )

# Summarize occurrences per Sea x Phylum
summarized_data <- merged_data %>%
  group_by(IHO_Sea, GTDB_Tk_Phylum) %>%
  summarise(pathway_present = sum(occurence, na.rm = TRUE), .groups = "drop") %>%
  mutate(pathway_present_log = log10(pathway_present + 1))

# --- Limit to TOP 50 phyla (frequency first, then total) ---------------------
N <- 50
phylum_scores <- summarized_data %>%
  group_by(GTDB_Tk_Phylum) %>%
  summarise(
    freq_across_seas = sum(pathway_present > 0, na.rm = TRUE),  # number of seas with ≥1 record
    total_occ        = sum(pathway_present, na.rm = TRUE),
    .groups = "drop"
  )

top_phyla <- phylum_scores %>%
  arrange(desc(freq_across_seas), desc(total_occ)) %>%
  slice_head(n = N) %>%
  pull(GTDB_Tk_Phylum)

summarized_data_top <- summarized_data %>%
  filter(GTDB_Tk_Phylum %in% top_phyla)

# --- Build wide matrix for clustering (dup-safe) ------------------------------
data_wide <- summarized_data_top %>%
  pivot_wider(
    names_from  = IHO_Sea,
    values_from = pathway_present_log,
    values_fill = 0
  ) %>%
  group_by(GTDB_Tk_Phylum) %>%   # collapse any duplicate phylum rows
  summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)), .groups = "drop")

stopifnot(nrow(data_wide) > 0)

mat <- as.matrix(data_wide[, -1, drop = FALSE])   # rows = phyla, cols = seas
rownames(mat) <- data_wide$GTDB_Tk_Phylum

# --- Row clustering (phyla) ---------------------------------------------------
if (cluster_rows && nrow(mat) > 1) {
  d_row  <- dist(mat, method = dist_method)
  hc_row <- hclust(d_row, method = hclust_method)
  ordered_phyla <- unique(rownames(mat)[hc_row$order])
} else {
  ordered_phyla <- rownames(mat)
}

# --- Column clustering (seas) -------------------------------------------------
if (cluster_cols && ncol(mat) > 1) {
  d_col  <- dist(t(mat), method = dist_method)
  hc_col <- hclust(d_col, method = hclust_method)
  ordered_seas <- unique(colnames(mat)[hc_col$order])
} else {
  ordered_seas <- colnames(mat)
}

# --- Apply clustered orders to long data --------------------------------------
summarized_data_top <- summarized_data_top %>%
  filter(GTDB_Tk_Phylum %in% ordered_phyla,
         IHO_Sea        %in% ordered_seas) %>%
  mutate(
    GTDB_Tk_Phylum = factor(GTDB_Tk_Phylum, levels = ordered_phyla),
    IHO_Sea        = factor(IHO_Sea,        levels = ordered_seas)
  )

# --- Plot (Sea on x, Phylum on y) — legend on RIGHT --------------------------
p <- ggplot(summarized_data_top, aes(x = IHO_Sea, y = GTDB_Tk_Phylum)) +
  geom_tile(
    aes(fill = pathway_present_log),
    color  = "white",
    size   = 0.3,
    width  = 0.98,
    height = 0.98
  ) +
  scale_x_discrete(expand = c(0, 0), drop = TRUE) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(
    option = "magma", direction = -1, na.value = "grey90",
    name = "Pathway presence (log)"
  ) +
  coord_fixed() +
  labs(title = NULL, x = "Sea", y = "Phylum") +
  theme_minimal(base_size = 22, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 42),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 33),
    axis.text.y = element_text(size = 33),
    axis.text   = element_text(color = "black"),
    axis.ticks.length = unit(3, "mm"),
    plot.margin = margin(20, 20, 20, 20, "pt"),
    legend.position  = "right",
    legend.direction = "vertical",
    legend.title = element_text(face = "bold", size = 33),
    legend.text  = element_text(size = 38),
    legend.key.height = unit(12, "mm"),
    legend.key.width  = unit(18, "mm")
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(140, "mm"),  # tall & thin for right side
      barwidth  = unit(12, "mm"),
      ticks = TRUE,
      label.theme = element_text(size = 22),
      title.theme = element_text(face = "bold", size = 26)
    )
  )

p

# --- Export ------------------------------------------------------------------
# Consider using vector (SVG/PDF) for large layouts; raster TIFF can be huge.
ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/03_general_biogeography/05_IHO_sea_phylum_top50.tiff",
  plot = p, width = 820, height = 1200, units = "mm",
  dpi = 600, compression = "lzw"
)

ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/03_general_biogeography/05_IHO_sea_phylum_top50.svg",
  plot = p, device = "svg",
  width = 920, height = 1200, units = "mm", dpi = 1200
)

#################
#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/01_pathway_taxa_region")


# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)


# Replace the file paths with your actual input file paths
pathway_data <- read.csv("CF_pathways_merged.csv")
coord_attr_01 <- read.csv("EEZ_output_with_full_names.csv")
coord_attr_02 <- read.csv("output_with_ocean_sea_names.csv")
cat_values_mtd <- read.csv("MAGs_Mtd_cat.csv")


# Categorize depth into zones using consistent logic
categories <- cat_values_mtd %>%
  mutate(
    Depth_Category = case_when(
      sample_depth <= 50 ~ "Shallow Zone",
      sample_depth > 50 & sample_depth <= 200 ~ "Epipelagic",
      sample_depth > 200 & sample_depth <= 1000 ~ "Mesopelagic",
      sample_depth > 1000 & sample_depth <= 4000 ~ "Bathypelagic",
      sample_depth > 4000 & sample_depth <= 6000 ~ "Abyssopelagic",
      sample_depth > 6000 ~ "Hadalpelagic",
      TRUE ~ "Unknown"  # Catch any unexpected cases
    )
  )


# Categorize depth into zones using consistent logic
categories <- categories %>%
  mutate(
    Salinity_Category = case_when(
      sample_salinity <= 10 ~ "Low Salinity",
      sample_salinity > 10 & sample_salinity <= 25 ~ "Moderate Salinity",
      sample_salinity > 25 & sample_salinity <= 35 ~ "High Salinity",
      sample_salinity > 35 ~ "Very High Salinity",
      TRUE ~ "Unknown"  # Catch any unexpected cases
    )
  )

categories <- categories %>%
  mutate(
    temp_Category = case_when(
      sample_temp <= 5 ~ "Cold Temperature",
      sample_temp > 5 & sample_temp <= 15 ~ "Moderate Temperature",
      sample_temp > 15 & sample_temp <= 25 ~ "Warm Temperature",
      sample_temp > 25 ~ "High Temperature",
      TRUE ~ "Unknown"  # Catch any unexpected cases
    )
  )

#write.csv(categories, file = "temp_depth_sal_cat.csv")

# Merge the datasets
biogeographic_data <- left_join(coord_attr_01, coord_attr_02, by = "Genome_ID")
biogeographic_data <- left_join(biogeographic_data, categories, by = "Genome_ID")


merged_data <- read.csv("merged_data.csv")

# Assuming your data is in a dataframe called `merged_data`
# Filter to keep only rows where CBB == 1 (presence)
filtered_data <- merged_data %>%
  filter(CBB == 1)

############################################################

library(dplyr)
library(forcats)
library(ggplot2)
library(grid)  # for unit()

# ---- TEXT KNOBS ---------------------------------------------------------------
BASE    <- 27  # general text (tick labels, legend text)  [was 16]
TITLE   <- 30  # axis/legend titles                       [was 16]
X_ANGLE <- 90  # x-axis label angle (set 90 for vertical)

# ---- Summarize ---------------------------------------------------------------
ihosea_df <- filtered_data %>%
  filter(!is.na(IHO_Sea), IHO_Sea != "",
         !is.na(GTDB_Tk_Phylum), GTDB_Tk_Phylum != "") %>%
  mutate(CBB = as.numeric(CBB)) %>%
  count(IHO_Sea, GTDB_Tk_Phylum, wt = CBB, name = "pathway_present") %>%
  mutate(pathway_present_log = log10(pathway_present + 1))

# Order phyla (columns) by median signal; order seas (rows) by total signal
phylum_order <- ihosea_df %>%
  group_by(GTDB_Tk_Phylum) %>%
  summarise(med = median(pathway_present_log, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>% pull(GTDB_Tk_Phylum)

sea_order <- ihosea_df %>%
  group_by(IHO_Sea) %>%
  summarise(total = sum(pathway_present, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(IHO_Sea)

# Shared color limits (clip at 99th percentile but keep non-zero range)
all_vals <- ihosea_df$pathway_present_log
lo <- 0
hi <- unname(quantile(all_vals, 0.99, na.rm = TRUE))
if (!is.finite(hi) || hi <= lo) hi <- max(all_vals, na.rm = TRUE)

# ---- Plot --------------------------------------------------------------------
p <- ihosea_df %>%
  mutate(
    GTDB_Tk_Phylum = factor(GTDB_Tk_Phylum, levels = phylum_order),
    IHO_Sea        = factor(IHO_Sea,        levels = sea_order)
  ) %>%
  ggplot(aes(x = GTDB_Tk_Phylum, y = IHO_Sea)) +
  geom_tile(aes(fill = pathway_present_log), color = "white", size = 0.25) +
  scale_x_discrete(expand = c(0, 0), drop = TRUE) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(
    option = "magma", direction = -1, na.value = "grey90",
    name = expression("CBB presence ("*log[10]*"(count+1))"),
    limits = c(lo, hi), oob = scales::squish,
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      direction = "vertical",
      barheight = unit(120, "mm"),   # taller, easier to read
      barwidth  = unit(6,   "mm"),
      ticks = TRUE
    )
  ) +
  coord_fixed() +
  labs(
    title = "",
    x = "Phylum", y = "IHO sea"
  ) +
  theme_minimal(base_size = BASE, base_family = "Helvetica") +
  theme(
    panel.grid    = element_blank(),
    plot.title    = element_text(face = "bold", size = TITLE),
    axis.title    = element_text(face = "bold", size = TITLE),
    axis.text     = element_text(color = "black", size = BASE),
    axis.text.x   = element_text(angle = X_ANGLE, hjust = 1, vjust = 1, size = BASE),
    legend.position = "right",
    legend.title  = element_text(face = "bold", size = TITLE),
    legend.text   = element_text(size = BASE),
    legend.key.height = unit(8, "mm"),
    legend.key.width  = unit(6, "mm"),
    plot.margin   = margin(8, 10, 8, 8, "pt")
  )

p

# ---- High-res export (dimensions unchanged) ----------------------------------
ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/01_pathway_taxa_region/01_CBB_IHOsea_phylum.svg",
  plot = p, width = 24, height = 26, units = "in"
)

ggsave("01_CBB_IHOsea_phylum.png",  p, width = 24, height = 26, units = "in", dpi = 600, bg = "white")
ggsave("01_CBB_IHOsea_phylum.tiff", p, width = 24, height = 26, units = "in", dpi = 600, compression = "lzw", bg = "white")

####################################################################
####################################################################

library(dplyr)
library(forcats)
library(ggplot2)
library(grid)  # for unit()

# ---- TEXT KNOBS ---------------------------------------------------------------
BASE  <- 25  # general text (tick labels, legend text)
TITLE <- 27  # axis/legend titles
X_ANGLE <- 90  # keep your 60° x labels (set to 90 if you prefer)

# ---- Summarize ---------------------------------------------------------------
ocean_df <- filtered_data %>%
  filter(!is.na(ocean_sea_name), ocean_sea_name != "",
         !is.na(GTDB_Tk_Phylum), GTDB_Tk_Phylum != "") %>%
  mutate(CBB = as.numeric(CBB)) %>%
  count(ocean_sea_name, GTDB_Tk_Phylum, wt = CBB, name = "pathway_present") %>%
  mutate(pathway_present_log = log10(pathway_present + 1))

# Order phyla (columns) by median signal; order oceans/seas (rows) by total signal
phylum_order <- ocean_df %>%
  group_by(GTDB_Tk_Phylum) %>%
  summarise(med = median(pathway_present_log, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>% pull(GTDB_Tk_Phylum)

ocean_order <- ocean_df %>%
  group_by(ocean_sea_name) %>%
  summarise(total = sum(pathway_present, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(ocean_sea_name)

# Robust color limits
all_vals <- ocean_df$pathway_present_log
lo <- 0
hi <- unname(quantile(all_vals, 0.99, na.rm = TRUE))
if (!is.finite(hi) || hi <= lo) hi <- max(all_vals, na.rm = TRUE)

# ---- Plot --------------------------------------------------------------------
p <- ocean_df %>%
  mutate(
    GTDB_Tk_Phylum = factor(GTDB_Tk_Phylum, levels = phylum_order),
    ocean_sea_name = factor(ocean_sea_name, levels = ocean_order)
  ) %>%
  ggplot(aes(x = GTDB_Tk_Phylum, y = ocean_sea_name)) +
  geom_tile(aes(fill = pathway_present_log), color = "white", size = 0.25) +
  scale_x_discrete(expand = c(0, 0), drop = TRUE) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_viridis_c(
    option = "magma", direction = -1, na.value = "grey90",
    name = expression("CBB presence ("*log[10]*"(count+1))"),
    limits = c(lo, hi), oob = scales::squish,
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      direction = "vertical",
      barheight = unit(120, "mm"),   # taller, easier to read
      barwidth  = unit(6,   "mm"),
      ticks = TRUE
    )
  ) +
  coord_fixed() +
  labs(
    title = "",
    x = "Phylum", y = "Ocean/Sea name"
  ) +
  theme_minimal(base_size = BASE, base_family = "Helvetica") +
  theme(
    panel.grid    = element_blank(),
    plot.title    = element_text(face = "bold", size = TITLE),
    axis.title    = element_text(face = "bold", size = TITLE),
    axis.text     = element_text(color = "black", size = BASE),
    axis.text.x   = element_text(angle = X_ANGLE, hjust = 1, vjust = 1, size = BASE),
    legend.position = "right",
    legend.title  = element_text(face = "bold", size = TITLE),
    legend.text   = element_text(size = BASE),
    plot.margin   = margin(8, 10, 8, 8, "pt")
  )

p

# ---- High-res export (dimensions unchanged) ----------------------------------
ggsave(
  filename = "/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/01_pathway_taxa_region/03_CBB_OceanSeaName_phylum.svg",
  plot = p, width = 22, height = 10, units = "in"
)
ggsave("03_CBB_OceanSeaName_phylum.png",  p, width = 22, height = 10, units = "in", dpi = 600, bg = "white")
ggsave("03_CBB_OceanSeaName_phylum.tiff", p, width = 22, height = 10, units = "in", dpi = 600, compression = "lzw", bg = "white")


# =========================
# One column (Depth / Salinity / Temperature), shared legend at RIGHT
# =========================
library(dplyr)
library(forcats)
library(ggplot2)
library(patchwork)

# If you need different subsets per panel, set them here:
filtered_data_depth    <- filtered_data
filtered_data_salinity <- filtered_data
filtered_data_temp     <- filtered_data

# ---- Helper to summarize one gradient (uses dplyr::count to avoid masking) ---
sum_by <- function(df, cat_col) {
  df %>%
    filter(
      !is.na({{cat_col}}), {{cat_col}} != "",
      !is.na(GTDB_Tk_Phylum), GTDB_Tk_Phylum != ""
    ) %>%
    mutate(CBB = as.numeric(CBB)) %>%  # ensure numeric for wt=
    count({{cat_col}}, GTDB_Tk_Phylum, wt = CBB, name = "pathway_present") %>%
    mutate(
      pathway_present_log = log10(pathway_present + 1),
      Category = as.character({{cat_col}})
    ) %>%
    select(Category, GTDB_Tk_Phylum, pathway_present_log)
}

# ---- Build datasets -----------------------------------------------------------
depth_df <- sum_by(filtered_data_depth,    Depth_Category)
sal_df   <- sum_by(filtered_data_salinity, Salinity_Category)
temp_df  <- sum_by(filtered_data_temp,     temp_Category)

# ---- Global phylum order (by median across ALL panels) ------------------------
phylum_order <- bind_rows(depth_df, sal_df, temp_df) %>%
  group_by(GTDB_Tk_Phylum) %>%
  summarise(med = median(pathway_present_log, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(GTDB_Tk_Phylum)

# ---- Shared color limits (robust to outliers) --------------------------------
all_vals <- bind_rows(depth_df, sal_df, temp_df) %>% pull(pathway_present_log)
lo <- 0
hi <- unname(quantile(all_vals, 0.99, na.rm = TRUE))
if (!is.finite(hi)) hi <- max(all_vals, na.rm = TRUE)

# ---- SIZE KNOBS ---------------------------------------------------------------
BASE  <- 20  # general text
TITLE <- 22  # axis/legend titles
TAG   <- 22  # (a)(b)(c) tags

# ---- Heatmap factory (VERTICAL colorbar for right-side legend) ----------------
make_heat <- function(df, y_lab) {
  df %>%
    mutate(GTDB_Tk_Phylum = factor(GTDB_Tk_Phylum, levels = phylum_order)) %>%
    ggplot(aes(x = GTDB_Tk_Phylum, y = Category)) +
    geom_tile(aes(fill = pathway_present_log), color = "white", size = 0.2) +
    scale_x_discrete(expand = c(0, 0), drop = TRUE) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_viridis_c(
      option = "magma", direction = -1, na.value = "grey90",
      name = expression("CBB presence ("*log[10]*"(count+1))"),
      limits = c(lo, hi), oob = scales::squish
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top", title.hjust = 0.5,
        direction = "vertical",
        barheight = unit(90, "mm"),
        barwidth  = unit(6,  "mm"),
        ticks = TRUE
      )
    ) +
    coord_fixed() +
    labs(x = "Phylum", y = y_lab) +
    theme_minimal(base_size = BASE, base_family = "Helvetica") +
    theme(
      panel.grid  = element_blank(),
      axis.title  = element_text(face = "bold", size = TITLE),
      axis.text   = element_text(color = "black", size = BASE),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE),
      legend.title = element_text(face = "bold", size = TITLE),
      legend.text  = element_text(size = BASE),
      plot.margin  = margin(6, 10, 14, 6, "pt")
    )
}

p_depth <- make_heat(depth_df, "Depth")
p_sal   <- make_heat(sal_df,   "Salinity")
p_temp  <- make_heat(temp_df,  "Temperature")

# ---- Combine: ONE COLUMN, THREE ROWS; shared legend at RIGHT -----------------
combined <- (p_depth / p_sal / p_temp) +
  plot_layout(ncol = 1, guides = "collect", heights = c(1, 1, 1)) +
  plot_annotation(
    tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
    theme = theme(plot.tag = element_text(size = TAG, face = "bold"))
  ) &
  theme(
    legend.position  = "right",
    legend.direction = "vertical",
    legend.box       = "vertical"
  )

# Preview
combined

# ---- High-resolution exports --------------------------------------------------
ggsave("CBB_gradients_triptych.png",  combined, width = 14, height = 16, units = "in", dpi = 600, bg = "white")
ggsave("CBB_gradients_triptych.tiff", combined, width = 14, height = 16, units = "in", dpi = 600, compression = "lzw", bg = "white")
ggsave("CBB_gradients_triptych.svg",  combined, width = 14, height = 16, units = "in")

#################################################
