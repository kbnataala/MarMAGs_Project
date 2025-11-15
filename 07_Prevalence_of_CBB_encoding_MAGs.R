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

###############################################
##### Plotting prevalence figures
###############################################


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
## Fig 5a

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
## Fig 5b

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
## Fig 5c

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
## Fig 5d

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
## Fig S12

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
## Fig S13

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
## Fig S11

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



