#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/04_library_figures")

library(dplyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(ggmap)
library(sp)
library(maps)
library(cowplot)

################################################
### Figure 2c


#loading the libraries edited data dataframe
libs <- read.csv(file = "/Thesis/work_package_3/Data_analysis/01_figures/04_library_figures/libs_data_edited.csv", header = TRUE, stringsAsFactors = FALSE)

# Function to create pie chart for a given category
create_pie_chart <- function(data, category_column, title) {
  # Generate frequency table
  freq_table <- data %>%
    count(!!sym(category_column)) %>%
    mutate(
      percentage = n / sum(n) * 100,
      labels = paste0(if_else(is.na(!!sym(category_column)), "NA", !!sym(category_column)), " (", round(percentage, 2), "%)"),
      Category = if_else(is.na(!!sym(category_column)), "NA", !!sym(category_column))
    ) %>%
    arrange(desc(percentage)) # Ensure order matches pie slices
  
  # Define a better color palette
  palette_colors <- brewer.pal(n = min(9, nrow(freq_table)), "Set2") # Use Set2 palette
  assigned_colors <- setNames(palette_colors[1:nrow(freq_table)], freq_table$labels)
  
  # Create pie chart
  ggplot(freq_table, aes(x = "", y = percentage, fill = labels)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(
      values = assigned_colors # Ensure colors match the legend labels
    ) +
    labs(
      title = title,
      fill = category_column,
      y = NULL,
      x = NULL
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
}


# Generate and plot for Salinity Categories
salinity_pie <- create_pie_chart(libs, "Salinity.Categories", "Salinity Categories") +
  theme(text = element_text(size = 20))  # Adjust the size as needed
print(salinity_pie)

# Generate and plot for Temperature Categories
temperature_pie <- create_pie_chart(libs, "Temperature.Categories", "Temperature Categories") +
  theme(text = element_text(size = 20))  # Adjust the size as needed
print(temperature_pie)

# Generate and plot for Depth Categories
depth_pie <- create_pie_chart(libs, "Depth.Categories", "Depth Categories") +
  theme(text = element_text(size = 20))  # Adjust the size as needed
print(depth_pie)

combined_pie <- grid.arrange(temperature_pie, salinity_pie, depth_pie, ncol = 3)

svg("/Thesis/work_package_3/Data_analysis/01_figures/04_library_figures/02_pie_charts.svg", width = 16, height = 6)

# Step 4: Plot the pie charts
grid.arrange(temperature_pie, salinity_pie, depth_pie, ncol = 3)

dev.off()

################################################
### Figure 2

# Function to create bar plots for frequencies with percentages
create_bar_plot <- function(data, column_name, title, top_n = 10) {
  # Frequency table with percentages
  freq_table <- data %>%
    count(!!sym(column_name)) %>%
    arrange(desc(n)) %>%
    mutate(
      percentage = (n / sum(n)) * 100, # Calculate percentage
      !!sym(column_name) := if_else(
        !!sym(column_name) == "" | is.na(!!sym(column_name)),
        "NA",
        !!sym(column_name)
      ),
      !!sym(column_name) := fct_reorder(!!sym(column_name), n) # For ordered bars
    ) %>%
    slice_max(n, n = top_n) # Get top n entries
  
  # Bar plot with percentage labels
  ggplot(freq_table, aes(x = !!sym(column_name), y = n, fill = !!sym(column_name))) +
    geom_col(color = "black", width = 0.7) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), 
              hjust = -0.2, size = 6) +  # Adjusted hjust for better placement
    coord_flip() + # Horizontal bars
    scale_fill_brewer(palette = "Set3") +
    expand_limits(y = max(freq_table$n) * 1.1) +  # Expand y-axis to avoid text cutoff
    labs(
      title = title,
      x = ".",
      y = "Frequency",
      fill = column_name
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title.x = element_text(size = 19),
      axis.title.y = element_text(size = 20),
      legend.position = "none" # Remove legend for cleaner display
    )
}

################
### Figure 2b

svg("/Thesis/work_package_3/Data_analysis/01_figures/04_library_figures/03_ocean_sea.svg", width = 14, height = 8)

# Generate and plot for ocean_sea_name
ocean_sea_plot <- create_bar_plot(libs, "Ocean_and_Sea", "")
print(ocean_sea_plot)

dev.off()

################
### Figure 2d

svg("/Thesis/work_package_3/Data_analysis/01_figures/04_library_figures/06_IHO_sea.svg", width = 12, height = 8)
# Generate and plot for IHO_Sea
iho_sea_plot <- create_bar_plot(libs, "IHO_Sea", "", top_n = 10)
print(iho_sea_plot)

dev.off()

################
### Figure S19

svg("/Thesis/work_package_3/Data_analysis/01_figures/04_library_figures/07_Marine_Region.svg", width = 18, height = 8)
# Generate and plot for IHO_Sea
marine_region_plot <- create_bar_plot(libs, "Marine_Region", "", top_n = 12)
print(marine_region_plot)

dev.off()

################
### Figure S20

svg("/Thesis/work_package_3/Data_analysis/01_figures/04_library_figures/05_territory.svg", width = 18, height = 8)
# Generate and plot for TERRITORY1
territory_plot <- create_bar_plot(libs, "Territory", "", top_n = 12)
print(territory_plot)

dev.off()


################
### Figure 2a



#loading the final marmags dataframe
marmags_dataframe <- read.csv(file = "01_MarMAGs_Dataframe_03062024_.csv", header = TRUE, stringsAsFactors = FALSE)

fig1a <- marmags_dataframe[, c("Genome_source", "sample_latitude", "sample_longitude")]


# Remove rows with missing latitude or longitude
fig1a <- fig1a[!is.na(fig1a$sample_latitude) & !is.na(fig1a$sample_longitude), ]

# Define the base world map
mapworld <- borders("world", colour = "gray85", fill = "gray80")

svg("/Thesis/work_package_3/Data_analysis/01_figures/03_biogeography_analysis/04_world_map.svg", width = 12, height = 8)

# Plot the data
sample_map <- ggplot(fig1a) +
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

