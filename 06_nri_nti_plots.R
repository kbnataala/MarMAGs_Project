##########################
#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI")


# Load necessary libraries
library(ggplot2)
library(dplyr)

#loading the data
data <- read.csv(file = "scatter_plot_all.csv", header = TRUE, stringsAsFactors = FALSE)


# # Ensure Metrics are factor variables with specific shapes
# data$Metrics <- factor(data$Metrics, levels = c("NRI", "NTI"))
# 
# svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/01_all_attributes_nri_nti.svg", width = 12, height = 8)
# # Scatter dot plot
# ggplot(data, aes(x = Attributes, y = NRI_and_NTI, color = Classes, shape = Metrics)) +
#   geom_point(size = 4, alpha = 0.8, position = position_jitter(width = 0.2)) +  # Jitter for better separation
#   scale_shape_manual(values = c(16, 17)) +  # Circle for NRI (16), Triangle for NTI (17)
#   theme_minimal() +
#   geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +  # Dashed reference lines
#   labs(title = "Dot Plot of NRI and NTI Grouped by Attributes",
#        x = "Attributes",
#        y = "NRI and NTI Values",
#        shape = "Metrics",
#        color = "Classes") +
#   guides(color = "none") +  # Suppress legend for colors
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 20),
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
#     axis.text.y = element_text(size = 25),
#     axis.title.x = element_text(size = 25),
#     axis.title.y = element_text(size = 25),
#     legend.key.width = unit(1, "cm"),
#     legend.key.height = unit(1, "cm"),
#     legend.position = "top",
#     panel.grid.major.x = element_blank(),
#     legend.title = element_text(size = 25),
#     legend.text = element_text(size = 25)
#   )
# 
# dev.off()


# Filter data for NRI only
data_NRI <- data %>% filter(Metrics == "NRI")

# Save NRI plot
svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/NRI_plot.svg", width = 12, height = 8)

# Scatter dot plot for NRI
ggplot(data_NRI, aes(x = Attributes, y = NRI_and_NTI, color = Classes)) +
  geom_point(size = 4, alpha = 0.8, position = position_jitter(width = 0.2)) +  # Jitter for better separation
  theme_minimal() +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +  # Dashed reference lines
  labs(title = "Dot Plot of NRI Grouped by Attributes",
       x = "Attributes",
       y = "NRI Values",
       color = "Classes") +
  guides(color = "none") +  # Suppress legend for colors
  ylim(c(-10, 65)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  )

dev.off()

#######################

# Filter data for NTI only
data_NTI <- data %>% filter(Metrics == "NTI")

# Save NTI plot
svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/NTI_plot.svg", width = 12, height = 8)

# Scatter dot plot for NTI with triangles
ggplot(data_NTI, aes(x = Attributes, y = NRI_and_NTI, color = Classes)) +
  geom_point(size = 4, alpha = 0.8, position = position_jitter(width = 0.2), shape = 17) +  # Triangle shape
  theme_minimal() +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +  # Dashed reference lines
  labs(title = "Dot Plot of NTI Grouped by Attributes",
       x = "Attributes",
       y = "NTI Values",
       color = "Classes") +
  ylim(c(-10, 65)) +
  guides(color = "none") +  # Suppress legend for colors
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    plot.tag = element_text(size = 30, face = "bold")  # Increase tag size
  )

dev.off()



############################################################################
library(ggplot2)
library(gridExtra)

# Function to standardize bar plots
standardize_bar_chart <- function(plot) {
  plot +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
      axis.text.y = element_text(size = 25),
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.position = "top",
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 25),
      plot.tag = element_text(size = 30, face = "bold")  # Increase tag size
    ) +
    scale_y_continuous(limits = c(-4, 30))  # Set a common Y-axis limit for consistency
}

#################################################
### Depth categories

data <- read.csv(file = "depth_categories.csv", header = TRUE, stringsAsFactors = FALSE)
data$Depth <- factor(data$Depth, levels = c("Shallow Zone", "Epipelagic", "Bathypelagic", "Mesopelagic", "Abyssopelagic"))

plot1 <- ggplot(data, aes(x = Depth, y = NRI_and_NTI, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Standardized bar width
  geom_errorbar(aes(ymin = NRI_and_NTI - Std_dev, ymax = NRI_and_NTI + Std_dev), 
                position = position_dodge(width = 0.6), width = 0.2) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +
  scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +
  coord_flip() +
  labs(title = "", x = "Depth", y = "", fill = "Metrics")

plot1 <- standardize_bar_chart(plot1)

svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/03_depth_nri_nti.svg", width = 12, height = 8)
print(plot1)
dev.off()

####################################################
##### Salinity categories

data <- read.csv(file = "salinity_categories.csv", header = TRUE, stringsAsFactors = FALSE)
# Ensure Salinity is a factor with proper ordering
data$Salinity <- factor(data$Salinity, levels = c("Very High", "High", "Moderate", "Low"))

# Create the horizontal grouped divergent bar plot
plot2 <- ggplot(data, aes(x = Salinity, y = NRI_and_NTI, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Grouped bars
  geom_errorbar(aes(ymin = NRI_and_NTI - Std_dev, ymax = NRI_and_NTI + Std_dev), 
                position = position_dodge(width = 0.6), width = 0.2) +  # Error bars
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +  # Dashed reference lines
  scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +  # Custom colors
  coord_flip() +  # Flip the plot to horizontal
  theme_minimal() +
  labs(title = "", x = "Salinity", y = "", fill = "Metrics")

plot2 <- standardize_bar_chart(plot2)  

svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/09_Salinity_nri_nti.svg", width = 8, height = 8)
print(plot2)
dev.off() 

#######################################################
### Salinity categories

data <- read.csv(file = "temperature_categories.csv", header = TRUE, stringsAsFactors = FALSE)
# Ensure Salinity is a factor with proper ordering
data$Temperature_Categories <- factor(data$Temperature_Categories, levels = c( "High", "Warm", "Moderate", "Cold"))

# Create the horizontal grouped divergent bar plot
plot3 <- ggplot(data, aes(x = Temperature_Categories, y = NRI_and_NTI, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Grouped bars
  geom_errorbar(aes(ymin = NRI_and_NTI - Std_dev, ymax = NRI_and_NTI + Std_dev), 
                position = position_dodge(width = 0.6), width = 0.2) +  # Error bars
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +  # Dashed reference lines
  scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +  # Custom colors
  coord_flip() +  # Flip the plot to horizontal
  theme_minimal() +
  labs(title = "", x = "Temperature", y = "NRI and NTI", fill = "Metrics")

plot3 <- standardize_bar_chart(plot3)

svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/11_temperature_nri_nti.svg", width = 8, height = 8)
print(plot3)  
dev.off() 

################################################################
################################################################

# Function to standardize bar plots
standardize_bar_chart <- function(plot) {
  plot +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
      axis.text.y = element_text(size = 25),
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.position = "top",
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 25),
      plot.tag = element_text(size = 30, face = "bold")  # Increase tag size
    ) +
    scale_y_continuous(limits = c(-4, 60))  # Set a common Y-axis limit for consistency
}

####################################

data <- read.csv(file = "phylum.csv", header = TRUE, stringsAsFactors = FALSE)
# Ensure Salinity is a factor with proper ordering
data$Phylum <- factor(data$Phylum, levels = c("Verrucomicrobiota", "Myxococcota", "Spirochaetota", "Patescibacteria",
                                              "Calditrichota", "Campylobacterota", "Chloroflexota", "Desulfobacterota",
                                              "Bacteroidota", "Actinomycetota", "SAR324", "Cyanobacteriota", "Pseudomonadota"))

# Create the horizontal grouped divergent bar plot
plot4 <- ggplot(data, aes(x = Phylum, y = NRI_and_NTI, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Grouped bars
  geom_errorbar(aes(ymin = NRI_and_NTI - Std_dev, ymax = NRI_and_NTI + Std_dev), 
                position = position_dodge(width = 0.6), width = 0.2) +  # Error bars
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +  # Dashed reference lines
  scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +  # Custom colors
  coord_flip() +  # Flip the plot to horizontal
  theme_minimal() +
  labs(title = "", x = "Phylum", y = "", fill = "Metrics")

plot4 <- standardize_bar_chart(plot4)

svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/08_Phylum_nri_nti.svg", width = 8, height = 10)  
print(plot4)
dev.off() 

####################################

# Function to standardize bar plots
standardize_bar_chart <- function(plot) {
  plot +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
      axis.text.y = element_text(size = 25),
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.position = "top",
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 25),
      plot.tag = element_text(size = 30, face = "bold")  # Increase tag size
    ) +
    scale_y_continuous(limits = c(-4, 30))  # Set a common Y-axis limit for consistency
}

####################################

data <- read.csv(file = "ocean.csv", header = TRUE, stringsAsFactors = FALSE)

# Ensure Salinity is a factor with proper ordering
data$Ocean <- factor(data$Ocean, levels = c("Southern Ocean", "S-China and E-A Seas",
                                            "Arctic Ocean", "Mediterranean Region", "Baltic Sea", 
                                            "South Atlantic Ocean", "South Pacific Ocean", "Indian Ocean",
                                             "North Atlantic Ocean", "North Pacific Ocean"))

# Create the horizontal grouped divergent bar plot
plot5 <- ggplot(data, aes(x = Ocean, y = NRI_and_NTI, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Grouped bars
  geom_errorbar(aes(ymin = NRI_and_NTI - Std_dev, ymax = NRI_and_NTI + Std_dev), 
                position = position_dodge(width = 0.6), width = 0.2) +  # Error bars
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +  # Dashed reference lines
  scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +  # Custom colors
  coord_flip() +  # Flip the plot to horizontal
  theme_minimal() +
  labs(title = "", x = "Ocean", y = "NRI and NTI", fill = "Metrics")

plot5 <- standardize_bar_chart(plot5)

svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/06_Ocean_nri_nti.svg", width = 14, height = 9)  
print(plot5)
dev.off() 

##################################################

library(ggplot2)
library(patchwork)

# Simulating 3 other smaller plots (replace with real plots)
plot1 <- plot1 + labs(tag = "(a)")
plot2 <- plot2 + labs(tag = "(b)")
plot3 <- plot3 + labs(tag = "(c)")
plot4 <- plot4 + labs(tag = "(d)")
plot5 <- plot5 + labs(tag = "(e)")

# Extract the common legend
get_legend <- function(a_plot) {
  legend <- ggplotGrob(a_plot + theme(legend.position = "bottom"))$grobs
  legend <- legend[[which(sapply(legend, function(x) x$name) == "guide-box")]]
  return(legend)
}

legend <- get_legend(ggplot(data, aes(x = Ocean, y = NRI_and_NTI, fill = Metrics)) +
                       geom_bar(stat = "identity", position = "dodge") +
                       scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +
                       theme(legend.position = "bottom", legend.text = element_text(size = 25), 
                             legend.title = element_text(size = 25)))

# Arrange the plots (3 smaller on left, 2 larger on right)
final_layout <- ((plot1 / plot2 / plot3) | (plot4 / plot5)) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Save as SVG
svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/12_merged_plots.svg", width = 14, height = 18)
print(final_layout)
#grid::grid.draw(legend)  # Draw the common legend
dev.off()

# Display the layout in RStudio
final_layout

###################################################
###################################################

# Function to standardize bar plots
standardize_bar_chart <- function(plot) {
  plot +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
      axis.text.y = element_text(size = 25),
      axis.title.x = element_text(size = 25),
      axis.title.y = element_text(size = 25),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.position = "top",
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 25),
      plot.tag = element_text(size = 30, face = "bold")  # Increase tag size
    ) +
    scale_y_continuous(limits = c(-4, 63))  # Set a common Y-axis limit for consistency
}

###################################################

# Load first dataset
data1 <- read.csv(file = "class.csv", header = TRUE, stringsAsFactors = FALSE)
data1$Class <- factor(data1$Class, levels = c("Acidimicrobiia", "Actinomycetes", "Alphaproteobacteria",
                                              "Anaerolineae", "Bacteroidia", "Calditrichia",
                                              "Campylobacteria", "Cyanobacteriia", "Desulfobulbia",
                                              "Dissulfuribacteria", "Gammaproteobacteria", "Gracilibacteria",
                                              "Ktedonobacteria", "Leptospirae", "SAR324", "Verrucomicrobiae",
                                              "Zetaproteobacteria"))

plot1 <- ggplot(data1, aes(x = Class, y = NRI_and_NTI, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Standardized bar width
  geom_errorbar(aes(ymin = NRI_and_NTI - Std_dev, ymax = NRI_and_NTI + Std_dev), 
                position = position_dodge(width = 0.6), width = 0.2) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +
  scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +
  coord_flip() +
  labs(title = "Class NRI and NTI", x = "Class", y = "NRI_and_NTI", fill = "Metrics")

plot1 <- standardize_bar_chart(plot1)



# Save plots
svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/02_class_nri_nti.svg", width = 12, height = 12)
print(plot1)
dev.off()



# Display plots side by side for comparison
grid.arrange(plot1, plot2, ncol = 2)

data <- read.csv(file = "eez.csv", header = TRUE, stringsAsFactors = FALSE)

# Ensure Salinity is a factor with proper ordering
data$EEZ <- factor(data$EEZ, levels = c("Algerian", "Australian", "British (Bermuda)",
                                        "British (Cayman Islands)", "British (Pitcairn)",
                                        "Canadian", "Cape Verdean", "Chilean", "Chilean (Easter Island)",
                                        "Chinese", "Croatian", "Cypriote", "Danish", "Danish (Greenland)",
                                        "French", "French (French Polynesia)", "German", "Greek", "Indian",
                                        "Israeli", "Italian", "Japanese", "Latvian", "Madagascan", "Maldivian",
                                        "Maltese", "Mexican", "Micronesian", "New Zealand", "Norwegian",
                                        "Norwegian (Svalbard)", "Not_assigned", "Overlapping claim Taiwan: Taiwan / China",
                                        "Panamanian", "Papua New Guinean", "Portuguese (Azores)", "Portuguese (Madeira)",
                                        "Russian", "Saudi Arabian", "South African", "South African (Prince Edward Islands)",
                                        "Spanish", "Swedish", "Tongan", "United States",
                                        "United States (Alaska)", "United States (Guam)",
                                        "United States (Hawaii)", "United States (Puerto Rico)"))

svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/04_eez_nri_nti.svg", width = 14, height = 17)

# Create the horizontal grouped divergent bar plot
ggplot(data, aes(x = EEZ, y = NRI_and_NTI, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Grouped bars
  geom_errorbar(aes(ymin = NRI_and_NTI - Std_dev, ymax = NRI_and_NTI + Std_dev), 
                position = position_dodge(width = 0.6), width = 0.2) +  # Error bars
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +  # Dashed reference lines
  scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +  # Custom colors
  coord_flip() +  # Flip the plot to horizontal
  theme_minimal() +
  labs(title = "",
       x = "Exclusive Economic Zones",
       y = "NRI_and_NTI",
       fill = "Metrics") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.position = "top",
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  )

dev.off() 

###########

data <- read.csv(file = "MarRegion.csv", header = TRUE, stringsAsFactors = FALSE)

# Ensure Salinity is a factor with proper ordering
data$MarRegion <- factor(data$MarRegion, levels = c("Algerian part of the Mediterranean Sea - Western Basin", "Australian part of the Indian Ocean",
                                        "Australian part of the Tasman Sea", "British part of the Caribbean Sea (Cayman Islands)",
                                        "British part of the North Atlantic Ocean (Bermuda)", "British part of the South Pacific Ocean (Pitcairn)",
                                        "Canadian part of the Baffin Bay", "Canadian part of the Coastal Waters of Southeast Alaska and British Columbia",
                                        "Canadian part of the North Pacific Ocean", "Canadian part of the Northwestern Passages",
                                        "Cape Verdean part of the North Atlantic Ocean", "Chilean part of the South Pacific Ocean",
                                        "Chilean part of the South Pacific Ocean (Easter Island)", "Chinese part of the Eastern China Sea",
                                        "Chinese part of the South China Sea", "Croatian part of the Adriatic Sea", 
                                        "Cypriote part of the Mediterranean Sea - Eastern Basin", "Danish part of the Baffin Bay (Greenland)",
                                        "Danish part of the Davis Strait (Greenland)", "Danish part of the Kattegat",
                                        "Danish part of the North Sea", "French part of the Mediterranean Sea - Western Basin",
                                        "French part of the South Pacific Ocean (French Polynesia)", "German part of the Baltic Sea",
                                        "German part of the North Sea", "Greek part of the Aegean Sea",
                                        "Greek part of the Ionian Sea", "High Seas of the Arabian Sea", "High Seas of the Arctic Ocean",
                                        "High Seas of the Bering Sea", "High Seas of the Greenland Sea", "High Seas of the Indian Ocean",
                                        "High Seas of the Labrador Sea", "High Seas of the North Atlantic Ocean",
                                        "High Seas of the North Pacific Ocean", "High Seas of the Norwegian Sea", 
                                        "High Seas of the South Atlantic Ocean", "High Seas of the South Pacific Ocean",
                                        "High Seas of the Southern Ocean", "Indian part of the Arabian Sea",
                                        "Israeli part of the Mediterranean Sea - Eastern Basin", "Italian part of the Ionian Sea",
                                        "Italian part of the Mediterranean Sea - Eastern Basin", "Italian part of the Tyrrhenian Sea",
                                        "Japanese part of the North Pacific Ocean", "Latvian part of the Baltic Sea",
                                        "Madagascan part of the Mozambique Channel", "Maldivian part of the Indian Ocean",
                                        "Maldivian part of the Laccadive Sea", "Maltese part of the Mediterranean Sea - Eastern Basin",
                                        "Mexican part of the Gulf of California", "Mexican part of the Gulf of Mexico",
                                        "Mexican part of the North Pacific Ocean", "Micronesian part of the North Pacific Ocean",
                                        "New Zealand part of the South Pacific Ocean", "Norwegian part of the Norwegian Sea",
                                        "Norwegian part of the Norwegian Sea (Svalbard)", "Overlapping claim Taiwan: Taiwan / China part of the Philippine Sea",
                                        "Papua New Guinean part of the Bismarck Sea", "Portuguese part of the North Atlantic Ocean (Azores)",
                                        "Portuguese part of the North Atlantic Ocean (Madeira)", "Russian part of the Barentsz Sea",
                                        "Russian part of the East Siberian Sea", "Russian part of the Kara Sea",
                                        "Russian part of the Laptev Sea", "Saudi Arabian part of the Red Sea", 
                                        "South African part of the Indian Ocean", "South African part of the Indian Ocean (Prince Edward Islands)",
                                        "South African part of the South Atlantic Ocean", "Spanish part of the Balearic (Iberian Sea)",
                                        "Spanish part of the Mediterranean Sea - Western Basin", "Spanish part of the North Atlantic Ocean",
                                        "Swedish part of the Baltic Sea", "Swedish part of the Gulf of Bothnia", "Tongan part of the South Pacific Ocean",
                                        "United States part of the Arctic Ocean (Alaska)", "United States part of the Caribbean Sea (Puerto Rico)",
                                        "United States part of the Gulf of Mexico", "United States part of the North Atlantic Ocean",
                                        "United States part of the North Pacific Ocean", "United States part of the North Pacific Ocean (Guam)",
                                        "United States part of the North Pacific Ocean (Hawaii)", "United States part of the Philippine Sea (Guam)"))

svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/05_MarRegion_nri_nti.svg", width = 18, height = 34)

# Create the horizontal grouped divergent bar plot
ggplot(data, aes(x = MarRegion, y = NRI_and_NTI, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Grouped bars
  geom_errorbar(aes(ymin = NRI_and_NTI - Std_dev, ymax = NRI_and_NTI + Std_dev), 
                position = position_dodge(width = 0.6), width = 0.2) +  # Error bars
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +  # Dashed reference lines
  scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +  # Custom colors
  coord_flip() +  # Flip the plot to horizontal
  theme_minimal() +
  labs(title = "",
       x = "Marine Regions",
       y = "NRI_and_NTI",
       fill = "Metrics") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.position = "top",
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  )

dev.off() 

###########


###########

data <- read.csv(file = "Order.csv", header = TRUE, stringsAsFactors = FALSE)

# Ensure Salinity is a factor with proper ordering
data$Order <- factor(data$Order, levels = c("21-64-14", "Acetobacterales", "Acidiferrobacterales", "Acidimicrobiales",
                                            "Actinomycetales", "AKS1", "Anaerolineales", "Arenicellales", "Beggiatoales",
                                            "BMS3Bbin11", "Burkholderiales", "BW-2", "CACIAM-69d", "CAIVPW01", "CAKZMO01",
                                            "CALZJG01", "Campylobacterales", "Caulobacterales", "Chitinophagales", "Chromatiales",
                                            "Cyanobacteriales", "Cytophagales", "DAOVJZ01", "Desulfobulbales", "Dissulfuribacterales",
                                            "DRMK01", "DSM-19610", "DSM-26407", "Ectothiorhodospirales", "Enterobacterales",
                                            "Flavobacteriales", "GCA-001735895", "GCA-2400775", "GCF-002020875", "Geminicoccales",
                                            "GRL18", "Halothiobacillales", "JAADGQ01", "JAADHQ01", "JAADHS01", "JAAOII01", "JACKNK01",
                                            "JAFLJP01", "JAJDOJ01", "Kiloniellales", "Ktedonobacterales", "Leptospirales",
                                            "Mariprofundales", "Methylococcales", "Nevskiales", "Nitrosococcales", "PCC-6307",
                                            "Phormidesmiales", "Promineifilales", "Propionibacteriales", "PS1", "Pseudomonadales",
                                            "QNFN01", "RBG-13-44-9", "Rhizobiales", "Rhodobacterales", "Rhodospirillales",
                                            "S012-40", "S36-B12", "SAR324", "SG8-11", "Sphingomonadales", "Steroidobacterales",
                                            "SZUA-140", "SZUA-152", "SZUA-229", "SZUA-76", "Tenderiales", "Thermostichales",
                                            "Thiomicrospirales", "Thiotrichales", "UBA1369", "UBA2966", "UBA5794", "UBA6429",
                                            "UBA8366", "UBA9214", "Unclassified_Gammaproteobacteria", "UWMA-0217", "Woeseiales",
                                            "Xanthomonadales", "XN24"))

svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/07_Order_nri_nti.svg", width = 15, height = 34)

# Create the horizontal grouped divergent bar plot
ggplot(data, aes(x = Order, y = NRI_and_NTI, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Grouped bars
  geom_errorbar(aes(ymin = NRI_and_NTI - Std_dev, ymax = NRI_and_NTI + Std_dev), 
                position = position_dodge(width = 0.6), width = 0.2) +  # Error bars
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +  # Dashed reference lines
  scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +  # Custom colors
  coord_flip() +  # Flip the plot to horizontal
  theme_minimal() +
  labs(title = "",
       x = "Order",
       y = "NRI_and_NTI",
       fill = "Metrics") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.position = "top",
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  )

dev.off() 

###########



###########

data <- read.csv(file = "sea.csv", header = TRUE, stringsAsFactors = FALSE)

# Ensure Salinity is a factor with proper ordering
data$IHO_sea <- factor(data$IHO_sea, levels = c("Adriatic Sea", "Aegean Sea", "Arabian Sea", "Arctic Ocean",
                                                  "Baffin Bay", "Balearic (Iberian Sea)", "Baltic Sea", "Barentsz Sea",
                                                  "Bering Sea", "Bismarck Sea", "Caribbean Sea", "Davis Strait",
                                                  "East Siberian Sea", "Eastern China Sea", "Greenland Sea",
                                                  "Gulf of Bothnia", "Gulf of California", "Gulf of Mexico",
                                                  "Indian Ocean", "Ionian Sea", "Kara Sea", "Kattegat", 
                                                  "Labrador Sea", "Laccadive Sea", "Laptev Sea", "Mediterranean Sea - Eastern Basin",
                                                  "Mediterranean Sea - Western Basin", "Mozambique Channel", "North Atlantic Ocean",
                                                  "North Pacific Ocean", "North Sea", "Norwegian Sea", "Philippine Sea",
                                                  "Red Sea", "South Atlantic Ocean", "South China Sea", "South Pacific Ocean",
                                                  "Southern Ocean", "Tasman Sea", "The Coastal Waters of Southeast Alaska and British Columbia",
                                                  "The Northwestern Passages", "Tyrrhenian Sea"))

svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/10_IHO_sea_nri_nti.svg", width = 16, height = 24)

# Create the horizontal grouped divergent bar plot
ggplot(data, aes(x = IHO_sea, y = NRI_and_NTI, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +  # Grouped bars
  geom_errorbar(aes(ymin = NRI_and_NTI - Std_dev, ymax = NRI_and_NTI + Std_dev), 
                position = position_dodge(width = 0.6), width = 0.2) +  # Error bars
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +  # Dashed reference lines
  scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +  # Custom colors
  coord_flip() +  # Flip the plot to horizontal
  theme_minimal() +
  labs(title = "",
       x = "Sea",
       y = "NRI_and_NTI",
       fill = "Metrics") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.position = "top",
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25)
  )

dev.off() 

##############################################################
##############################################################
##############################################################

#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis")

# Load your data
data <- read.csv("salinity_categories.csv")

# Ensure column names are clean
names(data) <- tolower(gsub("\\.", "_", names(data)))

# Separate NRI and NTI
nri_data <- subset(data, metrics == "NRI")
nti_data <- subset(data, metrics == "NTI")

# Function to perform independent t-test from means, SDs, and sample sizes
ind_ttest <- function(mean1, sd1, n1, mean2, sd2, n2) {
  se <- sqrt((sd1^2 / n1) + (sd2^2 / n2))
  t_value <- (mean1 - mean2) / se
  df_num <- (sd1^2 / n1 + sd2^2 / n2)^2
  df_den <- ((sd1^2 / n1)^2 / (n1 - 1)) + ((sd2^2 / n2)^2 / (n2 - 1))
  df <- df_num / df_den
  p_value <- 2 * pt(-abs(t_value), df)
  crit_t <- qt(0.975, df)
  return(c(t_value = t_value, p_value = p_value, df = df, critical_t = crit_t))
}

# Compare salinity groups for each metric
salinity_levels <- unique(data$salinity)
results <- data.frame()

for (metric_data in list(nri_data, nti_data)) {
  metric_name <- unique(metric_data$metrics)
  for (i in 1:(length(salinity_levels) - 1)) {
    for (j in (i + 1):length(salinity_levels)) {
      group1 <- subset(metric_data, salinity == salinity_levels[i])
      group2 <- subset(metric_data, salinity == salinity_levels[j])
      stats <- ind_ttest(group1$nri_and_nti, group1$std_dev, group1$number,
                         group2$nri_and_nti, group2$std_dev, group2$number)
      results <- rbind(results, data.frame(
        Metric = metric_name,
        Group1 = salinity_levels[i],
        Group2 = salinity_levels[j],
        t_value = round(stats["t_value"], 3),
        p_value = signif(stats["p_value"], 3),
        df = round(stats["df"]),
        critical_t = round(stats["critical_t"], 2)
      ))
    }
  }
}

print(results)

#########################
###student t/test independent

# Load your data
data <- read.csv("salinity_categories.csv")

# Clean column names for consistency
names(data) <- tolower(gsub("\\.", "_", names(data)))

# Separate into NRI and NTI data
nri_data <- subset(data, metrics == "NRI")
nti_data <- subset(data, metrics == "NTI")

# Define Student's t-test function (equal variances assumed)
student_ttest <- function(mean1, sd1, n1, mean2, sd2, n2) {
  # Pooled variance
  pooled_var <- (((n1 - 1) * sd1^2) + ((n2 - 1) * sd2^2)) / (n1 + n2 - 2)
  se <- sqrt(pooled_var * (1/n1 + 1/n2))
  t_val <- (mean1 - mean2) / se
  df <- n1 + n2 - 2
  p_val <- 2 * pt(-abs(t_val), df)
  crit_t <- qt(0.975, df)
  return(c(t_value = t_val, p_value = p_val, df = df, critical_t = crit_t))
}

# Perform pairwise comparisons
compare_groups <- function(metric_data, metric_name) {
  salinity_levels <- unique(metric_data$salinity)
  result_list <- list()
  
  for (i in 1:(length(salinity_levels) - 1)) {
    for (j in (i + 1):length(salinity_levels)) {
      g1 <- subset(metric_data, salinity == salinity_levels[i])
      g2 <- subset(metric_data, salinity == salinity_levels[j])
      
      stats <- student_ttest(
        g1$nri_and_nti, g1$std_dev, g1$number,
        g2$nri_and_nti, g2$std_dev, g2$number
      )
      
      result_list[[length(result_list) + 1]] <- data.frame(
        Metric = metric_name,
        Group1 = salinity_levels[i],
        Group2 = salinity_levels[j],
        t_value = round(stats["t_value"], 3),
        p_value = signif(stats["p_value"], 3),
        df = stats["df"],
        critical_t = round(stats["critical_t"], 2)
      )
    }
  }
  do.call(rbind, result_list)
}

# Run for both metrics
nri_results <- compare_groups(nri_data, "NRI")
nti_results <- compare_groups(nti_data, "NTI")

# Combine and print
final_results <- rbind(nri_results, nti_results)
print(final_results)

write.csv(final_results, file = "salinity_independent_t_test_results_df.csv")


################################
###student t/test dependent

# Load your data
data <- read.csv("salinity_categories.csv")

# Clean column names
names(data) <- tolower(gsub("\\.", "_", names(data)))

# Ensure Metrics column is properly interpreted
data$metrics <- toupper(data$metrics)

# Split NRI and NTI
nri_data <- subset(data, metrics == "NRI")
nti_data <- subset(data, metrics == "NTI")

# Sort to ensure matching order
nri_data <- nri_data[order(nri_data$salinity), ]
nti_data <- nti_data[order(nti_data$salinity), ]

# Paired t-test function (manual and built-in)
paired_ttest <- function(nri, nti, sd_nri, sd_nti, n, salinity) {
  # Difference
  diff <- nti - nri
  sd_diff <- sqrt(sd_nri^2 + sd_nti^2)  # Approximate
  se_diff <- sd_diff / sqrt(n)
  t_val <- diff / se_diff
  df <- n - 1
  p_val <- 2 * pt(-abs(t_val), df)
  crit_t <- qt(0.975, df)
  return(data.frame(
    Salinity = salinity,
    Mean_NRI = nri,
    Mean_NTI = nti,
    Mean_Difference = diff,
    t_value = round(t_val, 3),
    p_value = signif(p_val, 3),
    df = df,
    critical_t = round(crit_t, 2)
  ))
}

# Apply the test across rows
paired_results <- mapply(
  paired_ttest,
  nri = nri_data$nri_and_nti,
  nti = nti_data$nri_and_nti,
  sd_nri = nri_data$std_dev,
  sd_nti = nti_data$std_dev,
  n = nri_data$number,
  salinity = nri_data$salinity,
  SIMPLIFY = FALSE
)

# Combine all into one result table
paired_results_df <- do.call(rbind, paired_results)
print(paired_results_df)

write.csv(paired_results_df, file = "salinity_dependent_t_test_results_df.csv")

#################################
########### independent t-test

# Load your CSV
data <- read.csv("territory.csv")

# Clean column names
names(data) <- tolower(gsub("\\.", "_", names(data)))

# Convert required columns to numeric (in case they were read as factors/characters)
data$number <- as.numeric(as.character(data$number))
data$nri_and_nti <- as.numeric(as.character(data$nri_and_nti))
data$std_dev <- as.numeric(as.character(data$std_dev))


# Separate by metric
nri_data <- subset(data, metrics == "NRI")
nti_data <- subset(data, metrics == "NTI")

# Student's t-test function (enhanced)
student_ttest <- function(mean1, sd1, n1, mean2, sd2, n2, g1, g2, metric) {
  # Decide which group has higher mean
  if (mean1 >= mean2) {
    m_high <- mean1; s_high <- sd1; n_high <- n1; g_high <- g1
    m_low  <- mean2; s_low  <- sd2; n_low  <- n2; g_low  <- g2
  } else {
    m_high <- mean2; s_high <- sd2; n_high <- n2; g_high <- g2
    m_low  <- mean1; s_low  <- sd1; n_low  <- n1; g_low  <- g1
  }
  
  # Pooled variance and SE
  pooled_var <- (((n_high - 1) * s_high^2) + ((n_low - 1) * s_low^2)) / (n_high + n_low - 2)
  se <- sqrt(pooled_var * (1/n_high + 1/n_low))
  t_val <- (m_high - m_low) / se
  df <- n_high + n_low - 2
  p_val <- 2 * pt(-abs(t_val), df)
  crit_t <- qt(0.975, df)
  
  return(data.frame(
    Metric = metric,
    Group_High = g_high,
    Group_Low = g_low,
    Mean_High = round(m_high, 3),
    Mean_Low = round(m_low, 3),
    Mean_Difference = round(m_high - m_low, 3),
    t_value = round(t_val, 3),
    p_value = signif(p_val, 3),
    df = df,
    critical_t = round(crit_t, 2)
  ))
}

# Function to perform all pairwise comparisons for a given metric
compare_groups <- function(metric_data, metric_name) {
  salinity_levels <- unique(metric_data$salinity)
  result_list <- list()
  
  for (i in 1:(length(salinity_levels) - 1)) {
    for (j in (i + 1):length(salinity_levels)) {
      g1 <- subset(metric_data, salinity == salinity_levels[i])
      g2 <- subset(metric_data, salinity == salinity_levels[j])
      
      result_list[[length(result_list) + 1]] <- student_ttest(
        mean1 = g1$nri_and_nti, sd1 = g1$std_dev, n1 = g1$number,
        mean2 = g2$nri_and_nti, sd2 = g2$std_dev, n2 = g2$number,
        g1 = g1$salinity, g2 = g2$salinity,
        metric = metric_name
      )
    }
  }
  do.call(rbind, result_list)
}

# Run for both metrics
nri_results <- compare_groups(nri_data, "NRI")
nti_results <- compare_groups(nti_data, "NTI")

# Combine and export or print
final_results <- rbind(nri_results, nti_results)
print(final_results)

# Optional: write to CSV
write.csv(final_results, file = "territory_independent_t_test_results_df.csv")


###############################
##### t-test dependent 

# Load your data
data <- read.csv("territory.csv")

# Clean column names
names(data) <- tolower(gsub("\\.", "_", names(data)))

# Standardize the metrics column
data$metrics <- toupper(data$metrics)

# Split into NRI and NTI subsets
nri_data <- subset(data, metrics == "NRI")
nti_data <- subset(data, metrics == "NTI")

# Sort by salinity to align both datasets
nri_data <- nri_data[order(nri_data$salinity), ]
nti_data <- nti_data[order(nti_data$salinity), ]

# Ensure relevant columns are numeric
nri_data$nri_and_nti <- as.numeric(as.character(nri_data$nri_and_nti))
nti_data$nri_and_nti <- as.numeric(as.character(nti_data$nri_and_nti))
nri_data$std_dev <- as.numeric(as.character(nri_data$std_dev))
nti_data$std_dev <- as.numeric(as.character(nti_data$std_dev))
nri_data$number <- as.numeric(as.character(nri_data$number))

# Define paired t-test function
paired_ttest <- function(nri, nti, sd_nri, sd_nti, n, salinity) {
  diff <- nti - nri
  sd_diff <- sqrt(sd_nri^2 + sd_nti^2)  # Estimated SD of difference
  se_diff <- sd_diff / sqrt(n)
  t_val <- diff / se_diff
  df <- n - 1
  p_val <- 2 * pt(-abs(t_val), df)
  crit_t <- qt(0.975, df)
  
  return(data.frame(
    Salinity = salinity,
    Mean_NRI = round(nri, 3),
    Mean_NTI = round(nti, 3),
    Mean_Difference = round(diff, 3),
    t_value = round(t_val, 3),
    p_value = signif(p_val, 3),
    df = df,
    critical_t = round(crit_t, 2)
  ))
}

# Apply the test across all rows
paired_results <- mapply(
  paired_ttest,
  nri = nri_data$nri_and_nti,
  nti = nti_data$nri_and_nti,
  sd_nri = nri_data$std_dev,
  sd_nti = nti_data$std_dev,
  n = nri_data$number,
  salinity = nri_data$salinity,
  SIMPLIFY = FALSE
)

# Combine into a result data frame
paired_results_df <- do.call(rbind, paired_results)

# View results
print(paired_results_df)

# Optional: Export results
write.csv(paired_results_df, "territory_paired_ttest_nri_nti_results.csv", row.names = FALSE)































##################################
##################################
## New NRI and NTI plots
# ==============================
# Divergent grouped bars + stats
# ==============================
#install.packages("igraph")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(multcompView)
  library(igraph)
})

p_to_stars <- function(p) {
  ifelse(is.na(p), "ns",
         ifelse(p <= 0.001, "***",
                ifelse(p <= 0.01,  "**",
                       ifelse(p <= 0.05, "*", "ns")
                )
         )
  )
}

cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
  if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
  pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
  if (length(pcol) == 0) stop("Independent-test table has no p-value column.")
  pcol <- pcol[1]
  sub <- indep_df %>% filter(Metric == metric, Group_High %in% cats, Group_Low %in% cats)
  if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))
  pv <- sub[[pcol]]
  names(pv) <- paste(pmin(sub$Group_High, sub$Group_Low),
                     pmax(sub$Group_High, sub$Group_Low), sep = "-")
  cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
  out <- setNames(rep("", length(cats)), cats); out[names(cld)] <- cld; out
}

cld_single_letter <- function(indep_df, metric, cats, alpha = 0.05) {
  if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
  pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
  if (length(pcol) == 0) stop("Independent-test table has no p-value column.")
  pcol <- pcol[1]
  sub <- indep_df %>% filter(Metric == metric, Group_High %in% cats, Group_Low %in% cats)
  g <- igraph::make_empty_graph(n = length(cats), directed = FALSE)
  igraph::V(g)$name <- cats
  if (nrow(sub) > 0) {
    apply(sub, 1, function(r) {
      a <- r[["Group_High"]]; b <- r[["Group_Low"]]; p <- as.numeric(r[[pcol]])
      if (!is.na(p) && p > alpha) g <<- igraph::add_edges(g, c(as.character(a), as.character(b)))
    })
  }
  comp <- igraph::components(g)$membership
  labs <- LETTERS[seq_along(unique(comp))]
  setNames(labs[match(comp, sort(unique(comp)))], igraph::V(g)$name)
}

plot_divergent_nri_nti <- function(
    data_csv,
    paired_csv = NULL,
    indep_csv  = NULL,
    category_levels = c("High","Warm","Moderate","Cold"),
    cld_mode = c("true","single","none"),
    cld_add_spaces = TRUE,
    save_file = NULL
) {
  cld_mode <- match.arg(cld_mode)
  
  # ---- 1) Read data & normalize category column name ----
  data <- read.csv(data_csv, header = TRUE, stringsAsFactors = FALSE)
  cand <- intersect(c("Temperature_Categories","Salinity","salinity","Category"), names(data))
  if (length(cand) == 0) stop("Data table: no category column named Temperature_Categories/Salinity/salinity/Category.")
  names(data)[names(data) == cand[1]] <- "Temperature_Categories"
  data$Temperature_Categories <- factor(data$Temperature_Categories, levels = category_levels)
  
  # Optional stats tables
  paired <- if (!is.null(paired_csv) && file.exists(paired_csv)) read.csv(paired_csv, stringsAsFactors = FALSE) else NULL
  indep  <- if (!is.null(indep_csv)  && file.exists(indep_csv))  read.csv(indep_csv,  stringsAsFactors = FALSE) else NULL
  
  # Harmonize category col in paired
  if (!is.null(paired)) {
    pcand <- intersect(c("Temperature_Categories","Salinity","salinity","Category"), names(paired))
    if (length(pcand) == 0) stop("Paired table: no category column.")
    names(paired)[names(paired) == pcand[1]] <- "Temperature_Categories"
    pcol <- intersect(c("p_value","p.value","p-value"), names(paired))
    if (length(pcol) == 0) stop("Paired table: no p-value column.")
    paired <- paired %>%
      mutate(Temperature_Categories = factor(Temperature_Categories, levels = levels(data$Temperature_Categories)),
             stars = p_to_stars(.data[[pcol[1]]])) %>%
      select(Temperature_Categories, stars)
  }
  
  # ---- 2) Build base plot (your aesthetics) ----
  dodge_w <- 0.6
  plot3 <- ggplot(data, aes(x = Temperature_Categories, y = NRI_and_NTI, fill = Metrics)) +
    geom_bar(stat = "identity", position = position_dodge(width = dodge_w), width = 0.6) +
    geom_errorbar(aes(ymin = NRI_and_NTI - Std_dev, ymax = NRI_and_NTI + Std_dev),
                  position = position_dodge(width = dodge_w), width = 0.2) +
    geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "#CC79A7", linewidth = 1) +
    scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +
    coord_flip() +
    theme_minimal() +
    labs(title = "", x = "Depth", y = "NRI and NTI", fill = "Metrics")
  
  # ---- 3) Paired stars (if provided) ----
  if (!is.null(paired) && nrow(paired) > 0) {
    cat_pos <- data %>%
      group_by(Temperature_Categories) %>%
      summarise(
        top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, Std_dev, 0), na.rm = TRUE),
        bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, Std_dev, 0), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(y_pos = ifelse(top <= 0, bot - 0.15, top + 0.15))
    paired_ann <- left_join(paired, cat_pos, by = "Temperature_Categories")
    plot3 <- plot3 +
      geom_text(
        data = paired_ann,
        aes(x = Temperature_Categories, y = y_pos, label = stars),
        inherit.aes = FALSE,
        fontface = "bold", size = 4
      )
  }
  
  # ---- 4) Between-category letters (optional) ----
  if (!is.null(indep) && nrow(indep) > 0 && cld_mode != "none") {
    cats <- levels(data$Temperature_Categories)
    cld_fun <- switch(cld_mode, "true" = cld_letters_for_metric, "single" = cld_single_letter)
    letters_list <- lapply(c("NRI","NTI"), function(m) {
      lt <- cld_fun(indep, metric = m, cats = cats, alpha = 0.05)
      lbl <- unname(lt)
      if (cld_mode == "true" && isTRUE(cld_add_spaces)) lbl <- gsub("(?<=.)(?=.)", " ", lbl, perl = TRUE)
      tibble(Metrics = m, Temperature_Categories = names(lt), cld = lbl)
    })
    letters_df <- bind_rows(letters_list)
    
    letter_pos <- data %>%
      group_by(Temperature_Categories, Metrics) %>%
      summarise(
        y_top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, Std_dev, 0), na.rm = TRUE),
        y_bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, Std_dev, 0), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(y_pos = ifelse(y_top <= 0, y_bot - 0.10, y_top + 0.10))
    
    letters_ann <- left_join(letters_df, letter_pos,
                             by = c("Temperature_Categories","Metrics"))
    
    plot3 <- plot3 +
      geom_text(
        data = letters_ann,
        aes(x = Temperature_Categories, y = y_pos, label = cld, group = Metrics),
        position = position_dodge(width = dodge_w),
        size = 3.8
      )
  }
  
  if (exists("standardize_bar_chart")) plot3 <- standardize_bar_chart(plot3)
  
  if (!is.null(save_file)) {
    ext <- tolower(tools::file_ext(save_file))
    if (ext %in% c("png","jpg","jpeg","tiff","bmp")) ggsave(save_file, plot3, width = 8, height = 8, dpi = 300)
    else if (ext == "svg") ggsave(save_file, plot3, width = 8, height = 8)
    else if (ext == "pdf") ggsave(save_file, plot3, width = 8, height = 8,
                                  device = tryCatch(cairo_pdf, error = function(e) pdf))
    else ggsave(paste0(save_file, ".png"), plot3, width = 8, height = 8, dpi = 300)
  }
  
  return(plot3)
}


##############################

# If needed:
# install.packages(c("ggplot2","dplyr","tidyr","stringr","multcompView","igraph"))

# A) Paired stars only (no between-group letters)
temp_plot <- plot_divergent_nri_nti(
  data_csv   = "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/temperature_categories.csv",
  paired_csv = "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/temperature_paired_ttest_nri_nti_results.csv",
  indep_csv  = NULL,
  cld_mode   = "none",
#  save_file  = "/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/11_temperature_nri_nti.svg"
)
print(temp_plot)

# B) Paired + TRUE CLD (may show multiple letters like "A B")
temp_plot <- plot_divergent_nri_nti(
  data_csv   = "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/temperature_categories.csv",
  paired_csv = "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/temperature_paired_ttest_nri_nti_results.csv",
  indep_csv  = "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/temperature_independent_t_test_results_df.csv",
  cld_mode   = "true",      # real CLD
  cld_add_spaces = TRUE,    # render "AB" as "A B"
  save_file  = "/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/11_temperature_nri_nti.svg"
)
print(temp_plot)

# # C) Paired + SINGLE-LETTER grouping (exactly one letter per category)
# plot3 <- plot_divergent_nri_nti(
#   data_csv   = "temperature_categories.csv",
#   paired_csv = "temperature_paired_ttest_nri_nti_results.csv",
#   indep_csv  = "temperature_independent_t_test_results_df.csv",
#   cld_mode   = "single",    # one letter per category
#   save_file  = "/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/11_temperature_nri_nti.svg"
# )
# print(plot3)

depth_plot <- plot_divergent_nri_nti(
  data_csv   = "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/depth_categories.csv",
  paired_csv = "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/depth_paired_ttest_nri_nti_results.csv",
  indep_csv  = "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/depth_independent_t_test_results_df.csv",
  cld_mode   = "true",      # real CLD
  cld_add_spaces = TRUE,    # render "AB" as "A B"
  save_file  = "/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/11_temperature_nri_nti.svg"
)
print(depth_plot)
##################################################



































#############################################################################
# # ==============================
# # Temperature NRI/NTI — stars + CLD
# # ==============================
# 
# options(stringsAsFactors = FALSE)
# suppressPackageStartupMessages({
#   library(ggplot2)
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
#   library(tibble)
#   library(multcompView)
#   library(igraph)
#   library(tools)
# })
# 
# # ---------- PATHS (keeps your D: paths but falls back to /mnt/data) ----------
# choose_first_existing <- function(paths) {
#   ok <- paths[file.exists(paths)]
#   if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
#   ok[[1]]
# }
# data_csv   <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/temperature_categories.csv",
#   "/mnt/data/temperature_categories.csv"
# ))
# paired_csv <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/temperature_paired_ttest_nri_nti_results.csv",
#   "/mnt/data/temperature_paired_ttest_nri_nti_results.csv"
# ))
# indep_csv  <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/temperature_independent_t_test_results_df.csv",
#   "/mnt/data/temperature_independent_t_test_results_df.csv"
# ))
# 
# out_dir <- dirname(data_csv)
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # ---------- helpers ----------
# p_to_stars <- function(p) {
#   ifelse(is.na(p), "ns",
#          ifelse(p <= 0.001, "***",
#                 ifelse(p <= 0.01,  "**",
#                        ifelse(p <= 0.05, "*", "ns"))))
# }
# 
# cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
#   if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
#   # try multiple p-value spellings
#   pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
#   if (length(pcol) == 0) stop("Independent-test table has no p-value column named p_value/p.value/p-value.")
#   pcol <- pcol[1]
#   sub <- indep_df %>% filter(.data$Metric == metric,
#                              .data$Group_High %in% cats,
#                              .data$Group_Low %in% cats)
#   if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))
#   pv <- sub[[pcol]]
#   names(pv) <- paste(pmin(sub$Group_High, sub$Group_Low),
#                      pmax(sub$Group_High, sub$Group_Low), sep = "-")
#   cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
#   out <- setNames(rep("", length(cats)), cats); out[names(cld)] <- cld; out
# }
# 
# # ---------- data (read + normalize) ----------
# df <- read.csv(data_csv, check.names = FALSE)
# 
# # Category
# cat_col <- intersect(c("Temperature_Categories","Temperature","temperature",
#                        "Salinity","salinity","Category","Categories"),
#                      names(df))
# if (length(cat_col) == 0) stop("Data: no category column named Temperature_Categories/Temperature/Salinity/Category.")
# names(df)[names(df) == cat_col[1]] <- "Category"
# 
# # Metric
# mcol <- intersect(c("Metrics","Metric"), names(df))
# if (length(mcol) == 0) stop("Data: no 'Metrics'/'Metric' column.")
# names(df)[names(df) == mcol[1]] <- "Metrics"
# 
# # Value
# vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df))
# if (length(vcol) == 0) stop("Data: no value column (expected NRI_and_NTI/Value).")
# names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"
# 
# # Error (optional)
# ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
# if (length(ecol) == 0) {
#   warning("No error column found; using zeros (no error bars).")
#   df$.__err__ <- 0
# } else {
#   df$.__err__ <- df[[ecol[1]]]
# }
# 
# # Levels/order
# desired_levels <- c("High","Warm","Moderate","Cold")
# if (all(desired_levels %in% unique(df$Category))) {
#   df$Category <- factor(df$Category, levels = desired_levels)
# } else {
#   df$Category <- factor(df$Category, levels = unique(df$Category))
# }
# df$Metrics <- factor(df$Metrics, levels = c("NRI","NTI"))
# cat_levels <- levels(df$Category)  # <- use this inside mutate later (avoids df masking)
# 
# # ---------- paired stars ----------
# paired <- read.csv(paired_csv, check.names = FALSE)
# # normalize names
# names(paired) <- trimws(names(paired))
# pcat <- intersect(c("Category","Temperature_Categories","Temperature","Salinity","salinity"), names(paired))
# if (length(pcat) == 0) stop("Paired table: no category column.")
# names(paired)[names(paired) == pcat[1]] <- "Category"
# pcol <- intersect(c("p_value","p.value","p-value"), names(paired))
# if (length(pcol) == 0) stop("Paired table: no p-value column named p_value/p.value/p-value.")
# 
# paired <- paired %>%
#   mutate(Category = factor(trimws(as.character(.data$Category)), levels = cat_levels),
#          stars = p_to_stars(.data[[pcol[1]]])) %>%
#   select(Category, stars)
# 
# # Positions for stars
# cat_pos <- df %>%
#   group_by(Category) %>%
#   summarise(
#     top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#     bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(y_pos = ifelse(top <= 0, bot - 0.50, top + 0.50))
# 
# paired_ann <- left_join(paired, cat_pos, by = "Category")
# 
# # ---------- base plot ----------
# dodge_w <- 0.6
# fill_values <- c("NRI"="#332288","NTI"="#DDCC77")
# dash_at <- c(-2,2)
# 
# plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
#   geom_bar(stat = "identity", position = position_dodge(width = dodge_w), width = 0.6) +
#   geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
#                 position = position_dodge(width = dodge_w), width = 0.2) +
#   geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
#   scale_fill_manual(values = fill_values) +
#   coord_flip() +
#   theme_minimal() +
#   labs(title = "", x = "Temperature", y = "NRI and NTI", fill = "Metrics")
# 
# # ---------- (A) stars-only ----------
# plot_stars <- plot_base +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = 4
#   )
# 
# # ---------- (B) stars + CLD letters ----------
# # Read & sanitize independent tests
# indep <- read.csv(indep_csv, check.names = FALSE)
# nm <- names(indep)
# nm <- trimws(nm)
# nm[nm == "p-value"]  <- "p_value"
# nm[nm == "p.value"]  <- "p_value"
# nm[nm == "Group High"] <- "Group_High"
# nm[nm == "Group Low"]  <- "Group_Low"
# nm[nm == "Metric "]    <- "Metric"
# # fill any empty/NA colnames so dplyr won't error
# empty <- which(is.na(nm) | nm == "")
# if (length(empty) > 0) nm[empty] <- paste0("X", empty)
# names(indep) <- nm
# 
# # Build CLD labels with guard-rails (won't stop plotting if something's off)
# letters_df <- try({
#   cats <- levels(df$Category)
#   letters_list <- lapply(c("NRI","NTI"), function(m) {
#     lt <- cld_letters_for_metric(indep, metric = m, cats = cats, alpha = 0.05)
#     lbl <- gsub("(?<=.)(?=.)", " ", unname(lt), perl = TRUE)  # "AB" -> "A B"
#     tibble(Metrics = m, Category = names(lt), cld = lbl)
#   })
#   bind_rows(letters_list)
# }, silent = TRUE)
# 
# plot_cld <- plot_base +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = 4
#   )
# 
# if (!inherits(letters_df, "try-error") && is.data.frame(letters_df) && nrow(letters_df) > 0) {
#   # position per group for letters
#   letter_pos <- df %>%
#     group_by(Category, Metrics) %>%
#     summarise(
#       y_top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#       y_bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     mutate(y_pos = ifelse(y_top <= 0, y_bot - 0.65, y_top + 0.65))
#   letters_ann <- left_join(letters_df, letter_pos, by = c("Category","Metrics"))
#   
#   plot_cld <- plot_cld +
#     geom_text(
#       data = letters_ann,
#       aes(x = Category, y = y_pos, label = cld, group = Metrics),
#       position = position_dodge(width = dodge_w),
#       size = 3.8
#     )
# } else {
#   message("CLD letters were skipped due to column-name or data issues in the independent-tests file.")
# }
# 
# # ---------- save ----------
# file_stem_a <- file.path(out_dir, "temperature_nri_nti_stars")
# file_stem_b <- file.path(out_dir, "temperature_nri_nti_cld")
# 
# ggsave(paste0(file_stem_a, ".svg"), plot_stars, width = 8, height = 8)
# ggsave(paste0(file_stem_a, ".png"), plot_stars, width = 8, height = 8, dpi = 300)
# ggsave(paste0(file_stem_b, ".svg"), plot_cld,  width = 8, height = 8)
# ggsave(paste0(file_stem_b, ".png"), plot_cld,  width = 8, height = 8, dpi = 300)
# 
# message("Saved:\n  ", paste0(file_stem_a, c(".svg",".png"), collapse = "\n  "),
#         "\n  ", paste0(file_stem_b, c(".svg",".png"), collapse = "\n  "))
# 
# ################################
# 
# # (optional) small cleanup to avoid "Removed X rows" warnings
# paired_ann2  <- tryCatch(dplyr::filter(paired_ann, !is.na(y_pos), !is.na(stars)), error = function(e) paired_ann)
# letters_ann2 <- if (exists("letters_ann") && is.data.frame(letters_ann)) {
#   dplyr::filter(letters_ann, !is.na(y_pos), !is.na(cld), cld != "")
# } else NULL
# 
# plot_3 <- plot_base +
#   # paired t-test stars
#   geom_text(
#     data = paired_ann2,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold",
#     size = 3.8
#   ) +
#   # CLD letters (only if available)
#   {
#     if (!is.null(letters_ann2) && nrow(letters_ann2) > 0) {
#       geom_text(
#         data = letters_ann2,
#         aes(x = Category, y = y_pos, label = cld, group = Metrics),
#         position = position_dodge(width = dodge_w),
#         size = 3.4
#       )
#     } else {
#       NULL
#     }
#   }
# 
# # Use it later with: print(plot4)
# 
# # (optional) save to file when you want
# # ggsave(file.path(out_dir, "plot4.svg"), plot4, width = 10, height = 8)
# # ggsave(file.path(out_dir, "plot4.png"), plot4, width = 10, height = 8, dpi = 300)
# # (optional) save the R object:
# # saveRDS(plot4, file.path(out_dir, "plot4.rds"))
#######################################################

# ==============================
# Temperature NRI/NTI — stars + CLD (clean NA, colored CLD, Number>=5, ordering)
# ==============================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(multcompView)
  library(tools)
})

# ---------- paths ----------
choose_first_existing <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
  ok[[1]]
}
data_csv   <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/temperature_categories.csv",
  "/mnt/data/temperature_categories.csv"
))
paired_csv <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/temperature_paired_ttest_nri_nti_results.csv",
  "/mnt/data/temperature_paired_ttest_nri_nti_results.csv"
))
indep_csv  <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/temperature_independent_t_test_results_df.csv",
  "/mnt/data/temperature_independent_t_test_results_df.csv"
))

out_dir <- dirname(data_csv)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- knobs ----------
BASE      <- 22
TITLE     <- 24
STAR_SIZE <- 7
CLD_SIZE  <- 8.2
BAR_W     <- 0.65
DODGE_W   <- 0.60
dash_at   <- c(-2, 2)

fill_values <- c("NRI"="#332288","NTI"="#DDCC77")  # bar fills
cld_cols    <- c(NRI = "#332288", NTI = "#A68000") # CLD letters (colored by metric)

# Ordering: keep defined order if present, otherwise order by mean metric
KEEP_DEFINED_ORDER <- TRUE
defined_levels <- c("High","Warm","Moderate","Cold")
ORDER_METRIC   <- "NRI"   # or "NTI"
ORDER_DIR      <- "desc"  # "asc" or "desc"

# ---------- helpers ----------
p_to_stars <- function(p) ifelse(is.na(p), "ns",
                                 ifelse(p <= 0.001, "***",
                                        ifelse(p <= 0.01,  "**",
                                               ifelse(p <= 0.05, "*", "ns"))))

to_numeric <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\u2212", "-", x, fixed = TRUE)   # Unicode minus → ASCII
  x <- gsub("[^0-9eE+\\-\\.]", "", x)
  suppressWarnings(as.numeric(x))
}

safe_name <- function(x) { x <- trimws(as.character(x)); gsub("-", "–", x, fixed = TRUE) }

cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
  if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
  pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
  if (!length(pcol)) stop("Independent-tests table has no p-value column p_value/p.value/p-value.")
  pcol <- pcol[1]
  
  sub <- indep_df %>% filter(.data$Metric == metric,
                             .data$Group_High %in% cats,
                             .data$Group_Low  %in% cats)
  if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))
  
  n_hi <- safe_name(sub$Group_High)
  n_lo <- safe_name(sub$Group_Low)
  pair_names <- paste(pmin(n_hi, n_lo), pmax(n_hi, n_lo), sep = "-")
  pv <- sub[[pcol]]; names(pv) <- pair_names
  
  cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
  
  lookup <- setNames(cats, safe_name(cats))
  out <- setNames(rep("", length(cats)), cats)
  mapped <- intersect(names(cld), names(lookup))
  if (length(mapped)) out[lookup[mapped]] <- cld[mapped]
  out
}

# ---------- read + normalize ----------
df <- read.csv(data_csv, check.names = FALSE)
names(df) <- trimws(names(df))

# Category
cat_col <- intersect(c("Temperature_Categories","Temperature","temperature",
                       "Salinity","salinity","Category","Categories"),
                     names(df))
if (!length(cat_col)) stop("Data: no category column named Temperature/Category.")
names(df)[names(df) == cat_col[1]] <- "Category"

# Metric / Value / Error
mcol <- intersect(c("Metrics","Metric"), names(df)); if (!length(mcol)) stop("No 'Metrics' column.")
names(df)[names(df) == mcol[1]] <- "Metrics"
vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df)); if (!length(vcol)) stop("No value column (NRI_and_NTI/Value).")
names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"
ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]

# Filter categories with Number < 5 (auto-detect count column)
num_col_idx <- grep("(?i)^(n$|n_|^n\\.|^n\\s*$|number|num|count|samples|n_samples|n_samp|n\\.samples|n\\.samp)$|(?i)number", names(df), perl = TRUE)
if (length(num_col_idx) > 0) {
  num_col <- names(df)[num_col_idx[1]]
  df[[num_col]] <- to_numeric(df[[num_col]])
  df <- df %>% filter(.data[[num_col]] >= 5)
} else {
  message("No 'Number'/'Count' column found — skipping <5 filter.")
}

# Clean rows and coerce numerics; keep only NRI/NTI; drop NA/blank/"NA"
df <- df %>%
  mutate(
    Category    = trimws(as.character(Category)),
    Metrics     = trimws(as.character(Metrics)),
    NRI_and_NTI = to_numeric(NRI_and_NTI),
    .__err__    = to_numeric(.__err__)
  ) %>%
  filter(!is.na(Category), nzchar(Category), toupper(Category) != "NA",
         !is.na(Metrics), Metrics %in% c("NRI","NTI"))

# ---------- order categories ----------
if (KEEP_DEFINED_ORDER && all(defined_levels %in% unique(df$Category))) {
  ord_levels <- defined_levels
} else {
  ord_tbl <- df %>%
    group_by(Category) %>%
    summarise(order_mean = mean(ifelse(Metrics == ORDER_METRIC, NRI_and_NTI, NA), na.rm = TRUE),
              .groups = "drop") %>%
    arrange(if (ORDER_DIR == "asc") order_mean else desc(order_mean))
  ord_levels <- ord_tbl$Category
}
df$Category <- factor(df$Category, levels = ord_levels)
df$Metrics  <- factor(df$Metrics,  levels = c("NRI","NTI"))
df <- droplevels(df)

# ---------- paired stars ----------
paired <- read.csv(paired_csv, check.names = FALSE)
names(paired) <- trimws(names(paired))
pcat <- intersect(c("Category","Temperature_Categories","Temperature","Salinity","salinity"), names(paired))
if (!length(pcat)) stop("Paired table: no category column.")
names(paired)[names(paired) == pcat[1]] <- "Category"
pcol <- intersect(c("p_value","p.value","p-value"), names(paired)); if (!length(pcol)) stop("Paired table: no p-value column.")

paired <- paired %>%
  mutate(
    Category = trimws(as.character(.data$Category)),
    stars    = p_to_stars(.data[[pcol[1]]])
  ) %>%
  filter(Category %in% ord_levels, nzchar(Category), toupper(Category) != "NA") %>%
  mutate(Category = factor(Category, levels = ord_levels)) %>%
  select(Category, stars) %>%
  droplevels()

# ---------- dynamic positions & limits ----------
y_min_data <- min(df$NRI_and_NTI - ifelse(df$NRI_and_NTI < 0, abs(df$.__err__), 0), na.rm = TRUE)
y_max_data <- max(df$NRI_and_NTI + ifelse(df$NRI_and_NTI > 0, abs(df$.__err__), 0), na.rm = TRUE)
rng_base   <- range(c(y_min_data, y_max_data, dash_at), na.rm = TRUE)

pad <- max(diff(rng_base) * 0.06, 0.15)  # gap bars↔labels
tol <- max(diff(rng_base) * 0.02, 0.10)  # avoid dashed lines

cat_pos <- df %>%
  group_by(Category) %>%
  summarise(
    top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
    bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_star = ifelse(top <= 0, bot - pad, top + pad))

for (d in dash_at) {
  cat_pos$y_star <- ifelse(abs(cat_pos$y_star - d) < tol,
                           ifelse(cat_pos$y_star >= d, cat_pos$y_star + pad*0.5, cat_pos$y_star - pad*0.5),
                           cat_pos$y_star)
}
paired_ann <- left_join(paired, cat_pos, by = "Category")

# ---------- base plot ----------
plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
  geom_col(position = position_dodge(width = DODGE_W), width = BAR_W) +
  geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
                position = position_dodge(width = DODGE_W), width = 0.22, linewidth = 0.7) +
  geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
  scale_fill_manual(values = fill_values, name = "Metrics") +
  coord_flip() +
  labs(title = "", x = "Temperature", y = "NRI and NTI") +
  theme_minimal(base_size = BASE) +
  theme(
    axis.title.x = element_text(size = TITLE, face = "bold"),
    axis.title.y = element_text(size = TITLE, face = "bold"),
    axis.text.x  = element_text(size = BASE),
    axis.text.y  = element_text(size = BASE),
    legend.title = element_text(size = TITLE, face = "bold"),
    legend.text  = element_text(size = BASE),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 14, 10, 10)
  )

# ---------- (A) stars-only (with halo) ----------
plot_stars <- plot_base +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE + 1.2, color = "white"
  ) +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE, color = "black"
  )

# ---------- (B) stars + CLD letters (colored, with halo) ----------
indep <- read.csv(indep_csv, check.names = FALSE)
nm <- trimws(names(indep)); nm[!nzchar(nm) | is.na(nm)] <- paste0("X", which(!nzchar(nm) | is.na(nm)))
names(indep) <- nm
nm <- names(indep)
nm <- sub("^p[._-]?value$", "p_value", nm, ignore.case = TRUE)
nm <- sub("^Group[ _]High$", "Group_High", nm)
nm <- sub("^Group[ _]Low$",  "Group_Low",  nm)
nm <- sub("^Metric[ ]*$",    "Metric",     nm)
names(indep) <- nm
drop_cols <- grep("^(Unnamed: ?\\d+|X\\d+)$", names(indep))
if (length(drop_cols)) indep <- indep[, -drop_cols, drop = FALSE]

req <- c("Metric","Group_High","Group_Low","p_value")
missing <- setdiff(req, names(indep))
if (length(missing)) stop("Independent-tests table missing: ", paste(missing, collapse = ", "))

indep <- indep %>%
  mutate(
    Metric     = trimws(as.character(Metric)),
    Group_High = trimws(as.character(Group_High)),
    Group_Low  = trimws(as.character(Group_Low))
  )

cats <- levels(df$Category)
letters_df <- bind_rows(
  { lt <- cld_letters_for_metric(indep, "NRI", cats, alpha = 0.05)
  tibble(Metrics="NRI", Category=names(lt),
         cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) },
  { lt <- cld_letters_for_metric(indep, "NTI", cats, alpha = 0.05)
  tibble(Metrics="NTI", Category=names(lt),
         cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) }
) %>%
  filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA") %>%
  filter(Category %in% ord_levels) %>%
  mutate(Category = factor(Category, levels = ord_levels)) %>%
  droplevels()

plot_cld <- plot_stars  # includes halo stars already

if (nrow(letters_df) > 0) {
  letter_pos <- df %>%
    group_by(Category, Metrics) %>%
    summarise(
      y_top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
      y_bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_cld = ifelse(y_top <= 0, y_bot - pad*1.25, y_top + pad*1.25))
  
  for (d in dash_at) {
    letter_pos$y_cld <- ifelse(abs(letter_pos$y_cld - d) < tol,
                               ifelse(letter_pos$y_cld >= d, letter_pos$y_cld + pad*0.5, letter_pos$y_cld - pad*0.5),
                               letter_pos$y_cld)
  }
  star_map <- cat_pos %>% select(Category, y_star)
  letter_pos <- left_join(letter_pos, star_map, by = "Category") %>%
    mutate(y_cld = ifelse(!is.na(y_star) & abs(y_cld - y_star) < (pad*0.6),
                          ifelse(y_cld >= y_star, y_cld + pad*0.5, y_cld - pad*0.5),
                          y_cld)) %>%
    select(-y_star)
  
  letters_ann <- inner_join(letters_df, letter_pos, by = c("Category","Metrics")) %>%
    filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA")
  
  plot_cld <- plot_cld +
    geom_text(
      data = letters_ann,
      aes(x = Category, y = y_cld, label = cld, group = Metrics),
      position = position_dodge(width = DODGE_W),
      size = CLD_SIZE + 1.4, fontface = "bold", color = "white"
    ) +
    geom_text(
      data = letters_ann,
      aes(x = Category, y = y_cld, label = cld, group = Metrics, color = Metrics),
      position = position_dodge(width = DODGE_W),
      size = CLD_SIZE, fontface = "bold"
    ) +
    scale_color_manual(values = cld_cols, guide = "none")
} else {
  message("CLD letters are empty (no matching categories or p-values).")
}

# ---------- y-limits (bars + stars + letters + dashed lines) ----------
all_y <- c(df$NRI_and_NTI + df$.__err__, df$NRI_and_NTI - df$.__err__,
           paired_ann$y_star,
           if (exists("letters_ann")) letters_ann$y_cld else numeric(0),
           dash_at)
ymin <- floor((min(all_y, na.rm = TRUE) - pad*0.4) * 10) / 10
ymax <- ceiling((max(all_y, na.rm = TRUE) + pad*0.4) * 10) / 10

plot_stars <- plot_stars + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))
plot_cld   <- plot_cld   + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))

# ---------- save ----------
file_stem_a <- file.path(out_dir, sprintf("temperature_nri_nti_stars_orderBy%s_%s", ORDER_METRIC, ORDER_DIR))
file_stem_b <- file.path(out_dir, sprintf("temperature_nri_nti_cld_orderBy%s_%s",   ORDER_METRIC, ORDER_DIR))

ggsave(paste0(file_stem_a, ".svg"), plot_stars, width = 10, height = 10)
ggsave(paste0(file_stem_a, ".png"), plot_stars, width = 10, height = 10, dpi = 300)
ggsave(paste0(file_stem_b, ".svg"), plot_cld,  width = 10, height = 10)
ggsave(paste0(file_stem_b, ".png"), plot_cld,  width = 10, height = 10, dpi = 300)

message("Saved:\n  ", paste0(file_stem_a, c(".svg",".png"), collapse = "\n  "),
        "\n  ", paste0(file_stem_b, c(".svg",".png"), collapse = "\n  "))

# optional: keep final object
plot3 <- plot_cld

###############################
# # ==============================
# # Depth NRI/NTI — stars + CLD (clean NA, colored letters, Number>=5)
# # ==============================
# 
# options(stringsAsFactors = FALSE)
# suppressPackageStartupMessages({
#   library(ggplot2)
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
#   library(tibble)
#   library(multcompView)
#   library(tools)
# })
# 
# # ---------- paths ----------
# choose_first_existing <- function(paths) {
#   ok <- paths[file.exists(paths)]
#   if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
#   ok[[1]]
# }
# data_csv   <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/depth_categories.csv",
#   "/mnt/data/depth_categories.csv"
# ))
# paired_csv <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/depth_paired_ttest_nri_nti_results.csv",
#   "/mnt/data/depth_paired_ttest_nri_nti_results.csv"
# ))
# indep_csv  <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/depth_independent_t_test_results_df.csv",
#   "/mnt/data/depth_independent_t_test_results_df.csv"
# ))
# 
# out_dir <- dirname(data_csv)
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # ---------- knobs ----------
# BASE      <- 22
# TITLE     <- 24
# STAR_SIZE <- 7
# CLD_SIZE  <- 8.2
# BAR_W     <- 0.65
# DODGE_W   <- 0.60
# dash_at   <- c(-2, 2)
# 
# fill_values <- c("NRI"="#332288","NTI"="#DDCC77")  # bar fills
# cld_cols    <- c(NRI = "#332288", NTI = "#A68000") # CLD colors (with halo)
# 
# # Ordering: by mean NRI or NTI, or keep ecological order if present
# ORDER_METRIC         <- "NRI"   # or "NTI"
# ORDER_DIR            <- "desc"  # "asc" or "desc"
# KEEP_ECOLOGICAL_ORDER <- TRUE
# desired_levels <- c("Shallow Zone","Epipelagic","Mesopelagic","Bathypelagic","Abyssopelagic","Hadalpelagic")
# 
# # ---------- helpers ----------
# p_to_stars <- function(p) ifelse(is.na(p), "ns",
#                            ifelse(p <= 0.001, "***",
#                            ifelse(p <= 0.01,  "**",
#                            ifelse(p <= 0.05, "*", "ns"))))
# 
# safe_name <- function(x) { x <- trimws(as.character(x)); gsub("-", "–", x, fixed = TRUE) }
# 
# cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
#   if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
#   pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
#   if (length(pcol) == 0) stop("Independent-tests table has no p-value column p_value/p.value/p-value.")
#   pcol <- pcol[1]
# 
#   sub <- indep_df %>% filter(.data$Metric == metric,
#                              .data$Group_High %in% cats,
#                              .data$Group_Low  %in% cats)
#   if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))
# 
#   n_hi <- safe_name(sub$Group_High)
#   n_lo <- safe_name(sub$Group_Low)
#   pair_names <- paste(pmin(n_hi, n_lo), pmax(n_hi, n_lo), sep = "-")
#   pv <- sub[[pcol]]; names(pv) <- pair_names
# 
#   cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
# 
#   lookup <- setNames(cats, safe_name(cats))
#   out <- setNames(rep("", length(cats)), cats)
#   mapped <- intersect(names(cld), names(lookup))
#   if (length(mapped)) out[lookup[mapped]] <- cld[mapped]
#   out
# }
# 
# to_numeric <- function(x) {
#   x <- trimws(as.character(x))
#   x <- gsub("\u2212", "-", x, fixed = TRUE)
#   x <- gsub("[^0-9eE+\\-\\.]", "", x)
#   suppressWarnings(as.numeric(x))
# }
# 
# # ---------- read + normalize (DEPTH) ----------
# df <- read.csv(data_csv, check.names = FALSE)
# names(df) <- trimws(names(df))
# 
# # Category
# cat_col <- intersect(c("Depth_Categories","Depth_Category","Depth","depth","Category","Categories","salinity","Salinity"), names(df))
# if (length(cat_col) == 0) stop("Data: couldn't find a depth/category column. Available: ", paste(names(df), collapse = " | "))
# names(df)[names(df) == cat_col[1]] <- "Category"
# 
# # Metric, Value, Error
# mcol <- intersect(c("Metrics","Metric"), names(df)); if (!length(mcol)) stop("No 'Metrics' column.")
# names(df)[names(df) == mcol[1]] <- "Metrics"
# vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df)); if (!length(vcol)) stop("No value column (NRI_and_NTI/Value).")
# names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"
# ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
# df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]
# 
# # ---- filter categories with Number < 5 (auto-detect count column) ----
# num_col_idx <- grep("(?i)^(n$|n_|^n\\.|^n\\s*$|number|num|count|samples|n_samples|n_samp|n\\.samples|n\\.samp)$|(?i)number", names(df), perl = TRUE)
# if (length(num_col_idx) > 0) {
#   num_col <- names(df)[num_col_idx[1]]
#   df[[num_col]] <- to_numeric(df[[num_col]])
#   df <- df %>% filter(.data[[num_col]] >= 5)
# } else {
#   message("No 'Number'/'Count' column found — skipping <5 filter.")
# }
# 
# # Clean rows and keep only NRI/NTI
# df <- df %>%
#   mutate(
#     Category    = trimws(as.character(Category)),
#     Metrics     = trimws(as.character(Metrics)),
#     NRI_and_NTI = to_numeric(NRI_and_NTI),
#     .__err__    = to_numeric(.__err__)
#   ) %>%
#   filter(!is.na(Category), nzchar(Category), toupper(Category) != "NA",
#          !is.na(Metrics), Metrics %in% c("NRI","NTI"))
# 
# # ---------- ORDER depth categories ----------
# if (KEEP_ECOLOGICAL_ORDER && all(desired_levels %in% unique(df$Category))) {
#   ord_levels <- desired_levels
# } else {
#   ord_tbl <- df %>%
#     group_by(Category) %>%
#     summarise(order_mean = mean(ifelse(Metrics == ORDER_METRIC, NRI_and_NTI, NA), na.rm = TRUE),
#               .groups = "drop") %>%
#     arrange(if (ORDER_DIR == "asc") order_mean else desc(order_mean))
#   ord_levels <- ord_tbl$Category
# }
# 
# df$Category <- factor(df$Category, levels = ord_levels)
# df$Metrics  <- factor(df$Metrics,  levels = c("NRI","NTI"))
# df <- droplevels(df)
# 
# # ---------- paired stars ----------
# paired <- read.csv(paired_csv, check.names = FALSE)
# names(paired) <- trimws(names(paired))
# pcat <- intersect(c("Depth_Categories","Depth_Category","Depth","depth","Category","Salinity","salinity"), names(paired))
# if (length(pcat) == 0) stop("Paired table: no category column. Available: ", paste(names(paired), collapse = " | "))
# names(paired)[names(paired) == pcat[1]] <- "Category"
# 
# pcol <- intersect(c("p_value","p.value","p-value"), names(paired)); if (!length(pcol)) stop("Paired table: no p-value column.")
# paired <- paired %>%
#   mutate(
#     Category = trimws(as.character(.data$Category)),
#     stars    = p_to_stars(.data[[pcol[1]]])
#   ) %>%
#   filter(Category %in% ord_levels, nzchar(Category), toupper(Category) != "NA") %>%
#   mutate(Category = factor(Category, levels = ord_levels)) %>%
#   select(Category, stars) %>%
#   droplevels()
# 
# # ---------- dynamic positions & limits ----------
# y_min_data <- min(df$NRI_and_NTI - ifelse(df$NRI_and_NTI < 0, abs(df$.__err__), 0), na.rm = TRUE)
# y_max_data <- max(df$NRI_and_NTI + ifelse(df$NRI_and_NTI > 0, abs(df$.__err__), 0), na.rm = TRUE)
# rng_base   <- range(c(y_min_data, y_max_data, dash_at), na.rm = TRUE)
# 
# pad <- max(diff(rng_base) * 0.06, 0.15)  # gap bars↔labels
# tol <- max(diff(rng_base) * 0.02, 0.10)  # avoid dashed lines
# 
# cat_pos <- df %>%
#   group_by(Category) %>%
#   summarise(
#     top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#     bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(y_star = ifelse(top <= 0, bot - pad, top + pad))
# 
# for (d in dash_at) {
#   cat_pos$y_star <- ifelse(abs(cat_pos$y_star - d) < tol,
#                            ifelse(cat_pos$y_star >= d, cat_pos$y_star + pad*0.5, cat_pos$y_star - pad*0.5),
#                            cat_pos$y_star)
# }
# paired_ann <- left_join(paired, cat_pos, by = "Category")
# 
# # ---------- base plot ----------
# plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
#   geom_col(position = position_dodge(width = DODGE_W), width = BAR_W) +
#   geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
#                 position = position_dodge(width = DODGE_W), width = 0.22, linewidth = 0.7) +
#   geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
#   scale_fill_manual(values = fill_values, name = "Metrics") +
#   coord_flip() +
#   labs(title = "", x = "Depth", y = "NRI and NTI") +
#   theme_minimal(base_size = BASE) +
#   theme(
#     axis.title.x = element_text(size = TITLE, face = "bold"),
#     axis.title.y = element_text(size = TITLE, face = "bold"),
#     axis.text.x  = element_text(size = BASE),
#     axis.text.y  = element_text(size = BASE),
#     legend.title = element_text(size = TITLE, face = "bold"),
#     legend.text  = element_text(size = BASE),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_line(linetype = "dotted", linewidth = 0.3),
#     panel.grid.major.y = element_blank(),
#     plot.margin = margin(10, 14, 10, 10)
#   )
# 
# # ---------- (A) stars-only (with halo) ----------
# plot_stars <- plot_base +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_star, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = STAR_SIZE + 1.2, color = "white"
#   ) +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_star, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = STAR_SIZE, color = "black"
#   )
# 
# # ---------- (B) stars + CLD letters (colored by metric, halo) ----------
# # Read & sanitize independent tests
# indep <- read.csv(indep_csv, check.names = FALSE)
# nm <- trimws(names(indep)); nm[!nzchar(nm) | is.na(nm)] <- paste0("X", which(!nzchar(nm) | is.na(nm)))
# names(indep) <- nm
# nm <- names(indep)
# nm <- sub("^p[._-]?value$", "p_value", nm, ignore.case = TRUE)
# nm <- sub("^Group[ _]High$", "Group_High", nm)
# nm <- sub("^Group[ _]Low$",  "Group_Low",  nm)
# nm <- sub("^Metric[ ]*$",    "Metric",     nm)
# names(indep) <- nm
# drop_cols <- grep("^(Unnamed: ?\\d+|X\\d+)$", names(indep))
# if (length(drop_cols)) indep <- indep[, -drop_cols, drop = FALSE]
# 
# req <- c("Metric","Group_High","Group_Low","p_value")
# missing <- setdiff(req, names(indep))
# if (length(missing)) stop("Independent-tests table missing: ", paste(missing, collapse = ", "))
# 
# indep <- indep %>%
#   mutate(
#     Metric     = trimws(as.character(Metric)),
#     Group_High = trimws(as.character(Group_High)),
#     Group_Low  = trimws(as.character(Group_Low))
#   )
# 
# cats <- levels(df$Category)
# letters_df <- bind_rows(
#   { lt <- cld_letters_for_metric(indep, "NRI", cats, alpha = 0.05)
#     tibble(Metrics="NRI", Category=names(lt),
#            cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) },
#   { lt <- cld_letters_for_metric(indep, "NTI", cats, alpha = 0.05)
#     tibble(Metrics="NTI", Category=names(lt),
#            cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) }
# ) %>%
#   filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA") %>%
#   filter(Category %in% ord_levels) %>%
#   mutate(Category = factor(Category, levels = ord_levels)) %>%
#   droplevels()
# 
# plot_cld <- plot_stars  # includes halo stars already
# 
# if (nrow(letters_df) > 0) {
#   letter_pos <- df %>%
#     group_by(Category, Metrics) %>%
#     summarise(
#       y_top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#       y_bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     mutate(y_cld = ifelse(y_top <= 0, y_bot - pad*1.25, y_top + pad*1.25))
# 
#   # keep away from dashed lines & stars
#   for (d in dash_at) {
#     letter_pos$y_cld <- ifelse(abs(letter_pos$y_cld - d) < tol,
#                                ifelse(letter_pos$y_cld >= d, letter_pos$y_cld + pad*0.5, letter_pos$y_cld - pad*0.5),
#                                letter_pos$y_cld)
#   }
#   star_map <- cat_pos %>% select(Category, y_star)
#   letter_pos <- left_join(letter_pos, star_map, by = "Category") %>%
#     mutate(y_cld = ifelse(!is.na(y_star) & abs(y_cld - y_star) < (pad*0.6),
#                           ifelse(y_cld >= y_star, y_cld + pad*0.5, y_cld - pad*0.5),
#                           y_cld)) %>%
#     select(-y_star)
# 
#   letters_ann <- inner_join(letters_df, letter_pos, by = c("Category","Metrics")) %>%
#     filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA")
# 
#   plot_cld <- plot_cld +
#     geom_text(
#       data = letters_ann,
#       aes(x = Category, y = y_cld, label = cld, group = Metrics),
#       position = position_dodge(width = DODGE_W),
#       size = CLD_SIZE + 1.4, fontface = "bold", color = "white"
#     ) +
#     geom_text(
#       data = letters_ann,
#       aes(x = Category, y = y_cld, label = cld, group = Metrics, color = Metrics),
#       position = position_dodge(width = DODGE_W),
#       size = CLD_SIZE, fontface = "bold"
#     ) +
#     scale_color_manual(values = cld_cols, guide = "none")
# } else {
#   message("CLD letters are empty (no matching categories or p-values).")
# }
# 
# # ---------- y-limits (bars + stars + letters + dashed lines) ----------
# all_y <- c(df$NRI_and_NTI + df$.__err__, df$NRI_and_NTI - df$.__err__,
#            paired_ann$y_star,
#            if (exists("letters_ann")) letters_ann$y_cld else numeric(0),
#            dash_at)
# ymin <- floor((min(all_y, na.rm = TRUE) - pad*0.4) * 10) / 10
# ymax <- ceiling((max(all_y, na.rm = TRUE) + pad*0.4) * 10) / 10
# 
# plot_stars <- plot_stars + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))
# plot_cld   <- plot_cld   + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))
# 
# # ---------- save ----------
# file_stem_a <- file.path(out_dir, sprintf("depth_nri_nti_stars_orderBy%s_%s", ORDER_METRIC, ORDER_DIR))
# file_stem_b <- file.path(out_dir, sprintf("depth_nri_nti_cld_orderBy%s_%s",   ORDER_METRIC, ORDER_DIR))
# 
# ggsave(paste0(file_stem_a, ".svg"), plot_stars, width = 10, height = 10)
# ggsave(paste0(file_stem_a, ".png"), plot_stars, width = 10, height = 10, dpi = 300)
# ggsave(paste0(file_stem_b, ".svg"), plot_cld,  width = 10, height = 10)
# ggsave(paste0(file_stem_b, ".png"), plot_cld,  width = 10, height = 10, dpi = 300)
# 
# message("Saved:\n  ", paste0(file_stem_a, c(".svg",".png"), collapse = "\n  "),
#         "\n  ", paste0(file_stem_b, c(".svg",".png"), collapse = "\n  "))
# 
# # Use it later with: print(plot4)
# 
# # (optional) save to file when you want
# # ggsave(file.path(out_dir, "plot4.svg"), plot4, width = 10, height = 8)
# # ggsave(file.path(out_dir, "plot4.png"), plot4, width = 10, height = 8, dpi = 300)
# # (optional) save the R object:
# # saveRDS(plot4, file.path(out_dir, "plot4.rds"))

##############################################
# ==============================
# Depth NRI/NTI — stars + CLD (clean NA, colored letters, Number>=5)
# ==============================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(multcompView)
  library(tools)
})

# ---------- paths ----------
choose_first_existing <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
  ok[[1]]
}
data_csv   <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/depth_categories.csv",
  "/mnt/data/depth_categories.csv"
))
paired_csv <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/depth_paired_ttest_nri_nti_results.csv",
  "/mnt/data/depth_paired_ttest_nri_nti_results.csv"
))
indep_csv  <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/depth_independent_t_test_results_df.csv",
  "/mnt/data/depth_independent_t_test_results_df.csv"
))

out_dir <- dirname(data_csv)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- knobs ----------
BASE      <- 22
TITLE     <- 24
STAR_SIZE <- 7
CLD_SIZE  <- 8.2
BAR_W     <- 0.65
DODGE_W   <- 0.60
dash_at   <- c(-2, 2)

fill_values <- c("NRI"="#332288","NTI"="#DDCC77")  # bar fills
cld_cols    <- c(NRI = "#332288", NTI = "#A68000") # CLD colors (with halo)

# Ordering: by mean NRI or NTI, or keep ecological order if present
ORDER_METRIC         <- "NRI"   # or "NTI"
ORDER_DIR            <- "desc"  # "asc" or "desc"
KEEP_ECOLOGICAL_ORDER <- TRUE
desired_levels <- c("Shallow Zone","Epipelagic","Mesopelagic","Bathypelagic","Abyssopelagic","Hadalpelagic")

# ---------- helpers ----------
p_to_stars <- function(p) ifelse(is.na(p), "ns",
                                 ifelse(p <= 0.001, "***",
                                        ifelse(p <= 0.01,  "**",
                                               ifelse(p <= 0.05, "*", "ns"))))

safe_name <- function(x) { x <- trimws(as.character(x)); gsub("-", "–", x, fixed = TRUE) }

cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
  if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
  pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
  if (length(pcol) == 0) stop("Independent-tests table has no p-value column p_value/p.value/p-value.")
  pcol <- pcol[1]
  
  sub <- indep_df %>% filter(.data$Metric == metric,
                             .data$Group_High %in% cats,
                             .data$Group_Low  %in% cats)
  if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))
  
  n_hi <- safe_name(sub$Group_High)
  n_lo <- safe_name(sub$Group_Low)
  pair_names <- paste(pmin(n_hi, n_lo), pmax(n_hi, n_lo), sep = "-")
  pv <- sub[[pcol]]; names(pv) <- pair_names
  
  cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
  
  lookup <- setNames(cats, safe_name(cats))
  out <- setNames(rep("", length(cats)), cats)
  mapped <- intersect(names(cld), names(lookup))
  if (length(mapped)) out[lookup[mapped]] <- cld[mapped]
  out
}

to_numeric <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\u2212", "-", x, fixed = TRUE)
  x <- gsub("[^0-9eE+\\-\\.]", "", x)
  suppressWarnings(as.numeric(x))
}

# ---------- read + normalize (DEPTH) ----------
df <- read.csv(data_csv, check.names = FALSE)
names(df) <- trimws(names(df))

# Category
cat_col <- intersect(c("Depth_Categories","Depth_Category","Depth","depth","Category","Categories","salinity","Salinity"), names(df))
if (length(cat_col) == 0) stop("Data: couldn't find a depth/category column. Available: ", paste(names(df), collapse = " | "))
names(df)[names(df) == cat_col[1]] <- "Category"

# Metric, Value, Error
mcol <- intersect(c("Metrics","Metric"), names(df)); if (!length(mcol)) stop("No 'Metrics' column.")
names(df)[names(df) == mcol[1]] <- "Metrics"
vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df)); if (!length(vcol)) stop("No value column (NRI_and_NTI/Value).")
names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"
ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]

# ---- filter categories with Number < 5 (auto-detect count column) ----
num_col_idx <- grep("(?i)^(n$|n_|^n\\.|^n\\s*$|number|num|count|samples|n_samples|n_samp|n\\.samples|n\\.samp)$|(?i)number", names(df), perl = TRUE)
if (length(num_col_idx) > 0) {
  num_col <- names(df)[num_col_idx[1]]
  df[[num_col]] <- to_numeric(df[[num_col]])
  df <- df %>% filter(.data[[num_col]] >= 5)
} else {
  message("No 'Number'/'Count' column found — skipping <5 filter.")
}

# Clean rows and keep only NRI/NTI
df <- df %>%
  mutate(
    Category    = trimws(as.character(Category)),
    Metrics     = trimws(as.character(Metrics)),
    NRI_and_NTI = to_numeric(NRI_and_NTI),
    .__err__    = to_numeric(.__err__)
  ) %>%
  filter(!is.na(Category), nzchar(Category), toupper(Category) != "NA",
         !is.na(Metrics), Metrics %in% c("NRI","NTI"))

# ---------- ORDER depth categories ----------
if (KEEP_ECOLOGICAL_ORDER && all(desired_levels %in% unique(df$Category))) {
  ord_levels <- desired_levels
} else {
  ord_tbl <- df %>%
    group_by(Category) %>%
    summarise(order_mean = mean(ifelse(Metrics == ORDER_METRIC, NRI_and_NTI, NA), na.rm = TRUE),
              .groups = "drop") %>%
    arrange(if (ORDER_DIR == "asc") order_mean else desc(order_mean))
  ord_levels <- ord_tbl$Category
}

df$Category <- factor(df$Category, levels = ord_levels)
df$Metrics  <- factor(df$Metrics,  levels = c("NRI","NTI"))
df <- droplevels(df)

# ---------- paired stars ----------
paired <- read.csv(paired_csv, check.names = FALSE)
names(paired) <- trimws(names(paired))
pcat <- intersect(c("Depth_Categories","Depth_Category","Depth","depth","Category","Salinity","salinity"), names(paired))
if (length(pcat) == 0) stop("Paired table: no category column. Available: ", paste(names(paired), collapse = " | "))
names(paired)[names(paired) == pcat[1]] <- "Category"

pcol <- intersect(c("p_value","p.value","p-value"), names(paired)); if (!length(pcol)) stop("Paired table: no p-value column.")
paired <- paired %>%
  mutate(
    Category = trimws(as.character(.data$Category)),
    stars    = p_to_stars(.data[[pcol[1]]])
  ) %>%
  filter(Category %in% ord_levels, nzchar(Category), toupper(Category) != "NA") %>%
  mutate(Category = factor(Category, levels = ord_levels)) %>%
  select(Category, stars) %>%
  droplevels()

# ---------- dynamic positions & limits ----------
y_min_data <- min(df$NRI_and_NTI - ifelse(df$NRI_and_NTI < 0, abs(df$.__err__), 0), na.rm = TRUE)
y_max_data <- max(df$NRI_and_NTI + ifelse(df$NRI_and_NTI > 0, abs(df$.__err__), 0), na.rm = TRUE)
rng_base   <- range(c(y_min_data, y_max_data, dash_at), na.rm = TRUE)

pad <- max(diff(rng_base) * 0.06, 0.15)  # gap bars↔labels
tol <- max(diff(rng_base) * 0.02, 0.10)  # avoid dashed lines

cat_pos <- df %>%
  group_by(Category) %>%
  summarise(
    top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
    bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_star = ifelse(top <= 0, bot - pad, top + pad))

for (d in dash_at) {
  cat_pos$y_star <- ifelse(abs(cat_pos$y_star - d) < tol,
                           ifelse(cat_pos$y_star >= d, cat_pos$y_star + pad*0.5, cat_pos$y_star - pad*0.5),
                           cat_pos$y_star)
}
paired_ann <- left_join(paired, cat_pos, by = "Category")

# ---------- base plot ----------
plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
  geom_col(position = position_dodge(width = DODGE_W), width = BAR_W) +
  geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
                position = position_dodge(width = DODGE_W), width = 0.22, linewidth = 0.7) +
  geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
  scale_fill_manual(values = fill_values, name = "Metrics") +
  coord_flip() +
  labs(title = "", x = "Depth", y = "NRI and NTI") +
  theme_minimal(base_size = BASE) +
  theme(
    axis.title.x = element_text(size = TITLE, face = "bold"),
    axis.title.y = element_text(size = TITLE, face = "bold"),
    axis.text.x  = element_text(size = BASE),
    axis.text.y  = element_text(size = BASE),
    legend.title = element_text(size = TITLE, face = "bold"),
    legend.text  = element_text(size = BASE),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 14, 10, 10)
  )

# ---------- (A) stars-only (with halo) ----------
plot_stars <- plot_base +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE + 1.2, color = "white"
  ) +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE, color = "black"
  )

# ---------- (B) stars + CLD letters (colored by metric, halo) ----------
# Read & sanitize independent tests
indep <- read.csv(indep_csv, check.names = FALSE)
nm <- trimws(names(indep)); nm[!nzchar(nm) | is.na(nm)] <- paste0("X", which(!nzchar(nm) | is.na(nm)))
names(indep) <- nm
nm <- names(indep)
nm <- sub("^p[._-]?value$", "p_value", nm, ignore.case = TRUE)
nm <- sub("^Group[ _]High$", "Group_High", nm)
nm <- sub("^Group[ _]Low$",  "Group_Low",  nm)
nm <- sub("^Metric[ ]*$",    "Metric",     nm)
names(indep) <- nm
drop_cols <- grep("^(Unnamed: ?\\d+|X\\d+)$", names(indep))
if (length(drop_cols)) indep <- indep[, -drop_cols, drop = FALSE]

req <- c("Metric","Group_High","Group_Low","p_value")
missing <- setdiff(req, names(indep))
if (length(missing)) stop("Independent-tests table missing: ", paste(missing, collapse = ", "))

indep <- indep %>%
  mutate(
    Metric     = trimws(as.character(Metric)),
    Group_High = trimws(as.character(Group_High)),
    Group_Low  = trimws(as.character(Group_Low))
  )

cats <- levels(df$Category)
letters_df <- bind_rows(
  { lt <- cld_letters_for_metric(indep, "NRI", cats, alpha = 0.05)
  tibble(Metrics="NRI", Category=names(lt),
         cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) },
  { lt <- cld_letters_for_metric(indep, "NTI", cats, alpha = 0.05)
  tibble(Metrics="NTI", Category=names(lt),
         cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) }
) %>%
  filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA") %>%
  filter(Category %in% ord_levels) %>%
  mutate(Category = factor(Category, levels = ord_levels)) %>%
  droplevels()

plot_cld <- plot_stars  # includes halo stars already

if (nrow(letters_df) > 0) {
  letter_pos <- df %>%
    group_by(Category, Metrics) %>%
    summarise(
      y_top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
      y_bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_cld = ifelse(y_top <= 0, y_bot - pad*1.25, y_top + pad*1.25))
  
  # keep away from dashed lines & stars
  for (d in dash_at) {
    letter_pos$y_cld <- ifelse(abs(letter_pos$y_cld - d) < tol,
                               ifelse(letter_pos$y_cld >= d, letter_pos$y_cld + pad*0.5, letter_pos$y_cld - pad*0.5),
                               letter_pos$y_cld)
  }
  star_map <- cat_pos %>% select(Category, y_star)
  letter_pos <- left_join(letter_pos, star_map, by = "Category") %>%
    mutate(y_cld = ifelse(!is.na(y_star) & abs(y_cld - y_star) < (pad*0.6),
                          ifelse(y_cld >= y_star, y_cld + pad*0.5, y_cld - pad*0.5),
                          y_cld)) %>%
    select(-y_star)
  
  letters_ann <- inner_join(letters_df, letter_pos, by = c("Category","Metrics")) %>%
    filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA")
  
  plot_cld <- plot_cld +
    geom_text(
      data = letters_ann,
      aes(x = Category, y = y_cld, label = cld, group = Metrics),
      position = position_dodge(width = DODGE_W),
      size = CLD_SIZE + 1.4, fontface = "bold", color = "white"
    ) +
    geom_text(
      data = letters_ann,
      aes(x = Category, y = y_cld, label = cld, group = Metrics, color = Metrics),
      position = position_dodge(width = DODGE_W),
      size = CLD_SIZE, fontface = "bold"
    ) +
    scale_color_manual(values = cld_cols, guide = "none")
} else {
  message("CLD letters are empty (no matching categories or p-values).")
}

# ---------- y-limits (bars + stars + letters + dashed lines) ----------
all_y <- c(df$NRI_and_NTI + df$.__err__, df$NRI_and_NTI - df$.__err__,
           paired_ann$y_star,
           if (exists("letters_ann")) letters_ann$y_cld else numeric(0),
           dash_at)
ymin <- floor((min(all_y, na.rm = TRUE) - pad*0.4) * 10) / 10
ymax <- ceiling((max(all_y, na.rm = TRUE) + pad*0.4) * 10) / 10

plot_stars <- plot_stars + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))
plot_cld   <- plot_cld   + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))

# ---------- save ----------
file_stem_a <- file.path(out_dir, sprintf("depth_nri_nti_stars_orderBy%s_%s", ORDER_METRIC, ORDER_DIR))
file_stem_b <- file.path(out_dir, sprintf("depth_nri_nti_cld_orderBy%s_%s",   ORDER_METRIC, ORDER_DIR))

ggsave(paste0(file_stem_a, ".svg"), plot_stars, width = 10, height = 10)
ggsave(paste0(file_stem_a, ".png"), plot_stars, width = 10, height = 10, dpi = 300)
ggsave(paste0(file_stem_b, ".svg"), plot_cld,  width = 10, height = 10)
ggsave(paste0(file_stem_b, ".png"), plot_cld,  width = 10, height = 10, dpi = 300)

message("Saved:\n  ", paste0(file_stem_a, c(".svg",".png"), collapse = "\n  "),
        "\n  ", paste0(file_stem_b, c(".svg",".png"), collapse = "\n  "))

# save your final plot object
plot1 <- plot_cld   # or use plot_stars if you prefer that version

# (optional) show it
print(plot1)

# (optional) serialize to disk to reuse later
# saveRDS(plot_01, "plot_01.rds")
# plot_01 <- readRDS("plot_01.rds")


##############################################

# # ==============================
# # Salinity NRI/NTI — stars + CLD (robust numeric coercion)
# # ==============================
# 
# options(stringsAsFactors = FALSE)
# suppressPackageStartupMessages({
#   library(ggplot2)
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
#   library(tibble)
#   library(multcompView)
#   library(igraph)
#   library(tools)
# })
# 
# # ---------- paths ----------
# choose_first_existing <- function(paths) {
#   ok <- paths[file.exists(paths)]
#   if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
#   ok[[1]]
# }
# data_csv   <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/salinity_categories.csv",
#   "/mnt/data/salinity_categories.csv"
# ))
# paired_csv <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/salinity_paired_ttest_nri_nti_results.csv",
#   "/mnt/data/salinity_paired_ttest_nri_nti_results.csv"
# ))
# indep_csv  <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/salinity_independent_t_test_results_df.csv",
#   "/mnt/data/salinity_independent_t_test_results_df.csv"
# ))
# 
# out_dir <- dirname(data_csv)
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # ---------- helpers ----------
# p_to_stars <- function(p) {
#   ifelse(is.na(p), "ns",
#          ifelse(p <= 0.001, "***",
#                 ifelse(p <= 0.01,  "**",
#                        ifelse(p <= 0.05, "*", "ns"))))
# }
# 
# cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
#   if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
#   pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
#   if (length(pcol) == 0) stop("Independent-tests table has no p-value column named p_value/p.value/p-value.")
#   pcol <- pcol[1]
#   sub <- indep_df %>% filter(.data$Metric == metric,
#                              .data$Group_High %in% cats,
#                              .data$Group_Low  %in% cats)
#   if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))
#   pv <- sub[[pcol]]
#   names(pv) <- paste(pmin(sub$Group_High, sub$Group_Low),
#                      pmax(sub$Group_High, sub$Group_Low), sep = "-")
#   cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
#   out <- setNames(rep("", length(cats)), cats); out[names(cld)] <- cld; out
# }
# 
# # Coerce messy numerics safely (fix unicode minus, strip noise)
# to_numeric <- function(x) {
#   x <- trimws(as.character(x))
#   x <- gsub("\u2212", "-", x, fixed = TRUE)              # Unicode minus -> ASCII hyphen
#   x <- gsub("[^0-9eE+\\-\\.]", "", x)                    # keep digits/sign/decimal/exponent
#   suppressWarnings(as.numeric(x))
# }
# 
# # ---------- read + normalize (SALINITY main table) ----------
# df <- read.csv(data_csv, check.names = FALSE)
# names(df) <- trimws(names(df))
# 
# # Category
# cat_col <- intersect(c("Salinity","salinity","Category","Categories"), names(df))
# if (length(cat_col) == 0) stop("Data: couldn't find a category column. Available columns are: ", paste(names(df), collapse = ", "))
# names(df)[names(df) == cat_col[1]] <- "Category"
# 
# # Metric
# mcol <- intersect(c("Metrics","Metric"), names(df))
# if (length(mcol) == 0) stop("Data: no 'Metrics'/'Metric' column.")
# names(df)[names(df) == mcol[1]] <- "Metrics"
# 
# # Value
# vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df))
# if (length(vcol) == 0) stop("Data: no value column (expected NRI_and_NTI/Value).")
# names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"
# 
# # Error
# ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
# if (length(ecol) == 0) {
#   warning("No error column found; using zeros (no error bars).")
#   df$.__err__ <- 0
# } else {
#   df$.__err__ <- df[[ecol[1]]]
# }
# 
# # ---- CLEAN out footer/blank rows; coerce numerics AFTER filtering ----
# df <- df %>%
#   mutate(
#     Category = trimws(as.character(Category)),
#     Metrics  = trimws(as.character(Metrics))
#   ) %>%
#   filter(!is.na(Category), !is.na(Metrics), Metrics %in% c("NRI","NTI")) %>%
#   mutate(
#     NRI_and_NTI = to_numeric(NRI_and_NTI),
#     .__err__    = to_numeric(.__err__)
#   )
# 
# # Order salinity levels
# desired_levels <- c("Very High","High","Moderate","Low")
# if (all(desired_levels %in% unique(df$Category))) {
#   df$Category <- factor(df$Category, levels = desired_levels)
# } else {
#   df$Category <- factor(df$Category, levels = unique(df$Category))
# }
# df$Metrics  <- factor(df$Metrics, levels = c("NRI","NTI"))
# cat_levels  <- levels(df$Category)
# 
# # ---------- paired stars ----------
# paired <- read.csv(paired_csv, check.names = FALSE)
# names(paired) <- trimws(names(paired))
# pcat <- intersect(c("Salinity","salinity","Category"), names(paired))
# if (length(pcat) == 0) stop("Paired table: no category column. Available: ", paste(names(paired), collapse = ", "))
# names(paired)[names(paired) == pcat[1]] <- "Category"
# 
# pcol <- intersect(c("p_value","p.value","p-value"), names(paired))
# if (length(pcol) == 0) stop("Paired table: no p-value column named p_value/p.value/p-value.")
# 
# paired <- paired %>%
#   mutate(Category = factor(trimws(as.character(.data$Category)), levels = cat_levels),
#          stars    = p_to_stars(.data[[pcol[1]]])) %>%
#   select(Category, stars)
# 
# # Positions for stars (change 0.15 to move stars)
# cat_pos <- df %>%
#   group_by(Category) %>%
#   summarise(
#     top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#     bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(y_pos = ifelse(top <= 0, bot - 0.50, top + 0.50))
# 
# paired_ann <- left_join(paired, cat_pos, by = "Category")
# 
# # ---------- base plot ----------
# dodge_w <- 0.6
# fill_values <- c("NRI"="#332288","NTI"="#DDCC77")
# dash_at <- c(-2,2)
# 
# plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
#   geom_bar(stat = "identity", position = position_dodge(width = dodge_w), width = 0.6) +
#   geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
#                 position = position_dodge(width = dodge_w), width = 0.2) +
#   geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
#   scale_fill_manual(values = fill_values) +
#   coord_flip() +
#   theme_minimal() +
#   labs(title = "", x = "Salinity", y = "NRI and NTI", fill = "Metrics")
# 
# # ---------- (A) stars-only ----------
# plot_stars <- plot_base +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = 3.8
#   )
# 
# # ---------- (B) stars + CLD letters ----------
# indep <- read.csv(indep_csv, check.names = FALSE)
# nm <- trimws(names(indep))
# empty_idx <- which(!nzchar(nm) | is.na(nm))
# if (length(empty_idx) > 0) nm[empty_idx] <- paste0("X", empty_idx)
# names(indep) <- nm
# nm <- names(indep)
# nm <- sub("^p[._-]?value$", "p_value", nm, ignore.case = TRUE)
# nm <- sub("^Group[ _]High$", "Group_High", nm)
# nm <- sub("^Group[ _]Low$",  "Group_Low",  nm)
# nm <- sub("^Metric[ ]*$",    "Metric",     nm)
# names(indep) <- nm
# drop_cols <- grep("^(Unnamed: ?\\d+|X\\d+)$", names(indep))
# if (length(drop_cols)) indep <- indep[, -drop_cols, drop = FALSE]
# req <- c("Metric","Group_High","Group_Low","p_value")
# missing <- setdiff(req, names(indep))
# if (length(missing)) stop("Independent-tests table is missing column(s): ", paste(missing, collapse = ", "))
# indep <- indep %>%
#   mutate(
#     Metric     = trimws(as.character(Metric)),
#     Group_High = trimws(as.character(Group_High)),
#     Group_Low  = trimws(as.character(Group_Low))
#   )
# 
# cats <- levels(df$Category)
# letters_list <- lapply(c("NRI","NTI"), function(m) {
#   lt  <- cld_letters_for_metric(indep, metric = m, cats = cats, alpha = 0.05)
#   tibble(
#     Metrics  = m,
#     Category = names(lt),
#     cld      = gsub("(?<=.)(?=.)", " ", unname(lt), perl = TRUE)
#   )
# })
# letters_df <- bind_rows(letters_list)
# 
# plot_cld <- plot_base +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = 3.8
#   )
# 
# if (nrow(letters_df) > 0) {
#   letter_pos <- df %>%
#     group_by(Category, Metrics) %>%
#     summarise(
#       y_top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#       y_bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     mutate(y_pos = ifelse(y_top <= 0, y_bot - 0.90, y_top + 0.90))
#   
#   letters_ann <- left_join(letters_df, letter_pos, by = c("Category","Metrics"))
#   
#   plot_cld <- plot_cld +
#     geom_text(
#       data = letters_ann,
#       aes(x = Category, y = y_pos, label = cld, group = Metrics),
#       position = position_dodge(width = dodge_w),
#       size = 3.4
#     )
# } else {
#   message("CLD letters computed empty (no matching categories or missing p-values).")
# }
# 
# # ---------- save ----------
# file_stem_a <- file.path(out_dir, "salinity_nri_nti_stars")
# file_stem_b <- file.path(out_dir, "salinity_nri_nti_cld")
# 
# ggsave(paste0(file_stem_a, ".svg"), plot_stars, width = 8, height = 8)
# ggsave(paste0(file_stem_a, ".png"), plot_stars, width = 8, height = 8, dpi = 300)
# ggsave(paste0(file_stem_b, ".svg"), plot_cld,  width = 8, height = 8)
# ggsave(paste0(file_stem_b, ".png"), plot_cld,  width = 8, height = 8, dpi = 300)
# 
# message("Saved:\n  ", paste0(file_stem_a, c(".svg",".png"), collapse = "\n  "),
#         "\n  ", paste0(file_stem_b, c(".svg",".png"), collapse = "\n  "))
# ######################
# 
# # (optional) small cleanup to avoid "Removed X rows" warnings
# paired_ann2  <- tryCatch(dplyr::filter(paired_ann, !is.na(y_pos), !is.na(stars)), error = function(e) paired_ann)
# letters_ann2 <- if (exists("letters_ann") && is.data.frame(letters_ann)) {
#   dplyr::filter(letters_ann, !is.na(y_pos), !is.na(cld), cld != "")
# } else NULL
# 
# plot_2 <- plot_base +
#   # paired t-test stars
#   geom_text(
#     data = paired_ann2,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold",
#     size = 3.8
#   ) +
#   # CLD letters (only if available)
#   {
#     if (!is.null(letters_ann2) && nrow(letters_ann2) > 0) {
#       geom_text(
#         data = letters_ann2,
#         aes(x = Category, y = y_pos, label = cld, group = Metrics),
#         position = position_dodge(width = dodge_w),
#         size = 3.4
#       )
#     } else {
#       NULL
#     }
#   }
# 
# # Use it later with: print(plot4)
# 
# # (optional) save to file when you want
# # ggsave(file.path(out_dir, "plot4.svg"), plot4, width = 10, height = 8)
# # ggsave(file.path(out_dir, "plot4.png"), plot4, width = 10, height = 8, dpi = 300)
# # (optional) save the R object:
# # saveRDS(plot4, file.path(out_dir, "plot4.rds"))
#####################################

# ==============================
# Salinity NRI/NTI — stars + CLD (clean NA, colored CLD, Number>=5, ordering)
# ==============================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(multcompView)
  library(tools)
})

# ---------- paths ----------
choose_first_existing <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
  ok[[1]]
}
data_csv   <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/salinity_categories.csv",
  "/mnt/data/salinity_categories.csv"
))
paired_csv <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/salinity_paired_ttest_nri_nti_results.csv",
  "/mnt/data/salinity_paired_ttest_nri_nti_results.csv"
))
indep_csv  <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/salinity_independent_t_test_results_df.csv",
  "/mnt/data/salinity_independent_t_test_results_df.csv"
))

out_dir <- dirname(data_csv)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- knobs ----------
BASE      <- 22
TITLE     <- 24
STAR_SIZE <- 7
CLD_SIZE  <- 8.2
BAR_W     <- 0.65
DODGE_W   <- 0.60
dash_at   <- c(-2, 2)

fill_values <- c("NRI"="#332288","NTI"="#DDCC77")  # bar fills
cld_cols    <- c(NRI = "#332288", NTI = "#A68000") # CLD letters (colored by metric)

# Ordering: keep defined order if present, otherwise order by mean metric
KEEP_DEFINED_ORDER <- TRUE
defined_levels <- c("Very High","High","Moderate","Low")
ORDER_METRIC   <- "NRI"   # or "NTI"
ORDER_DIR      <- "desc"  # "asc" or "desc"

# ---------- helpers ----------
p_to_stars <- function(p) ifelse(is.na(p), "ns",
                                 ifelse(p <= 0.001, "***",
                                        ifelse(p <= 0.01,  "**",
                                               ifelse(p <= 0.05, "*", "ns"))))

# Coerce messy numerics safely
to_numeric <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\u2212", "-", x, fixed = TRUE)   # Unicode minus → ASCII
  x <- gsub("[^0-9eE+\\-\\.]", "", x)
  suppressWarnings(as.numeric(x))
}

# hyphen-safe names for multcompView
safe_name <- function(x) { x <- trimws(as.character(x)); gsub("-", "–", x, fixed = TRUE) }

# CLD letters from independent tests (metric-specific)
cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
  if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
  pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
  if (!length(pcol)) stop("Independent-tests table has no p-value column p_value/p.value/p-value.")
  pcol <- pcol[1]
  
  sub <- indep_df %>% filter(.data$Metric == metric,
                             .data$Group_High %in% cats,
                             .data$Group_Low  %in% cats)
  if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))
  
  n_hi <- safe_name(sub$Group_High)
  n_lo <- safe_name(sub$Group_Low)
  pair_names <- paste(pmin(n_hi, n_lo), pmax(n_hi, n_lo), sep = "-")
  pv <- sub[[pcol]]; names(pv) <- pair_names
  
  cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
  
  lookup <- setNames(cats, safe_name(cats))
  out <- setNames(rep("", length(cats)), cats)
  mapped <- intersect(names(cld), names(lookup))
  if (length(mapped)) out[lookup[mapped]] <- cld[mapped]
  out
}

# ---------- read + normalize ----------
df <- read.csv(data_csv, check.names = FALSE)
names(df) <- trimws(names(df))

# Category
cat_col <- intersect(c("Salinity","salinity","Category","Categories"), names(df))
if (!length(cat_col)) stop("Data: couldn't find a category column. Available: ", paste(names(df), collapse = " | "))
names(df)[names(df) == cat_col[1]] <- "Category"

# Metric / Value / Error
mcol <- intersect(c("Metrics","Metric"), names(df)); if (!length(mcol)) stop("No 'Metrics' column.")
names(df)[names(df) == mcol[1]] <- "Metrics"
vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df)); if (!length(vcol)) stop("No value column (NRI_and_NTI/Value).")
names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"
ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]

# Filter categories with Number < 5 (auto-detect count column)
num_col_idx <- grep("(?i)^(n$|n_|^n\\.|^n\\s*$|number|num|count|samples|n_samples|n_samp|n\\.samples|n\\.samp)$|(?i)number", names(df), perl = TRUE)
if (length(num_col_idx) > 0) {
  num_col <- names(df)[num_col_idx[1]]
  df[[num_col]] <- to_numeric(df[[num_col]])
  df <- df %>% filter(.data[[num_col]] >= 5)
} else {
  message("No 'Number'/'Count' column found — skipping <5 filter.")
}

# Clean rows and coerce numerics; keep only NRI/NTI; drop NA/blank/"NA"
df <- df %>%
  mutate(
    Category    = trimws(as.character(Category)),
    Metrics     = trimws(as.character(Metrics)),
    NRI_and_NTI = to_numeric(NRI_and_NTI),
    .__err__    = to_numeric(.__err__)
  ) %>%
  filter(!is.na(Category), nzchar(Category), toupper(Category) != "NA",
         !is.na(Metrics), Metrics %in% c("NRI","NTI"))

# ---------- order salinity categories ----------
if (KEEP_DEFINED_ORDER && all(defined_levels %in% unique(df$Category))) {
  ord_levels <- defined_levels
} else {
  ord_tbl <- df %>%
    group_by(Category) %>%
    summarise(order_mean = mean(ifelse(Metrics == ORDER_METRIC, NRI_and_NTI, NA), na.rm = TRUE),
              .groups = "drop") %>%
    arrange(if (ORDER_DIR == "asc") order_mean else desc(order_mean))
  ord_levels <- ord_tbl$Category
}
df$Category <- factor(df$Category, levels = ord_levels)
df$Metrics  <- factor(df$Metrics,  levels = c("NRI","NTI"))
df <- droplevels(df)

# ---------- paired stars ----------
paired <- read.csv(paired_csv, check.names = FALSE)
names(paired) <- trimws(names(paired))
pcat <- intersect(c("Salinity","salinity","Category"), names(paired))
if (!length(pcat)) stop("Paired table: no category column. Available: ", paste(names(paired), collapse = " | "))
names(paired)[names(paired) == pcat[1]] <- "Category"

pcol <- intersect(c("p_value","p.value","p-value"), names(paired)); if (!length(pcol)) stop("Paired table: no p-value column.")
paired <- paired %>%
  mutate(
    Category = trimws(as.character(.data$Category)),
    stars    = p_to_stars(.data[[pcol[1]]])
  ) %>%
  filter(Category %in% ord_levels, nzchar(Category), toupper(Category) != "NA") %>%
  mutate(Category = factor(Category, levels = ord_levels)) %>%
  select(Category, stars) %>%
  droplevels()

# ---------- dynamic positions & limits ----------
y_min_data <- min(df$NRI_and_NTI - ifelse(df$NRI_and_NTI < 0, abs(df$.__err__), 0), na.rm = TRUE)
y_max_data <- max(df$NRI_and_NTI + ifelse(df$NRI_and_NTI > 0, abs(df$.__err__), 0), na.rm = TRUE)
rng_base   <- range(c(y_min_data, y_max_data, dash_at), na.rm = TRUE)

pad <- max(diff(rng_base) * 0.06, 0.15)  # gap bars↔labels
tol <- max(diff(rng_base) * 0.02, 0.10)  # avoid dashed lines

cat_pos <- df %>%
  group_by(Category) %>%
  summarise(
    top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
    bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_star = ifelse(top <= 0, bot - pad, top + pad))

# nudge stars away from dashed lines
for (d in dash_at) {
  cat_pos$y_star <- ifelse(abs(cat_pos$y_star - d) < tol,
                           ifelse(cat_pos$y_star >= d, cat_pos$y_star + pad*0.5, cat_pos$y_star - pad*0.5),
                           cat_pos$y_star)
}
paired_ann <- left_join(paired, cat_pos, by = "Category")

# ---------- base plot ----------
plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
  geom_col(position = position_dodge(width = DODGE_W), width = BAR_W) +
  geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
                position = position_dodge(width = DODGE_W), width = 0.22, linewidth = 0.7) +
  geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
  scale_fill_manual(values = fill_values, name = "Metrics") +
  coord_flip() +
  labs(title = "", x = "Salinity", y = "NRI and NTI") +
  theme_minimal(base_size = BASE) +
  theme(
    axis.title.x = element_text(size = TITLE, face = "bold"),
    axis.title.y = element_text(size = TITLE, face = "bold"),
    axis.text.x  = element_text(size = BASE),
    axis.text.y  = element_text(size = BASE),
    legend.title = element_text(size = TITLE, face = "bold"),
    legend.text  = element_text(size = BASE),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 14, 10, 10)
  )

# ---------- (A) stars-only (with halo) ----------
plot_stars <- plot_base +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE + 1.2, color = "white"
  ) +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE, color = "black"
  )

# ---------- (B) stars + CLD letters (colored, with halo) ----------
indep <- read.csv(indep_csv, check.names = FALSE)
nm <- trimws(names(indep)); nm[!nzchar(nm) | is.na(nm)] <- paste0("X", which(!nzchar(nm) | is.na(nm)))
names(indep) <- nm
nm <- names(indep)
nm <- sub("^p[._-]?value$", "p_value", nm, ignore.case = TRUE)
nm <- sub("^Group[ _]High$", "Group_High", nm)
nm <- sub("^Group[ _]Low$",  "Group_Low",  nm)
nm <- sub("^Metric[ ]*$",    "Metric",     nm)
names(indep) <- nm
drop_cols <- grep("^(Unnamed: ?\\d+|X\\d+)$", names(indep))
if (length(drop_cols)) indep <- indep[, -drop_cols, drop = FALSE]

req <- c("Metric","Group_High","Group_Low","p_value")
missing <- setdiff(req, names(indep))
if (length(missing)) stop("Independent-tests table missing: ", paste(missing, collapse = ", "))

indep <- indep %>%
  mutate(
    Metric     = trimws(as.character(Metric)),
    Group_High = trimws(as.character(Group_High)),
    Group_Low  = trimws(as.character(Group_Low))
  )

cats <- levels(df$Category)
letters_df <- bind_rows(
  { lt <- cld_letters_for_metric(indep, "NRI", cats, alpha = 0.05)
  tibble(Metrics="NRI", Category=names(lt),
         cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) },
  { lt <- cld_letters_for_metric(indep, "NTI", cats, alpha = 0.05)
  tibble(Metrics="NTI", Category=names(lt),
         cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) }
) %>%
  filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA") %>%
  filter(Category %in% ord_levels) %>%
  mutate(Category = factor(Category, levels = ord_levels)) %>%
  droplevels()

plot_cld <- plot_stars  # start from stars (halo stars already drawn)

if (nrow(letters_df) > 0) {
  letter_pos <- df %>%
    group_by(Category, Metrics) %>%
    summarise(
      y_top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
      y_bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_cld = ifelse(y_top <= 0, y_bot - pad*1.25, y_top + pad*1.25))
  
  # keep away from dashed lines & stars
  for (d in dash_at) {
    letter_pos$y_cld <- ifelse(abs(letter_pos$y_cld - d) < tol,
                               ifelse(letter_pos$y_cld >= d, letter_pos$y_cld + pad*0.5, letter_pos$y_cld - pad*0.5),
                               letter_pos$y_cld)
  }
  star_map <- cat_pos %>% select(Category, y_star)
  letter_pos <- left_join(letter_pos, star_map, by = "Category") %>%
    mutate(y_cld = ifelse(!is.na(y_star) & abs(y_cld - y_star) < (pad*0.6),
                          ifelse(y_cld >= y_star, y_cld + pad*0.5, y_cld - pad*0.5),
                          y_cld)) %>%
    select(-y_star)
  
  letters_ann <- inner_join(letters_df, letter_pos, by = c("Category","Metrics")) %>%
    filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA")
  
  plot_cld <- plot_cld +
    geom_text(
      data = letters_ann,
      aes(x = Category, y = y_cld, label = cld, group = Metrics),
      position = position_dodge(width = DODGE_W),
      size = CLD_SIZE + 1.4, fontface = "bold", color = "white"
    ) +
    geom_text(
      data = letters_ann,
      aes(x = Category, y = y_cld, label = cld, group = Metrics, color = Metrics),
      position = position_dodge(width = DODGE_W),
      size = CLD_SIZE, fontface = "bold"
    ) +
    scale_color_manual(values = cld_cols, guide = "none")
} else {
  message("CLD letters are empty (no matching categories or p-values).")
}

# ---------- y-limits (bars + stars + letters + dashed lines) ----------
all_y <- c(df$NRI_and_NTI + df$.__err__, df$NRI_and_NTI - df$.__err__,
           paired_ann$y_star,
           if (exists("letters_ann")) letters_ann$y_cld else numeric(0),
           dash_at)
ymin <- floor((min(all_y, na.rm = TRUE) - pad*0.4) * 10) / 10
ymax <- ceiling((max(all_y, na.rm = TRUE) + pad*0.4) * 10) / 10

plot_stars <- plot_stars + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))
plot_cld   <- plot_cld   + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))

# ---------- save ----------
file_stem_a <- file.path(out_dir, sprintf("salinity_nri_nti_stars_orderBy%s_%s", ORDER_METRIC, ORDER_DIR))
file_stem_b <- file.path(out_dir, sprintf("salinity_nri_nti_cld_orderBy%s_%s",   ORDER_METRIC, ORDER_DIR))

ggsave(paste0(file_stem_a, ".svg"), plot_stars, width = 10, height = 10)
ggsave(paste0(file_stem_a, ".png"), plot_stars, width = 10, height = 10, dpi = 300)
ggsave(paste0(file_stem_b, ".svg"), plot_cld,  width = 10, height = 10)
ggsave(paste0(file_stem_b, ".png"), plot_cld,  width = 10, height = 10, dpi = 300)

message("Saved:\n  ", paste0(file_stem_a, c(".svg",".png"), collapse = "\n  "),
        "\n  ", paste0(file_stem_b, c(".svg",".png"), collapse = "\n  "))

# (optional) keep a plot object
plot2 <- plot_cld



#####################################

# # ==============================
# # Ocean NRI/NTI — stars + CLD (hyphen-safe CLD + no extra categories)
# # ==============================
# 
# options(stringsAsFactors = FALSE)
# suppressPackageStartupMessages({
#   library(ggplot2)
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
#   library(tibble)
#   library(multcompView)
#   library(igraph)
#   library(tools)
# })
# 
# # ---------- paths ----------
# choose_first_existing <- function(paths) {
#   ok <- paths[file.exists(paths)]
#   if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
#   ok[[1]]
# }
# data_csv   <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/ocean.csv",
#   "/mnt/data/ocean.csv"
# ))
# paired_csv <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/ocean_paired_ttest_nri_nti_results.csv",
#   "/mnt/data/ocean_paired_ttest_nri_nti_results.csv"
# ))
# indep_csv  <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/ocean_independent_t_test_results_df.csv",
#   "/mnt/data/ocean_independent_t_test_results_df.csv"
# ))
# 
# out_dir <- dirname(data_csv)
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # ---------- helpers ----------
# p_to_stars <- function(p) {
#   ifelse(is.na(p), "ns",
#          ifelse(p <= 0.001, "***",
#                 ifelse(p <= 0.01,  "**",
#                        ifelse(p <= 0.05, "*", "ns"))))
# }
# 
# # Replace internal ASCII hyphens in labels by en dash “–” so there's exactly
# # ONE ASCII hyphen left (the pair separator) in names fed to multcompView.
# safe_name <- function(x) {
#   x <- trimws(as.character(x))
#   gsub("-", "–", x, fixed = TRUE)
# }
# 
# cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
#   if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
#   pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
#   if (length(pcol) == 0) stop("Independent-test table has no p-value column named p_value/p.value/p-value.")
#   pcol <- pcol[1]
#   
#   sub <- indep_df %>%
#     filter(.data$Metric == metric,
#            .data$Group_High %in% cats,
#            .data$Group_Low  %in% cats)
#   
#   if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))
#   
#   # Build names with exactly ONE ASCII hyphen as separator
#   n_hi <- safe_name(sub$Group_High)
#   n_lo <- safe_name(sub$Group_Low)
#   pair_names <- paste(pmin(n_hi, n_lo), pmax(n_hi, n_lo), sep = "-")
#   
#   pv <- sub[[pcol]]
#   names(pv) <- pair_names
#   
#   cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
#   
#   # Map back to original cats (not the safe-name strings)
#   # Build a lookup: safe -> original
#   lookup <- setNames(cats, safe_name(cats))
#   out <- setNames(rep("", length(cats)), cats)
#   names(out) <- cats
#   # cld is named by safe names; remap names to original where possible
#   mapped <- intersect(names(cld), names(lookup))
#   if (length(mapped)) {
#     out[lookup[mapped]] <- cld[mapped]
#   }
#   out
# }
# 
# # ---------- read + normalize (OCEAN main table) ----------
# df <- read.csv(data_csv, check.names = FALSE)
# names(df) <- trimws(names(df))
# 
# # Category (your file uses 'salinity')
# cat_col <- intersect(c("salinity","Salinity","Category","Categories"), names(df))
# if (length(cat_col) == 0) stop("Data: couldn't find a category column. Available columns are: ", paste(names(df), collapse = ", "))
# names(df)[names(df) == cat_col[1]] <- "Category"
# 
# # Metric
# mcol <- intersect(c("Metrics","Metric"), names(df))
# if (length(mcol) == 0) stop("Data: no 'Metrics'/'Metric' column.")
# names(df)[names(df) == mcol[1]] <- "Metrics"
# 
# # Value
# vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df))
# if (length(vcol) == 0) stop("Data: no value column (expected NRI_and_NTI/Value).")
# names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"
# 
# # Error
# ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
# df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]
# 
# # Keep order of appearance; drop unused levels to avoid “extra” blank bars
# df$Category <- factor(df$Category, levels = unique(df$Category))
# df$Metrics  <- factor(df$Metrics,  levels = c("NRI","NTI"))
# df <- droplevels(df)
# cat_levels  <- levels(df$Category)
# 
# # ---------- paired stars ----------
# paired <- read.csv(paired_csv, check.names = FALSE)
# names(paired) <- trimws(names(paired))
# pcat <- intersect(c("Salinity","salinity","Category"), names(paired))
# if (length(pcat) == 0) stop("Paired table: no category column. Available: ", paste(names(paired), collapse = ", "))
# names(paired)[names(paired) == pcat[1]] <- "Category"
# 
# pcol <- intersect(c("p_value","p.value","p-value"), names(paired))
# if (length(pcol) == 0) stop("Paired table: no p-value column named p_value/p.value/p-value.")
# 
# paired <- paired %>%
#   mutate(Category = factor(trimws(as.character(.data$Category)), levels = cat_levels),
#          stars    = p_to_stars(.data[[pcol[1]]])) %>%
#   select(Category, stars)
# 
# # Positions for stars (tweak 0.15 to move stars)
# cat_pos <- df %>%
#   group_by(Category) %>%
#   summarise(
#     top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#     bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(y_pos = ifelse(top <= 0, bot - 0.50, top + 0.50))
# 
# paired_ann <- left_join(paired, cat_pos, by = "Category")
# 
# # ---------- base plot ----------
# dodge_w <- 0.6
# fill_values <- c("NRI"="#332288","NTI"="#DDCC77")
# dash_at <- c(-2,2)
# 
# plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
#   geom_bar(stat = "identity", position = position_dodge(width = dodge_w), width = 0.6) +
#   geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
#                 position = position_dodge(width = dodge_w), width = 0.2) +
#   geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
#   scale_fill_manual(values = fill_values) +
#   coord_flip() +
#   theme_minimal() +
#   labs(title = "", x = "Ocean / Sea", y = "NRI and NTI", fill = "Metrics")
# 
# # ---------- (A) stars-only ----------
# plot_stars <- plot_base +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = 3.8
#   )
# 
# # ---------- (B) stars + CLD letters (hyphen-safe) ----------
# indep <- read.csv(indep_csv, check.names = FALSE)
# 
# # Fix empty/NA names first, then normalize
# nm <- trimws(names(indep))
# empty_idx <- which(!nzchar(nm) | is.na(nm))
# if (length(empty_idx) > 0) nm[empty_idx] <- paste0("X", empty_idx)
# names(indep) <- nm
# nm <- names(indep)
# nm <- sub("^p[._-]?value$", "p_value", nm, ignore.case = TRUE)
# nm <- sub("^Group[ _]High$", "Group_High", nm)
# nm <- sub("^Group[ _]Low$",  "Group_Low",  nm)
# nm <- sub("^Metric[ ]*$",    "Metric",     nm)
# names(indep) <- nm
# 
# # Drop index-like columns safely
# drop_cols <- grep("^(Unnamed: ?\\d+|X\\d+)$", names(indep))
# if (length(drop_cols)) indep <- indep[, -drop_cols, drop = FALSE]
# 
# # Check required columns
# req <- c("Metric","Group_High","Group_Low","p_value")
# missing <- setdiff(req, names(indep))
# if (length(missing)) stop("Independent-tests table is missing column(s): ", paste(missing, collapse = ", "))
# 
# # Tidy values to match Category factor exactly
# indep <- indep %>%
#   mutate(
#     Metric     = trimws(as.character(Metric)),
#     Group_High = trimws(as.character(Group_High)),
#     Group_Low  = trimws(as.character(Group_Low))
#   )
# 
# # Build CLD labels (safe to hyphens)
# cats <- levels(df$Category)
# letters_list <- lapply(c("NRI","NTI"), function(m) {
#   lt  <- cld_letters_for_metric(indep, metric = m, cats = cats, alpha = 0.05)
#   tibble(
#     Metrics  = m,
#     Category = names(lt),
#     cld      = gsub("(?<=.)(?=.)", " ", unname(lt), perl = TRUE)  # "AB" -> "A B"
#   )
# })
# letters_df <- bind_rows(letters_list)
# 
# plot_cld <- plot_base +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = 3.8
#   )
# 
# if (nrow(letters_df) > 0) {
#   letter_pos <- df %>%
#     group_by(Category, Metrics) %>%
#     summarise(
#       y_top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#       y_bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     mutate(y_pos = ifelse(y_top <= 0, y_bot - 1.10, y_top + 1.10))
#   
#   # inner_join to drop rows that lack a y position (prevents "Removed X rows..." warning)
#   letters_ann <- inner_join(letters_df, letter_pos, by = c("Category","Metrics")) %>%
#     filter(!is.na(cld), nzchar(cld))
#   
#   plot_cld <- plot_cld +
#     geom_text(
#       data = letters_ann,
#       aes(x = Category, y = y_pos, label = cld, group = Metrics),
#       position = position_dodge(width = dodge_w),
#       size = 3.4
#     )
# } else {
#   message("CLD letters computed empty (no matching categories or missing p-values).")
# }
# 
# # ---------- save ----------
# file_stem_a <- file.path(out_dir, "ocean_nri_nti_stars")
# file_stem_b <- file.path(out_dir, "ocean_nri_nti_cld")
# 
# ggsave(paste0(file_stem_a, ".svg"), plot_stars, width = 10, height = 8)
# ggsave(paste0(file_stem_a, ".png"), plot_stars, width = 10, height = 8, dpi = 300)
# ggsave(paste0(file_stem_b, ".svg"), plot_cld,  width = 10, height = 8)
# ggsave(paste0(file_stem_b, ".png"), plot_cld,  width = 10, height = 8, dpi = 300)
# 
# message("Saved:\n  ", paste0(file_stem_a, c(".svg",".png"), collapse = "\n  "),
#         "\n  ", paste0(file_stem_b, c(".svg",".png"), collapse = "\n  "))
# ####################################
# 
# # --- EXTENSION: single figure with BOTH stars and CLD, stored as `plot4` ---
# 
# # (optional) small cleanup to avoid "Removed X rows" warnings
# paired_ann2  <- tryCatch(dplyr::filter(paired_ann, !is.na(y_pos), !is.na(stars)), error = function(e) paired_ann)
# letters_ann2 <- if (exists("letters_ann") && is.data.frame(letters_ann)) {
#   dplyr::filter(letters_ann, !is.na(y_pos), !is.na(cld), cld != "")
# } else NULL
# 
# plot_5 <- plot_base +
#   # paired t-test stars
#   geom_text(
#     data = paired_ann2,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold",
#     size = 3.8
#   ) +
#   # CLD letters (only if available)
#   {
#     if (!is.null(letters_ann2) && nrow(letters_ann2) > 0) {
#       geom_text(
#         data = letters_ann2,
#         aes(x = Category, y = y_pos, label = cld, group = Metrics),
#         position = position_dodge(width = dodge_w),
#         size = 3.4
#       )
#     } else {
#       NULL
#     }
#   }
# 
# # Use it later with: print(plot4)
# 
# # (optional) save to file when you want
# # ggsave(file.path(out_dir, "plot4.svg"), plot4, width = 10, height = 8)
# # ggsave(file.path(out_dir, "plot4.png"), plot4, width = 10, height = 8, dpi = 300)
# # (optional) save the R object:
# # saveRDS(plot4, file.path(out_dir, "plot4.rds"))
# 
# #dev.off()
# ######################################

# ==============================
# Ocean NRI/NTI — stars + CLD (clean NA, Number>=5, colored CLD, hyphen-safe)
# ==============================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(multcompView)
  library(tools)
})

# ---------- paths ----------
choose_first_existing <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
  ok[[1]]
}
data_csv   <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/ocean.csv",
  "/mnt/data/ocean.csv"
))
paired_csv <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/ocean_paired_ttest_nri_nti_results.csv",
  "/mnt/data/ocean_paired_ttest_nri_nti_results.csv"
))
indep_csv  <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/ocean_independent_t_test_results_df.csv",
  "/mnt/data/ocean_independent_t_test_results_df.csv"
))

out_dir <- dirname(data_csv)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- knobs (match your other figures) ----------
BASE      <- 22
TITLE     <- 24
STAR_SIZE <- 7
CLD_SIZE  <- 8.2
BAR_W     <- 0.65
DODGE_W   <- 0.60
dash_at   <- c(-2, 2)

fill_values <- c("NRI"="#332288","NTI"="#DDCC77")    # bar fills
cld_cols    <- c(NRI = "#332288", NTI = "#A68000")   # CLD letters by metric

# Ordering: oceans don't have a canonical order, so default to order by a metric
KEEP_DEFINED_ORDER <- FALSE
defined_levels <- character(0)   # put a custom vector here if you have one
ORDER_METRIC   <- "NRI"          # or "NTI"
ORDER_DIR      <- "desc"         # "asc" or "desc"

# ---------- helpers ----------
p_to_stars <- function(p) ifelse(is.na(p), "ns",
                                 ifelse(p <= 0.001, "***",
                                        ifelse(p <= 0.01,  "**",
                                               ifelse(p <= 0.05, "*", "ns"))))

to_numeric <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\u2212", "-", x, fixed = TRUE)   # Unicode minus → ASCII
  x <- gsub("[^0-9eE+\\-\\.]", "", x)
  suppressWarnings(as.numeric(x))
}

# Replace ASCII hyphens inside labels by en dash “–” (so CLD pair names use a single ASCII hyphen)
safe_name <- function(x) { x <- trimws(as.character(x)); gsub("-", "–", x, fixed = TRUE) }

cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
  if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
  pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
  if (!length(pcol)) stop("Independent-tests table has no p-value column p_value/p.value/p-value.")
  pcol <- pcol[1]
  
  sub <- indep_df %>% filter(.data$Metric == metric,
                             .data$Group_High %in% cats,
                             .data$Group_Low  %in% cats)
  if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))
  
  n_hi <- safe_name(sub$Group_High)
  n_lo <- safe_name(sub$Group_Low)
  pair_names <- paste(pmin(n_hi, n_lo), pmax(n_hi, n_lo), sep = "-")
  pv <- sub[[pcol]]; names(pv) <- pair_names
  
  cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
  
  lookup <- setNames(cats, safe_name(cats))     # safe -> original
  out <- setNames(rep("", length(cats)), cats)
  mapped <- intersect(names(cld), names(lookup))
  if (length(mapped)) out[lookup[mapped]] <- cld[mapped]
  out
}

# ---------- read + normalize ----------
df <- read.csv(data_csv, check.names = FALSE)
names(df) <- trimws(names(df))

# Category / Metric / Value / Error
cat_col <- intersect(c("Ocean","ocean","Sea","sea","Ocean_Sea","Ocean/Sea",
                       "Salinity","salinity","Category","Categories"), names(df))
if (!length(cat_col)) stop("Data: couldn't find a category column. Available: ", paste(names(df), collapse = " | "))
names(df)[names(df) == cat_col[1]] <- "Category"

mcol <- intersect(c("Metrics","Metric"), names(df)); if (!length(mcol)) stop("No 'Metrics' column.")
names(df)[names(df) == mcol[1]] <- "Metrics"

vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df)); if (!length(vcol)) stop("No value column (NRI_and_NTI/Value).")
names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"

ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]

# Filter categories with Number < 5 (auto-detect count column)
num_col_idx <- grep("(?i)^(n$|n_|^n\\.|^n\\s*$|number|num|count|samples|n_samples|n_samp|n\\.samples|n\\.samp)$|(?i)number", names(df), perl = TRUE)
if (length(num_col_idx) > 0) {
  num_col <- names(df)[num_col_idx[1]]
  df[[num_col]] <- to_numeric(df[[num_col]])
  df <- df %>% filter(.data[[num_col]] >= 5)
} else {
  message("No 'Number'/'Count' column found — skipping <5 filter.")
}

# Clean rows & coerce numerics; keep only NRI/NTI; drop NA/blank/"NA"
df <- df %>%
  mutate(
    Category    = trimws(as.character(Category)),
    Metrics     = trimws(as.character(Metrics)),
    NRI_and_NTI = to_numeric(NRI_and_NTI),
    .__err__    = to_numeric(.__err__)
  ) %>%
  filter(!is.na(Category), nzchar(Category), toupper(Category) != "NA",
         !is.na(Metrics), Metrics %in% c("NRI","NTI"))

# ---------- order categories ----------
if (KEEP_DEFINED_ORDER && length(defined_levels) && all(defined_levels %in% unique(df$Category))) {
  ord_levels <- defined_levels
} else {
  ord_tbl <- df %>%
    group_by(Category) %>%
    summarise(order_mean = mean(ifelse(Metrics == ORDER_METRIC, NRI_and_NTI, NA), na.rm = TRUE),
              .groups = "drop") %>%
    arrange(if (ORDER_DIR == "asc") order_mean else desc(order_mean))
  ord_levels <- ord_tbl$Category
}
df$Category <- factor(df$Category, levels = ord_levels)
df$Metrics  <- factor(df$Metrics,  levels = c("NRI","NTI"))
df <- droplevels(df)

# ---------- paired stars ----------
paired <- read.csv(paired_csv, check.names = FALSE)
names(paired) <- trimws(names(paired))
pcat <- intersect(c("Category","Ocean","ocean","Sea","sea","Salinity","salinity"), names(paired))
if (!length(pcat)) stop("Paired table: no category column.")
names(paired)[names(paired) == pcat[1]] <- "Category"
pcol <- intersect(c("p_value","p.value","p-value"), names(paired)); if (!length(pcol)) stop("Paired table: no p-value column.")

paired <- paired %>%
  mutate(
    Category = trimws(as.character(.data$Category)),
    stars    = p_to_stars(.data[[pcol[1]]])
  ) %>%
  filter(Category %in% ord_levels, nzchar(Category), toupper(Category) != "NA") %>%
  mutate(Category = factor(Category, levels = ord_levels)) %>%
  select(Category, stars) %>%
  droplevels()

# ---------- dynamic positions & limits ----------
y_min_data <- min(df$NRI_and_NTI - ifelse(df$NRI_and_NTI < 0, abs(df$.__err__), 0), na.rm = TRUE)
y_max_data <- max(df$NRI_and_NTI + ifelse(df$NRI_and_NTI > 0, abs(df$.__err__), 0), na.rm = TRUE)
rng_base   <- range(c(y_min_data, y_max_data, dash_at), na.rm = TRUE)

pad <- max(diff(rng_base) * 0.06, 0.15)  # gap bars↔labels
tol <- max(diff(rng_base) * 0.02, 0.10)  # avoid dashed lines

cat_pos <- df %>%
  group_by(Category) %>%
  summarise(
    top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
    bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_star = ifelse(top <= 0, bot - pad, top + pad))

for (d in dash_at) {
  cat_pos$y_star <- ifelse(abs(cat_pos$y_star - d) < tol,
                           ifelse(cat_pos$y_star >= d, cat_pos$y_star + pad*0.5, cat_pos$y_star - pad*0.5),
                           cat_pos$y_star)
}
paired_ann <- left_join(paired, cat_pos, by = "Category")

# ---------- base plot ----------
plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
  geom_col(position = position_dodge(width = DODGE_W), width = BAR_W) +
  geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
                position = position_dodge(width = DODGE_W), width = 0.22, linewidth = 0.7) +
  geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
  scale_fill_manual(values = fill_values, name = "Metrics") +
  coord_flip() +
  labs(title = "", x = "Ocean / Sea", y = "NRI and NTI") +
  theme_minimal(base_size = BASE) +
  theme(
    axis.title.x = element_text(size = TITLE, face = "bold"),
    axis.title.y = element_text(size = TITLE, face = "bold"),
    axis.text.x  = element_text(size = BASE),
    axis.text.y  = element_text(size = BASE),
    legend.title = element_text(size = TITLE, face = "bold"),
    legend.text  = element_text(size = BASE),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 14, 10, 10)
  )

# ---------- (A) stars-only (with halo) ----------
plot_stars <- plot_base +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE + 1.2, color = "white"
  ) +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE, color = "black"
  )

# ---------- (B) stars + CLD letters (colored, hyphen-safe, with halo) ----------
indep <- read.csv(indep_csv, check.names = FALSE)
nm <- trimws(names(indep)); nm[!nzchar(nm) | is.na(nm)] <- paste0("X", which(!nzchar(nm) | is.na(nm)))
names(indep) <- nm
nm <- names(indep)
nm <- sub("^p[._-]?value$", "p_value", nm, ignore.case = TRUE)
nm <- sub("^Group[ _]High$", "Group_High", nm)
nm <- sub("^Group[ _]Low$",  "Group_Low",  nm)
nm <- sub("^Metric[ ]*$",    "Metric",     nm)
names(indep) <- nm
drop_cols <- grep("^(Unnamed: ?\\d+|X\\d+)$", names(indep))
if (length(drop_cols)) indep <- indep[, -drop_cols, drop = FALSE]

req <- c("Metric","Group_High","Group_Low","p_value")
missing <- setdiff(req, names(indep))
if (length(missing)) stop("Independent-tests table missing: ", paste(missing, collapse = ", "))

indep <- indep %>%
  mutate(
    Metric     = trimws(as.character(Metric)),
    Group_High = trimws(as.character(Group_High)),
    Group_Low  = trimws(as.character(Group_Low))
  )

cats <- levels(df$Category)
letters_df <- bind_rows(
  { lt <- cld_letters_for_metric(indep, "NRI", cats, alpha = 0.05)
  tibble(Metrics="NRI", Category=names(lt),
         cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) },
  { lt <- cld_letters_for_metric(indep, "NTI", cats, alpha = 0.05)
  tibble(Metrics="NTI", Category=names(lt),
         cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) }
) %>%
  filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA") %>%
  filter(Category %in% ord_levels) %>%
  mutate(Category = factor(Category, levels = ord_levels)) %>%
  droplevels()

plot_cld <- plot_stars  # includes halo stars already

if (nrow(letters_df) > 0) {
  letter_pos <- df %>%
    group_by(Category, Metrics) %>%
    summarise(
      y_top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
      y_bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_cld = ifelse(y_top <= 0, y_bot - pad*1.25, y_top + pad*1.25))
  
  # keep letters away from dashed lines & stars
  for (d in dash_at) {
    letter_pos$y_cld <- ifelse(abs(letter_pos$y_cld - d) < tol,
                               ifelse(letter_pos$y_cld >= d, letter_pos$y_cld + pad*0.5, letter_pos$y_cld - pad*0.5),
                               letter_pos$y_cld)
  }
  star_map <- cat_pos %>% select(Category, y_star)
  letter_pos <- left_join(letter_pos, star_map, by = "Category") %>%
    mutate(y_cld = ifelse(!is.na(y_star) & abs(y_cld - y_star) < (pad*0.6),
                          ifelse(y_cld >= y_star, y_cld + pad*0.5, y_cld - pad*0.5),
                          y_cld)) %>%
    select(-y_star)
  
  letters_ann <- inner_join(letters_df, letter_pos, by = c("Category","Metrics")) %>%
    filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA")
  
  plot_cld <- plot_cld +
    geom_text(
      data = letters_ann,
      aes(x = Category, y = y_cld, label = cld, group = Metrics),
      position = position_dodge(width = DODGE_W),
      size = CLD_SIZE + 1.4, fontface = "bold", color = "white"
    ) +
    geom_text(
      data = letters_ann,
      aes(x = Category, y = y_cld, label = cld, group = Metrics, color = Metrics),
      position = position_dodge(width = DODGE_W),
      size = CLD_SIZE, fontface = "bold"
    ) +
    scale_color_manual(values = cld_cols, guide = "none")
} else {
  message("CLD letters are empty (no matching categories or p-values).")
}

# ---------- y-limits (bars + stars + letters + dashed lines) ----------
all_y <- c(df$NRI_and_NTI + df$.__err__, df$NRI_and_NTI - df$.__err__,
           paired_ann$y_star,
           if (exists("letters_ann")) letters_ann$y_cld else numeric(0),
           dash_at)
ymin <- floor((min(all_y, na.rm = TRUE) - pad*0.4) * 10) / 10
ymax <- ceiling((max(all_y, na.rm = TRUE) + pad*0.4) * 10) / 10

plot_stars <- plot_stars + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))
plot_cld   <- plot_cld   + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))

# ---------- save ----------
file_stem_a <- file.path(out_dir, sprintf("ocean_nri_nti_stars_orderBy%s_%s", ORDER_METRIC, ORDER_DIR))
file_stem_b <- file.path(out_dir, sprintf("ocean_nri_nti_cld_orderBy%s_%s",   ORDER_METRIC, ORDER_DIR))

ggsave(paste0(file_stem_a, ".svg"), plot_stars, width = 10, height = 10)
ggsave(paste0(file_stem_a, ".png"), plot_stars, width = 10, height = 10, dpi = 300)
ggsave(paste0(file_stem_b, ".svg"), plot_cld,  width = 10, height = 10)
ggsave(paste0(file_stem_b, ".png"), plot_cld,  width = 10, height = 10, dpi = 300)

message("Saved:\n  ", paste0(file_stem_a, c(".svg",".png"), collapse = "\n  "),
        "\n  ", paste0(file_stem_b, c(".svg",".png"), collapse = "\n  "))

# optional: keep final object
plot5 <- plot_cld



########################################
# # ==============================
# # Sea NRI/NTI — stars + CLD (adapted to your files)
# # ==============================
# 
# options(stringsAsFactors = FALSE)
# suppressPackageStartupMessages({
#   library(ggplot2)
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
#   library(tibble)
#   library(multcompView)
#   library(igraph)
#   library(tools)
# })
# 
# # ---------- paths ----------
# choose_first_existing <- function(paths) {
#   ok <- paths[file.exists(paths)]
#   if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
#   ok[[1]]
# }
# data_csv   <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/sea.csv",
#   "/mnt/data/sea.csv"
# ))
# paired_csv <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/sea_paired_ttest_nri_nti_results.csv",
#   "/mnt/data/sea_paired_ttest_nri_nti_results.csv"
# ))
# indep_csv  <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/sea_independent_t_test_results_df.csv",
#   "/mnt/data/sea_independent_t_test_results_df.csv"
# ))
# 
# out_dir <- dirname(data_csv)
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # ---------- helpers ----------
# p_to_stars <- function(p) {
#   ifelse(is.na(p), "ns",
#          ifelse(p <= 0.001, "***",
#                 ifelse(p <= 0.01,  "**",
#                        ifelse(p <= 0.05, "*", "ns"))))
# }
# 
# # Replace internal ASCII hyphens in labels by en dash “–” so there's exactly
# # ONE ASCII hyphen (the pair separator) in names fed to multcompView.
# safe_name <- function(x) {
#   x <- trimws(as.character(x))
#   gsub("-", "–", x, fixed = TRUE)
# }
# 
# cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
#   if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
#   pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
#   if (length(pcol) == 0) stop("Independent-test table has no p-value column named p_value/p.value/p-value.")
#   pcol <- pcol[1]
#   
#   sub <- indep_df %>% filter(.data$Metric == metric,
#                              .data$Group_High %in% cats,
#                              .data$Group_Low  %in% cats)
#   if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))
#   
#   n_hi <- safe_name(sub$Group_High)
#   n_lo <- safe_name(sub$Group_Low)
#   pair_names <- paste(pmin(n_hi, n_lo), pmax(n_hi, n_lo), sep = "-")
#   
#   pv <- sub[[pcol]]
#   names(pv) <- pair_names
#   
#   cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
#   
#   # Map back to original category names
#   lookup <- setNames(cats, safe_name(cats))
#   out <- setNames(rep("", length(cats)), cats)
#   mapped <- intersect(names(cld), names(lookup))
#   if (length(mapped)) out[lookup[mapped]] <- cld[mapped]
#   out
# }
# 
# # ---------- read + normalize (SEA main table) ----------
# df <- read.csv(data_csv, check.names = FALSE)
# names(df) <- trimws(names(df))
# 
# # Category
# cat_col <- intersect(c("salinity","Salinity","Category","Categories"), names(df))
# if (length(cat_col) == 0) stop("Data: couldn't find a category column. Available columns are: ", paste(names(df), collapse = ", "))
# names(df)[names(df) == cat_col[1]] <- "Category"
# 
# # Metric
# mcol <- intersect(c("Metrics","Metric"), names(df))
# if (length(mcol) == 0) stop("Data: no 'Metrics'/'Metric' column.")
# names(df)[names(df) == mcol[1]] <- "Metrics"
# 
# # Value
# vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df))
# if (length(vcol) == 0) stop("Data: no value column (expected NRI_and_NTI/Value).")
# names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"
# 
# # Error
# ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
# df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]
# 
# # Keep order of appearance and drop unused levels to avoid “extra” empty bars
# df$Category <- factor(df$Category, levels = unique(df$Category))
# df$Metrics  <- factor(df$Metrics, levels = c("NRI","NTI"))
# df <- droplevels(df)
# cat_levels  <- levels(df$Category)
# 
# # ---------- paired stars ----------
# paired <- read.csv(paired_csv, check.names = FALSE)
# names(paired) <- trimws(names(paired))
# pcat <- intersect(c("Salinity","salinity","Category"), names(paired))
# if (length(pcat) == 0) stop("Paired table: no category column. Available: ", paste(names(paired), collapse = ", "))
# names(paired)[names(paired) == pcat[1]] <- "Category"
# 
# pcol <- intersect(c("p_value","p.value","p-value"), names(paired))
# if (length(pcol) == 0) stop("Paired table: no p-value column named p_value/p.value/p-value.")
# 
# paired <- paired %>%
#   mutate(Category = factor(trimws(as.character(.data$Category)), levels = cat_levels),
#          stars    = p_to_stars(.data[[pcol[1]]])) %>%
#   select(Category, stars)
# 
# # Positions for stars (change 0.15 to move stars)
# cat_pos <- df %>%
#   group_by(Category) %>%
#   summarise(
#     top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#     bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(y_pos = ifelse(top <= 0, bot - 0.45, top + 0.45))
# 
# paired_ann <- left_join(paired, cat_pos, by = "Category")
# 
# # ---------- base plot ----------
# dodge_w <- 0.6
# fill_values <- c("NRI"="#332288","NTI"="#DDCC77")
# dash_at <- c(-2,2)
# 
# plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
#   geom_bar(stat = "identity", position = position_dodge(width = dodge_w), width = 0.6) +
#   geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
#                 position = position_dodge(width = dodge_w), width = 0.2) +
#   geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
#   scale_fill_manual(values = fill_values) +
#   coord_flip() +
#   theme_minimal() +
#   labs(title = "", x = "Sea", y = "NRI and NTI", fill = "Metrics")
# 
# # ---------- (A) stars-only ----------
# plot_stars <- plot_base +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = 3.8
#   )
# 
# # ---------- (B) stars + CLD letters (hyphen-safe) ----------
# # Read & sanitize independent tests
# indep <- read.csv(indep_csv, check.names = FALSE)
# 
# # Fix empty/NA names first, then normalize
# nm <- trimws(names(indep))
# empty_idx <- which(!nzchar(nm) | is.na(nm))
# if (length(empty_idx) > 0) nm[empty_idx] <- paste0("X", empty_idx)
# names(indep) <- nm
# nm <- names(indep)
# nm <- sub("^p[._-]?value$", "p_value", nm, ignore.case = TRUE)
# nm <- sub("^Group[ _]High$", "Group_High", nm)
# nm <- sub("^Group[ _]Low$",  "Group_Low",  nm)
# nm <- sub("^Metric[ ]*$",    "Metric",     nm)
# names(indep) <- nm
# 
# # Drop index-like columns safely
# drop_cols <- grep("^(Unnamed: ?\\d+|X\\d+)$", names(indep))
# if (length(drop_cols)) indep <- indep[, -drop_cols, drop = FALSE]
# 
# # Check required columns
# req <- c("Metric","Group_High","Group_Low","p_value")
# missing <- setdiff(req, names(indep))
# if (length(missing)) stop("Independent-tests table is missing column(s): ", paste(missing, collapse = ", "))
# 
# # Tidy values to match Category factor exactly
# indep <- indep %>%
#   mutate(
#     Metric     = trimws(as.character(Metric)),
#     Group_High = trimws(as.character(Group_High)),
#     Group_Low  = trimws(as.character(Group_Low))
#   )
# 
# # Build CLD labels (safe to hyphens)
# cats <- levels(df$Category)
# letters_list <- lapply(c("NRI","NTI"), function(m) {
#   lt  <- cld_letters_for_metric(indep, metric = m, cats = cats, alpha = 0.05)
#   tibble(
#     Metrics  = m,
#     Category = names(lt),
#     cld      = gsub("(?<=.)(?=.)", " ", unname(lt), perl = TRUE)  # "AB" -> "A B"
#   )
# })
# letters_df <- bind_rows(letters_list)
# 
# plot_cld <- plot_base +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = 3.8
#   )
# 
# if (nrow(letters_df) > 0) {
#   letter_pos <- df %>%
#     group_by(Category, Metrics) %>%
#     summarise(
#       y_top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#       y_bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     mutate(y_pos = ifelse(y_top <= 0, y_bot - 1.90, y_top + 1.90))  # <-- adjust offset here
#   
#   # inner_join to drop rows lacking y positions (avoids “Removed X rows...” warnings)
#   letters_ann <- inner_join(letters_df, letter_pos, by = c("Category","Metrics")) %>%
#     filter(!is.na(cld), nzchar(cld))
#   
#   plot_cld <- plot_cld +
#     geom_text(
#       data = letters_ann,
#       aes(x = Category, y = y_pos, label = cld, group = Metrics),
#       position = position_dodge(width = dodge_w),
#       size = 3.4
#     )
# } else {
#   message("CLD letters computed empty (no matching categories or missing p-values).")
# }
# 
# # ---------- save ----------
# file_stem_a <- file.path(out_dir, "sea_nri_nti_stars")
# file_stem_b <- file.path(out_dir, "sea_nri_nti_cld")
# 
# ggsave(paste0(file_stem_a, ".svg"), plot_stars, width = 10, height = 10)
# ggsave(paste0(file_stem_a, ".png"), plot_stars, width = 10, height = 10, dpi = 300)
# ggsave(paste0(file_stem_b, ".svg"), plot_cld,  width = 10, height = 10)
# ggsave(paste0(file_stem_b, ".png"), plot_cld,  width = 10, height = 10, dpi = 300)
# 
# message("Saved:\n  ", paste0(file_stem_a, c(".svg",".png"), collapse = "\n  "),
#         "\n  ", paste0(file_stem_b, c(".svg",".png"), collapse = "\n  "))
#########################################################

# ==============================
# Sea NRI/NTI — stars + CLD (clean NA, colored letters by metric)
# ==============================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(multcompView)
  library(tools)
})

# ---------- paths ----------
choose_first_existing <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
  ok[[1]]
}
data_csv   <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/sea.csv",
  "/mnt/data/sea.csv"
))
paired_csv <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/sea_paired_ttest_nri_nti_results.csv",
  "/mnt/data/sea_paired_ttest_nri_nti_results.csv"
))
indep_csv  <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/sea_independent_t_test_results_df.csv",
  "/mnt/data/sea_independent_t_test_results_df.csv"
))

out_dir <- dirname(data_csv)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- knobs ----------
BASE      <- 22
TITLE     <- 24
STAR_SIZE <- 7
CLD_SIZE  <- 8.2      # larger, clearer
BAR_W     <- 0.65
DODGE_W   <- 0.60
dash_at   <- c(-2, 2)
fill_values <- c("NRI"="#332288","NTI"="#DDCC77")

# CLD letter colors by metric
cld_cols <- c(NRI = "#332288", NTI = "#A68000")  # deep blue vs. darker gold

# Order seas by mean NRI or NTI (“NRI” or “NTI”), ascending/descending
ORDER_METRIC <- "NRI"    # or "NTI"
ORDER_DIR    <- "desc"   # "asc" or "desc"

# ---------- helpers ----------
p_to_stars <- function(p) ifelse(is.na(p), "ns",
                                 ifelse(p <= 0.001, "***",
                                        ifelse(p <= 0.01,  "**",
                                               ifelse(p <= 0.05, "*", "ns"))))

# safe hyphen handling for multcompView
safe_name <- function(x) { x <- trimws(as.character(x)); gsub("-", "–", x, fixed = TRUE) }

# CLD letters from independent t-tests (metric-specific)
cld_letters_for_metric <- function(indep_df, metric, cats, alpha = 0.05) {
  if (is.null(indep_df) || nrow(indep_df) == 0) return(setNames(rep("", length(cats)), cats))
  pcol <- intersect(c("p_value","p.value","p-value"), names(indep_df))
  if (length(pcol) == 0) stop("Independent-test table has no p-value column p_value/p.value/p-value.")
  pcol <- pcol[1]
  
  sub <- indep_df %>% filter(.data$Metric == metric,
                             .data$Group_High %in% cats,
                             .data$Group_Low  %in% cats)
  if (nrow(sub) == 0) return(setNames(rep("", length(cats)), cats))
  
  n_hi <- safe_name(sub$Group_High)
  n_lo <- safe_name(sub$Group_Low)
  pair_names <- paste(pmin(n_hi, n_lo), pmax(n_hi, n_lo), sep = "-")
  
  pv <- sub[[pcol]]; names(pv) <- pair_names
  cld <- multcompView::multcompLetters(pv, threshold = alpha)$Letters
  
  lookup <- setNames(cats, safe_name(cats))
  out <- setNames(rep("", length(cats)), cats)
  mapped <- intersect(names(cld), names(lookup))
  if (length(mapped)) out[lookup[mapped]] <- cld[mapped]
  out
}

# Coerce numerics safely (handles Unicode minus)
to_numeric <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\u2212", "-", x, fixed = TRUE)
  x <- gsub("[^0-9eE+\\-\\.]", "", x)
  suppressWarnings(as.numeric(x))
}

# ---------- read + normalize (SEA) ----------
df <- read.csv(data_csv, check.names = FALSE)
names(df) <- trimws(names(df))

# Standardize core columns
cat_col <- intersect(c("Sea","sea","salinity","Salinity","Category","Categories"), names(df))
if (length(cat_col) == 0) stop("Data: couldn't find a category column. Available: ", paste(names(df), collapse = " | "))
names(df)[names(df) == cat_col[1]] <- "Category"

mcol <- intersect(c("Metrics","Metric"), names(df))
if (length(mcol) == 0) stop("Data: no 'Metrics'/'Metric' column.")
names(df)[names(df) == mcol[1]] <- "Metrics"

vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df))
if (length(vcol) == 0) stop("Data: no value column (expected NRI_and_NTI/Value).")
names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"

ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]

# ---- filter categories with Number < 5 (robust detection of column) ----
num_col_idx <- grep("(?i)^(n$|n_|^n\\.|^n\\s*$|number|num|count|samples|n_samples|n_samp|n\\.samples|n\\.samp)$|(?i)number", names(df), perl = TRUE)
if (length(num_col_idx) > 0) {
  num_col <- names(df)[num_col_idx[1]]
  df[[num_col]] <- to_numeric(df[[num_col]])
  df <- df %>% filter(.data[[num_col]] >= 5)
} else {
  message("No 'Number'/'Count' column found — skipping <5 filter.")
}

# Clean rows (also drop NA/blank/'NA' categories defensively)
df <- df %>%
  mutate(
    Category    = trimws(as.character(Category)),
    Metrics     = trimws(as.character(Metrics)),
    NRI_and_NTI = to_numeric(NRI_and_NTI),
    .__err__    = to_numeric(.__err__)
  ) %>%
  filter(!is.na(Category), nzchar(Category), toupper(Category) != "NA",
         !is.na(Metrics), Metrics %in% c("NRI","NTI"))

# ---------- ORDER seas by mean NRI or NTI ----------
ord_tbl <- df %>%
  group_by(Category) %>%
  summarise(
    order_mean = mean(ifelse(Metrics == ORDER_METRIC, NRI_and_NTI, NA), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(if (ORDER_DIR == "asc") order_mean else desc(order_mean))

ord_levels <- ord_tbl$Category
df$Category <- factor(df$Category, levels = ord_levels)
df$Metrics  <- factor(df$Metrics,  levels = c("NRI","NTI"))
df <- droplevels(df)

# ---------- paired stars ----------
paired <- read.csv(paired_csv, check.names = FALSE)
names(paired) <- trimws(names(paired))
pcat <- intersect(c("Sea","sea","Salinity","salinity","Category"), names(paired))
if (length(pcat) == 0) stop("Paired table: no category column. Available: ", paste(names(paired), collapse = " | "))
names(paired)[names(paired) == pcat[1]] <- "Category"

pcol <- intersect(c("p_value","p.value","p-value"), names(paired))
if (length(pcol) == 0) stop("Paired table: no p-value column named p_value/p.value/p-value.")

paired <- paired %>%
  mutate(
    Category = trimws(as.character(.data$Category)),
    stars    = p_to_stars(.data[[pcol[1]]])
  ) %>%
  filter(Category %in% ord_levels, nzchar(Category), toupper(Category) != "NA") %>%
  mutate(Category = factor(Category, levels = ord_levels)) %>%
  select(Category, stars) %>%
  droplevels()

# ---------- dynamic positions & limits (include bars, stars, letters, dashed lines) ----------
y_min_data <- min(df$NRI_and_NTI - ifelse(df$NRI_and_NTI < 0, abs(df$.__err__), 0), na.rm = TRUE)
y_max_data <- max(df$NRI_and_NTI + ifelse(df$NRI_and_NTI > 0, abs(df$.__err__), 0), na.rm = TRUE)
rng_base   <- range(c(y_min_data, y_max_data, dash_at), na.rm = TRUE)

pad <- max(diff(rng_base) * 0.06, 0.15)  # gap bars↔labels
tol <- max(diff(rng_base) * 0.02, 0.10)  # avoid dashed lines

cat_pos <- df %>%
  group_by(Category) %>%
  summarise(
    top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
    bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_star = ifelse(top <= 0, bot - pad, top + pad))

# nudge stars away from dashed lines
for (d in dash_at) {
  cat_pos$y_star <- ifelse(abs(cat_pos$y_star - d) < tol,
                           ifelse(cat_pos$y_star >= d, cat_pos$y_star + pad*0.5, cat_pos$y_star - pad*0.5),
                           cat_pos$y_star)
}
paired_ann <- left_join(paired, cat_pos, by = "Category")

# ---------- base plot ----------
plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
  geom_col(position = position_dodge(width = DODGE_W), width = BAR_W) +
  geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
                position = position_dodge(width = DODGE_W), width = 0.22, linewidth = 0.7) +
  geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
  scale_fill_manual(values = fill_values, name = "Metrics") +
  coord_flip() +
  labs(title = "", x = "Sea", y = "NRI and NTI") +
  theme_minimal(base_size = BASE) +
  theme(
    axis.title.x = element_text(size = TITLE, face = "bold"),
    axis.title.y = element_text(size = TITLE, face = "bold"),
    axis.text.x  = element_text(size = BASE),
    axis.text.y  = element_text(size = BASE),
    legend.title = element_text(size = TITLE, face = "bold"),
    legend.text  = element_text(size = BASE),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 14, 10, 10)
  )

# ---------- (A) stars-only (with halo) ----------
plot_stars <- plot_base +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE + 1.2, color = "white"
  ) +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE, color = "black"
  )

# ---------- (B) stars + CLD letters (colored by metric, halo) ----------
# Read, sanitize independent tests
indep <- read.csv(indep_csv, check.names = FALSE)
nm <- trimws(names(indep)); nm[!nzchar(nm) | is.na(nm)] <- paste0("X", which(!nzchar(nm) | is.na(nm)))
names(indep) <- nm
nm <- names(indep)
nm <- sub("^p[._-]?value$", "p_value", nm, ignore.case = TRUE)
nm <- sub("^Group[ _]High$", "Group_High", nm)
nm <- sub("^Group[ _]Low$",  "Group_Low",  nm)
nm <- sub("^Metric[ ]*$",    "Metric",     nm)
names(indep) <- nm
drop_cols <- grep("^(Unnamed: ?\\d+|X\\d+)$", names(indep))
if (length(drop_cols)) indep <- indep[, -drop_cols, drop = FALSE]

req <- c("Metric","Group_High","Group_Low","p_value")
missing <- setdiff(req, names(indep))
if (length(missing)) stop("Independent-tests table missing: ", paste(missing, collapse = ", "))

indep <- indep %>%
  mutate(
    Metric     = trimws(as.character(Metric)),
    Group_High = trimws(as.character(Group_High)),
    Group_Low  = trimws(as.character(Group_Low))
  )

cats <- levels(df$Category)
letters_df <- bind_rows(
  { lt <- cld_letters_for_metric(indep, "NRI", cats, alpha = 0.05)
  tibble(Metrics="NRI", Category=names(lt),
         cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) },
  { lt <- cld_letters_for_metric(indep, "NTI", cats, alpha = 0.05)
  tibble(Metrics="NTI", Category=names(lt),
         cld = gsub("(?<=.)(?=.)", " ", toupper(unname(lt)), perl = TRUE)) }
) %>%
  # remove any NA/blank/"NA" labels, and keep only categories we actually plot
  filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA") %>%
  filter(Category %in% ord_levels) %>%
  mutate(Category = factor(Category, levels = ord_levels)) %>%
  droplevels()

plot_cld <- plot_stars  # start from stars (already has halo stars)

if (nrow(letters_df) > 0) {
  # per-metric label positions (avoid overlap with bars & stars)
  letter_pos <- df %>%
    group_by(Category, Metrics) %>%
    summarise(
      y_top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
      y_bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_cld = ifelse(y_top <= 0, y_bot - pad*1.25, y_top + pad*1.25))
  
  # keep away from dashed lines & stars
  for (d in dash_at) {
    letter_pos$y_cld <- ifelse(abs(letter_pos$y_cld - d) < tol,
                               ifelse(letter_pos$y_cld >= d, letter_pos$y_cld + pad*0.5, letter_pos$y_cld - pad*0.5),
                               letter_pos$y_cld)
  }
  star_map <- cat_pos %>% select(Category, y_star)
  letter_pos <- left_join(letter_pos, star_map, by = "Category") %>%
    mutate(y_cld = ifelse(!is.na(y_star) & abs(y_cld - y_star) < (pad*0.6),
                          ifelse(y_cld >= y_star, y_cld + pad*0.5, y_cld - pad*0.5),
                          y_cld)) %>%
    select(-y_star)
  
  letters_ann <- inner_join(letters_df, letter_pos, by = c("Category","Metrics")) %>%
    filter(!is.na(cld), nzchar(trimws(cld)), toupper(trimws(cld)) != "NA")
  
  # draw CLD letters with white halo, then colored by metric
  plot_cld <- plot_cld +
    geom_text(
      data = letters_ann,
      aes(x = Category, y = y_cld, label = cld, group = Metrics),
      position = position_dodge(width = DODGE_W),
      size = CLD_SIZE + 1.4, fontface = "bold", color = "white"
    ) +
    geom_text(
      data = letters_ann,
      aes(x = Category, y = y_cld, label = cld, group = Metrics, color = Metrics),
      position = position_dodge(width = DODGE_W),
      size = CLD_SIZE, fontface = "bold"
    ) +
    scale_color_manual(values = cld_cols, guide = "none")
} else {
  message("CLD letters are empty (no matching categories or p-values).")
}

# ---------- y-limits (bars + stars + letters + dashed lines) ----------
all_y <- c(df$NRI_and_NTI + df$.__err__, df$NRI_and_NTI - df$.__err__,
           paired_ann$y_star,
           if (exists("letters_ann")) letters_ann$y_cld else numeric(0),
           dash_at)
ymin <- floor((min(all_y, na.rm = TRUE) - pad*0.4) * 10) / 10
ymax <- ceiling((max(all_y, na.rm = TRUE) + pad*0.4) * 10) / 10

plot_stars <- plot_stars + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))
plot_cld   <- plot_cld   + scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02)))

# ---------- save ----------
file_stem_a <- file.path(out_dir, sprintf("sea_nri_nti_stars_orderBy%s_%s", ORDER_METRIC, ORDER_DIR))
file_stem_b <- file.path(out_dir, sprintf("sea_nri_nti_cld_orderBy%s_%s",   ORDER_METRIC, ORDER_DIR))

ggsave(paste0(file_stem_a, ".svg"), plot_stars, width = 20, height = 16)
ggsave(paste0(file_stem_a, ".png"), plot_stars, width = 20, height = 16, dpi = 300)
ggsave(paste0(file_stem_b, ".svg"), plot_cld,  width = 20, height = 16)
ggsave(paste0(file_stem_b, ".png"), plot_cld,  width = 20, height = 16, dpi = 300)

message("Saved:\n  ", paste0(file_stem_a, c(".svg",".png"), collapse = "\n  "),
        "\n  ", paste0(file_stem_b, c(".svg",".png"), collapse = "\n  "))




########################################################
# # ==============================
# # Phylum NRI/NTI — stars only (paired t-test)
# # ==============================
# 
# options(stringsAsFactors = FALSE)
# suppressPackageStartupMessages({
#   library(ggplot2)
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
#   library(tibble)
#   library(tools)
# })
# 
# # ---------- paths ----------
# choose_first_existing <- function(paths) {
#   ok <- paths[file.exists(paths)]
#   if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
#   ok[[1]]
# }
# data_csv   <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/phylum.csv",
#   "/mnt/data/phylum.csv"
# ))
# paired_csv <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/phylum_paired_ttest_nri_nti_results.csv",
#   "/mnt/data/phylum_paired_ttest_nri_nti_results.csv"
# ))
# 
# out_dir <- dirname(data_csv)
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # ---------- helpers ----------
# p_to_stars <- function(p) {
#   ifelse(is.na(p), "ns",
#          ifelse(p <= 0.001, "***",
#                 ifelse(p <= 0.01,  "**",
#                        ifelse(p <= 0.05, "*", "ns"))))
# }
# 
# # Coerce messy numerics safely (fix unicode minus, strip noise)
# to_numeric <- function(x) {
#   x <- trimws(as.character(x))
#   x <- gsub("\u2212", "-", x, fixed = TRUE)         # Unicode minus -> ASCII hyphen
#   x <- gsub("[^0-9eE+\\-\\.]", "", x)               # keep digits/sign/decimal/exponent
#   suppressWarnings(as.numeric(x))
# }
# 
# # ---------- read + normalize (PHYLUM main table) ----------
# df <- read.csv(data_csv, check.names = FALSE)
# names(df) <- trimws(names(df))
# 
# # Category (your file uses 'salinity' for the phylum names)
# cat_col <- intersect(c("salinity","Salinity","Category","Categories"), names(df))
# if (length(cat_col) == 0) stop("Data: couldn't find a category column. Available columns are: ", paste(names(df), collapse = ", "))
# names(df)[names(df) == cat_col[1]] <- "Category"
# 
# # Metric
# mcol <- intersect(c("Metrics","Metric"), names(df))
# if (length(mcol) == 0) stop("Data: no 'Metrics'/'Metric' column.")
# names(df)[names(df) == mcol[1]] <- "Metrics"
# 
# # Value
# vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df))
# if (length(vcol) == 0) stop("Data: no value column (expected NRI_and_NTI/Value).")
# names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"
# 
# # Error (Std_dev, SD, ...)
# ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
# df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]
# 
# # Clean rows, coerce numerics
# df <- df %>%
#   mutate(
#     Category    = trimws(as.character(Category)),
#     Metrics     = trimws(as.character(Metrics)),
#     NRI_and_NTI = to_numeric(NRI_and_NTI),
#     .__err__    = to_numeric(.__err__)
#   ) %>%
#   filter(!is.na(Category), !is.na(Metrics), Metrics %in% c("NRI","NTI"))
# 
# # Order phylum levels: keep order of appearance (change to sort(unique(.)) to alphabetize)
# df$Category <- factor(df$Category, levels = unique(df$Category))
# df$Metrics  <- factor(df$Metrics,  levels = c("NRI","NTI"))
# df <- droplevels(df)
# cat_levels <- levels(df$Category)
# 
# # ---------- paired stars ----------
# paired <- read.csv(paired_csv, check.names = FALSE)
# names(paired) <- trimws(names(paired))
# 
# pcat <- intersect(c("Salinity","salinity","Category"), names(paired))
# if (length(pcat) == 0) stop("Paired table: no category column. Available: ", paste(names(paired), collapse = ", "))
# names(paired)[names(paired) == pcat[1]] <- "Category"
# 
# pcol <- intersect(c("p_value","p.value","p-value"), names(paired))
# if (length(pcol) == 0) stop("Paired table: no p-value column named p_value/p.value/p-value.")
# 
# paired <- paired %>%
#   mutate(
#     Category = factor(trimws(as.character(.data$Category)), levels = cat_levels),
#     stars    = p_to_stars(.data[[pcol[1]]])
#   ) %>%
#   select(Category, stars)
# 
# # Position for stars (change 0.15 to move stars)
# cat_pos <- df %>%
#   group_by(Category) %>%
#   summarise(
#     top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#     bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(y_pos = ifelse(top <= 0, bot - 0.55, top + 0.55))  # <-- adjust offset here
# 
# paired_ann <- left_join(paired, cat_pos, by = "Category")
# 
# # ---------- base plot (explicitly define; no stars yet) ----------
# dodge_w <- 0.6
# fill_values <- c("NRI"="#332288","NTI"="#DDCC77")
# dash_at <- c(-2, 2)
# 
# plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
#   geom_bar(stat = "identity", position = position_dodge(width = dodge_w), width = 0.6) +
#   geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
#                 position = position_dodge(width = dodge_w), width = 0.2) +
#   geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
#   scale_fill_manual(values = fill_values) +
#   coord_flip() +
#   theme_minimal() +
#   labs(title = "", x = "Phylum", y = "NRI and NTI", fill = "Metrics")
# 
# # ---------- final single figure (stars only) stored as `plot4` ----------
# plot_4 <- plot_base +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = 3.8
#   )
# 
# # ---------- save ONLY plot4 ----------
# file_stem <- file.path(out_dir, "phylum_nri_nti_plot4")
# ggsave(paste0(file_stem, ".svg"), plot4, width = 10, height = 10)
# ggsave(paste0(file_stem, ".png"), plot4, width = 10, height = 10, dpi = 300)
# message("Saved:\n  ", paste0(file_stem, c(".svg",".png"), collapse = "\n  "))
# 
# # ---------- OPTIONAL: if you ever compute CLD later ----------
# # if (exists("letters_ann") && is.data.frame(letters_ann) && nrow(letters_ann) > 0) {
# #   plot4 <- plot4 +
# #     geom_text(
# #       data = letters_ann,
# #       aes(x = Category, y = y_pos, label = cld, group = Metrics),
# #       position = position_dodge(width = dodge_w),
# #       size = 3.4
# #     )
# ggsave(paste0(file_stem, "_with_cld.svg"), plot_4, width = 10, height = 10)
# ggsave(paste0(file_stem, "_with_cld.png"), plot_4, width = 10, height = 10, dpi = 300)
# # }
# 
# # ---------- save ----------
# file_stem <- file.path(out_dir, "phylum_nri_nti_stars")
# ggsave(paste0(file_stem, ".svg"), plot_stars, width = 10, height = 10)
# ggsave(paste0(file_stem, ".png"), plot_stars, width = 10, height = 10, dpi = 300)
# message("Saved:\n  ", paste0(file_stem, c(".svg",".png"), collapse = "\n  "))
#######################################

# ==============================
# Phylum NRI/NTI — paired t-test ONLY (stars, no CLD)
# ==============================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(tools)
})

# ---------- paths ----------
choose_first_existing <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (!length(ok)) stop("None of these exist: ", paste(paths, collapse = " | "))
  ok[[1]]
}
data_csv   <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/phylum.csv",
  "/mnt/data/phylum.csv"
))
paired_csv <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/phylum_paired_ttest_nri_nti_results.csv",
  "/mnt/data/phylum_paired_ttest_nri_nti_results.csv"
))
out_dir <- dirname(data_csv)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- knobs (match your style) ----------
BASE      <- 22
TITLE     <- 24
STAR_SIZE <- 7
BAR_W     <- 0.65
DODGE_W   <- 0.60
dash_at   <- c(-2, 2)

fill_values <- c("NRI"="#332288","NTI"="#DDCC77")

# Ordering: choose which metric to sort phyla by
KEEP_DEFINED_ORDER <- FALSE
defined_levels <- character(0)   # put fixed order here if you want to lock it
ORDER_METRIC   <- "NRI"          # or "NTI"
ORDER_DIR      <- "desc"         # "asc" or "desc"

# ---------- helpers ----------
p_to_stars <- function(p) ifelse(is.na(p), "ns",
                                 ifelse(p <= 0.001, "***",
                                        ifelse(p <= 0.01,  "**",
                                               ifelse(p <= 0.05, "*", "ns"))))
to_numeric <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\u2212", "-", x, fixed = TRUE)   # Unicode minus → ASCII
  x <- gsub("[^0-9eE+\\-\\.]", "", x)
  suppressWarnings(as.numeric(x))
}

# ---------- read + normalize (PHYLUM main table) ----------
df <- read.csv(data_csv, check.names = FALSE)
names(df) <- trimws(names(df))

# Category (phylum name)
cat_col <- intersect(c("Phylum","phylum","Category","Categories","salinity","Salinity"), names(df))
if (!length(cat_col)) stop("Couldn't find a phylum/category column.")
names(df)[names(df) == cat_col[1]] <- "Category"

# Metric / Value / Error
mcol <- intersect(c("Metrics","Metric"), names(df)); if (!length(mcol)) stop("No 'Metrics' column.")
names(df)[names(df) == mcol[1]] <- "Metrics"
vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df)); if (!length(vcol)) stop("No value column.")
names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"
ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]

# Optional: filter out small sample sizes (<5) if a count column exists
num_col_idx <- grep("(?i)(^n$|^n[_ .]?samples?$|^count$|counts?|number)", names(df), perl = TRUE)
if (length(num_col_idx) > 0) {
  num_col <- names(df)[num_col_idx[1]]
  df[[num_col]] <- to_numeric(df[[num_col]])
  df <- df %>% filter(.data[[num_col]] >= 5)
} else {
  message("No 'Number/Count' column found — skipping <5 filter.")
}

# Clean & coerce
df <- df %>%
  mutate(
    Category    = trimws(as.character(Category)),
    Metrics     = trimws(as.character(Metrics)),
    NRI_and_NTI = to_numeric(NRI_and_NTI),
    .__err__    = to_numeric(.__err__)
  ) %>%
  filter(!is.na(Category), nzchar(Category), toupper(Category) != "NA",
         !is.na(Metrics),  Metrics %in% c("NRI","NTI"))

# ---------- order phyla ----------
if (KEEP_DEFINED_ORDER && length(defined_levels) && all(defined_levels %in% unique(df$Category))) {
  ord_levels <- defined_levels
} else {
  ord_tbl <- df %>%
    group_by(Category) %>%
    summarise(order_mean = mean(ifelse(Metrics == ORDER_METRIC, NRI_and_NTI, NA), na.rm = TRUE),
              .groups = "drop") %>%
    arrange(if (ORDER_DIR == "asc") order_mean else desc(order_mean))
  ord_levels <- ord_tbl$Category
}
df$Category <- factor(df$Category, levels = ord_levels)
df$Metrics  <- factor(df$Metrics,  levels = c("NRI","NTI"))
df <- droplevels(df)

# ---------- paired stars ----------
paired <- read.csv(paired_csv, check.names = FALSE)
names(paired) <- trimws(names(paired))
pcat <- intersect(c("Phylum","phylum","Category","Salinity","salinity"), names(paired))
if (!length(pcat)) stop("Paired table: no category column.")
names(paired)[names(paired) == pcat[1]] <- "Category"
pcol <- intersect(c("p_value","p.value","p-value"), names(paired)); if (!length(pcol)) stop("Paired table: no p-value column.")

paired <- paired %>%
  mutate(
    Category = trimws(as.character(.data$Category)),
    stars    = p_to_stars(.data[[pcol[1]]])
  ) %>%
  filter(Category %in% ord_levels, nzchar(Category), toupper(Category) != "NA") %>%
  mutate(Category = factor(Category, levels = ord_levels)) %>%
  select(Category, stars) %>%
  droplevels()

# ---------- star positions & y-limits ----------
y_min_data <- min(df$NRI_and_NTI - ifelse(df$NRI_and_NTI < 0, abs(df$.__err__), 0), na.rm = TRUE)
y_max_data <- max(df$NRI_and_NTI + ifelse(df$NRI_and_NTI > 0, abs(df$.__err__), 0), na.rm = TRUE)
rng_base   <- range(c(y_min_data, y_max_data, dash_at), na.rm = TRUE)
pad <- max(diff(rng_base) * 0.06, 0.15)
tol <- max(diff(rng_base) * 0.02, 0.10)

cat_pos <- df %>%
  group_by(Category) %>%
  summarise(
    top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
    bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_star = ifelse(top <= 0, bot - pad, top + pad))

for (d in dash_at) {
  cat_pos$y_star <- ifelse(abs(cat_pos$y_star - d) < tol,
                           ifelse(cat_pos$y_star >= d, cat_pos$y_star + pad*0.5, cat_pos$y_star - pad*0.5),
                           cat_pos$y_star)
}
paired_ann <- left_join(paired, cat_pos, by = "Category")

# ---------- base plot ----------
plot_base <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
  geom_col(position = position_dodge(width = DODGE_W), width = BAR_W) +
  geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
                position = position_dodge(width = DODGE_W), width = 0.22, linewidth = 0.7) +
  geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
  scale_fill_manual(values = fill_values, name = "Metrics") +
  coord_flip() +
  labs(title = "", x = "Phylum", y = "NRI and NTI") +
  theme_minimal(base_size = BASE) +
  theme(
    axis.title.x = element_text(size = TITLE, face = "bold"),
    axis.title.y = element_text(size = TITLE, face = "bold"),
    axis.text.x  = element_text(size = BASE),
    axis.text.y  = element_text(size = BASE),
    legend.title = element_text(size = TITLE, face = "bold"),
    legend.text  = element_text(size = BASE),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 14, 10, 10)
  )

# ---------- final (stars only, haloed) ----------
plot_stars <- plot_base +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE + 1.2, color = "white"
  ) +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_star, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE, color = "black"
  )

# dynamic y-limits to include bars, stars, and dashed lines
all_y <- c(df$NRI_and_NTI + df$.__err__, df$NRI_and_NTI - df$.__err__,
           paired_ann$y_star, dash_at)
ymin <- floor((min(all_y, na.rm = TRUE) - 0.1) * 10) / 10
ymax <- ceiling((max(all_y, na.rm = TRUE) + 0.1) * 10) / 10
plot_stars <- plot_stars + scale_y_continuous(limits = c(ymin, ymax),
                                              expand = expansion(mult = c(0.02, 0.02)))

# ---------- save ----------
file_stem <- file.path(out_dir, sprintf("phylum_nri_nti_stars_orderBy%s_%s", ORDER_METRIC, ORDER_DIR))
ggsave(paste0(file_stem, ".svg"), plot_stars, width = 10, height = 10)
ggsave(paste0(file_stem, ".png"), plot_stars, width = 10, height = 10, dpi = 300)

# for your patchwork assembly
plot4 <- plot_stars



#######################################


# # --- EXTENSION: single figure with BOTH stars and CLD, stored as `plot4` ---
# 
# # (optional) small cleanup to avoid "Removed X rows" warnings
# paired_ann2  <- tryCatch(dplyr::filter(paired_ann, !is.na(y_pos), !is.na(stars)), error = function(e) paired_ann)
# letters_ann2 <- if (exists("letters_ann") && is.data.frame(letters_ann)) {
#   dplyr::filter(letters_ann, !is.na(y_pos), !is.na(cld), cld != "")
# } else NULL
# 
# plot_4 <- plot_base +
#   # paired t-test stars
#   geom_text(
#     data = paired_ann2,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold",
#     size = 3.8
#   ) +
#   # CLD letters (only if available)
#   {
#     if (!is.null(letters_ann2) && nrow(letters_ann2) > 0) {
#       geom_text(
#         data = letters_ann2,
#         aes(x = Category, y = y_pos, label = cld, group = Metrics),
#         position = position_dodge(width = dodge_w),
#         size = 3.4
#       )
#     } else {
#       NULL
#     }
#   }
# 
# # Use it later with: print(plot4)
# 
# # (optional) save to file when you want
# # ggsave(file.path(out_dir, "plot4.svg"), plot4, width = 10, height = 8)
# # ggsave(file.path(out_dir, "plot4.png"), plot4, width = 10, height = 8, dpi = 300)
# # (optional) save the R object:
# # saveRDS(plot4, file.path(out_dir, "plot4.rds"))
# 
# #dev.off()
# 
# 
# ##########################################

# # ==============================
# # Class NRI/NTI — stars only (paired t-test)
# # ==============================
# 
# options(stringsAsFactors = FALSE)
# suppressPackageStartupMessages({
#   library(ggplot2)
#   library(dplyr)
#   library(tidyr)
#   library(stringr)
#   library(tibble)
#   library(tools)
# })
# 
# # ---------- paths ----------
# choose_first_existing <- function(paths) {
#   ok <- paths[file.exists(paths)]
#   if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
#   ok[[1]]
# }
# data_csv   <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/class.csv",
#   "/mnt/data/class.csv"
# ))
# paired_csv <- choose_first_existing(c(
#   "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/class_paired_ttest_nri_nti_results.csv",
#   "/mnt/data/class_paired_ttest_nri_nti_results.csv"
# ))
# 
# out_dir <- dirname(data_csv)
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # ---------- helpers ----------
# p_to_stars <- function(p) {
#   ifelse(is.na(p), "ns",
#          ifelse(p <= 0.001, "***",
#                 ifelse(p <= 0.01,  "**",
#                        ifelse(p <= 0.05, "*", "ns"))))
# }
# 
# # Coerce messy numerics safely (handles Unicode minus etc.)
# to_numeric <- function(x) {
#   x <- trimws(as.character(x))
#   x <- gsub("\u2212", "-", x, fixed = TRUE)     # Unicode minus → ASCII hyphen
#   x <- gsub("[^0-9eE+\\-\\.]", "", x)           # keep digits/sign/decimal/exponent
#   suppressWarnings(as.numeric(x))
# }
# 
# # ---------- read + normalize (CLASS main table) ----------
# df <- read.csv(data_csv, check.names = FALSE)
# names(df) <- trimws(names(df))
# 
# # Category (your file uses 'salinity' with class names)
# cat_col <- intersect(c("salinity","Salinity","Category","Categories"), names(df))
# if (length(cat_col) == 0) stop("Data: couldn't find a category column. Available: ", paste(names(df), collapse = ", "))
# names(df)[names(df) == cat_col[1]] <- "Category"
# 
# # Metric
# mcol <- intersect(c("Metrics","Metric"), names(df))
# if (length(mcol) == 0) stop("Data: no 'Metrics'/'Metric' column.")
# names(df)[names(df) == mcol[1]] <- "Metrics"
# 
# # Value
# vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df))
# if (length(vcol) == 0) stop("Data: no value column (expected NRI_and_NTI/Value).")
# names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"
# 
# # Error
# ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
# df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]
# 
# # Clean + coerce numerics; keep only NRI/NTI rows
# df <- df %>%
#   mutate(
#     Category    = trimws(as.character(Category)),
#     Metrics     = trimws(as.character(Metrics)),
#     NRI_and_NTI = to_numeric(NRI_and_NTI),
#     .__err__    = to_numeric(.__err__)
#   ) %>%
#   filter(!is.na(Category), !is.na(Metrics), Metrics %in% c("NRI","NTI"))
# 
# # Order class levels: keep order of appearance (change to sort(unique(.)) to alphabetize)
# df$Category <- factor(df$Category, levels = unique(df$Category))
# df$Metrics  <- factor(df$Metrics,  levels = c("NRI","NTI"))
# df <- droplevels(df)
# cat_levels  <- levels(df$Category)
# 
# # ---------- paired stars only ----------
# paired <- read.csv(paired_csv, check.names = FALSE)
# names(paired) <- trimws(names(paired))
# 
# pcat <- intersect(c("Salinity","salinity","Category"), names(paired))
# if (length(pcat) == 0) stop("Paired table: no category column. Available: ", paste(names(paired), collapse = ", "))
# names(paired)[names(paired) == pcat[1]] <- "Category"
# 
# pcol <- intersect(c("p_value","p.value","p-value"), names(paired))
# if (length(pcol) == 0) stop("Paired table: no p-value column named p_value/p.value/p-value.")
# 
# paired <- paired %>%
#   mutate(
#     Category = factor(trimws(as.character(.data$Category)), levels = cat_levels),
#     stars    = p_to_stars(.data[[pcol[1]]])
#   ) %>%
#   select(Category, stars)
# 
# # Star positions (tweak 0.15 to move stars up/down)
# cat_pos <- df %>%
#   group_by(Category) %>%
#   summarise(
#     top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
#     bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(y_pos = ifelse(top <= 0, bot - 0.75, top + 0.75))
# 
# paired_ann <- left_join(paired, cat_pos, by = "Category")
# 
# # ---------- plot ----------
# dodge_w <- 0.6
# fill_values <- c("NRI"="#332288","NTI"="#DDCC77")
# dash_at <- c(-2, 2)
# 
# plot_stars <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
#   geom_bar(stat = "identity", position = position_dodge(width = dodge_w), width = 0.6) +
#   geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
#                 position = position_dodge(width = dodge_w), width = 0.2) +
#   geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
#   scale_fill_manual(values = fill_values) +
#   coord_flip() +
#   theme_minimal() +
#   labs(title = "", x = "Class", y = "NRI and NTI", fill = "Metrics") +
#   geom_text(
#     data = paired_ann,
#     aes(x = Category, y = y_pos, label = stars),
#     inherit.aes = FALSE,
#     fontface = "bold", size = 5.5
#   )
# 
# # ---------- save ----------
# file_stem <- file.path(out_dir, "class_nri_nti_stars")
# ggsave(paste0(file_stem, ".svg"), plot_stars, width = 10, height = 10)
# ggsave(paste0(file_stem, ".png"), plot_stars, width = 10, height = 10, dpi = 300)
# message("Saved:\n  ", paste0(file_stem, c(".svg",".png"), collapse = "\n  "))
##########################################
# ==============================
# Class NRI/NTI — stars only (paired t-test) — improved & robust
# ==============================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(tools)
})

# ---------- paths ----------
choose_first_existing <- function(paths) {
  ok <- paths[file.exists(paths)]
  if (length(ok) == 0) stop("None of these exist: ", paste(paths, collapse = " | "))
  ok[[1]]
}
data_csv   <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/class.csv",
  "/mnt/data/class.csv"
))
paired_csv <- choose_first_existing(c(
  "D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/new_nri_nti_analysis/class_paired_ttest_nri_nti_results.csv",
  "/mnt/data/class_paired_ttest_nri_nti_results.csv"
))

out_dir <- dirname(data_csv)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- helpers ----------
p_to_stars <- function(p) {
  ifelse(is.na(p), "ns",
         ifelse(p <= 0.001, "***",
                ifelse(p <= 0.01,  "**",
                       ifelse(p <= 0.05, "*", "ns"))))
}

# Coerce messy numerics safely (handles Unicode minus etc.)
to_numeric <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\u2212", "-", x, fixed = TRUE)     # Unicode minus → ASCII hyphen
  x <- gsub("[^0-9eE+\\-\\.]", "", x)           # keep digits/sign/decimal/exponent
  suppressWarnings(as.numeric(x))
}

# ---------- read + normalize (CLASS main table) ----------
df <- read.csv(data_csv, check.names = FALSE)
names(df) <- trimws(names(df))

# Category (e.g., salinity class)
cat_col <- intersect(c("salinity","Salinity","Category","Categories"), names(df))
if (length(cat_col) == 0) stop("Data: couldn't find a category column. Available: ", paste(names(df), collapse = ", "))
names(df)[names(df) == cat_col[1]] <- "Category"

# Metric
mcol <- intersect(c("Metrics","Metric"), names(df))
if (length(mcol) == 0) stop("Data: no 'Metrics'/'Metric' column.")
names(df)[names(df) == mcol[1]] <- "Metrics"

# Value
vcol <- intersect(c("NRI_and_NTI","Value","value"), names(df))
if (length(vcol) == 0) stop("Data: no value column (expected NRI_and_NTI/Value).")
names(df)[names(df) == vcol[1]] <- "NRI_and_NTI"

# Error (optional)
ecol <- intersect(c("Std_dev","SD","sd","SE","se","Std error","Std_error"), names(df))
df$.__err__ <- if (length(ecol) == 0) 0 else df[[ecol[1]]]

# Clean + keep only NRI/NTI rows
df <- df %>%
  mutate(
    Category    = trimws(as.character(Category)),
    Metrics     = trimws(as.character(Metrics)),
    NRI_and_NTI = to_numeric(NRI_and_NTI),
    .__err__    = to_numeric(.__err__)
  ) %>%
  filter(!is.na(Category), !is.na(Metrics), Metrics %in% c("NRI","NTI"))

# ---------- paired stars only ----------
paired <- read.csv(paired_csv, check.names = FALSE)
names(paired) <- trimws(names(paired))

pcat <- intersect(c("Salinity","salinity","Category"), names(paired))
if (length(pcat) == 0) stop("Paired table: no category column. Available: ", paste(names(paired), collapse = ", "))
names(paired)[names(paired) == pcat[1]] <- "Category"

pcol <- intersect(c("p_value","p.value","p-value"), names(paired))
if (length(pcol) == 0) stop("Paired table: no p-value column named p_value/p.value/p-value.")

paired <- paired %>%
  mutate(
    Category = trimws(as.character(.data$Category)),
    stars    = p_to_stars(.data[[pcol[1]]])
  ) %>%
  select(Category, stars)

# ---------- ordering options & style knobs ----------
ORDER_BY  <- "abs_mean"  # choose: "mean", "abs_mean", "NRI", "NTI"
ORDER_DIR <- "desc"      # "asc" or "desc"

BASE      <- 22          # general text size
TITLE     <- 24          # axis/legend title size
STAR_SIZE <- 7           # "ns/*/**/***" size
BAR_W     <- 0.65
DODGE_W   <- 0.60
dash_at   <- c(-2, 2)    # reference lines to KEEP

fill_values <- c("NRI"="#332288","NTI"="#DDCC77")

# ---------- compute ordering ----------
ord_tbl <- df %>%
  group_by(Category) %>%
  summarize(
    mean     = mean(NRI_and_NTI, na.rm = TRUE),
    abs_mean = mean(abs(NRI_and_NTI), na.rm = TRUE),
    NRI      = mean(ifelse(Metrics == "NRI", NRI_and_NTI, NA), na.rm = TRUE),
    NTI      = mean(ifelse(Metrics == "NTI", NRI_and_NTI, NA), na.rm = TRUE),
    .groups = "drop"
  )

ord_vec <- ord_tbl[[ORDER_BY]]
ord_tbl <- ord_tbl %>% arrange(if (ORDER_DIR == "asc") ord_vec else desc(ord_vec))
ord_levels <- ord_tbl$Category

df$Category     <- factor(df$Category, levels = ord_levels)
paired$Category <- factor(paired$Category, levels = ord_levels)
df$Metrics      <- factor(df$Metrics, levels = c("NRI","NTI"))
df              <- droplevels(df)
paired          <- droplevels(paired)

# ---------- star positions: dynamic, non-overlapping with bars ----------
# data range including asymmetric error bars
y_min_data <- min(df$NRI_and_NTI - ifelse(df$NRI_and_NTI < 0, abs(df$.__err__), 0), na.rm = TRUE)
y_max_data <- max(df$NRI_and_NTI + ifelse(df$NRI_and_NTI > 0, abs(df$.__err__), 0), na.rm = TRUE)

# include the dashed reference lines in the baseline range
rng_base <- range(c(y_min_data, y_max_data, dash_at), na.rm = TRUE)

pad <- max(diff(rng_base) * 0.06, 0.15)   # gap between bars and star text
tol <- max(diff(rng_base) * 0.02, 0.10)   # avoid dashed lines

cat_pos <- df %>%
  group_by(Category) %>%
  summarize(
    top = max(NRI_and_NTI + ifelse(NRI_and_NTI >= 0, .__err__, 0), na.rm = TRUE),
    bot = min(NRI_and_NTI - ifelse(NRI_and_NTI <  0, .__err__, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(y_pos = ifelse(top <= 0, bot - pad, top + pad))

# nudge stars away from dashed lines if too close
for (d in dash_at) {
  cat_pos$y_pos <- ifelse(abs(cat_pos$y_pos - d) < tol,
                          ifelse(cat_pos$y_pos >= d, cat_pos$y_pos + pad*0.5, cat_pos$y_pos - pad*0.5),
                          cat_pos$y_pos)
}

paired_ann <- left_join(paired, cat_pos, by = "Category")

# y-limits include bars, stars, AND dashed lines (robust fallback)
ymin <- min(rng_base[1], min(cat_pos$y_pos, na.rm = TRUE)) - pad*0.3
ymax <- max(rng_base[2], max(cat_pos$y_pos, na.rm = TRUE)) + pad*0.3
if (!is.finite(ymin) | !is.finite(ymax) | (ymax <= ymin)) {
  ymin <- rng_base[1] - abs(rng_base[1])*0.1 - 0.5
  ymax <- rng_base[2] + abs(rng_base[2])*0.1 + 0.5
}
# round a bit so ticks look nice
ymin <- floor(ymin * 10) / 10
ymax <- ceiling(ymax * 10) / 10

# ---------- plot ----------
plot_stars <- ggplot(df, aes(x = Category, y = NRI_and_NTI, fill = Metrics)) +
  geom_col(position = position_dodge(width = DODGE_W), width = BAR_W) +
  geom_errorbar(aes(ymin = NRI_and_NTI - .__err__, ymax = NRI_and_NTI + .__err__),
                position = position_dodge(width = DODGE_W), width = 0.2, linewidth = 0.6) +
  geom_hline(yintercept = dash_at, linetype = "dashed", color = "#CC79A7", linewidth = 1) +
  scale_fill_manual(values = fill_values, name = "Metrics") +
  coord_flip() +
  labs(title = "", x = "Class", y = "NRI and NTI") +
  geom_text(
    data = paired_ann,
    aes(x = Category, y = y_pos, label = stars),
    inherit.aes = FALSE,
    fontface = "bold", size = STAR_SIZE
  ) +
  scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02, 0.02))) +
  theme_minimal(base_size = BASE) +
  theme(
    axis.title.x = element_text(size = TITLE, face = "bold"),
    axis.title.y = element_text(size = TITLE, face = "bold"),
    axis.text.x  = element_text(size = BASE),
    axis.text.y  = element_text(size = BASE),
    legend.title = element_text(size = TITLE, face = "bold"),
    legend.text  = element_text(size = BASE),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linetype = "dotted", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 14, 10, 10)
  )

# ---------- save ----------
file_stem <- file.path(out_dir, sprintf("class_nri_nti_stars_%s_%s", ORDER_BY, ORDER_DIR))
ggsave(paste0(file_stem, ".svg"), plot_stars, width = 10, height = 10)
ggsave(paste0(file_stem, ".png"), plot_stars, width = 10, height = 10, dpi = 300)
message("Saved:\n  ", paste0(file_stem, c(".svg",".png"), collapse = "\n  "))


##########################################
# library(patchwork)
# 
# # (Optional) add tags
# plot1 <- plot1 + labs(tag = "(a)")
# plot2 <- plot2 + labs(tag = "(b)")
# plot3 <- plot3 + labs(tag = "(c)")
# plot4 <- plot4 + labs(tag = "(d)")
# plot5 <- plot5 + labs(tag = "(e)")
# 
# # helper: hide legend for a plot
# p_nolegend <- function(p) p + theme(legend.position = "none")
# 
# # build columns: left has plot1–3; right has plot4–5 (keep legend on plot4)
# left_col  <- p_nolegend(plot1) / p_nolegend(plot2) / p_nolegend(plot3)
# right_col <- plot4 / p_nolegend(plot5) / plot_spacer()
# 
# final_plot <- (left_col | right_col) +
#   plot_layout(guides = "collect", widths = c(1, 1)) &
#   theme(
#     legend.position = "bottom",
#     legend.text  = element_text(size = 25),
#     legend.title = element_text(size = 25)
#   )
# 
# # save if you want
# ggsave(file.path(out_dir, "13_merged.png"), final_plot, width = 14, height = 20, dpi = 600)
# #ggsave(file.path(out_dir, "13_merged.svg"), final_plot, width = 16, height = 12)

library(patchwork)

# tag panels (kept from your code)
plot1 <- plot1 + labs(tag = "(a)")
plot2 <- plot2 + labs(tag = "(b)")
plot3 <- plot3 + labs(tag = "(c)")
plot4 <- plot4 + labs(tag = "(d)")
plot5 <- plot5 + labs(tag = "(e)")

# helpers
strip_legend <- function(p) p + theme(legend.position = "none")

# reclaim a bit of space on the right column so areas look proportional
tight_right <- function(p, left = -10, right = 6, top = 6, bottom = 6) {
  p + theme(
    plot.margin = margin(top, right, bottom, left),
    axis.title.y = element_text(margin = margin(r = 8))
  )
}

# layout: (d) spans two rows on the right; (e) sits below it
design <- "
AB
CB
DE
"

final_layout <-
  wrap_plots(
    A = strip_legend(plot1),            # row1, col1
    B = tight_right(plot4),             # row1–2, col2 (bigger)
    C = strip_legend(plot2),            # row2, col1
    D = strip_legend(plot3),            # row3, col1
    E = tight_right(strip_legend(plot5)), # row3, col2
    design = design
  ) +
  plot_layout(
    guides  = "collect",                # collect one legend from any plot that has it
    widths  = c(0.9, 1.1),              # give the right column a touch more room
    heights = c(1, 1, 1)                # equal row heights; panel (d) == rows 1+2
  ) &
  theme(
    legend.position = "bottom",
    legend.text  = element_text(size = 25),
    legend.title = element_text(size = 25),
    plot.tag     = element_text(size = 18, face = "bold")
  )

# Save
ggsave(file.path(out_dir, "13_merged_plots.svg"), final_layout, width = 14, height = 18)
ggsave(file.path(out_dir, "13_merged_plots.png"), final_layout, width = 14, height = 18, dpi = 600)

########################





library(ggplot2)
library(patchwork)

# Simulating 3 other smaller plots (replace with real plots)
plot1 <- plot1 + labs(tag = "(a)")
plot2 <- plot2 + labs(tag = "(b)")
plot3 <- plot3 + labs(tag = "(c)")
plot4 <- plot4 + labs(tag = "(d)")
plot5 <- plot5 + labs(tag = "(e)")

# Extract the common legend
get_legend <- function(a_plot) {
  legend <- ggplotGrob(a_plot + theme(legend.position = "bottom"))$grobs
  legend <- legend[[which(sapply(legend, function(x) x$name) == "guide-box")]]
  return(legend)
}

legend <- get_legend(ggplot(data, aes(x = Ocean, y = NRI_and_NTI, fill = Metrics)) +
                       geom_bar(stat = "identity", position = "dodge") +
                       scale_fill_manual(values = c("NRI" = "#332288", "NTI" = "#DDCC77")) +
                       theme(legend.position = "bottom", legend.text = element_text(size = 25), 
                             legend.title = element_text(size = 25)))

# Arrange the plots (3 smaller on left, 2 larger on right)
final_layout <- ((plot1 / plot2 / plot3) | (plot4 / plot5)) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Save as SVG
svg("/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI/13_merged_plots.svg", width = 18, height = 22)
print(final_layout)
#grid::grid.draw(legend)  # Draw the common legend
dev.off()





