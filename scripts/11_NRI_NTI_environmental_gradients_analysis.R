##########################
#set working directory
setwd("D:/Thesis/work_package_3/Data_analysis/01_figures/06_NRI_NTI")


# Load necessary libraries
library(ggplot2)
library(dplyr)

#loading the data
data <- read.csv(file = "scatter_plot_all.csv", header = TRUE, stringsAsFactors = FALSE)

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
## NRI and NTI plots
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

####################################################################

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





