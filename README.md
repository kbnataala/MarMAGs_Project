# MarMAGs_Project

This repository contains the scripts used for analysis in the manuscript "An analyses of almost 135 thousand metagenome-assembled genomes indicate that habitat filtering constrained by evolutionary history drives the assemblage of carbon-fixing species in marine systems".

Authors: Muhammad Kabiru Nata’ala, Anderson P. Avila Santos, Jonas Coelho Kasmanas, Eriik Borchert, Jan Zarzycki, Newton C. M. Gomes, Rodrigo Costa, Tina Keller‑Costa, Mirjam perner, Marcel Niklausz, Sabine Kleinsteuber, Ivan A Berg, Tobias J. Erb, Peter F. Stadler, Ulisses Rocha


## 01_Retrieving_ocean_and_sea_names.py

Scripts for assigning ocean/sea names and maritime regions to sample coordinates using spatial joins with global marine shapefiles (GOaS and EEZ–IHO). Enables standardized geographic annotation for downstream analyses.

## 02_Distribution_of_metagenomes.R

R script for generating publication-ready figures summarizing marine MAG libraries, including pie charts of salinity/temperature/depth categories, bar plots of ocean/sea/IHO/marine region/territory distributions, and a global map of sample locations colored by genome source (MarMAGs, GEM catalog, OceanDNA).

## 03_Percentage_of_missingness_metadata.R

R script for quantifying and visualizing metadata missingness in the MarMAGs dataset, generating a field-wise bar plot (percentage of missing values per metadata field) and a flipped sample × field missingness heatmap with automatic sizing for large tables.

## 04_Metagenome-assembled-genomes_quality_metrics.R

R script for comparing genome quality metrics (CheckM completeness, contamination, strain heterogeneity, and quality score) across genome sources (MarMAGs, GEM catalog, OceanDNA), using Levene’s tests, Welch ANOVA, pairwise Welch t-tests, and generating violin/box plots with significance annotations plus a combined multi-panel figure.

## 05_Multivariate_biogeography_of_all_and_carbon_fixing_phyla.R

R script for testing and visualizing the biogeographic structure of marine prokaryotes and carbon-fixing phyla across environmental gradients and regions. The workflow builds phylum-by-sample presence matrices, fits mvabund manyGLM models for depth, salinity, temperature, oceans and IHO seas, and exports multivariate and univariate ANOVA summaries. It also produces publication grade clustered heatmaps of pathway and CBB presence by phylum and region/gradient to summarize community turnover in environmental and geographic space.

## 06_Correlation_analysis.R

R workflow for gradient-based co-occurrence analysis of all MAGs accross the different phyla and CBB encoding MAGs in each phyla across temperature, depth, and salinity, building presence–absence matrices, testing monotonicity, computing smoothed (σ = 8) curve correlations with circular-shift permutation + FDR, and generating heatmaps / multi-panel figures for significant phylum–phylum associations and top-50 networks.

## 07_Prevalence_of_CBB_encoding_MAGs.R

R workflow for quantifying and comparing carbon-fixation pathway prevalence across taxonomic groups and environmental categories (phyla, classes, depth, temperature, salinity, oceans/seas), using proportion tests for pairwise differences and generating publication-ready barplots with prevalence labels and significance brackets, plus combined multi-panel figures.

## 08_Bubble_plot_showing_the_presence_and_absence_of_key_carbon_fixation_genes.R

R script for generating a bubble plot of carbon-fixation genes across the top 30 phyla, merging taxonomy and pathway presence/absence tables, aggregating genome counts, scaling bubble size by % genomes per phylum–gene–pathway, annotating phylum sample sizes, and exporting high-resolution PNG/PDF/SVG figures.

## 09_Phylogenetic_structure_nri_nti.R

R script for computing phylogenetic structure metrics (NRI and NTI) for each community using an abundance matrix and a GTDB phylogenetic tree. The workflow calculates ses.mpd and ses.mntd with 999 null permutations, merges NRI/NTI outputs by phylum, and exports the combined results to CSV.

## 10_Tree_reduction_sensitivity_nri_nti.R

This script benchmarks how sensitive NRI and NTI are to using reduced phylogenetic trees. It reads precomputed NRI/NTI values across multiple regions for different percentages of the tree retained, fits regressions of metric vs. tree size, and standardizes slopes to quantify effect size. It then (i) produces a forest plot of standardized slopes with an equivalence band, (ii) computes and plots rank preservation (Spearman ρ) relative to the full tree, and (iii) generates faceted regression panels for NRI and NTI. Finally, it exports these three publication-grade figures plus minimal regression and rank-preservation summary tables.

## 11_NRI_NTI_environmental_gradients_analysis.R

R script for analyzing and visualizing NRI and NTI across environmental categories (e.g. salinity, temperature, depth). The workflow computes independent and paired t-tests, derives significance stars and compact letter displays, and generates publication-grade barplots and scatter plots with error bars and reference lines, together with minimal summary tables.















