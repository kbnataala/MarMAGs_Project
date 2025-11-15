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

## 05_05_Multivariate_biogeography_of_all_and_carbon_fixing_phyla.R

R script for testing and visualizing the biogeographic structure of marine prokaryotes and carbon-fixing phyla across environmental gradients and regions. The workflow builds phylum-by-sample presence matrices, fits mvabund manyGLM models for depth, salinity, temperature, oceans and IHO seas, and exports multivariate and univariate ANOVA summaries. It also produces publication grade clustered heatmaps of pathway and CBB presence by phylum and region/gradient to summarize community turnover in environmental and geographic space.

## 06_Correlation_analysis.R

R workflow for gradient-based co-occurrence analysis of all MAGs accross the different phyla and CBB encoding MAGs in each phyla across temperature, depth, and salinity, building presence–absence matrices, testing monotonicity, computing smoothed (σ = 8) curve correlations with circular-shift permutation + FDR, and generating heatmaps / multi-panel figures for significant phylum–phylum associations and top-50 networks.



























