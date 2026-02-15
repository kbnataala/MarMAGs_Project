# MarMAGs Project: Biogeography and Phylogenetic Analysis of Marine Carbon-Fixing Metagenome-Assembled Genomes

---

## Overview

This repository contains all analysis code, curated processed datasets, and figure-generation workflows supporting the manuscript:

**Title**: Habitat filtering constrained by evolutionary history drives the assemblage of carbon-fixing species in marine systems  
**Authors**: Muhammad Kabiru Nata’ala et al.  
**Corresponding Author**: Ulisses Rocha  
**Journal**: [To be updated upon publication]  
**DOI**: [To be updated upon publication]

This study investigates the biogeographic distribution and phylogenetic structure of marine metagenome-assembled genomes (MAGs) encoding carbon fixation pathways. The project integrates large-scale genome recovery, functional annotation, phylogenomics, and ecological modeling to evaluate how evolutionary history and environmental gradients shape marine carbon-fixing microbial communities.

This repository reproduces all downstream analyses and figures derived from curated processed datasets included here. Upstream high-performance computing (HPC) steps are documented for transparency but are not fully executable from this repository due to data size and infrastructure constraints.

---

## Abstract

### Background  
Microbial carbon fixation is central to marine food webs and global biogeochemical cycles. Yet, the drivers of its distribution—particularly among uncultured lineages—remain unclear. The “everything is everywhere, but the environment selects” (EiE-BES) hypothesis suggests broad microbial dispersal constrained by environmental selection but often overlooks phylogenetic constraints on key functional traits such as carbon fixation.

### Results  
We analysed 134,644 metagenome-assembled genomes (MAGs) derived from 6,561 marine metagenomes, including 73,944 newly recovered genomes and 60,700 from public datasets. Species-level operational taxonomic units were delineated and placed within a phylogenomic framework using GTDB-Tk (release 220).

Carbon fixation capacity was widespread but dominated by taxa encoding the Calvin–Benson–Bassham (CBB) cycle (4,253 MAGs). Geography, depth, salinity, and temperature significantly influenced the distributions of both total and carbon-fixing MAGs (p < 0.001). CBB-capable MAGs were phylogenetically clustered within specific taxonomic groups and associated with distinct environmental regimes.

Phylogenetic structure metrics (Net Relatedness Index and Nearest Taxon Index) consistently indicated strong clustering across environmental gradients. These findings suggest that carbon fixation traits are evolutionarily conserved and that marine microbial community assembly is shaped by habitat filtering constrained by evolutionary history.

An interactive web application enables exploration and download of MAG characteristics and carbon fixation annotations.

### Conclusion  
Habitat filtering constrained by evolutionary history shapes the assembly of marine carbon-fixing prokaryotes. This extends the EiE-BES framework by incorporating evolutionary limits to environmental selection and provides a foundation for predicting microbial responses to environmental change.

---

## Key Findings

- Carbon fixation pathways exhibit strong phylogenetic constraint across marine prokaryotes.
- Environmental gradients significantly structure both total and carbon-fixing MAG distributions.
- CBB-capable MAGs are phylogenetically clustered and associated with distinct environmental regimes.
- Phylogenetic clustering patterns remain robust under tree-size reduction.
- Habitat filtering dominates over neutral assembly processes across marine environments.

---

## Reproducibility Scope

###Included in this repository
 - Final curated processed tables used for manuscript analyses
 - Statistical analysis templates
 - Figure-generation scripts
 - Documentation of analytical steps
 - Command templates for upstream genome characterization

###Not included
 - Raw sequencing reads
 - MAG assembly and binning workflows
 - Large intermediate HPC outputs
 - Phylogenetic tree reconstruction from raw genomes
 - GTDB reference databases

Upstream workflows are described in pipeline_reconstruction/ but require external data and computational resources.

---

## Repository Structure

```
MarMAGs_Project/
├── README.md
├── results/
│   ├── tables/                 # Final curated manuscript tables
│   ├── supplementary_tables/   # Supplementary tables (if included)
│   └── figures/                # Final manuscript figures
├── scripts/
│   ├── figures/                # Figure-generation templates
│   └── analysis_templates/     # Statistical analysis examples
├── pipeline_reconstruction/
│   └── mag_characterization/   # CheckM / GTDB-Tk / BBTools command templates
├── docs/
│   ├── DATA_OVERVIEW.md
│   ├── SOFTWARE_VERSIONS.md
│   ├── MANUAL_CURATION_LOG.md
│   └── WORKFLOW_SUMMARY.md
└── environment/
    └── R_packages.txt

```

## Processed Data
All manuscript analyses begin from curated processed datasets used to generate final statistical outputs and figures. These are included in:

results/tables/

Some tables were manually harmonised (e.g., column renaming, label standardisation, matrix formatting). All modifications are documented in:

docs/MANUAL_CURATION_LOG.md

---

## Analysis Documentation

This repository starts from curated processed tables located in `data/processed/`.  
Each script operates independently on its own input table(s) and performs local input checks before running.

Scripts are organised logically according to the manuscript structure:

1. `01_mag_characterization.R`  
2. `02_carbon_fixation.R`  
3. `03_phylogenetic_analysis.R`  
4. `04_biogeography.R`  
5. `05_statistical_analysis.R`  
6. `06_figures_main.R`  

These scripts document the analytical approach used in the manuscript.
They may require adaptation depending on local paths and software environments.

## Upstream Genome Characterization

The directory:

pipeline_reconstruction/mag_characterization/

contains command templates reconstructing:
 - CheckM quality assessment
 - GTDB-Tk taxonomy assignment (release 220)
 - BBTools genome metrics
 - Quality Score calculation (QS ≥ 50 threshold)

These scripts are provided for methodological transparency and require external data and HPC resources.

## Software Environment

Analyses were conducted using:
 - R ≥ 4.3
 - Key packages: tidyverse, vegan, ape, phangorn, ggtree, mvabund, picante
 - Python (for geospatial processing and data integration)
 - CheckM
 - GTDB-Tk v2.4.0
 - BBTools

Exact software versions are documented in:


docs/SOFTWARE_VERSIONS.md


## Web Application

The MarMAGs interactive web application allows users to:
 - Filter MAGs by taxonomy, quality metrics, and environmental attributes
 - Query carbon fixation pathways and gene presence
 - Visualise geographic distributions
 - Download filtered datasets

Access:

https://webapp.ufz.de/marmags/

## Citation

If you use this repository, please cite:

Muhammad, K. N. et al.  
Habitat filtering constrained by evolutionary history drives the assemblage of carbon-fixing species in marine systems.  
[Journal, Year, DOI]

---

## License

See `LICENSE` for licensing information.

---

## Contact

Muhammad Kabiru Nata’ala  
Email: kmuhammad@atb-potsdam.de  

Corresponding author:  
Ulisses Rocha – ulisses.rocha@ufz.de
