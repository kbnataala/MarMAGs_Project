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

### Reproducible from this repository

- Statistical analyses (manyGLM, prevalence tests, correlation analyses, t-tests)
- NRI and NTI calculations using curated phylogenetic subsets
- Environmental gradient binning and prevalence calculations
- All manuscript figures
- Final result tables

### Not reproducible from this repository alone

- Raw sequencing read download
- MAG assembly and binning workflows
- Large intermediate HPC outputs
- Full phylogenetic reconstruction from raw genomes

Upstream workflows are documented in the `pipeline_reconstruction/` directory but require external data and computational resources.

---

## Repository Structure

```
MarMAGs_Project/
├── README.md
├── data/
│   ├── raw/                  # Pointer directory (raw data not hosted)
│   ├── processed/            # Curated analytical datasets used in analysis
│   └── metadata/             # Harmonised environmental/sample metadata
├── scripts/                  # Sequential analysis scripts
├── figures/                  # Generated figures (main + supplementary)
├── results/                  # Statistical outputs and summary tables
├── environment/              # R and Python dependency documentation
├── docs/                     # Manual curation log + data dictionary
└── pipeline_reconstruction/  # Documented upstream workflow from Methods
```

---

## Data Description

### Raw Data  
Raw sequencing reads originate from public marine metagenomes (e.g., SRA). These are not redistributed here due to size and infrastructure constraints.

See `data/raw/README.md` for acquisition guidance.

### Processed Data  
All analyses begin from curated processed datasets located in:

```
data/processed/
```

Some processed tables were manually harmonised (e.g., column renaming, label standardisation, matrix formatting). All modifications are documented in:

```
docs/MANUAL_CURATION_LOG.md
```

This ensures transparency and traceability.

---

## Analysis Workflow

Scripts are organised sequentially:

1. `01_data_preprocessing.R`
2. `02_mag_characterization.R`
3. `03_carbon_fixation.R`
4. `04_phylogenetic_analysis.R`
5. `05_biogeography.R`
6. `06_statistical_analysis.R`
7. `07_figures_main.R`

Each script writes outputs to the `results/` and `figures/` directories.

---

## Quick Start

Clone the repository:

```bash
git clone https://github.com/kbnataala/MarMAGs_Project.git
cd MarMAGs_Project
```

Install R dependencies:

```r
source("environment/install_packages.R")
```

Run the reproducible subset:

```r
source("run_all.R")
```

Outputs will be written to:

- `results/`
- `figures/`

---

## Software Environment

Analyses were conducted using:

- R ≥ 4.3
- Key packages: tidyverse, vegan, ape, phangorn, ggtree, mvabund, picante
- Python (for geospatial processing and web application components)

Exact package versions are listed in:

```
environment/R_packages.txt
```

---

## Web Application

The MarMAGs interactive web application allows users to:

- Filter MAGs by taxonomy, quality metrics, and environmental attributes
- Query carbon fixation pathways and gene presence
- Visualise geographic distributions
- Download filtered datasets

Access: https://webapp.ufz.de/marmags/

---

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
