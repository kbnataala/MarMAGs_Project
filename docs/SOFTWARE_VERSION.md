# Software Versions

## Purpose

This document records the main software tools, databases, and computational environments used in the MarMAGs study.

---

## Genome Recovery and Characterization

### MuDoGeR
- **Tool:** MuDoGeR
- **Version:** v1.0.1

### CheckM
- **Tool:** CheckM
- **Version:** v1.1.6

### GTDB-Tk
- **Tool:** GTDB-Tk
- **Version:** v2.4.0
- **Reference database:** GTDB Release 220

### BBTools
- **Tool:** BBTools / statswrapper.sh
- **Version:** [add if known]

### dRep
- **Tool:** dRep
- **Version:** v3.2.2

---

## Functional Annotation

### METABOLIC
- **Tool:** METABOLIC
- **Version:** v4.0

### Prodigal
- **Tool:** Prodigal
- **Version:** v2.6.3

### HMMER
- **Tool:** HMMER
- **Version:** v3.3.2

### Annotation databases
- **KOfam database:** July 2019 release
- **Pfam:** version 32.0
- **TIGRfam:** version 15.0

---

## Phylogenetic Analysis

### FastTree
- **Tool:** FastTree
- **Version:** v2.1.11

### Treemmer
- **Tool:** Treemmer
- **Version:** v0.3

### PICANTE
- **Tool/package:** picante
- **Language:** R

---

## Statistical Analysis

### R
- **Language:** R
- **Version:** 4.3 or later

### Key R packages
- tidyverse
- vegan
- ape
- phangorn
- ggtree
- mvabund
- picante
- ggplot2
- patchwork
- cowplot
- pheatmap
- gridExtra
- writexl

> Exact package versions can be added from `sessionInfo()` if available.

---

## Geospatial Processing and Web Application

### Python
- **Language:** Python
- **Version:** 3.10

### Python libraries
- geopandas
- shapely
- pandas
- streamlit
- plotly express
- folium
- matplotlib
- st-aggrid

---

## Computational Environment

### Operating environment
- Linux-based high-performance computing (HPC) environment for upstream genome recovery and characterization
- Local R/Python environments for downstream figure generation and table processing

---

## Notes

Some exact versions may depend on the original HPC environment and installed package state at the time of analysis. Where exact patch versions are unavailable, the documented versions reflect the versions stated in the manuscript methods and supplementary methods.
