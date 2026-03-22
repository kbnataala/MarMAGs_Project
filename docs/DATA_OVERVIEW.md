# Data Overview

## Purpose

This document summarizes the datasets used in the MarMAGs project and clarifies which components are included in this repository and which are not.

---

## Study Scope

This study analyzed marine metagenome-assembled genomes (MAGs) to investigate the biogeographic distribution, phylogenetic structure, and ecological assembly of carbon-fixing prokaryotes in marine systems.

### Summary of scale

- **Marine metagenomes analyzed:** 6,561
- **Total MAGs analyzed:** 134,644
- **Newly recovered MAGs:** 73,944
- **MAGs from previous public datasets:** 60,700
- **Dominant carbon-fixation pathway:** Calvin–Benson–Bassham (CBB)
- **CBB-encoding MAGs:** 4,253

---

## Data Included in This Repository

This repository includes selected final outputs and curated materials used to support the manuscript.

### Included content

- Final curated result tables used in the manuscript
- Supplementary tables (where file size permits)
- Final manuscript figures
- Analysis and plotting scripts provided as documentation templates
- Upstream command templates describing MAG characterization and related processing
- Documentation of software versions and manual table curation

---

## Data Not Included in This Repository

The following are **not redistributed** in this repository:

- Raw sequencing reads
- Full metagenomic assemblies
- Intermediate binning outputs
- Large phylogenetic intermediate files
- GTDB reference databases
- METABOLIC databases
- Full HPC working directories

These files are omitted due to their size and the computational infrastructure required to generate and process them.

---

## Raw Data Provenance

Raw sequencing reads were obtained from public marine metagenomic repositories, including sources accessible through the Sequence Read Archive (SRA) and public marine genome resources.

The upstream workflow included:

1. Selection of marine metagenomes
2. MAG recovery using MuDoGeR
3. Quality filtering with CheckM
4. Taxonomic classification using GTDB-Tk
5. Functional annotation of carbon fixation pathways
6. Phylogenetic and ecological analyses

These upstream steps are documented in `pipeline_reconstruction/`.

---

## Processed Data Used for Analysis

Downstream analyses reported in the manuscript were conducted using curated processed tables derived from upstream outputs.

These processed tables contain information such as:

- Genome identifiers
- Taxonomic classification
- Sample coordinates
- Environmental metadata
- Carbon fixation pathway presence/absence
- NRI and NTI values
- Prevalence summaries
- Pairwise statistical comparisons

---

## Repository Data Philosophy

This repository is structured as a **documentation and transparency companion** rather than a complete reproducible pipeline.

The goal is to provide:

- transparency of analytical logic
- access to final curated outputs
- command templates for upstream methods
- documentation of manual curation decisions

rather than a full end-to-end re-executable workflow.

---

## Notes on Manual Curation

Some final tables were manually harmonized prior to inclusion in the repository. These edits included operations such as:

- renaming column headers
- standardizing categorical labels
- reformatting tables for consistency
- harmonizing matrix structures

Details are recorded in:

`docs/MANUAL_CURATION_LOG.md`
