# Workflow Summary

## Purpose

This document provides a concise overview of the analytical workflow used in the MarMAGs project, from raw metagenomic data to ecological and phylogenetic interpretation.

---

## Overview of Workflow

The study combined public marine metagenomic data, genome-resolved analysis, functional annotation, and eco-phylogenetic modeling to investigate the assembly of carbon-fixing marine microorganisms.

---

## Step 1 — Marine Metagenome Selection

Marine metagenomes were selected from public repositories based on study scope and data availability. The selection focused on marine systems and large-scale metagenomic datasets suitable for MAG recovery.

---

## Step 2 — MAG Recovery

Metagenome-assembled genomes (MAGs) were recovered using MuDoGeR. This included:

- quality control of reads
- assembly
- binning
- MAG extraction

---

## Step 3 — Quality Assessment

Recovered genomes were assessed using CheckM.

Key metrics included:

- completeness
- contamination
- strain heterogeneity

A quality score (QS) was computed as:

\[
QS = completeness - 5 \times contamination
\]

MAGs were defined as genomes with:

\[
QS \geq 50
\]

---

## Step 4 — Taxonomic Assignment

MAG taxonomy was assigned using GTDB-Tk with GTDB release 220.

This provided standardized taxonomic labels across ranks, including:

- domain
- phylum
- class
- order
- family
- genus
- species

---

## Step 5 — Genome Metrics

Additional genome statistics were summarized using BBTools, including:

- genome size
- number of scaffolds
- N50
- GC content

---

## Step 6 — Functional Annotation

Carbon fixation pathways were identified using METABOLIC and associated databases.

The main pathways considered included:

- Calvin–Benson–Bassham cycle (CBB)
- reverse tricarboxylic acid cycle (rTCA)
- 3-hydroxypropionate bicycle (3-HPB)
- dicarboxylate/hydroxybutyrate cycle (DCHB)
- hydroxypropionate/hydroxybutyrate cycle (HPHB)
- bacterial Wood–Ljungdahl pathway (bWLJ)
- archaeal Wood–Ljungdahl pathway (aWLJ)

---

## Step 7 — Dereplication and Species-Level Clustering

MAGs were grouped into species-level operational taxonomic units using dRep with a 95% ANI threshold.

---

## Step 8 — Phylogenetic Reconstruction

Phylogenetic reconstruction was carried out using GTDB-Tk alignments and FastTree. Tree-size reduction sensitivity analyses were further performed using Treemmer.

---

## Step 9 — Environmental and Biogeographic Analysis

Environmental metadata and geospatial annotations were used to evaluate associations with:

- depth
- salinity
- temperature
- marine regions
- seas and oceans

These analyses included multivariate modeling, prevalence analysis, and spatial summaries.

---

## Step 10 — Phylogenetic Structure Analysis

Community phylogenetic structure was assessed using:

- Net Relatedness Index (NRI)
- Nearest Taxon Index (NTI)

These analyses were used to test for clustering and infer the role of habitat filtering versus neutral assembly.

---

## Step 11 — Figure and Table Generation

Final figures and manuscript tables were generated from curated processed outputs. In this repository, scripts are retained as documentation templates to illustrate how these figures and summaries were produced.

---

## Repository Role in This Workflow

This repository does **not** reproduce the full workflow end-to-end.

Instead, it serves as a:

- documentation companion
- transparency archive
- collection of final curated outputs
- set of command and script templates

for the thesis and manuscript.
