---
title: "MAG Characterization Pipeline (Reconstruction)"

---

# Overview

This directory reconstructs the MAG characterization steps used in the study:

> *Habitat filtering constrained by evolutionary history drives the assemblage of carbon-fixing species in marine systems*

The pipeline mirrors the MuDoGeR Module 2 workflow and includes:

- CheckM (genome quality metrics)
- GTDB-Tk (taxonomy assignment; GTDB release 220)
- BBTools (genome statistics)
- Quality Score filtering (QS)

All MAG FASTA files are assumed to be located in a single directory.

---

# Required Software

- **CheckM**
- **GTDB-Tk v2.4.0**
- **BBTools**
- Python ≥ 3.8
- R (optional, for rendering)

GTDB reference database (release 220) must be installed and the environment variable set:

```bash
export GTDBTK_DATA_PATH=/path/to/gtdb/release220
```

---

# Input

MAG genome files:

```
path/to/mags/folder/*.fa
```

All MAG FASTA files must share the same extension (default: `.fa`).

---

# Step 1 — Genome Quality Assessment (CheckM)

```bash
bash 01_run_checkm.sh input/mags/folder checkm/output/result/folder 32 fa
```

Output:

```
results/folder/mag_characterization/checkm/output_checkm.tsv
```

Key columns used:

- Completeness
- Contamination
- Strain heterogeneity

---

# Step 2 — Taxonomic Classification (GTDB-Tk)

```bash
bash 02_run_gtdbtk.sh data/mags results/mag_characterization/gtdbtk 32 fa
```

Output:

```
results/mag_characterization/gtdbtk/gtdbtk_result.tsv
```

Taxonomy is assigned according to GTDB release 220.

---

# Step 3 — Genome Metrics (BBTools)

```bash
bash 03_run_bbtools.sh data/mags results/mag_characterization/bbtools fa
```

Output:

```
results/mag_characterization/bbtools/bbtools.tsv
```

Metrics include:

- Genome size (bp)
- GC content
- N50
- Number of scaffolds

---

# Step 4 — Merge Results and Compute Quality Score

```bash
python3 04_merge_metrics.py \
  --checkm results/mag_characterization/checkm/output_checkm.tsv \
  --bbtools results/mag_characterization/bbtools/bbtools.tsv \
  --gtdb results/mag_characterization/gtdbtk/gtdbtk_result.tsv \
  --outdir results/mag_characterization/merged \
  --qs_cutoff 50
```

---

# Quality Score Definition

Quality Score (QS) is defined as:

\[
QS = \text{Completeness} - 5 \times \text{Contamination}
\]

MAGs are defined as genomes with:

\[
QS \geq 50
\]

This follows the same definition used in the manuscript and MuDoGeR module outputs.

---

# Final Outputs

```
results/mag_characterization/merged/
├── mag_metadata.tsv
├── mag_taxonomy.tsv
└── mags_qs50.tsv
```

- **mag_metadata.tsv** — combined genome metrics, taxonomy, and QS  
- **mag_taxonomy.tsv** — taxonomy split into ranks  
- **mags_qs50.tsv** — MAGs passing QS ≥ 50  

---

# Notes

- This pipeline reconstructs upstream genome characterization.
- It does not include assembly, binning, or dereplication.
- Large intermediate files are not tracked in Git.
- Shapefiles and reference databases must be installed locally.

---

# Reproducibility Statement

These scripts replicate the MAG characterization strategy described in the Methods section of the manuscript and correspond to the MuDoGeR Module 2 workflow.

