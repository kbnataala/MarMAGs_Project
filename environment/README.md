# Environment

This directory documents the downstream software environment used for analysis and figure generation in the MarMAGs project.

## Contents

- `R_packages.txt` — list of R packages and versions used in the analysis and visualization workflows

## Notes

The package list was compiled from the analysis environment used during the study and is provided for transparency. It is not intended to function as a lockfile or exact reproducible environment specification.

Upstream analyses performed on HPC infrastructure additionally relied on external tools such as:

- MuDoGeR
- CheckM
- GTDB-Tk
- BBTools
- METABOLIC
- FastTree
- Treemmer

These tools and versions are summarized in `docs/SOFTWARE_VERSIONS.md`.
