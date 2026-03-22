# Manual Curation Log

## Purpose

This document records manual changes made to curated tables included in the repository.

Because this repository is intended as a documentation and transparency companion, it is important to explicitly note where outputs were manually harmonized rather than regenerated directly from scripts.

---

## Types of Manual Changes Recorded Here

Examples of manual curation include:

- renaming column headers
- standardizing spelling or category labels
- reordering columns
- merging or splitting table fields
- removing empty or redundant columns
- adapting table titles for consistency with the manuscript
- formatting supplementary tables for presentation

---

## Curation Table

| File | Source / Upstream origin | Manual change | Reason |
|------|---------------------------|---------------|--------|
| results/tables/[example_file].csv | HPC-derived output / processed table | Renamed columns to more descriptive labels | Improve clarity and consistency |
| results/supplementary_tables/[example_file].xlsx | Final formatted supplementary table | Adjusted title and reordered columns | Match manuscript presentation |
| results/tables/[example_matrix].csv | Derived analysis matrix | Standardized category labels | Ensure consistency across figures |

---

## Notes for Updating

For each manually edited file, record:

1. **File name**
2. **Where it came from**
3. **What was changed**
4. **Why the change was made**

This helps distinguish:
- direct computational outputs
- manually curated presentation-ready tables

---

## Current Status

This log should be updated whenever:

- a processed table is renamed
- a supplementary table is reformatted
- a category label is standardized manually
- a final manuscript table differs from the raw exported analysis table
