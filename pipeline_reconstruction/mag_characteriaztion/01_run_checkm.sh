#!/usr/bin/env bash
set -euo pipefail

MAG_DIR="${1:-data/mags}"
OUT_DIR="${2:-results/mag_characterization/checkm}"
CORES="${3:-16}"
EXT="${4:-fa}"

mkdir -p "$OUT_DIR"

# Output format consistent with MuDoGeR: tab table with header
checkm lineage_wf \
  -t "$CORES" \
  --reduced_tree \
  --tab_table \
  -x "$EXT" \
  -f "$OUT_DIR/output_checkm.tsv" \
  "$MAG_DIR" \
  "$OUT_DIR"

echo "Done: $OUT_DIR/output_checkm.tsv"
