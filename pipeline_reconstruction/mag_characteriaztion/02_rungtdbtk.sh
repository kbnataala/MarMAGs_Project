#!/usr/bin/env bash
set -euo pipefail

MAG_DIR="${1:-data/mags}"
OUT_DIR="${2:-results/mag_characterization/gtdbtk}"
CORES="${3:-16}"
EXT="${4:-fa}"

# Must be set by the user (or module): export GTDBTK_DATA_PATH=/path/to/release220
: "${GTDBTK_DATA_PATH:?Need GTDBTK_DATA_PATH set (e.g., GTDB release 220)}"

mkdir -p "$OUT_DIR"

gtdbtk classify_wf \
  --skip_ani_screen \
  --extension "$EXT" \
  --cpus "$CORES" \
  --genome_dir "$MAG_DIR" \
  --out_dir "$OUT_DIR"

# Merge bac120 + ar53 (MuDoGeR-style)
if ls "$OUT_DIR"/gtdbtk.*summary.tsv >/dev/null 2>&1; then
  awk 'FNR==1 && NR!=1 {next} {print}' "$OUT_DIR"/gtdbtk.*summary.tsv > "$OUT_DIR/gtdbtk_result.tsv"
  echo "Done: $OUT_DIR/gtdbtk_result.tsv"
else
  echo "ERROR: no GTDB-Tk summary files found in $OUT_DIR" >&2
  exit 1
fi
