#!/usr/bin/env bash
set -euo pipefail

MAG_DIR="${1:-data/mags}"
OUT_DIR="${2:-results/mag_characterization/bbtools}"
EXT="${3:-fa}"

mkdir -p "$OUT_DIR"

# Requires BBTools in PATH (statswrapper.sh)
statswrapper.sh "$MAG_DIR"/*."$EXT" > "$OUT_DIR/bbtools.tsv"

echo "Done: $OUT_DIR/bbtools.tsv"
