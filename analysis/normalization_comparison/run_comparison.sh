#!/usr/bin/env bash
# Regenerate the normalization comparison in figures/normalization_comparison.
#
# Runs the example pipeline with the original code (branch `master`) and the
# corrected code (branch `proposed-fixes`), estimates the detection limit by
# simulation, re-renders the corrected heatmaps with the colour floor at that
# limit, and assembles the comparison figures.
#
# Requirements: the htpamdaenv conda environment (environment.yml), local
# branches `master` and `proposed-fixes`, and macOS (`sips` is used to
# rasterize the pipeline's PDF heatmaps).
#
# Usage, from anywhere in the repository:
#   conda activate htpamdaenv
#   bash analysis/normalization_comparison/run_comparison.sh [work_dir]
set -euo pipefail

REPO=$(git rev-parse --show-toplevel)
HERE="$REPO/analysis/normalization_comparison"
OUT="$REPO/figures/normalization_comparison"
if [ -n "${1:-}" ]; then
    mkdir -p "$1"
    WORK=$(cd "$1" && pwd)
else
    WORK=$(mktemp -d)
fi
echo "Working directory: $WORK"

cleanup() {
    git -C "$REPO" worktree remove --force "$WORK/original" 2>/dev/null || true
    git -C "$REPO" worktree remove --force "$WORK/corrected" 2>/dev/null || true
}
trap cleanup EXIT

git -C "$REPO" worktree add --detach "$WORK/original" master
git -C "$REPO" worktree add --detach "$WORK/corrected" proposed-fixes

for version in original corrected; do
    echo "== Running pipeline: $version"
    (cd "$WORK/$version/code" && python library_QC.py && python PAMDA_complete.py)
done

echo "== Simulating null background"
python "$HERE/simulate_null_background.py" "$WORK/original" "$WORK/corrected" "$WORK/null_simulation"
LOD=$(cat "$WORK/null_simulation/detection_limit.txt")

echo "== Rendering corrected heatmaps with colour floor $LOD"
(cd "$WORK/corrected" && python - "$LOD" <<'EOF'
import sys
sys.path.insert(0, 'code')
import PAMDA
PAMDA.rate2heatmap('lod_floor', 'barcode_csv/example_PAMDA_barcodes.csv', 4, 0,
                   heatmap_fixed_min=float(sys.argv[1]), heatmap_fixed_max=-1.5,
                   input_csv='output/example_PAMDA_data/PAM_start_0_length_4/PAMDA_3_rates.csv')
EOF
)

echo "== Rasterizing heatmaps"
render() {
    mkdir -p "$WORK/png/$2"
    for pdf in "$1"/PAMDA_HEATMAP_*.pdf; do
        sips -s format png -s dpiHeight 150 -s dpiWidth 150 "$pdf" \
            --out "$WORK/png/$2/$(basename "${pdf%.pdf}").png" >/dev/null
    done
}
render "$WORK/original/figures/example_PAMDA_data/PAM_start_0_length_4" original
render "$WORK/corrected/figures/example_PAMDA_data/PAM_start_0_length_4" corrected
render "$WORK/corrected/figures/lod_floor/PAM_start_0_length_4" corrected_lod

echo "== Assembling comparison figures"
python "$HERE/make_comparison_figures.py" "$REPO" "$WORK/png" "$OUT" "$LOD"
cp "$WORK/original/output/example_PAMDA_data/PAM_start_0_length_4/PAMDA_3_rates.csv" "$OUT/rates_original.csv"
cp "$WORK/corrected/output/example_PAMDA_data/PAM_start_0_length_4/PAMDA_3_rates.csv" "$OUT/rates_corrected.csv"
cp "$WORK/null_simulation/summary.txt" "$OUT/null_simulation_summary.txt"
echo "Done. Figures written to $OUT"
