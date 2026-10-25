#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
# ==============================================================================
# VERVE module 03c: merge Telescope reports (single job)
#
#   - locus-level and family-level count matrices, plus a per-sample stats
#     table (scripts/merge_telescope.py)
#   - MultiQC over the bowtie2 logs, for the ERV-arm alignment rates
#
# Counts are RAW. Normalisation uses the module 02 size factors and happens in
# R, so that a change of model does not mean rerunning the pipeline.
# ==============================================================================

source "${SLURM_SUBMIT_DIR:-.}/config.sh"

module use "$MODULE_PATH_QIB"
module load "$PKG_MULTIQC"

echo "Job started"; hostname; date
echo ""

OUT_STEP="$OUT_DIR/03_erv"
COUNTS_DIR="$OUT_STEP/counts"
MULTIQC_DIR="$OUT_STEP/multiqc"
for D in "$COUNTS_DIR" "$MULTIQC_DIR"; do
    [[ ! -d "$D" ]] && mkdir -p "$D"
done

MERGE_SCRIPT="${SLURM_SUBMIT_DIR:-.}/scripts/merge_telescope.py"
if [[ ! -f "$MERGE_SCRIPT" ]]; then
    echo "ERROR: $MERGE_SCRIPT not found"
    exit 1
fi

echo "=== Merging Telescope reports ==="
python3 "$MERGE_SCRIPT" \
    --sample-sheet "$SAMPLE_SHEET" \
    --telescope-dir "$OUT_STEP/telescope" \
    --out-prefix "$COUNTS_DIR/erv_counts"

if [[ $? -ne 0 ]]; then
    echo "ERROR: merge failed"
    exit 1
fi
echo ""

echo "=== MultiQC (bowtie2 logs) ==="
multiqc "$OUT_STEP/align" \
    --outdir "$MULTIQC_DIR" \
    --filename "multiqc_report.html" \
    --title "VERVE module 03 ERV alignment" \
    --force \
    $MULTIQC_EXTRA_ARGS

if [[ $? -ne 0 ]]; then
    echo "ERROR: MultiQC failed"
    exit 1
fi

echo ""
echo "Job done"; date
