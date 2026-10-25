#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
# ==============================================================================
# VERVE module 04b: merge viral metrics (single job)
#
#   - per-accession and per-taxid tables, segments summed, with RPM against
#     both denominators (scripts/merge_viral.py)
#   - MultiQC over the bowtie2 logs, for the viral-arm alignment rates
#
# No positivity call is made here. Counts, breadth and covered windows are
# written out; the positivity rule is set in the report once the pilot numbers
# are in hand, the same way the module 01 thresholds were.
# ==============================================================================

source "${SLURM_SUBMIT_DIR:-.}/config.sh"

module use "$MODULE_PATH_QIB"
module load "$PKG_MULTIQC"

echo "Job started"; hostname; date
echo ""

OUT_STEP="$OUT_DIR/04_viral"
COUNTS_DIR="$OUT_STEP/counts"
MULTIQC_DIR="$OUT_STEP/multiqc"
for D in "$COUNTS_DIR" "$MULTIQC_DIR"; do
    [[ ! -d "$D" ]] && mkdir -p "$D"
done

MERGE_SCRIPT="${SLURM_SUBMIT_DIR:-.}/scripts/merge_viral.py"
if [[ ! -f "$MERGE_SCRIPT" ]]; then
    echo "ERROR: $MERGE_SCRIPT not found"
    exit 1
fi

echo "=== Merging per-sample viral metrics ==="
python3 "$MERGE_SCRIPT" \
    --sample-sheet "$SAMPLE_SHEET" \
    --out-dir "$OUT_DIR" \
    --taxa "$VIRAL_TAXA" \
    --out-prefix "$COUNTS_DIR/viral" \
    --window-size "$VIRAL_WINDOW_SIZE"

if [[ $? -ne 0 ]]; then
    echo "ERROR: merge failed"
    exit 1
fi
echo ""

echo "=== MultiQC (bowtie2 logs) ==="
multiqc "$OUT_STEP/align" \
    --outdir "$MULTIQC_DIR" \
    --filename "multiqc_report.html" \
    --title "VERVE module 04 viral mapping" \
    --force \
    $MULTIQC_EXTRA_ARGS

if [[ $? -ne 0 ]]; then
    echo "ERROR: MultiQC failed"
    exit 1
fi

echo ""
echo "Job done"; date
