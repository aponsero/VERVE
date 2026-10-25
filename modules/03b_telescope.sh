#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
# ==============================================================================
# VERVE module 03b: Telescope EM reassignment (array job)
#
# telescope assign reads the permissive alignment from 03a and reassigns
# ambiguously mapped fragments to the most probable locus, rather than
# discarding them.
#
# Single core: Telescope states that multiple cores are not supported. Memory,
# not CPU, is the binding constraint.
# ==============================================================================

source "${SLURM_SUBMIT_DIR:-.}/config.sh"

module use "$MODULE_PATH_QIB"
module load "$PKG_TELESCOPE"

echo "Job started"; hostname; date
echo ""

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_SHEET" | tr -d '\r')
SAMPLE_ID=$(echo "$LINE" | awk '{print $1}')

if [[ -z "$SAMPLE_ID" ]]; then
    echo "ERROR: could not parse line $((SLURM_ARRAY_TASK_ID + 1)) of $SAMPLE_SHEET"
    exit 1
fi

BAM="$OUT_DIR/03_erv/align/${SAMPLE_ID}.erv.bam"
if [[ ! -s "$BAM" ]]; then
    echo "ERROR: alignment not found for $SAMPLE_ID: $BAM"
    echo "       has 03a completed for this sample?"
    exit 1
fi

SAMPLE_OUT="$OUT_DIR/03_erv/telescope/${SAMPLE_ID}"
[[ ! -d "$SAMPLE_OUT" ]] && mkdir -p "$SAMPLE_OUT"

# Telescope names its outputs <exp_tag>-<file>. The exact report name differs
# between versions (TE_counts.tsv in recent docs, telescope_report.tsv in
# earlier ones), so the check and merge steps glob for it rather than assume.
if ls "$SAMPLE_OUT"/${SAMPLE_ID}-*counts.tsv "$SAMPLE_OUT"/${SAMPLE_ID}-*report.tsv > /dev/null 2>&1; then
    echo "Telescope output already exists for $SAMPLE_ID, skipping"
    echo "Job done"; date
    exit 0
fi

echo "Running telescope assign : $SAMPLE_ID"
echo "  annotation             : $ERV_GTF"
echo "  stranded mode          : $TELESCOPE_STRANDED_MODE"
echo "  attribute              : $TELESCOPE_ATTRIBUTE"
echo ""

telescope assign \
    "$BAM" \
    "$ERV_GTF" \
    --outdir "$SAMPLE_OUT" \
    --exp_tag "$SAMPLE_ID" \
    --attribute "$TELESCOPE_ATTRIBUTE" \
    --stranded_mode "$TELESCOPE_STRANDED_MODE" \
    --ncpu 1 \
    $TELESCOPE_EXTRA_ARGS

STATUS=$?
if [[ "$STATUS" -ne 0 ]]; then
    echo "ERROR: telescope assign exited with status $STATUS for $SAMPLE_ID"
    exit "$STATUS"
fi

echo ""
echo "Outputs in $SAMPLE_OUT:"
ls -la "$SAMPLE_OUT"

if [[ "$ERV_DELETE_BAM_AFTER_ASSIGN" == true ]]; then
    echo ""
    echo "Deleting $BAM (ERV_DELETE_BAM_AFTER_ASSIGN is true)"
    rm -f "$BAM"
fi

echo ""
echo "Job done"; date
