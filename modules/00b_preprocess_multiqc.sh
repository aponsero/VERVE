#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
# ==============================================================================
# VERVE module 00b: aggregate preprocessing report (single job)
#
# MultiQC over the fastp JSON reports and, if produced, the FastQC output.
# Submitted by 00.run_preprocess.sh with --dependency=afterok on the array job.
# ==============================================================================

source "${SLURM_SUBMIT_DIR:-.}/config.sh"

module use "$MODULE_PATH_QIB"
module load "$PKG_MULTIQC"

echo "Job started"; hostname; date
echo ""

OUT_STEP="$OUT_DIR/00_preprocess"
MULTIQC_DIR="$OUT_STEP/multiqc"
[[ ! -d "$MULTIQC_DIR" ]] && mkdir -p "$MULTIQC_DIR"

# ── Sanity check on expected inputs ───────────────────────────────────────────
N_SAMPLES=$(tail -n +2 "$SAMPLE_SHEET" | grep -c . )
N_JSON=$(find "$OUT_STEP/fastp" -name "*.fastp.json" 2>/dev/null | wc -l)

echo "Samples in sample sheet : $N_SAMPLES"
echo "fastp JSON reports found: $N_JSON"

if [[ "$N_JSON" -eq 0 ]]; then
    echo "ERROR: no fastp JSON reports found in $OUT_STEP/fastp"
    echo "       the array job did not produce output; nothing to aggregate"
    exit 1
fi

if [[ "$N_JSON" -lt "$N_SAMPLES" ]]; then
    echo "WARNING: fewer reports than samples. MultiQC will run on what exists,"
    echo "         but the report will be incomplete. Check the array job logs"
    echo "         in $LOG_DIR and rerun the missing tasks."
fi
echo ""

# ── MultiQC ───────────────────────────────────────────────────────────────────
echo "Running MultiQC"
multiqc "$OUT_STEP" \
    --outdir "$MULTIQC_DIR" \
    --filename "multiqc_report.html" \
    --title "$MULTIQC_TITLE" \
    --force \
    $MULTIQC_EXTRA_ARGS

MULTIQC_STATUS=$?
if [[ "$MULTIQC_STATUS" -ne 0 ]]; then
    echo "ERROR: MultiQC exited with status $MULTIQC_STATUS"
    exit "$MULTIQC_STATUS"
fi

echo ""
echo "Report written to: $MULTIQC_DIR/multiqc_report.html"
echo "Job done"; date
