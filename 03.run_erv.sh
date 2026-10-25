#!/bin/bash -l
# ==============================================================================
# VERVE step 03: submit ERV quantification
#
# Usage: bash 03.run_erv.sh
# Run from the pipeline directory (the one containing config.sh).
#
# 03a bowtie2 (array) -> 03b telescope assign (array) -> 03c merge.
# Both arrays skip samples that are already done, so the script is safe to
# resubmit after a partial failure.
# ==============================================================================

source config.sh

ERR=false

for VAR in SAMPLE_SHEET OUT_DIR MODULE_PATH_QIB \
           PKG_BOWTIE2 PKG_TELESCOPE PKG_SAMTOOLS PKG_MULTIQC \
           BOWTIE2_INDEX ERV_GTF \
           PARTITION_ERV_ALIGN PARTITION_ERV_ASSIGN PARTITION_ERV_MERGE; do
    if [[ -z "${!VAR}" ]]; then
        echo "ERROR: $VAR is empty in config.sh"; ERR=true
    fi
done

[[ -n "$SAMPLE_SHEET" && ! -f "$SAMPLE_SHEET" ]] && { echo "ERROR: sample sheet not found: $SAMPLE_SHEET"; ERR=true; }
[[ -n "$ERV_GTF" && ! -f "$ERV_GTF" ]] && { echo "ERROR: ERV annotation not found: $ERV_GTF"; ERR=true; }
[[ -n "$BOWTIE2_INDEX" && ! -s "${BOWTIE2_INDEX}.1.bt2" && ! -s "${BOWTIE2_INDEX}.1.bt2l" ]] && \
    { echo "ERROR: no bowtie2 index at prefix: $BOWTIE2_INDEX"; ERR=true; }

case "$TELESCOPE_STRANDED_MODE" in
    RF|FR|R|F|None) ;;
    *) echo "ERROR: TELESCOPE_STRANDED_MODE must be RF, FR, R, F or None (got '$TELESCOPE_STRANDED_MODE')"; ERR=true ;;
esac

for S in 03a_erv_align 03b_telescope 03c_erv_merge; do
    [[ ! -f "modules/${S}.sh" ]] && { echo "ERROR: modules/${S}.sh not found. Run from the pipeline directory."; ERR=true; }
done
[[ ! -f "scripts/merge_telescope.py" ]] && { echo "ERROR: scripts/merge_telescope.py not found"; ERR=true; }

[[ "$ERR" == true ]] && exit 1

sed -i 's/\r//' "$SAMPLE_SHEET"
[[ ! -d "$LOG_DIR" ]] && mkdir -p "$LOG_DIR"
[[ ! -d "$OUT_DIR/03_erv" ]] && mkdir -p "$OUT_DIR/03_erv"

N_SAMPLES=$(tail -n +2 "$SAMPLE_SHEET" | grep -c . )
if [[ "$N_SAMPLES" -lt 1 ]]; then
    echo "ERROR: no samples found in $SAMPLE_SHEET"
    exit 1
fi

echo "=== VERVE step 03: ERV quantification ==="
echo "  Sample sheet  : $SAMPLE_SHEET"
echo "  Samples       : $N_SAMPLES"
echo "  Output        : $OUT_DIR/03_erv"
echo "  Annotation    : $ERV_GTF"
echo "  bowtie2       : $PKG_BOWTIE2 ($BOWTIE2_PRESET -k $BOWTIE2_K --score-min $BOWTIE2_SCORE_MIN)"
echo "  Telescope     : $PKG_TELESCOPE (stranded mode $TELESCOPE_STRANDED_MODE)"
echo "  Delete BAMs   : $ERV_DELETE_BAM_AFTER_ASSIGN"
echo ""

ALIGN_ID=$(sbatch \
    --job-name=03a_erv_align \
    --partition="$PARTITION_ERV_ALIGN" \
    --time="$TIME_ERV_ALIGN" \
    --cpus-per-task="$CPUS_ERV_ALIGN" \
    --mem-per-cpu="$MEM_PER_CPU_ERV_ALIGN" \
    --array=1-${N_SAMPLES} \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    modules/03a_erv_align.sh | awk '{print $NF}')
echo "03a bowtie2 array submitted   : $ALIGN_ID"

ASSIGN_ID=$(sbatch \
    --job-name=03b_telescope \
    --partition="$PARTITION_ERV_ASSIGN" \
    --time="$TIME_ERV_ASSIGN" \
    --cpus-per-task="$CPUS_ERV_ASSIGN" \
    --mem="$MEM_ERV_ASSIGN" \
    --array=1-${N_SAMPLES} \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    --dependency=afterok:${ALIGN_ID} \
    modules/03b_telescope.sh | awk '{print $NF}')
echo "03b telescope array submitted : $ASSIGN_ID (after $ALIGN_ID)"

MERGE_ID=$(sbatch \
    --job-name=03c_erv_merge \
    --partition="$PARTITION_ERV_MERGE" \
    --time="$TIME_ERV_MERGE" \
    --cpus-per-task="$CPUS_ERV_MERGE" \
    --mem-per-cpu="$MEM_PER_CPU_ERV_MERGE" \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    --dependency=afterok:${ASSIGN_ID} \
    modules/03c_erv_merge.sh | awk '{print $NF}')
echo "03c merge submitted           : $MERGE_ID (after $ASSIGN_ID)"

echo ""
echo "Next step: bash 03.check_erv.sh once all jobs complete"
