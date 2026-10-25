#!/bin/bash -l
# ==============================================================================
# VERVE step 04: submit RNA viral detection
#
# Usage: bash 04.run_viral.sh
# Run from the pipeline directory (the one containing config.sh).
#
# 04a map and count (array) -> 04b merge. The array skips samples that already
# have metrics, so the script is safe to resubmit after a partial failure.
# ==============================================================================

source config.sh

ERR=false

for VAR in SAMPLE_SHEET OUT_DIR MODULE_PATH_QIB \
           PKG_BOWTIE2 PKG_SAMTOOLS PKG_MULTIQC \
           VIRAL_INDEX VIRAL_TAXA \
           PARTITION_VIRAL PARTITION_VIRAL_MERGE; do
    if [[ -z "${!VAR}" ]]; then
        echo "ERROR: $VAR is empty in config.sh"; ERR=true
    fi
done

[[ -n "$SAMPLE_SHEET" && ! -f "$SAMPLE_SHEET" ]] && { echo "ERROR: sample sheet not found: $SAMPLE_SHEET"; ERR=true; }
[[ -n "$VIRAL_TAXA" && ! -f "$VIRAL_TAXA" ]] && { echo "ERROR: taxon table not found: $VIRAL_TAXA"; ERR=true; }
[[ -n "$VIRAL_INDEX" && ! -s "${VIRAL_INDEX}.1.bt2" && ! -s "${VIRAL_INDEX}.1.bt2l" ]] && \
    { echo "ERROR: no bowtie2 index at prefix: $VIRAL_INDEX"; ERR=true; }

for S in 04a_viral_map 04b_viral_merge; do
    [[ ! -f "modules/${S}.sh" ]] && { echo "ERROR: modules/${S}.sh not found. Run from the pipeline directory."; ERR=true; }
done
[[ ! -f "scripts/merge_viral.py" ]] && { echo "ERROR: scripts/merge_viral.py not found"; ERR=true; }

[[ "$ERR" == true ]] && exit 1

sed -i 's/\r//' "$SAMPLE_SHEET"
[[ ! -d "$LOG_DIR" ]] && mkdir -p "$LOG_DIR"
[[ ! -d "$OUT_DIR/04_viral" ]] && mkdir -p "$OUT_DIR/04_viral"

N_SAMPLES=$(tail -n +2 "$SAMPLE_SHEET" | grep -c . )
if [[ "$N_SAMPLES" -lt 1 ]]; then
    echo "ERROR: no samples found in $SAMPLE_SHEET"
    exit 1
fi

N_VIRAL_SEQS=$(( $(wc -l < "$VIRAL_TAXA") - 1 ))

echo "=== VERVE step 04: RNA viral detection ==="
echo "  Sample sheet   : $SAMPLE_SHEET"
echo "  Samples        : $N_SAMPLES"
echo "  Output         : $OUT_DIR/04_viral"
echo "  Index          : $VIRAL_INDEX"
echo "  Taxon table    : $VIRAL_TAXA ($N_VIRAL_SEQS sequences)"
echo "  bowtie2        : $PKG_BOWTIE2 ($BOWTIE2_VIRAL_PRESET, best hit)"
echo "  High-conf MAPQ : $VIRAL_MIN_MAPQ"
echo "  Window size    : $VIRAL_WINDOW_SIZE bp"
echo "  Keep BAMs      : $VIRAL_KEEP_BAM"
echo ""

MAP_ID=$(sbatch \
    --job-name=04a_viral_map \
    --partition="$PARTITION_VIRAL" \
    --time="$TIME_VIRAL" \
    --cpus-per-task="$CPUS_VIRAL" \
    --mem-per-cpu="$MEM_PER_CPU_VIRAL" \
    --array=1-${N_SAMPLES} \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    modules/04a_viral_map.sh | awk '{print $NF}')
echo "04a map array submitted : $MAP_ID"

MERGE_ID=$(sbatch \
    --job-name=04b_viral_merge \
    --partition="$PARTITION_VIRAL_MERGE" \
    --time="$TIME_VIRAL_MERGE" \
    --cpus-per-task="$CPUS_VIRAL_MERGE" \
    --mem-per-cpu="$MEM_PER_CPU_VIRAL_MERGE" \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    --dependency=afterok:${MAP_ID} \
    modules/04b_viral_merge.sh | awk '{print $NF}')
echo "04b merge submitted     : $MERGE_ID (after $MAP_ID)"

echo ""
echo "Next step: bash 04.check_viral.sh once all jobs complete"

