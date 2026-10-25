#!/bin/bash -l
# ==============================================================================
# VERVE step 00: submit preprocessing
#
# Usage: bash 00.run_preprocess.sh
# Run from the pipeline directory (the one containing config.sh).
# ==============================================================================

source config.sh

# ── Sanity checks ─────────────────────────────────────────────────────────────
if [[ -z "$SAMPLE_SHEET" || -z "$OUT_DIR" ]]; then
    echo "ERROR: SAMPLE_SHEET and OUT_DIR must be set in config.sh"
    exit 1
fi

if [[ ! -f "$SAMPLE_SHEET" ]]; then
    echo "ERROR: Sample sheet not found: $SAMPLE_SHEET"
    exit 1
fi

for VAR in MODULE_PATH_QIB PKG_FASTP PKG_MULTIQC; do
    if [[ -z "${!VAR}" ]]; then
        echo "ERROR: $VAR is empty in config.sh"
        exit 1
    fi
done

if [[ "$RUN_FASTQC" == true && -z "$PKG_FASTQC" ]]; then
    echo "ERROR: RUN_FASTQC=true requires PKG_FASTQC to be set in config.sh"
    exit 1
fi

if [[ ! -f "modules/00a_preprocess.sh" || ! -f "modules/00b_preprocess_multiqc.sh" ]]; then
    echo "ERROR: module scripts not found. Run this from the pipeline directory."
    exit 1
fi

# Strip Windows carriage returns once, here, so every downstream script inherits
# a clean sample sheet
sed -i 's/\r//' "$SAMPLE_SHEET"

[[ ! -d "$LOG_DIR" ]] && mkdir -p "$LOG_DIR"
[[ ! -d "$OUT_DIR/00_preprocess" ]] && mkdir -p "$OUT_DIR/00_preprocess"

# Count samples, excluding the header and any trailing blank lines
N_SAMPLES=$(tail -n +2 "$SAMPLE_SHEET" | grep -c . )
if [[ "$N_SAMPLES" -lt 1 ]]; then
    echo "ERROR: no samples found in $SAMPLE_SHEET"
    exit 1
fi

echo "=== VERVE step 00: preprocessing ==="
echo "  Sample sheet : $SAMPLE_SHEET"
echo "  Samples      : $N_SAMPLES"
echo "  Output       : $OUT_DIR/00_preprocess"
echo "  Logs         : $LOG_DIR"
echo "  fastp        : $PKG_FASTP"
echo "  FastQC       : $([[ "$RUN_FASTQC" == true ]] && echo "$PKG_FASTQC" || echo "skipped")"
echo "  MultiQC      : $PKG_MULTIQC"
echo ""

# ── Submit array job ──────────────────────────────────────────────────────────
JOB_ID=$(sbatch \
    --job-name=00a_preprocess \
    --partition="$PARTITION_PREPROCESS" \
    --time="$TIME_PREPROCESS" \
    --cpus-per-task="$CPUS_PREPROCESS" \
    --mem-per-cpu="$MEM_PER_CPU_PREPROCESS" \
    --array=1-${N_SAMPLES} \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    modules/00a_preprocess.sh | awk '{print $NF}')
echo "Array job submitted : $JOB_ID"

# ── Submit MultiQC merge job ──────────────────────────────────────────────────
MERGE_ID=$(sbatch \
    --job-name=00b_preprocess_multiqc \
    --partition="$PARTITION_PREPROCESS_MERGE" \
    --time="$TIME_PREPROCESS_MERGE" \
    --cpus-per-task="$CPUS_PREPROCESS_MERGE" \
    --mem-per-cpu="$MEM_PER_CPU_PREPROCESS_MERGE" \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    --dependency=afterok:${JOB_ID} \
    modules/00b_preprocess_multiqc.sh | awk '{print $NF}')
echo "MultiQC job submitted : $MERGE_ID (runs after $JOB_ID completes)"

echo ""
echo "Next step: bash 00.check_preprocess.sh once all jobs complete"
