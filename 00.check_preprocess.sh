#!/bin/bash -l
# ==============================================================================
# VERVE step 00: check preprocessing outputs
#
# Usage: bash 00.check_preprocess.sh
# Run from the pipeline directory once the array and MultiQC jobs have finished.
# ==============================================================================

source config.sh

OUT_STEP="$OUT_DIR/00_preprocess"

echo "=== Preprocessing check ==="
echo ""
echo "Configuration"
echo "  Sample sheet : $SAMPLE_SHEET"
echo "  Output       : $OUT_STEP"
echo "  Logs         : $LOG_DIR"
echo "  FastQC       : $RUN_FASTQC"
echo "  fastp        : $PKG_FASTP"
echo "  MultiQC      : $PKG_MULTIQC"
echo "  fastp params : Q>=$FASTP_QUALIFIED_QUALITY, unqualified<=${FASTP_UNQUALIFIED_PERCENT}%, N<=$FASTP_N_BASE_LIMIT, length>=$FASTP_LENGTH_REQUIRED"
echo ""

PASS=0; FAIL=0; MISSING=0; EMPTY=0; BLANK_EMPTY=0

while IFS=$'\t' read -r SAMPLE_ID READ1 READ2 BODY_SITE GROUP BATCH; do
    [[ "$SAMPLE_ID" == "sample_id" ]] && continue
    [[ -z "$SAMPLE_ID" ]] && continue
    SAMPLE_ID=$(echo "$SAMPLE_ID" | tr -d '\r')
    GROUP=$(echo "$GROUP" | tr -d '\r')
    SAMPLE_FAIL=false
    SAMPLE_NOTE=""

    # A negative control with no amplifiable material can legitimately yield an
    # empty trimmed FASTQ. That is reported, not counted as a pipeline failure.
    IS_BLANK=false
    [[ "$GROUP" == blank* ]] && IS_BLANK=true

    FILES=( \
        "$OUT_STEP/trimmed/${SAMPLE_ID}_R1.trimmed.fastq.gz" \
        "$OUT_STEP/trimmed/${SAMPLE_ID}_R2.trimmed.fastq.gz" \
        "$OUT_STEP/fastp/${SAMPLE_ID}.fastp.json" \
        "$OUT_STEP/fastp/${SAMPLE_ID}.fastp.html" )

    if [[ "$RUN_FASTQC" == true ]]; then
        FILES+=( \
            "$OUT_STEP/fastqc/${SAMPLE_ID}_R1.trimmed_fastqc.zip" \
            "$OUT_STEP/fastqc/${SAMPLE_ID}_R2.trimmed_fastqc.zip" )
    fi

    for FILE in "${FILES[@]}"; do
        if [[ ! -f "$FILE" ]]; then
            echo "MISSING    : $SAMPLE_ID -- $(basename "$FILE")"
            ((MISSING++)); SAMPLE_FAIL=true
        elif [[ ! -s "$FILE" ]]; then
            if [[ "$IS_BLANK" == true ]]; then
                echo "NO READS   : $SAMPLE_ID -- $(basename "$FILE") (blank control, expected)"
                ((BLANK_EMPTY++)); SAMPLE_NOTE=" (blank, no reads)"
            else
                echo "EMPTY      : $SAMPLE_ID -- $(basename "$FILE")"
                ((EMPTY++)); SAMPLE_FAIL=true
            fi
        fi
    done

    if [[ "$SAMPLE_FAIL" == false ]]; then
        echo "OK         : $SAMPLE_ID$SAMPLE_NOTE"
        ((PASS++))
    else
        ((FAIL++))
    fi

done < "$SAMPLE_SHEET"

# ── Aggregated output ─────────────────────────────────────────────────────────
echo ""
echo "=== Aggregated output ==="
AGG_FAIL=false

FILE="$OUT_STEP/multiqc/multiqc_report.html"
if [[ ! -f "$FILE" ]]; then
    echo "MISSING    : multiqc_report.html -- check whether the MultiQC job ran"
    AGG_FAIL=true
elif [[ ! -s "$FILE" ]]; then
    echo "EMPTY      : multiqc_report.html"
    AGG_FAIL=true
else
    echo "OK         : multiqc_report.html"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Summary ==="
echo "  Total : $(( PASS + FAIL )) | Passed : $PASS | Failed : $FAIL"
echo "  (missing: $MISSING | empty: $EMPTY | blank controls with no reads: $BLANK_EMPTY)"
echo ""

if [[ "$FAIL" -gt 0 || "$AGG_FAIL" == true ]]; then
    echo "WARNING: some outputs are missing or empty."
    echo "         Check the logs in $LOG_DIR before proceeding."
    echo "         To resubmit a single failed task:"
    echo "           sbatch --array=<task_id> ... modules/00a_preprocess.sh"
else
    echo "All checks passed."
    echo "Open $OUT_STEP/multiqc/multiqc_report.html and inspect the reports"
    echo "before moving on to module 01."
fi
