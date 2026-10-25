#!/bin/bash -l
# ==============================================================================
# VERVE step 01: check host alignment outputs
#
# Usage: bash 01.check_host.sh
# Run from the pipeline directory once the submitted jobs have finished.
#
# Always checks 01a and 01b outputs and prints the strandedness table.
# Checks 01c and 01d outputs only if STRANDEDNESS is set in config.sh.
# ==============================================================================

source config.sh

OUT_STEP="$OUT_DIR/01_host"

echo "=== Host alignment check ==="
echo ""
echo "Configuration"
echo "  Sample sheet : $SAMPLE_SHEET"
echo "  Output       : $OUT_STEP"
echo "  Logs         : $LOG_DIR"
echo "  STAR index   : $STAR_INDEX"
echo "  GTF          : $GTF"
echo "  Strandedness : ${STRANDEDNESS:-not set}"
echo "  STAR         : $PKG_STAR"
echo "  SortMeRNA    : $PKG_SORTMERNA"
echo "  Picard       : $PKG_PICARD"
echo "  featureCounts: $PKG_SUBREAD"
echo ""

PHASE2=false
[[ -n "$STRANDEDNESS" ]] && PHASE2=true

PASS=0; FAIL=0; MISSING=0; EMPTY=0; BLANK_EMPTY=0

check_file () {
    # $1 file, $2 sample, $3 is_blank. Sets SAMPLE_FAIL and SAMPLE_NOTE.
    local FILE="$1" SID="$2" BLANK="$3"
    if [[ ! -f "$FILE" ]]; then
        echo "MISSING    : $SID -- $(basename "$FILE")"
        ((MISSING++)); SAMPLE_FAIL=true
    elif [[ ! -s "$FILE" ]]; then
        if [[ "$BLANK" == true ]]; then
            echo "NO READS   : $SID -- $(basename "$FILE") (blank control, expected)"
            ((BLANK_EMPTY++)); SAMPLE_NOTE=" (blank, no reads)"
        else
            echo "EMPTY      : $SID -- $(basename "$FILE")"
            ((EMPTY++)); SAMPLE_FAIL=true
        fi
    fi
}

STRAND_TABLE=""

while IFS=$'\t' read -r SAMPLE_ID READ1 READ2 BODY_SITE GROUP BATCH; do
    [[ "$SAMPLE_ID" == "sample_id" || -z "$SAMPLE_ID" ]] && continue
    SAMPLE_ID=$(echo "$SAMPLE_ID" | tr -d '\r')
    GROUP=$(echo "$GROUP" | tr -d '\r')
    SAMPLE_FAIL=false
    SAMPLE_NOTE=""
    IS_BLANK=false
    [[ "$GROUP" == blank* ]] && IS_BLANK=true

    PREFIX="$OUT_STEP/align/${SAMPLE_ID}."
    FILES=( \
        "${PREFIX}Aligned.sortedByCoord.out.bam" \
        "${PREFIX}Aligned.sortedByCoord.out.bam.bai" \
        "${PREFIX}Log.final.out" \
        "${PREFIX}ReadsPerGene.out.tab" \
        "$OUT_STEP/unmapped/${SAMPLE_ID}_R1.unmapped.fastq.gz" \
        "$OUT_STEP/unmapped/${SAMPLE_ID}_R2.unmapped.fastq.gz" \
        "$OUT_STEP/rrna/${SAMPLE_ID}.log" )

    if [[ "$PHASE2" == true ]]; then
        FILES+=( \
            "$OUT_STEP/markdup/${SAMPLE_ID}.markdup.bam" \
            "$OUT_STEP/bamqc/${SAMPLE_ID}.markdup_metrics.txt" \
            "$OUT_STEP/bamqc/${SAMPLE_ID}.rnaseq_metrics.txt" \
            "$OUT_STEP/bamqc/${SAMPLE_ID}.flagstat" \
            "$OUT_STEP/bamqc/${SAMPLE_ID}.idxstats" )
    fi

    for FILE in "${FILES[@]}"; do
        check_file "$FILE" "$SAMPLE_ID" "$IS_BLANK"
    done

    if [[ "$SAMPLE_FAIL" == false ]]; then
        echo "OK         : $SAMPLE_ID$SAMPLE_NOTE"
        ((PASS++))
    else
        ((FAIL++))
    fi

    # Strandedness totals: gene rows only (the first 4 lines are N_unmapped,
    # N_multimapping, N_noFeature, N_ambiguous)
    RPG="${PREFIX}ReadsPerGene.out.tab"
    if [[ -s "$RPG" ]]; then
        ROW=$(awk -v s="$SAMPLE_ID" 'NR>4 {u+=$2; f+=$3; r+=$4}
              END {fr = (f+r > 0) ? sprintf("%.3f", f/(f+r)) : "NA";
                   printf "%-24s %14d %14d %14d %10s", s, u, f, r, fr}' "$RPG")
        STRAND_TABLE+="$ROW"$'\n'
    fi

done < "$SAMPLE_SHEET"

# ── Strandedness table ────────────────────────────────────────────────────────
echo ""
echo "=== Strandedness (STAR ReadsPerGene totals, gene rows only) ==="
printf "%-24s %14s %14s %14s %10s\n" "sample" "col2_unstr" "col3_fwd" "col4_rev" "col3_share"
printf "%s" "$STRAND_TABLE"
echo ""
echo "col3_share = col3 / (col3 + col4). Reading guide, not a rule: close to 0.5"
echo "suggests unstranded, close to 1 forward, close to 0 reverse. Set"
echo "STRANDEDNESS in config.sh yourself; the pipeline does not infer it."

# ── Aggregated outputs ────────────────────────────────────────────────────────
AGG_FAIL=false
if [[ "$PHASE2" == true ]]; then
    echo ""
    echo "=== Aggregated output ==="
    for FILE in "$OUT_STEP/counts/gene_counts.txt" \
                "$OUT_STEP/counts/gene_counts.txt.summary" \
                "$OUT_STEP/multiqc/multiqc_report.html"; do
        if [[ ! -f "$FILE" ]]; then
            echo "MISSING    : $(basename "$FILE")"; AGG_FAIL=true
        elif [[ ! -s "$FILE" ]]; then
            echo "EMPTY      : $(basename "$FILE")"; AGG_FAIL=true
        else
            echo "OK         : $(basename "$FILE")"
        fi
    done
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "=== Summary ==="
echo "  Phase checked : $([[ "$PHASE2" == true ]] && echo "1 and 2 (01a to 01d)" || echo "1 only (01a, 01b)")"
echo "  Total : $(( PASS + FAIL )) | Passed : $PASS | Failed : $FAIL"
echo "  (missing: $MISSING | empty: $EMPTY | blank controls with no reads: $BLANK_EMPTY)"
echo ""

if [[ "$FAIL" -gt 0 || "$AGG_FAIL" == true ]]; then
    echo "WARNING: some outputs are missing or empty."
    echo "         Check the logs in $LOG_DIR before proceeding."
elif [[ "$PHASE2" == false ]]; then
    echo "Phase 1 complete. Inspect the strandedness table, set STRANDEDNESS in"
    echo "config.sh, then rerun: bash 01.run_host.sh"
else
    echo "All checks passed."
    echo "Open $OUT_STEP/multiqc/multiqc_report.html for the feasibility readout."
fi
