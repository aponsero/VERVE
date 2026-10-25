#!/bin/bash -l
# ==============================================================================
# VERVE step 04: check viral detection outputs
#
# Usage: bash 04.check_viral.sh
# Run from the pipeline directory once the submitted jobs have finished.
# ==============================================================================

source config.sh

OUT_STEP="$OUT_DIR/04_viral"

echo "=== Viral detection check ==="
echo ""
echo "Configuration"
echo "  Sample sheet   : $SAMPLE_SHEET"
echo "  Output         : $OUT_STEP"
echo "  Index          : $VIRAL_INDEX"
echo "  Taxon table    : $VIRAL_TAXA"
echo "  bowtie2        : $PKG_BOWTIE2 ($BOWTIE2_VIRAL_PRESET, best hit)"
echo "  High-conf MAPQ : $VIRAL_MIN_MAPQ"
echo "  Keep BAMs      : $VIRAL_KEEP_BAM"
echo ""

PASS=0; FAIL=0

while IFS=$'\t' read -r SAMPLE_ID READ1 READ2 BODY_SITE GROUP BATCH; do
    [[ "$SAMPLE_ID" == "sample_id" || -z "$SAMPLE_ID" ]] && continue
    SAMPLE_ID=$(echo "$SAMPLE_ID" | tr -d '\r')
    SAMPLE_FAIL=false

    LOG="$OUT_STEP/align/${SAMPLE_ID}.bowtie2.log"
    METRICS="$OUT_STEP/counts_per_sample/${SAMPLE_ID}.viral_metrics.tsv"
    BAM="$OUT_STEP/align/${SAMPLE_ID}.viral.bam"

    if [[ ! -s "$LOG" ]]; then
        echo "MISSING    : $SAMPLE_ID -- bowtie2 log"; SAMPLE_FAIL=true
    fi
    if [[ ! -s "$METRICS" ]]; then
        echo "MISSING    : $SAMPLE_ID -- metrics table"; SAMPLE_FAIL=true
    fi
    if [[ ! -s "$BAM" && "$VIRAL_KEEP_BAM" == true ]]; then
        echo "MISSING    : $SAMPLE_ID -- viral BAM"; SAMPLE_FAIL=true
    fi

    if [[ "$SAMPLE_FAIL" == false ]]; then
        # Metrics hold one row per viral sequence with at least one read; an
        # empty table is a legitimate result, not a failure
        N_REFS=$(( $(wc -l < "$METRICS") - 1 ))
        READS=$(awk -F'\t' 'NR > 1 { total += $4 } END { print total + 0 }' "$METRICS")
        RATE=$(grep "overall alignment rate" "$LOG" | awk '{print $1}')
        printf "OK         : %-14s viral refs hit %-6s viral reads %-9s align rate %s\n" \
               "$SAMPLE_ID" "$N_REFS" "$READS" "${RATE:-NA}"
        ((PASS++))
    else
        ((FAIL++))
    fi

done < "$SAMPLE_SHEET"

echo ""
echo "=== Aggregated output ==="
AGG_FAIL=false
for FILE in "$OUT_STEP/counts/viral_per_accession.tsv" \
            "$OUT_STEP/counts/viral_per_taxid.tsv" \
            "$OUT_STEP/counts/viral_denominators.tsv" \
            "$OUT_STEP/counts/viral_counts_taxid.tsv" \
            "$OUT_STEP/counts/viral_counts_family.tsv" \
            "$OUT_STEP/multiqc/multiqc_report.html"; do
    if [[ ! -s "$FILE" ]]; then
        echo "MISSING    : $(basename "$FILE")"; AGG_FAIL=true
    else
        echo "OK         : $(basename "$FILE")"
    fi
done

if [[ -s "$OUT_STEP/counts/viral_per_taxid.tsv" ]]; then
    echo ""
    echo "Taxa with reads in at least one sample: $(( $(wc -l < "$OUT_STEP/counts/viral_counts_taxid.tsv") - 1 ))"
    echo ""
    echo "Flagged taxa detected (contamination-prone; not findings without scrutiny):"
    awk -F'\t' 'NR > 1 && $7 != "" { printf "  %-14s %-40s reads %s\n", $1, $3, $11 }' \
        "$OUT_STEP/counts/viral_per_taxid.tsv" | head -20
    if ! awk -F'\t' 'NR > 1 && $7 != "" { found = 1 } END { exit !found }' \
        "$OUT_STEP/counts/viral_per_taxid.tsv"; then
        echo "  none"
    fi
fi

echo ""
echo "=== Summary ==="
echo "  Total : $(( PASS + FAIL )) | Passed : $PASS | Failed : $FAIL"
echo ""

if [[ "$FAIL" -gt 0 || "$AGG_FAIL" == true ]]; then
    echo "WARNING: some outputs are missing. Check the logs in $LOG_DIR."
else
    echo "All checks passed."
    echo "No positivity call has been made: counts, breadth and covered windows"
    echo "are written out, and the rule is set in the module 04 report."
fi

if [[ -d "$OUT_STEP/align" ]]; then
    echo ""
    echo "Disk used by viral alignments: $(du -sh "$OUT_STEP/align" 2>/dev/null | cut -f1)"
fi

