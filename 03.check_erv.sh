#!/bin/bash -l
# ==============================================================================
# VERVE step 03: check ERV quantification outputs
#
# Usage: bash 03.check_erv.sh
# Run from the pipeline directory once the submitted jobs have finished.
# ==============================================================================

source config.sh

OUT_STEP="$OUT_DIR/03_erv"

echo "=== ERV quantification check ==="
echo ""
echo "Configuration"
echo "  Sample sheet  : $SAMPLE_SHEET"
echo "  Output        : $OUT_STEP"
echo "  Annotation    : $ERV_GTF"
echo "  bowtie2       : $PKG_BOWTIE2 ($BOWTIE2_PRESET -k $BOWTIE2_K --score-min $BOWTIE2_SCORE_MIN)"
echo "  Telescope     : $PKG_TELESCOPE (stranded mode $TELESCOPE_STRANDED_MODE)"
echo "  BAMs deleted  : $ERV_DELETE_BAM_AFTER_ASSIGN"
echo ""

PASS=0; FAIL=0

while IFS=$'\t' read -r SAMPLE_ID READ1 READ2 BODY_SITE GROUP BATCH; do
    [[ "$SAMPLE_ID" == "sample_id" || -z "$SAMPLE_ID" ]] && continue
    SAMPLE_ID=$(echo "$SAMPLE_ID" | tr -d '\r')
    SAMPLE_FAIL=false

    BAM="$OUT_STEP/align/${SAMPLE_ID}.erv.bam"
    LOG="$OUT_STEP/align/${SAMPLE_ID}.bowtie2.log"
    REPORT=$(ls "$OUT_STEP/telescope/${SAMPLE_ID}/${SAMPLE_ID}-"*counts.tsv \
                "$OUT_STEP/telescope/${SAMPLE_ID}/${SAMPLE_ID}-"*report.tsv 2>/dev/null | head -1)

    # The BAM is optional once Telescope has run, if deletion is switched on
    if [[ ! -s "$BAM" && "$ERV_DELETE_BAM_AFTER_ASSIGN" != true ]]; then
        echo "MISSING    : $SAMPLE_ID -- alignment BAM"
        SAMPLE_FAIL=true
    fi
    if [[ ! -s "$LOG" ]]; then
        echo "MISSING    : $SAMPLE_ID -- bowtie2 log"
        SAMPLE_FAIL=true
    fi
    if [[ -z "$REPORT" ]]; then
        echo "MISSING    : $SAMPLE_ID -- Telescope report"
        SAMPLE_FAIL=true
    fi

    if [[ "$SAMPLE_FAIL" == false ]]; then
        RATE=$(grep "overall alignment rate" "$LOG" | awk '{print $1}')
        # Skip the comment line (absent in some Telescope versions) and the
        # header, then count loci with a non-zero final_count (column 3)
        LOCI=$(awk '/^#/ { next }
                    !seen_header { seen_header = 1; next }
                    $1 != "__no_feature" && $3 > 0 { n++ }
                    END { print n+0 }' "$REPORT")
        SIZE=$([[ -s "$BAM" ]] && du -h "$BAM" | cut -f1 || echo "-")
        printf "OK         : %-14s align %-7s loci detected %-7s BAM %s\n" \
               "$SAMPLE_ID" "${RATE:-NA}" "$LOCI" "$SIZE"
        ((PASS++))
    else
        ((FAIL++))
    fi

done < "$SAMPLE_SHEET"

echo ""
echo "=== Aggregated output ==="
AGG_FAIL=false
for FILE in "$OUT_STEP/counts/erv_counts_locus.tsv" \
            "$OUT_STEP/counts/erv_counts_family.tsv" \
            "$OUT_STEP/counts/erv_counts_stats.tsv" \
            "$OUT_STEP/multiqc/multiqc_report.html"; do
    if [[ ! -s "$FILE" ]]; then
        echo "MISSING    : $(basename "$FILE")"; AGG_FAIL=true
    else
        echo "OK         : $(basename "$FILE")"
    fi
done

if [[ -s "$OUT_STEP/counts/erv_counts_locus.tsv" ]]; then
    echo ""
    echo "Loci in matrix     : $(( $(wc -l < "$OUT_STEP/counts/erv_counts_locus.tsv") - 1 ))"
    echo "Families in matrix : $(( $(wc -l < "$OUT_STEP/counts/erv_counts_family.tsv") - 1 ))"
fi

echo ""
echo "=== Summary ==="
echo "  Total : $(( PASS + FAIL )) | Passed : $PASS | Failed : $FAIL"
echo ""

if [[ "$FAIL" -gt 0 || "$AGG_FAIL" == true ]]; then
    echo "WARNING: some outputs are missing. Check the logs in $LOG_DIR."
else
    echo "All checks passed."
    echo "Counts are RAW. Normalise with the module 02 size factors in R."
fi

if [[ -d "$OUT_STEP/align" ]]; then
    echo ""
    echo "Disk used by ERV alignments: $(du -sh "$OUT_STEP/align" 2>/dev/null | cut -f1)"
fi
