#!/bin/bash -l
# ==============================================================================
# Diagnose one ERV locus in one sample.
#
# Usage, from the pipeline directory:
#   bash scripts/erv_locus_diagnostic.sh <SAMPLE_ID> <LOCUS>
#   bash scripts/erv_locus_diagnostic.sh CMEH_Vag HERVH_20p11.23b
#
# Coordinates are read from the annotation, so nothing needs filling in by hand.
# The sorted BAM is kept in 03_erv/diagnostics so repeat runs skip the sort.
#
# Reports:
#   1. the locus as annotated (span, strand, repName, geneRegion)
#   2. coverage across the locus in 20 bins: even coverage suggests a
#      transcript, a single spike suggests reads piling on a short motif
#   3. MAPQ distribution: bowtie2 gives multi-mapping reads MAPQ 0 or 1, so a
#      wall of MAPQ 0 means the EM model had nothing unambiguous to work with
#   4. whether the locus sits near annotated rRNA
# ==============================================================================

source config.sh

SAMPLE_ID="$1"
LOCUS="$2"

if [[ -z "$SAMPLE_ID" || -z "$LOCUS" ]]; then
    echo "Usage: bash scripts/erv_locus_diagnostic.sh <SAMPLE_ID> <LOCUS>"
    exit 1
fi

module use "$MODULE_PATH_QIB"
module load "$PKG_SAMTOOLS"

BAM="$OUT_DIR/03_erv/align/${SAMPLE_ID}.erv.bam"
DIAG_DIR="$OUT_DIR/03_erv/diagnostics"
SORTED="$DIAG_DIR/${SAMPLE_ID}.sorted.bam"
mkdir -p "$DIAG_DIR"

if [[ ! -s "$BAM" ]]; then
    echo "ERROR: $BAM not found"
    exit 1
fi

# ── 1. The locus as annotated ─────────────────────────────────────────────────
echo "=== $LOCUS in $SAMPLE_ID ==="
echo ""

COORDS=$(awk -v locus="$LOCUS" -F'\t' '
    $9 ~ "locus \"" locus "\";" {
        chrom = $1
        if (min == 0 || $4 < min) { min = $4 }
        if ($5 > max) { max = $5 }
        strand = $7
        n++
    }
    END { if (n > 0) print chrom, min, max, strand, n }' "$ERV_GTF")

if [[ -z "$COORDS" ]]; then
    echo "ERROR: $LOCUS not found in $ERV_GTF"
    echo "       check the exact spelling against the counts matrix"
    exit 1
fi

read -r CHROM START END STRAND N_EXONS <<< "$COORDS"
SPAN=$(( END - START + 1 ))

echo "Chromosome : $CHROM"
echo "Span       : ${START}-${END} (${SPAN} bp, strand $STRAND, $N_EXONS exon features)"
echo ""
echo "Annotation attributes:"
awk -v locus="$LOCUS" -F'\t' '$9 ~ "locus \"" locus "\";" { print "  " $9 }' "$ERV_GTF" | head -5
echo ""

# ── 2. Coverage ───────────────────────────────────────────────────────────────
if [[ ! -s "${SORTED}.bai" ]]; then
    echo "Sorting and indexing (once per sample, kept in $DIAG_DIR)"
    samtools sort -@ 8 -o "$SORTED" "$BAM" || exit 1
    samtools index "$SORTED" || exit 1
    echo ""
fi

echo "=== Coverage across the locus (20 bins) ==="
samtools depth -a -r "${CHROM}:${START}-${END}" "$SORTED" | \
awk -v start="$START" -v span="$SPAN" '
    {
        bin = int(($2 - start) * 20 / span)
        if (bin > 19) { bin = 19 }
        total[bin] += $3
        n[bin]++
        if ($3 > peak) { peak = $3; peak_pos = $2 }
        sum += $3
        positions++
    }
    END {
        if (positions == 0) { print "  no coverage"; exit }
        mean_overall = sum / positions
        for (i = 0; i < 20; i++) {
            mean = (n[i] > 0) ? total[i] / n[i] : 0
            bars = int(mean * 50 / (peak > 0 ? peak : 1))
            printf "  bin %2d  mean depth %10.1f  ", i + 1, mean
            for (j = 0; j < bars; j++) { printf "#" }
            printf "\n"
        }
        printf "\n  mean depth over locus : %.1f\n", mean_overall
        printf "  peak depth            : %d at position %d\n", peak, peak_pos
        printf "  peak / mean           : %.1f\n", (mean_overall > 0 ? peak / mean_overall : 0)
    }'
echo ""

# ── 3. MAPQ ───────────────────────────────────────────────────────────────────
echo "=== MAPQ distribution of reads overlapping the locus ==="
samtools view "$SORTED" "${CHROM}:${START}-${END}" | awk '
    { mapq[$5]++; total++ }
    END {
        if (total == 0) { print "  no reads"; exit }
        for (q in mapq) { printf "  MAPQ %-4s %8d  %5.1f%%\n", q, mapq[q], 100 * mapq[q] / total }
        printf "  total %12d\n", total
    }' | sort -k2 -n
echo ""

# ── 4. rRNA proximity ─────────────────────────────────────────────────────────
if [[ -n "$RRNA_INTERVAL_LIST" && -s "$RRNA_INTERVAL_LIST" ]]; then
    echo "=== Annotated rRNA genes within 10 kb ==="
    awk -v chrom="$CHROM" -v start="$START" -v end="$END" -F'\t' '
        !/^@/ && $1 == chrom && $3 > start - 10000 && $2 < end + 10000 {
            print "  " $1 ":" $2 "-" $3 "  " $5
            found++
        }
        END { if (found == 0) { print "  none" } }' "$RRNA_INTERVAL_LIST"
else
    echo "=== rRNA proximity: skipped (RRNA_INTERVAL_LIST not set) ==="
fi

echo ""
echo "Done."
