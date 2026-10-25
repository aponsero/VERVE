#!/bin/bash -l
# ==============================================================================
# Diagnose one viral reference in one sample.
#
# Usage, from the pipeline directory:
#   bash scripts/viral_taxon_diagnostic.sh <SAMPLE_ID> <ACCESSION or name>
#   bash scripts/viral_taxon_diagnostic.sh CMWZ_Nas NC_001617
#   bash scripts/viral_taxon_diagnostic.sh CMEH_Sal "Torque teno midi virus 9"
#
# The second argument is matched against the accession first, then against the
# virus name in the taxon table. A name matching several accessions (a
# segmented virus, or several strains) lists them and stops, so the choice is
# explicit.
#
# This is what separates the two readings of a low-breadth hit:
#   reads concentrated in a few separated regions  -> transcription of a large
#       DNA virus, or a conserved region of a divergent strain
#   one narrow spike with very high depth          -> cross-mapping
# ==============================================================================

source config.sh

SAMPLE_ID="$1"
QUERY="$2"

if [[ -z "$SAMPLE_ID" || -z "$QUERY" ]]; then
    echo "Usage: bash scripts/viral_taxon_diagnostic.sh <SAMPLE_ID> <ACCESSION or name>"
    exit 1
fi

module use "$MODULE_PATH_QIB"
module load "$PKG_SAMTOOLS"

BAM="$OUT_DIR/04_viral/align/${SAMPLE_ID}.viral.bam"
if [[ ! -s "$BAM" ]]; then
    echo "ERROR: $BAM not found."
    echo "       VIRAL_KEEP_BAM must have been true when module 04 ran."
    exit 1
fi

# ── Resolve the accession ─────────────────────────────────────────────────────
MATCHES=$(awk -F'\t' -v query="$QUERY" '
    NR == 1 { next }
    $1 == query { print; exact = 1 }
    END { if (!exact) { exit 1 } }' "$VIRAL_TAXA")

if [[ -z "$MATCHES" ]]; then
    MATCHES=$(awk -F'\t' -v query="$QUERY" '
        NR == 1 { next }
        index(tolower($3), tolower(query)) { print }' "$VIRAL_TAXA")
fi

N_MATCHES=$(echo "$MATCHES" | grep -c . )

if [[ "$N_MATCHES" -eq 0 ]]; then
    echo "ERROR: '$QUERY' matched no accession or virus name in $VIRAL_TAXA"
    exit 1
fi

if [[ "$N_MATCHES" -gt 1 ]]; then
    echo "'$QUERY' matches $N_MATCHES sequences. Pick one accession:"
    echo "$MATCHES" | awk -F'\t' '{ printf "  %-14s %-50s segment %s of %s\n", $1, $3, ($4 == "" ? "-" : $4), $5 }'
    exit 1
fi

ACCESSION=$(echo "$MATCHES" | cut -f1)
VIRUS_NAME=$(echo "$MATCHES" | cut -f3)
GENUS=$(echo "$MATCHES" | cut -f6)
FAMILY=$(echo "$MATCHES" | cut -f7)
FLAGS=$(echo "$MATCHES" | cut -f9)
REF_LEN=$(echo "$MATCHES" | cut -f10)

echo "=== $VIRUS_NAME in $SAMPLE_ID ==="
echo ""
echo "Accession  : $ACCESSION"
echo "Genus      : ${GENUS:-(none)}"
echo "Family     : ${FAMILY:-(none)}"
echo "Length     : $REF_LEN bp"
[[ -n "$FLAGS" ]] && echo "Flags      : $FLAGS"
echo ""

# ── Summary from samtools ─────────────────────────────────────────────────────
echo "=== samtools coverage ==="
samtools coverage -r "$ACCESSION" "$BAM"
echo ""

N_READS=$(samtools view -c "$BAM" "$ACCESSION")
if [[ "$N_READS" -eq 0 ]]; then
    echo "No reads on $ACCESSION in this sample. Nothing further to show."
    exit 0
fi

# ── Coverage profile ──────────────────────────────────────────────────────────
echo "=== Coverage across the reference (40 bins) ==="
samtools depth -a -r "$ACCESSION" "$BAM" | awk -v len="$REF_LEN" '
    {
        bin = int(($2 - 1) * 40 / len)
        if (bin > 39) { bin = 39 }
        total[bin] += $3
        n[bin]++
        if ($3 > peak) { peak = $3; peak_pos = $2 }
        if ($3 > 0) { covered++ }
        sum += $3
        positions++
    }
    END {
        if (positions == 0) { print "  no coverage"; exit }
        mean_overall = sum / positions
        for (i = 0; i < 40; i++) {
            mean = (n[i] > 0) ? total[i] / n[i] : 0
            bars = int(mean * 50 / (peak > 0 ? peak : 1))
            printf "  bin %2d  mean depth %10.1f  ", i + 1, mean
            for (j = 0; j < bars; j++) { printf "#" }
            printf "\n"
        }
        printf "\n  covered bases         : %d of %d (%.2f%%)\n", covered, positions, 100 * covered / positions
        printf "  mean depth over ref   : %.1f\n", mean_overall
        printf "  peak depth            : %d at position %d\n", peak, peak_pos
        printf "  peak / mean           : %.1f\n", (mean_overall > 0 ? peak / mean_overall : 0)
    }'
echo ""

# ── Covered blocks ────────────────────────────────────────────────────────────
# How many separated stretches the coverage falls into. One block is a spike;
# several well-separated blocks on a large genome is what transcription looks
# like.
echo "=== Covered blocks (gaps of more than 100 bp start a new block) ==="
samtools depth "$BAM" -r "$ACCESSION" | awk '
    $3 > 0 {
        if (prev == 0 || $2 - prev > 100) {
            if (prev != 0) { printf "  %8d - %-8d  %6d bp  mean depth %8.1f\n", start, prev, prev - start + 1, block_sum / (prev - start + 1) }
            start = $2
            block_sum = 0
            blocks++
        }
        prev = $2
        block_sum += $3
    }
    END {
        if (prev != 0) { printf "  %8d - %-8d  %6d bp  mean depth %8.1f\n", start, prev, prev - start + 1, block_sum / (prev - start + 1) }
        printf "\n  blocks: %d\n", blocks + 0
    }' | tail -25
echo ""

# ── MAPQ ──────────────────────────────────────────────────────────────────────
echo "=== MAPQ distribution (primary alignments) ==="
samtools view -F 0x100 "$BAM" "$ACCESSION" | awk '
    { mapq[$5]++; total++ }
    END {
        if (total == 0) { print "  no primary alignments"; exit }
        for (q in mapq) { printf "  MAPQ %-4s %8d  %5.1f%%\n", q, mapq[q], 100 * mapq[q] / total }
        printf "  total %12d\n", total
    }' | sort -k2 -n
echo ""
echo "bowtie2 gives a read it could not place confidently a MAPQ of 0 or 1, so"
echo "a wall of those means the assignment to this reference is not supported."
echo ""
echo "Done."
