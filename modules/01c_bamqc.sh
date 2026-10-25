#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
# ==============================================================================
# VERVE module 01c: alignment QC (array job, one task per sample)
#
#   - Picard MarkDuplicates: duplicates MARKED, never removed. The marked BAM is
#     the one used downstream (01d featureCounts, module 02).
#   - Picard CollectRnaSeqMetrics: coding, UTR, intronic, intergenic fractions;
#     MEDIAN_CV_COVERAGE (comparable to the IMPACC exclusion criterion);
#     MEDIAN_5PRIME_TO_3PRIME_BIAS.
#   - samtools flagstat and idxstats.
#
# Requires STRANDEDNESS to be set in config.sh.
# ==============================================================================

source "${SLURM_SUBMIT_DIR:-.}/config.sh"

module use "$MODULE_PATH_QIB"
module load "$PKG_PICARD"
module load "$PKG_SAMTOOLS"

echo "Job started"; hostname; date
echo ""

# ── Strandedness mapping ──────────────────────────────────────────────────────
case "$STRANDEDNESS" in
    unstranded) PICARD_STRAND="NONE" ;;
    forward)    PICARD_STRAND="FIRST_READ_TRANSCRIPTION_STRAND" ;;
    reverse)    PICARD_STRAND="SECOND_READ_TRANSCRIPTION_STRAND" ;;
    *)
        echo "ERROR: STRANDEDNESS must be unstranded, forward or reverse (got '$STRANDEDNESS')"
        exit 1 ;;
esac

# ── Parse sample sheet ────────────────────────────────────────────────────────
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_SHEET" | tr -d '\r')
SAMPLE_ID=$(echo "$LINE" | awk '{print $1}')

if [[ -z "$SAMPLE_ID" ]]; then
    echo "ERROR: could not parse line $((SLURM_ARRAY_TASK_ID + 1)) of $SAMPLE_SHEET"
    exit 1
fi

BAM="$OUT_DIR/01_host/align/${SAMPLE_ID}.Aligned.sortedByCoord.out.bam"
if [[ ! -s "$BAM" ]]; then
    echo "ERROR: STAR BAM not found or empty for $SAMPLE_ID: $BAM"
    exit 1
fi

# MarkDuplicates fails with a NullPointerException on reads without a read
# group. Catch it here with a clear message instead.
if ! samtools view -H "$BAM" | grep -q '^@RG'; then
    echo "ERROR: $BAM has no read group (@RG). Rerun 01a for this sample:"
    echo "       01a now adds one and realigns BAMs that lack it."
    exit 1
fi

# ── Outputs ───────────────────────────────────────────────────────────────────
MD_DIR="$OUT_DIR/01_host/markdup"
QC_DIR="$OUT_DIR/01_host/bamqc"
TMP_DIR="$OUT_DIR/01_host/tmp/${SAMPLE_ID}"
for D in "$MD_DIR" "$QC_DIR" "$TMP_DIR"; do
    [[ ! -d "$D" ]] && mkdir -p "$D"
done

MD_BAM="$MD_DIR/${SAMPLE_ID}.markdup.bam"
MD_METRICS="$QC_DIR/${SAMPLE_ID}.markdup_metrics.txt"
RNA_METRICS="$QC_DIR/${SAMPLE_ID}.rnaseq_metrics.txt"

# Java heap for Picard, applies to every java call in this job
export _JAVA_OPTIONS="-Xmx${PICARD_XMX}"

echo "QC for sample   : $SAMPLE_ID"
echo "  strandedness  : $STRANDEDNESS ($PICARD_STRAND)"
echo ""

# ── MarkDuplicates ────────────────────────────────────────────────────────────
# REMOVE_DUPLICATES is written out explicitly (it is also Picard's default) so
# that the behaviour is visible in the script. Do not change it: removing
# duplicates destroys the complexity measurement this module exists to make.
echo "Running Picard MarkDuplicates"
picard MarkDuplicates \
    --INPUT "$BAM" \
    --OUTPUT "$MD_BAM" \
    --METRICS_FILE "$MD_METRICS" \
    --REMOVE_DUPLICATES false \
    --CREATE_INDEX true \
    --TMP_DIR "$TMP_DIR"

if [[ $? -ne 0 ]]; then
    echo "ERROR: MarkDuplicates failed for $SAMPLE_ID"
    exit 1
fi
echo ""

# ── CollectRnaSeqMetrics ──────────────────────────────────────────────────────
RRNA_ARG=""
if [[ -n "$RRNA_INTERVAL_LIST" ]]; then
    RRNA_ARG="--RIBOSOMAL_INTERVALS $RRNA_INTERVAL_LIST"
fi

echo "Running Picard CollectRnaSeqMetrics"
picard CollectRnaSeqMetrics \
    --INPUT "$MD_BAM" \
    --OUTPUT "$RNA_METRICS" \
    --REF_FLAT "$REFFLAT" \
    --STRAND_SPECIFICITY "$PICARD_STRAND" \
    $RRNA_ARG \
    --TMP_DIR "$TMP_DIR"

if [[ $? -ne 0 ]]; then
    echo "ERROR: CollectRnaSeqMetrics failed for $SAMPLE_ID"
    exit 1
fi
echo ""

# ── samtools ──────────────────────────────────────────────────────────────────
echo "Running samtools flagstat and idxstats"
samtools flagstat -@ "$CPUS_BAMQC" "$MD_BAM" > "$QC_DIR/${SAMPLE_ID}.flagstat" \
    || { echo "ERROR: flagstat failed for $SAMPLE_ID"; exit 1; }
samtools idxstats "$MD_BAM" > "$QC_DIR/${SAMPLE_ID}.idxstats" \
    || { echo "ERROR: idxstats failed for $SAMPLE_ID"; exit 1; }

rm -rf "$TMP_DIR"

echo ""
echo "Job done"; date
