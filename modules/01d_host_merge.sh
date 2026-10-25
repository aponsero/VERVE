#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
# ==============================================================================
# VERVE module 01d: gene counts and aggregate report (single job)
#
#   - featureCounts on all marked BAMs in one call, giving the gene count matrix
#     directly (samples in sample sheet order). Parameters are featureCounts
#     defaults apart from paired-end fragment counting and strandedness.
#   - MultiQC over all module 01 outputs: the feasibility readout.
#
# Submitted after 01b and 01c with --dependency=afterok on both.
# ==============================================================================

source "${SLURM_SUBMIT_DIR:-.}/config.sh"

module use "$MODULE_PATH_QIB"
module load "$PKG_SUBREAD"
module load "$PKG_MULTIQC"

echo "Job started"; hostname; date
echo ""

case "$STRANDEDNESS" in
    unstranded) FC_STRAND=0 ;;
    forward)    FC_STRAND=1 ;;
    reverse)    FC_STRAND=2 ;;
    *)
        echo "ERROR: STRANDEDNESS must be unstranded, forward or reverse (got '$STRANDEDNESS')"
        exit 1 ;;
esac

OUT_STEP="$OUT_DIR/01_host"
COUNTS_DIR="$OUT_STEP/counts"
MULTIQC_DIR="$OUT_STEP/multiqc"
for D in "$COUNTS_DIR" "$MULTIQC_DIR"; do
    [[ ! -d "$D" ]] && mkdir -p "$D"
done

# ── Collect BAMs in sample sheet order ────────────────────────────────────────
BAMS=()
N_SAMPLES=0
while IFS=$'\t' read -r SAMPLE_ID _; do
    [[ "$SAMPLE_ID" == "sample_id" || -z "$SAMPLE_ID" ]] && continue
    SAMPLE_ID=$(echo "$SAMPLE_ID" | tr -d '\r')
    ((N_SAMPLES++))
    MD_BAM="$OUT_STEP/markdup/${SAMPLE_ID}.markdup.bam"
    if [[ -s "$MD_BAM" ]]; then
        BAMS+=("$MD_BAM")
    else
        echo "WARNING: no marked BAM for $SAMPLE_ID, excluded from the count matrix"
    fi
done < "$SAMPLE_SHEET"

echo "Samples in sample sheet : $N_SAMPLES"
echo "BAMs to count           : ${#BAMS[@]}"
echo ""

if [[ "${#BAMS[@]}" -eq 0 ]]; then
    echo "ERROR: no marked BAMs found in $OUT_STEP/markdup"
    exit 1
fi

# ── featureCounts ─────────────────────────────────────────────────────────────
# subread 2.0.1: -p alone counts fragments (read pairs). From 2.0.2 onwards
# --countReadPairs must be added as well, or reads are counted instead.
echo "Running featureCounts (strandedness: $STRANDEDNESS, -s $FC_STRAND)"
featureCounts \
    -T "$CPUS_HOST_MERGE" \
    -p \
    -s "$FC_STRAND" \
    -a "$GTF" \
    -o "$COUNTS_DIR/gene_counts.txt" \
    $FEATURECOUNTS_EXTRA_ARGS \
    "${BAMS[@]}"

if [[ $? -ne 0 ]]; then
    echo "ERROR: featureCounts failed"
    exit 1
fi
echo "Gene count matrix: $COUNTS_DIR/gene_counts.txt"
echo ""

# ── MultiQC ───────────────────────────────────────────────────────────────────
# Scans STAR logs and ReadsPerGene, SortMeRNA logs, Picard metrics, samtools
# flagstat/idxstats and the featureCounts summary.
echo "Running MultiQC"
multiqc "$OUT_STEP" \
    --outdir "$MULTIQC_DIR" \
    --filename "multiqc_report.html" \
    --title "$MULTIQC_TITLE_HOST" \
    --force \
    $MULTIQC_EXTRA_ARGS

if [[ $? -ne 0 ]]; then
    echo "ERROR: MultiQC failed"
    exit 1
fi

echo ""
echo "Report written to: $MULTIQC_DIR/multiqc_report.html"
echo "Job done"; date
