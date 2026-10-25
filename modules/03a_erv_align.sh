#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
# ==============================================================================
# VERVE module 03a: permissive alignment for the ERV arm (array job)
#
# bowtie2 against GRCh38 with multi-mapping allowed, piped straight into a BAM.
#
# NOT coordinate-sorted, deliberately. Telescope requires all alignments of a
# read pair to be adjacent in the file, which is how bowtie2 writes them;
# sorting by position would break that and cost a lot of I/O for nothing.
#
# This is a separate alignment from module 01: STAR there answers "where does
# this read come from in the transcriptome", bowtie2 here answers "everywhere
# this read could have come from", which is what the EM model needs.
# ==============================================================================

source "${SLURM_SUBMIT_DIR:-.}/config.sh"

module load "$PKG_BOWTIE2"
module use "$MODULE_PATH_QIB"
module load "$PKG_SAMTOOLS"

echo "Job started"; hostname; date
echo ""

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_SHEET" | tr -d '\r')
SAMPLE_ID=$(echo "$LINE" | awk '{print $1}')

if [[ -z "$SAMPLE_ID" ]]; then
    echo "ERROR: could not parse line $((SLURM_ARRAY_TASK_ID + 1)) of $SAMPLE_SHEET"
    exit 1
fi

TRIM1="$OUT_DIR/00_preprocess/trimmed/${SAMPLE_ID}_R1.trimmed.fastq.gz"
TRIM2="$OUT_DIR/00_preprocess/trimmed/${SAMPLE_ID}_R2.trimmed.fastq.gz"

for F in "$TRIM1" "$TRIM2"; do
    if [[ ! -f "$F" ]]; then
        echo "ERROR: trimmed input not found for $SAMPLE_ID: $F"
        exit 1
    fi
done

ALIGN_DIR="$OUT_DIR/03_erv/align"
[[ ! -d "$ALIGN_DIR" ]] && mkdir -p "$ALIGN_DIR"
BAM="$ALIGN_DIR/${SAMPLE_ID}.erv.bam"
LOG="$ALIGN_DIR/${SAMPLE_ID}.bowtie2.log"

if [[ -s "$BAM" ]]; then
    echo "Alignment already exists for $SAMPLE_ID, skipping"
    echo "Job done"; date
    exit 0
fi

echo "Aligning sample : $SAMPLE_ID"
echo "  index         : $BOWTIE2_INDEX"
echo "  settings      : $BOWTIE2_PRESET -k $BOWTIE2_K --score-min $BOWTIE2_SCORE_MIN"
echo ""

# bowtie2 writes its summary to stderr, so it is teed to a log for MultiQC and
# for the check script.
bowtie2 \
    $BOWTIE2_PRESET \
    -k "$BOWTIE2_K" \
    --score-min "$BOWTIE2_SCORE_MIN" \
    -p "$CPUS_ERV_ALIGN" \
    -x "$BOWTIE2_INDEX" \
    -1 "$TRIM1" \
    -2 "$TRIM2" \
    $BOWTIE2_EXTRA_ARGS \
    2> "$LOG" \
    | samtools view -b -o "${BAM}.tmp" -

STATUS=("${PIPESTATUS[@]}")
cat "$LOG"

if [[ "${STATUS[0]}" -ne 0 || "${STATUS[1]}" -ne 0 ]]; then
    echo "ERROR: alignment failed for $SAMPLE_ID (bowtie2 ${STATUS[0]}, samtools ${STATUS[1]})"
    rm -f "${BAM}.tmp"
    exit 1
fi

# Rename only on success, so a partial BAM is never mistaken for a finished one
mv "${BAM}.tmp" "$BAM"

echo ""
echo "BAM: $BAM ($(du -h "$BAM" | cut -f1))"
echo "Job done"; date
