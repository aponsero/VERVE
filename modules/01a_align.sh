#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
# ==============================================================================
# VERVE module 01a: host alignment (array job, one task per sample)
#
# STAR against GRCh38. Produces:
#   - coordinate-sorted BAM + index (input to 01c, 01d, module 02)
#   - ReadsPerGene table (--quantMode GeneCounts), used to decide strandedness
#   - unmapped read pairs as gzipped FASTQ (input to module 04, so STAR never
#     has to be rerun for the viral arm)
#
# Alignment parameters are STAR defaults. Only output options are set here.
# A read group (ID, SM, LB = sample_id) is written to the header and to every
# read: Picard MarkDuplicates requires it, and uses LB to decide which reads can
# be duplicates of each other (one library per sample). It does not change the
# alignment.
#
# Samples whose outputs already exist AND carry a read group are skipped, so the
# step can be resubmitted safely. BAMs from before the read group was added are
# realigned automatically.
# ==============================================================================

source "${SLURM_SUBMIT_DIR:-.}/config.sh"

module use "$MODULE_PATH_QIB"
module load "$PKG_STAR"
module load "$PKG_SAMTOOLS"

echo "Job started"; hostname; date
echo ""

# ── Parse sample sheet ────────────────────────────────────────────────────────
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_SHEET" | tr -d '\r')
SAMPLE_ID=$(echo "$LINE" | awk '{print $1}')

if [[ -z "$SAMPLE_ID" ]]; then
    echo "ERROR: could not parse line $((SLURM_ARRAY_TASK_ID + 1)) of $SAMPLE_SHEET"
    exit 1
fi

# ── Inputs from module 00 ─────────────────────────────────────────────────────
TRIM1="$OUT_DIR/00_preprocess/trimmed/${SAMPLE_ID}_R1.trimmed.fastq.gz"
TRIM2="$OUT_DIR/00_preprocess/trimmed/${SAMPLE_ID}_R2.trimmed.fastq.gz"

for F in "$TRIM1" "$TRIM2"; do
    if [[ ! -f "$F" ]]; then
        echo "ERROR: trimmed input not found for $SAMPLE_ID: $F"
        echo "       has module 00 completed for this sample?"
        exit 1
    fi
done

# ── Outputs ───────────────────────────────────────────────────────────────────
ALIGN_DIR="$OUT_DIR/01_host/align"
UNMAPPED_DIR="$OUT_DIR/01_host/unmapped"
for D in "$ALIGN_DIR" "$UNMAPPED_DIR"; do
    [[ ! -d "$D" ]] && mkdir -p "$D"
done

PREFIX="$ALIGN_DIR/${SAMPLE_ID}."
BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
RPG="${PREFIX}ReadsPerGene.out.tab"
UNMAP1="$UNMAPPED_DIR/${SAMPLE_ID}_R1.unmapped.fastq.gz"
UNMAP2="$UNMAPPED_DIR/${SAMPLE_ID}_R2.unmapped.fastq.gz"

# ── Skip if already done ──────────────────────────────────────────────────────
if [[ -s "$BAM" && -s "${BAM}.bai" && -s "$RPG" && -f "$UNMAP1" && -f "$UNMAP2" ]]; then
    if samtools view -H "$BAM" | grep -q '^@RG'; then
        echo "Outputs already exist for $SAMPLE_ID, skipping alignment"
        echo "Job done"; date
        exit 0
    fi
    echo "Existing BAM for $SAMPLE_ID has no read group: realigning"
    echo ""
fi

echo "Aligning sample : $SAMPLE_ID"
echo "  STAR index    : $STAR_INDEX"
echo ""

# STAR refuses to start if its temporary directory exists (left over from a
# failed previous attempt)
rm -rf "${PREFIX}_STARtmp"

# ── STAR ──────────────────────────────────────────────────────────────────────
STAR \
    --runThreadN "$CPUS_ALIGN" \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$TRIM1" "$TRIM2" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$PREFIX" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outReadsUnmapped Fastx \
    --outSAMattrRGline "ID:${SAMPLE_ID}" "SM:${SAMPLE_ID}" "LB:${SAMPLE_ID}" \
    $STAR_EXTRA_ARGS

STAR_STATUS=$?
if [[ "$STAR_STATUS" -ne 0 ]]; then
    echo "ERROR: STAR exited with status $STAR_STATUS for $SAMPLE_ID"
    exit "$STAR_STATUS"
fi
echo "STAR complete"
echo ""

# ── Index BAM ─────────────────────────────────────────────────────────────────
samtools index -@ "$CPUS_ALIGN" "$BAM"
if [[ $? -ne 0 ]]; then
    echo "ERROR: samtools index failed for $SAMPLE_ID"
    exit 1
fi

# ── Compress and move unmapped reads for module 04 ────────────────────────────
for MATE in 1 2; do
    SRC="${PREFIX}Unmapped.out.mate${MATE}"
    DEST="$UNMAPPED_DIR/${SAMPLE_ID}_R${MATE}.unmapped.fastq.gz"
    if [[ ! -f "$SRC" ]]; then
        echo "ERROR: STAR did not write $SRC"
        exit 1
    fi
    gzip -c "$SRC" > "$DEST" && rm "$SRC"
    if [[ $? -ne 0 ]]; then
        echo "ERROR: failed to compress unmapped mate $MATE for $SAMPLE_ID"
        exit 1
    fi
done
echo "Unmapped reads written to $UNMAPPED_DIR"

echo ""
echo "Job done"; date
