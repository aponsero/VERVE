#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
# ==============================================================================
# VERVE module 01b: rRNA measurement (array job, one task per sample)
#
# SortMeRNA against two databases kept as separate files, so the log reports
# them separately:
#   - prokaryotic 16S + 23S (bacteria and archaea)
#   - eukaryotic 18S + 28S, a PROXY for human rRNA (fungal rRNA also lands
#     here, so it may be inflated at sites with fungal load, vaginal in
#     particular)
#
# MEASUREMENT ONLY. No reads are removed and no read files are written
# (--fastx is deliberately not passed). STAR in 01a receives the full module 00
# output regardless of what this step finds.
#
# Runs on all reads, no subsampling.
# ==============================================================================

source "${SLURM_SUBMIT_DIR:-.}/config.sh"

module use "$MODULE_PATH_QIB"
module load "$PKG_SORTMERNA"

echo "Job started"; hostname; date
echo ""

# ── Parse sample sheet ────────────────────────────────────────────────────────
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

# ── Outputs ───────────────────────────────────────────────────────────────────
RRNA_DIR="$OUT_DIR/01_host/rrna"
[[ ! -d "$RRNA_DIR" ]] && mkdir -p "$RRNA_DIR"

# --aligned sets the output prefix: the log is written as ${SAMPLE_ID}.log
ALIGNED_PREFIX="$RRNA_DIR/${SAMPLE_ID}"
LOG_FILE="${ALIGNED_PREFIX}.log"

if [[ -s "$LOG_FILE" ]]; then
    echo "SortMeRNA log already exists for $SAMPLE_ID, skipping"
    echo "Job done"; date
    exit 0
fi

# Per-sample working directory. SortMeRNA refuses to run if its key-value
# store (workdir/kvdb) is not empty, so it is cleared before and after.
WORKDIR="$RRNA_DIR/work/${SAMPLE_ID}"
rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"

echo "Measuring rRNA : $SAMPLE_ID"
echo "  prokaryote   : $SORTMERNA_DB_PROKARYOTE"
echo "  eukaryote    : $SORTMERNA_DB_EUKARYOTE"
echo "  index        : $SORTMERNA_IDX_DIR"
echo ""

# --index 0: use the prebuilt index only, never build it inside an array task
sortmerna \
    --ref "$SORTMERNA_DB_PROKARYOTE" \
    --ref "$SORTMERNA_DB_EUKARYOTE" \
    --reads "$TRIM1" \
    --reads "$TRIM2" \
    --idx-dir "$SORTMERNA_IDX_DIR" \
    --index 0 \
    --workdir "$WORKDIR" \
    --aligned "$ALIGNED_PREFIX" \
    --threads "$CPUS_RRNA"

SMR_STATUS=$?
rm -rf "$WORKDIR"

if [[ "$SMR_STATUS" -ne 0 ]]; then
    echo "ERROR: SortMeRNA exited with status $SMR_STATUS for $SAMPLE_ID"
    # Remove any partial log so that a resubmission does not skip this sample
    rm -f "$LOG_FILE"
    exit "$SMR_STATUS"
fi

echo "SortMeRNA complete: $LOG_FILE"
echo ""
echo "Job done"; date
