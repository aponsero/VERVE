#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
# ==============================================================================
# VERVE module 00a: preprocessing (array job, one task per sample)
#
# fastp adapter and quality trimming, then FastQC on the trimmed reads.
# Submitted by 00.run_preprocess.sh; do not call sbatch on this file directly.
#
# NOTE ON DUPLICATES
# fastp is run WITHOUT --dedup, deliberately and permanently. Module 01 needs
# the library complexity and duplication measurement, and removing duplicates
# here destroys it. Duplicates are marked after alignment (module 01), never
# removed at the read level. --dont_eval_duplication is also not passed, so the
# duplication rate estimate still appears in the fastp report.
# This is not exposed as a config switch: it is a design constraint of the
# pipeline, not a user preference.
# ==============================================================================

source "${SLURM_SUBMIT_DIR:-.}/config.sh"

module use "$MODULE_PATH_QIB"
module load "$PKG_FASTP"
if [[ "$RUN_FASTQC" == true ]]; then
    module load "$PKG_FASTQC"
fi

echo "Job started"; hostname; date
echo ""

# ── Parse sample sheet ────────────────────────────────────────────────────────
# tr -d '\r' guards against carriage returns from Excel-saved sample sheets
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_SHEET" | tr -d '\r')
SAMPLE_ID=$(echo "$LINE" | awk '{print $1}')
READ1=$(echo "$LINE"    | awk '{print $2}')
READ2=$(echo "$LINE"    | awk '{print $3}')

if [[ -z "$SAMPLE_ID" || -z "$READ1" || -z "$READ2" ]]; then
    echo "ERROR: could not parse line $((SLURM_ARRAY_TASK_ID + 1)) of $SAMPLE_SHEET"
    echo "       got: sample_id='$SAMPLE_ID' read1='$READ1' read2='$READ2'"
    exit 1
fi

for F in "$READ1" "$READ2"; do
    if [[ ! -f "$F" ]]; then
        echo "ERROR: input file not found for $SAMPLE_ID: $F"
        exit 1
    fi
done

echo "Processing sample : $SAMPLE_ID"
echo "  read1 : $READ1"
echo "  read2 : $READ2"
echo ""

# ── Output directories ────────────────────────────────────────────────────────
OUT_STEP="$OUT_DIR/00_preprocess"
TRIM_DIR="$OUT_STEP/trimmed"
FASTP_DIR="$OUT_STEP/fastp"
FASTQC_DIR="$OUT_STEP/fastqc"

for D in "$TRIM_DIR" "$FASTP_DIR" "$FASTQC_DIR"; do
    [[ ! -d "$D" ]] && mkdir -p "$D"
done

TRIM1="$TRIM_DIR/${SAMPLE_ID}_R1.trimmed.fastq.gz"
TRIM2="$TRIM_DIR/${SAMPLE_ID}_R2.trimmed.fastq.gz"

# ── Adapter handling ──────────────────────────────────────────────────────────
# Priority: explicit adapter sequences from config > auto-detection for PE.
# fastp always performs overlap-based adapter trimming for paired-end reads in
# addition to whichever of these applies.
ADAPTER_ARGS=""
if [[ -n "$FASTP_ADAPTER_R1" ]]; then
    ADAPTER_ARGS="--adapter_sequence $FASTP_ADAPTER_R1"
    if [[ -n "$FASTP_ADAPTER_R2" ]]; then
        ADAPTER_ARGS="$ADAPTER_ARGS --adapter_sequence_r2 $FASTP_ADAPTER_R2"
    fi
    echo "Adapters : explicit sequences from config.sh"
elif [[ "$FASTP_DETECT_ADAPTER_PE" == true ]]; then
    ADAPTER_ARGS="--detect_adapter_for_pe"
    echo "Adapters : auto-detection enabled (--detect_adapter_for_pe)"
else
    echo "Adapters : overlap-based trimming only"
fi
echo ""

# ── fastp ─────────────────────────────────────────────────────────────────────
echo "Running fastp"
fastp \
    --in1 "$READ1" \
    --in2 "$READ2" \
    --out1 "$TRIM1" \
    --out2 "$TRIM2" \
    $ADAPTER_ARGS \
    --qualified_quality_phred "$FASTP_QUALIFIED_QUALITY" \
    --unqualified_percent_limit "$FASTP_UNQUALIFIED_PERCENT" \
    --n_base_limit "$FASTP_N_BASE_LIMIT" \
    --length_required "$FASTP_LENGTH_REQUIRED" \
    --thread "$CPUS_PREPROCESS" \
    --json "$FASTP_DIR/${SAMPLE_ID}.fastp.json" \
    --html "$FASTP_DIR/${SAMPLE_ID}.fastp.html" \
    --report_title "$SAMPLE_ID" \
    $FASTP_EXTRA_ARGS

FASTP_STATUS=$?
if [[ "$FASTP_STATUS" -ne 0 ]]; then
    echo "ERROR: fastp exited with status $FASTP_STATUS for $SAMPLE_ID"
    exit "$FASTP_STATUS"
fi
echo "fastp complete"
echo ""

# ── FastQC on trimmed reads ───────────────────────────────────────────────────
if [[ "$RUN_FASTQC" == true ]]; then
    echo "Running FastQC"
    fastqc \
        --threads "$CPUS_PREPROCESS" \
        --outdir "$FASTQC_DIR" \
        "$TRIM1" "$TRIM2"

    FASTQC_STATUS=$?
    if [[ "$FASTQC_STATUS" -ne 0 ]]; then
        echo "ERROR: FastQC exited with status $FASTQC_STATUS for $SAMPLE_ID"
        exit "$FASTQC_STATUS"
    fi
    echo "FastQC complete"
else
    echo "FastQC skipped (RUN_FASTQC=false)"
fi

echo ""
echo "Job done"; date
