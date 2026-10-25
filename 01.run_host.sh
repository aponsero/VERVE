#!/bin/bash -l
# ==============================================================================
# VERVE step 01: submit host alignment and feasibility QC
#
# Usage: bash 01.run_host.sh
# Run from the pipeline directory (the one containing config.sh).
#
# Two phases, controlled by STRANDEDNESS in config.sh:
#   STRANDEDNESS empty -> submits 01a (STAR) and 01b (SortMeRNA) only.
#                         Then run 01.check_host.sh, decide strandedness from
#                         the ReadsPerGene totals, set it in config.sh, and
#                         rerun this script.
#   STRANDEDNESS set   -> submits 01a and 01b (both skip completed samples),
#                         then 01c after 01a, then 01d after 01b and 01c.
# ==============================================================================

source config.sh

# ── Sanity checks ─────────────────────────────────────────────────────────────
ERR=false

for VAR in SAMPLE_SHEET OUT_DIR MODULE_PATH_QIB \
           PKG_STAR PKG_SAMTOOLS PKG_SORTMERNA \
           STAR_INDEX SORTMERNA_DB_PROKARYOTE SORTMERNA_DB_EUKARYOTE SORTMERNA_IDX_DIR \
           PARTITION_ALIGN PARTITION_RRNA; do
    if [[ -z "${!VAR}" ]]; then
        echo "ERROR: $VAR is empty in config.sh"; ERR=true
    fi
done

[[ -n "$SAMPLE_SHEET" && ! -f "$SAMPLE_SHEET" ]] && { echo "ERROR: sample sheet not found: $SAMPLE_SHEET"; ERR=true; }
[[ -n "$STAR_INDEX" && ! -d "$STAR_INDEX" ]]     && { echo "ERROR: STAR index directory not found: $STAR_INDEX"; ERR=true; }
[[ -n "$SORTMERNA_IDX_DIR" && ! -d "$SORTMERNA_IDX_DIR" ]] && { echo "ERROR: SortMeRNA index directory not found: $SORTMERNA_IDX_DIR"; ERR=true; }
for F in "$SORTMERNA_DB_PROKARYOTE" "$SORTMERNA_DB_EUKARYOTE"; do
    [[ -n "$F" && ! -f "$F" ]] && { echo "ERROR: SortMeRNA database not found: $F"; ERR=true; }
done

PHASE2=false
if [[ -n "$STRANDEDNESS" ]]; then
    case "$STRANDEDNESS" in
        unstranded|forward|reverse) PHASE2=true ;;
        *) echo "ERROR: STRANDEDNESS must be empty, unstranded, forward or reverse (got '$STRANDEDNESS')"; ERR=true ;;
    esac
fi

if [[ "$PHASE2" == true ]]; then
    for VAR in PKG_PICARD PKG_SUBREAD PKG_MULTIQC GTF REFFLAT \
               PARTITION_BAMQC PARTITION_HOST_MERGE; do
        if [[ -z "${!VAR}" ]]; then
            echo "ERROR: $VAR is empty in config.sh (needed once STRANDEDNESS is set)"; ERR=true
        fi
    done
    [[ -n "$GTF" && ! -f "$GTF" ]]         && { echo "ERROR: GTF not found: $GTF"; ERR=true; }
    [[ -n "$REFFLAT" && ! -f "$REFFLAT" ]] && { echo "ERROR: refFlat not found: $REFFLAT"; ERR=true; }
    [[ -n "$RRNA_INTERVAL_LIST" && ! -f "$RRNA_INTERVAL_LIST" ]] && { echo "ERROR: rRNA interval list not found: $RRNA_INTERVAL_LIST"; ERR=true; }
fi

for S in 01a_align 01b_rrna 01c_bamqc 01d_host_merge; do
    [[ ! -f "modules/${S}.sh" ]] && { echo "ERROR: modules/${S}.sh not found. Run from the pipeline directory."; ERR=true; }
done

[[ "$ERR" == true ]] && exit 1

sed -i 's/\r//' "$SAMPLE_SHEET"
[[ ! -d "$LOG_DIR" ]] && mkdir -p "$LOG_DIR"
[[ ! -d "$OUT_DIR/01_host" ]] && mkdir -p "$OUT_DIR/01_host"

N_SAMPLES=$(tail -n +2 "$SAMPLE_SHEET" | grep -c . )
if [[ "$N_SAMPLES" -lt 1 ]]; then
    echo "ERROR: no samples found in $SAMPLE_SHEET"
    exit 1
fi

echo "=== VERVE step 01: host alignment and feasibility QC ==="
echo "  Sample sheet : $SAMPLE_SHEET"
echo "  Samples      : $N_SAMPLES"
echo "  Output       : $OUT_DIR/01_host"
echo "  Logs         : $LOG_DIR"
echo "  Strandedness : ${STRANDEDNESS:-not set (phase 1 only)}"
echo ""

# ── 01a: STAR (array) ─────────────────────────────────────────────────────────
ALIGN_ID=$(sbatch \
    --job-name=01a_align \
    --partition="$PARTITION_ALIGN" \
    --time="$TIME_ALIGN" \
    --cpus-per-task="$CPUS_ALIGN" \
    --mem-per-cpu="$MEM_PER_CPU_ALIGN" \
    --array=1-${N_SAMPLES} \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    modules/01a_align.sh | awk '{print $NF}')
echo "01a STAR array submitted      : $ALIGN_ID"

# ── 01b: SortMeRNA (array, independent of 01a) ────────────────────────────────
RRNA_ID=$(sbatch \
    --job-name=01b_rrna \
    --partition="$PARTITION_RRNA" \
    --time="$TIME_RRNA" \
    --cpus-per-task="$CPUS_RRNA" \
    --mem-per-cpu="$MEM_PER_CPU_RRNA" \
    --array=1-${N_SAMPLES} \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    modules/01b_rrna.sh | awk '{print $NF}')
echo "01b SortMeRNA array submitted : $RRNA_ID"

if [[ "$PHASE2" == false ]]; then
    echo ""
    echo "STRANDEDNESS is not set: 01c and 01d were NOT submitted."
    echo "When these jobs finish, run: bash 01.check_host.sh"
    echo "It prints the ReadsPerGene column totals. Set STRANDEDNESS in"
    echo "config.sh, then rerun: bash 01.run_host.sh"
    exit 0
fi

# ── 01c: BAM QC (array, after 01a) ────────────────────────────────────────────
BAMQC_ID=$(sbatch \
    --job-name=01c_bamqc \
    --partition="$PARTITION_BAMQC" \
    --time="$TIME_BAMQC" \
    --cpus-per-task="$CPUS_BAMQC" \
    --mem-per-cpu="$MEM_PER_CPU_BAMQC" \
    --array=1-${N_SAMPLES} \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    --dependency=afterok:${ALIGN_ID} \
    modules/01c_bamqc.sh | awk '{print $NF}')
echo "01c BAM QC array submitted    : $BAMQC_ID (after $ALIGN_ID)"

# ── 01d: featureCounts + MultiQC (after 01b and 01c) ──────────────────────────
MERGE_ID=$(sbatch \
    --job-name=01d_host_merge \
    --partition="$PARTITION_HOST_MERGE" \
    --time="$TIME_HOST_MERGE" \
    --cpus-per-task="$CPUS_HOST_MERGE" \
    --mem-per-cpu="$MEM_PER_CPU_HOST_MERGE" \
    --output="$LOG_DIR/%x_%A_%a.out" \
    --error="$LOG_DIR/%x_%A_%a.err" \
    --dependency=afterok:${RRNA_ID}:${BAMQC_ID} \
    modules/01d_host_merge.sh | awk '{print $NF}')
echo "01d merge submitted           : $MERGE_ID (after $RRNA_ID and $BAMQC_ID)"

echo ""
echo "Next step: bash 01.check_host.sh once all jobs complete"
