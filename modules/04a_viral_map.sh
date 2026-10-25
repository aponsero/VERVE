#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
# ==============================================================================
# VERVE module 04a: viral mapping and counting (array job, one task per sample)
#
# Maps the reads STAR could not place (written by 01a) against GRCh38 plus the
# viral set in a single index. GRCh38 is there as a decoy: a read that is
# really human prefers its human source over a spurious viral hit. Only the
# alignments that landed on a viral sequence are kept, so the BAM stays small.
#
# Best-hit mapping, not -k: one alignment per read, MAPQ meaningful.
#
# Per-sample metrics written for 04b:
#   reads per accession, all primary alignments and those at or above
#   VIRAL_MIN_MAPQ; covered bases; mean depth; mean MAPQ; and the number of
#   distinct covered windows, which is the part of the breadth question that
#   survives the difference between a 7 kb and a 150 kb genome.
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

UNMAP1="$OUT_DIR/01_host/unmapped/${SAMPLE_ID}_R1.unmapped.fastq.gz"
UNMAP2="$OUT_DIR/01_host/unmapped/${SAMPLE_ID}_R2.unmapped.fastq.gz"

for F in "$UNMAP1" "$UNMAP2"; do
    if [[ ! -f "$F" ]]; then
        echo "ERROR: non-host reads not found for $SAMPLE_ID: $F"
        echo "       has module 01 completed for this sample?"
        exit 1
    fi
done

OUT_STEP="$OUT_DIR/04_viral"
BAM_DIR="$OUT_STEP/align"
COUNT_DIR="$OUT_STEP/counts_per_sample"
for D in "$BAM_DIR" "$COUNT_DIR"; do
    [[ ! -d "$D" ]] && mkdir -p "$D"
done

BAM="$BAM_DIR/${SAMPLE_ID}.viral.bam"
LOG="$BAM_DIR/${SAMPLE_ID}.bowtie2.log"
METRICS="$COUNT_DIR/${SAMPLE_ID}.viral_metrics.tsv"

if [[ -s "$METRICS" ]]; then
    echo "Metrics already exist for $SAMPLE_ID, skipping"
    echo "Job done"; date
    exit 0
fi

# Viral sequence names, so human alignments can be dropped from the stream
VIRAL_NAMES="$COUNT_DIR/.${SAMPLE_ID}.viral_names"
tail -n +2 "$VIRAL_TAXA" | cut -f1 | sort -u > "$VIRAL_NAMES"
echo "Viral sequences in reference: $(wc -l < "$VIRAL_NAMES")"
echo ""

echo "Mapping sample : $SAMPLE_ID"
echo "  index        : $VIRAL_INDEX"
echo "  settings     : $BOWTIE2_VIRAL_PRESET (best hit, no -k)"
echo ""

# ── Map, keep viral alignments only, sort ─────────────────────────────────────
bowtie2 \
    $BOWTIE2_VIRAL_PRESET \
    -p "$CPUS_VIRAL" \
    -x "$VIRAL_INDEX" \
    -1 "$UNMAP1" \
    -2 "$UNMAP2" \
    $BOWTIE2_VIRAL_EXTRA_ARGS \
    2> "$LOG" \
  | samtools view -h -F 4 - \
  | awk -v names="$VIRAL_NAMES" '
        BEGIN { while ((getline line < names) > 0) { viral[line] = 1 } }
        /^@/ { print; next }
        $3 in viral { print }
    ' \
  | samtools sort -@ "$CPUS_VIRAL" -o "${BAM}.tmp" -

STATUS=("${PIPESTATUS[@]}")
cat "$LOG"

for i in "${!STATUS[@]}"; do
    if [[ "${STATUS[$i]}" -ne 0 ]]; then
        echo "ERROR: step $i of the mapping pipeline failed for $SAMPLE_ID (status ${STATUS[$i]})"
        rm -f "${BAM}.tmp"
        exit 1
    fi
done

mv "${BAM}.tmp" "$BAM"
samtools index "$BAM" || exit 1

echo ""
echo "Viral BAM: $BAM ($(du -h "$BAM" | cut -f1))"
echo ""

# ── Per-accession metrics ─────────────────────────────────────────────────────
# samtools coverage gives reads, covered bases, mean depth and mean MAPQ in one
# pass; a second pass with --min-MQ gives the high-confidence read count.
COV_ALL="$COUNT_DIR/.${SAMPLE_ID}.coverage_all"
COV_HIQ="$COUNT_DIR/.${SAMPLE_ID}.coverage_hiq"
WINDOWS="$COUNT_DIR/.${SAMPLE_ID}.windows"

samtools coverage "$BAM" > "$COV_ALL" || exit 1
samtools coverage --min-MQ "$VIRAL_MIN_MAPQ" "$BAM" > "$COV_HIQ" || exit 1

# Distinct covered windows per reference, from the covered positions only
samtools depth "$BAM" | awk -v size="$VIRAL_WINDOW_SIZE" '
    { seen[$1 "\t" int($2 / size)] = 1 }
    END { for (key in seen) { split(key, parts, "\t"); n[parts[1]]++ }
          for (ref in n) { print ref "\t" n[ref] } }
' > "$WINDOWS" || exit 1

# Join the three into one table, keeping only references with at least one read
awk -F'\t' -v OFS='\t' -v sample="$SAMPLE_ID" -v hiq="$COV_HIQ" -v win="$WINDOWS" '
    BEGIN {
        while ((getline line < hiq) > 0) {
            split(line, f, "\t")
            if (f[1] == "#rname") { continue }
            reads_hiq[f[1]] = f[4]
        }
        while ((getline line < win) > 0) {
            split(line, f, "\t")
            windows[f[1]] = f[2]
        }
        print "sample_id", "accession", "ref_length", "reads", "reads_hiq",
              "covered_bases", "breadth", "mean_depth", "mean_mapq", "windows_covered"
    }
    /^#/ { next }
    $4 > 0 {
        length_bp = $3
        breadth = (length_bp > 0) ? $5 / length_bp : 0
        print sample, $1, length_bp, $4, (($1 in reads_hiq) ? reads_hiq[$1] : 0),
              $5, sprintf("%.6f", breadth), $7, $9,
              (($1 in windows) ? windows[$1] : 0)
    }
' "$COV_ALL" > "$METRICS" || exit 1

N_REFS=$(( $(wc -l < "$METRICS") - 1 ))
echo "Viral sequences with at least one read: $N_REFS"
echo ""

rm -f "$VIRAL_NAMES" "$COV_ALL" "$COV_HIQ" "$WINDOWS"

if [[ "$VIRAL_KEEP_BAM" != true ]]; then
    echo "Deleting $BAM (VIRAL_KEEP_BAM is not true)"
    rm -f "$BAM" "${BAM}.bai"
fi

echo "Metrics: $METRICS"
echo "Job done"; date

