#!/bin/bash -l
#SBATCH --job-name=ref02_build
#SBATCH --partition=qib-compute
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=ref02_build_%j.out
#SBATCH --error=ref02_build_%j.err
# ==============================================================================
# Reference step 2: build derived references (module 01 references for VERVE)
#
# Usage: sbatch ref.02_build.sh      (from the references/ folder)
# Set --partition above before submitting. Run after ref.01_download.sh.
#
# Standalone: reads ref.config.sh only, never the pipeline's config.sh.
# At the end it prints the lines to paste into the pipeline's config.sh.
#
# Builds, in order (each step is skipped if its output already exists):
#   1. samtools faidx + Picard sequence dictionary
#   2. rRNA interval list for Picard (gene_biotype rRNA and Mt_rRNA)
#   3. refFlat for Picard (gtfToGenePred)
#   4. STAR index, built WITH the GTF (required by --quantMode GeneCounts)
#   5. SortMeRNA index for the prokaryotic and eukaryotic databases
# ==============================================================================

# Under sbatch, $0 points to a spooled copy, so the submission directory is used
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/ref.config.sh"

module use "$MODULE_PATH_QIB"
module load "$PKG_SAMTOOLS"
module load "$PKG_PICARD"
module load "$PKG_STAR"
module load "$PKG_SORTMERNA"
module load "$PKG_GTFTOGENEPRED"

echo "Job started"; hostname; date
echo ""

THREADS="${SLURM_CPUS_PER_TASK:-4}"
TODAY=$(date +%Y-%m-%d)

# ── Sanity checks ─────────────────────────────────────────────────────────────
for VAR in REF_DIR ENSEMBL_RELEASE; do
    if [[ -z "${!VAR}" ]]; then
        echo "ERROR: $VAR is empty in ref.config.sh"; exit 1
    fi
done
for F in "$GENOME_FASTA" "$GTF" "$SORTMERNA_DB_PROKARYOTE" "$SORTMERNA_DB_EUKARYOTE" "$MANIFEST"; do
    if [[ ! -s "$F" ]]; then
        echo "ERROR: $F not found. Has ref.01_download.sh completed?"; exit 1
    fi
done

manifest_add () {
    # $1 file, $2 source, $3 note. md5 only for single files.
    local MD5="NA"
    [[ -f "$1" ]] && MD5=$(md5sum "$1" | awk '{print $1}')
    printf "%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$TODAY" "$MD5" "$3" >> "$MANIFEST"
}

fail () { echo "ERROR: $1"; exit 1; }

# ── 1. FASTA index and sequence dictionary ────────────────────────────────────
echo "=== 1. FASTA index and sequence dictionary ==="
if [[ ! -s "${GENOME_FASTA}.fai" ]]; then
    samtools faidx "$GENOME_FASTA" || fail "samtools faidx"
    manifest_add "${GENOME_FASTA}.fai" "built" "$PKG_SAMTOOLS faidx"
else
    echo "exists, skipping"
fi
if [[ ! -s "$GENOME_DICT" ]]; then
    picard CreateSequenceDictionary --REFERENCE "$GENOME_FASTA" --OUTPUT "$GENOME_DICT" \
        || fail "CreateSequenceDictionary"
    manifest_add "$GENOME_DICT" "built" "$PKG_PICARD CreateSequenceDictionary"
else
    echo "exists, skipping"
fi
echo ""

# ── 2. rRNA interval list ─────────────────────────────────────────────────────
# Gene records with gene_biotype "rRNA" or "Mt_rRNA". rRNA_pseudogene is NOT
# included. GTF and interval list coordinates are both 1-based closed, so no
# conversion is needed.
echo "=== 2. rRNA interval list ==="
if [[ ! -s "$RRNA_INTERVAL_LIST" ]]; then
    {
        cat "$GENOME_DICT"
        awk -F'\t' 'BEGIN{OFS="\t"}
            $0 !~ /^#/ && $3 == "gene" && $9 ~ /gene_biotype "(rRNA|Mt_rRNA)"/ {
                match($9, /gene_id "[^"]+"/)
                id = substr($9, RSTART + 9, RLENGTH - 10)
                print $1, $4, $5, $7, id
            }' "$GTF"
    } > "$RRNA_INTERVAL_LIST" || fail "rRNA interval list"
    N_RRNA=$(grep -vc '^@' "$RRNA_INTERVAL_LIST")
    echo "rRNA intervals: $N_RRNA"
    [[ "$N_RRNA" -eq 0 ]] && fail "no rRNA genes found in GTF"
    manifest_add "$RRNA_INTERVAL_LIST" "built from GTF" "gene_biotype rRNA + Mt_rRNA; $N_RRNA intervals"
else
    echo "exists, skipping"
fi
echo ""

# ── 3. refFlat ────────────────────────────────────────────────────────────────
# genePredExt column 12 (name2 = gene name) is moved to the front to give the
# 11-column refFlat layout Picard expects.
echo "=== 3. refFlat ==="
if [[ ! -s "$REFFLAT" ]]; then
    TMP_GP="${REFFLAT}.genePred.tmp"
    gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons \
        "$GTF" "$TMP_GP" || fail "gtfToGenePred"
    awk 'BEGIN{OFS="\t"} {print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' "$TMP_GP" > "$REFFLAT" \
        || fail "refFlat conversion"
    rm -f "$TMP_GP"
    manifest_add "$REFFLAT" "built from GTF" "$PKG_GTFTOGENEPRED gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons"
else
    echo "exists, skipping"
fi
echo ""

# ── 4. STAR index ─────────────────────────────────────────────────────────────
# --limitGenomeGenerateRAM is a resource ceiling only; it does not change the
# index. Raised from the STAR default because the human genome with GTF can
# exceed it. Keep it below the --mem requested above.
echo "=== 4. STAR index ($STAR_INDEX) ==="
if [[ -s "$STAR_INDEX/SAindex" && -s "$STAR_INDEX/sjdbList.out.tab" ]]; then
    echo "exists, skipping"
else
    rm -rf "$STAR_INDEX"
    mkdir -p "$STAR_INDEX"
    STAR \
        --runMode genomeGenerate \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$GENOME_FASTA" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang "$SJDB_OVERHANG" \
        --limitGenomeGenerateRAM 60000000000 \
        --outFileNamePrefix "$STAR_INDEX/"
    if [[ $? -ne 0 || ! -s "$STAR_INDEX/SAindex" ]]; then
        rm -rf "$STAR_INDEX"
        fail "STAR genomeGenerate (partial index removed)"
    fi
    manifest_add "$STAR_INDEX" "built" "$PKG_STAR genomeGenerate; sjdbOverhang $SJDB_OVERHANG; GTF release $ENSEMBL_RELEASE"
fi
echo ""

# ── 5. SortMeRNA index ────────────────────────────────────────────────────────
# Built once here so that pipeline array tasks only read it (--index 0).
# The pipeline must pass exactly these --ref paths when using the index.
echo "=== 5. SortMeRNA index ($SORTMERNA_IDX_DIR) ==="
if [[ -d "$SORTMERNA_IDX_DIR" && -n "$(ls -A "$SORTMERNA_IDX_DIR" 2>/dev/null)" ]]; then
    echo "exists, skipping"
else
    SMR_WORK="$SORTMERNA_DIR/work_index"
    rm -rf "$SMR_WORK" "$SORTMERNA_IDX_DIR"
    mkdir -p "$SMR_WORK" "$SORTMERNA_IDX_DIR"
    sortmerna \
        --ref "$SORTMERNA_DB_PROKARYOTE" \
        --ref "$SORTMERNA_DB_EUKARYOTE" \
        --index 1 \
        --idx-dir "$SORTMERNA_IDX_DIR" \
        --workdir "$SMR_WORK" \
        --threads "$THREADS"
    STATUS=$?
    rm -rf "$SMR_WORK"
    if [[ "$STATUS" -ne 0 ]]; then
        rm -rf "$SORTMERNA_IDX_DIR"
        fail "SortMeRNA indexing (partial index removed)"
    fi
    manifest_add "$SORTMERNA_IDX_DIR" "built" "$PKG_SORTMERNA --index 1"
fi
echo ""

# ── Summary ───────────────────────────────────────────────────────────────────
echo "=== Summary ==="
for F in "${GENOME_FASTA}.fai" "$GENOME_DICT" "$RRNA_INTERVAL_LIST" "$REFFLAT" \
         "$STAR_INDEX/SAindex" "$SORTMERNA_IDX_DIR"; do
    if [[ -e "$F" ]]; then echo "OK      : $F"; else echo "MISSING : $F"; fi
done
echo ""
echo "Chromosome naming check (first 3 names in each; they must match):"
echo "  FASTA   : $(cut -f1 "${GENOME_FASTA}.fai" | head -3 | tr '\n' ' ')"
echo "  GTF     : $(grep -v '^#' "$GTF" | cut -f1 | uniq | head -3 | tr '\n' ' ')"
echo "  refFlat : $(cut -f3 "$REFFLAT" | uniq | head -3 | tr '\n' ' ')"
echo ""
echo "=== Lines for the pipeline config.sh (REFERENCES AND DATABASES) ==="
echo "STAR_INDEX=\"$STAR_INDEX\""
echo "GTF=\"$GTF\""
echo "REFFLAT=\"$REFFLAT\""
echo "RRNA_INTERVAL_LIST=\"$RRNA_INTERVAL_LIST\""
echo "SORTMERNA_DB_PROKARYOTE=\"$SORTMERNA_DB_PROKARYOTE\""
echo "SORTMERNA_DB_EUKARYOTE=\"$SORTMERNA_DB_EUKARYOTE\""
echo "SORTMERNA_IDX_DIR=\"$SORTMERNA_IDX_DIR\""
echo ""
echo "Job done"; date
