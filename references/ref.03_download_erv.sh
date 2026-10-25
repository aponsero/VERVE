#!/bin/bash
# ==============================================================================
# Reference step 3: Telescope ERV annotation (module 03)
#
# Usage: bash ref.03_download_erv.sh      (from the references/ folder)
# Run INTERACTIVELY on the download node. Needs ref.02_build.sh to have run
# first, because the renaming step reads the genome .fai.
#
# Downloads the retro.hg38.v1 annotation (HERV_rmsk.hg38.v2 plus L1Base.hg38.v1)
# and renames its contigs from UCSC to Ensembl, so that the ERV arm uses the
# same genome as the rest of the pipeline.
#
# Renaming: chr1 to 1, chrX to X, chrM to MT. Loci on unplaced and alt contigs
# (chr1_KI270706v1_random and similar) are DROPPED rather than mapped, since
# that would need an NCBI assembly report; the count dropped is reported and
# recorded in the manifest.
# ==============================================================================

# No pipefail: the preview commands below pipe a large file into head, and head
# closing the pipe early would otherwise be treated as a failure and exit.
set -eu

source "$(dirname "$0")/ref.config.sh"

for VAR in REF_DIR TELESCOPE_BUILD TELESCOPE_REPO_REF; do
    if [[ -z "${!VAR}" ]]; then
        echo "ERROR: $VAR is empty in ref.config.sh"; exit 1
    fi
done

if [[ ! -s "${GENOME_FASTA}.fai" ]]; then
    echo "ERROR: ${GENOME_FASTA}.fai not found. Run ref.02_build.sh first."
    exit 1
fi

mkdir -p "$ERV_DIR"
TODAY=$(date +%Y-%m-%d)

manifest_add () {
    local MD5
    MD5=$(md5sum "$1" | awk '{print $1}')
    printf "%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$TODAY" "$MD5" "$3" >> "$MANIFEST"
}

# ── Download ──────────────────────────────────────────────────────────────────
BASE="https://github.com/mlbendall/telescope_annotation_db/raw/${TELESCOPE_REPO_REF}/builds/${TELESCOPE_BUILD}"

echo "=== Telescope annotation: $TELESCOPE_BUILD ==="
wget -q -O "$ERV_GTF_UCSC" "${BASE}/transcripts.gtf"
wget -q -O "$ERV_GENES_TSV" "${BASE}/genes.tsv.gz"

# These files are stored with Git LFS. If the server hands back the pointer
# stub instead of the content, the GTF is a few hundred bytes of text starting
# with "version https://git-lfs".
if head -c 60 "$ERV_GTF_UCSC" | grep -q "git-lfs"; then
    echo "ERROR: got a Git LFS pointer, not the annotation."
    echo "       Clone the repository with git lfs instead:"
    echo "         git clone https://github.com/mlbendall/telescope_annotation_db.git"
    echo "         cd telescope_annotation_db && git lfs pull"
    exit 1
fi

manifest_add "$ERV_GTF_UCSC" "${BASE}/transcripts.gtf" "UCSC contig names, as downloaded; ref $TELESCOPE_REPO_REF"
manifest_add "$ERV_GENES_TSV" "${BASE}/genes.tsv.gz" "locus-level summary table; ref $TELESCOPE_REPO_REF"

echo "Lines in downloaded GTF: $(grep -vc '^#' "$ERV_GTF_UCSC")"
echo ""
echo "First annotation line, to confirm the attribute keys:"
awk '!/^#/ { print; exit }' "$ERV_GTF_UCSC"
echo ""
echo "Columns in the locus summary table:"
gzip -cd "$ERV_GENES_TSV" | awk 'NR == 1 { print }' || true
echo ""

# ── Rename contigs to Ensembl ─────────────────────────────────────────────────
# Contigs present in the genome are kept; anything else is dropped.
echo "=== Renaming contigs ==="
cut -f1 "${GENOME_FASTA}.fai" > "$ERV_DIR/genome_contigs.txt"

awk -v contigs="$ERV_DIR/genome_contigs.txt" '
BEGIN {
    FS = OFS = "\t"
    while ((getline line < contigs) > 0) { known[line] = 1 }
}
/^#/ { next }
{
    name = $1
    sub(/^chr/, "", name)
    if (name == "M") { name = "MT" }
    if (name in known) { $1 = name; print; kept++ }
    else { dropped++; dropped_contigs[$1] = 1 }
}
END {
    printf "kept %d lines, dropped %d lines\n", kept, dropped > "/dev/stderr"
    n = 0
    for (c in dropped_contigs) { n++ }
    printf "contigs dropped: %d\n", n > "/dev/stderr"
}' "$ERV_GTF_UCSC" > "$ERV_GTF"

N_LOCI=$(grep -o 'locus "[^"]*"' "$ERV_GTF" | sort -u | wc -l)
echo "Distinct loci in the renamed annotation: $N_LOCI"

manifest_add "$ERV_GTF" "renamed from $(basename "$ERV_GTF_UCSC")" \
    "UCSC to Ensembl contigs; loci on unplaced and alt contigs dropped; $N_LOCI loci"

echo ""
echo "Next step, on the cluster, from the references/ folder: sbatch ref.04_build_erv.sh"
