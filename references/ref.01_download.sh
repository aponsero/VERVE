#!/bin/bash
# ==============================================================================
# Reference step 1: download (module 01 references for VERVE)
#
# Usage: bash ref.01_download.sh      (from the references/ folder)
# Run INTERACTIVELY on the download node (internet access, no modules needed).
# Uses only wget, gunzip, sum, md5sum and standard shell tools.
#
# Standalone: reads ref.config.sh only, never the pipeline's config.sh.
#
# Downloads:
#   - Ensembl GRCh38 primary assembly FASTA and matching GTF, verified against
#     Ensembl's CHECKSUMS files
#   - SortMeRNA rRNA FASTA files from the sortmerna repository, concatenated
#     into a prokaryotic and a eukaryotic database (5S and 5.8S not used)
#
# Every file is recorded in $MANIFEST with source URL, date and md5.
# ==============================================================================

set -euo pipefail

source "$(dirname "$0")/ref.config.sh"

# ── Sanity checks ─────────────────────────────────────────────────────────────
for VAR in REF_DIR ENSEMBL_RELEASE ENSEMBL_ASSEMBLY SORTMERNA_DB_TAG; do
    if [[ -z "${!VAR}" ]]; then
        echo "ERROR: $VAR is empty in ref.config.sh"
        exit 1
    fi
done

mkdir -p "$ENSEMBL_DIR" "$SORTMERNA_DIR/raw"

TODAY=$(date +%Y-%m-%d)
if [[ ! -f "$MANIFEST" ]]; then
    printf "file\tsource\tdate\tmd5\tnote\n" > "$MANIFEST"
fi

manifest_add () {
    # $1 file, $2 source, $3 note
    local MD5
    MD5=$(md5sum "$1" | awk '{print $1}')
    printf "%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$TODAY" "$MD5" "$3" >> "$MANIFEST"
}

# ══════════════════════════════════════════════════════════════════════════════
# Ensembl
# ══════════════════════════════════════════════════════════════════════════════
ENS_BASE="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}"
FA_URL="$ENS_BASE/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_URL="$ENS_BASE/gtf/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf.gz"

echo "=== Ensembl release $ENSEMBL_RELEASE ==="

cd "$ENSEMBL_DIR"

wget -q -O CHECKSUMS.fasta "$ENS_BASE/fasta/homo_sapiens/dna/CHECKSUMS"
wget -q -O CHECKSUMS.gtf   "$ENS_BASE/gtf/homo_sapiens/CHECKSUMS"

download_ensembl () {
    # $1 url, $2 checksum file. Verified with the BSD 'sum' algorithm used by
    # Ensembl's CHECKSUMS files.
    local URL="$1" CK="$2" GZ
    GZ=$(basename "$URL")
    echo "Downloading $GZ"
    wget -c -q --show-progress -O "$GZ" "$URL"

    local EXPECTED OBSERVED
    EXPECTED=$(awk -v f="$GZ" '$3 == f {print $1, $2}' "$CK")
    OBSERVED=$(sum "$GZ" | awk '{print $1, $2}')
    if [[ -z "$EXPECTED" ]]; then
        echo "ERROR: $GZ not listed in $CK"; exit 1
    fi
    if [[ "$EXPECTED" != "$OBSERVED" ]]; then
        echo "ERROR: checksum mismatch for $GZ (expected '$EXPECTED', got '$OBSERVED')"
        echo "       delete the file and rerun"
        exit 1
    fi
    echo "  checksum OK"
    manifest_add "$ENSEMBL_DIR/$GZ" "$URL" "Ensembl CHECKSUMS verified (sum: $OBSERVED)"
}

download_ensembl "$FA_URL"  CHECKSUMS.fasta
download_ensembl "$GTF_URL" CHECKSUMS.gtf

echo "Decompressing (compressed copies are kept)"
gunzip -c "$(basename "$FA_URL")"  > "$GENOME_FASTA"
gunzip -c "$(basename "$GTF_URL")" > "$GTF"

# ── Assembly check against the GTF header ─────────────────────────────────────
GTF_BUILD=$(grep -m1 '^#!genome-build ' "$GTF" | awk '{print $2}')
echo "GTF genome-build: $GTF_BUILD (expected $ENSEMBL_ASSEMBLY)"
if [[ "$GTF_BUILD" != "$ENSEMBL_ASSEMBLY" ]]; then
    echo "ERROR: release $ENSEMBL_RELEASE is $GTF_BUILD, not $ENSEMBL_ASSEMBLY"
    exit 1
fi

GTF_HEADER=$(grep '^#!' "$GTF" | tr '\n' ' ')
manifest_add "$GENOME_FASTA" "$FA_URL" "decompressed; assembly $GTF_BUILD"
manifest_add "$GTF" "$GTF_URL" "decompressed; header: $GTF_HEADER"

# ══════════════════════════════════════════════════════════════════════════════
# SortMeRNA rRNA databases
# ══════════════════════════════════════════════════════════════════════════════
# Separate per-kingdom files from the sortmerna repository. Per the bundled
# README these are SILVA release 119 derived (clustered representatives).
echo ""
echo "=== SortMeRNA databases ($SORTMERNA_DB_TAG) ==="

SMR_BASE="https://raw.githubusercontent.com/sortmerna/sortmerna/${SORTMERNA_DB_TAG}/data/rRNA_databases"
PROK_FILES=(silva-bac-16s-id90 silva-bac-23s-id98 silva-arc-16s-id95 silva-arc-23s-id98)
EUK_FILES=(silva-euk-18s-id95 silva-euk-28s-id98)

cd "$SORTMERNA_DIR/raw"
for F in "${PROK_FILES[@]}" "${EUK_FILES[@]}"; do
    echo "Downloading ${F}.fasta"
    wget -q -O "${F}.fasta" "$SMR_BASE/${F}.fasta"
    if [[ ! -s "${F}.fasta" ]] || ! head -c 1 "${F}.fasta" | grep -q '>'; then
        echo "ERROR: ${F}.fasta is empty or not FASTA"; exit 1
    fi
    manifest_add "$SORTMERNA_DIR/raw/${F}.fasta" "$SMR_BASE/${F}.fasta" \
        "sortmerna repo tag $SORTMERNA_DB_TAG; SILVA 119 derived; no upstream checksum"
done
wget -q -O README.txt "$SMR_BASE/README.txt"

echo "Concatenating"
cat $(printf "%s.fasta " "${PROK_FILES[@]}") > "$SORTMERNA_DB_PROKARYOTE"
cat $(printf "%s.fasta " "${EUK_FILES[@]}")  > "$SORTMERNA_DB_EUKARYOTE"

manifest_add "$SORTMERNA_DB_PROKARYOTE" "concatenation" "${PROK_FILES[*]}"
manifest_add "$SORTMERNA_DB_EUKARYOTE"  "concatenation" "${EUK_FILES[*]}"

# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "=== Summary ==="
echo "  Genome FASTA : $GENOME_FASTA ($(grep -c '^>' "$GENOME_FASTA") sequences)"
echo "  GTF          : $GTF"
echo "  Prokaryote   : $SORTMERNA_DB_PROKARYOTE ($(grep -c '^>' "$SORTMERNA_DB_PROKARYOTE") sequences)"
echo "  Eukaryote    : $SORTMERNA_DB_EUKARYOTE ($(grep -c '^>' "$SORTMERNA_DB_EUKARYOTE") sequences)"
echo "  Manifest     : $MANIFEST"
echo ""
echo "Next step, on the cluster, from the references/ folder: sbatch ref.02_build.sh"
