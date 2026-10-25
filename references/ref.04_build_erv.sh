#!/bin/bash -l
#SBATCH --job-name=ref04_bowtie2
#SBATCH --partition=qib-compute
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=ref04_bowtie2_%j.out
#SBATCH --error=ref04_bowtie2_%j.err
# ==============================================================================
# Reference step 4: bowtie2 index for the ERV arm (module 03)
#
# Usage: sbatch ref.04_build_erv.sh      (from the references/ folder)
# Set --partition above before submitting. Run after ref.03_download_erv.sh.
#
# Built from the same primary assembly FASTA as the STAR index, so the ERV and
# host arms share one coordinate system.
# ==============================================================================

source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/ref.config.sh"

module use /qib/research-projects/bioboxes/lua
module load "$PKG_BOWTIE2"

echo "Job started"; hostname; date
echo ""

THREADS="${SLURM_CPUS_PER_TASK:-4}"
TODAY=$(date +%Y-%m-%d)

if [[ ! -s "$GENOME_FASTA" ]]; then
    echo "ERROR: $GENOME_FASTA not found. Has ref.01_download.sh completed?"
    exit 1
fi
if [[ ! -s "$ERV_GTF" ]]; then
    echo "ERROR: $ERV_GTF not found. Has ref.03_download_erv.sh completed?"
    exit 1
fi

echo "=== bowtie2 index ($BOWTIE2_INDEX) ==="
if [[ -s "${BOWTIE2_INDEX}.1.bt2" || -s "${BOWTIE2_INDEX}.1.bt2l" ]]; then
    echo "exists, skipping"
else
    mkdir -p "$(dirname "$BOWTIE2_INDEX")"
    bowtie2-build --threads "$THREADS" "$GENOME_FASTA" "$BOWTIE2_INDEX"
    if [[ $? -ne 0 ]]; then
        echo "ERROR: bowtie2-build failed (partial index left in place for inspection)"
        exit 1
    fi
    printf "%s\t%s\t%s\t%s\t%s\n" "$BOWTIE2_INDEX" "built" "$TODAY" "NA" \
        "$PKG_BOWTIE2 bowtie2-build on $(basename "$GENOME_FASTA")" >> "$MANIFEST"
fi

echo ""
echo "=== Contig naming check (ERV annotation against genome) ==="
echo "Annotation contigs, first 5: $(cut -f1 "$ERV_GTF" | uniq | head -5 | tr '\n' ' ')"
echo "Genome contigs, first 5    : $(cut -f1 "${GENOME_FASTA}.fai" | head -5 | tr '\n' ' ')"

MISSING=$(cut -f1 "$ERV_GTF" | sort -u | grep -vxF -f <(cut -f1 "${GENOME_FASTA}.fai") | wc -l)
echo "Annotation contigs absent from the genome: $MISSING (should be 0)"

echo ""
echo "=== Lines for the pipeline config.sh ==="
echo "BOWTIE2_INDEX=\"$BOWTIE2_INDEX\""
echo "ERV_GTF=\"$ERV_GTF\""
echo "ERV_GENES_TSV=\"$ERV_GENES_TSV\""
echo ""
echo "Job done"; date
