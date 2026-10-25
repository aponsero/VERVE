#!/bin/bash
# ==============================================================================
# Reference configuration (module 01 references for VERVE)
#
# Standalone: sourced only by ref.01_download.sh and ref.02_build.sh. It does
# not read, and is not read by, the pipeline's config.sh. Both reference
# scripts share this file so that download and build cannot disagree on paths
# or versions.
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# SETTINGS
# ══════════════════════════════════════════════════════════════════════════════
REF_DIR="/qib/platforms/Informatics/AlisePonsero/0.Databases/VERVE_0.0.0/databases"
                                 # absolute path where references are stored
ENSEMBL_RELEASE="116"            # Ensembl release NUMBER (the N in release-N).
                                 # Several releases share GRCh38.p14, so the
                                 # number is needed to fix the GTF version
ENSEMBL_ASSEMBLY="GRCh38.p14"    # checked against the downloaded GTF header
SORTMERNA_DB_TAG="v4.3.7"        # sortmerna git tag the rRNA FASTA files come
                                 # from (SILVA release 119 derived)
SJDB_OVERHANG=100                # STAR default. Ideal is read length - 1; the
                                 # STAR manual states 100 works nearly as well.
                                 # Changing it requires rebuilding the index

TELESCOPE_BUILD="retro.hg38.v1"  # annotation build in mlbendall/telescope_annotation_db
                                 # (HERV_rmsk.hg38.v2 plus L1Base.hg38.v1)
TELESCOPE_REPO_REF="master"      # git ref the annotation is downloaded from. No
                                 # release tags exist, so the date and md5 in
                                 # MANIFEST.tsv are what pin this version


# ══════════════════════════════════════════════════════════════════════════════
# PACKAGES (Lmod, used by ref.02_build.sh only)
# ══════════════════════════════════════════════════════════════════════════════
# Must match the versions the pipeline runs with: a STAR index can only be used
# by the STAR version that built it.
MODULE_PATH_QIB="/qib/research-projects/bioboxes/lua"
PKG_STAR="star/2.7.10a"
PKG_SAMTOOLS="samtools/1.23.1"
PKG_PICARD="picard/3.4.0"
PKG_SORTMERNA="sortmerna/4.3.7"
PKG_GTFTOGENEPRED="ucsc-gtftogenepred/482-0"   # provides gtfToGenePred
PKG_BOWTIE2="bowtie2/2.5.4"      # for the ERV arm. Lives in the NBI module tree,
                                 # not in bioboxes, so no module use is needed


# ══════════════════════════════════════════════════════════════════════════════
# DERIVED PATHS (do not edit; change the settings above)
# ══════════════════════════════════════════════════════════════════════════════
ENSEMBL_DIR="$REF_DIR/ensembl_${ENSEMBL_RELEASE}"
GENOME_FASTA="$ENSEMBL_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GENOME_DICT="$ENSEMBL_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.dict"
GTF="$ENSEMBL_DIR/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf"
REFFLAT="$ENSEMBL_DIR/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.refFlat.txt"
RRNA_INTERVAL_LIST="$ENSEMBL_DIR/Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.rRNA.interval_list"
STAR_INDEX="$ENSEMBL_DIR/${PKG_STAR//\//_}_sjdb${SJDB_OVERHANG}"

SORTMERNA_DIR="$REF_DIR/sortmerna_db_${SORTMERNA_DB_TAG}"
SORTMERNA_DB_PROKARYOTE="$SORTMERNA_DIR/prokaryote_16s_23s.fasta"
SORTMERNA_DB_EUKARYOTE="$SORTMERNA_DIR/eukaryote_18s_28s.fasta"
SORTMERNA_IDX_DIR="$SORTMERNA_DIR/idx_${PKG_SORTMERNA//\//_}"

ERV_DIR="$REF_DIR/telescope_${TELESCOPE_BUILD}"
ERV_GTF_UCSC="$ERV_DIR/transcripts.ucsc.gtf"      # as downloaded
ERV_GTF="$ERV_DIR/transcripts.ensembl.gtf"        # renamed to Ensembl contigs
ERV_GENES_TSV="$ERV_DIR/genes.tsv.gz"             # locus-level summary table
BOWTIE2_INDEX="$ENSEMBL_DIR/${PKG_BOWTIE2//\//_}/genome"
                                 # index PREFIX, not a directory

MANIFEST="$REF_DIR/MANIFEST.tsv"
