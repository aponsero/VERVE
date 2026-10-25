#!/bin/bash
# ==============================================================================
# VERVE pipeline configuration
# Viral Expression and Retroviral Element pipeline
#
# This file is sourced by every submission script, module script and check
# script. Nothing is hardcoded elsewhere.
#
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# INPUT / OUTPUT
# ══════════════════════════════════════════════════════════════════════════════
# All paths must be absolute.

SAMPLE_SHEET="samples.tsv"                  # absolute path to samples.tsv
OUT_DIR="VERVE_test"                       # absolute path to the results directory
LOG_DIR="$OUT_DIR/logs"          # override with an absolute path if you want
                                 # logs somewhere other than inside OUT_DIR


# ══════════════════════════════════════════════════════════════════════════════
# MODULE SWITCHES
# ══════════════════════════════════════════════════════════════════════════════
RUN_FASTQC=true                  # FastQC on trimmed reads, in addition to the
                                 # fastp report. Adds per-base metrics that
                                 # fastp does not report. MultiQC aggregates
                                 # both; it runs regardless of this setting.


# ══════════════════════════════════════════════════════════════════════════════
# MODULE 00a -- PREPROCESSING PARAMETERS (fastp)
# ══════════════════════════════════════════════════════════════════════════════

FASTP_QUALIFIED_QUALITY=15       # fastp default. Phred score at or above which
                                 # a base is called qualified
FASTP_UNQUALIFIED_PERCENT=40     # fastp default. Read discarded if more than
                                 # this percent of bases are unqualified
FASTP_N_BASE_LIMIT=5             # fastp default. Read discarded above this
                                 # many N bases
FASTP_LENGTH_REQUIRED=15         # fastp default. Reads shorter than this after
                                 # trimming are discarded.
                                 # REVISIT BEFORE MODULE 04: short reads are a
                                 # known false-positive source in competitive
                                 # viral mapping. Raising this is a design
                                 # decision, not a default to change quietly
FASTP_DETECT_ADAPTER_PE=true     # NOT a fastp default. Enables adapter
                                 # auto-detection for paired-end data in
                                 # addition to overlap-based detection. Set to
                                 # false and use FASTP_ADAPTER_R1/R2 if the
                                 # adapter sequences are known from the
                                 # sequencing provider
FASTP_ADAPTER_R1=""              # leave empty unless known; overrides detection
FASTP_ADAPTER_R2=""              # leave empty unless known; overrides detection
FASTP_EXTRA_ARGS=""              # any additional fastp arguments, appended last


# ══════════════════════════════════════════════════════════════════════════════
# MODULE 00b -- AGGREGATE REPORTING (MultiQC)
# ══════════════════════════════════════════════════════════════════════════════
MULTIQC_TITLE="VERVE module 00 preprocessing"
MULTIQC_EXTRA_ARGS=""


# ══════════════════════════════════════════════════════════════════════════════
# MODULE 01 -- HOST ALIGNMENT AND FEASIBILITY QC
# ══════════════════════════════════════════════════════════════════════════════
# Library strandedness. Set explicitly by the user; never inferred by the
# pipeline. Leave EMPTY on the first run of new data: 01.run_host.sh then
# submits only alignment (01a) and rRNA measurement (01b), and
# 01.check_host.sh prints the STAR ReadsPerGene column totals to decide from.
# Once set, rerun 01.run_host.sh to submit 01c and 01d (01a and 01b skip
# samples that are already done).
#
# Allowed values:
#   unstranded : featureCounts -s 0 | Picard NONE
#   forward    : featureCounts -s 1 | Picard FIRST_READ_TRANSCRIPTION_STRAND
#                (ReadsPerGene column 3 dominant)
#   reverse    : featureCounts -s 2 | Picard SECOND_READ_TRANSCRIPTION_STRAND
#                (ReadsPerGene column 4 dominant)
STRANDEDNESS="reverse"

# STAR: alignment parameters are left at STAR defaults. The module script sets
# output options only (sorted BAM, GeneCounts, unmapped reads to FASTQ).
STAR_EXTRA_ARGS=""

# featureCounts: left at defaults (exon features, gene_id meta-features,
# multi-mapping and multi-overlapping reads not counted, duplicates counted).
FEATURECOUNTS_EXTRA_ARGS=""

# Java heap for Picard. Keep below the SLURM memory request for 01c, to leave
# room for JVM overhead.
PICARD_XMX="12g"

MULTIQC_TITLE_HOST="VERVE module 01 host alignment and feasibility QC"


# ══════════════════════════════════════════════════════════════════════════════
# PACKAGES (Lmod)
# ══════════════════════════════════════════════════════════════════════════════
# Tools are loaded with:
#     module use $MODULE_PATH_QIB
#     module load $PKG_<TOOL>

MODULE_PATH_QIB="/qib/research-projects/bioboxes/lua"

PKG_FASTP="fastp/1.3.1"
PKG_FASTQC="fastqc/0.12.1"
PKG_MULTIQC="multiqc/1.33"
PKG_STAR="star/2.7.10a"
PKG_SAMTOOLS="samtools/1.23.1"
PKG_PICARD="picard/3.4.0"
PKG_SUBREAD="subread/2.0.1"          # featureCounts. NOTE: in 2.0.1, -p alone
                                     # counts fragments. From 2.0.2 onwards
                                     # --countReadPairs is also required; the
                                     # 01d script must change if this is upgraded
PKG_SORTMERNA="sortmerna/4.3.7"
PKG_BOWTIE2="bowtie2/2.5.4"          # ERV arm. Lives in the NBI module tree,
                                     # not in bioboxes, so no module use needed
PKG_TELESCOPE="telescope/1.0.4.1"    # the bioconda build of the hanalysis fork,
                                     # which is the one with --stranded_mode


# ══════════════════════════════════════════════════════════════════════════════
# REFERENCES AND DATABASES
# ══════════════════════════════════════════════════════════════════════════════
# Produced by the independent reference scripts (references/ folder), which do
# not read this file. ref.02_build.sh prints the exact lines to paste here.
# Provenance (URLs, dates, checksums) is in the MANIFEST.tsv of the reference
# directory.
#
# Consistency requirements for module 01:
#   - STAR_INDEX, GTF and REFFLAT from the same Ensembl release
#   - identical chromosome naming across all three
#   - STAR index built WITH the GTF (required by --quantMode GeneCounts)
#   - SortMeRNA database paths identical to those used to build the index

STAR_INDEX="VERVE_0.0.0/databases/ensembl_116/star_2.7.10a_sjdb100"                    # directory
GTF="VERVE_0.0.0/databases/ensembl_116/Homo_sapiens.GRCh38.116.gtf"                           # Ensembl GRCh38 GTF (uncompressed)
REFFLAT="VERVE_0.0.0/databases/ensembl_116/Homo_sapiens.GRCh38.116.refFlat.txt"                       # refFlat for Picard CollectRnaSeqMetrics
RRNA_INTERVAL_LIST="VERVE_0.0.0/databases/ensembl_116/Homo_sapiens.GRCh38.116.rRNA.interval_list"            # OPTIONAL. Picard interval list of rRNA and
                                 # Mt_rRNA genes. Leave empty to skip. 45S rDNA
                                 # is poorly represented in GRCh38, so this
                                 # underestimates human rRNA

SORTMERNA_DB_PROKARYOTE="VERVE_0.0.0/databases/sortmerna_db_v4.3.7/prokaryote_16s_23s.fasta"       # SILVA bacteria + archaea 16S and 23S
SORTMERNA_DB_EUKARYOTE="VERVE_0.0.0/databases/sortmerna_db_v4.3.7/eukaryote_18s_28s.fasta"        # SILVA eukaryote 18S + 28S. PROXY for human
                                 # rRNA: fungal rRNA also lands here, so the
                                 # fraction may be inflated at sites with fungal
                                 # load (vaginal in particular)
SORTMERNA_IDX_DIR="VERVE_0.0.0/databases/sortmerna_db_v4.3.7/idx_sortmerna_4.3.7"             # prebuilt SortMeRNA index for the two
                                 # databases above

BOWTIE2_INDEX="VERVE_0.0.0/databases/ensembl_116/bowtie2_2.5.4/genome"                 # index PREFIX (not a directory), built from
                                 # the same genome FASTA as the STAR index
ERV_GTF="VERVE_0.0.0/databases/telescope_retro.hg38.v1/transcripts.ensembl.gtf"                       # Telescope annotation, contigs renamed to
                                 # Ensembl (transcripts.ensembl.gtf)

VIRAL_INDEX="VERVE_0.0.0/databases/viral_virushostdb_22.10.2026/bowtie2_2.5.4/grch38_plus_viral"                   # index PREFIX for GRCh38 plus the viral set in
                                 # ONE index, so mapping is competitive and a
                                 # read that is really human prefers human.
                                 # Separate from BOWTIE2_INDEX above
VIRAL_TAXA="VERVE_0.0.0/databases/viral_virushostdb_22.10.2026/taxa.tsv"                    # taxa.tsv: accession to taxid, genus, family,
                                 # segment number, category and flags
# ══════════════════════════════════════════════════════════════════════════════
# MODULE 03 -- ERV QUANTIFICATION (Telescope)
# ══════════════════════════════════════════════════════════════════════════════
# bowtie2 settings follow the Telescope tutorial. They are deliberately
# permissive: Telescope's EM model reassigns ambiguous fragments, so it needs
# those alignments to exist rather than being discarded at the aligner.
#   --very-sensitive-local  thorough local alignment
#   -k 100                  report up to 100 alignments per read, instead of 1
#   --score-min L,0,1.6     accept alignments of roughly 95% identity
# Raising -k increases BAM size and Telescope runtime; lowering it discards
# the ambiguity the model exists to resolve.
BOWTIE2_K=100
BOWTIE2_SCORE_MIN="L,0,1.6"
BOWTIE2_PRESET="--very-sensitive-local"
BOWTIE2_EXTRA_ARGS=""

# Telescope. RF matches the reverse-stranded libraries (read 1 on the reverse
# strand), the same convention as featureCounts -s 2. The default is None,
# which reports both strands per feature and counts antisense overlap.
TELESCOPE_STRANDED_MODE="RF"
TELESCOPE_ATTRIBUTE="locus"          # Telescope default; the attribute key in
                                     # the retro.hg38.v1 annotation
TELESCOPE_EXTRA_ARGS=""              # reassign mode, priors and overlap
                                     # settings are left at Telescope defaults;
                                     # all reassignment modes appear in the
                                     # report regardless

# The -k 100 BAMs are large. Set to true to delete a sample's BAM once its
# Telescope run has succeeded; rerunning 03b then needs a realignment.
ERV_DELETE_BAM_AFTER_ASSIGN=false


# ══════════════════════════════════════════════════════════════════════════════
# MODULE 04 -- RNA VIRAL DETECTION
# ══════════════════════════════════════════════════════════════════════════════
# Input is the reads STAR could not place on the human genome (written by 01a).
#
# Unlike module 03 this is BEST-HIT mapping: no -k, so bowtie2 reports one
# alignment per read and MAPQ stays meaningful. Local alignment is kept for
# sensitivity to divergent strains, since the reference carries one genome per
# taxon and a strain some percent away still needs to be found.
BOWTIE2_VIRAL_PRESET="--very-sensitive-local"
BOWTIE2_VIRAL_EXTRA_ARGS=""

# Reads are counted twice: all primary alignments, and those at or above this
# MAPQ. The unfiltered count is the primary number, because near-identical
# strains of the same virus legitimately split reads and lower their MAPQ;
# that ambiguity dissolves on aggregation to genus and family. The filtered
# count is reported alongside so the difference is visible.
VIRAL_MIN_MAPQ=10

# Window size for the spread-of-coverage metric. Breadth alone is not
# comparable across genome sizes (ten reads cover a lot of a 7 kb enterovirus
# and almost none of a 150 kb herpesvirus), so the number of distinct covered
# windows is reported next to it.
VIRAL_WINDOW_SIZE=1000

# The BAM holds only reads that mapped to a viral sequence, so it is small and
# kept by default; it is what any positive call has to be checked against.
VIRAL_KEEP_BAM=true

# ══════════════════════════════════════════════════════════════════════════════
# SLURM -- MODULE 00a (fastp + FastQC, array)
# ══════════════════════════════════════════════════════════════════════════════
PARTITION_PREPROCESS="qib-compute"
TIME_PREPROCESS="02:00:00"
CPUS_PREPROCESS=8                # fastp does not scale usefully beyond 16
MEM_PER_CPU_PREPROCESS="4G"


# ══════════════════════════════════════════════════════════════════════════════
# SLURM -- MODULE 00b (MultiQC, single job)
# ══════════════════════════════════════════════════════════════════════════════
PARTITION_PREPROCESS_MERGE="qib-compute"
TIME_PREPROCESS_MERGE="01:00:00"
CPUS_PREPROCESS_MERGE=1
MEM_PER_CPU_PREPROCESS_MERGE="8G"


# ══════════════════════════════════════════════════════════════════════════════
# SLURM -- MODULE 01
# ══════════════════════════════════════════════════════════════════════════════

# 01a: STAR alignment (array). GRCh38 index needs about 32 GB resident.
PARTITION_ALIGN="qib-compute"
TIME_ALIGN="06:00:00"
CPUS_ALIGN=8
MEM_PER_CPU_ALIGN="5G"

# 01b: SortMeRNA on all reads (array). Slow on full libraries.
PARTITION_RRNA="qib-compute"
TIME_RRNA="24:00:00"
CPUS_RRNA=8
MEM_PER_CPU_RRNA="4G"

# 01c: Picard + samtools (array)
PARTITION_BAMQC="qib-compute"
TIME_BAMQC="04:00:00"
CPUS_BAMQC=2
MEM_PER_CPU_BAMQC="8G"

# 01d: featureCounts on all samples + MultiQC (single job)
PARTITION_HOST_MERGE="qib-compute"
TIME_HOST_MERGE="04:00:00"
CPUS_HOST_MERGE=8
MEM_PER_CPU_HOST_MERGE="4G"


# ══════════════════════════════════════════════════════════════════════════════
# SLURM -- MODULE 03
# ══════════════════════════════════════════════════════════════════════════════
# Starting values, to revise on the pilot. The vaginal samples (50 to 70%
# multi-mapping) are the ones to size on.

# 03a: permissive bowtie2 alignment (array)
PARTITION_ERV_ALIGN="qib-compute"
TIME_ERV_ALIGN="24:00:00"
CPUS_ERV_ALIGN=8
MEM_PER_CPU_ERV_ALIGN="4G"

# 03b: telescope assign (array). Telescope states multiple cores are not
# supported, so this is single-core and slow, with memory the binding
# constraint rather than CPU.
PARTITION_ERV_ASSIGN="qib-compute"
TIME_ERV_ASSIGN="48:00:00"
CPUS_ERV_ASSIGN=1
MEM_ERV_ASSIGN="64G"

# 03c: merge reports into count matrices (single job)
PARTITION_ERV_MERGE="qib-compute"
TIME_ERV_MERGE="02:00:00"
CPUS_ERV_MERGE=1
MEM_PER_CPU_ERV_MERGE="8G"


# ══════════════════════════════════════════════════════════════════════════════
# SLURM -- MODULE 04
# ══════════════════════════════════════════════════════════════════════════════
# 04a: bowtie2 against the combined index, then per-sample counting (array)
PARTITION_VIRAL="qib-compute"
TIME_VIRAL="12:00:00"
CPUS_VIRAL=8
MEM_PER_CPU_VIRAL="4G"

# 04b: merge into matrices (single job)
PARTITION_VIRAL_MERGE="qib-compute"
TIME_VIRAL_MERGE="02:00:00"
CPUS_VIRAL_MERGE=1
MEM_PER_CPU_VIRAL_MERGE="8G"
