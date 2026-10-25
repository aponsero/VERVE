# VERVE

**V**iral **E**xpression and **R**etroviral **V**ariant **E**valuation: a SLURM
pipeline for quantifying endogenous retrovirus expression and detecting RNA
viruses in human swab RNA-seq libraries.

Built for the ME virome and ERV study at QIB. Version 0.0.0.

---

## Overview

VERVE takes paired-end, stranded total RNA-seq from human swabs and produces
two independent readouts from the same libraries:

- **ERV expression**: locus- and family-level counts for 27,748 endogenous
  retrovirus and LINE-1 elements, normalised with size factors from cellular
  genes.
- **RNA viral detection**: reads, coverage breadth and RPM for 1,873 human-host
  virus sequences, with GRCh38 in the same index as a decoy.

It also produces the host gene expression and QC that both arms depend on, and
a feasibility assessment of whether the libraries support locus-level ERV work
at all.

The pipeline stops at module 04. Contamination and batch control (module 05 in
the original design) and the statistics (module 06) are downstream analyses
whose shape depends on the study design, and they are not part of the
pipeline.

### What VERVE does not do

- No positivity calls, exclusions or statistical tests. Every report labels
  and flags; nothing is filtered out of a pipeline output.
- No de novo virus discovery.
- No deduplication of ERV or viral alignments.
- No cell-type deconvolution (a coarse epithelial versus immune marker score
  is computed instead).

---

## Pipeline structure

```
config.sh                      all parameters, sourced by every script
samples.tsv                    the sample sheet

00.run_preprocess.sh           submit module 00
00.check_preprocess.sh         verify module 00 outputs
01.run_host.sh                 submit module 01
01.check_host.sh               verify module 01 outputs
01.qc_report.Rmd               feasibility QC report and flags
02.make_markers.py             build the marker panel (run once)
02.expression_report.Rmd       size factors and composition scores
03.run_erv.sh                  submit module 03
03.check_erv.sh                verify module 03 outputs
03.erv_report.Rmd              ERV report
04.run_viral.sh                submit module 04
04.check_viral.sh              verify module 04 outputs
04.viral_report.Rmd            viral detection report

modules/                       the SLURM job scripts, not called directly
  00a_preprocess.sh            fastp and FastQC (array)
  00b_preprocess_multiqc.sh    MultiQC (merge)
  01a_align.sh                 STAR (array)
  01b_rrna.sh                  SortMeRNA (array)
  01c_bamqc.sh                 Picard and samtools (array)
  01d_host_merge.sh            featureCounts and MultiQC (merge)
  03a_erv_align.sh             permissive bowtie2 (array)
  03b_telescope.sh             telescope assign (array)
  03c_erv_merge.sh             merge and MultiQC
  04a_viral_map.sh             bowtie2 and counting (array)
  04b_viral_merge.sh           merge and MultiQC

scripts/                       helpers
  merge_telescope.py           Telescope reports to count matrices
  merge_viral.py               viral metrics to count matrices
  erv_locus_diagnostic.sh      inspect one ERV locus in one sample
  viral_taxon_diagnostic.sh    inspect one viral reference in one sample

references/                    reference building, independent of the pipeline
  ref.config.sh                its own settings; does not read config.sh
  ref.01_download.sh           Ensembl genome, GTF, SortMeRNA databases
  ref.02_build.sh              STAR index, refFlat, rRNA intervals
  ref.03_download_erv.sh       Telescope annotation, contigs renamed
  ref.04_build_erv.sh          bowtie2 index for the ERV arm
  ref.05_make_viral_set.py     filter VirusHostDB to the viral set
  ref.06_build_viral_index.sh  bowtie2 index for the viral arm
```

Each step follows the same pattern: `NN.run_*.sh` submits, `NN.check_*.sh`
verifies, and where interpretation is needed an `.Rmd` report is knitted in
RStudio. Array jobs skip samples that are already done, so any `run` script is
safe to resubmit after a partial failure.

---

## Tools and versions

Loaded through Lmod from `/qib/research-projects/bioboxes/lua`, except bowtie2
which lives in the NBI tree.

| Tool | Version | Used for |
| --- | --- | --- |
| fastp | 1.3.1 | adapter and quality trimming (00a) |
| FastQC | 0.12.1 | per-base QC of trimmed reads (00a) |
| MultiQC | 1.33 | aggregate reports (00b, 01d, 03c, 04b) |
| STAR | 2.7.10a | spliced alignment to GRCh38 (01a) |
| SortMeRNA | 4.3.7 | rRNA fraction, measurement only (01b) |
| Picard | 3.4.0 | MarkDuplicates, CollectRnaSeqMetrics (01c) |
| samtools | 1.23.1 | BAM handling, coverage, depth (throughout) |
| subread (featureCounts) | 2.0.1 | gene counts (01d) |
| bowtie2 | 2.5.4 | permissive ERV alignment (03a), viral mapping (04a) |
| Telescope | 1.0.4.1 | EM reassignment of ERV reads (03b) |
| ucsc-gtftogenepred | 482-0 | refFlat for Picard (reference build) |

R packages for the reports: tidyverse, knitr, DESeq2, GSVA. Python for the
helper scripts: standard library plus openpyxl for `ref.05`.

### References

| Reference | Version | Notes |
| --- | --- | --- |
| Human genome | Ensembl 116, GRCh38.p14 | primary assembly, not toplevel |
| Gene annotation | Ensembl 116 GTF | same release as the STAR index |
| ERV annotation | `retro.hg38.v1` | HERV_rmsk.hg38.v2 plus L1Base.hg38.v1 |
| rRNA databases | sortmerna v4.3.7 | SILVA 119 derived, 5S and 5.8S excluded |
| Viral genomes | VirusHostDB, human subset | 1,873 sequences, 1,389 taxa |

Every downloaded and built file is recorded in `MANIFEST.tsv` in the reference
directory, with its source URL, date and md5. The compute nodes are
air-gapped, so a database cannot be re-pulled later to check what it was.

---

## Setting up

### 1. Build the references

Run once per installation. These scripts are standalone: they read
`references/ref.config.sh` and never touch the pipeline's `config.sh`.

Fill in `REF_DIR` and `ENSEMBL_RELEASE` in `references/ref.config.sh`, then:

```bash
cd references

# On the download node (internet access, no modules available)
bash ref.01_download.sh

# On the cluster. Set --partition in each script first.
sbatch ref.02_build.sh

# Back on the download node
bash ref.03_download_erv.sh

# On the cluster
sbatch ref.04_build_erv.sh
```

For the viral arm, download the VirusHostDB human-host subset table and the
complete genome FASTA yourself, then:

```bash
python3 ref.05_make_viral_set.py \
    --subset virushostdb_humanonly.xlsx \
    --fasta virushostdb.genomic.fna.gz \
    --out-dir $REF_DIR/viral_virushostdb_<date>

sbatch ref.06_build_viral_index.sh
```

`ref.02`, `ref.04` and `ref.06` each print the exact lines to paste into the
pipeline's `config.sh` when they finish.

Expect a few hours for each bowtie2 index and for the STAR index.

### 2. Configure the pipeline

Edit `config.sh`. The sections that must be filled in:

- `SAMPLE_SHEET` and `OUT_DIR`, both absolute paths. `OUT_DIR` can be
  anywhere; it is deliberately not tied to the script location.
- The reference paths, pasted from the `ref.*` output.
- The `PARTITION_*` values. Several are deliberately left empty so a run
  cannot start against the wrong queue by accident.

### 3. Write the sample sheet

Tab-separated, with a header:

```
sample_id    read1              read2              body_site    group      batch
CM28_Nas     /path/R1.fastq.gz  /path/R2.fastq.gz  nasal        control    1
Neg          /path/R1.fastq.gz  /path/R2.fastq.gz  nasal        blank      1
```

Rules:

- No spaces in `sample_id`.
- `body_site` is used to group samples for normalisation. Keep it consistent:
  `Nasal` and `nasal` are two different sites to the check scripts, and a
  group of one sample silently receives a size factor of exactly 1.
- `group` beginning with `blank` marks a negative control. Blanks run through
  the whole pipeline identically to samples, but are excluded from
  normalisation and from the plots, and reported separately.
- Columns beyond the first three are not read by modules 00 to 04 except as
  described above; they are there for the downstream analysis.
- Excel on Windows adds carriage returns. Every `run` script strips them from
  the sample sheet before submitting anything.

---

## Running

Each step in turn, waiting for the jobs to finish before moving on.

```bash
bash 00.run_preprocess.sh
bash 00.check_preprocess.sh

bash 01.run_host.sh          # phase 1: alignment and rRNA only
bash 01.check_host.sh        # prints the strandedness table
# set STRANDEDNESS in config.sh, then
bash 01.run_host.sh          # phase 2: adds QC, counts and the merge
bash 01.check_host.sh
```

**Module 01 runs in two phases.** The pipeline never infers strandedness. On
the first pass only alignment and rRNA measurement are submitted, and the
check script prints the STAR `ReadsPerGene` column totals with a column-3
share. You read that, set `STRANDEDNESS` to `unstranded`, `forward` or
`reverse` in `config.sh`, and rerun. The already-finished samples are skipped.

Then knit the two reports in RStudio, setting `out_dir` and `sample_sheet` in
each YAML header:

```
01.qc_report.Rmd             feasibility flags; writes qc_flags.tsv
02.expression_report.Rmd     size factors and marker scores
```

`02.expression_report.Rmd` needs a marker panel. Build it once:

```bash
python3 02.make_markers.py \
    --cellmarker Cell_marker_Human.xlsx \
    --gtf <GTF> \
    --out markers.tsv
```

`markers.tsv` should be version-controlled with the pipeline. The script
prints its full selection to stderr; keep that output alongside the file.

Then the two arms, which are independent of each other:

```bash
bash 03.run_erv.sh
bash 03.check_erv.sh
# knit 03.erv_report.Rmd

bash 04.run_viral.sh
bash 04.check_viral.sh
# knit 04.viral_report.Rmd
```

Module 03 requires `size_factors.tsv` from the module 02 report, so knit that
first.

---

## Output structure

```
$OUT_DIR/
  00_preprocess/
    trimmed/                  <sample>_R{1,2}.trimmed.fastq.gz
    fastp/                    per-sample JSON and HTML
    fastqc/                   per-sample zip and HTML
    multiqc/                  multiqc_report.html
  01_host/
    align/                    STAR BAMs, indexes, logs, ReadsPerGene
    unmapped/                 non-host reads, the input to module 04
    rrna/                     SortMeRNA logs
    markdup/                  duplicate-marked BAMs
    bamqc/                    Picard and samtools metrics
    counts/                   gene_counts.txt and its summary
    multiqc/                  multiqc_report.html
    qc_flags.tsv              written by 01.qc_report.Rmd
  02_expression/
    size_factors.tsv          consumed by module 03
    marker_scores.tsv         epithelial and immune ssGSEA scores
    markers_used.tsv          the panel, with a found-in-matrix flag
  03_erv/
    align/                    permissive bowtie2 BAMs (large)
    telescope/<sample>/        per-sample Telescope reports
    counts/                   erv_counts_{locus,family,stats}.tsv
    diagnostics/              sorted BAMs, written on demand
    multiqc/
  04_viral/
    align/                    viral-only BAMs (small) and bowtie2 logs
    counts_per_sample/        per-sample metric tables
    counts/                   viral_per_{accession,taxid}.tsv,
                              viral_counts_{taxid,family}.tsv,
                              viral_denominators.tsv, viral_calls.tsv
    multiqc/
  logs/                       <jobname>_<arrayid>_<taskid>.{out,err}
```

The handoff to downstream analysis is four files:

| File | Contents |
| --- | --- |
| `01_host/counts/gene_counts.txt` | host gene counts |
| `02_expression/size_factors.tsv` | normalisation for both arms |
| `03_erv/counts/erv_counts_locus.tsv` and `_family.tsv` | raw ERV counts |
| `04_viral/counts/viral_calls.tsv` | every viral taxon-sample row with its tier |

All counts are raw. Normalisation happens in R, so a change of model does not
mean rerunning the pipeline.

---

## Design decisions worth knowing

These are choices that would be reasonable to make differently, and that
affect interpretation.

**Duplicates are marked, never removed.** Removing them would destroy the
library-complexity measurement module 01 exists to make. The ERV and viral
alignments are not deduplicated at all, which matters for samples with high
duplication: a few original fragments amplified heavily will count as many.

**The ERV and viral arms each use a whole-genome index.** Mapping to a small
ERV-only or virus-only reference forces every partially similar read onto it.
The pilot showed this empirically for the viral arm: up to 71% of the reads
STAR could not place map back to human under local alignment, and without
GRCh38 in the same index they would have been assigned to viral genomes
instead.

**Module 03 is deliberately permissive, module 04 deliberately is not.**
Telescope needs ambiguous alignments to redistribute, so 03a uses `-k 100` at
about 95% identity. Module 04 wants one best hit per read so that MAPQ stays
meaningful, so 04a passes no `-k`.

**ERV counts are normalised by cellular genes, never by ERV counts.**
Self-normalisation would divide away a global ERV induction, which is one of
the effects the study is looking for.

**Size factors are estimated within each body site.** The median-of-ratios
estimator assumes most genes are unchanged between the samples compared, and
across these sites that fails: they differ in epithelial type, bacterial load
and recoverable human transcriptome. Pooling would absorb site differences
into the factors. The cost is that each factor rests on few libraries, so this
is a judgement that should be revisited with cohort sizes in hand.

**ssGSEA runs with `normalize = FALSE`.** The default rescales each score by
the range across all samples in the run, which would make a sample's value
depend on which others were run with it and would shift when cohort samples
are added. Raw scores are stable across runs.

**No positivity rule lives in the pipeline.** Counts, breadth, covered windows
and both RPM denominators are written out; the rule is set in the report, with
thresholds in an editable chunk at the top.

---

## Known limitations

**LINE-1 is roughly half the ERV annotation.** `retro.hg38.v1` combines
HERV_rmsk with L1Base, and 13,438 of its 27,748 loci are LINE-1. They are kept
in the annotation deliberately: removing them would push their reads onto HERV
loci during the EM step. They are reported separately throughout, and a
combined "ERV total" would be misleading. LINE-1 also resolves to very few
families, so family-level aggregation is not meaningful on that side.

**`unique_count` is confounded with family copy number.** A fragment counts as
unique only if all its alignments fall inside one annotated locus. In a family
with over a thousand near-identical copies, almost nothing qualifies however
real the expression, so a low unique fraction says the family is repetitive,
not that the signal is unreliable. It is reported next to family size for that
reason and must not be used as a quality filter on its own.

**Family totals are better supported than individual ERV loci.** Where
ambiguity is high, the family-level matrix is the defensible unit and
locus-level results are exploratory. The pilot's analysable locus sets (loci
detected in every sample of a site) were 82 to 254 HERV loci per site, not
14,000.

**One viral genome per taxon means divergent strains may be missed**, and
conversely that a real infection spreads its reads across many near-identical
references of the same family. Genus and family totals are the firm numbers.

**Viral breadth is not comparable across genome sizes.** The same covered-window
count is most of a 7 kb enterovirus and a few percent of a 472 kb
herpesvirus, which is why `windows_fraction` is reported and used in the rule.

**DNA viruses are detected through transcription, not presence.** A
herpesvirus hit is lytic transcription, which is a different claim from the
virus being present, and coverage concentrates on transcribed genes so breadth
stays low however real the signal. These land in the report's suggestive tier
and need a coverage profile before being called either way.

**Anelloviruses are not a reliable positive control here.** They are DNA
viruses detected via transcription, and among the most diverse families known
against a one-genome-per-taxon reference. In the pilot their reads were
present but stacked on short conserved stretches, which neither confirms nor
refutes them.

**The Picard ribosomal metric is meaningless on this genome.** The 45S rDNA
arrays are poorly represented in GRCh38, so `PCT_RIBOSOMAL_BASES` reads near
zero regardless of true rRNA content. SortMeRNA carries that metric.

**The SortMeRNA eukaryotic database is a proxy for human rRNA.** Fungal rRNA
matches it too, so the fraction may be inflated at sites with fungal load.

**Excluding bacteriophage removed a third of the viral accessions**, not the
2% the row count suggests, because those rows carry many accessions each. 46
further accessions in the table were absent from the VirusHostDB FASTA; they
are listed in `missing.tsv` in the viral set directory.

**Telescope `final_count` is an expected value under a model**, not a tally of
unambiguous reads. For any locus that matters, `unique_count` in the
per-sample report gives the unambiguous subset.

---

## Troubleshooting

**Reading the logs.** `$LOG_DIR/<jobname>_<arrayid>_<taskid>.err` for a failed
array task. The task ID is the sample's row in the sample sheet, excluding the
header: task 3 is the third sample.

**Resubmitting one sample.** All array steps skip samples whose outputs exist,
so rerunning the `run` script is usually enough. To force one task:

```bash
sbatch --array=3 --partition=<queue> ... modules/01a_align.sh
```

**A sample sheet that looks fine but produces odd groupings.** Check
capitalisation in `body_site`. The pilot had both `Nasal` and `nasal`, which
split three samples into six groups and gave each singleton a size factor of
exactly 1.

**MarkDuplicates fails with a NullPointerException about read groups.** The
BAM predates the read group that 01a now adds. Rerun 01a; it detects BAMs
without `@RG` and realigns them.

**A container-based module cannot see a file.** Several tools load as
Singularity containers, which only see host paths that are bound into them. If
a tool reports a missing file that exists, check whether `REF_DIR` or
`OUT_DIR` is on a path the site's Singularity configuration binds.

**Investigating one ERV locus or one viral taxon.** Two diagnostic scripts
report coverage across the feature, MAPQ distribution and, for the viral one,
how many separated blocks the coverage falls into:

```bash
bash scripts/erv_locus_diagnostic.sh CMEH_Vag HERVH_20p11.23b
bash scripts/viral_taxon_diagnostic.sh CMWZ_Nas "Human betaherpesvirus 6B"
```

Both need the arm's BAMs. The ERV one sorts a large BAM on first use for a
given sample and keeps it in `03_erv/diagnostics/`; run it as a batch job with
8 cores. The viral one needs `VIRAL_KEEP_BAM=true` to have been set.

**Running out of disk.** The module 03 BAMs are the largest thing the pipeline
writes, because `-k 100` reports up to 100 alignments per read. Both check
scripts print the disk used by their alignments.
`ERV_DELETE_BAM_AFTER_ASSIGN=true` deletes each ERV BAM once its Telescope run
succeeds, at the cost of needing a realignment to rerun 03b.

---

## Extending the pipeline

**Adding a reference, database or parameter option:** add the value to the
relevant conditional in the module script, update the error message that lists
valid values, update the comment in `config.sh`, update this README everywhere
the valid values appear, and test on one sample before releasing.

**Changing a tool version:** an index is only valid for the tool version that
built it. The STAR and bowtie2 index directory names carry their tool version
for this reason, so a new version gets a new index rather than silently
reusing an incompatible one. Update both `config.sh` and
`references/ref.config.sh`, which hold the versions independently.

**A note on featureCounts.** In subread 2.0.1, `-p` alone counts fragments.
From 2.0.2 onwards `--countReadPairs` is also required, or reads are counted
instead. Upgrading subread without changing `03c`'s call would silently change
every count.

---

## Downstream of the pipeline

Contamination and batch control, and the statistics, are not part of VERVE.
They depend on the study design, which is why they are left out rather than
built against assumptions. What they will need:

- A metadata file supplied separately: sample, participant, group, age, sex,
  batch, plate and ideally well position. It cannot be auto-generated.
- Blank controls per batch. The number available determines whether only a
  ratio-to-blank rule is possible or a prevalence-based method as well.
- Cases and controls randomised across plates. If they are not, batch and
  group may be confounded in a way no downstream correction fixes cleanly.
- A gene set for the interferon confounder check.

