#!/usr/bin/env python3
"""
Merge per-sample viral metrics into matrices, and aggregate to taxon level.

Accessions are aggregated to virus taxid first, which is what joins the
segments of a segmented virus into one number, then to genus and family.
Breadth after aggregation is covered bases over the summed length of all the
taxon's sequences, so a rotavirus is judged on its whole genome rather than
per segment.

windows_fraction is the covered windows over the total the reference divides
into. It is reported because the raw window count is not comparable across
genome sizes: 15 covered windows is most of a 7 kb enterovirus and 3% of a
472 kb herpesvirus.

RPM is reported against two denominators, because neither is right on its own:
  rpm_total   per million reads surviving trimming (STAR input)
  rpm_nonhost per million reads that STAR could not place, which is what was
              actually offered to this module
Samples differ enormously in rRNA content and in host fraction, so the two
can disagree by an order of magnitude. Both are written; the methods should
say which one a given figure uses.

Counts are READS, not fragments. Denominators are converted to reads by
doubling the pair counts that STAR and bowtie2 report.

Usage:
  python3 merge_viral.py --sample-sheet samples.tsv --out-dir <OUT_DIR> \\
      --taxa taxa.tsv --out-prefix <OUT_DIR>/04_viral/counts/viral
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--sample-sheet", required=True)
    parser.add_argument("--out-dir", required=True,
                        help="pipeline OUT_DIR, for finding the per-sample files")
    parser.add_argument("--taxa", required=True, help="taxa.tsv from ref.05")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--window-size", type=int, default=1000,
                        help="window size 04a used, for windows_fraction "
                             "(default 1000)")
    return parser.parse_args()


def read_samples(path):
    samples = []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample_id = (row.get("sample_id") or "").strip()
            if sample_id:
                samples.append(sample_id)
    return samples


def read_taxa(path):
    taxa = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            taxa[row["accession"]] = row
    return taxa


def read_metrics(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def first_number(path, pattern):
    """First capture group of pattern in a file, as a float, or None."""
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            found = re.search(pattern, line)
            if found:
                return float(found.group(1).replace(",", ""))
    return None


def denominators(out_dir, sample_id):
    """Read pairs after trimming (STAR input) and pairs offered to bowtie2."""
    star_log = os.path.join(out_dir, "01_host", "align",
                            f"{sample_id}.Log.final.out")
    total_pairs = first_number(star_log,
                               r"Number of input reads\s*\|\s*([0-9]+)")

    bowtie_log = os.path.join(out_dir, "04_viral", "align",
                              f"{sample_id}.bowtie2.log")
    nonhost_pairs = first_number(bowtie_log, r"^\s*([0-9]+) reads; of these:")

    return total_pairs, nonhost_pairs


def windows_total(length_bp, window_size):
    """Number of windows a reference of this length is divided into."""
    if not length_bp or window_size <= 0:
        return 0
    return max(1, -(-int(length_bp) // window_size))


def rpm(reads, pairs):
    if pairs is None or pairs <= 0:
        return ""
    return round(reads / (2.0 * pairs) * 1e6, 4)


def write_tsv(path, columns, rows):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_matrix(path, row_label, row_names, samples, values):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([row_label] + samples)
        for name in row_names:
            row = [name]
            for sample in samples:
                value = values.get(sample, {}).get(name, 0)
                row.append(int(value) if float(value).is_integer() else value)
            writer.writerow(row)


def main():
    args = parse_args()
    samples = read_samples(args.sample_sheet)
    taxa = read_taxa(args.taxa)
    print(f"Accessions in taxon table: {len(taxa)}", file=sys.stderr)

    # Total sequence length per taxid, for breadth after aggregation
    taxid_length = defaultdict(float)
    taxid_meta = {}
    for accession, row in taxa.items():
        taxid = row["virus_taxid"]
        try:
            taxid_length[taxid] += float(row.get("genome_length") or 0)
        except ValueError:
            pass
        taxid_meta.setdefault(taxid, {
            "virus_taxid": taxid,
            "virus_name": row.get("virus_name", ""),
            "genus": row.get("genus", ""),
            "family": row.get("family", ""),
            "category": row.get("category", ""),
            "flags": row.get("flags", ""),
            "segments_total": row.get("segments_total", ""),
        })

    per_accession = []
    per_taxid = []
    denom_rows = []
    counts_taxid = defaultdict(dict)
    counts_family = defaultdict(dict)
    unknown_accessions = set()
    present = []
    missing_metrics = []

    for sample_id in samples:
        metrics_path = os.path.join(args.out_dir, "04_viral",
                                    "counts_per_sample",
                                    f"{sample_id}.viral_metrics.tsv")
        if not os.path.exists(metrics_path):
            missing_metrics.append(sample_id)
            continue
        present.append(sample_id)

        total_pairs, nonhost_pairs = denominators(args.out_dir, sample_id)
        denom_rows.append({
            "sample_id": sample_id,
            "trimmed_pairs": int(total_pairs) if total_pairs else "",
            "nonhost_pairs": int(nonhost_pairs) if nonhost_pairs else "",
            "nonhost_fraction": (round(nonhost_pairs / total_pairs, 4)
                                 if total_pairs and nonhost_pairs else ""),
        })

        rows = read_metrics(metrics_path)

        # Per accession, with taxon annotation attached
        by_taxid = defaultdict(lambda: {"reads": 0.0, "reads_hiq": 0.0,
                                        "covered_bases": 0.0, "windows": 0.0,
                                        "accessions_hit": 0})
        for row in rows:
            accession = row["accession"]
            taxon = taxa.get(accession)
            if taxon is None:
                unknown_accessions.add(accession)
                continue

            reads = float(row["reads"])
            reads_hiq = float(row["reads_hiq"])

            per_accession.append({
                "sample_id": sample_id,
                "accession": accession,
                "virus_taxid": taxon["virus_taxid"],
                "virus_name": taxon.get("virus_name", ""),
                "genus": taxon.get("genus", ""),
                "family": taxon.get("family", ""),
                "segment_n": taxon.get("segment_n", ""),
                "ref_length": row["ref_length"],
                "reads": int(reads),
                "reads_hiq": int(reads_hiq),
                "covered_bases": int(float(row["covered_bases"])),
                "breadth": row["breadth"],
                "mean_depth": row["mean_depth"],
                "mean_mapq": row["mean_mapq"],
                "windows_covered": row["windows_covered"],
                "windows_total": windows_total(float(row["ref_length"]),
                                               args.window_size),
                "windows_fraction": round(
                    float(row["windows_covered"]) /
                    windows_total(float(row["ref_length"]), args.window_size), 4)
                    if windows_total(float(row["ref_length"]),
                                     args.window_size) > 0 else "",
                "flags": taxon.get("flags", ""),
            })

            entry = by_taxid[taxon["virus_taxid"]]
            entry["reads"] += reads
            entry["reads_hiq"] += reads_hiq
            entry["covered_bases"] += float(row["covered_bases"])
            entry["windows"] += float(row["windows_covered"])
            entry["accessions_hit"] += 1

        # Per taxid, segments summed
        for taxid, entry in by_taxid.items():
            meta = taxid_meta.get(taxid, {})
            total_length = taxid_length.get(taxid, 0.0)
            breadth = (entry["covered_bases"] / total_length
                       if total_length > 0 else 0.0)
            per_taxid.append({
                "sample_id": sample_id,
                "virus_taxid": taxid,
                "virus_name": meta.get("virus_name", ""),
                "genus": meta.get("genus", ""),
                "family": meta.get("family", ""),
                "category": meta.get("category", ""),
                "flags": meta.get("flags", ""),
                "segments_total": meta.get("segments_total", ""),
                "accessions_hit": entry["accessions_hit"],
                "total_ref_length": int(total_length),
                "reads": int(entry["reads"]),
                "reads_hiq": int(entry["reads_hiq"]),
                "covered_bases": int(entry["covered_bases"]),
                "breadth": round(breadth, 6),
                "windows_covered": int(entry["windows"]),
                "windows_total": windows_total(total_length, args.window_size),
                "windows_fraction": (
                    round(entry["windows"] /
                          windows_total(total_length, args.window_size), 4)
                    if windows_total(total_length, args.window_size) > 0 else ""),
                "rpm_total": rpm(entry["reads"], total_pairs),
                "rpm_nonhost": rpm(entry["reads"], nonhost_pairs),
            })
            counts_taxid[sample_id][taxid] = entry["reads"]
            family = meta.get("family") or "(no family)"
            counts_family[sample_id][family] = \
                counts_family[sample_id].get(family, 0) + entry["reads"]

        print(f"{sample_id}: {len(rows)} accessions with reads, "
              f"{len(by_taxid)} taxa", file=sys.stderr)

    if missing_metrics:
        print("WARNING: no metrics file for: " + ", ".join(missing_metrics),
              file=sys.stderr)
    if unknown_accessions:
        print(f"WARNING: {len(unknown_accessions)} accessions had reads but no "
              f"row in the taxon table, and were dropped: "
              f"{', '.join(sorted(unknown_accessions)[:5])}", file=sys.stderr)
    if not present:
        raise SystemExit("ERROR: no per-sample metrics found, nothing to merge")

    # ── Write ────────────────────────────────────────────────────────────────
    write_tsv(f"{args.out_prefix}_per_accession.tsv",
              ["sample_id", "accession", "virus_taxid", "virus_name", "genus",
               "family", "segment_n", "ref_length", "reads", "reads_hiq",
               "covered_bases", "breadth", "mean_depth", "mean_mapq",
               "windows_covered", "windows_total", "windows_fraction", "flags"],
              per_accession)

    write_tsv(f"{args.out_prefix}_per_taxid.tsv",
              ["sample_id", "virus_taxid", "virus_name", "genus", "family",
               "category", "flags", "segments_total", "accessions_hit",
               "total_ref_length", "reads", "reads_hiq", "covered_bases",
               "breadth", "windows_covered", "windows_total",
               "windows_fraction", "rpm_total", "rpm_nonhost"],
              per_taxid)

    write_tsv(f"{args.out_prefix}_denominators.tsv",
              ["sample_id", "trimmed_pairs", "nonhost_pairs",
               "nonhost_fraction"],
              denom_rows)

    taxid_rows = sorted({t for s in counts_taxid.values() for t in s})
    write_matrix(f"{args.out_prefix}_counts_taxid.tsv", "virus_taxid",
                 taxid_rows, present, counts_taxid)

    family_rows = sorted({f for s in counts_family.values() for f in s})
    write_matrix(f"{args.out_prefix}_counts_family.tsv", "family",
                 family_rows, present, counts_family)

    print("", file=sys.stderr)
    print(f"Samples merged : {len(present)}", file=sys.stderr)
    print(f"Taxa with reads: {len(taxid_rows)}", file=sys.stderr)
    print(f"Families       : {len(family_rows)}", file=sys.stderr)
    print(f"Written with prefix {args.out_prefix}", file=sys.stderr)


if __name__ == "__main__":
    main()
