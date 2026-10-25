#!/usr/bin/env python3
"""
Merge per-sample Telescope reports into count matrices.

Writes three files:
  <prefix>_locus.tsv    final_count per locus per sample
  <prefix>_family.tsv   the same counts summed per family
  <prefix>_stats.tsv    the run statistics from each report's comment line

Family is taken from the locus name: everything before the final underscore,
which is the cytoband (ERV316A3_1p36.33 gives ERV316A3, L1FLnI_1p36.13 gives
L1FLnI). The same rule covers the HERV and L1 components of retro.hg38.v1.

The __no_feature row is kept out of the matrices and reported in the stats
table instead.
"""

import argparse
import csv
import glob
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample-sheet", required=True)
    p.add_argument("--telescope-dir", required=True,
                   help="directory holding one subdirectory per sample")
    p.add_argument("--out-prefix", required=True,
                   help="output path prefix, e.g. .../erv_counts")
    p.add_argument("--count-column", default="final_count",
                   help="report column to use (default final_count, which is "
                        "the one for negative-binomial models such as DESeq2)")
    return p.parse_args()


def read_samples(path):
    samples = []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample_id = (row.get("sample_id") or "").strip()
            if sample_id:
                samples.append(sample_id)
    return samples


def find_report(telescope_dir, sample_id):
    """Telescope's report name differs between versions, so glob for it."""
    patterns = [
        os.path.join(telescope_dir, sample_id, f"{sample_id}-TE_counts.tsv"),
        os.path.join(telescope_dir, sample_id, f"{sample_id}-telescope_report.tsv"),
        os.path.join(telescope_dir, sample_id, f"{sample_id}-*counts.tsv"),
        os.path.join(telescope_dir, sample_id, f"{sample_id}-*report.tsv"),
    ]
    for pattern in patterns:
        hits = sorted(glob.glob(pattern))
        if hits:
            return hits[0]
    return None


def parse_report(path, count_column):
    """Return (counts by locus, run statistics from the comment line)."""
    counts = {}
    stats = {}
    with open(path, encoding="utf-8") as handle:
        first = handle.readline()
        if first.startswith("#"):
            # e.g. "## RunInfo\tkey:value\tkey:value"
            for field in first.lstrip("#").strip().split("\t"):
                if ":" in field:
                    key, value = field.split(":", 1)
                    stats[key.strip()] = value.strip()
            header_line = handle.readline()
        else:
            header_line = first

        header = header_line.rstrip("\n").split("\t")
        if count_column not in header:
            raise SystemExit(f"ERROR: column {count_column} not in {path}\n"
                             f"       columns are: {', '.join(header)}")
        locus_i = 0
        count_i = header.index(count_column)

        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= count_i:
                continue
            locus = fields[locus_i]
            try:
                value = float(fields[count_i])
            except ValueError:
                continue
            counts[locus] = value
    return counts, stats


def family_of(locus):
    return locus.rsplit("_", 1)[0] if "_" in locus else locus


def write_matrix(path, row_names, samples, values, row_label):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([row_label] + samples)
        for name in row_names:
            row = [name]
            for sample in samples:
                value = values[sample].get(name, 0.0)
                row.append(int(value) if float(value).is_integer() else value)
            writer.writerow(row)


def main():
    args = parse_args()
    samples = read_samples(args.sample_sheet)

    out_dir = os.path.dirname(os.path.abspath(args.out_prefix))
    os.makedirs(out_dir, exist_ok=True)

    counts_by_sample = {}
    stats_by_sample = {}
    no_feature = {}
    missing = []

    for sample_id in samples:
        report = find_report(args.telescope_dir, sample_id)
        if report is None:
            missing.append(sample_id)
            continue
        counts, stats = parse_report(report, args.count_column)
        no_feature[sample_id] = counts.pop("__no_feature", 0.0)
        counts_by_sample[sample_id] = counts
        stats_by_sample[sample_id] = stats
        print(f"{sample_id}: {len(counts)} loci from {os.path.basename(report)}",
              file=sys.stderr)

    if missing:
        print("WARNING: no Telescope report found for: " + ", ".join(missing),
              file=sys.stderr)
    if not counts_by_sample:
        raise SystemExit("ERROR: no reports found, nothing to merge")

    present = [s for s in samples if s in counts_by_sample]

    # Locus matrix
    loci = sorted({locus for counts in counts_by_sample.values() for locus in counts})
    write_matrix(f"{args.out_prefix}_locus.tsv", loci, present, counts_by_sample, "locus")

    # Family matrix
    family_counts = {}
    for sample_id, counts in counts_by_sample.items():
        per_family = {}
        for locus, value in counts.items():
            family = family_of(locus)
            per_family[family] = per_family.get(family, 0.0) + value
        family_counts[sample_id] = per_family
    families = sorted({f for per in family_counts.values() for f in per})
    write_matrix(f"{args.out_prefix}_family.tsv", families, present, family_counts, "family")

    # Run statistics
    stat_keys = sorted({k for stats in stats_by_sample.values() for k in stats})
    with open(f"{args.out_prefix}_stats.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "loci_detected", "no_feature"] + stat_keys)
        for sample_id in present:
            counts = counts_by_sample[sample_id]
            detected = sum(1 for v in counts.values() if v > 0)
            unassigned = no_feature.get(sample_id, 0.0)
            if float(unassigned).is_integer():
                unassigned = int(unassigned)
            writer.writerow([sample_id, detected, unassigned]
                            + [stats_by_sample[sample_id].get(k, "") for k in stat_keys])

    print(f"\nSamples merged : {len(present)}", file=sys.stderr)
    print(f"Loci           : {len(loci)}", file=sys.stderr)
    print(f"Families       : {len(families)}", file=sys.stderr)
    print(f"Written        : {args.out_prefix}_locus.tsv, "
          f"{args.out_prefix}_family.tsv, {args.out_prefix}_stats.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
