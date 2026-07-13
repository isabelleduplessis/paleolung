#!/usr/bin/env python3
"""
assign_features.py

Assigns each mapped read to exactly one genomic feature category based on
a priority-ordered list of BED files. Reads are assigned to the first BED
they overlap (highest priority), and each subsequent category excludes all
higher-priority reads. Unassigned reads are counted as "intergenic".

Outputs:
  - A TSV summary of counts and percentages per category (written to outdir)
  - One subset BAM per BED file (plus an intergenic BAM)

Requirements:
    bedtools, samtools on PATH
    pip install pysam

Usage:
    python assign_features.py \\
        -b reads.bam \\
        -f ribo.bed ss_exons.bed ds_exons.bed ss_introns.bed ds_introns.bed te.bed \\
        -n ribo ss_exon ds_exon ss_intron ds_intron te \\
        -s {0,1,2} \\
        [-o|--outdir results/] \\
        [--tsv output.tsv] \\
        [--sample SAMPLE_NAME]

Notes:
  - BED files and names must be given in priority order (highest first).
  - outdir will be created automatically if it does not exist.
  - Strandedness (-s):
      0 = unstranded (no -s/-S flag passed to bedtools intersect)
      1 = stranded, same strand (-s)
      2 = stranded, opposite strand (-S)
"""

import argparse
import os
import subprocess
import sys
import tempfile


# ── Helpers ──────────────────────────────────────────────────────────────────

def run(cmd, stdin=None, capture_stdout=False):
    """Run a shell command, optionally piping stdin and capturing stdout."""
    result = subprocess.run(
        cmd, shell=True,
        stdin=stdin,
        stdout=subprocess.PIPE if capture_stdout else None,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        print(f"Error running:\n  {cmd}", file=sys.stderr)
        print(result.stderr.decode(), file=sys.stderr)
        sys.exit(1)
    return result


def count_mapped(bam):
    """Count total mapped reads in BAM."""
    r = run(f"samtools view -c -F 4 {bam}", capture_stdout=True)
    return int(r.stdout.strip())


def samtools_index(bam):
    run(f"samtools index {bam}")


def strand_flag(strandedness):
    """Return the bedtools intersect strand flag for the given strandedness code."""
    if strandedness == 1:
        return "-s"
    elif strandedness == 2:
        return "-S"
    else:  # 0
        return ""


# ── Core logic ────────────────────────────────────────────────────────────────

def intersect_priority(bam, bed_files, names, outdir, sample, strandedness):
    """
    For each feature in priority order:
      1. Intersect BAM with this BED (-u: any overlap, report read once)
      2. Subtract all higher-priority BEDs (-v: keep only non-overlapping)
      3. Count and write subset BAM

    Returns list of (name, count, bam_path) tuples.
    """
    results = []
    higher_priority_beds = []
    sflag = strand_flag(strandedness)

    for bed, name in zip(bed_files, names):
        print(f"  Processing: {name} ({bed})", file=sys.stderr)

        # First: intersect with this feature's BED (-u = report read once if any overlap)
        cmd = f"bedtools intersect -a {bam} -b {bed} -u {sflag}".strip()

        # Then subtract each higher-priority BED in sequence
        for prev_bed in higher_priority_beds:
            cmd += f" | bedtools intersect -a stdin -b {prev_bed} -v {sflag}".rstrip()

        # Write to output BAM
        out_bam = os.path.join(outdir, f"{sample}.{name}.bam")

        # Pipe through samtools to write BAM (preserving header)
        full_cmd = f"{cmd} | samtools view -b -o {out_bam}"
        run(full_cmd)
        samtools_index(out_bam)

        # Count reads in the output BAM
        count_r = run(f"samtools view -c -F 4 {out_bam}", capture_stdout=True)
        count = int(count_r.stdout.strip())

        results.append((name, count, out_bam))
        higher_priority_beds.append(bed)

    return results


def write_intergenic_bam(bam, bed_files, outdir, sample, strandedness):
    """Extract reads that don't overlap ANY BED file."""
    print("  Processing: intergenic", file=sys.stderr)
    sflag = strand_flag(strandedness)
    cmd = f"samtools view -b -F 4 {bam}"
    for bed in bed_files:
        cmd += f" | bedtools intersect -a stdin -b {bed} -v {sflag}".rstrip()
    out_bam = os.path.join(outdir, f"{sample}.intergenic.bam")
    full_cmd = f"{cmd} | samtools view -b -o {out_bam}"
    run(full_cmd)
    samtools_index(out_bam)
    count_r = run(f"samtools view -c -F 4 {out_bam}", capture_stdout=True)
    count = int(count_r.stdout.strip())
    return count, out_bam


def write_tsv(tsv_path, sample, total, results, intergenic_count, intergenic_bam):
    """Write summary TSV with counts and percentages."""
    all_results = results + [("intergenic", intergenic_count, intergenic_bam)]

    header_counts = "\t".join(name for name, _, _ in all_results)
    header_pcts   = "\t".join(f"{name}_pct" for name, _, _ in all_results)
    header_bams   = "\t".join(f"{name}_bam" for name, _, _ in all_results)

    counts = "\t".join(str(c) for _, c, _ in all_results)
    pcts   = "\t".join(
        f"{c/total*100:.2f}" if total > 0 else "0.00"
        for _, c, _ in all_results
    )
    bams   = "\t".join(b for _, _, b in all_results)

    write_header = not os.path.exists(tsv_path)
    with open(tsv_path, "a") as f:
        if write_header:
            f.write(f"sample\ttotal\t{header_counts}\t{header_pcts}\t{header_bams}\n")
        f.write(f"{sample}\t{total}\t{counts}\t{pcts}\t{bams}\n")

    print(f"\nSummary TSV → {tsv_path}")


def print_summary(sample, total, results, intergenic_count, intergenic_bam):
    """Pretty-print results to stderr."""
    all_results = results + [("intergenic", intergenic_count, intergenic_bam)]
    width = max(len(n) for n, _, _ in all_results) + 2
    print(f"\n{'─'*50}", file=sys.stderr)
    print(f"  Sample : {sample}", file=sys.stderr)
    print(f"  Total  : {total:,}", file=sys.stderr)
    print(f"{'─'*50}", file=sys.stderr)
    assigned = 0
    for name, count, bam_path in all_results:
        pct = count / total * 100 if total > 0 else 0
        print(f"  {name:<{width}} {count:>10,}  ({pct:5.1f}%)  → {bam_path}", file=sys.stderr)
        assigned += count
    print(f"{'─'*50}", file=sys.stderr)
    if assigned != total:
        print(f"  WARNING: assigned {assigned:,} but total is {total:,} — check for duplicates or flags", file=sys.stderr)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Priority-based read assignment to genomic features")
    ap.add_argument("-b", "--bam",    required=True,
                    help="Input BAM file (must be indexed)")
    ap.add_argument("-f", "--beds",   required=True, nargs="+",
                    help="BED files in priority order (highest first)")
    ap.add_argument("-n", "--names",  nargs="+", default=None,
                    help="Feature names matching --beds order. Defaults to BED filenames.")
    ap.add_argument("-s", "--strandedness", type=int, choices=[0, 1, 2], default=0,
                    help="Strandedness for bedtools intersect: "
                         "0 = unstranded (no flag), "
                         "1 = same strand (-s), "
                         "2 = opposite strand (-S). "
                         "Default: 0.")
    ap.add_argument("-o", "--outdir", default=".",
                    help="Output directory for subset BAMs and TSV (default: current dir). "
                         "Created automatically if it does not exist.")
    ap.add_argument("--tsv",          default=None,
                    help="Output TSV path (default: <outdir>/feature_assignment.tsv)")
    ap.add_argument("--sample",       default=None,
                    help="Sample name for TSV output. Defaults to BAM filename stem.")
    args = ap.parse_args()

    # ── Validate inputs ───────────────────────────────────────────────────────
    for bed in args.beds:
        if not os.path.exists(bed):
            sys.exit(f"Error: BED file not found: {bed}")
    if not os.path.exists(args.bam):
        sys.exit(f"Error: BAM file not found: {args.bam}")
    if not os.path.exists(args.bam + ".bai"):
        print(f"Index not found for {args.bam}, indexing...", file=sys.stderr)
        samtools_index(args.bam)

    names = args.names or [os.path.splitext(os.path.basename(b))[0] for b in args.beds]
    if len(names) != len(args.beds):
        sys.exit(f"Error: --names has {len(names)} entries but --beds has {len(args.beds)}")

    sample = args.sample or os.path.splitext(os.path.basename(args.bam))[0]
    os.makedirs(args.outdir, exist_ok=True)

    # Default TSV to outdir if not specified
    tsv_path = args.tsv if args.tsv else os.path.join(args.outdir, "feature_assignment.tsv")

    strand_label = {0: "unstranded", 1: "stranded (same, -s)", 2: "stranded (opposite, -S)"}
    print(f"Strandedness: {strand_label[args.strandedness]}", file=sys.stderr)

    # ── Run ───────────────────────────────────────────────────────────────────
    print(f"Counting total mapped reads...", file=sys.stderr)
    total = count_mapped(args.bam)
    print(f"  Total mapped: {total:,}", file=sys.stderr)

    print(f"\nAssigning reads by priority...", file=sys.stderr)
    results = intersect_priority(
        args.bam, args.beds, names, args.outdir, sample, args.strandedness)

    intergenic_count, intergenic_bam = write_intergenic_bam(
        args.bam, args.beds, args.outdir, sample, args.strandedness)

    print_summary(sample, total, results, intergenic_count, intergenic_bam)
    write_tsv(tsv_path, sample, total, results, intergenic_count, intergenic_bam)


if __name__ == "__main__":
    main()
