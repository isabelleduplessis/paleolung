#!/usr/bin/env python3
"""
gene_strand_quant.py

For each gene in a GFF3 file, count reads (from a BAM file) that overlap the
gene's coordinates, and classify each read as mapping to the same strand as
the gene or the opposite strand (forward-stranded library: read strand ==
RNA sense strand).

Usage:
    python gene_strand_quant.py -b reads.bam -g genes.gff3 -o out.tsv \
        [--feature-type gene]
"""

import argparse
import re
import sys

import pysam


def parse_gff3_attributes(attr_string):
    """Parse the GFF3 column 9 attribute string into a dict."""
    attrs = {}
    for field in attr_string.strip().split(";"):
        field = field.strip()
        if not field or "=" not in field:
            continue
        key, value = field.split("=", 1)
        attrs[key.strip()] = value.strip()
    return attrs


def get_locus_id(attrs):
    """Prefer locus_tag, then ID, then Name, as the identifier for the gene."""
    for key in ("locus_tag", "ID", "Name"):
        if key in attrs:
            return attrs[key]
    return "unknown"


def load_genes(gff3_path, feature_type):
    """Yield (locus_id, seqid, start0, end, strand) for each matching feature.

    GFF3 coordinates are 1-based inclusive; we convert start to 0-based for
    use with pysam's fetch (which takes 0-based, half-open coordinates).
    """
    genes = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr_str = fields[:9]
            if ftype != feature_type:
                continue
            if strand not in ("+", "-"):
                sys.stderr.write(
                    f"Warning: skipping feature with unknown strand '{strand}': {line.strip()}\n"
                )
                continue
            attrs = parse_gff3_attributes(attr_str)
            locus_id = get_locus_id(attrs)
            start0 = int(start) - 1  # GFF3 1-based -> 0-based
            end = int(end)           # end stays as-is (0-based half-open == 1-based inclusive end)
            genes.append((locus_id, seqid, start0, end, strand))
    return genes


def classify_reads(bam, seqid, start0, end, gene_strand):
    """Fetch all reads overlapping [start0, end) on seqid and classify strand.

    Returns (n_overlap, n_same, n_opposite).
    Forward-stranded library: read maps to '+' -> read is on '+' strand of
    genome (is_reverse == False); read strand is compared directly to the
    gene's strand.
    """
    n_overlap = 0
    n_same = 0
    n_opposite = 0

    for read in bam.fetch(seqid, start0, end):
        n_overlap += 1
        read_strand = "-" if read.is_reverse else "+"
        if read_strand == gene_strand:
            n_same += 1
        else:
            n_opposite += 1

    return n_overlap, n_same, n_opposite


def main():
    parser = argparse.ArgumentParser(
        description="Quantify per-gene same-strand vs opposite-strand read overlap."
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file (indexed).")
    parser.add_argument("-g", "--gff3", required=True, help="Input GFF3 file.")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path.")
    parser.add_argument(
        "--feature-type",
        default="gene",
        help="GFF3 column 3 feature type to use as 'gene' (default: gene).",
    )
    args = parser.parse_args()

    genes = load_genes(args.gff3, args.feature_type)
    if not genes:
        sys.stderr.write(
            f"No features of type '{args.feature_type}' found in {args.gff3}\n"
        )
        sys.exit(1)

    # Ensure BAM is indexed; build an index if missing.
    bam_path = args.bam
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        if not bam.has_index():
            raise ValueError("no index")
    except (ValueError, OSError):
        sys.stderr.write(f"Indexing {bam_path}...\n")
        pysam.index(bam_path)
        bam = pysam.AlignmentFile(bam_path, "rb")

    bam_contigs = set(bam.references)

    with open(args.output, "w") as out:
        out.write(
            "locus_id\tseqid\tstart\tend\tstrand\tn_overlapping_reads\t"
            "n_same_strand\tn_opposite_strand\tpercent_same_strand\t"
            "percent_opposite_strand\n"
        )
        for locus_id, seqid, start0, end, strand in genes:
            if seqid not in bam_contigs:
                sys.stderr.write(
                    f"Warning: {seqid} (gene {locus_id}) not found in BAM contigs, skipping.\n"
                )
                continue

            n_overlap, n_same, n_opposite = classify_reads(bam, seqid, start0, end, strand)

            if n_overlap > 0:
                pct_same = 100.0 * n_same / n_overlap
                pct_opposite = 100.0 * n_opposite / n_overlap
            else:
                pct_same = 0.0
                pct_opposite = 0.0

            out.write(
                f"{locus_id}\t{seqid}\t{start0 + 1}\t{end}\t{strand}\t"
                f"{n_overlap}\t{n_same}\t{n_opposite}\t"
                f"{pct_same:.2f}\t{pct_opposite:.2f}\n"
            )

    bam.close()
    sys.stderr.write(f"Done. Wrote {len(genes)} gene rows to {args.output}\n")


if __name__ == "__main__":
    main()
