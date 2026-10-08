process FEATURE_COUNTS {
    tag "$meta.sample"
    label 'demand_1'
    container 'community.wave.seqera.io/library/subread:2.1.1--ab75f84c2fa06990'

    input: 
    tuple val(meta), path(dedup_bam), path(dedup_bai)
    path gtf
    val strand
    val multi

    output:
    tuple val(meta), path("${meta.sample}_exon_fc.txt"), emit: exon_counts
    tuple val(meta), path("${meta.sample}_exon_fc.txt.summary"), emit: exon_summary
    tuple val(meta), path("${meta.sample}_gene_fc.txt"), emit: gene_counts
    tuple val(meta), path("${meta.sample}_gene_fc.txt.summary"), emit: gene_summary

    script:
    def multi_args = multi ? '-M --fraction' : ''

    """
    featureCounts \\
        -t exon \\
        -g gene_id \\
        -s ${strand} \\
        ${multi_args} \\
        -a ${gtf} \\
        -o ${meta.sample}_exon_fc.txt \\
        ${dedup_bam}

    featureCounts \\
        -t gene \\
        -g gene_id \\
        -s ${strand} \\
        ${multi_args} \\
        -a ${gtf} \\
        -o ${meta.sample}_gene_fc.txt \\
        ${dedup_bam}
    """
}


process FEATURECOUNTS_SUMMARY {
    label 'demand_1'
    container 'community.wave.seqera.io/library/pip_pandas:ab7edebb26d10817'

    input:
    path fc_files   
    path gtf
    val level

    output:
    path "featurecounts_summary_${level}.csv", emit: csv

    script:
    """
    #!/usr/bin/env python3
    import glob
    import re
    import sys

    import pandas as pd

    GTF = "${gtf}"
    TAB = chr(9)

    # 1. Build gene_id -> label mapping from GTF
    gene_label = {}
    with open(GTF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split(TAB)
            if len(cols) < 9 or cols[2] != "exon":
                continue

            chrom = cols[0]
            gid = re.search(r'gene_id "([^"]+)"', line)
            gt  = re.search(r'gene_type "([^"]+)"', line)

            if gid:
                if gt:
                    label = gt.group(1)
                    # Reclassify chrM protein-coding genes
                    if chrom == "chrM" and label == "protein_coding":
                        label = "Mt_protein_coding"
                else:
                    label = "unknown"
                gene_label[gid.group(1)] = label

    print(f"{len(gene_label)} unique gene_ids mapped")

    # 2. Process each fc.txt file
    fc_paths = sorted(p for p in glob.glob("*.txt") if not p.endswith(".summary"))
    if not fc_paths:
        sys.exit("No featureCounts .txt files found")

    records = []
    for fc_path in fc_paths:
        sample = (
            fc_path
            .replace("_gene_fc.txt", "")
            .replace("_exon_fc.txt", "")
        )
        summary_path = fc_path + ".summary"

        # Parse .summary for alignment stats
        total_alignments = None
        assigned = None
        try:
            status_counts = {}
            with open(summary_path) as f:
                next(f)  # skip header line
                for line in f:
                    parts = line.strip().split(TAB)
                    if len(parts) >= 2:
                        status_counts[parts[0]] = int(parts[1])
            assigned = status_counts.get("Assigned", 0)
            total_alignments = sum(status_counts.values())
        except FileNotFoundError:
            pass

        # Parse counts
        fc = pd.read_csv(fc_path, sep=TAB, comment="#")
        fc.columns = [c.strip() for c in fc.columns]
        fc = fc[[fc.columns[0], fc.columns[-1]]].copy()
        fc.columns = ["gene_id", "count"]
        fc["label"] = fc["gene_id"].map(gene_label).fillna("unknown")
        label_counts = fc.groupby("label")["count"].sum()

        row = {
            "sample":              sample,
            "total_alignments":    total_alignments,
            "assigned_alignments": assigned,
            "pct_assigned":        round(100 * assigned / total_alignments, 2)
                                   if total_alignments and assigned else None,
        }
        row.update(label_counts.to_dict())
        records.append(row)

    # 3. Build final table
    df = pd.DataFrame(records)
    meta_cols  = ["sample", "total_alignments", "assigned_alignments", "pct_assigned"]
    label_cols = sorted(c for c in df.columns if c not in meta_cols)
    df = df[meta_cols + label_cols].fillna(0)
    for col in label_cols:
        df[col] = df[col].astype(int)

    df.to_csv("featurecounts_summary_${level}.csv", index=False)
    print(df[meta_cols + label_cols[:5]].to_string(index=False))
    """
}