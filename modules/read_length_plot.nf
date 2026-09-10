process READ_LENGTH_DISTRIBUTION {
    tag "$meta.sample"

    container 'community.wave.seqera.io/library/samtools_pip_matplotlib_pandas:83512509914d9566'

    input:
    tuple val(meta), path(dedup_bam), path(dedup_bai)

    output:
    tuple val(meta), path("${meta.sample}.read_length_distribution.*"), emit: plot

    script:
    """
    samtools view -F 0x904 ${dedup_bam} > ${meta.sample}.filtered.sam

    python3 <<EOF
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from collections import Counter

    counts = Counter()
    with open("${meta.sample}.filtered.sam") as f:
        for line in f:
            fields = line.rstrip("\\n").split("\\t")
            flag = int(fields[1])
            seq = fields[9]
            length = len(seq)
            strand = "reverse" if (flag & 16) else "forward"
            counts[(length, strand)] += 1

    df = pd.DataFrame(
        [(length, strand, n) for (length, strand), n in counts.items()],
        columns=["length", "strand", "count"]
    ).sort_values(["strand", "length"])
    df.to_csv("${meta.sample}.read_length_distribution.tsv", sep="\\t", index=False)

    pivot = df.pivot_table(index="length", columns="strand", values="count", fill_value=0)
    for col in ["forward", "reverse"]:
        if col not in pivot.columns:
            pivot[col] = 0
    pivot = pivot.sort_index()

    lengths = pivot.index.to_numpy()
    bar_width = 0.4

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(lengths - bar_width / 2, pivot["forward"], width=bar_width,
        label="forward", color="#ff3637")
    ax.bar(lengths + bar_width / 2, pivot["reverse"], width=bar_width,
        label="reverse", color="#3bafff")

    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Read count")
    ax.set_title("Host Mapped Read Length Distribution: ${meta.sample}")
    ax.legend(title="Strand")
    plt.tight_layout()
    plt.savefig("${meta.sample}.read_length_distribution.png", dpi=150)
    EOF

    rm ${meta.sample}.filtered.sam
    """
}