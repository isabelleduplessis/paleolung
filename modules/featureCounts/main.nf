process featureCountsGenes {
    input:
    tuple val(meta), path(bam)
    path(annotation)

    output:
    tuple val(meta), path("*_fc.txt"), emit: counts
    tuple val(meta), path("*_fc.txt.summary"), emit: summary

    script:
    def prefix = "${meta.sample_id}"
    """
    featureCounts \
        -t gene \
        -g gene_id \
        -s 1 \
        -a ${annotation} \
        -o ${prefix}_fc.txt \
        ${bam}
    """
}