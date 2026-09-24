process FILTER_RIBO {
    tag "${meta.sample}"
    label 'demand_2'
    container 'community.wave.seqera.io/library/bedtools_samtools:650829d3c536249b'

    input:
    tuple val(meta), path(sorted_bam)
    path rrna_bed

    output:
    tuple val(meta), path("*.filtered.bam"), emit: filtered_bam
    tuple val(meta), path("*.ribo.bam"), emit: ribo_bam

    script:
    """
    samtools index ${sorted_bam}

    bedtools intersect -v -abam ${sorted_bam} -b ${rrna_bed} > ${meta.sample}.filtered.bam

    bedtools intersect -abam ${sorted_bam} -b ${rrna_bed} > ${meta.sample}.ribo.bam

    """
}