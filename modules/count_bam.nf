process COUNT_BAM_READS {
    tag "${meta.sample}"
    label 'demand_1'
    container 'community.wave.seqera.io/library/samtools:1.21--1929a413d00e1d7a'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), env('COUNT'), emit: count

    script:
    """
    COUNT=\$(samtools view -c -F 0x900 ${bam})
    """
}