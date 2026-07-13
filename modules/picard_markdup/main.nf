process rmDup {

    container 'community.wave.seqera.io/library/picard_samtools:42878316047c4b3f'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.rmdup.bam"), emit: bam
    path "*.metrics.txt",                 emit: metrics

    script:
    def prefix = "${meta.sample_id}"
    """
    picard MarkDuplicates \
        -I ${bam} \
        -O ${prefix}.rmdup.bam \
        -REMOVE_DUPLICATES true \
        -M ${prefix}.metrics.txt

    samtools index ${prefix}.rmdup.bam
    """
}