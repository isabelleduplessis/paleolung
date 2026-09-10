process MARK_DUPLICATES {
    tag "$meta.sample"
    label 'demand_2'
    container 'community.wave.seqera.io/library/picard_samtools:42878316047c4b3f'

    input:
    tuple val(meta), path(sorted_bam)

    output:
    tuple val(meta), path("${meta.sample}.dedup.bam"), path("${meta.sample}.dedup.bam.bai"), emit: dedup_bam
    tuple val(meta), path("${meta.sample}.dedup_read_count.txt"),                             emit: dedup_read_count
    tuple val(meta), path("${meta.sample}.metrics.txt"),                                      emit: dedup_metrics

    // still want to make a better way of making sure there is enough memory

    script:
    """
    picard -Xmx${(task.memory.toGiga() - 2)}g MarkDuplicates \
        -I ${sorted_bam} \
        -O ${meta.sample}.dedup.bam \
        -REMOVE_DUPLICATES true \
        -M ${meta.sample}.metrics.txt

    samtools index ${meta.sample}.dedup.bam
    samtools view -c ${meta.sample}.dedup.bam > ${meta.sample}.dedup_read_count.txt
    """
}