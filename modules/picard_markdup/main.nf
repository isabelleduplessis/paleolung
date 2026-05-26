process rmDup {

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.rmdup.bam"), emit: bam
    path "*.metrics.txt",                 emit: metrics

    script:
    def prefix = "${meta.sample_id}"
    """
    java -jar \$EBROOTPICARD/picard.jar MarkDuplicates \
        -I ${bam} \
        -O ${prefix}.rmdup.bam \
        -REMOVE_DUPLICATES true \
        -M ${prefix}.metrics.txt

    samtools index ${prefix}.rmdup.bam
    """
}