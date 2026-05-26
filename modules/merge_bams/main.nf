process mergeBams {
    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.sort.bam"),       emit: sorted_bam
    tuple val(meta), path("*.unmapped.fastq"), emit: unmapped_reads

    script:
    def prefix = "${meta.sample_id}"

    """
    samtools merge -f -@ ${task.cpus} ${prefix}.merged.bam ${bams}

    samtools view -bh -F4 -q 20 ${prefix}.merged.bam \
        | samtools sort -o ${prefix}.sort.bam

    samtools view -bh -F4 -U ${prefix}.unmapped.bam ${prefix}.merged.bam > /dev/null

    samtools bam2fq ${prefix}.unmapped.bam > ${prefix}.unmapped.fastq

    rm ${prefix}.merged.bam

    """
}
