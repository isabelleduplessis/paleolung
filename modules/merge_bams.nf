process MERGE_BAMS {
    tag "$meta.sample"
    label 'demand_3'
    container 'community.wave.seqera.io/library/samtools:1.24--a1a83d6551ff8faa'

    input:
    tuple val(meta), path(bams)
    val PE
    val mapq
    val multi

    output:
    tuple val(meta), path("${meta.sample}.sort.bam"),            emit: sorted_bam
    tuple val(meta), path("*unmapped*.fastq.gz", arity: '1..*'), emit: unmapped_reads

    script:
    def mapq_filter = multi ? '' : "-q ${mapq}"

    """
    samtools merge -f -@ ${task.cpus} ${meta.sample}.merged.bam ${bams}

    samtools view -bh -F 4 ${mapq_filter} ${meta.sample}.merged.bam | samtools sort -@ ${task.cpus} -o ${meta.sample}.sort.bam

    samtools view -bh -f 4 -F 0x100 -o ${meta.sample}.unmapped.bam ${meta.sample}.merged.bam

    if [ "${PE}" = "true" ]; then
        samtools sort -n -@ ${task.cpus} ${meta.sample}.unmapped.bam \\
            | samtools bam2fq \\
                -1 ${meta.sample}_unmapped_R1.fastq \\
                -2 ${meta.sample}_unmapped_R2.fastq \\
                -0 /dev/null \\
                -s /dev/null \\
                --threads ${task.cpus} \\
                -
        gzip ${meta.sample}_unmapped_R1.fastq ${meta.sample}_unmapped_R2.fastq
    else
        samtools bam2fq ${meta.sample}.unmapped.bam | gzip > ${meta.sample}.unmapped.fastq.gz
    fi
    """
}