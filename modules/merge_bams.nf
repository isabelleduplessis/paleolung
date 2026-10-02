process MERGE_BAMS {
    tag "$meta.sample"
    label 'demand_3'
    container 'community.wave.seqera.io/library/samtools:1.24--a1a83d6551ff8faa'

    input:
    tuple val(meta), path(bams)
    val PE
    val mapq

    output:
    tuple val(meta), path("${meta.sample}.sort.bam"),                     emit: sorted_bam
    tuple val(meta), path("*unmapped*.fastq.gz", arity: '1..*'),          emit: unmapped_reads
    tuple val(meta), path("${meta.sample}.multi.mapped.bam"), path("${meta.sample}.multi.mapped.bam.bai"), emit: multimapped_bam, optional: true
    tuple val(meta), path("*.multimapped*.fastq"),                       emit: multimapped_reads, optional: true
    script:
    def merge_and_extract = """
        set -o pipefail
        samtools merge -f -@ ${task.cpus} ${meta.sample}.merged.bam ${bams}

        samtools view -bh -F 4 -o ${meta.sample}.mapped.bam ${meta.sample}.merged.bam
        samtools sort -@ ${task.cpus} -o ${meta.sample}.sort.bam ${meta.sample}.mapped.bam

        samtools view -bh -f 4 -F 0x100 \\
            -o ${meta.sample}.unmapped.bam ${meta.sample}.merged.bam
    """
    def unmapped_to_fastq = PE
        ? """
            samtools sort -n -@ ${task.cpus} ${meta.sample}.unmapped.bam \\
                | samtools bam2fq \\
                    -1 ${meta.sample}_unmapped_R1.fastq \\
                    -2 ${meta.sample}_unmapped_R2.fastq \\
                    -0 /dev/null \\
                    -s /dev/null \\
                    --threads ${task.cpus} \\
                    -
                    gzip ${meta.sample}_unmapped_R1.fastq ${meta.sample}_unmapped_R2.fastq
        """
        : """
            samtools bam2fq ${meta.sample}.unmapped.bam \\
                | gzip > ${meta.sample}.unmapped.fastq.gz
        """

    if (params.sample_type == 'DNA') {
        """
        ${merge_and_extract}
        ${unmapped_to_fastq}
        rm ${meta.sample}.mapped.bam
        """
    } else {
        """
        ${merge_and_extract}
        ${unmapped_to_fastq}

        samtools view -bh -q ${mapq} -F 0x100 \\
            -U ${meta.sample}.multi.raw.bam \\
            -o ${meta.sample}.unique.primary.bam ${meta.sample}.mapped.bam
        samtools sort -@ ${task.cpus} \\
            -o ${meta.sample}.multi.mapped.bam ${meta.sample}.multi.raw.bam
        samtools index -@ ${task.cpus} ${meta.sample}.multi.mapped.bam

        if [ "${PE}" = "true" ]; then
            samtools collate -O -u ${meta.sample}.multi.mapped.bam \\
                | samtools bam2fq \\
                    -1 ${meta.sample}.multimapped_R1.fastq \\
                    -2 ${meta.sample}.multimapped_R2.fastq \\
                    -0 ${meta.sample}.multimapped_other.fastq \\
                    -s ${meta.sample}.multimapped_single.fastq \\
                    -n -
        else
            samtools bam2fq ${meta.sample}.multi.mapped.bam > ${meta.sample}.multimapped.fastq
        fi
        rm ${meta.sample}.mapped.bam ${meta.sample}.multi.raw.bam ${meta.sample}.unique.primary.bam
        """
    }
}