process RNA_BWA_MAPPING {
    tag "${meta.sample}"
    label 'demand_6'
    container 'community.wave.seqera.io/library/bwa_samtools:22527ddb3cc80fb3'

    input:
    tuple val(meta), path(reads)
    tuple path(host_reference), path(host_index_files)
    val seed_length
    val map_missing_prob
    val gap_fraction
    val PE
    val mapq

    output:
    tuple val(meta), path("${meta.sample}.bwa.bam"), emit: sorted_bam
    tuple val(meta), path("*.unmapped*.fastq"), emit: unmapped_reads
    tuple val(meta), path("*.multimapped*.fastq"), emit: multimapped_reads

    script:
    def rg = "@RG\\tID:${meta.sample}\\tSM:${meta.sample}"
    if (PE) {
        """
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${reads[0]} > R1.sai
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${reads[1]} > R2.sai
        bwa sampe -r "${rg}" ${host_reference} R1.sai R2.sai ${reads[0]} ${reads[1]} \
            | samtools view -bSh > ${meta.sample}.sampe.bam

        samtools view -bh -F4 \
            -U ${meta.sample}.unmapped.bam \
            -o ${meta.sample}.mapped.bam ${meta.sample}.sampe.bam

        samtools view -bh -q ${mapq} \
            -U ${meta.sample}.multi.mapped.bam \
            -o ${meta.sample}.unique.mapped.bam ${meta.sample}.mapped.bam

        samtools sort -o ${meta.sample}.bwa.bam ${meta.sample}.unique.mapped.bam

        samtools collate -O -u ${meta.sample}.unmapped.bam | samtools bam2fq -1 ${meta.sample}.unmapped_1.fastq -2 ${meta.sample}.unmapped_2.fastq -0 /dev/null -s /dev/null -n -
        samtools collate -O -u ${meta.sample}.multi.mapped.bam | samtools bam2fq -1 ${meta.sample}.multimapped_1.fastq -2 ${meta.sample}.multimapped_2.fastq -0 /dev/null -s /dev/null -n -
        """
    } else {
        """
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${reads} > reads.sai
        bwa samse -r "${rg}" ${host_reference} reads.sai ${reads} \
            | samtools view -bSh > ${meta.sample}.samse.bam

        samtools view -bh -F4 \
            -U ${meta.sample}.unmapped.bam \
            -o ${meta.sample}.mapped.bam ${meta.sample}.samse.bam

        samtools view -bh -q ${mapq} \
            -U ${meta.sample}.multi.mapped.bam \
            -o ${meta.sample}.unique.mapped.bam ${meta.sample}.mapped.bam

        samtools sort -o ${meta.sample}.bwa.bam ${meta.sample}.unique.mapped.bam

        samtools bam2fq ${meta.sample}.unmapped.bam > ${meta.sample}.unmapped.fastq
        samtools bam2fq ${meta.sample}.multi.mapped.bam > ${meta.sample}.multimapped.fastq
        """
    }
}