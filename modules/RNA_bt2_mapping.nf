process RNA_BT2_MAPPING {
    tag "${meta.sample}"
    label  'demand_6'
    container 'community.wave.seqera.io/library/bowtie2_samtools:6243cc3671a7ee48'

    input:
    tuple val(meta), path(reads)
    path index_files
    val index_prefix
    val PE
    val mapq

    output:
    tuple val(meta), path("*.bt2.bam"), emit: sorted_bam
    tuple val(meta), path("*.unmapped*.fastq"), emit: unmapped_reads
    tuple val(meta), path("*.multimapped*.fastq"), emit: multimapped_reads
    tuple val(meta), path("${meta.sample}.multi.mapped.bam"), emit: multimapped_bam

    script:
    if (PE) {
        """
        bowtie2 -k 10 -x ${index_prefix} -1 ${reads[0]} -2 ${reads[1]} \
            --end-to-end --sensitive \
            -S ${meta.sample}.sam \
            --threads ${task.cpus}

        samtools addreplacerg \
            -r ID:BT2 \
            -r SM:${meta.sample} \
            -o ${meta.sample}.rg.sam ${meta.sample}.sam

        samtools view -bh -F4 \
            -U ${meta.sample}.unmapped.bam \
            -o ${meta.sample}.mapped.bam ${meta.sample}.rg.sam

        samtools view -bh -q ${mapq} -F 0x100 \
            -U ${meta.sample}.multi.mapped.bam \
            -o ${meta.sample}.unique.mapped.bam ${meta.sample}.mapped.bam

        samtools sort -o ${meta.sample}.bt2.bam ${meta.sample}.unique.mapped.bam

        samtools collate -O -u ${meta.sample}.unmapped.bam | samtools bam2fq -1 ${meta.sample}.unmapped_1.fastq -2 ${meta.sample}.unmapped_2.fastq -0 /dev/null -s /dev/null -n -
        samtools collate -O -u ${meta.sample}.multi.mapped.bam | samtools bam2fq -1 ${meta.sample}.multimapped_1.fastq -2 ${meta.sample}.multimapped_2.fastq -0 /dev/null -s /dev/null -n -

        rm *.sam
        """
    } else {
        """
        bowtie2 -k 10 -x ${index_prefix} -U ${reads} \
            --end-to-end --sensitive \
            -S ${meta.sample}.sam \
            --threads ${task.cpus}

        samtools addreplacerg \
            -r ID:BT2 \
            -r SM:${meta.sample} \
            -o ${meta.sample}.rg.sam ${meta.sample}.sam

        samtools view -bh -F4 \
            -U ${meta.sample}.unmapped.bam \
            -o ${meta.sample}.mapped.bam ${meta.sample}.rg.sam

        samtools view -bh -q ${mapq} -F 0x100 \
            -U ${meta.sample}.multi.mapped.bam \
            -o ${meta.sample}.unique.mapped.bam ${meta.sample}.mapped.bam

        samtools sort -o ${meta.sample}.bt2.bam ${meta.sample}.unique.mapped.bam

        samtools bam2fq ${meta.sample}.unmapped.bam > ${meta.sample}.unmapped.fastq
        samtools bam2fq ${meta.sample}.multi.mapped.bam > ${meta.sample}.multimapped.fastq

        rm *.sam
        """
    }
}