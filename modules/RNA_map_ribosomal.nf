process MAP_RIBO {
    tag "${meta.sample}"
    label  'demand_6'
    container 'community.wave.seqera.io/library/bowtie2_samtools:6243cc3671a7ee48'

    input:
    tuple val(meta), path(reads)
    path index_files
    val index_prefix
    val PE

    output:
    tuple val(meta), path("*.sort.bam"), emit: sorted_bam
    tuple val(meta), path("*.ribo.unmapped*.fastq"), emit: unmapped_reads

    script:
    if (PE) {
        """
        bowtie2 -x ${index_prefix} -1 ${reads[0]} -2 ${reads[1]} \
            --end-to-end --sensitive \
            -S ${meta.sample}.sam -a \
            --threads ${task.cpus}

        samtools addreplacerg \
            -r ID:BT2 \
            -r SM:${meta.sample} \
            -o ${meta.sample}.rg.sam ${meta.sample}.sam

        samtools view -bh -F4 \
            -U ${meta.sample}.unmapped.bam \
            ${meta.sample}.rg.sam \
            | samtools sort -o ${meta.sample}.sort.bam

        rm *.sam

        samtools collate -O -u ${meta.sample}.unmapped.bam  | samtools bam2fq -1 ${meta.sample}.ribo.unmapped_1.fastq -2 ${meta.sample}.ribo.unmapped_2.fastq -0 /dev/null -s /dev/null -n -
        """
    } else {
        """
        bowtie2 -x ${index_prefix} -U ${reads} \
            --end-to-end --sensitive \
            -S ${meta.sample}.sam -a \
            --threads ${task.cpus}

        samtools addreplacerg \
            -r ID:BT2 \
            -r SM:${meta.sample} \
            -o ${meta.sample}.rg.sam ${meta.sample}.sam

        samtools view -bh -F4 \
            -U ${meta.sample}.unmapped.bam \
            ${meta.sample}.rg.sam \
            | samtools sort -o ${meta.sample}.sort.bam

        rm *.sam

        samtools bam2fq ${meta.sample}.unmapped.bam > ${meta.sample}.ribo.unmapped.fastq
        """
    }
}