process RNA_BT2_NCRNA_MAPPING {
    tag "${meta.sample} ${target}"
    label 'demand_6'
    container 'community.wave.seqera.io/library/bowtie2_samtools:6243cc3671a7ee48'

    input:
    tuple val(meta), path(reads)
    path index_files
    val index_prefix
    val target
    val PE

    output:
    tuple val(meta), path("${meta.sample}.${target}.sort.bam"), emit: sorted_bam
    tuple val(meta), path("${meta.sample}.${target}.unmapped*.fastq"), emit: unmapped_reads

    script:
    if (PE) {
        """
        bowtie2 -x ${index_prefix} -1 ${reads[0]} -2 ${reads[1]} \\
            --end-to-end --sensitive \\
            -S ${meta.sample}.${target}.sam -a \\
            --threads ${task.cpus}

        samtools addreplacerg \\
            -r ID:BT2 \\
            -r SM:${meta.sample} \\
            -o ${meta.sample}.${target}.rg.sam ${meta.sample}.${target}.sam

        samtools view -bh -F4 \\
            -U ${meta.sample}.${target}.unmapped.bam \\
            ${meta.sample}.${target}.rg.sam \\
            | samtools sort -o ${meta.sample}.${target}.sort.bam

        rm -f *.sam

        samtools collate -O -u ${meta.sample}.${target}.unmapped.bam | samtools bam2fq -1 ${meta.sample}.${target}.unmapped_1.fastq -2 ${meta.sample}.${target}.unmapped_2.fastq -0 /dev/null -s /dev/null -n -
        """
    } else {
        """
        bowtie2 -x ${index_prefix} -U ${reads} \\
            --end-to-end --sensitive \\
            -S ${meta.sample}.${target}.sam -a \\
            --threads ${task.cpus}

        samtools addreplacerg \\
            -r ID:BT2 \\
            -r SM:${meta.sample} \\
            -o ${meta.sample}.${target}.rg.sam ${meta.sample}.${target}.sam

        samtools view -bh -F4 \\
            -U ${meta.sample}.${target}.unmapped.bam \\
            ${meta.sample}.${target}.rg.sam \\
            | samtools sort -o ${meta.sample}.${target}.sort.bam

        rm -f *.sam

        samtools bam2fq ${meta.sample}.${target}.unmapped.bam > ${meta.sample}.${target}.unmapped.fastq
        """
    }
}
