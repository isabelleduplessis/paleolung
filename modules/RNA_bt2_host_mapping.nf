process RNA_BT2_MAPPING {
    tag "${meta.sample}"
    label  'demand_6'
    container 'community.wave.seqera.io/library/bowtie2_samtools:6243cc3671a7ee48'

    input:
    tuple val(meta), path(reads)
    path index_files
    val index_prefix
    val PE
    val max_hits

    output:
    tuple val(meta), path("*.bt2.bam"), emit: sorted_bam

    script:
    if (PE) {
        """
        bowtie2 -k ${max_hits} -x ${index_prefix} -1 ${reads[0]} -2 ${reads[1]} \
            --end-to-end --sensitive \
            -S ${meta.sample}.sam \
            --threads ${task.cpus}

        samtools addreplacerg \
            -r ID:BT2 \
            -r SM:${meta.sample} \
            -o ${meta.sample}.rg.sam ${meta.sample}.sam

        samtools sort -@ ${task.cpus} -o ${meta.sample}.bt2.bam ${meta.sample}.rg.sam

        rm *.sam
        """
    } else {
        """
        bowtie2 -k ${max_hits} -x ${index_prefix} -U ${reads} \
            --end-to-end --sensitive \
            -S ${meta.sample}.sam \
            --threads ${task.cpus}

        samtools addreplacerg \
            -r ID:BT2 \
            -r SM:${meta.sample} \
            -o ${meta.sample}.rg.sam ${meta.sample}.sam

        samtools sort -@ ${task.cpus} -o ${meta.sample}.bt2.bam ${meta.sample}.rg.sam

        rm *.sam
        """
    }
}