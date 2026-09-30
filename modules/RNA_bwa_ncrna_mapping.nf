process RNA_BWA_NCRNA_MAPPING {
    tag "${meta.sample} ${target}"
    label 'demand_6'
    container 'community.wave.seqera.io/library/bwa_samtools:22527ddb3cc80fb3'

    input:
    tuple val(meta), path(reads)
    tuple path(reference), path(index_files)
    val target
    val PE
    val seed_length
    val map_missing_prob
    val gap_fraction

    output:
    tuple val(meta), path("${meta.sample}.${target}.sort.bam"), emit: sorted_bam
    tuple val(meta), path("${meta.sample}.${target}.unmapped*.fastq"), emit: unmapped_reads

    script:
    def rg = "@RG\\tID:${meta.sample}\\tSM:${meta.sample}"
    if (PE) {
        """
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${reference} ${reads[0]} > R1.sai
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${reference} ${reads[1]} > R2.sai
        bwa sampe -r "${rg}" ${reference} R1.sai R2.sai ${reads[0]} ${reads[1]} \\
            | samtools view -bSh > ${meta.sample}.${target}.bam

        samtools view -bh -F4 \\
            -U ${meta.sample}.${target}.unmapped.bam \\
            ${meta.sample}.${target}.bam \\
            | samtools sort -o ${meta.sample}.${target}.sort.bam

        samtools collate -O -u ${meta.sample}.${target}.unmapped.bam | samtools bam2fq -1 ${meta.sample}.${target}.unmapped_1.fastq -2 ${meta.sample}.${target}.unmapped_2.fastq -0 /dev/null -s /dev/null -n -
        rm -f *.sai ${meta.sample}.${target}.bam
        """
    } else {
        """
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${reference} ${reads} > reads.sai
        bwa samse -r "${rg}" ${reference} reads.sai ${reads} \\
            | samtools view -bSh > ${meta.sample}.${target}.bam

        samtools view -bh -F4 \\
            -U ${meta.sample}.${target}.unmapped.bam \\
            ${meta.sample}.${target}.bam \\
            | samtools sort -o ${meta.sample}.${target}.sort.bam

        samtools bam2fq ${meta.sample}.${target}.unmapped.bam > ${meta.sample}.${target}.unmapped.fastq
        rm -f reads.sai ${meta.sample}.${target}.bam
        """
    }
}
