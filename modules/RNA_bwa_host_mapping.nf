process RNA_BWA_MAPPING {
    tag "${meta.sample}.chunk${meta.chunk}"
    label 'demand_6'
    container 'community.wave.seqera.io/library/bwa_samtools:22527ddb3cc80fb3'

    input:
    tuple val(meta), path(chunk_reads)
    tuple path(host_reference), path(host_index_files)
    val seed_length
    val map_missing_prob
    val gap_fraction
    val PE

    output:
    tuple val(meta), path("${meta.sample}.chunk${meta.chunk}.bam"), emit: bam

    script:
    def out_bam = "${meta.sample}.chunk${meta.chunk}.bam"
    def rg      = "@RG\\tID:${meta.sample}.chunk${meta.chunk}\\tSM:${meta.sample}"
    if (PE) {
        """
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${chunk_reads[0]} > R1.sai
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${chunk_reads[1]} > R2.sai
        bwa sampe -r "${rg}" ${host_reference} R1.sai R2.sai ${chunk_reads[0]} ${chunk_reads[1]} \\
            | samtools view -bSh > ${out_bam}
        """
    } else {
        """
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${chunk_reads[0]} > reads.sai
        bwa samse -r "${rg}" ${host_reference} reads.sai ${chunk_reads[0]} \\
            | samtools view -bSh > ${out_bam}
        """
    }
}

process RNA_BWA_MAPPING_MERGE {
    tag "${meta.sample}"
    label 'demand_6'
    container 'community.wave.seqera.io/library/bwa_samtools:22527ddb3cc80fb3'

    input:
    tuple val(meta), path(bams)
    val PE
    val mapq
    val max_hits

    output:
    tuple val(meta), path("${meta.sample}.bwa.bam"), emit: sorted_bam
    tuple val(meta), path("${meta.sample}.multi.mapped.bam"), path("${meta.sample}.multi.mapped.bam.bai"), emit: multimapped_bam
    tuple val(meta), path("*.unmapped*.fastq"), emit: unmapped_reads
    tuple val(meta), path("*.multimapped*.fastq"), emit: multimapped_reads

    script:
    def xa_n = max_hits.toInteger() - 1   // XA holds the alternatives; the primary record is the +1
    if (PE) {
        """
        samtools merge -f -@ ${task.cpus} ${meta.sample}.sampe.bam ${bams}

        samtools view -bh -F4 \\
            -U ${meta.sample}.unmapped.bam \\
            -o ${meta.sample}.mapped.bam ${meta.sample}.sampe.bam

        samtools view -bh -q ${mapq} \\
            -U ${meta.sample}.multi.raw.bam \\
            -o ${meta.sample}.unique.mapped.bam ${meta.sample}.mapped.bam

        samtools sort -o ${meta.sample}.bwa.bam ${meta.sample}.unique.mapped.bam

        # expand XA tags into separate alignment records
        samtools view -h ${meta.sample}.multi.raw.bam \\
            | xa2multi.pl \\
            | samtools sort -o ${meta.sample}.multi.mapped.bam -
        samtools index ${meta.sample}.multi.mapped.bam

        samtools collate -O -u ${meta.sample}.unmapped.bam | samtools bam2fq -1 ${meta.sample}.unmapped_1.fastq -2 ${meta.sample}.unmapped_2.fastq -0 /dev/null -s /dev/null -n -
        samtools collate -O -u ${meta.sample}.multi.raw.bam | samtools bam2fq -1 ${meta.sample}.multimapped_1.fastq -2 ${meta.sample}.multimapped_2.fastq -0 /dev/null -s /dev/null -n -
        """
    } else {
        """
        samtools merge -f -@ ${task.cpus} ${meta.sample}.samse.bam ${bams}

        samtools view -bh -F4 \\
            -U ${meta.sample}.unmapped.bam \\
            -o ${meta.sample}.mapped.bam ${meta.sample}.samse.bam

        samtools view -bh -q ${mapq} \\
            -U ${meta.sample}.multi.raw.bam \\
            -o ${meta.sample}.unique.mapped.bam ${meta.sample}.mapped.bam

        samtools sort -o ${meta.sample}.bwa.bam ${meta.sample}.unique.mapped.bam

        # expand XA tags into separate alignment records
        samtools view -h ${meta.sample}.multi.raw.bam \\
            | xa2multi.pl \\
            | samtools sort -o ${meta.sample}.multi.mapped.bam -
        samtools index ${meta.sample}.multi.mapped.bam

        samtools bam2fq ${meta.sample}.unmapped.bam > ${meta.sample}.unmapped.fastq
        samtools bam2fq ${meta.sample}.multi.raw.bam > ${meta.sample}.multimapped.fastq
        """
    }
}