process mapBWAHumanRibo {
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.sort.bam"), emit: sorted_bam
    tuple val(meta), path("*.unmapped.fastq"), emit: unmapped_reads

    script:
    def prefix = "${meta.sample_id}"
    """
    bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} ${params.ribo_bwa_ref} ${reads} \
    | bwa samse ${params.ribo_bwa_ref} -r "@RG\\tID:${prefix}\\tSM:${prefix}" - ${reads} \
    | samtools view -bSh > ${prefix}.bam

    samtools view -bh -F4 -q 30 \
        -U ${prefix}.unmapped.bam \
        ${prefix}.bam \
        | samtools sort -o ${prefix}.sort.bam

    samtools bam2fq ${prefix}.unmapped.bam > ${prefix}.unmapped.fastq

    """
}

process mapBWAChunks {
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    def prefix = "${reads.baseName}"

    """
    bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} ${params.human_bwa_ref} ${reads} \
    | bwa samse ${params.human_bwa_ref} -r "@RG\\tID:${prefix}\\tSM:${prefix}" - ${reads} \
    | samtools view -bSh > ${prefix}.bam

    """
}

