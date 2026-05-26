process mapBowtie2HumanRibo {
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.sort.bam"), emit: sorted_bam
    tuple val(meta), path ("*.unmapped.fastq"), emit: unmapped_reads

    script:
    def prefix = "${meta.sample_id}"
    """
    bowtie2 -x ${params.ribo_bt2_ref} -U ${reads} \
        --end-to-end --sensitive \
        -S ${prefix}.sam \
        --threads ${task.cpus}

    samtools addreplacerg \
        -r ID:BT2 \
        -r SM:${prefix} \
        -o ${prefix}.rg.sam ${prefix}.sam

    samtools view -bh -F4 -q 20 \
        -U ${prefix}.unmapped.bam \
        ${prefix}.rg.sam \
        | samtools sort -o ${prefix}.sort.bam

    rm *.sam

    samtools bam2fq ${prefix}.unmapped.bam > ${prefix}.unmapped.fastq

    """
}

process mapBowtie2Chunks {
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    def prefix = "${reads.baseName}"

    """
    bowtie2 -x ${params.human_bt2_ref} -U ${reads} \
        --end-to-end --sensitive \
        -S ${prefix}.sam \
        --threads ${task.cpus}

    samtools addreplacerg \
        -r ID:BT2 \
        -r SM:${prefix} \
        -o ${prefix}.rg.sam ${prefix}.sam

    samtools view -bSh -o ${prefix}.bam ${prefix}.rg.sam

    rm *.sam

    """
}

