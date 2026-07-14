process mapBWAHumanRibo {

    container 'community.wave.seqera.io/library/bwa_samtools:22527ddb3cc80fb3'

    input:
    tuple val(meta), path(reads)
    path ref_files
    path ref_fasta

    output:
    tuple val(meta), path("*.sort.bam"), emit: sorted_bam
    tuple val(meta), path("*.unmapped.fastq"), emit: unmapped_reads

    script:
    def prefix = "${meta.sample_id}"
    """
    bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} ${ref_fasta} ${reads} \
    | bwa samse ${ref_fasta} -r "@RG\\tID:${prefix}\\tSM:${prefix}" - ${reads} \
    | samtools view -bSh > ${prefix}.bam

    samtools view -bh -F4 -q 30 \
        -U ${prefix}.unmapped.bam \
        ${prefix}.bam \
        | samtools sort -o ${prefix}.sort.bam

    samtools bam2fq ${prefix}.unmapped.bam > ${prefix}.unmapped.fastq

    """
}

process mapBWAChunks {

    container 'community.wave.seqera.io/library/bwa_samtools:22527ddb3cc80fb3'

    input:
    tuple val(meta), path(reads)
    path ref_files
    path ref_fasta

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    def prefix = "${reads.baseName}"

    """
    bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} ${ref_fasta} ${reads} \
    | bwa samse ${ref_fasta} -r "@RG\\tID:${prefix}\\tSM:${prefix}" - ${reads} \
    | samtools view -bSh -q 30 > ${prefix}.bam
    # added mapq30 dna bwa filter here, later on it filters by 20 for RNA, but this will protect dna

    """
}

