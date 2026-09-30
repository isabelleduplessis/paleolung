process FASTP {
    tag "$meta.sample"
    label 'demand_3'
    container 'community.wave.seqera.io/library/fastp:1.3.6--5a6797673f0eb245'

    input:
    tuple val(meta), path(reads)
    val PE
    val trim_length
    val trim_quality
    val merge_PE
    val min_overlap

    output:
    tuple val(meta), path("${meta.sample}*.trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(meta), path("${meta.sample}.merged.fastq.gz"), emit: merged_reads, optional: true
    path "*.fastp.json", emit: json
    path "*.fastp.html", emit: html

    script:
    def read1 = reads[0]
    def read2 = PE ? "-I ${reads[1]}" : ''
    def output1 = PE ? "${meta.sample}_R1.trimmed.fastq.gz" : "${meta.sample}.trimmed.fastq.gz"
    def output2 = PE ? "-O ${meta.sample}_R2.trimmed.fastq.gz" : ''
    def prefix = "${meta.sample}"
    def merge_options = merge_PE
        ? "--merge --merged_out ${prefix}.merged.fastq.gz --include_unmerged --overlap_len_require ${min_overlap}"
        : ''
    """
    fastp \
        -i ${read1} \
        ${read2} \
        -o ${output1} \
        ${output2} \
        -j ${prefix}.fastp.json \
        -h ${prefix}.fastp.html \
        -q ${trim_quality} \
        -l ${trim_length} \
        --trim_poly_g \
        --trim_poly_x \
        ${merge_options} \
        --thread ${task.cpus}
    """
}