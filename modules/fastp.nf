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
    tuple val(meta), path("${meta.sample}*.trimmed.fastq.gz"), emit: trimmed_reads, optional: true
    tuple val(meta), path("${meta.sample}.merged.fastq.gz"), emit: merged_reads, optional: true
    path "*.fastp.json", emit: json
    path "*.fastp.html", emit: html

    script:
    def prefix = "${meta.sample}"
    def read1 = reads[0]
    def read2 = PE ? "-I ${reads[1]}" : ''
    
    def outputs
    if (merge_PE) {
        outputs = "--merge --merged_out ${prefix}.merged.fastq.gz --include_unmerged --overlap_len_require ${min_overlap}"
    } else if (PE) {
        outputs = "-o ${prefix}_R1.trimmed.fastq.gz -O ${prefix}_R2.trimmed.fastq.gz"
    } else {
        outputs = "-o ${prefix}.trimmed.fastq.gz"
    }
    """
    fastp \\
        -i ${read1} \\
        ${read2} \\
        ${outputs} \\
        -j ${prefix}.fastp.json \\
        -h ${prefix}.fastp.html \\
        -q ${trim_quality} \\
        -l ${trim_length} \\
        --trim_poly_g \\
        --trim_poly_x \\
        --thread ${task.cpus}
    """
}