process TRIM_GALORE {
    tag "$meta.sample"
    label 'demand_3'
    container 'community.wave.seqera.io/library/trim-galore:2.3.0--2b255ff43100a747'

    input:
    tuple val(meta), path(reads)
    val PE

    output:
    tuple val(meta), path("${meta.sample}*.{fq,fastq}.gz"), emit: trimmed_reads
    tuple val(meta), path("*trimming_report.txt"),          emit: trim_report, optional: true

    script:
    def paired = PE ? "--paired" : ""
    """
    trim_galore ${paired} \
        --length ${params.trim_length} \
        --quality ${params.trim_quality} \
        --basename ${meta.sample} \
        --cores ${task.cpus} \
        ${reads}
    """
}