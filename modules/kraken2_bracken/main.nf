process kraken {
    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.kraken2"),        emit: kraken_out
    tuple val(meta), path("*.kraken2.report"), emit: report
    tuple val(meta), path("*.bracken"),        emit: bracken

    script:
    def prefix = "${meta.sample_id}"
    """
    kraken2 --db ${db} ${reads} --threads ${task.cpus} --minimum-hit-groups 3 \
        --report-minimizer-data --report ${prefix}.kraken2.report > ${prefix}.kraken2

    bracken -d ${db} -i ${prefix}.kraken2.report -o ${prefix}.bracken -r 100 -l S
    """
}