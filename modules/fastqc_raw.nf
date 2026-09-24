process FASTQC_RAW {
    tag "$meta.sample"
    label 'demand_2'
    container 'community.wave.seqera.io/library/fastqc:0.12.1--9971ea336a9eddae'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.{html,zip}"),                      emit: qc

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}