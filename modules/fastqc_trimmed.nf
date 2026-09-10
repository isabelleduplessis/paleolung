process FASTQC_TRIMMED {
    tag "$meta.sample"
    label 'demand_2'
    container 'community.wave.seqera.io/library/fastqc:0.12.1--9971ea336a9eddae'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.{html,zip}"),                          emit: qc
    tuple val(meta), path("${meta.sample}.trimmed_read_count.txt"), emit: read_count

    script:
    """
    fastqc -t ${task.cpus} ${reads}

    total=0
    for f in ${reads}; do
        n=\$(zcat -f "\$f" | wc -l)
        total=\$((total + n / 4))
    done
    echo \$total > ${meta.sample}.trimmed_read_count.txt
    """
}