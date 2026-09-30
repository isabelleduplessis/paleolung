process SPADES {
    tag "${meta.sample}"
    label 'demand_6'
    container 'community.wave.seqera.io/library/spades:4.3.0--fefe8aa442d71a4d'

    input:
    tuple val(meta), path(reads)
    val PE

    output:
    tuple val(meta), path("${meta.sample}_spades_assembly"), emit: assembly

    script:
    def prefix = "${meta.sample}"
    if (PE) {
        """
        spades.py \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o ${prefix}_spades_assembly \\
            --threads ${task.cpus} \\
            --careful
        """
    } else {
        """
        spades.py \\
            -s ${reads[0]} \\
            -o ${prefix}_spades_assembly \\
            --threads ${task.cpus} \\
            --careful
        """
    }
}
