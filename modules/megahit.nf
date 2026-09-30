process MEGAHIT {
    tag "${meta.sample}"
    label 'demand_6'
    container 'community.wave.seqera.io/library/megahit:1.2.9--90a5f977c5aa9d92'

    input:
    tuple val(meta), path(reads)
    val PE

    output:
    tuple val(meta), path("${meta.sample}_megahit_assembly"), emit: assembly

    script:
    def prefix = "${meta.sample}"
    if (PE) {
        """
        megahit \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -o ${prefix}_megahit_assembly \\
            -t ${task.cpus}
        """
    } else {
        """
        megahit \\
            -r ${reads[0]} \\
            -o ${prefix}_megahit_assembly \\
            -t ${task.cpus}
        """
    }
}
