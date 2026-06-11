process mapDamageHuman {
    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam.baseName}_damage/"), emit: damage_outdir

    script:
    """
    mapDamage \
        -i ${bam} \
        -r ${params.human_bwa_ref} \
        -d ${bam.baseName}_damage \
        --merge-reference-sequences \
        --no-stats
    """
}