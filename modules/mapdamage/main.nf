process mapDamageHuman {
    input:
    path bam

    output:
    path "${bam.baseName}_damage/", emit: damage_outdir

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