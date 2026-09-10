process mapDamageHuman {

    container 'community.wave.seqera.io/library/mapdamage2:2.2.3--65457d747ad55400'

    input:
    tuple val(meta), path(bam)
    path hg_ref

    output:
    tuple val(meta), path("${bam.baseName}_damage/"), emit: damage_outdir

    script:
    """
    mapDamage \
        -i ${bam} \
        -r ${hg_ref} \
        -d ${bam.baseName}_damage \
        --merge-reference-sequences \
        --no-stats
    """
}