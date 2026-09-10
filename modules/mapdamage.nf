process MAPDAMAGE {
    tag "$meta.sample"
    label 'demand_1'
    container 'community.wave.seqera.io/library/mapdamage2:2.2.3--65457d747ad55400'

    input:
    tuple val(meta), path(dedup_bam), path(dedup_bai)
    path host_reference

    output:
    tuple val(meta), path("${meta.sample}_mapdamage"), emit: results

    script:
    """
    mapDamage \
        -i ${dedup_bam} \
        -r ${host_reference} \
        -d ${meta.sample}_mapdamage \
        --merge-reference-sequences \
        --no-stats
    """
}