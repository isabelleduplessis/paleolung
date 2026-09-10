process INDEX_HOST_REFERENCE {
    label 'demand_1'
    container 'community.wave.seqera.io/library/bwa:0.7.19--5f6c9601d4d7a45a'

    input:
    path host_reference

    output:
    tuple path(host_reference), path("${host_reference}.*"), emit: index

    script:
    """
    bwa index ${host_reference}
    """
}