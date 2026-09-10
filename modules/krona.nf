process KRONA {
    tag "$meta.sample"
    container 'community.wave.seqera.io/library/krakentools_krona_pandas:1411638929ab7a28'
    label 'demand_1'

    input:
    tuple val(meta), path(bracken_kreport)
    val min_fraction

    output:
    tuple val(meta), path("*.filtered.kreport"), emit: filtered_kreport
    tuple val(meta), path("*.krona.txt"),        emit: krona_txt
    tuple val(meta), path("*.krona.html"),       emit: html

    script:
    """
    filter_bracken.py \\
        --input ${bracken_kreport} \\
        --threshold ${min_fraction} \\
        --output ${meta.sample}.filtered.kreport

    kreport2krona.py \\
        -r ${meta.sample}.filtered.kreport \\
        -o ${meta.sample}.krona.txt

    ktImportText \\
        ${meta.sample}.krona.txt \\
        -o ${meta.sample}.krona.html
    """
}