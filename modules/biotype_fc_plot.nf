process BIOTYPE_COMPOSITION_PLOT {
    label 'demand_1'
    container 'community.wave.seqera.io/library/pip_matplotlib_pandas:4c06d17aa1014f89'

    input:
    path csv
    val level

    output:
    path "*.pdf", emit: pdf
    path "*.png", emit: png

    script:
    def prefix = csv.baseName.replace("featurecounts_summary", "biotype_composition")
    """
    plot_biotypes.py \\
        --input ${csv} --level ${level} \\
        --output-prefix ${prefix}
    """
}