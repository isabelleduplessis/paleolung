process trim {
    input:
    tuple val(meta), path(raw_reads)

    output:
    tuple val(meta), path("*.trimmed.fastq.gz"), emit: trimmed_reads
    path "*.fastp.json",                          emit: json
    path "*.fastp.html",                          emit: html

    script:
    def prefix = "${meta.sample_id}"
    """
    fastp -i ${raw_reads} -o ${prefix}.trimmed.fastq.gz \\
    -j ${prefix}.fastp.json \\
    -h ${prefix}.fastp.html \\
    -q 33 \\
    -l ${params.trim_length} \\
    --trim_poly_g \\
    --trim_poly_x \\
    """
}
