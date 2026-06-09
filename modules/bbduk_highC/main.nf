process highC {
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.highC.unmapped.fastq"), emit: filtered_reads

    script:
    def prefix = "${meta.sample_id}"
    """
    bbduk.sh in=${reads} out=${prefix}.highC.unmapped.fastq \
        entropy=0.90 entropywindow=20 entropyk=5
    """
}