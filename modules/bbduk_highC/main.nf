process highC {

    container 'community.wave.seqera.io/library/bbmap:39.81--6de2cad4dc0eb52c'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.highC.unmapped.fastq"), emit: filtered_reads

    script:
    def prefix = "${meta.sample_id}"
    """
    bbduk.sh in=${reads} out=${prefix}.highC.unmapped.fastq \
        entropy=0.90 entropywindow=20 entropyk=5 threads=${task.cpus}
    """
}