process seqkit_stats_highC {

    container 'community.wave.seqera.io/library/seqkit:2.13.0--4cc296e17e90b043'

    input:
    path fastq_files

    output:
    path "highC_stats.txt", emit: statsfile

    script:
    """
    seqkit stats *.highC.unmapped.fastq > highC_stats.txt
    """
}