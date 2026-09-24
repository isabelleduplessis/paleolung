process COUNT_FASTQ_READS {
    tag "${meta.sample}"
    label 'demand_1'

    input:
    tuple val(meta), path(reads)
    val PE

    output:
    tuple val(meta), env('COUNT'), emit: count

    script:
    if (PE) {
        """
        COUNT_1=\$(( \$(zcat -f ${reads[0]} | wc -l) / 4 ))
        COUNT_2=\$(( \$(zcat -f ${reads[1]} | wc -l) / 4 ))
        COUNT=\$(( COUNT_1 + COUNT_2 ))
        """
    } else {
        """
        COUNT=\$(( \$(zcat -f ${reads} | wc -l) / 4 ))
        """
    }
}