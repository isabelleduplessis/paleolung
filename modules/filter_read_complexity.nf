process FILTER_COMPLEXITY {
    tag "$meta.sample"
    label 'demand_2'
    container 'community.wave.seqera.io/library/bbmap:39.81--6de2cad4dc0eb52c'

    input:
    tuple val(meta), path(reads)
    val PE

    output:
    tuple val(meta), path("${meta.sample}.highC*.fastq.gz", arity: '1..*'), emit: filtered_fastq

    script:
    if (PE) {
        """
        bbduk.sh \
            in1=${reads[0]} \
            in2=${reads[1]} \
            out1=${meta.sample}.highC_R1.fastq.gz \
            out2=${meta.sample}.highC_R2.fastq.gz \
            entropy=0.90 \
            entropywindow=20 \
            entropyk=5 \
            threads=${task.cpus}

        """
    } else {
        """
        bbduk.sh \
            in=${reads} \
            out=${meta.sample}.highC.fastq.gz \
            entropy=0.90 \
            entropywindow=20 \
            entropyk=5 \
            threads=${task.cpus}

        """
    }
}