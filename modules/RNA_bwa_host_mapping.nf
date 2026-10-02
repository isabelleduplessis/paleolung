process RNA_BWA_MAPPING {
    tag "${meta.sample}.${meta.chunk}"
    label 'demand_6'
    container 'community.wave.seqera.io/library/bwa_samtools:22527ddb3cc80fb3'

    input:
    tuple val(meta), path(chunk_reads)
    tuple path(host_reference), path(host_index_files)
    val seed_length
    val map_missing_prob
    val gap_fraction
    val PE

    output:
    tuple val(meta), path("${meta.sample}.chunk${meta.chunk}.bam"), emit: sorted_bam

    script:
    def out_bam = "${meta.sample}.chunk${meta.chunk}.bam"
    def sorted_bam = "${meta.sample}.chunk${meta.chunk}.sorted.bam"
    def rg      = "@RG\\tID:${meta.sample}.chunk${meta.chunk}\\tSM:${meta.sample}"
    if (PE) {
        """
        set -o pipefail
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${chunk_reads[0]} > R1.sai
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${chunk_reads[1]} > R2.sai
        bwa sampe -r "${rg}" ${host_reference} R1.sai R2.sai ${chunk_reads[0]} ${chunk_reads[1]} \\
            | samtools view -bSh > ${out_bam}

        samtools view -h ${out_bam} \\
            | xa2multi.pl \\
            | samtools sort -@ ${task.cpus} -o ${sorted_bam} -
        mv ${sorted_bam} ${out_bam}
        """
    } else {
        """
        set -o pipefail
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${chunk_reads[0]} > reads.sai
        bwa samse -r "${rg}" ${host_reference} reads.sai ${chunk_reads[0]} \\
            | samtools view -bSh > ${out_bam}

        samtools view -h ${out_bam} \\
            | xa2multi.pl \\
            | samtools sort -@ ${task.cpus} -o ${sorted_bam} -
        mv ${sorted_bam} ${out_bam}
        """
    }
}
