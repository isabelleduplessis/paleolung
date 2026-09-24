process HOST_DNA_MAPPING {
    tag "${meta.sample}.chunk${meta.chunk}"
    label 'demand_6'
    container 'community.wave.seqera.io/library/bwa_samtools:22527ddb3cc80fb3'

    input:
    tuple val(meta), path(chunk_reads)
    tuple path(host_reference), path(host_index_files)
    val PE
    val seed_length
    val map_missing_prob
    val gap_fraction

    output:
    tuple val(meta), path("${meta.sample}.chunk${meta.chunk}.bam"), emit: bam

    script:
    def out_bam = "${meta.sample}.chunk${meta.chunk}.bam"
    def rg      = "@RG\\tID:${meta.sample}.chunk${meta.chunk}\\tSM:${meta.sample}"
    if (PE) {
        """
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${chunk_reads[0]} > R1.sai
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${chunk_reads[1]} > R2.sai
        bwa sampe -r "${rg}" ${host_reference} R1.sai R2.sai ${chunk_reads[0]} ${chunk_reads[1]} \
            | samtools view -bSh > ${out_bam}
        """
    } else {
        """
        bwa aln -l ${seed_length} -n ${map_missing_prob} -o ${gap_fraction} -t ${task.cpus} ${host_reference} ${chunk_reads} > reads.sai
        bwa samse -r "${rg}" ${host_reference} reads.sai ${chunk_reads} \
            | samtools view -bSh > ${out_bam}
        """
    }
}