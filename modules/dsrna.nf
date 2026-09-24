process DSRNA_INTERSECT {

    tag "$meta.sample"
    label 'demand_1'
    container 'community.wave.seqera.io/library/bedtools_mapdamage2:a199a90b418cb941'

    input:
    tuple val(meta), path(dedup_bam), path(dedup_bai)
    path host_reference
    path dsrna_bed
    val strand

    output:
    tuple val(meta), path("${meta.sample}_dsRNA_mapdamage"), emit: dsRNA_mapdamage
    tuple val(meta), path("${meta.sample}_ssRNA_mapdamage"), emit: ssRNA_mapdamage

    script:
    def sflag = strand == 1 ? '-s' : strand == 2 ? '-S' : ''
    """
    bedtools intersect -a ${dedup_bam} -b ${dsrna_bed} -u ${sflag} > ${meta.sample}.dsRNA.bam
    bedtools intersect -a ${dedup_bam} -b ${dsrna_bed} -v ${sflag} > ${meta.sample}.ssRNA.bam

    mapDamage \
        -i ${meta.sample}.dsRNA.bam \
        -r ${host_reference} \
        -d ${meta.sample}_dsRNA_mapdamage \
        --merge-reference-sequences \
        --no-stats

    mapDamage \
        -i ${meta.sample}.ssRNA.bam \
        -r ${host_reference} \
        -d ${meta.sample}_ssRNA_mapdamage \
        --merge-reference-sequences \
        --no-stats


    mv ${meta.sample}.dsRNA.bam ${meta.sample}_dsRNA_mapdamage
    mv ${meta.sample}.ssRNA.bam ${meta.sample}_ssRNA_mapdamage

    """

}