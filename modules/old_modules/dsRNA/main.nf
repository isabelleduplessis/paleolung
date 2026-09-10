process dsRNA_arms_intersect {

    container 'community.wave.seqera.io/library/bedtools_pysam_samtools:a95b7b995eeeb173'

    input:
    tuple val(meta), path(bam)
    path dsRNAome_bed

    output:
    tuple val(meta), path("${meta.sample_id}.feature_assignment.tsv"), emit: feature_tsv
    tuple val(meta), path("${meta.sample_id}"),                        emit: sample_outdir

    script:
    def sample = "${meta.sample_id}"
    """
    assign_features.py \
        -b ${bam} \
        -f ${dsRNAome_bed} \
        -s 0 \
        --outdir ${sample}

    mv ${sample}/feature_assignment.tsv ${sample}.feature_assignment.tsv
    """
}


process dsRNA_arms_summary {
    input:
    path tsv_files

    output:
    path "dsRNAarms_summary.tsv", emit: summary

    script:
    """
    awk 'FNR==1 && NR!=1 {next} 1' *.tsv > dsRNAarms_summary.tsv
    """
}

process dsRNA_mapdamage {

    container 'community.wave.seqera.io/library/mapdamage2:2.2.3--65457d747ad55400'

    input:
    tuple val(meta), path(feat_bam)
    path hg_ref

    output:
    tuple val(meta), path("*.damage"), emit: damage

    script:
    def sample = "${meta.sample_id}"
    def feat = feat_bam.baseName.replace("${sample}.", "")
    """

    mapDamage \
        -i ${feat_bam} \
        -r ${hg_ref} \
        -d ${feat}.damage \
        --merge-reference-sequences \
        --no-stats
    """
}