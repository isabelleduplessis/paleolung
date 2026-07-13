process dsRNA_arms_intersect {

    container 'community.wave.seqera.io/library/bedtools_pysam_samtools:a95b7b995eeeb173'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample_id}.feature_assignment.tsv"), emit: feature_tsv
    tuple val(meta), path("${meta.sample_id}"),                        emit: sample_outdir

    script:
    def sample = "${meta.sample_id}"
    """
    python assign_features.py \
        -b ${bam} \
        -f ${params.dsRNAome} \
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

    container 'community.wave.seqera.io/library/pip_mapdamage:596dc5f355e3d192'

    input:
    tuple val(meta), path(feat_bam)

    output:
    tuple val(meta), path("*.damage"), emit: damage

    script:
    def sample = "${meta.sample_id}"
    def feat = feat_bam.baseName.replace("${sample}.", "")
    """

    mapDamage \
        -i ${feat_bam} \
        -r ${params.human_bwa_ref} \
        -d ${feat}.damage \
        --merge-reference-sequences \
        --no-stats
    """
}