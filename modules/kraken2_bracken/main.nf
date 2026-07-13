process kraken {

    container 'community.wave.seqera.io/library/bracken_kraken2:f7d86e30f284602f'

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.kraken2"),        emit: kraken_out
    tuple val(meta), path("*.kraken2.report"), emit: report
    tuple val(meta), path("*.bracken"),        emit: bracken

    script:
    def prefix = "${meta.sample_id}"
    """
    kraken2 --db ${db} ${reads} --threads ${task.cpus} --minimum-hit-groups 3 \
        --report-minimizer-data --report ${prefix}.kraken2.report > ${prefix}.kraken2

    bracken -d ${db} -i ${prefix}.kraken2.report -o ${prefix}.bracken -r 100 -l S
    """
}

process bracken_summary {

    container 'community.wave.seqera.io/library/bracken_kraken2:f7d86e30f284602f'

    input:
    path bracken_files

    output:
    path "bracken_summary.tsv", emit: summaryfile

    script:
    """
    (
    echo -e "sample\tname\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"
    for f in *.bracken; do
        sample=\$(basename "\$f" .bracken)
        awk -v s="\$sample" 'NR>1 {print s"\t"\$0}' "\$f"
    done
    ) > bracken_summary.tsv
    """
}