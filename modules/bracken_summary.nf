process BRACKEN_SUMMARY {
    container 'community.wave.seqera.io/library/bracken_kraken2:f7d86e30f284602f'
    label 'demand_1'

    input:
    path bracken_files

    output:
    path 'bracken_summary.tsv', emit: summary

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
