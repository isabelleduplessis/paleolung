process TAXONOMIC_PROFILING {
    tag "$meta.sample"
    container 'community.wave.seqera.io/library/bracken_kraken2:f7d86e30f284602f'
    label 'demand_4'

    input:
    tuple val(meta), path(reads)
    path db
    val bracken_read_len
    val bracken_tax_level

    output:
    tuple val(meta), path("*.kraken2"),        emit: kraken_out
    tuple val(meta), path("*.kraken2.report"), emit: report
    tuple val(meta), path("*.bracken"),        emit: bracken
    tuple val (meta), path ("*.bracken.kreport"), emit: bracken_kreport
    tuple val (meta), path ("*_bracken_summary.tsv"), emit: summaryfile

    script:
    """
    kraken2 --db ${db} ${reads} --threads ${task.cpus} --minimum-hit-groups 3 \
        --report-minimizer-data --report ${meta.sample}.kraken2.report > ${meta.sample}.kraken2

    bracken -d ${db} -i ${meta.sample}.kraken2.report -o ${meta.sample}.bracken -w ${meta.sample}.bracken.kreport -r ${bracken_read_len} -l ${bracken_tax_level} 

    (
    echo -e "sample\tname\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"
    for f in *.bracken; do
        sample=\$(basename "\$f" .bracken)
        awk -v s="\$sample" 'NR>1 {print s"\t"\$0}' "\$f"
    done
    ) > ${meta.sample}_bracken_summary.tsv

    """
}