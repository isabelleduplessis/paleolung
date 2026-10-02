#!/usr/bin/env nextflow
include { ANCIENT_DNA_HOST_MAPPING } from './aDNA.nf'
include { ANCIENT_RNA_HOST_MAPPING } from './aRNA.nf'
include { MARK_DUPLICATES } from '../modules/mark_duplicates.nf'
include { MAPDAMAGE } from '../modules/mapdamage.nf'
include { FILTER_COMPLEXITY } from '../modules/filter_read_complexity.nf'
include { READ_LENGTH_DISTRIBUTION } from '../modules/read_length_plot.nf'
include { COUNT_FASTQ_READS as COUNT_FASTQ_READS_TRIMMED } from '../modules/count_fastq.nf'
include { COUNT_FASTQ_READS as COUNT_FASTQ_READS_REMAINDER } from '../modules/count_fastq.nf'
include { COUNT_BAM_READS } from '../modules/count_bam.nf'
include { DSRNA_INTERSECT } from '../modules/dsrna.nf'
include { FEATURE_COUNTS } from '../modules/featurecounts.nf'
include { FEATURECOUNTS_SUMMARY as FEATURECOUNTS_SUMMARY_EXON } from '../modules/featurecounts.nf'
include { FEATURECOUNTS_SUMMARY as FEATURECOUNTS_SUMMARY_GENE } from '../modules/featurecounts.nf'
include { BIOTYPE_COMPOSITION_PLOT as BIOTYPE_COMPOSITION_PLOT_EXON } from '../modules/biotype_fc_plot.nf'
include { BIOTYPE_COMPOSITION_PLOT as BIOTYPE_COMPOSITION_PLOT_GENE } from '../modules/biotype_fc_plot.nf'

workflow HOST_MAPPING {
    take:
    ch_raw_reads
    ch_trimmed_reads
    PE

    main:
    def filter_ribo_value = params.filter_ribo == null ? 'true' : params.filter_ribo.toString().toLowerCase()
    if (!(filter_ribo_value in ['true', 'false'])) {
        error "ERROR: --filter_ribo must be 'true' or 'false' (got: ${params.filter_ribo})"
    }
    def filter_ribo_enabled = filter_ribo_value == 'true'

    // Count trimmed reads once; the DNA chunk mapper also uses these counts
    // to determine when each sample's chunk BAMs can be merged.
    ch_trimmed_for_count = ch_trimmed_reads
        .map { meta, reads -> tuple(meta + [stage:'trimmed'], reads) }
    COUNT_FASTQ_READS_TRIMMED(ch_trimmed_for_count, PE)

    // perform mapping depending on sample type
    if (params.sample_type == 'DNA') {
        ANCIENT_DNA_HOST_MAPPING(ch_trimmed_reads, COUNT_FASTQ_READS_TRIMMED.out.count, PE)
        ch_unmapped_reads = ANCIENT_DNA_HOST_MAPPING.out.unmapped_reads
        ch_host_bam = ANCIENT_DNA_HOST_MAPPING.out.host_bam_with_dups
        ch_ribo_mapped = channel.empty()
        ch_trna_mapped = channel.empty()
        ch_ribo_intersect = channel.empty()
        ch_multimapped_reads = channel.empty()
        ch_multimapped_bam = channel.empty()
    } else {
        ANCIENT_RNA_HOST_MAPPING(ch_trimmed_reads, params.RNA_host_mapper, params.ncRNA_mapper, PE)
        ch_unmapped_reads = ANCIENT_RNA_HOST_MAPPING.out.unmapped_reads
        ch_host_bam = ANCIENT_RNA_HOST_MAPPING.out.host_bam_with_dups
        ch_ribo_mapped = ANCIENT_RNA_HOST_MAPPING.out.ribo_mapped_bam
        ch_trna_mapped = ANCIENT_RNA_HOST_MAPPING.out.trna_mapped_bam
        ch_ribo_intersect = ANCIENT_RNA_HOST_MAPPING.out.ribo_intersect_bam
        ch_multimapped_reads = ANCIENT_RNA_HOST_MAPPING.out.multimapped_reads
        ch_multimapped_bam = ANCIENT_RNA_HOST_MAPPING.out.multimapped_bam
    }

    // host mapping post processing
    MARK_DUPLICATES(ch_host_bam)
    MAPDAMAGE(MARK_DUPLICATES.out.dedup_bam, file(params.host_reference))
    READ_LENGTH_DISTRIBUTION(MARK_DUPLICATES.out.dedup_bam)
    FILTER_COMPLEXITY(ch_unmapped_reads, PE)
    ch_host_bam_publish = params.publish_host_bam
        ? MARK_DUPLICATES.out.dedup_bam
        : channel.empty()

    ch_dsrna_mapdamage = channel.empty()
    ch_ssrna_mapdamage = channel.empty()
    ch_fc_genes = channel.empty()
    ch_fc_exons = channel.empty()
    ch_biotype_genes = channel.empty()
    ch_biotype_exons = channel.empty()

    if (params.sample_type == 'RNA') {
        // RNA specific post processing
        DSRNA_INTERSECT(MARK_DUPLICATES.out.dedup_bam, file(params.host_reference), file(params.dsRNA_bed), params.strand)
        FEATURE_COUNTS(MARK_DUPLICATES.out.dedup_bam, file(params.hg_gtf), params.strand)

        ch_dsrna_mapdamage = DSRNA_INTERSECT.out.dsRNA_mapdamage
        ch_ssrna_mapdamage = DSRNA_INTERSECT.out.ssRNA_mapdamage

        // Collect exon counts + summaries across all samples
        exon_files = FEATURE_COUNTS.out.exon_counts
            .map { _meta, txt -> txt }
            .mix(FEATURE_COUNTS.out.exon_summary.map { _meta, summ -> summ })
            .collect()

        // Collect gene counts + summaries across all samples
        gene_files = FEATURE_COUNTS.out.gene_counts
            .map { _meta, txt -> txt }
            .mix(FEATURE_COUNTS.out.gene_summary.map { _meta, summ -> summ })
            .collect()

        FEATURECOUNTS_SUMMARY_EXON(exon_files, file(params.hg_gtf), 'exon')
        FEATURECOUNTS_SUMMARY_GENE(gene_files, file(params.hg_gtf), 'gene')

        BIOTYPE_COMPOSITION_PLOT_EXON(FEATURECOUNTS_SUMMARY_EXON.out.csv, 'exon')
        BIOTYPE_COMPOSITION_PLOT_GENE(FEATURECOUNTS_SUMMARY_GENE.out.csv, 'gene')

        ch_fc_genes = FEATURECOUNTS_SUMMARY_EXON.out.csv
        ch_fc_exons = FEATURECOUNTS_SUMMARY_GENE.out.csv
        ch_biotype_genes = BIOTYPE_COMPOSITION_PLOT_GENE.out.pdf
        ch_biotype_exons = BIOTYPE_COMPOSITION_PLOT_EXON.out.pdf
    }

    // Generate the read-count report across all samples.
    ch_to_count = ch_raw_reads.map { meta, reads -> tuple(meta + [stage:'raw'], reads) }
        .mix(ch_unmapped_reads.map { meta, reads -> tuple(meta + [stage:'unmapped'], reads) })
        .mix(FILTER_COMPLEXITY.out.filtered_fastq.map { meta, reads -> tuple(meta + [stage:'complexity_filtered'], reads) })
        .mix(ch_multimapped_reads.map { meta, reads -> tuple(meta + [stage:'multimapped'], reads) })

    COUNT_FASTQ_READS_REMAINDER(ch_to_count, PE)
    ch_all_counts = COUNT_FASTQ_READS_REMAINDER.out.count
        .map { meta, count -> tuple(meta.sample, meta.stage, count as Integer) }
        .mix(COUNT_FASTQ_READS_TRIMMED.out.count
            .map { meta, count -> tuple(meta.sample, meta.stage, count as Integer) })

    ch_bam_to_count = ch_host_bam.map { meta, bam -> tuple(meta + [stage:'host_mapped'], bam) }
        .mix(MARK_DUPLICATES.out.dedup_bam.map { meta, bam, _bai -> tuple(meta + [stage:'host_mapped_deduped'], bam) })
        .mix(ch_ribo_mapped.map { meta, bam -> tuple(meta + [stage:'ribo_mapped'], bam) })
        .mix(ch_trna_mapped.map { meta, bam -> tuple(meta + [stage:'trna_mapped'], bam) })
        .mix(ch_ribo_intersect.map { meta, bam -> tuple(meta + [stage:'ribo_intersect'], bam) })

    COUNT_BAM_READS(ch_bam_to_count)
    ch_all_bam_counts = COUNT_BAM_READS.out.count
        .map { meta, count -> tuple(meta.sample, meta.stage, count as Integer) }

    ch_combined_counts = ch_all_counts.mix(ch_all_bam_counts)
    ch_counts_wide = ch_combined_counts
        .map { sample, stage, count -> tuple(sample, [(stage): count]) }
        .groupTuple()
        .map { sample, maps -> tuple(sample, maps.inject([:]) { acc, m -> acc + m }) }

    ch_report_rows = ch_counts_wide.map { sample, counts ->
        def row = [
            sample,
            params.sample_type,
            counts.raw ?: 'NA',
            counts.trimmed ?: 'NA'
        ]
        if (params.sample_type == 'RNA') {
            row << (counts.ribo_mapped ?: 'NA')
            row << (counts.trna_mapped ?: 'NA')
            if (filter_ribo_enabled) {
                row << (counts.ribo_intersect ?: 'NA')
            }
        }
        row << (counts.host_mapped ?: 'NA')
        row << (counts.host_mapped_deduped ?: 'NA')
        if (params.sample_type == 'RNA') {
            row << (counts.multimapped ?: 'NA')
        }
        row << (counts.unmapped ?: 'NA')
        row << (counts.complexity_filtered ?: 'NA')
        row.join('\t')
    }

    def header = ['sample','sample_type','raw_reads','trimmed_reads']
    if (params.sample_type == 'RNA') {
        header << 'ribo_mapped'
        header << 'trna_mapped'
        if (filter_ribo_enabled) {
            header << 'ribo_intersect'
        }
    }
    header.addAll(['host_mapped','host_mapped_deduped'])
    if (params.sample_type == 'RNA') {
        header << 'host_multimapped'
    }
    header.addAll(['host_unmapped','host_unmapped_high_complexity'])

    ch_report_file = ch_report_rows
        .collectFile(
            name: 'read_count_report.tsv',
            newLine: true,
            seed: header.join('\t')
        )

    emit:
    host_bam = ch_host_bam
    host_bam_deduped = ch_host_bam_publish
    mapdamage = MAPDAMAGE.out.results
    read_length_distribution = READ_LENGTH_DISTRIBUTION.out.plot
    unmapped_reads = ch_unmapped_reads
    high_complexity_reads = FILTER_COMPLEXITY.out.filtered_fastq
    summary_stats = ch_report_file
    multimapped_bam = ch_multimapped_bam
    dsRNA_dir = ch_dsrna_mapdamage
    ssRNA_dir = ch_ssrna_mapdamage
    fc_genes = ch_fc_genes
    fc_exons = ch_fc_exons
    biotype_genes = ch_biotype_genes
    biotype_exons = ch_biotype_exons
}
