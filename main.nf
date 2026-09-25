#!/user/bin/env nextflow
include { FASTQC_RAW } from './modules/fastqc_raw.nf'
include { TRIM_GALORE } from './modules/trim.nf'
include { FASTQC_TRIMMED } from './modules/fastqc_trimmed.nf'
include { ANCIENT_DNA_HOST_MAPPING } from './workflows/aDNA'
include { ANCIENT_RNA_HOST_MAPPING } from './workflows/aRNA'
include { MARK_DUPLICATES } from './modules/mark_duplicates.nf'
include { MAPDAMAGE } from './modules/mapdamage.nf'
include { FILTER_COMPLEXITY } from './modules/filter_read_complexity.nf'
include { READ_LENGTH_DISTRIBUTION } from './modules/read_length_plot.nf'
include { COUNT_FASTQ_READS } from './modules/count_fastq.nf'
include { COUNT_BAM_READS } from './modules/count_bam.nf'
include { DSRNA_INTERSECT } from './modules/dsrna.nf'
include { FEATURE_COUNTS } from './modules/featurecounts.nf'
include { FEATURECOUNTS_SUMMARY as FEATURECOUNTS_SUMMARY_EXON } from './modules/featurecounts.nf'
include { FEATURECOUNTS_SUMMARY as FEATURECOUNTS_SUMMARY_GENE } from './modules/featurecounts.nf'
include { BIOTYPE_COMPOSITION_PLOT as BIOTYPE_COMPOSITION_PLOT_EXON } from './modules/biotype_fc_plot.nf'
include { BIOTYPE_COMPOSITION_PLOT as BIOTYPE_COMPOSITION_PLOT_GENE } from './modules/biotype_fc_plot.nf'



// eventually make it so that you can skip some of the mapping processes if user provides fastq files for tax porfiling for example or contig assembly

def parseSamplesheet(String path, boolean isPE) {

    def sheetFile = file(path)
    if (!sheetFile.exists()) {
        error "ERROR: samplesheet not found at ${path}"
    }

    def lines = sheetFile.readLines().findAll { it.trim() }
    if (lines.isEmpty()) {
        error "ERROR: samplesheet ${path} is empty"
    }

    def header = lines[0].split(',', -1).collect { it.trim() }
    def dataLines = lines.drop(1)

    def sampleIdx = header.indexOf('sample')
    def r1Idx     = header.indexOf('R1')
    def r2Idx     = header.indexOf('R2')
    def hasR2Col  = (r2Idx != -1)

    if (sampleIdx == -1 || r1Idx == -1) {
        error "ERROR: samplesheet must contain at minimum 'sample' and 'R1' columns"
    }

    if (isPE && !hasR2Col) {
        error "ERROR: --PE was specified but the samplesheet has no 'R2' column"
    }
    if (!isPE && hasR2Col) {
        error "ERROR: samplesheet has an 'R2' column but --PE was not specified. " +
              "Either remove the column or re-run with --PE."
    }

    def parsedRows = []

    dataLines.eachWithIndex { line, i ->
        def rowNum = i + 2
        def fields = line.split(',', -1).collect { it.trim() }

        def sampleName = fields.size() > sampleIdx ? fields[sampleIdx] : ''
        def r1Path     = fields.size() > r1Idx ? fields[r1Idx] : ''
        def r2Path     = hasR2Col && fields.size() > r2Idx ? fields[r2Idx] : ''

        if (!sampleName) { error "ERROR: row ${rowNum} of samplesheet is missing a 'sample' value" }
        if (!r1Path)     { error "ERROR: row ${rowNum} of samplesheet is missing an 'R1' value" }

        // uniform PE/SE across every row
        if (isPE && !r2Path) {
            error "ERROR: row ${rowNum} — --PE was specified but 'R2' is blank for this row. " +
                  "All rows must be paired-end when --PE is set."
        }
        if (!isPE && r2Path) {
            error "ERROR: row ${rowNum} — 'R2' has a value but --PE was not specified."
        }

        parsedRows << [
            sample : sampleName,
            r1     : r1Path,
            r2     : r2Path
        ]
    }

    def countsBySample = parsedRows.countBy { it.sample }
    def dupSamples = countsBySample.findAll { k, v -> v > 1 }.keySet()
    if (dupSamples && !params.merge_inputs) {
        error "ERROR: samplesheet has duplicate 'sample' names (${dupSamples.join(', ')}) " +
              "but --merge_inputs was not set. Pass --merge_inputs to treat same-named rows " +
              "as separate inputs to be merged after mapping, or fix the samplesheet."
    }

    return parsedRows
}


workflow {
	main:

	sheetRows = parseSamplesheet(params.samplesheet, params.PE)

	Channel
        .fromList(sheetRows)
        .map { row ->
            def meta = [sample: row.sample]
            def reads = params.PE ? [file(row.r1), file(row.r2)] : [file(row.r1)]
            tuple(meta, reads)
        }
        .set { ch_raw_reads }

    // preprocess reads
    FASTQC_RAW(ch_raw_reads)
    TRIM_GALORE(ch_raw_reads, params.PE, params.trim_length, params.trim_quality)
    FASTQC_TRIMMED(TRIM_GALORE.out.trimmed_reads)

    // perform mapping depending on sample type
    if (params.sample_type == 'DNA') {
        ANCIENT_DNA_HOST_MAPPING(TRIM_GALORE.out.trimmed_reads)
        ch_unmapped_reads    = ANCIENT_DNA_HOST_MAPPING.out.unmapped_reads
        ch_host_bam = ANCIENT_DNA_HOST_MAPPING.out.host_bam_with_dups
        ch_ribo_mapped    = Channel.empty()
        ch_ribo_intersect = Channel.empty()
        ch_multimapped_reads = Channel.empty()
        ch_multimapped_bam = Channel.empty()

    } else {
        ANCIENT_RNA_HOST_MAPPING(TRIM_GALORE.out.trimmed_reads, params.RNA_host_mapper)
        ch_unmapped_reads    = ANCIENT_RNA_HOST_MAPPING.out.unmapped_reads
        ch_host_bam = ANCIENT_RNA_HOST_MAPPING.out.host_bam_with_dups
        ch_ribo_mapped = ANCIENT_RNA_HOST_MAPPING.out.ribo_mapped_bam
        ch_ribo_intersect = ANCIENT_RNA_HOST_MAPPING.out.ribo_intersect_bam
        ch_multimapped_reads = ANCIENT_RNA_HOST_MAPPING.out.multimapped_reads
        ch_multimapped_bam = ANCIENT_RNA_HOST_MAPPING.out.multimapped_bam
    }


    // host mapping post processing
    MARK_DUPLICATES(ch_host_bam)
    MAPDAMAGE(MARK_DUPLICATES.out.dedup_bam, file(params.host_reference))
    READ_LENGTH_DISTRIBUTION(MARK_DUPLICATES.out.dedup_bam)
    FILTER_COMPLEXITY(ch_unmapped_reads, params.PE)
    ch_host_bam_publish = params.publish_host_bam
        ? MARK_DUPLICATES.out.dedup_bam
        : Channel.empty()

    if (params.sample_type == 'RNA') {
        // RNA specific post processing
        DSRNA_INTERSECT(MARK_DUPLICATES.out.dedup_bam, file(params.host_reference), file(params.dsRNA_bed), params.strand)
        FEATURE_COUNTS(MARK_DUPLICATES.out.dedup_bam, file(params.hg_gtf), params.strand)

        // Collect exon counts + summaries across all samples
        exon_files = FEATURE_COUNTS.out.exon_counts
            .map { meta, txt -> txt }
            .mix(FEATURE_COUNTS.out.exon_summary.map { meta, summ -> summ })
            .collect()

        // Collect gene counts + summaries across all samples
        gene_files = FEATURE_COUNTS.out.gene_counts
            .map { meta, txt -> txt }
            .mix(FEATURE_COUNTS.out.gene_summary.map { meta, summ -> summ })
            .collect()

        FEATURECOUNTS_SUMMARY_EXON(exon_files, file(params.hg_gtf), "exon")
        FEATURECOUNTS_SUMMARY_GENE(gene_files, file(params.hg_gtf), "gene")

        BIOTYPE_COMPOSITION_PLOT_EXON(FEATURECOUNTS_SUMMARY_EXON.out.csv, "exon")
        BIOTYPE_COMPOSITION_PLOT_GENE(FEATURECOUNTS_SUMMARY_GENE.out.csv, "gene")
    }

    // GENERATING COUNT REPORT ACROSS ALL SAMPLES

    // get fastq read counts
    ch_to_count = ch_raw_reads.map { meta, f -> tuple(meta + [stage:'raw'], f) }
        .mix( TRIM_GALORE.out.trimmed_reads.map { meta, f -> tuple(meta + [stage:'trimmed'], f) } )
        .mix( ch_unmapped_reads.map { meta, f -> tuple(meta + [stage:'unmapped'], f) } )
        .mix( FILTER_COMPLEXITY.out.filtered_fastq.map { meta, f -> tuple(meta + [stage:'complexity_filtered'], f) } )
        .mix( ch_multimapped_reads.map { meta, f -> tuple(meta + [stage:'multimapped'], f) } )

    COUNT_FASTQ_READS(ch_to_count, params.PE)
    ch_all_counts = COUNT_FASTQ_READS.out.count
        .map { meta, count -> tuple(meta.sample, meta.stage, count as Integer) }


    // get bam read counts
    ch_bam_to_count = ch_host_bam.map { meta, bam -> tuple(meta + [stage:'host_mapped'], bam) }
        .mix( MARK_DUPLICATES.out.dedup_bam.map { meta, bam, bai -> tuple(meta + [stage:'host_mapped_deduped'], bam) } )
        .mix( ch_ribo_mapped.map { meta, bam -> tuple(meta + [stage:'ribo_mapped'], bam) } )
        .mix( ch_ribo_intersect.map { meta, bam -> tuple(meta + [stage:'ribo_intersect'], bam) } )

    COUNT_BAM_READS(ch_bam_to_count)

    ch_all_bam_counts = COUNT_BAM_READS.out.count
        .map { meta, count -> tuple(meta.sample, meta.stage, count as Integer) }


    // combine counts into stats

    ch_combined_counts = ch_all_counts.mix(ch_all_bam_counts)

    ch_counts_wide = ch_combined_counts
        .map { sample, stage, count -> tuple(sample, [(stage): count]) }
        .groupTuple()
        .map { sample, maps -> tuple(sample, maps.inject([:]) { acc, m -> acc + m }) }

    ch_report_rows = ch_counts_wide.map { sample, counts ->
        def row = [
            sample,
            params.sample_type,
            counts.raw     ?: 'NA',
            counts.trimmed ?: 'NA'
        ]
        if (params.sample_type == 'RNA') {
            row << (counts.ribo_mapped    ?: 'NA')
            row << (counts.ribo_intersect ?: 'NA')
        }
        row << (counts.host_mapped         ?: 'NA')
        row << (counts.host_mapped_deduped ?: 'NA')
        if (params.sample_type == 'RNA') {
            row << (counts.multimapped ?: 'NA')
        }
        row << (counts.unmapped            ?: 'NA')
        row << (counts.complexity_filtered ?: 'NA')
        row.join('\t')
    }

    def header = params.sample_type == 'RNA' ?
        ['sample','sample_type','raw_reads','trimmed_reads',
         'ribo_mapped','ribo_intersect',
         'host_mapped','host_mapped_deduped',
         'host_multimapped','host_unmapped','host_unmapped_high_complexity'] :
        ['sample','sample_type','raw_reads','trimmed_reads',
         'host_mapped','host_mapped_deduped',
         'host_unmapped','host_unmapped_high_complexity']

    ch_report_file = ch_report_rows
        .collectFile(
            name: 'read_count_report.tsv',
            newLine: true,
            seed: header.join('\t')
        )


    publish:
    fastqc_raw = FASTQC_RAW.out.qc
    fastqc_trimmed = FASTQC_TRIMMED.out.qc
    trimmed_reads = TRIM_GALORE.out.trimmed_reads
    host_bam = ch_host_bam
    host_bam_deduped = ch_host_bam_publish
    mapdamage = MAPDAMAGE.out.results
    read_length_distribution = READ_LENGTH_DISTRIBUTION.out.plot
    high_complexity_reads = FILTER_COMPLEXITY.out.filtered_fastq
    summary_stats = ch_report_file
    multimapped_bam = ch_multimapped_bam
    dsRNA_dir = DSRNA_INTERSECT.out.dsRNA_mapdamage
    ssRNA_dir = DSRNA_INTERSECT.out.ssRNA_mapdamage
    fc_genes = FEATURECOUNTS_SUMMARY_EXON.out.csv
    fc_exons = FEATURECOUNTS_SUMMARY_GENE.out.csv
    biotype_genes = BIOTYPE_COMPOSITION_PLOT_GENE.out.pdf
    biotype_exons = BIOTYPE_COMPOSITION_PLOT_EXON.out.pdf

}

output {

    fastqc_raw {
        path { meta, qc_files -> "reads/fastqc/raw/${meta.sample}" }
    }

    fastqc_trimmed {
        path { meta, qc_files -> "reads/fastqc/trimmed/${meta.sample}" }
    }

    trimmed_reads {
        path { meta, reads -> "reads/trimmed/${meta.sample}" }
    }

    host_bam {
        path { meta, bam -> "host_mapping/${params.RNA_host_mapper}/${meta.sample}" }
    }

    host_bam_deduped {
        path { meta, bam, bai -> "host_mapping/${params.RNA_host_mapper}/${meta.sample}" }
    }

    multimapped_bam {
        path { meta, bam -> "multimapped_bam/${meta.sample}" }
    }

    mapdamage {
        path { meta, dir -> "damage_patterns/${params.RNA_host_mapper}/host_mapped" }
    }

    dsRNA_dir {
        path { meta, dir -> "damage_patterns/${params.RNA_host_mapper}/dsRNA" }
    }

    ssRNA_dir {
        path { meta, dir -> "damage_patterns/${params.RNA_host_mapper}/ssRNA" }
    }

    read_length_distribution {
        path { meta, plot -> "plots/${params.RNA_host_mapper}/read_length_distribution" }
    }

    high_complexity_reads {
        path { meta, fastq -> "reads/high_complexity_reads" }
    }

    summary_stats {
        path "stats/${params.RNA_host_mapper}"
    }

    fc_genes {
        path "stats/${params.RNA_host_mapper}"
    }

    fc_exons {
        path "stats/${params.RNA_host_mapper}"
    }

    biotype_genes {
        path "plots/${params.RNA_host_mapper}/biotype_composition"
    }

    biotype_exons {
        path "plots/${params.RNA_host_mapper}/biotype_composition"
    }
}