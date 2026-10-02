#!/user/bin/env nextflow
include { PREPROCESSING } from './workflows/PREPROCESSING.nf'
include { HOST_MAPPING } from './workflows/HOST_MAPPING.nf'
include { TAXONOMIC_PROFILING } from './workflows/TAXONOMIC_PROFILING.nf'
include { CONTIG_ASSEMBLY } from './workflows/CONTIG_ASSEMBLY.nf'

def parseSamplesheet(String path, boolean pairedInput) {

    def sheetFile = file(path)
    if (!sheetFile.exists()) {
        error "ERROR: samplesheet not found at ${path}"
    }

    def lines = sheetFile.readLines().findAll { line -> line.trim() }
    if (lines.isEmpty()) {
        error "ERROR: samplesheet ${path} is empty"
    }

    def header = lines[0].split(',', -1).collect { col -> col.trim() }
    def dataLines = lines.drop(1)

    def sampleIdx = header.indexOf('sample')
    def r1Idx     = header.indexOf('R1')
    def r2Idx     = header.indexOf('R2')
    def hasR2Col  = (r2Idx != -1)

    if (sampleIdx == -1 || r1Idx == -1) {
        error "ERROR: samplesheet must contain at minimum 'sample' and 'R1' columns"
    }

    if (pairedInput && !hasR2Col) {
        error "ERROR: --read_type '${params.read_type}' requires an 'R2' column in the samplesheet"
    }
    if (!pairedInput && hasR2Col) {
        error "ERROR: samplesheet has an 'R2' column but --read_type '${params.read_type ?: 'se'}' is single-end. " +
              "Set --read_type pe or --read_type merge_pe."
    }

    def parsedRows = []

    dataLines.eachWithIndex { line, i ->
        def rowNum = i + 2
        def fields = line.split(',', -1).collect { field -> field.trim() }

        def sampleName = fields.size() > sampleIdx ? fields[sampleIdx] : ''
        def r1Path     = fields.size() > r1Idx ? fields[r1Idx] : ''
        def r2Path     = hasR2Col && fields.size() > r2Idx ? fields[r2Idx] : ''

        if (!sampleName) { error "ERROR: row ${rowNum} of samplesheet is missing a 'sample' value" }
        if (!r1Path)     { error "ERROR: row ${rowNum} of samplesheet is missing an 'R1' value" }

        // All rows must use the same paired-end or single-end layout.
        if (pairedInput && !r2Path) {
            error "ERROR: row ${rowNum} — paired-end read_type requires an 'R2' value. " +
                  "All rows must contain both R1 and R2."
        }
        if (!pairedInput && r2Path) {
            error "ERROR: row ${rowNum} — 'R2' has a value but read_type is single-end."
        }

        parsedRows << [
            sample : sampleName,
            r1     : r1Path,
            r2     : r2Path
        ]
    }

    def countsBySample = parsedRows.countBy { row -> row.sample }
    def dupSamples = countsBySample.findAll { _k, v -> v > 1 }.keySet()
    if (dupSamples && !params.merge_inputs) {
        error "ERROR: samplesheet has duplicate 'sample' names (${dupSamples.join(', ')}) " +
              "but --merge_inputs was not set. Pass --merge_inputs to treat same-named rows " +
              "as separate inputs to be merged after mapping, or fix the samplesheet."
    }

    return parsedRows
}


workflow {
    main:

    // warnings based on parameters
    def trimmer = params.trimmer?.toString()?.toLowerCase() ?: 'trim_galore'
    if (!(trimmer in ['fastp', 'trim_galore'])) {
        error "ERROR: --trimmer must be 'fastp' or 'trim_galore', got '${params.trimmer}'"
    }
    def read_type = params.read_type?.toString()?.toLowerCase() ?: 'se'
    if (!(read_type in ['se', 'pe', 'merge_pe'])) {
        error "ERROR: --read_type must be 'se', 'pe', or 'merge_pe', got '${params.read_type}'"
    }
    def input_PE = read_type in ['pe', 'merge_pe']
    def merge_PE = read_type == 'merge_pe'
    if (merge_PE && trimmer != 'fastp') {
        error "ERROR: --read_type merge_pe requires --trimmer fastp (got '${trimmer}')"
    }
    if (merge_PE && params.skip_host_mapping) {
        error 'ERROR: --read_type merge_pe cannot be used with --skip_host_mapping because preprocessing is skipped'
    }
    def min_overlap = params.min_overlap == null ? 30 : params.min_overlap as Integer
    if (min_overlap < 1) {
        error "ERROR: --min_overlap must be a positive integer, got '${params.min_overlap}'"
    }
    def downstream_PE = read_type == 'pe'

    def assembly_tool = params.assembly_tool?.toString()?.toLowerCase()
    if (assembly_tool && !(assembly_tool in ['spades', 'megahit'])) {
        error "ERROR: --assembly_tool must be 'spades' or 'megahit', got '${params.assembly_tool}'"
    }

    sheetRows = parseSamplesheet(params.samplesheet, input_PE)
    inputs_per_sample = sheetRows.countBy { it.sample }

    channel
        .fromList(sheetRows)
        .map { row ->
            def meta = [sample: row.sample, num_inputs: inputs_per_sample[row.sample]]
            def reads = input_PE ? [file(row.r1), file(row.r2)] : [file(row.r1)]
            tuple(meta, reads)
        }
        .set { ch_raw_reads }

    if (params.skip_host_mapping) {
        // In skip_host_mapping mode, R1/R2 already point to the high-complexity reads.
        // Skip trimming and host mapping and profile these files directly.
        TAXONOMIC_PROFILING(ch_raw_reads, downstream_PE)

        ch_fastqc_raw = channel.empty()
        ch_fastqc_trimmed = channel.empty()
        ch_trimmed_reads = channel.empty()
        ch_host_bam = channel.empty()
        ch_host_bam_deduped = channel.empty()
        ch_mapdamage = channel.empty()
        ch_read_length_distribution = channel.empty()
        ch_high_complexity_reads = ch_raw_reads
        ch_summary_stats = channel.empty()
        ch_multimapped_bam = channel.empty()
        ch_dsRNA_dir = channel.empty()
        ch_ssRNA_dir = channel.empty()
        ch_fc_genes = channel.empty()
        ch_fc_exons = channel.empty()
        ch_biotype_genes = channel.empty()
        ch_biotype_exons = channel.empty()
    } else {
        PREPROCESSING(ch_raw_reads, input_PE, params.trim_length, params.trim_quality, trimmer, merge_PE, min_overlap)
        HOST_MAPPING(ch_raw_reads, PREPROCESSING.out.trimmed_reads, downstream_PE)
        TAXONOMIC_PROFILING(HOST_MAPPING.out.high_complexity_reads, downstream_PE)

        ch_fastqc_raw = PREPROCESSING.out.fastqc_raw
        ch_fastqc_trimmed = PREPROCESSING.out.fastqc_trimmed
        ch_trimmed_reads = PREPROCESSING.out.trimmed_reads
        ch_host_bam = HOST_MAPPING.out.host_bam
        ch_host_bam_deduped = HOST_MAPPING.out.host_bam_deduped
        ch_mapdamage = HOST_MAPPING.out.mapdamage
        ch_read_length_distribution = HOST_MAPPING.out.read_length_distribution
        ch_high_complexity_reads = HOST_MAPPING.out.high_complexity_reads
        ch_summary_stats = HOST_MAPPING.out.summary_stats
        ch_multimapped_bam = HOST_MAPPING.out.multimapped_bam
        ch_dsRNA_dir = HOST_MAPPING.out.dsRNA_dir
        ch_ssRNA_dir = HOST_MAPPING.out.ssRNA_dir
        ch_fc_genes = HOST_MAPPING.out.fc_genes
        ch_fc_exons = HOST_MAPPING.out.fc_exons
        ch_biotype_genes = HOST_MAPPING.out.biotype_genes
        ch_biotype_exons = HOST_MAPPING.out.biotype_exons
    }

    ch_contig_assembly = channel.empty()
    if (assembly_tool) {
        CONTIG_ASSEMBLY(ch_high_complexity_reads, downstream_PE)
        ch_contig_assembly = CONTIG_ASSEMBLY.out.assembly
    }

    publish:
    fastqc_raw = ch_fastqc_raw
    fastqc_trimmed = ch_fastqc_trimmed
    trimmed_reads = ch_trimmed_reads
    host_bam = ch_host_bam
    host_bam_deduped = ch_host_bam_deduped
    mapdamage = ch_mapdamage
    read_length_distribution = ch_read_length_distribution
    high_complexity_reads = ch_high_complexity_reads
    summary_stats = ch_summary_stats
    multimapped_bam = ch_multimapped_bam
    dsRNA_dir = ch_dsRNA_dir
    ssRNA_dir = ch_ssRNA_dir
    fc_genes = ch_fc_genes
    fc_exons = ch_fc_exons
    biotype_genes = ch_biotype_genes
    biotype_exons = ch_biotype_exons
    contig_assembly = ch_contig_assembly
    bracken_summary = TAXONOMIC_PROFILING.out.bracken_summary
    krona_html = TAXONOMIC_PROFILING.out.krona_html
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
        path { meta, bam -> params.sample_type == 'RNA'
            ? "host_mapping/${params.RNA_host_mapper}/${meta.sample}"
            : "host_mapping/${meta.sample}" }
    }

    host_bam_deduped {
        path { meta, bam, bai -> params.sample_type == 'RNA'
            ? "host_mapping/${params.RNA_host_mapper}/${meta.sample}"
            : "host_mapping/${meta.sample}" }
    }

    multimapped_bam {
        path { meta, bam -> "multimapped_bam/${meta.sample}" }
    }

    mapdamage {
        path { meta, dir -> params.sample_type == 'RNA'
            ? "damage_patterns/${params.RNA_host_mapper}/host_mapped"
            : "damage_patterns/host_mapped" }
    }

    dsRNA_dir {
        path { meta, dir -> params.sample_type == 'RNA'
            ? "damage_patterns/${params.RNA_host_mapper}/dsRNA"
            : "damage_patterns/dsRNA" }
    }

    ssRNA_dir {
        path { meta, dir -> params.sample_type == 'RNA'
            ? "damage_patterns/${params.RNA_host_mapper}/ssRNA"
            : "damage_patterns/ssRNA" }
    }

    read_length_distribution {
        path { meta, plot -> params.sample_type == 'RNA'
            ? "plots/${params.RNA_host_mapper}/read_length_distribution"
            : "plots/read_length_distribution" }
    }

    high_complexity_reads {
        path { meta, fastq -> "reads/high_complexity_reads" }
    }

    summary_stats {
        path "stats${params.sample_type == 'RNA' ? '/' + params.RNA_host_mapper : ''}"
    }

    fc_genes {
        path "stats${params.sample_type == 'RNA' ? '/' + params.RNA_host_mapper : ''}"
    }

    fc_exons {
        path "stats${params.sample_type == 'RNA' ? '/' + params.RNA_host_mapper : ''}"
    }

    biotype_genes {
        path "plots${params.sample_type == 'RNA' ? '/' + params.RNA_host_mapper : ''}/biotype_composition"
    }

    biotype_exons {
        path "plots${params.sample_type == 'RNA' ? '/' + params.RNA_host_mapper : ''}/biotype_composition"
    }

    contig_assembly {
        path { meta, assembly -> "assemblies/${params.assembly_tool?.toString()?.toLowerCase()}/${meta.sample}" }
    }

    bracken_summary {
        path 'taxonomic_profiling'
    }

    krona_html {
        path { meta, html -> "taxonomic_profiling/krona/${meta.sample}.html" }
    }
}