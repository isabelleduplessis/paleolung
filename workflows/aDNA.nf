#!/usr/bin/env nextflow

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

    if (!params.samplesheet)    { error "ERROR: --samplesheet is required" }
    if (!params.host_reference) { error "ERROR: --host_reference is required" }

    sheetRows = parseSamplesheet(params.samplesheet, params.PE)

    Channel
        .fromList(sheetRows)
        .map { row ->
            def meta = [sample: row.sample]
            def reads = params.PE ? [file(row.r1), file(row.r2)] : [file(row.r1)]
            tuple(meta, reads)
        }
        .set { ch_raw_reads }

    FASTQC_RAW(ch_raw_reads)
    TRIM_GALORE(ch_raw_reads, params.PE)
    FASTQC_TRIMMED(TRIM_GALORE.out.trimmed_reads)


    bwaIdxExts = ['amb', 'ann', 'bwt', 'pac', 'sa']
    existingIdxFiles = bwaIdxExts.collect { file("${params.host_reference}.${it}") }
    hostIdxExists = existingIdxFiles.every { it.exists() }

    if (hostIdxExists) {
        ch_host_index = Channel.value( tuple(file(params.host_reference), existingIdxFiles) )
    } else {
        INDEX_HOST_REFERENCE(file(params.host_reference))
        ch_host_index = INDEX_HOST_REFERENCE.out.index.first()
    }

    base_chunks = params.PE
        ? TRIM_GALORE.out.trimmed_reads
            // flatten [meta, [R1, R2]] -> [meta, R1, R2] so pe:true can sync-split
            .map { meta, reads -> tuple(meta, reads[0], reads[1]) }
            .splitFastq(by: params.chunk_size, pe: true, file: true)
            // normalize back to [meta, [chunkR1, chunkR2]]
            .map { meta, r1, r2 -> tuple(meta, [r1, r2]) }
        : TRIM_GALORE.out.trimmed_reads
            .splitFastq(by: params.chunk_size, file: true)
            .map { meta, r -> tuple(meta, [r]) }

    ch_chunks = base_chunks
        .map { meta, reads -> tuple(meta.id, meta, reads) }
        .groupTuple()
        .flatMap { id, metas, reads_list ->
            reads_list.indices.collect { i ->
                tuple(metas[i] + [id: "${id}_chunk${i + 1}", chunk: i + 1], reads_list[i])
            }
        }

    HOST_DNA_MAPPING(ch_chunks, ch_host_index, params.PE, params.seed_length, params.map_missing_prob, params.gap_fraction)
    ch_bams_by_sample = HOST_DNA_MAPPING.out.bam
        .map { meta, bam -> tuple(meta.sample, bam) }
        .groupTuple()
        .map { sample, bams -> tuple([sample: sample], bams) }

    MERGE_BAMS(ch_bams_by_sample, params.PE, params.map_quality)

    MARK_DUPLICATES(MERGE_BAMS.out.sorted_bam)

    MAPDAMAGE(MARK_DUPLICATES.out.dedup_bam, file(params.host_reference))

    FILTER_COMPLEXITY(MERGE_BAMS.out.unmapped_reads, params.PE)

    READ_LENGTH_DISTRIBUTION(MARK_DUPLICATES.out.dedup_bam)


    // if the user has provided a database for kraken2 and bracken, run it, else continue without it
    if (params.kraken2_db) {
        db_ch = file(params.kraken2_db, checkIfExists: true)
        TAXONOMIC_PROFILING(FILTER_COMPLEXITY.out.filtered_fastq, db_ch, params.bracken_read_len, params.bracken_tax_level)
        kraken_out_ch  = TAXONOMIC_PROFILING.out.kraken_out
        report_ch      = TAXONOMIC_PROFILING.out.report
        bracken_ch     = TAXONOMIC_PROFILING.out.bracken
        summary_ch     = TAXONOMIC_PROFILING.out.summaryfile

        KRONA(TAXONOMIC_PROFILING.out.bracken_kreport, params.krona_min_fraction)
        krona_ch = KRONA.out.html

    } else {
        kraken_out_ch  = Channel.empty()
        report_ch      = Channel.empty()
        bracken_ch     = Channel.empty()
        summary_ch     = Channel.empty()
        krona_ch = Channel.empty()
    }


    // Summary stats
    ch_raw_counts = FASTQC_RAW.out.read_count
        .map { meta, f -> tuple(meta.sample, f.text.trim().toInteger()) }
        .groupTuple()
        .map { sample, counts -> tuple(sample, counts.sum()) }

    ch_trimmed_counts = FASTQC_TRIMMED.out.read_count
        .map { meta, f -> tuple(meta.sample, f.text.trim().toInteger()) }
        .groupTuple()
        .map { sample, counts -> tuple(sample, counts.sum()) }

    ch_mapped_counts = MERGE_BAMS.out.mapped_read_count
        .map { meta, f -> tuple(meta.sample, f.text.trim().toInteger()) }

    ch_dedup_counts = MARK_DUPLICATES.out.dedup_read_count
        .map { meta, f -> tuple(meta.sample, f.text.trim().toInteger()) }

    ch_unmapped_counts = MERGE_BAMS.out.unmapped_read_count
        .map { meta, f -> tuple(meta.sample, f.text.trim().toInteger()) }

    ch_highc_counts = FILTER_COMPLEXITY.out.read_count
        .map { meta, f -> tuple(meta.sample, f.text.trim().toInteger()) }

    ch_stats_rows = ch_raw_counts
        .join(ch_trimmed_counts)
        .join(ch_mapped_counts)
        .join(ch_dedup_counts)
        .join(ch_unmapped_counts)
        .join(ch_highc_counts)
        .map { sample, raw, trimmed, mapped, dedup, unmapped, highc ->
            "${sample}\t${raw}\t${trimmed}\t${mapped}\t${dedup}\t${unmapped}\t${highc}"
        }

    ch_summary_stats = ch_stats_rows.collectFile(
        name: 'read_stats.tsv',
        newLine: true,
        sort: true,
        seed: [
            'sample',
            'raw_reads',
            'trimmed_reads',
            'host_mapped_reads',
            'host_mapped_dedup_reads',
            'unmapped_reads',
            'high_complexity_unmapped_reads'
        ].join('\t')
    )


    ch_host_bam_publish = params.publish_host_bam
        ? MARK_DUPLICATES.out.dedup_bam
        : Channel.empty()

    publish:

    fastqc_raw                = FASTQC_RAW.out.qc
    fastqc_trimmed            = FASTQC_TRIMMED.out.qc
    trimmed_reads             = TRIM_GALORE.out.trimmed_reads
    host_reference_index      = ch_host_index
    host_bam                  = ch_host_bam_publish
    mapdamage                 = MAPDAMAGE.out.results
    read_length_distribution  = READ_LENGTH_DISTRIBUTION.out.plot
    unmapped_reads            = MERGE_BAMS.out.unmapped_reads
    high_complexity_reads     = FILTER_COMPLEXITY.out.filtered_fastq
    summary_stats             = ch_summary_stats
    krona_html = krona_ch    

}


output {

    fastqc_raw {
        path { meta, qc_files -> "fastqc/raw/${meta.sample}" }
    }

    fastqc_trimmed {
        path { meta, qc_files -> "fastqc/trimmed/${meta.sample}" }
    }

    trimmed_reads {
        path { meta, reads -> "trimmed_reads/${meta.sample}" }
    }

    host_reference_index {
        path { ref, idx_files -> "host_reference_index" }
    }

    host_bam {
        path { meta, bam, bai -> "host_bam/${meta.sample}" }
    }

    mapdamage {
        path { meta, dir -> "mapdamage" }
    }

    read_length_distribution {
        path { meta, plot -> "read_length_distribution" }
    }

    unmapped_reads {
        path { meta, fastq -> "unmapped_reads" }
    }

    high_complexity_reads {
        path { meta, fastq -> "high_complexity_reads" }
    }

    krona_html {
        path { meta, plot -> "krona" }
    }

    summary_stats {
        path 'summary'
    }
}