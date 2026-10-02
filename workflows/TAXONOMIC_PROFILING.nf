#!/usr/bin/env nextflow
include { KRAKEN2_PROFILING } from '../modules/kraken2_profiling.nf'
include { KRONA } from '../modules/krona.nf'
include { BRACKEN_SUMMARY } from '../modules/bracken_summary.nf'

workflow TAXONOMIC_PROFILING {
    take:
    ch_high_complexity_reads
    PE

    main:
    if (!params.kraken2_db) {
        error 'ERROR: --kraken2_db must point to a Kraken2 database'
    }

    KRAKEN2_PROFILING(
        ch_high_complexity_reads,
        file(params.kraken2_db),
        PE,
        params.bracken_read_len,
        params.bracken_tax_level
    )

    KRONA(
        KRAKEN2_PROFILING.out.bracken_kreport,
        params.krona_min_fraction
    )

    ch_bracken_summaries = KRAKEN2_PROFILING.out.bracken
        .map { _meta, summary_file -> summary_file }
        .collect()
    BRACKEN_SUMMARY(ch_bracken_summaries)

    emit:
    bracken_summary = BRACKEN_SUMMARY.out.summary
    krona_html = KRONA.out.html
}
