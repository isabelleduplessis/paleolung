#!/usr/bin/env nextflow
include { FASTQC_RAW } from '../modules/fastqc_raw.nf'
include { TRIM_GALORE } from '../modules/trim.nf'
include { FASTP } from '../modules/fastp.nf'
include { FASTQC_TRIMMED } from '../modules/fastqc_trimmed.nf'

workflow PREPROCESSING {
    take:
    ch_raw_reads
    PE
    trim_length
    trim_quality
    trimmer
    merge_PE
    min_overlap

    main:
    FASTQC_RAW(ch_raw_reads)

    if (trimmer == 'fastp') {
        FASTP(ch_raw_reads, PE, trim_length, trim_quality, merge_PE, min_overlap)
        if (merge_PE) {
            ch_trimmed_reads = FASTP.out.merged_reads
        } else {
            ch_trimmed_reads = FASTP.out.trimmed_reads
        }
    } else {
        TRIM_GALORE(ch_raw_reads, PE, trim_length, trim_quality)
        ch_trimmed_reads = TRIM_GALORE.out.trimmed_reads
    }

    FASTQC_TRIMMED(ch_trimmed_reads)

    emit:
    fastqc_raw = FASTQC_RAW.out.qc
    fastqc_trimmed = FASTQC_TRIMMED.out.qc
    trimmed_reads = ch_trimmed_reads
}
