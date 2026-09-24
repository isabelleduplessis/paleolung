#!/usr/bin/env nextflow
include { MAP_RIBO } from '../modules/RNA_map_ribosomal.nf'
include { STAR_MAPPING } from '../modules/RNA_star_mapping.nf'
include { RNA_BT2_MAPPING } from '../modules/RNA_bt2_mapping.nf'
include { RNA_BWA_MAPPING } from '../modules/RNA_bwa_mapping.nf'
include { FILTER_RIBO } from '../modules/filter_ribo.nf'

workflow ANCIENT_RNA_HOST_MAPPING {

    take:
    ch_trimmed_reads
    RNA_host_mapper

    main:

    // prepare ribo reference files
    ribo_bt2_index_ch = Channel.fromPath("${params.ribo_bt2_ref}*.bt2*").collect()
    ribo_bt2_prefix   = file(params.ribo_bt2_ref).name

    // map to ribosomal
    MAP_RIBO(ch_trimmed_reads, ribo_bt2_index_ch, ribo_bt2_prefix, params.PE)

    // map to host - can test different mappers
    if (RNA_host_mapper == 'STAR') {
        hg_star_index_dir = Channel.fromPath(params.hg_star_ref, type: 'dir').first()
        STAR_MAPPING(MAP_RIBO.out.unmapped_reads, hg_star_index_dir, params.PE, params.STAR_mapq)
        ch_host_mapped        = STAR_MAPPING.out.sorted_bam
        ch_host_unmapped      = STAR_MAPPING.out.unmapped_reads
        ch_multimapped_reads  = STAR_MAPPING.out.multimapped_reads
        ch_multimapped_bam = STAR_MAPPING.out.multimapped_bam
        
    }
    else if (RNA_host_mapper == 'bowtie2') {
        bt2_prefix            = file(params.bt2_ref).name
        hg_bowtie2_index_dir  = Channel.fromPath("${params.bt2_ref}*.bt2*").collect()
        RNA_BT2_MAPPING(MAP_RIBO.out.unmapped_reads, hg_bowtie2_index_dir, bt2_prefix, params.PE, params.bt2_mapq)
        ch_host_mapped        = RNA_BT2_MAPPING.out.sorted_bam
        ch_host_unmapped      = RNA_BT2_MAPPING.out.unmapped_reads
        ch_multimapped_reads  = RNA_BT2_MAPPING.out.multimapped_reads
        ch_multimapped_bam = RNA_BT2_MAPPING.out.multimapped_bam
    }
    else if (RNA_host_mapper == 'BWA') {
        def bwa_index_files = files("${params.host_reference}.*")
        hg_bwa_index_ch = Channel.value( tuple(file(params.host_reference), bwa_index_files) )
        RNA_BWA_MAPPING(MAP_RIBO.out.unmapped_reads, hg_bwa_index_ch, params.PE, params.seed_length, params.map_missing_prob, params.gap_fraction, params.bwa_mapq)
        ch_host_mapped        = RNA_BWA_MAPPING.out.sorted_bam
        ch_host_unmapped      = RNA_BWA_MAPPING.out.unmapped_reads
        ch_multimapped_reads  = RNA_BWA_MAPPING.out.multimapped_reads
        ch_multimapped_bam = Channel.empty()
    }
    else {
        error "ERROR: --RNA_host_mapper must be one of 'STAR', 'bowtie2', 'BWA' (got: ${params.RNA_host_mapper})"
    }

    // filter out ribosomal again
    FILTER_RIBO(ch_host_mapped, file(params.rrna_bed))
    
    emit:
    ribo_mapped_bam = MAP_RIBO.out.sorted_bam
    host_bam_with_dups = FILTER_RIBO.out.filtered_bam
    ribo_intersect_bam = FILTER_RIBO.out.ribo_bam
    unmapped_reads            = ch_host_unmapped
    multimapped_reads = ch_multimapped_reads
    multimapped_bam = ch_multimapped_bam
}
