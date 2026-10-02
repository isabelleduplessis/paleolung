#!/usr/bin/env nextflow
include { RNA_BT2_NCRNA_MAPPING as RNA_BT2_RIBO } from '../modules/RNA_bt2_ncrna_mapping.nf'
include { RNA_BT2_NCRNA_MAPPING as RNA_BT2_TRNA } from '../modules/RNA_bt2_ncrna_mapping.nf'
include { RNA_BWA_NCRNA_MAPPING as RNA_BWA_RIBO } from '../modules/RNA_bwa_ncrna_mapping.nf'
include { RNA_BWA_NCRNA_MAPPING as RNA_BWA_TRNA } from '../modules/RNA_bwa_ncrna_mapping.nf'
include { COUNT_FASTQ_READS } from '../modules/count_fastq.nf'
include { STAR_MAPPING } from '../modules/RNA_star_host_mapping.nf'
include { RNA_BT2_MAPPING } from '../modules/RNA_bt2_host_mapping.nf'
include { RNA_BWA_MAPPING } from '../modules/RNA_bwa_host_mapping.nf'
include { MERGE_BAMS } from '../modules/merge_bams.nf'
include { FILTER_RIBO } from '../modules/filter_ribo.nf'

workflow ANCIENT_RNA_HOST_MAPPING {

    take:
    ch_trimmed_reads
    RNA_host_mapper
    ncRNA_mapper
    PE

    main:

    def ncRNA_mapper_selected = ncRNA_mapper?.toString() ?: 'bowtie2'
    if (!(ncRNA_mapper_selected in ['bowtie2', 'BWA'])) {
        error "ERROR: --ncRNA_mapper must be one of 'bowtie2' or 'BWA' (got: ${ncRNA_mapper_selected})"
    }

    def filter_ribo_value = params.filter_ribo == null ? 'true' : params.filter_ribo.toString().toLowerCase()
    if (!(filter_ribo_value in ['true', 'false'])) {
        error "ERROR: --filter_ribo must be 'true' or 'false' (got: ${params.filter_ribo})"
    }
    def filter_ribo_enabled = filter_ribo_value == 'true'

    if (ncRNA_mapper_selected == 'bowtie2') {
        // map first to ribosomal RNA, then map the remainder to tRNA
        ribo_bt2_index_ch = channel.fromPath("${params.ribo_bt2_ref}*.bt2*").collect()
        ribo_bt2_prefix   = file(params.ribo_bt2_ref).name
        RNA_BT2_RIBO(ch_trimmed_reads, ribo_bt2_index_ch, ribo_bt2_prefix, 'ribo', PE)

        trna_bt2_index_ch = channel.fromPath("${params.trna_bt2_ref}*.bt2*").collect()
        trna_bt2_prefix   = file(params.trna_bt2_ref).name
        RNA_BT2_TRNA(RNA_BT2_RIBO.out.unmapped_reads, trna_bt2_index_ch, trna_bt2_prefix, 'trna', PE)

        ch_ribo_mapped_bam = RNA_BT2_RIBO.out.sorted_bam
        ch_trna_mapped_bam = RNA_BT2_TRNA.out.sorted_bam
        ch_trna_unmapped_reads = RNA_BT2_TRNA.out.unmapped_reads
    } else {
        // BWA indexes are supplied as prebuilt reference prefixes and sidecar files
        if (!params.ribo_bwa_ref || !params.trna_bwa_ref) {
            error "ERROR: --ncRNA_mapper BWA requires both --ribo_bwa_ref and --trna_bwa_ref"
        }
        ribo_bwa_index_ch = channel.value(tuple(file(params.ribo_bwa_ref), files("${params.ribo_bwa_ref}.*")))
        RNA_BWA_RIBO(ch_trimmed_reads, ribo_bwa_index_ch, 'ribo', PE, params.seed_length, params.map_missing_prob, params.gap_fraction)

        trna_bwa_index_ch = channel.value(tuple(file(params.trna_bwa_ref), files("${params.trna_bwa_ref}.*")))
        RNA_BWA_TRNA(RNA_BWA_RIBO.out.unmapped_reads, trna_bwa_index_ch, 'trna', PE, params.seed_length, params.map_missing_prob, params.gap_fraction)

        ch_ribo_mapped_bam = RNA_BWA_RIBO.out.sorted_bam
        ch_trna_mapped_bam = RNA_BWA_TRNA.out.sorted_bam
        ch_trna_unmapped_reads = RNA_BWA_TRNA.out.unmapped_reads
    }

    // map to host - can test different mappers
    if (RNA_host_mapper == 'STAR') {
        hg_star_index_dir = channel.fromPath(params.hg_star_ref, type: 'dir').first()
        STAR_MAPPING(ch_trna_unmapped_reads, hg_star_index_dir, PE, params.max_hits)

        ch_bams_by_sample = STAR_MAPPING.out.sorted_bam
            .map { meta, bam -> tuple(groupKey(meta.sample, (meta.num_inputs ?: 1) as int), bam) }
            .groupTuple()
            .map { sample_key, bams -> tuple([sample: sample_key.getGroupTarget()], bams) }
    }
    else if (RNA_host_mapper == 'bowtie2') {
        bt2_prefix            = file(params.bt2_ref).name
        hg_bowtie2_index_dir  = channel.fromPath("${params.bt2_ref}*.bt2*").collect()
        RNA_BT2_MAPPING(ch_trna_unmapped_reads, hg_bowtie2_index_dir, bt2_prefix, PE, params.max_hits)

        ch_bams_by_sample = RNA_BT2_MAPPING.out.sorted_bam
            .map { meta, bam -> tuple(groupKey(meta.sample, (meta.num_inputs ?: 1) as int), bam) }
            .groupTuple()
            .map { sample_key, bams -> tuple([sample: sample_key.getGroupTarget()], bams) }
    }
    else if (RNA_host_mapper == 'BWA') {
        def bwa_index_files = files("${params.host_reference}.*")
        hg_bwa_index_ch = channel.value( tuple(file(params.host_reference), bwa_index_files) )

        // Count the reads remaining after ncRNA depletion, since this is the
        // input that will be split and mapped to the host in chunks.
        COUNT_FASTQ_READS(ch_trna_unmapped_reads, PE)
        ch_expected_chunks = COUNT_FASTQ_READS.out.count
            .map { meta, count -> tuple(groupKey(meta.sample, (meta.num_inputs ?: 1) as int), count as long) }
            .groupTuple()
            .map { sample_key, counts ->
                def total_reads = Math.floorDiv(counts.sum() as long, (PE ? 2L : 1L))
                def reads_per_chunk = params.chunk_size as long
                def chunk_count = (total_reads + reads_per_chunk - 1L).intdiv(reads_per_chunk)
                tuple(sample_key.getGroupTarget(), chunk_count)
            }

        base_chunks = PE
            ? ch_trna_unmapped_reads
                .map { meta, reads -> tuple(meta, reads[0], reads[1]) }
                .splitFastq(by: params.chunk_size, pe: true, file: true)
                .map { meta, r1, r2 -> tuple(meta, [r1, r2]) }
            : ch_trna_unmapped_reads
                .splitFastq(by: params.chunk_size, file: true)
                .map { meta, r -> tuple(meta, [r]) }

        chunk_counts = [:].withDefault { 0 }
        ch_chunks = base_chunks
            .map { meta, reads ->
                def chunk = chunk_counts[meta.sample] + 1
                chunk_counts[meta.sample] = chunk
                tuple(meta.sample, meta + [id: "${meta.sample}_chunk${chunk}", chunk: chunk], reads)
            }
            .combine(ch_expected_chunks, by: 0)
            .map { _sample, meta, reads, num_chunks ->
                tuple(meta + [num_chunks: num_chunks], reads)
            }

        RNA_BWA_MAPPING(ch_chunks, hg_bwa_index_ch, params.seed_length, params.map_missing_prob, params.gap_fraction, PE)
        ch_bams_by_sample = RNA_BWA_MAPPING.out.sorted_bam
            .map { meta, bam -> tuple(groupKey(meta.sample, meta.num_chunks as int), bam) }
            .groupTuple()
            .map { sample_key, bams -> tuple([sample: sample_key.getGroupTarget()], bams) }
    }
    else {
        error "ERROR: --RNA_host_mapper must be one of 'STAR', 'bowtie2', 'BWA' (got: ${params.RNA_host_mapper})"
    }

    def merge_mapq = RNA_host_mapper == 'STAR'    ? params.STAR_mapq :
                 RNA_host_mapper == 'bowtie2' ? params.bt2_mapq :
                 params.bwa_mapq

    MERGE_BAMS(ch_bams_by_sample, PE, merge_mapq)
    ch_host_mapped   = MERGE_BAMS.out.sorted_bam
    ch_host_unmapped = MERGE_BAMS.out.unmapped_reads
    ch_multimapped_reads = MERGE_BAMS.out.multimapped_reads
    ch_multimapped_bam = MERGE_BAMS.out.multimapped_bam

    // Filter host-mapped BAMs for rRNA only when requested. The tRNA-unmapped
    // reads always go directly into host mapping above.
    if (filter_ribo_enabled) {
        FILTER_RIBO(ch_host_mapped, file(params.rrna_bed))
        ch_host_bam_with_dups = FILTER_RIBO.out.filtered_bam
        ch_ribo_intersect_bam = FILTER_RIBO.out.ribo_bam
    } else {
        ch_host_bam_with_dups = ch_host_mapped
        ch_ribo_intersect_bam = channel.empty()
    }
    
    emit:
    ribo_mapped_bam = ch_ribo_mapped_bam
    trna_mapped_bam = ch_trna_mapped_bam
    host_bam_with_dups = ch_host_bam_with_dups
    ribo_intersect_bam = ch_ribo_intersect_bam
    unmapped_reads            = ch_host_unmapped
    multimapped_reads = ch_multimapped_reads
    multimapped_bam = ch_multimapped_bam
}
