#!/usr/bin/env nextflow
include { INDEX_HOST_REFERENCE } from '../modules/bwa_index.nf'
include { HOST_DNA_MAPPING } from '../modules/host_dna_mapping.nf'
include { MERGE_BAMS } from '../modules/merge_bams.nf'

workflow ANCIENT_DNA_HOST_MAPPING {

    take:
    ch_trimmed_reads
    ch_trimmed_read_counts
    PE

    main:

    // prepare reference files
    bwaIdxExts = ['amb', 'ann', 'bwt', 'pac', 'sa']
    existingIdxFiles = bwaIdxExts.collect { ext -> file("${params.host_reference}.${ext}") }
    hostIdxExists = existingIdxFiles.every { f -> f.exists() }

    if (hostIdxExists) {
        ch_host_index = channel.value( tuple(file(params.host_reference), existingIdxFiles) )
    } else {
        INDEX_HOST_REFERENCE(file(params.host_reference))
        ch_host_index = INDEX_HOST_REFERENCE.out.index
    }

    // Use counts supplied by main.nf to determine each sample's expected
    ch_expected_chunks = ch_trimmed_read_counts
        .map { meta, count ->
            def total_reads = Math.floorDiv(count as long, (PE ? 2L : 1L))
            def reads_per_chunk = params.chunk_size as long
            def chunk_count = (total_reads + reads_per_chunk - 1L).intdiv(reads_per_chunk)
            tuple(meta.sample, chunk_count)
        }

    // prepare read mapping in chunks
    base_chunks = PE
        ? ch_trimmed_reads
            // flatten [meta, [R1, R2]] -> [meta, R1, R2] so pe:true can sync-split
            .map { meta, reads -> tuple(meta, reads[0], reads[1]) }
            .splitFastq(by: params.chunk_size, pe: true, file: true)
            // normalize back to [meta, [chunkR1, chunkR2]]
            .map { meta, r1, r2 -> tuple(meta, [r1, r2]) }
        : ch_trimmed_reads
            .splitFastq(by: params.chunk_size, file: true)
            .map { meta, r -> tuple(meta, [r]) }

    // Number chunks as they arrive so mapping can start without waiting for
    // groupTuple() to see the end of the entire input channel.
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

    // map to host
    HOST_DNA_MAPPING(ch_chunks, ch_host_index, PE, params.seed_length, params.map_missing_prob, params.gap_fraction)
    ch_bams_by_sample = HOST_DNA_MAPPING.out.bam
        .map { meta, bam -> tuple(groupKey(meta.sample, meta.num_chunks as int), bam) }
        .groupTuple()
        .map { sample_key, bams -> tuple([sample: sample_key.getGroupTarget()], bams) }

    MERGE_BAMS(ch_bams_by_sample, PE, params.map_quality, false)
    ch_host_unmapped = MERGE_BAMS.out.unmapped_reads

    emit:
    host_reference_index      = ch_host_index
    host_bam_with_dups = MERGE_BAMS.out.sorted_bam
    unmapped_reads = ch_host_unmapped
}


