#!/usr/bin/env nextflow
include { INDEX_HOST_REFERENCE } from '../modules/bwa_index.nf'
include { HOST_DNA_MAPPING } from '../modules/host_dna_mapping.nf'
include { MERGE_BAMS } from '../modules/merge_bams.nf'

workflow ANCIENT_DNA_HOST_MAPPING {

    take:
    ch_trimmed_reads

    main:

    // prepare reference files
    bwaIdxExts = ['amb', 'ann', 'bwt', 'pac', 'sa']
    existingIdxFiles = bwaIdxExts.collect { file("${params.host_reference}.${it}") }
    hostIdxExists = existingIdxFiles.every { it.exists() }

    if (hostIdxExists) {
        ch_host_index = Channel.value( tuple(file(params.host_reference), existingIdxFiles) )
    } else {
        INDEX_HOST_REFERENCE(file(params.host_reference))
        ch_host_index = INDEX_HOST_REFERENCE.out.index.first()
    }

    // prepare read mapping in chunks
    base_chunks = params.PE
        ? ch_trimmed_reads
            // flatten [meta, [R1, R2]] -> [meta, R1, R2] so pe:true can sync-split
            .map { meta, reads -> tuple(meta, reads[0], reads[1]) }
            .splitFastq(by: params.chunk_size, pe: true, file: true)
            // normalize back to [meta, [chunkR1, chunkR2]]
            .map { meta, r1, r2 -> tuple(meta, [r1, r2]) }
        : ch_trimmed_reads
            .splitFastq(by: params.chunk_size, file: true)
            .map { meta, r -> tuple(meta, [r]) }

    ch_chunks = base_chunks
        .map { meta, reads -> tuple(meta.sample, meta, reads) }
        .groupTuple()
        .flatMap { sample, metas, reads_list ->
            reads_list.indices.collect { i ->
                tuple(metas[i] + [id: "${sample}_chunk${i + 1}", chunk: i + 1], reads_list[i])
            }
        }

    // map to host
    HOST_DNA_MAPPING(ch_chunks, ch_host_index, params.PE, params.seed_length, params.map_missing_prob, params.gap_fraction)
    ch_bams_by_sample = HOST_DNA_MAPPING.out.bam
        .map { meta, bam -> tuple(meta.sample, bam) }
        .groupTuple()
        .map { sample, bams -> tuple([sample: sample], bams) }

    MERGE_BAMS(ch_bams_by_sample, params.PE, params.map_quality)
    ch_host_unmapped = MERGE_BAMS.out.unmapped_reads

    emit:
    host_reference_index      = ch_host_index
    host_bam_with_dups = MERGE_BAMS.out.sorted_bam
    unmapped_reads = ch_host_unmapped
}


