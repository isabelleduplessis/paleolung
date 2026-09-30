#!/usr/bin/env nextflow
include { SPADES } from '../modules/spades.nf'
include { MEGAHIT } from '../modules/megahit.nf'

workflow CONTIG_ASSEMBLY {
    take:
    ch_high_complexity_reads
    PE

    main:
    def assembly_tool = params.assembly_tool?.toString()?.toLowerCase()
    if (assembly_tool == 'spades') {
        SPADES(ch_high_complexity_reads, PE)
        ch_assembly = SPADES.out.assembly
    } else if (assembly_tool == 'megahit') {
        MEGAHIT(ch_high_complexity_reads, PE)
        ch_assembly = MEGAHIT.out.assembly
    } else {
        error "ERROR: assembly_tool must be 'spades' or 'megahit', got '${assembly_tool}'"
    }

    emit:
    assembly = ch_assembly
}
