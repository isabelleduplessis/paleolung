include { trim } from './modules/fastp_trim/main.nf'
include { mapBowtie2HumanRibo } from './modules/map_bowtie2/main.nf'
include { mapBowtie2Chunks } from './modules/map_bowtie2/main.nf'
include { mergeBams } from './modules/merge_bams/main.nf'

workflow {
    main:

    reads_ch = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
        def meta = [
            sample_id: row.sample_id
        ]
        tuple(meta, file(row.read_path))
        }

    // trim reads
    trim(reads_ch)

    // map to ribo
    mapBowtie2HumanRibo(trim.out.trimmed_reads)

    // map to human
    chunks_ch = mapBowtie2HumanRibo.out.unmapped_reads
        .splitFastq(by: 20000000, file: true)

    mapBowtie2Chunks(chunks_ch)
    grouped_chunks_ch = mapBowtie2Chunks.out.bam.groupTuple()

    mergeBams(grouped_chunks_ch)

    // human mapping post processing



    // PUBLISH OUTPUTS 

    publish:
    trimmed_reads = trim.out.trimmed_reads
    ribo_sorted_bam = mapBowtie2HumanRibo.out.sorted_bam
    ribo_unmapped = mapBowtie2HumanRibo.out.unmapped_reads

}

output {
    trimmed_reads {
        path { "${params.output_dir}/${params.run_name}/trimlength_${params.trim_length}" }
    }
    ribo_sorted_bam {
        path { "${params.output_dir}/${params.run_name}/ribo_mapping" }
    }
    ribo_unmapped {
        path { "${params.output_dir}/${params.run_name}/ribo_mapping" }
    }
}