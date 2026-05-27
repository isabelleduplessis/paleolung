include { trim } from './modules/fastp_trim/main.nf'
include { mapBowtie2HumanRibo } from './modules/map_bowtie2/main.nf'
include { mapBowtie2Chunks } from './modules/map_bowtie2/main.nf'
include { mergeBams } from './modules/merge_bams/main.nf'
include { rmDup } from './modules/picard_markdup/main.nf'
include { mapDamageHuman } from './modules/mapdamage/main.nf'
include { mapBWAChunks } from './modules/map_BWA/main.nf'
include { mapBWAHumanRibo } from './modules/map_BWA/main.nf'

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

    // Declare channels that are branch-specific
    def ribo_sorted_bam_ch = channel.empty()
    def ribo_unmapped_ch = channel.empty()
    def grouped_chunks_ch = channel.empty()

    if (params.sample_type == 'RNA') {
        // map to ribo
        mapBowtie2HumanRibo(trim.out.trimmed_reads)

        ribo_sorted_bam_ch = mapBowtie2HumanRibo.out.sorted_bam
        ribo_unmapped_ch = mapBowtie2HumanRibo.out.unmapped_reads

        // map to human
        chunks_ch = mapBowtie2HumanRibo.out.unmapped_reads
            .splitFastq(by: 20000000, file: true)

        mapBowtie2Chunks(chunks_ch)

        grouped_chunks_ch = mapBowtie2Chunks.out.bam.groupTuple()

    } else if (params.sample_type == 'DNA') {

        mapBWAHumanRibo(trim.out.trimmed_reads)

        ribo_sorted_bam_ch = mapBWAHumanRibo.out.sorted_bam
        ribo_unmapped_ch = mapBWAHumanRibo.out.unmapped_reads

        chunks_ch = mapBWAHumanRibo.out.unmapped_reads.splitFastq(by: 20000000, file: true)

        mapBWAChunks(chunks_ch)
        
        grouped_chunks_ch = mapBWAChunks.out.bam.groupTuple()

    }

    mergeBams(grouped_chunks_ch)

    // human mapping post processing
    rmDup(mergeBams.out.sorted_bam)

    // mapDamage
    mapDamageHuman(rmDup.out.bam.map { meta, bam -> bam })

    publish:
    trimmed_reads       = trim.out.trimmed_reads
    ribo_sorted_bam     = ribo_sorted_bam_ch
    ribo_unmapped       = ribo_unmapped_ch
    human_mapped        = mergeBams.out.sorted_bam
    human_mapped_rmdup  = rmDup.out.bam
    human_rmdup_report  = rmDup.out.metrics
    human_mapdamage     = mapDamageHuman.out.damage_outdir
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
    human_mapped {
        path { "${params.output_dir}/${params.run_name}/hg_mapping" }
    }
    human_mapped_rmdup {
        path { "${params.output_dir}/${params.run_name}/hg_mapping" }
    }
    human_rmdup_report {
        path { "${params.output_dir}/${params.run_name}/hg_mapping" }
    }
    human_mapdamage {
        path { "${params.output_dir}/${params.run_name}/hg_mapping" }
    }
}