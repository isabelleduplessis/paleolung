include { trim } from './modules/fastp_trim/main.nf'

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

    trim(reads_ch) // it takes in a channel with (metadata, file)

    publish:
    trimmed_reads = trim.out.trimmed_reads

}

output {
    trimmed_reads {
        path { "${params.trim_length}" }
    }
}