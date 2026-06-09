include { rmDup } from './modules/picard_markdup/main.nf'
include { mapDamageHuman } from './modules/mapdamage/main.nf'

process trim {
    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("*_1.trimmed.fastq.gz"), path("*_2.trimmed.fastq.gz"), emit: trimmed_reads
    path "*.fastp.json",                                                          emit: json
    path "*.fastp.html",                                                          emit: html

    script:
    def prefix = "${meta.sample_id}"
    """
    fastp \
        -i ${reads1} -I ${reads2} \
        -o ${prefix}_1.trimmed.fastq.gz -O ${prefix}_2.trimmed.fastq.gz \
        -j ${prefix}.fastp.json \
        -h ${prefix}.fastp.html \
        -q 33 \
        -l ${params.trim_length} \
        --trim_poly_g \
        --trim_poly_x \
        --thread ${task.cpus}
    """
}

process mapBWAChunks {
    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    def prefix = "${meta.sample_id}_${reads1.baseName}"
    """
    bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} ${params.human_bwa_ref} ${reads1} > R1.sai
    bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} ${params.human_bwa_ref} ${reads2} > R2.sai

    bwa sampe ${params.human_bwa_ref} \
        -r "@RG\\tID:${meta.sample_id}\\tSM:${meta.sample_id}" \
        R1.sai R2.sai ${reads1} ${reads2} \
        | samtools view -bSh - > ${prefix}.bam
    """
}

process mergeBams {
    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.sort.bam"),                                          emit: sorted_bam
    tuple val(meta), path("*_unmapped_R1.fastq.gz"), path("*_unmapped_R2.fastq.gz"), emit: unmapped_reads

    script:
    def prefix = "${meta.sample_id}"
    """
    samtools merge -f -@ ${task.cpus} ${prefix}.merged.bam ${bams}

    samtools view -bh -F4 -q 30 ${prefix}.merged.bam \
        | samtools sort -o ${prefix}.sort.bam

    samtools view -bh -f4 ${prefix}.merged.bam > ${prefix}.unmapped.bam

    samtools sort -n -@ ${task.cpus} ${prefix}.unmapped.bam \
        | samtools bam2fq \
            -1 ${prefix}_unmapped_R1.fastq.gz \
            -2 ${prefix}_unmapped_R2.fastq.gz \
            -0 /dev/null \
            -s /dev/null \
            --threads ${task.cpus} \
            -

    rm ${prefix}.merged.bam ${prefix}.unmapped.bam
    """
}

process spades {
    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.sample_id}_spades_assembly"), emit: assembly

    script:
    def prefix = "${meta.sample_id}"
    """
    spades.py \
        -1 ${reads1} \
        -2 ${reads2} \
        -o ${prefix}_spades_assembly \
        --threads ${task.cpus} \
        --careful
    """
}

process megahit {
    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.sample_id}_megahit_assembly"), emit: assembly

    script:
    def prefix = "${meta.sample_id}"
    """
    megahit -1 ${reads1} -2 ${reads2} -o ${prefix}_megahit_assembly -t ${task.cpus}
    """
}

workflow {
    main:

    reads_ch = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample_id: row.sample_id
            ]
            tuple(meta, file(row.read1), file(row.read2))
        }

    // trim reads
    trim(reads_ch)

    // split into chunks (pe: true keeps pairs synchronized)
    chunks_ch = trim.out.trimmed_reads
        .splitFastq(by: 20000000, pe: true, file: true)

    mapBWAChunks(chunks_ch)

    grouped_chunks_ch = mapBWAChunks.out.bam.groupTuple()

    mergeBams(grouped_chunks_ch)

    // human mapping post processing
    rmDup(mergeBams.out.sorted_bam)

    // mapDamage
    mapDamageHuman(rmDup.out.bam.map { meta, bam -> bam })

    spades(mergeBams.out.unmapped_reads)
    megahit(mergeBams.out.unmapped_reads)

    publish:
    trimmed_reads       = trim.out.trimmed_reads
    human_mapped        = mergeBams.out.sorted_bam
    human_mapped_rmdup  = rmDup.out.bam
    human_rmdup_report  = rmDup.out.metrics
    human_mapdamage     = mapDamageHuman.out.damage_outdir
    human_unmapped =  mergeBams.out.unmapped_reads
    spades_contigs = spades.out.assembly
    megahit_contigs = megahit.out.assembly
}

output {
    trimmed_reads {
        path { "${params.output_dir}/${params.run_name}/trimlength_${params.trim_length}" }
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
    human_unmapped {
        path { "${params.output_dir}/${params.run_name}/hg_mapping" }
    }
    spades_contigs {
        path { "${params.output_dir}/${params.run_name}/hg_mapping" }
    }
    megahit_contigs {
        path { "${params.output_dir}/${params.run_name}/hg_mapping" }
    }
}
