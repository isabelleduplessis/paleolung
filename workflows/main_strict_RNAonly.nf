include { trim } from './modules/fastp_trim/main.nf'
include { mapBowtie2HumanRiboStrict } from './modules/map_bowtie2/main.nf'

params.hg_star_ref = '/fh/fast/blancomelo_d/user/ebarnett/reference_genomes/homo_sapiens/wholegenome/STAR'


process STAR_MAPPING {
    container 'community.wave.seqera.io/library/samtools_star:22cc46d5e973aa9c'

    input:
    tuple val(meta), path(reads)
    path index_dir

    output:
    tuple val(meta), path("*.star.bam"), emit: sorted_bam

    script:
    def prefix = "${meta.sample_id}"
    """
    STAR --runThreadN ${task.cpus} \
                  --genomeDir ${index_dir} \
                  --readFilesIn ${reads} \
                  --outSAMtype BAM SortedByCoordinate \
                  --outSAMattributes NH HI AS NM MD \
                  --alignEndsType EndToEnd \
                  --outFilterMismatchNoverReadLmax 0.1 \
                  --outFilterScoreMinOverLread 0.80 \
                  --outFilterMatchNminOverLread 0.80 \
                  --alignSJoverhangMin 10 \
                  --alignSJDBoverhangMin 5 \
                  --seedSearchStartLmax 15 \
                  --alignIntronMin 20 \
                  --alignIntronMax 500000 \
                  --alignMatesGapMax 500000 \
                  --outFilterMultimapNmax 20 \
                  --outFileNamePrefix ${prefix}_

                samtools view -bh -F4 ${prefix}_Aligned.sortedByCoord.out.bam | samtools sort -o ${prefix}.star.bam
    """
}

workflow {
    main:
    // prepare all input parameter file paths
    ribo_bt2_index_ch  = Channel.fromPath("${params.ribo_bt2_ref}*.bt2*").collect()
    ribo_bt2_prefix    = file(params.ribo_bt2_ref).name

    ribo_hg_star_index_dir = Channel.fromPath(params.hg_star_ref, type: 'dir').first() // this does the same thing as collect apparently


    reads_ch = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample_id:   row.sample_id,
                sample_type: row.sample_type
            ]
            tuple(meta, file(row.read_path))
        }
        .filter { meta, reads -> meta.sample_type == 'RNA' }

    // trim reads
    trim(reads_ch)

    // ===== RNA branch: bowtie2 =====
    mapBowtie2HumanRiboStrict(trim.out.trimmed_reads, ribo_bt2_index_ch, ribo_bt2_prefix)


    STAR_MAPPING(mapBowtie2HumanRiboStrict.out.unmapped_reads, ribo_hg_star_index_dir)


    publish:
    trimmed_reads       = trim.out.trimmed_reads
    ribo_sorted_bam     = mapBowtie2HumanRiboStrict.out.sorted_bam
    ribo_unmapped       = mapBowtie2HumanRiboStrict.out.unmapped_reads
    human_mapped        = STAR_MAPPING.out.sorted_bam

}

output {
    trimmed_reads {
        path { meta, reads ->
            "RNAstrict/${meta.sample_id}/trimlength_${params.trim_length}/trimmed_reads/"
        }
    }

    ribo_sorted_bam {
        path { meta, bam ->
            "RNAstrict/${meta.sample_id}/trimlength_${params.trim_length}/ribo_mapping/"
        }
    }

    ribo_unmapped {
        path { meta, reads ->
            "RNAstrict/${meta.sample_id}/trimlength_${params.trim_length}/ribo_mapping/"
        }
    }

    human_mapped {
        path { meta, bam ->
            "RNAstrict/${meta.sample_id}/trimlength_${params.trim_length}/hg_mapping/"
        }
    }
}


