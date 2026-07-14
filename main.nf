include { trim } from './modules/fastp_trim/main.nf'
include { mapBowtie2HumanRibo } from './modules/map_bowtie2/main.nf'
include { mapBowtie2Chunks } from './modules/map_bowtie2/main.nf'
include { mergeBams } from './modules/merge_bams/main.nf'
include { rmDup } from './modules/picard_markdup/main.nf'
include { mapDamageHuman } from './modules/mapdamage/main.nf'
include { mapBWAChunks } from './modules/map_BWA/main.nf'
include { mapBWAHumanRibo } from './modules/map_BWA/main.nf'
include { featureCountsGenes } from './modules/featureCounts/main.nf'
include { featureCountsGenes as featureCountsGenesNotRmdup } from './modules/featureCounts/main.nf'
include { kraken } from './modules/kraken2_bracken/main.nf'
include { highC } from './modules/bbduk_highC/main.nf'
include { seqkit_stats_highC } from './modules/seqkit_stats/main.nf'
include { bracken_summary } from './modules/kraken2_bracken/main.nf'
include { dsRNA_arms_intersect } from './modules/dsRNA/main.nf'
include { dsRNA_arms_summary } from './modules/dsRNA/main.nf'
include { dsRNA_mapdamage } from './modules/dsRNA/main.nf'
include { exon_intron_coverage } from './modules/exon_intron/main.nf'
include { exon_intron_summary } from './modules/exon_intron/main.nf'

// process prePullContainers {

//     output:
//     val true, emit: ready

//     script:
//     """
//     apptainer pull docker://community.wave.seqera.io/library/bbmap:39.81--6de2cad4dc0eb52c
//     apptainer pull docker://community.wave.seqera.io/library/bedtools_pysam_samtools:a95b7b995eeeb173
//     apptainer pull docker://community.wave.seqera.io/library/bowtie2_samtools:6243cc3671a7ee48
//     apptainer pull docker://community.wave.seqera.io/library/bracken_kraken2:f7d86e30f284602f
//     apptainer pull docker://community.wave.seqera.io/library/bwa_samtools:22527ddb3cc80fb3
//     apptainer pull docker://community.wave.seqera.io/library/fastp:1.3.6--5a6797673f0eb245
//     apptainer pull docker://community.wave.seqera.io/library/mapdamage2:2.2.3--65457d747ad55400
//     apptainer pull docker://community.wave.seqera.io/library/picard_samtools:42878316047c4b3f
//     apptainer pull docker://community.wave.seqera.io/library/samtools:1.24--a1a83d6551ff8faa
//     apptainer pull docker://community.wave.seqera.io/library/seqkit:2.13.0--4cc296e17e90b043
//     apptainer pull docker://community.wave.seqera.io/library/subread:2.1.1--ab75f84c2fa06990
//     """
// }

workflow {
    main:

    // if (workflow.containerEngine == 'apptainer') {
    //     prePullContainers()
    //     ready_ch = prePullContainers.out.ready.collect()
    // } else {
    //     ready_ch = Channel.value(true)
    // }

    // prepare all input parameter file paths
    ribo_bt2_index_ch  = Channel.fromPath("${params.ribo_bt2_ref}*.bt2*").collect()
    ribo_bt2_prefix    = file(params.ribo_bt2_ref).name

    human_bt2_index_ch = Channel.fromPath("${params.human_bt2_ref}*.bt2*").collect()
    human_bt2_prefix   = file(params.human_bt2_ref).name

    ribo_bwa_ref_ch    = Channel.fromPath("${params.ribo_bwa_ref}*").collect()
    human_bwa_ref_ch   = Channel.fromPath("${params.human_bwa_ref}*").collect()

    ribo_bwa_fasta_ch  = Channel.fromPath(params.ribo_bwa_ref).collect()
    human_bwa_fasta_ch = Channel.fromPath(params.human_bwa_ref).collect()

    dsrnaome_ch        = Channel.fromPath(params.dsRNAome).collect()

    // done, move on to workflow. Next, add these as inputs to the processes that use them and update the processes. 


    reads_ch = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample_id:   row.sample_id,
                sample_type: row.sample_type
            ]
            tuple(meta, file(row.read_path))
        }

    // reads_ch = reads_ch.combine(prePullContainers.out.ready)
    //     .map { meta, reads, ready -> tuple(meta, reads) }

    // trim all reads regardless of type
    trim(reads_ch)

    // branch trimmed reads by sample_type
    trim.out.trimmed_reads.branch { meta, reads ->
        RNA: meta.sample_type == 'RNA'
        DNA: meta.sample_type == 'DNA'
    }.set { branched_trimmed }

    // ===== RNA branch: bowtie2 =====
    mapBowtie2HumanRibo(branched_trimmed.RNA, ribo_bt2_index_ch, ribo_bt2_prefix)

    rna_chunks_ch = mapBowtie2HumanRibo.out.unmapped_reads
        .splitFastq(by: 20000000, file: true)
        .groupTuple()
        .flatMap { meta, reads_list ->
            def n = reads_list.size()
            def new_meta = meta + [num_chunks: n]
            reads_list.collect { r -> tuple(new_meta, r) }
        }

    mapBowtie2Chunks(rna_chunks_ch, human_bt2_index_ch, human_bt2_prefix)

    rna_grouped_ch = mapBowtie2Chunks.out.bam
        .map { meta, bam -> tuple(groupKey(meta, meta.num_chunks), bam) }
        .groupTuple()
        .map { meta, bams -> tuple(meta.subMap('sample_id', 'sample_type'), bams) }

    // ===== DNA branch: BWA =====
    mapBWAHumanRibo(branched_trimmed.DNA, ribo_bwa_ref_ch, ribo_bwa_fasta_ch)

    dna_chunks_ch = mapBWAHumanRibo.out.unmapped_reads
        .splitFastq(by: 15000000, file: true)
        .groupTuple()
        .flatMap { meta, reads_list ->
            def n = reads_list.size()
            def new_meta = meta + [num_chunks: n]
            reads_list.collect { r -> tuple(new_meta, r) }
        }

    mapBWAChunks(dna_chunks_ch, human_bwa_ref_ch, human_bwa_fasta_ch)

    dna_grouped_ch = mapBWAChunks.out.bam
        .map { meta, bam -> tuple(groupKey(meta, meta.num_chunks), bam) }
        .groupTuple()
        .map { meta, bams -> tuple(meta.subMap('sample_id', 'sample_type'), bams) }

    // ===== Merge both branches back together =====
    grouped_chunks_ch = rna_grouped_ch.mix(dna_grouped_ch)

    mergeBams(grouped_chunks_ch)

    // featureCounts on non-rmduped BAM
    annotation_ch = channel.fromPath(params.gene_annotation)
    featureCountsGenesNotRmdup(mergeBams.out.sorted_bam, annotation_ch.collect())

    // rmDup
    rmDup(mergeBams.out.sorted_bam)

    // featureCounts on rmduped BAM
    featureCountsGenes(rmDup.out.bam, annotation_ch.collect())

    // mapDamage
    mapDamageHuman(rmDup.out.bam, human_bwa_fasta_ch)
    

    // combine ribo outputs from both branches
    ribo_sorted_bam_ch = mapBowtie2HumanRibo.out.sorted_bam
        .mix(mapBWAHumanRibo.out.sorted_bam)
    ribo_unmapped_ch = mapBowtie2HumanRibo.out.unmapped_reads
        .mix(mapBWAHumanRibo.out.unmapped_reads)

    // get highC unmapped
    highC(mergeBams.out.unmapped_reads)

    // run kraken2/bracken on highC unmapped

    kraken_db_ch = channel.fromPath(params.kraken_db, type: 'dir')

    kraken(highC.out.filtered_reads, kraken_db_ch.collect())

    // stats files
    seqkit_stats_highC(
        highC.out.filtered_reads
            .map { meta, fastq -> fastq }
            .collect()
    )

    bracken_summary(
        kraken.out.bracken
            .map { meta, bracken -> bracken }
            .collect()
    )

    // dsRNAome analysis

    dsRNA_arms_intersect(rmDup.out.bam, dsrnaome_ch)

    dsRNA_arms_summary(
        dsRNA_arms_intersect.out.feature_tsv
            .map { meta, tsv -> tsv }
            .collect()
    )

    ch_feat_bams = dsRNA_arms_intersect.out.sample_outdir
        .flatMap { meta, dir ->
            dir.listFiles()
            .findAll { it.name.endsWith('.bam') }
            .collect { bam -> [meta, bam] }
        }

    dsRNA_mapdamage(ch_feat_bams, human_bwa_fasta_ch)

    // exon intron analysis

    exon_bed_ch   = channel.fromPath(params.exon_bed)
    intron_bed_ch = channel.fromPath(params.intron_bed)

    exon_intron_coverage(rmDup.out.bam, exon_bed_ch.collect(), intron_bed_ch.collect())

    exon_intron_summary(
        exon_intron_coverage.out.coverage
            .map { meta, tsv -> tsv }
            .collect()
    )

    // done


    publish:
    trimmed_reads       = trim.out.trimmed_reads
    ribo_sorted_bam     = ribo_sorted_bam_ch
    ribo_unmapped       = ribo_unmapped_ch
    human_mapped        = mergeBams.out.sorted_bam
    human_mapped_rmdup  = rmDup.out.bam
    human_rmdup_report  = rmDup.out.metrics
    human_mapdamage     = mapDamageHuman.out.damage_outdir
    fc_notrmdup         = featureCountsGenesNotRmdup.out.counts
    fc_notrmdup_summary = featureCountsGenesNotRmdup.out.summary
    fc_rmdup            = featureCountsGenes.out.counts
    fc_rmdup_summary    = featureCountsGenes.out.summary
    highC_reads = highC.out.filtered_reads
    kraken_out     = kraken.out.kraken_out
    kraken_report  = kraken.out.report
    kraken_bracken = kraken.out.bracken
    seqkit_highC = seqkit_stats_highC.out.statsfile
    bracken_summary = bracken_summary.out.summaryfile
    dsrna_summary      = dsRNA_arms_summary.out.summary
    dsrna_mapdamage    = dsRNA_mapdamage.out.damage
    dsrna_outdirs      = dsRNA_arms_intersect.out.sample_outdir
    dsrna_feature_tsvs = dsRNA_arms_intersect.out.feature_tsv
    exon_intron_coverage = exon_intron_coverage.out.coverage
    exon_intron_summary  = exon_intron_summary.out.summary
}

output {
    trimmed_reads {
        path { meta, reads ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/trimmed_reads/"
        }
    }

    ribo_sorted_bam {
        path { meta, bam ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/ribo_mapping/"
        }
    }

    ribo_unmapped {
        path { meta, reads ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/ribo_mapping/"
        }
    }

    human_mapped {
        path { meta, bam ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/hg_mapping/"
        }
    }

    human_mapped_rmdup {
        path { meta, bam ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/hg_mapping/"
        }
    }

    human_rmdup_report {
        path { meta, metrics ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/hg_mapping/"
        }
    }

    human_mapdamage {
        path { meta, outdir ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/hg_mapping/mapdamage/"
        }
    }

    fc_notrmdup {
        path { meta, counts ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/hg_mapping/featurecounts_notrmdup/"
        }
    }

    fc_notrmdup_summary {
        path { meta, summary ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/hg_mapping/featurecounts_notrmdup/"
        }
    }

    fc_rmdup {
        path { meta, counts ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/hg_mapping/featurecounts_rmdup/"
        }
    }

    fc_rmdup_summary {
        path { meta, summary ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/hg_mapping/featurecounts_rmdup/"
        }
    }
    highC_reads {
        path { meta, f ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/trimmed_reads/"
        }
    }
    kraken_out {
        path { meta, f ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/kraken/"
        }
    }

    kraken_report {
        path { meta, f ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/kraken/"
        }
    }

    kraken_bracken {
        path { meta, f ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/kraken/"
        }
    }
    seqkit_highC {
        path { 
            "reports/"
        }
    }

    bracken_summary {
        path { 
            "reports/"
        }
    }
    dsrna_summary {
        path {
            "reports/"
        }
    }

    dsrna_outdirs {
        path { meta, dir ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/dsrna/"
        }
    }

    dsrna_feature_tsvs {
        path { meta, tsv ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/dsrna/"
        }
    }

    dsrna_mapdamage {
        path { meta, damage ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/dsrna/mapdamage/"
        }
    }
    exon_intron_coverage {
        path { meta, tsv ->
            "${meta.sample_type}/${meta.sample_id}/trimlength_${params.trim_length}/exon_intron_coverage/"
        }
    }

    exon_intron_summary {
        path {
            "reports/"
        }
    }
}


