#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PATHOGEN COMPETITIVE MAPPING PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Maps single-end DNA and RNA reads competitively to strains of a pathogen,
    selects the winner by highest depth + breadth of coverage, then performs
    downstream analyses on the winning strain reference.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// ============================================================================
// PROCESSES 
// ============================================================================

/*
 * Combine all strain .fna files from the pathogen directory into one multi-FASTA
 * and build a lookup mapping sequence IDs → strain subdirectory names (for GFF linking).
 *
 * Input: the pathogen directory (e.g., spneumo_test/) containing subdirectories
 *        each with a .fna and .gff file.
 *
 * NOTE on FASTA headers: each .fna has headers like:
 *   >NC_003028.3 Streptococcus pneumoniae TIGR4, complete genome
 * The first whitespace-delimited field (NC_003028.3) is the sequence ID that
 * samtools coverage will report in its output. The lookup file maps these
 * sequence IDs back to their source subdirectory so we can find the matching GFF.
 */
process COMBINE_STRAIN_FASTAS {
    tag "${pathogen_dir.baseName}"

    input:
    tuple path(pathogen_dir), val(pathogen)

    output:
    tuple path(pathogen_dir), path("combined_strains.fa"), emit: combined_fasta
    tuple path(pathogen_dir), path("seqid_to_strain.tsv"), emit: strain_lookup

    script:
    """
    cat ${pathogen_dir}/*/*.fna > combined_strains.fa

    > seqid_to_strain.tsv
    for fna in ${pathogen_dir}/*/*.fna; do
        strain_dir=\$(basename \$(dirname "\$fna"))
        grep '^>' "\$fna" | while read -r header; do
            seq_id=\$(echo "\$header" | sed 's/^>//' | awk '{print \$1}')
            echo -e "\${seq_id}\\t\${strain_dir}" >> seqid_to_strain.tsv
        done
    done
    """
}

/*
 * Index the combined multi-FASTA with both bwa and bowtie2.
 * Outputs the FASTA alongside all index files so downstream processes
 * can stage them together in the same work directory.
 */
process INDEX_COMBINED_FASTA {
    tag "${pathogen_dir.baseName}"

    input:
    tuple path(pathogen_dir), path(combined_fasta)

    output:
    tuple path(pathogen_dir), path(combined_fasta), path("${combined_fasta}.*"), emit: bwa_indexed
    tuple path(pathogen_dir), path(combined_fasta), path("bt2_index.*"),         emit: bt2_indexed

    script:
    """
    bwa index ${combined_fasta}
    bowtie2-build ${combined_fasta} bt2_index
    """
}

/*
 * Competitive mapping: align DNA reads to the combined multi-FASTA so reads
 * compete for the best alignment across all strains simultaneously.
 * Uses bwa aln/samse (single-end), followed by Picard MarkDuplicates.
 *
 * samtools coverage output columns include:
 *   #rname  startpos  endpos  numreads  covbases  coverage  meandepth  ...
 * where #rname is the sequence ID (first field from FASTA header).
 */
process COMPETITIVE_MAP_DNA {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(dna_reads), path(combined_fasta), path(bwa_index)

    output:
    tuple val(meta), path("competitive.dedup.sorted.bam"), path("competitive.dedup.sorted.bam.bai"), emit: competitive_bam
    tuple val(meta), path("per_strain_coverage.txt"), emit: coverage_stats

    script:
    def prefix = "${meta.sample}"
    """
    # Map reads — bwa finds index files next to the fasta automatically
    bwa aln -n 0.01 -l 1000 -t ${task.cpus} ${combined_fasta} ${dna_reads} \
    | bwa samse ${combined_fasta} -r "@RG\\tID:${prefix}\\tSM:${prefix}" - ${dna_reads} \
    | samtools view -bSh | samtools sort -@ ${task.cpus} -o competitive.sorted.bam

    samtools index competitive.sorted.bam

    # Mark and remove duplicates with Picard
    java -jar \$EBROOTPICARD/picard.jar MarkDuplicates \
        -I competitive.sorted.bam \
        -O competitive.dedup.sorted.bam \
        -M competitive_dedup_metrics.txt \
        --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT

    samtools index competitive.dedup.sorted.bam

    # Per-sequence coverage stats with -q 30 map quality filter
    samtools coverage -q 30 competitive.dedup.sorted.bam > per_strain_coverage.txt
    """
}

/*
 * Select the winning strain: highest depth and breadth of coverage.
 * Parses samtools coverage output and picks the sequence with the best
 * combined depth (meandepth, col 7) and breadth (coverage%, col 6).
 *
 * Outputs the winning SEQUENCE ID (e.g., NC_003028.3) — NOT the strain dir name.
 * The next process uses the lookup table to find the corresponding GFF.
 */
process SELECT_WINNER {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(coverage_stats)

    output:
    tuple val(meta), env('WINNER_SEQID'), emit: winner

    script:
    """
    # Pick winner by: sort by breadth (col6) desc, then meandepth (col7) desc, take first
    # Skip header line (starts with #)
    WINNER_SEQID=\$(grep -v '^#' ${coverage_stats} \\
        | sort -t\$'\\t' -k6,6nr -k7,7nr \\
        | head -n1 \\
        | cut -f1)

    export WINNER_SEQID
    """
}

/*
 * Extract the winner strain FASTA from the combined file and find its
 * corresponding GFF using the sequence ID → strain directory lookup.
 */
process EXTRACT_WINNER_REF {
    tag "${meta.sample}_${meta.pathogen}_${winner_seqid}"

    input:
    tuple val(meta), val(winner_seqid), path(combined_fasta), path(strain_lookup), path(pathogen_dir)

    output:
    tuple val(meta), val(winner_seqid), path("winner.fa"), path("winner.gff"), emit: winner_ref

    script:
    """
    # Extract the winning sequence from the combined FASTA
    samtools faidx ${combined_fasta}
    samtools faidx ${combined_fasta} ${winner_seqid} > winner.fa

    # Look up which strain subdirectory this sequence ID belongs to
    STRAIN_DIR=\$(grep -w "${winner_seqid}" ${strain_lookup} | cut -f2)

    # Copy the matching GFF
    cp ${pathogen_dir}/\${STRAIN_DIR}/*.gff winner.gff
    """
}

/*
 * Index the winner FASTA with both bwa and bowtie2.
 * The FASTA + all index files are emitted together so downstream mapping
 * processes receive them in the same work directory.
 */
process INDEX_WINNER {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), val(winner_seqid), path(winner_fa), path(winner_gff)

    output:
    tuple val(meta), path(winner_fa), path("${winner_fa}.*"),  emit: bwa_indexed
    tuple val(meta), path(winner_fa), path("bt2_winner.*"),    emit: bt2_indexed
    tuple val(meta), path(winner_gff),                         emit: winner_gff

    script:
    """
    # bwa index — creates winner.fa.bwt, .pac, .ann, .amb, .sa
    bwa index ${winner_fa}

    # bowtie2 index
    bowtie2-build ${winner_fa} bt2_winner
    """
}

/*
 * Full mapping of DNA reads to the winner strain + duplicate removal + coverage stats.
 * bwa finds index files collocated with the FASTA in the work directory.
 */
process MAP_DNA_TO_WINNER {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(winner_fa), path(bwa_index), path(dna_reads)

    output:
    tuple val(meta), path("winner.dedup.sorted.bam"), path("winner.dedup.sorted.bam.bai"), emit: winner_bam
    tuple val(meta), path("coverage_stats.txt"), emit: coverage_stats

    script:
    def prefix = "${meta.sample}"
    """
    bwa aln -n 0.1 -l 1000 -t ${task.cpus} ${winner_fa} ${dna_reads} \
    | bwa samse ${winner_fa} -r "@RG\\tID:${prefix}\\tSM:${prefix}" - ${dna_reads} \
    | samtools view -bSh | samtools sort -@ ${task.cpus} -o winner.sorted.bam

    samtools index winner.sorted.bam

    # Mark and remove duplicates with Picard
    java -jar \$EBROOTPICARD/picard.jar MarkDuplicates \
        -I winner.sorted.bam \
        -O winner.dedup.sorted.bam \
        -M winner_dedup_metrics.txt \
        --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT

    samtools index winner.dedup.sorted.bam

    # Coverage stats with -q 30 map quality filter
    samtools coverage -q 30 winner.dedup.sorted.bam > coverage_stats.txt
    """
}

/*
 * mapDamage analysis on the DNA alignment to the winner.
 */
process MAPDAMAGE_DNA {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(bam), path(bai), path(winner_fa)

    output:
    tuple val(meta), path("mapDamage/*"), emit: results

    script:
    """
    mapDamage \
        -i ${bam} \
        -r ${winner_fa} \
        -d mapDamage \
        --merge-reference-sequences \
        --no-stats
    """
}

/*
 * Mosdepth high-coverage region detection on the winner BAM.

process MOSDEPTH_HIGHCOV {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("mosdepth_out/*"), emit: results

    script:
    """
    // TODO: run mosdepth to identify high-coverage regions
    """
}
 */

/*
 * Generate a consensus genome from the winner mapping.
 */
process CONSENSUS_GENOME {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(bam), path(bai), path(winner_fa)

    output:
    tuple val(meta), path("consensus.fa"), emit: consensus_fa
    tuple val(meta), path("consensus.sorted.bam"), path("consensus.sorted.bam.bai"), emit: consensus_bam

    script:
    """
    // TODO: call variants, generate consensus FASTA
    // TODO: re-map reads to consensus, sort, index
    """
}

/*
 * Annotate the consensus genome (liftoff or PGAP).
 */
process ANNOTATE_CONSENSUS {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(consensus_fa), path(winner_fa), path(winner_gff)

    output:
    tuple val(meta), path("consensus_annotated.gff"), emit: consensus_gff

    script:
    """
    // TODO: run liftoff or PGAP to transfer/create annotations on consensus
    """
}

/*
 * AMRFinder on the consensus genome.
 */
process AMRFINDER {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(consensus_fa), path(consensus_gff)

    output:
    tuple val(meta), path("amrfinder_results/*"), emit: results

    script:
    """
    // TODO: run AMRFinderPlus on consensus genome
    """
}

/*
 * Map DNA reads to virulence factor database.
 */
process VIRULENCE_FACTOR_DNA {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(dna_reads)

    output:
    tuple val(meta), path("vf_dna_results/*"), emit: results

    script:
    """
    // TODO: map DNA reads to virulence factor database
    """
}

// ============================================================================
// RNA PROCESSES
// ============================================================================

/*
 * Map RNA reads to the winner strain reference using bowtie2 with default
 * parameters. Reads are filtered at map quality -q 20.
 */
process MAP_RNA_TO_WINNER {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(winner_fa), path(bt2_index), path(rna_reads)

    output:
    tuple val(meta), path("rna_winner.sorted.bam"), path("rna_winner.sorted.bam.bai"), emit: rna_winner_bam

    script:
    def prefix = "${meta.sample_id}"
    """
    bowtie2 -x bt2_winner -U ${rna_reads} \
        --end-to-end --sensitive \
        -S ${prefix}.sam \
        --threads ${task.cpus}

    samtools addreplacerg \
        -r ID:BT2 \
        -r SM:${prefix} \
        -o ${prefix}.rg.sam ${prefix}.sam

    samtools view -bSh -q 20 -o ${prefix}.bam ${prefix}.rg.sam
    samtools sort ${prefix}.bam -@ ${task.cpus} -o rna_winner.sorted.bam

    samtools index rna_winner.sorted.bam
    samtools coverage rna_winner.sorted.bam
    """
}

/*
 * mapDamage on the RNA alignment to the winner.
 */
process MAPDAMAGE_RNA {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(bam), path(bai), path(winner_fa)

    output:
    tuple val(meta), path("mapDamage_out/*"), emit: results

    script:
    """
    // TODO: run mapDamage2 on RNA BAM against winner_fa
    """
}

/*
 * Map RNA reads to the consensus genome.
 */
process MAP_RNA_TO_CONSENSUS {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(consensus_fa), path(rna_reads)

    output:
    tuple val(meta), path("rna_consensus.sorted.bam"), path("rna_consensus.sorted.bam.bai"), emit: rna_consensus_bam

    script:
    """
    // TODO: map RNA reads to consensus genome, sort, index
    """
}

/*
 * featureCounts on the RNA consensus alignment.
 */
process FEATURECOUNTS {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(bam), path(bai), path(consensus_gff)

    output:
    tuple val(meta), path("featurecounts_results/*"), emit: results

    script:
    """
    // TODO: run featureCounts using consensus GFF annotation
    """
}

/*
 * RNA QC — infer_experiment to check strand specificity.
 */
process RNA_QC {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(bam), path(bai), path(winner_gff)

    output:
    tuple val(meta), path("infer_experiment_results.txt"), emit: results

    script:
    """
    // TODO: run RSeQC infer_experiment.py
    """
}

/*
 * Map RNA reads to virulence factor database.
 */
process VIRULENCE_FACTOR_RNA {
    tag "${meta.sample}_${meta.pathogen}"

    input:
    tuple val(meta), path(rna_reads)

    output:
    tuple val(meta), path("vf_rna_results/*"), emit: results

    script:
    """
    // TODO: map RNA reads to virulence factor database
    """
}


workflow {

    main:
    // ------------------------------------------------------------------------
    // Parameter setup
    // ------------------------------------------------------------------------
    ch_samplesheet = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample:       row.sample,
                pathogen:     row.pathogen_name,
                pathogen_dir: row.pathogen_dir   // ← added
            ]
            def dna_reads    = file(row.dna_trimmed_reads, checkIfExists: true)
            def rna_reads    = file(row.rna_trimmed_reads, checkIfExists: true)
            def pathogen_dir = file(row.pathogen_dir, checkIfExists: true)

            [ meta, dna_reads, rna_reads, pathogen_dir ]
        }

    // Derive focused channels from the samplesheet
    ch_dna_reads    = ch_samplesheet.map { meta, dna, rna, pathogen_dir -> [ meta, dna ] }
    ch_rna_reads    = ch_samplesheet.map { meta, dna, rna, pathogen_dir -> [ meta, rna ] }
    ch_pathogen_dir = ch_samplesheet.map { meta, dna, rna, pathogen_dir -> [ meta, pathogen_dir ] }

    // ========================================================================
    // DNA BRANCH — Competitive Mapping & Downstream
    // ========================================================================

    /*
     * Step 1: Combine all strain FASTAs and build seq_id → strain lookup
     *         Deduplicated by full pathogen_dir path so shared references are only built once
     */
    ch_pathogen_dir_unique = ch_pathogen_dir
        .map { meta, dir -> [dir, meta.pathogen] }
        .unique { dir, pathogen -> dir.toRealPath() }

    COMBINE_STRAIN_FASTAS(ch_pathogen_dir_unique)

    /*
     * Step 2: Index the combined multi-FASTA with bwa and bowtie2
     */
    INDEX_COMBINED_FASTA(COMBINE_STRAIN_FASTAS.out.combined_fasta)

    /*
     * Step 3: Competitive mapping — join samples back to their index via pathogen_dir path
     */


    ch_competitive_input = ch_dna_reads
        .map { meta, reads -> [file(meta.pathogen_dir).toRealPath().toString(), meta, reads] }
        .combine(
            INDEX_COMBINED_FASTA.out.bwa_indexed
                .map { pathogen_dir, fasta, index -> 
                    [pathogen_dir.toRealPath().toString(), fasta, index] 
                },
            by: 0
        )
        .map { dir, meta, reads, fasta, index -> [meta, reads, fasta, index] }

    COMPETITIVE_MAP_DNA(ch_competitive_input)

    /*
     * Step 4: Determine per-strain coverage and pick the winner
     */
    SELECT_WINNER(COMPETITIVE_MAP_DNA.out.coverage_stats)

    /*
     * Step 5: Extract the winner strain FASTA + GFF using the lookup
     */
    ch_winner_with_refs = SELECT_WINNER.out.winner
        .combine(
            COMBINE_STRAIN_FASTAS.out.combined_fasta
                .map { pathogen_dir, fasta -> [pathogen_dir.toString(), fasta] },
            by: 0  // ❌ won't work — SELECT_WINNER.out.winner has [meta, winner_seqid], not pathogen_dir
        )

    // We need to key by pathogen_dir string from meta
    ch_winner_with_refs = SELECT_WINNER.out.winner
        .map { meta, winner_seqid -> [file(meta.pathogen_dir).toRealPath().toString(), meta, winner_seqid] }
        .combine(
            COMBINE_STRAIN_FASTAS.out.combined_fasta
                .map { pathogen_dir, fasta -> [pathogen_dir.toRealPath().toString(), fasta] },
            by: 0
        )
        .combine(
            COMBINE_STRAIN_FASTAS.out.strain_lookup
                .map { pathogen_dir, lookup -> [pathogen_dir.toRealPath().toString(), lookup] },
            by: 0
        )
        .combine(
            ch_pathogen_dir
                .map { meta, dir -> [file(meta.pathogen_dir).toRealPath().toString(), dir] },
            by: 0
        )
        // shape: [pathogen_dir_str, meta, winner_seqid, combined_fasta, strain_lookup, pathogen_dir_path]
        .map { dir, meta, winner_seqid, fasta, lookup, pathogen_dir ->
            [meta, winner_seqid, fasta, lookup, pathogen_dir]
        }

    EXTRACT_WINNER_REF(ch_winner_with_refs)

    /*
     * Step 6: Index the winner FASTA with bwa and bowtie2
     */
    INDEX_WINNER(EXTRACT_WINNER_REF.out.winner_ref)

    /*
     * Step 7: Full DNA mapping to the winner strain (with dedup)
     */
    ch_dna_winner_input = INDEX_WINNER.out.bwa_indexed
        .combine(ch_dna_reads, by: 0)
        // shape: [meta, winner_fa, bwa_index_files, dna_reads]

    MAP_DNA_TO_WINNER(ch_dna_winner_input)

    /*
     * Step 8: mapDamage on DNA winner alignment
     */
    ch_mapdamage_dna_input = MAP_DNA_TO_WINNER.out.winner_bam
        .combine(
            INDEX_WINNER.out.bwa_indexed.map { meta, fa, idx -> [meta, fa] },
            by: 0
        )
        // shape: [meta, bam, bai, winner_fa]

    MAPDAMAGE_DNA(ch_mapdamage_dna_input)

    /*
     * Step 9: RNA mapping to winner
     */
    ch_rna_winner_input = INDEX_WINNER.out.bt2_indexed
        .combine(ch_rna_reads, by: 0)
        // shape: [meta, winner_fa, bt2_index, rna_reads]

    MAP_RNA_TO_WINNER(ch_rna_winner_input)

    // ========================================================================
    // PUBLISH
    // ========================================================================

    publish:
    competitive_mapping        = COMPETITIVE_MAP_DNA.out.competitive_bam
    competitive_coverage       = COMPETITIVE_MAP_DNA.out.coverage_stats
    initial_mapping_to_winner  = MAP_DNA_TO_WINNER.out.winner_bam
    dna_coverage_stats         = MAP_DNA_TO_WINNER.out.coverage_stats
    mapdamage_dna              = MAPDAMAGE_DNA.out.results
    rna_mapping_to_winner      = MAP_RNA_TO_WINNER.out.rna_winner_bam
}


// ============================================================================
// OUTPUT BLOCK — defines where each published channel lands on disk
// ============================================================================

output {

    // --- DNA: Competitive Mapping ---
    competitive_mapping {
        path { meta, bam, bai ->
            bam >> "${meta.sample}/${meta.pathogen}/DNA/competitive_mapping/"
            bai >> "${meta.sample}/${meta.pathogen}/DNA/competitive_mapping/"
        }
    }

    competitive_coverage {
        path { meta, stats ->
            stats >> "${meta.sample}/${meta.pathogen}/DNA/competitive_mapping/"
        }
    }

    // --- DNA: Initial Mapping to Winner ---
    initial_mapping_to_winner {
        path { meta, bam, bai ->
            bam >> "${meta.sample}/${meta.pathogen}/DNA/initial_mapping_to_winner/"
            bai >> "${meta.sample}/${meta.pathogen}/DNA/initial_mapping_to_winner/"
        }
    }

    dna_coverage_stats {
        path { meta, stats ->
            stats >> "${meta.sample}/${meta.pathogen}/DNA/initial_mapping_to_winner/"
        }
    }

    mapdamage_dna {
        path { meta, files ->
            files >> "${meta.sample}/${meta.pathogen}/DNA/initial_mapping_to_winner/mapDamage_results/"
        }
    }

    // --- RNA: Mapping to Winner ---
    rna_mapping_to_winner {
        path { meta, bam, bai ->
            bam >> "${meta.sample}/${meta.pathogen}/RNA/mapping_to_winner/"
            bai >> "${meta.sample}/${meta.pathogen}/RNA/mapping_to_winner/"
        }
    }

/*     mosdepth_highcov {
        path { meta, files ->
            files >> "${meta.sample}/${meta.pathogen}/DNA/initial_mapping_to_winner/mosdepth_highcov_regions/"
        }
    } */
/* 
    // --- DNA: Consensus Genome ---
    consensus_genome {
        path { meta, fa ->
            fa >> "${meta.sample}/${meta.pathogen}/DNA/consensus_genome_creation_and_mapping/"
        }
    }

    consensus_mapping {
        path { meta, bam, bai ->
            bam >> "${meta.sample}/${meta.pathogen}/DNA/consensus_genome_creation_and_mapping/"
            bai >> "${meta.sample}/${meta.pathogen}/DNA/consensus_genome_creation_and_mapping/"
        }
    }

    consensus_annotation {
        path { meta, gff ->
            gff >> "${meta.sample}/${meta.pathogen}/DNA/consensus_genome_creation_and_mapping/liftoff_or_pgap_annot/"
        }
    }

    // --- DNA: AMRFinder ---
    amrfinder {
        path { meta, files ->
            files >> "${meta.sample}/${meta.pathogen}/DNA/amrfinder/"
        }
    }

    // --- DNA: Virulence Factor Mapping ---
    virulence_factor_dna {
        path { meta, files ->
            files >> "${meta.sample}/${meta.pathogen}/DNA/virulence_factor_mapping/"
        }
    }

    

    mapdamage_rna {
        path { meta, files ->
            files >> "${meta.sample}/${meta.pathogen}/RNA/mapping_to_winner/mapDamage_results/"
        }
    }

    // --- RNA: Map to Consensus Genome ---
    rna_mapping_to_consensus {
        path { meta, bam, bai ->
            bam >> "${meta.sample}/${meta.pathogen}/RNA/map_to_consensus_genome/"
            bai >> "${meta.sample}/${meta.pathogen}/RNA/map_to_consensus_genome/"
        }
    }

    // --- RNA: featureCounts ---
    featurecounts {
        path { meta, files ->
            files >> "${meta.sample}/${meta.pathogen}/RNA/featurecounts/"
        }
    }

    // --- RNA: QC ---
    rna_qc {
        path { meta, files ->
            files >> "${meta.sample}/${meta.pathogen}/RNA/rnaQC/"
        }
    }

    // --- RNA: Virulence Factor Mapping ---
    virulence_factor_rna {
        path { meta, files ->
            files >> "${meta.sample}/${meta.pathogen}/RNA/virulence_factor_mapping/"
        }
    } */
}