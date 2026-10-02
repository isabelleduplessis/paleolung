process STAR_MAPPING {
    tag "${meta.sample}"
    label  'demand_6'
    container 'community.wave.seqera.io/library/samtools_star:22cc46d5e973aa9c'

    input:
    tuple val(meta), path(reads)
    path index_dir
    val PE
    val max_hits

    output:
    tuple val(meta), path("*.star.bam"), emit: sorted_bam

    script:
    if (PE) {
        """
        STAR --runThreadN ${task.cpus} \
                    --genomeDir ${index_dir} \
                    --readFilesIn ${reads[0]} ${reads[1]} \
                    --outSAMtype BAM SortedByCoordinate \
                    --outSAMunmapped Within \
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
                    --outFilterMultimapNmax ${max_hits} \
                    --outFileNamePrefix ${meta.sample}_

        mv ${meta.sample}_Aligned.sortedByCoord.out.bam ${meta.sample}.star.bam
        """
    } else {
        """
        STAR --runThreadN ${task.cpus} \
                    --genomeDir ${index_dir} \
                    --readFilesIn ${reads} \
                    --outSAMtype BAM SortedByCoordinate \
                    --outSAMunmapped Within \
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
                    --outFilterMultimapNmax ${max_hits} \
                    --outFileNamePrefix ${meta.sample}_

        mv ${meta.sample}_Aligned.sortedByCoord.out.bam ${meta.sample}.star.bam
        """
    }
}