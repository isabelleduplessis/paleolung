process STAR_MAPPING {
    tag "${meta.sample}"
    label  'demand_6'
    container 'community.wave.seqera.io/library/samtools_star:22cc46d5e973aa9c'

    input:
    tuple val(meta), path(reads)
    path index_dir
    val PE
    val mapq

    output:
    tuple val(meta), path("*.star.bam"), emit: sorted_bam
    tuple val(meta), path("*.unmapped*.fastq"), emit: unmapped_reads
    tuple val(meta), path("*.multimapped*.fastq"), emit: multimapped_reads
    tuple val(meta), path("${meta.sample}.multi.mapped.bam"), emit: multimapped_bam

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
                    --outFilterMultimapNmax 20 \
                    --outFileNamePrefix ${meta.sample}_

        samtools view -bh -F4 \
            -U ${meta.sample}.unmapped.bam \
            -o ${meta.sample}.mapped.bam \
            ${meta.sample}_Aligned.sortedByCoord.out.bam

        samtools view -bh -q ${mapq} -F 0x100 \
            -U ${meta.sample}.multi.mapped.bam \
            -o ${meta.sample}.unique.mapped.bam \
            ${meta.sample}.mapped.bam

        samtools sort -o ${meta.sample}.star.bam ${meta.sample}.unique.mapped.bam

        samtools collate -O -u ${meta.sample}.unmapped.bam | samtools bam2fq -1 ${meta.sample}.unmapped_1.fastq -2 ${meta.sample}.unmapped_2.fastq -0 /dev/null -s /dev/null -n -
        samtools collate -O -u ${meta.sample}.multi.mapped.bam | samtools bam2fq -1 ${meta.sample}.multimapped_1.fastq -2 ${meta.sample}.multimapped_2.fastq -0 /dev/null -s /dev/null -n -
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
                    --outFilterMultimapNmax 10 \
                    --outFileNamePrefix ${meta.sample}_

        samtools view -bh -F4 \
            -U ${meta.sample}.unmapped.bam \
            -o ${meta.sample}.mapped.bam \
            ${meta.sample}_Aligned.sortedByCoord.out.bam

        samtools view -bh -q ${mapq} -F 0x100 \
            -U ${meta.sample}.multi.mapped.bam \
            -o ${meta.sample}.unique.mapped.bam \
            ${meta.sample}.mapped.bam

        samtools sort -o ${meta.sample}.star.bam ${meta.sample}.unique.mapped.bam

        samtools bam2fq ${meta.sample}.unmapped.bam > ${meta.sample}.unmapped.fastq
        samtools bam2fq ${meta.sample}.multi.mapped.bam > ${meta.sample}.multimapped.fastq
        """
    }
}