process MERGE_BAMS {
    tag "$meta.sample"
    label 'demand_3'
    container 'community.wave.seqera.io/library/samtools:1.24--a1a83d6551ff8faa'

    input:
    tuple val(meta), path(bams)
    val PE
    val mapq

    output:
    tuple val(meta), path("${meta.sample}.sort.bam"),                     emit: sorted_bam
    tuple val(meta), path("*unmapped*.fastq.gz", arity: '1..*'),          emit: unmapped_reads
    tuple val(meta), path("${meta.sample}.mapped_read_count.txt"),        emit: mapped_read_count
    tuple val(meta), path("${meta.sample}.unmapped_read_count.txt"),      emit: unmapped_read_count
    script:
    if (PE) {
        """
        samtools cat -o ${meta.sample}.merged.bam ${bams}

        samtools view -bh -F4 -q ${mapq} ${meta.sample}.merged.bam \
            | samtools sort -@ ${task.cpus} -o ${meta.sample}.sort.bam -

        samtools view -c ${meta.sample}.sort.bam > ${meta.sample}.mapped_read_count.txt

        samtools view -bh -F4 -U ${meta.sample}.unmapped.bam ${meta.sample}.merged.bam > /dev/null

        samtools sort -n -@ ${task.cpus} ${meta.sample}.unmapped.bam \
            | samtools bam2fq \
                -1 ${meta.sample}_unmapped_R1.fastq.gz \
                -2 ${meta.sample}_unmapped_R2.fastq.gz \
                -0 /dev/null \
                -s /dev/null \
                --threads ${task.cpus} \
                -

        total=0
        for f in ${meta.sample}_unmapped_R1.fastq.gz ${meta.sample}_unmapped_R2.fastq.gz; do
            n=\$(zcat -f "\$f" | wc -l)
            total=\$((total + n / 4))
        done
        echo \$total > ${meta.sample}.unmapped_read_count.txt
        """
    } else {
        """
        samtools cat -o ${meta.sample}.merged.bam ${bams}

        samtools view -bh -F4 -q ${mapq} ${meta.sample}.merged.bam \
            | samtools sort -@ ${task.cpus} -o ${meta.sample}.sort.bam -

        samtools view -c ${meta.sample}.sort.bam > ${meta.sample}.mapped_read_count.txt

        samtools view -bh -F4 -U ${meta.sample}.unmapped.bam ${meta.sample}.merged.bam > /dev/null

        samtools bam2fq ${meta.sample}.unmapped.bam | gzip > ${meta.sample}.unmapped.fastq.gz

        n=\$(zcat -f ${meta.sample}.unmapped.fastq.gz | wc -l)
        echo \$((n / 4)) > ${meta.sample}.unmapped_read_count.txt
        """
    }
}