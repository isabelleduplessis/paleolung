process exon_intron_coverage {
    input:
    tuple val(meta), path(bam)
    path exon_bed
    path intron_bed

    output:
    tuple val(meta), path("${meta.sample_id}.exon_intron_coverage.tsv"), emit: coverage

    script:
    def sample = "${meta.sample_id}"
    """
    ml SAMtools
    ml BEDTools
    ml fhPython

    samtools index ${bam}

    # breadth (cols: chr, start, end, name, score, strand, reads, bases_covered, length, fraction)
    bedtools coverage -a ${exon_bed}   -b ${bam} -s > ${sample}.exon_breadth.tsv
    bedtools coverage -a ${intron_bed} -b ${bam} -s > ${sample}.intron_breadth.tsv

    # mean depth per interval
    bedtools coverage -a ${exon_bed}   -b ${bam} -s -mean > ${sample}.exon_depth.tsv
    bedtools coverage -a ${intron_bed} -b ${bam} -s -mean > ${sample}.intron_depth.tsv

    python << 'EOF'
import csv
from collections import defaultdict

def parse_breadth(filepath):
    data = defaultdict(lambda: {'bases_covered': 0, 'length': 0})
    with open(filepath) as f:
        for line in f:
            cols = line.strip().split('\\t')
            gene_id      = cols[3].split('|')[0]
            bases_covered = int(cols[7])
            length        = int(cols[8])
            data[gene_id]['bases_covered'] += bases_covered
            data[gene_id]['length']        += length
    return data

def parse_depth(filepath):
    data = defaultdict(lambda: {'depth_sum': 0.0, 'length': 0})
    with open(filepath) as f:
        for line in f:
            cols    = line.strip().split('\\t')
            gene_id = cols[3].split('|')[0]
            length  = int(cols[2]) - int(cols[1])
            depth   = float(cols[6]) if cols[6] != '.' else 0.0
            data[gene_id]['depth_sum'] += depth * length
            data[gene_id]['length']    += length
    return data

exon_breadth   = parse_breadth('${sample}.exon_breadth.tsv')
intron_breadth = parse_breadth('${sample}.intron_breadth.tsv')
exon_depth     = parse_depth('${sample}.exon_depth.tsv')
intron_depth   = parse_depth('${sample}.intron_depth.tsv')

all_genes = set(exon_breadth.keys()) | set(intron_breadth.keys())

with open('${sample}.exon_intron_coverage.tsv', 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\\t')
    writer.writerow([
        'sample', 'gene_id',
        'exon_bases_covered', 'exon_length', 'exon_breadth', 'exon_mean_depth',
        'intron_bases_covered', 'intron_length', 'intron_breadth', 'intron_mean_depth'
    ])
    for gene in sorted(all_genes):
        e_b = exon_breadth.get(gene,   {'bases_covered': 0, 'length': 0})
        i_b = intron_breadth.get(gene, {'bases_covered': 0, 'length': 0})
        e_d = exon_depth.get(gene,     {'depth_sum': 0.0, 'length': 0})
        i_d = intron_depth.get(gene,   {'depth_sum': 0.0, 'length': 0})

        exon_breadth_val   = e_b['bases_covered'] / e_b['length'] if e_b['length'] > 0 else 0
        intron_breadth_val = i_b['bases_covered'] / i_b['length'] if i_b['length'] > 0 else 0
        exon_depth_val     = e_d['depth_sum'] / e_d['length']     if e_d['length'] > 0 else 0
        intron_depth_val   = i_d['depth_sum'] / i_d['length']     if i_d['length'] > 0 else 0

        writer.writerow([
            '${sample}', gene,
            e_b['bases_covered'], e_b['length'], round(exon_breadth_val,   4), round(exon_depth_val,   4),
            i_b['bases_covered'], i_b['length'], round(intron_breadth_val, 4), round(intron_depth_val, 4)
        ])
EOF
    """
}

process exon_intron_summary {
    input:
    path tsv_files

    output:
    path "exon_intron_summary.tsv", emit: summary

    script:
    """
    awk 'FNR==1 && NR!=1 {next} 1' *.exon_intron_coverage.tsv > exon_intron_summary.tsv
    """
}