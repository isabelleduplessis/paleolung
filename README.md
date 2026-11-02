# Paleolung

## Getting Started

### Installing and setting up Nextflow version 26

Our cluster doesn't have the most recent version of nextflow so we have to install it:

```
ml Java/17.0.6
cd ~/.local/bin
export NXF_VER=26.04.1
curl -s https://get.nextflow.io | bash
chmod +x nextflow
export PATH="$PATH:$HOME/.local/bin"
```
Add the export path line to ~/.bash_profile and ~/.bashrc as well if not already:
```
echo 'export PATH="$PATH:$HOME/.local/bin"' >> ~/.bashrc && echo 'export PATH="$PATH:$HOME/.local/bin"' >> ~/.bash_profile
```


To run nextflow v26, Java version 17 or later is needed, which we will load when running the pipeline. Nextflow 26 uses v2 strict syntax parser.


### Set up configuration file


See how to configure here: https://sciwiki.fredhutch.org/datademos/on_gizmo/

You might already have this if you've run nextflow on our cluster before. I would store this in a safe directory that will be easy to access all the time from anywhere for any pipeline. 


Example gizmo config file `nextflow.gizmo.config`:
```
workDir = '/hpc/temp/path/to/tmp/dir/work/'

apptainer {
    enabled = true
    autoMounts = true
    cacheDir = '/hpc/temp/path/to/tmp/dir/cache/'
    runOptions = '--containall -B \$TMPDIR'
}

env {
    TMPDIR="\$TMPDIR"
}

process {
    executor = 'slurm'
    queue = 'campus-new'

    // below address apptainer pull error, allegedly
    errorStrategy = 'retry'
    maxRetries = 3
}

// hard copy all result files, so that if tmp dir is cleared you still have results. 

workflow {
    output.mode = 'copy'
}
```

### Get pipeline
```
git clone https://github.com/isabelleduplessis/paleolung.git
cd paleolung
```

### Run 

Example nf-run script on HPC cluster after installation (run as sbatch job and add sbatch params as needed):

```
#!/bin/bash
#SBATCH --job-name=nf-PALEO

set -euo pipefail

# Nextflow Version
NXF_VER=26.04.1

# Nextflow Configuration File
NXF_CONFIG=/path/to/nextflow.gizmo.config

# Load java
ml Java/17.0.6

# Load the Apptainer module
ml Apptainer
# Make sure that the apptainer executables are in the PATH
export PATH=$APPTAINERROOT/bin/:$PATH

nf_samplesheet=/path/to/samplesheet.csv
nf_outdir=outdir_name

# Run the workflow
NXF_VER=$NXF_VER \
nextflow \
    run \
    -c ${NXF_CONFIG} \
    main.nf \
    --samplesheet ${nf_samplesheet} --sample_type RNA --read_type se \
    -o ${nf_outdir} -resume

```

The first time you run the pipeline, it will take a while because it has to download all the docker images

## Parameters

### Required

`--sample_type` (required): 
- `DNA` or `RNA`

`--samplesheet` (required): 
- CSV file with a header containing `sample` and `R1`
columns. Add an `R2` column when you have paired end reads. Each row must have a unique sample name.



### More Options

All default values are specified in nextflow.config. You can change the parameter values in the config file or by adding them to the `nextflow run` command. 

For any input files where default is in `${projectDir}/references`, I can give you the path to copy them over.

| Parameter | Description | Input Options | Default |
|---|---|---|---|
| `--read_type` | Read layout. `se` is single-end, `pe` is paired-end, and `merge_pe` merges paired-end reads with fastp (requires `--trimmer fastp`) | `se`, `pe`, `merge_pe` | `se` |
| `--trimmer` | Trimming tool | `trim_galore`, `fastp` | `trim_galore` |
| `--trim_length` | Minimum read length kept after trimming | Integer | `30` |
| `--trim_quality` | Quality score cutoff used in trimming | Integer | `33` |
| `--min_overlap` | Minimum overlap length used by fastp to merge paired-end reads | Integer | `30` |
| `--host_reference` | Host fasta file used for DNA host mapping, mapDamage, and dsRNA analysis | Path to fasta file | `${projectDir}/references/GRCh38.p14.genome.fa` |
| `--map_quality` | Minimum mapping quality for BWA host alignment for DNA | Integer, 0–37 | `30` |
| `--map_missing_prob` | `-n` in BWA alignment | Number: fraction between 0 and 1, or integer edit distance | `0.01` |
| `--seed_length` | `-l` in BWA alignment | Integer | `1024` |
| `--gap_fraction` | `-o` in BWA alignment | Integer| `2` |
| `--chunk_size` | Reads per mapping job in BWA host alignment | Integer | `20000000` |
| `--skip_host_mapping` | Set to `true` when the samplesheet already contains trimmed reads that do not need host reads removed. Runs taxonomic profiling and optionally contig assembly on the reads provided in the samplesheet | `true`, `false` | `false` |
| `--publish_host_bam` | Output deduplicated host BAM files. Set to `false` if you don't need the large BAM files containing host-mapped reads | `true`, `false` | `true` |
| `--kraken2_db` | Kraken2 and Bracken database files | Path to database directory | `/shared/biodata/microbiome/kraken2/k2_standard_20230605` |
| `--bracken_read_len` | Read length for Bracken | Integer, (should match the available database k-mer distribution file) | `100` |
| `--bracken_tax_level` | Taxonomic level to re-estimate read counts to in Bracken | `D` (domain), `P` (phylum), `C` (class), `O` (order), `F` (family), `G` (genus), `S` (species) | `S` |
| `--krona_min_fraction` | Minimum abundance (in percent) for a taxon to appear in Krona output | Number, 0–100 | `0.5` |
| `--assembly_tool` | Contig assembly tool | `spades`, `megahit`, `null` (skips assembly) | `megahit` |


### RNA-Specific Parameters

| Parameter | Description | Input Options | Default |
|---|---|---|---|
| `--ncRNA_mapper` | Noncoding RNA mapper (for tRNA and rRNA) | `bowtie2`, `BWA` | `bowtie2` |
| `--ribo_bt2_ref` | bowtie2 index for rRNA (used if `--ncRNA_mapper bowtie2`) | Path to bowtie2 index prefix | `${projectDir}/references/rRNA.bt2` |
| `--trna_bt2_ref` | bowtie2 index for tRNA (used if `--ncRNA_mapper bowtie2`) | Path to bowtie2 index prefix | `${projectDir}/references/hg38-tRNAs` |
| `--ribo_bwa_ref` | Pre-indexed rRNA fasta file (used if `--ncRNA_mapper BWA`) | Path to BWA-indexed fasta | `${projectDir}/references/rRNA.fa` |
| `--trna_bwa_ref` | Pre-indexed tRNA fasta file (used if `--ncRNA_mapper BWA`) | Path to BWA-indexed fasta | `${projectDir}/references/hg38-tRNAs.fa` |
| `--RNA_host_mapper` | Host RNA mapper | `STAR`, `bowtie2`, `BWA` | `bowtie2` |
| `--bt2_ref` | bowtie2 index of host reference (used if `--RNA_host_mapper bowtie2`) | Path to bowtie2 index prefix | `${projectDir}/references/GRCh38.p14.genome.bt2` |
| `--hg_star_ref` | STAR index of host reference (used if `--RNA_host_mapper STAR`) | Path to STAR index directory | `${projectDir}/references/STAR_index` |
| `--host_reference` | Host reference fasta used as the reference when `--RNA_host_mapper BWA` | Path to fasta file | `${projectDir}/references/GRCh38.p14.genome.fa` |
| `--bwa_mapq` | Minimum mapping quality for the BWA host mapper | Integer, 0–37 | `20` |
| `--bt2_mapq` | Minimum mapping quality for the bowtie2 host mapper | Integer, 0–42 (or 40 for paired end) | `20` |
| `--STAR_mapq` | Minimum mapping quality for the STAR host mapper | `0`, `1`, `3`, or `255` | `20` |
| `--max_hits` | Maximum reported alignment hits for RNA host mapping if allowing multimappers | Integer, > 1 | `100` |
| `--multi` | If `true`, include multimapped reads fractionally in RNA featureCounts and retain multimapped host alignments during BAM merging | `true`, `false` | `false` |
| `--hg_gtf` | Host genome GTF used for feature counting | Path to GTF file | `${projectDir}/references/gencode.v49.chr_patch_hapl_scaff.annotation.withrmsk.gtf` |
| `--dsRNA_bed` | dsRNA annotation BED | Path to BED file | `${projectDir}/references/dsRNA_RNAduplex_sort.bed` |
| `--strand` | Strandedness setting for featureCounts and dsRNA analysis | `0` (unstranded), `1` (forward stranded), `2` (reverse stranded) | `0` (unstranded)|
| `--filter_ribo` | If `true`, filter host-mapped RNA reads overlapping with `rrna_bed` and publish the intersecting rRNA BAM | `true`, `false` | `false` |
| `--rrna_bed` | Bedfile containing rRNA annotations | Path to bed file | `${projectDir}/references/rmsk_annot_rRNA.bed` |

`--ncRNA_mapper`: noncoding RNA mapper (for tRNA and rRNA), `bowtie2` or `BWA`, default `bowtie2`

- If `bowtie2`:
    - `--ribo_bt2_ref`: bowtie2 index for rRNA
        - default
`${projectDir}/references/rRNA.bt2`
    - `--trna_bt2_ref`: bowtie2 index for tRNA
        - default
`${projectDir}/references/hg38-tRNAs`
- If `BWA`:
    - `--ribo_bwa_ref`: pre indexed rRNA fasta file
        - default `${projectDir}/references/rRNA.fa`
    - `--trna_bwa_ref`: pre indexed tRNA fasta file
        - default `${projectDir}/references/hg38-tRNAs.fa`


## Outputs

The workflow publishes the following named outputs under the output directory
specified with `-o`

`read_count_report.tsv` contains stats for all samples including: sample, sample_type, raw_reads, trimmed_reads, ribo_mapped, trna_mapped, ribo_intersect, host_mapped, host_mapped_deduped, host_multimapped, host_unmapped, host_unmapped_high_complexity

Example output structure for RNA:
```
<outdir>/
├── reads/
│   ├── fastqc/
│   │   ├── raw/
│   │   │   └── {sample}/
│   │   │       ├── {sample}_fastqc.html
│   │   │       └── {sample}_fastqc.zip
│   │   └── trimmed/
│   │       └── {sample}/
│   │           ├── {sample}_trimmed_fastqc.html
│   │           └── {sample}_trimmed_fastqc.zip
│   ├── trimmed/
│   │   └── {sample}/
│   │       └── {sample}_trimmed.fq.gz
│   └── high_complexity_reads/
│       └── {sample}.highC.fastq.gz
├── host_mapping/
│   └── {RNA_host_mapper}/
│       └── {sample}/
│           ├── {sample}.filtered.bam
│           ├── {sample}.dedup.bam
│           └── {sample}.dedup.bam.bai
├── damage_patterns/
│   └── {RNA_host_mapper}/
│       ├── host_mapped/
│       │   └── {sample}_mapdamage/
│       ├── dsRNA/
│       │   └── {sample}_dsRNA_mapdamage/
│       └── ssRNA/
│           └── {sample}_ssRNA_mapdamage/
├── plots/
│   └── {RNA_host_mapper}/
│       ├── read_length_distribution/
│       │   ├── {sample}.read_length_distribution.png
│       │   └── {sample}.read_length_distribution.tsv
│       └── biotype_composition/
│           ├── biotype_composition_gene.pdf
│           └── biotype_composition_exon.pdf
├── stats/
│   └── {RNA_host_mapper}/
│       ├── read_count_report.tsv
│       ├── featurecounts_summary_gene.csv
│       └── featurecounts_summary_exon.csv
├── contigs/
│   └── megahit/
│       └── {sample}_megahit_assembly/
└── taxonomic_profiling/
    ├── bracken_summary.tsv
    └── krona/
        └── {sample}.krona.html
```

Example output structure for DNA:
```
<outdir>/
├── reads/
│   ├── fastqc/
│   │   ├── raw/
│   │   │   └── {sample}/
│   │   │       ├── {sample}_fastqc.html
│   │   │       └── {sample}_fastqc.zip
│   │   └── trimmed/
│   │       └── {sample}/
│   │           ├── {sample}_trimmed_fastqc.html
│   │           └── {sample}_trimmed_fastqc.zip
│   ├── trimmed/
│   │   └── {sample}/
│   │       └── {sample}_trimmed.fq.gz
│   └── high_complexity_reads/
│       └── {sample}.highC.fastq.gz
├── host_mapping/
│   └── {sample}/
│       ├── {sample}.sort.bam
│       ├── {sample}.dedup.bam
│       └── {sample}.dedup.bam.bai
├── damage_patterns/
│   └── host_mapped/
│       └── {sample}_mapdamage/
├── plots/
│   └── read_length_distribution/
│       ├── {sample}.read_length_distribution.png
│       └── {sample}.read_length_distribution.tsv
├── stats/
│   └── read_count_report.tsv
├── contigs/
│   └── megahit/
│       └── {sample}_megahit_assembly/
└── taxonomic_profiling/
    ├── bracken_summary.tsv
    └── krona/
        └── {sample}.krona.html
```


## Errors
If you get a "Failed to pull singularity image" error on FH HPC, just rerun with -resume