# Paleolung

## Installing and setting up Nextflow version 26

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


## Set up configuration
Example gizmo config file: nextflow.gizmo.config

See how to configure here: https://sciwiki.fredhutch.org/datademos/on_gizmo/

You might already have this if you've run nextflow on our cluster before. I would store this in a safe directory that will be easy to access all the time from anywhere for any pipeline. 
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
// does double amount of storage being used until tmp dir clears
workflow {
    output.mode = 'copy'
}
```

## Get pipeline
```
git clone https://github.com/isabelleduplessis/paleolung.git
cd paleolung
```

## Run 

Example nf-run script on HPC cluster after installation (run as sbatch job and add sbatch params as needed):

```
#!/bin/bash

# Nextflow run
set -euo pipefail

# Nextflow Version
NXF_VER=26.04.1

# Nextflow Configuration File
NXF_CONFIG=nextflow.gizmo.config

# Load the Nextflow module (if running on rhino/gizmo)
ml Java/17.0.6

# Load the Apptainer module
ml Apptainer
# Make sure that the apptainer executables are in the PATH
export PATH=$APPTAINERROOT/bin/:$PATH

# Run the workflow
NXF_VER=$NXF_VER \
nextflow \
    run \
    -c ${NXF_CONFIG} \
    main.nf \
    --samplesheet ${nf_samplesheet} --sample_type RNA \
    -o ${nf_outdir} --RNA_host_mapper bowtie2 --trim_length 30

# Log commands

tmpl=$(mktemp)
trap 'rm -f "$tmpl"' EXIT
cat > "$tmpl" << 'EOF'
# ${process}
${script}
EOF

s=$(nextflow log last -f tag | grep -v '^-$' | sed -n 1p)
nextflow log last -F "tag == \"$s\"" -t "$tmpl" \
    | sed "s/$s/\${sample}/g" > "$nf_outdir/commands.txt"

sed -i -E 's/^[ \t]+//; s/[ \t]+$//; s/[ \t]{2,}/ /g' "$nf_outdir/commands.txt"

echo "Done!"
```

The first time you run it, it will take a while because it has to download all the docker images

## Required input parameters

Pass parameters on the command line with `--parameter_name parameter_value`, or set
them in the Nextflow configuration file.

The defaults below are those in `nextflow.config`. I have path to references/ folder if you'd like to use those.

### Read input

`samplesheet` (required): CSV file with a header containing `sample` and `R1`
columns. Add an `R2` column when `read_type` is `pe` or `merge_pe`. All rows
must use the same read layout, and each row must have a unique `sample` name.

`read_type` (default `se`): select `se` for single-end reads, `pe` for paired-end
reads, or `merge_pe` to merge paired-end reads with fastp. Both paired layouts
require `R2` values. `merge_pe` requires `--trimmer fastp`; merged reads and
unmerged singletons are treated as single-end downstream. It cannot be combined
with `--skip_host_mapping`.

`sample_type` (required): `DNA` or `RNA`.

`skip_host_mapping` (default `false`): set to `true` when the samplesheet already
contains high-complexity reads. This skips preprocessing and host mapping and runs
taxonomic profiling (and optional assembly) directly on the supplied reads.

### Common processing parameters

`trimmer` (default `trim_galore`): choose `trim_galore` or `fastp`.

`trim_length` (default `30`): minimum read length retained by the selected trimmer.

`trim_quality` (default `33`): quality threshold used by the selected trimmer.

`min_overlap` (default `30`): minimum overlap length used by fastp with
`read_type=merge_pe`.

`host_reference` (default `${projectDir}/references/GRCh38.p14.genome.fa`): host
reference FASTA. Used for DNA host mapping, mapDamage, and dsRNA analysis; RNA
BWA host mapping also uses this FASTA and its BWA index files.

`map_quality` (default `30`): minimum mapping quality for DNA BAM merging.

`map_missing_prob` (default `0.01`), `seed_length` (default `1024`), and
`gap_fraction` (default `2`): BWA alignment settings.

`chunk_size` (default `20000000`): reads per mapping chunk for chunked BWA mapping.

`publish_host_bam` (default `true`): publish deduplicated host BAM files as
`host_bam_deduped`.

`assembly_tool` (default `megahit`): choose `megahit` or `spades` to assemble
high-complexity reads. Set `assembly_tool = null` in the config file to skip
assembly. Both tools support single-end and paired-end inputs.

### RNA parameters

`ncRNA_mapper` (default `bowtie2`): choose `bowtie2` or `BWA` for sequential
rRNA and tRNA mapping. For Bowtie2, set `ribo_bt2_ref` (default
`${projectDir}/references/rRNA.bt2`) and `trna_bt2_ref` (default
`${projectDir}/references/hg38-tRNAs`). For BWA, set `ribo_bwa_ref` (default
`${projectDir}/references/rRNA.fa`) and `trna_bwa_ref` (default
`${projectDir}/references/hg38-tRNAs.fa`); prebuilt BWA index sidecars must be
available for each prefix.

`RNA_host_mapper` (default `bowtie2`): choose `STAR`, `bowtie2`, or `BWA` for
host mapping. STAR requires `hg_star_ref` (default
`${projectDir}/references/STAR_index`). Bowtie2 requires `bt2_ref` (default
`${projectDir}/references/GRCh38.p14.genome.bt2`). BWA uses `host_reference`
and its index sidecars.

`bwa_mapq`, `bt2_mapq`, and `STAR_mapq` (each default `20`): minimum mapping
quality for the selected RNA host mapper.

`max_hits` (default `100`): maximum reported alignment hits for RNA host mapping.

`multi` (default `false`): include multimapped reads fractionally in RNA
featureCounts and retain multimapped host alignments during BAM merging when set
to `true`.

`hg_gtf` (default
`${projectDir}/references/gencode.v49.chr_patch_hapl_scaff.annotation.withrmsk.gtf`):
host genome GTF used for feature counting.

`dsRNA_bed` (default `${projectDir}/references/dsRNA_RNAduplex_sort.bed`):
dsRNA annotation BED used by dsRNA analysis.

`strand` (default `0`): strandedness setting for feature counting and dsRNA
analysis.

`filter_ribo` (default `false` in `nextflow.config`): set to `true` to filter
host-mapped RNA reads against `rrna_bed` (default
`${projectDir}/references/rmsk_annot_rRNA.bed`) and publish the intersecting
rRNA BAM. In that mode, `ribo_intersect` is also added to the RNA summary
statistics.

### Taxonomic profiling parameters

`kraken2_db` (default `/shared/biodata/microbiome/kraken2/k2_standard_20230605`):
directory containing the Kraken2 database and Bracken database files; required
for taxonomic profiling. Set this to an accessible database path if the configured
Fred Hutch path is unavailable.

Taxonomic profiling runs Kraken2/Bracken and Krona on high-complexity reads
remaining after host mapping, or directly on the supplied reads with
`--skip_host_mapping`. `bracken_read_len` (default `100`) is the Bracken read
length, `bracken_tax_level` (default `S`) is the Bracken taxonomic level, and
`krona_min_fraction` (default `0.5` percent) is the minimum Krona input fraction.

## Outputs

The workflow publishes the following named outputs under the output directory
specified with `-o`. Outputs are conditional on the selected sample type and
options; outputs for skipped workflow stages are not created.

- `fastqc_raw`: FastQC results for raw reads; absent with `--skip_host_mapping`.
- `fastqc_trimmed`: FastQC results for trimmed reads; absent with
    `--skip_host_mapping`.
- `trimmed_reads`: reads retained after trimming; absent with
    `--skip_host_mapping`.
- `host_bam`: host-mapped BAM files before duplicate removal.
- `host_bam_deduped`: deduplicated host BAM files when
    `publish_host_bam=true` (the configured default).
- `mapdamage`: mapDamage results for deduplicated host mappings.
- `read_length_distribution`: read-length plots for deduplicated host mappings.
- `high_complexity_reads`: reads retained after read-complexity filtering, or
    the input reads supplied with `--skip_host_mapping`.
- `summary_stats`: `read_count_report.tsv`, with read counts by sample and
    processing stage; absent with `--skip_host_mapping`. RNA summaries include
    `ribo_mapped` and `trna_mapped`; `ribo_intersect` is included when
    `--filter_ribo true`.
- `dsRNA_dir` and `ssRNA_dir`: dsRNA- and ssRNA-associated mapping and mapDamage
    output directories for RNA runs.
- `fc_genes` and `fc_exons`: combined featureCounts summary CSVs for RNA runs.
    In the current workflow wiring, the exon summary is assigned to `fc_genes` and
    the gene summary to `fc_exons`.
- `biotype_genes` and `biotype_exons`: gene- and exon-level biotype composition
    PDF plots for RNA runs.
- `contig_assembly`: per-sample SPAdes or MEGAHIT assembly output when assembly
    is enabled (MEGAHIT is the configured default), under `contigs/<assembler>`.
- `bracken_summary`: combined Bracken abundance matrix with one row per sample
    and one column per detected taxon.
- `krona_html`: per-sample interactive Krona taxonomy reports.


## Errors
If you get a "Failed to pull singularity image" error on FH HPC, just rerun with -resume