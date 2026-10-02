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


To run nextflow v26, Java version 17 or later is needed, which we will load when running the pipeline. Nextflow 26 uses v2 syntax parser which is strict and quite different than previous versions. Extra care is needed when writing code with v26 as most existing documentation and training data uses the v1 parser


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
    -o ${nf_outdir} --RNA_host_mapper bowtie2 -resume --trim_length 30

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

### Read input

`samplesheet`: CSV file with a header containing `sample` and `R1` columns. Add
an `R2` column when `read_type` is `pe` or `merge_pe`. All rows must use the
same read layout.

`read_type`: select `se` for single-end reads, `pe` for paired-end reads, or
`merge_pe` to merge paired-end reads with fastp. Defaults to `se`. Both `pe`
and `merge_pe` require `R2` values in the samplesheet. With `merge_pe`, the
merged reads are treated as single-end by all downstream steps. For example,
use `--read_type merge_pe --trimmer fastp` to enable fastp merging.

`sample_type`: `DNA` or `RNA`.

Every samplesheet row must have a unique `sample` name. To process separate
libraries or read sets independently, give each row its own sample name.

`skip_host_mapping`: set to `true` when `R1` (and `R2` for paired-end data) already points
to high-complexity reads produced by a previous host-mapping run. This skips
preprocessing and host mapping and runs taxonomic profiling directly on those
reads. The samplesheet format is unchanged.

### Common processing parameters

`trim_length`: minimum read length retained by the selected trimmer.

`trim_quality`: quality threshold used by the selected trimmer.

`trimmer`: read trimmer to use; set to `fastp` or `trim_galore`. Defaults to
`trim_galore` when omitted.

`min_overlap`: minimum paired-read overlap length used by fastp when
`read_type=merge_pe`. Defaults to `30` bases. The merged FASTQ also includes
reads that could not be paired into an overlap, as singletons. `merge_pe`
cannot be used with `--skip_host_mapping`.

`host_reference`: host reference used for BWA mapping/indexing and downstream mapDamage and dsRNA analysis.

`hg_gtf`: host genome GTF used for feature counting.

`dsRNA_bed`: dsRNA annotation BED file.

`strand`: feature-counting and dsRNA strandedness setting.

`chunk_size`: number of reads per mapping chunk.

`seed_length`, `map_missing_prob`, `gap_fraction`: BWA alignment settings.

`map_quality`: minimum mapping quality used when merging BAM files.

`publish_host_bam`: set to `true` to publish deduplicated host BAM files.

`assembly_tool`: optionally set to `spades` or `megahit` to assemble the
high-complexity reads. Assembly is skipped when this parameter is omitted.
The selected assembler supports both single-end and paired-end inputs, as
selected by `PE` and the samplesheet layout.

### RNA parameters

`ribo_bt2_ref`: Bowtie2 index prefix for the ribosomal RNA reference; required
when `ncRNA_mapper=bowtie2`.

`trna_bt2_ref`: Bowtie2 index prefix for the tRNA reference; required when
`ncRNA_mapper=bowtie2`. Reads that do not map to the ribosomal reference are
mapped to tRNA before host mapping.

`ribo_bwa_ref` and `trna_bwa_ref`: BWA reference/index prefixes for ribosomal
RNA and tRNA; both are required when `ncRNA_mapper=BWA`. Supply prebuilt BWA
index sidecars for each prefix.

`ncRNA_mapper`: mapper used for ribosomal RNA and tRNA; choose `bowtie2` or
`BWA`. Defaults to `bowtie2`.

`rrna_bed`: ribosomal RNA annotation BED file.

`filter_ribo`: set to `true` to remove host-mapped reads overlapping the
rRNA annotation and emit the intersecting rRNA BAM separately. Defaults to
`true`; set `--filter_ribo false` to skip this post-mapping filter. Reads
unmapped from tRNA mapping still proceed directly to host mapping either way.

`RNA_host_mapper`: one of `STAR`, `bowtie2`, or `BWA`.

For `RNA_host_mapper=STAR`, provide `hg_star_ref` and `STAR_mapq`.
For `RNA_host_mapper=bowtie2`, provide `bt2_ref` and `bt2_mapq`.
For `RNA_host_mapper=BWA`, provide `bwa_mapq`.

### Taxonomic profiling parameters

By default, taxonomic profiling runs Kraken2/Bracken and Krona on the
high-complexity reads remaining after host mapping. With `--skip_host_mapping`, it instead
uses the high-complexity reads supplied in the samplesheet. Provide
`kraken2_db`, a directory containing the Kraken2 database and Bracken database
files. `bracken_read_len` defaults to `100`, `bracken_tax_level` defaults to
`S`, and `krona_min_fraction` defaults to `0.5` percent for the Krona input
filter.

## Outputs

The workflow publishes the following named outputs under the output directory
specified with `-o`:

- `fastqc_raw`: FastQC results for the raw reads.
- `fastqc_trimmed`: FastQC results for the trimmed reads.
- `trimmed_reads`: reads retained after trimming with the selected `trimmer`.
- `host_bam`: host-mapped BAM files before duplicate removal.
- `host_bam_deduped`: deduplicated host BAM files when `publish_host_bam=true`.
- `mapdamage`: mapDamage results for deduplicated host mappings.
- `read_length_distribution`: read-length plots for deduplicated host mappings.
- `high_complexity_reads`: reads retained after read-complexity filtering, or
    the input reads supplied when running with `--skip_host_mapping`.
- `contig_assembly`: per-sample SPAdes or MEGAHIT assembly directory, published
    under `assemblies/<assembler>/<sample>` when `assembly_tool` is selected.
- `summary_stats`: `read_count_report.tsv` containing read counts by sample and processing stage. RNA summaries include `ribo_mapped` and `trna_mapped`; `ribo_intersect` is included only when `--filter_ribo true`.
- `multimapped_bam`: multimapped host BAM files for RNA mapping methods that produce them.
- `dsRNA_dir`: dsRNA-associated mapping and mapDamage output directories.
- `ssRNA_dir`: ssRNA-associated mapDamage output directories.
- `fc_genes`: combined feature-count summary CSV published under the `fc_genes` name.
- `fc_exons`: combined feature-count summary CSV published under the `fc_exons` name.
- `biotype_genes`: gene-level biotype composition PDF plot.
- `biotype_exons`: exon-level biotype composition PDF plot.
- `bracken_summary`: combined Bracken abundance matrix with one row per sample
    and one column per detected taxon.
- `krona_html`: per-sample interactive Krona taxonomy reports.


## Errors
If you get a "Failed to pull singularity image" error on FH HPC, just rerun with -resume