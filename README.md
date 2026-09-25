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
an `R2` column for paired-end data. All rows must use the same single-end or
paired-end layout.

`PE`: set to `true` when the samplesheet contains paired-end reads.

`sample_type`: `DNA` or `RNA`.

`merge_inputs`: set to `true` to allow duplicate `sample` names in the
samplesheet and merge their inputs after mapping.

### Common processing parameters

`trim_length`: minimum read length retained by Trim Galore.

`trim_quality`: quality threshold used by Trim Galore.

`host_reference`: host reference used for BWA mapping/indexing and downstream mapDamage and dsRNA analysis.

`hg_gtf`: host genome GTF used for feature counting.

`dsRNA_bed`: dsRNA annotation BED file.

`strand`: feature-counting and dsRNA strandedness setting.

`chunk_size`: number of reads per mapping chunk.

`seed_length`, `map_missing_prob`, `gap_fraction`: BWA alignment settings.

`map_quality`: minimum mapping quality used when merging BAM files.

`publish_host_bam`: set to `true` to publish deduplicated host BAM files.

### RNA parameters

`ribo_bt2_ref`: Bowtie2 index prefix for the ribosomal RNA reference.

`rrna_bed`: ribosomal RNA annotation BED file.

`RNA_host_mapper`: one of `STAR`, `bowtie2`, or `BWA`.

For `RNA_host_mapper=STAR`, provide `hg_star_ref` and `STAR_mapq`.
For `RNA_host_mapper=bowtie2`, provide `bt2_ref` and `bt2_mapq`.
For `RNA_host_mapper=BWA`, provide `bwa_mapq`.

## Outputs

The workflow publishes the following named outputs under the output directory
specified with `-o`:

- `fastqc_raw`: FastQC results for the raw reads.
- `fastqc_trimmed`: FastQC results for the trimmed reads.
- `trimmed_reads`: Trim Galore reads retained after adapter and quality trimming.
- `host_bam`: host-mapped BAM files before duplicate removal.
- `host_bam_deduped`: deduplicated host BAM files when `publish_host_bam=true`.
- `mapdamage`: mapDamage results for deduplicated host mappings.
- `read_length_distribution`: read-length plots for deduplicated host mappings.
- `high_complexity_reads`: reads retained after read-complexity filtering.
- `summary_stats`: `read_count_report.tsv` containing read counts by sample and processing stage.
- `multimapped_bam`: multimapped host BAM files for RNA mapping methods that produce them.
- `dsRNA_dir`: dsRNA-associated mapping and mapDamage output directories.
- `ssRNA_dir`: ssRNA-associated mapDamage output directories.
- `fc_genes`: combined feature-count summary CSV published under the `fc_genes` name.
- `fc_exons`: combined feature-count summary CSV published under the `fc_exons` name.
- `biotype_genes`: gene-level biotype composition PDF plot.
- `biotype_exons`: exon-level biotype composition PDF plot.


## Errors
If you get a "Failed to pull singularity image" error on FH HPC, just rerun with -resume