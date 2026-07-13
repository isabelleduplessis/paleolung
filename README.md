# paleolung
Pipeline for ancient pathogen metagenomics analysis

Pipeline for ancient RNA authentication and pathogen metagenomics


## Installing and setting up Nextflow version 26
```
cd ~/.local/bin
export NXF_VER=26.04.1
curl -s https://get.nextflow.io | bash
chmod +x nextflow
export PATH="$PATH:$HOME/.local/bin"
```
Add the export path line to ~/.bash_profile and ~/.bashrc as well if not already
```
echo 'export PATH="$PATH:$HOME/.local/bin"' >> ~/.bashrc && echo 'export PATH="$PATH:$HOME/.local/bin"' >> ~/.bash_profile
```


To run nextflow v26, Java version 17 at minimum is needed.


## Set up configuration
Example gizmo config file: nextflow.gizmo.config
```
workDir = '/path/to/tmp/dir/work/'

apptainer {
    enabled = true
    autoMounts = true
    cacheDir = '/path/to/tmp/dir/cache/'
    runOptions = '--containall -B \$TMPDIR'
}

env {
    TMPDIR="\$TMPDIR"
}

process {
    executor = 'slurm'
    queue = 'campus-new'
    // optional fields: 
    // errorStrategy = 'retry'
    // maxRetries = 3
}

// hard copy all result files, so that if tmp dir is cleared you still have results. Does double amount of storage being used until tmp dir clears
workflow {
    output.mode = 'copy'
}
```



## Run 

Example nf-run script on HPC cluster after installation:

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
    -output-dir results/test_results \
    --samplesheet './data/test/samplesheet.csv'
```

## Required input parameters

(call with --parameter_name parameter_value)

`trim_length`: minimum trim length for fastp

`ribo_bt2_ref`: bowtie2 index for human ribosomal RNA reference

`ribo_bwa_ref`: human ribosomal RNA reference fasta

`human_bt2_ref`: bowtie2 index for human genome reference

`human_bwa_ref`: human genome reference fasta

`kraken_db`: kraken standard database directory

`exon_bed`: exon bed6 file

`intron_bed`: intron bed6 file

`gene_annotation`: human genome gtf, genes only

`dsRNAome`: dsRNAome arms annotation


## Provided input files

`exon_bed`: `./assets/exons.bed`

`intron_bed`: `./assets/introns.bed`

`gene_annotation`: `./assets/hg38_v49_genes.gtf`

`dsRNAome`: `./assets/dsRNAome.bed`

