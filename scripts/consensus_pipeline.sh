#!/usr/bin/env bash
#SBATCH --job-name=consensus_genomes
#SBATCH --cpus-per-task=8
#SBATCH --output=pipeline_out/slurm_%j.out
#SBATCH --error=pipeline_out/slurm_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=idupless@fredhutch.org

ml SAMtools
ml BWA
ml picard
ml BCFtools
ml BEDTools
ml Bowtie2

SAMPLE_NAME="$1"
PATHOGEN_NAME="$2"
REFERENCE_FASTA="$3"
REFERENCE_GFF="$4"
DNA_TRIMMED_READS="$5"
RNA_TRIMMED_READS="$6"

OUTPUT_DIR="./pipeline_out"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ─────────────────────────────────────────────
# Output dirs and logging
# ─────────────────────────────────────────────
SAMPLE_OUT="$OUTPUT_DIR/$SAMPLE_NAME/$PATHOGEN_NAME"
DNA_OUT="$SAMPLE_OUT/dna"
RNA_OUT="$SAMPLE_OUT/rna"
mkdir -p "$DNA_OUT" "$RNA_OUT"

LOG_FILE="$SAMPLE_OUT/pipeline.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO]  $*" | tee -a "$LOG_FILE"; }
err() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" | tee -a "$LOG_FILE"; exit 1; }

PREFIX="${DNA_OUT}/${SAMPLE_NAME}.${PATHOGEN_NAME}"
RNA_PREFIX="${RNA_OUT}/${SAMPLE_NAME}.${PATHOGEN_NAME}"

log "Sample:    $SAMPLE_NAME"
log "Pathogen:  $PATHOGEN_NAME"
log "Reference: $REFERENCE_FASTA"
log "DNA reads: $DNA_TRIMMED_READS"
log "RNA reads: $RNA_TRIMMED_READS"
log "Threads:   $THREADS"
log "Prefix:    $PREFIX"

# ─────────────────────────────────────────────
# Validate inputs
# ─────────────────────────────────────────────
[[ ! -f "$REFERENCE_FASTA" ]]   && err "FASTA not found: $REFERENCE_FASTA"
[[ ! -f "$DNA_TRIMMED_READS" ]] && err "DNA reads not found: $DNA_TRIMMED_READS"
[[ ! -f "$RNA_TRIMMED_READS" ]] && err "RNA reads not found: $RNA_TRIMMED_READS"

# ─────────────────────────────────────────────
# Step 1 — Index reference
# ─────────────────────────────────────────────
if [[ ! -f "${REFERENCE_FASTA}.bwt" ]]; then
    log "Indexing reference for BWA"
    bwa index "$REFERENCE_FASTA" 2>> "$LOG_FILE"
    log "BWA indexing complete"
else
    log "BWA index already exists, skipping"
fi

# ─────────────────────────────────────────────
# Step 2 — Map DNA reads (single-end, ancient DNA params)
# ─────────────────────────────────────────────
BAM="${PREFIX}.bam"
DEDUP_BAM="${PREFIX}.dedup.bam"
DEDUP_METRICS="${PREFIX}.dedup_metrics.txt"
RESCALED_BAM="${PREFIX}.dedup.rescaled.bam"
MAPDAMAGE_DIR="${DNA_OUT}/${SAMPLE_NAME}.${PATHOGEN_NAME}.mapdamage"

# log "Mapping DNA reads with bwa aln/samse"
# bwa aln \
#     -n 0.1 \
#     -l 1000 \
#     -t "$THREADS" \
#     "$REFERENCE_FASTA" \
#     "$DNA_TRIMMED_READS" \
#     2>> "$LOG_FILE" \
# | bwa samse \
#     -r "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}" \
#     "$REFERENCE_FASTA" \
#     - \
#     "$DNA_TRIMMED_READS" \
#     2>> "$LOG_FILE" \
# | samtools view -F 4 -q 20 -bSh \
# | samtools sort -@ "$THREADS" -o "$BAM" \
#     2>> "$LOG_FILE"

# samtools index "$BAM" 2>> "$LOG_FILE"
# log "BAM ready: $BAM"

# ─────────────────────────────────────────────
# Step 3 — Mark and remove duplicates with Picard
# ─────────────────────────────────────────────
# log "Removing duplicates with Picard"
# java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
#     -I "$BAM" \
#     -O "$DEDUP_BAM" \
#     -M "$DEDUP_METRICS" \
#     --REMOVE_DUPLICATES true \
#     2>> "$LOG_FILE"

# samtools index "$DEDUP_BAM" 2>> "$LOG_FILE"
# log "Dedup BAM ready: $DEDUP_BAM"

# ─────────────────────────────────────────────
# Step 4 — mapDamage rescaling
# ─────────────────────────────────────────────
log "Running mapDamage"
set --
ml Miniforge3
source activate
conda activate mapdamage

# mapDamage \
#     -i "$DEDUP_BAM" \
#     -r "$REFERENCE_FASTA" \
#     -d "$MAPDAMAGE_DIR" \
#     --rescale \
#     --rescale-out="$RESCALED_BAM" --merge-reference-sequences \
#     2>> "$LOG_FILE"

# samtools index "$RESCALED_BAM" 2>> "$LOG_FILE"
# log "Rescaled BAM ready: $RESCALED_BAM"

conda deactivate

# ─────────────────────────────────────────────
# Step 5 — Consensus genome
# ─────────────────────────────────────────────
conda activate bgzip

DRAFT_VCF="${PREFIX}.draft.vcf"
DRAFT_VCF_GZ="${DRAFT_VCF}.gz"
NOCOV_BED="${PREFIX}.nocov.bed"
CONSENSUS_FA="${PREFIX}.draftgenome.fa"

# log "Calling variants with bcftools mpileup"
# bcftools mpileup -Q 20 \
#     -Ou \
#     -f "$REFERENCE_FASTA" \
#     "$RESCALED_BAM" \
#     2>> "$LOG_FILE" \
# | bcftools call \
#     -mv \
#     -Ov \
#     -o "$DRAFT_VCF" \
#     2>> "$LOG_FILE"

# bgzip "$DRAFT_VCF" 2>> "$LOG_FILE"
# bcftools index "$DRAFT_VCF_GZ" 2>> "$LOG_FILE"

# log "Generating no-coverage mask"
# bedtools genomecov \
#     -ibam "$RESCALED_BAM" \
#     -bga \
#     2>> "$LOG_FILE" \
# | awk '$4 < 1' > "$NOCOV_BED"

# log "Generating consensus"
# bcftools consensus \
#     "$DRAFT_VCF_GZ" \
#     -f "$REFERENCE_FASTA" \
#     -m "$NOCOV_BED" \
#     > "$CONSENSUS_FA" \
#     2>> "$LOG_FILE"

# log "Consensus written: $CONSENSUS_FA"
conda deactivate

# ─────────────────────────────────────────────
# Step 6 — Index consensus for Bowtie2 (RNA)
# ─────────────────────────────────────────────
BT2_INDEX="${CONSENSUS_FA%.*}"

# if [[ ! -f "${BT2_INDEX}.1.bt2" && ! -f "${BT2_INDEX}.1.bt2l" ]]; then
#     log "Indexing consensus FASTA for Bowtie2"
#     bowtie2-build \
#         --threads "$THREADS" \
#         "$CONSENSUS_FA" \
#         "$BT2_INDEX" \
#         2>> "$LOG_FILE" || err "Bowtie2 indexing failed"
#     log "Bowtie2 indexing complete"
# else
#     log "Bowtie2 index already exists, skipping"
# fi

# ─────────────────────────────────────────────
# Step 7 — Map RNA reads to consensus (Bowtie2, single-end)
# ─────────────────────────────────────────────
RNA_BAM="${RNA_PREFIX}.bam"
RNA_DEDUP_BAM="${RNA_PREFIX}.dedup.bam"
RNA_DEDUP_METRICS="${RNA_PREFIX}.dedup_metrics.txt"
RNA_SAM="${RNA_PREFIX}.sam"
RNA_RG_SAM="${RNA_PREFIX}.rg.sam"

# log "Mapping RNA reads with Bowtie2"
# bowtie2 \
#     -x "$BT2_INDEX" \
#     -U "$RNA_TRIMMED_READS" \
#     --end-to-end --sensitive \
#     --no-unal \
#     --threads "$THREADS" \
#     -S "$RNA_SAM" \
#     2>> "$LOG_FILE" \
#     || err "RNA mapping failed"

# log "Adding read group to RNA SAM"
# samtools addreplacerg \
#     -r "ID:${SAMPLE_NAME}" \
#     -r "SM:${SAMPLE_NAME}" \
#     -o "$RNA_RG_SAM" \
#     "$RNA_SAM" \
#     2>> "$LOG_FILE" \
#     || err "samtools addreplacerg failed"

# log "Filtering and sorting RNA BAM"
# samtools view -bh -F 4 -q 20 "$RNA_RG_SAM" \
#     | samtools sort -@ "$THREADS" -o "$RNA_BAM" \
#     2>> "$LOG_FILE" \
#     || err "samtools view/sort failed on RNA"

# rm -f "$RNA_SAM" "$RNA_RG_SAM"
# log "Removed intermediate SAM files"

# samtools index "$RNA_BAM" 2>> "$LOG_FILE"
# log "RNA BAM ready: $RNA_BAM"

# ─────────────────────────────────────────────
# Step 8 — Dedup RNA BAM
# ─────────────────────────────────────────────
# log "Removing duplicates from RNA BAM with Picard"
# java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
#     -I "$RNA_BAM" \
#     -O "$RNA_DEDUP_BAM" \
#     -M "$RNA_DEDUP_METRICS" \
#     --REMOVE_DUPLICATES true \
#     2>> "$LOG_FILE" \
#     || err "Picard dedup failed on RNA BAM"

# samtools index "$RNA_DEDUP_BAM" 2>> "$LOG_FILE"
# log "RNA dedup BAM ready: $RNA_DEDUP_BAM"


# ─────────────────────────────────────────────
# Step 9 — Liftoff: annotate consensus from reference GFF
# ─────────────────────────────────────────────
LIFTED_GFF="${PREFIX}.consensus.gff3"
UNMAPPED_FEATURES="${PREFIX}.unmapped_features.txt"

log "Running Liftoff to transfer annotations to consensus genome"
conda activate liftoff   # adjust to your env name

liftoff \
    -g "$REFERENCE_GFF" \
    -o "$LIFTED_GFF" \
    -u "$UNMAPPED_FEATURES" \
    "$CONSENSUS_FA" \
    "$REFERENCE_FASTA" \
    2>> "$LOG_FILE" \
    || err "Liftoff failed"

log "Lifted annotation written: $LIFTED_GFF"
log "Unmapped features written: $UNMAPPED_FEATURES"

conda deactivate

# ─────────────────────────────────────────────
# Step 10 — featureCounts: quantify RNA reads
# ─────────────────────────────────────────────
COUNTS_OUT="${RNA_OUT}/${SAMPLE_NAME}.${PATHOGEN_NAME}.featurecounts.txt"

log "Running featureCounts on RNA dedup BAM"
conda activate subread   

featureCounts \
    -a "$LIFTED_GFF" \
    -o "$COUNTS_OUT" \
    -F GFF \
    -t CDS \
    -g locus_tag \
    -s 1 \
    "$RNA_DEDUP_BAM" \
    2>> "$LOG_FILE" \
    || err "featureCounts failed"

log "featureCounts output: $COUNTS_OUT"

conda deactivate

# ─────────────────────────────────────────────
# Build locus_tag → gene name mapping from lifted GFF
# ─────────────────────────────────────────────
LOCUS_MAP="${RNA_OUT}/${SAMPLE_NAME}.${PATHOGEN_NAME}.locus_to_name.tsv"

log "Building locus_tag to gene name mapping"
grep -v "^#" "$LIFTED_GFF" \
    | awk -F'\t' '$3=="gene" {
        match($9, /locus_tag=([^;]+)/, l)
        match($9, /Name=([^;]+)/, n)
        if (l[1] != "") print l[1] "\t" n[1]
    }' \
    > "$LOCUS_MAP" \
    || err "Failed to build locus_tag mapping"

log "Locus mapping written: $LOCUS_MAP"

# ─────────────────────────────────────────────
# Merge mapping into featureCounts output
# ─────────────────────────────────────────────
COUNTS_NAMED="${COUNTS_OUT%.txt}.named.txt"
ml fhPython

log "Merging gene names into featureCounts output"
python3 - <<EOF
import pandas as pd

counts = pd.read_csv("$COUNTS_OUT", sep="\t", comment="#")
mapping = pd.read_csv("$LOCUS_MAP", sep="\t", header=None, names=["locus_tag", "gene_name"])

merged = counts.merge(mapping, left_on="Geneid", right_on="locus_tag", how="left")
merged.insert(1, "gene_name", merged.pop("gene_name"))
merged.drop(columns=["locus_tag"], inplace=True)

merged.to_csv("$COUNTS_NAMED", sep="\t", index=False)
print(f"Written: $COUNTS_NAMED")
EOF

log "Named counts written: $COUNTS_NAMED"

log "Done. Outputs in: $SAMPLE_OUT"


