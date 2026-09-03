#!/usr/bin/env bash
#SBATCH --job-name=vfdbmapping
#SBATCH --cpus-per-task=8
#SBATCH --output=vfdb_mapping_out/slurm_%j.out
#SBATCH --error=vfdb_mapping_out/slurm_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=idupless@fredhutch.org

ml SAMtools
ml BWA
ml picard
ml CD-HIT
ml Bowtie2


# run script:
# tail -n +2 samples.tsv | while IFS=$'\t' read -r sample pathogen dna_reads rna_reads vf_reference_fasta; do
#     sbatch vfdb_mapping.sh "$sample" "$pathogen" "$dna_reads" "$rna_reads" "$vf_reference_fasta"
# done

SAMPLE_NAME="$1"
PATHOGEN_NAME="$2"
DNA_TRIMMED_READS="$3"
RNA_TRIMMED_READS="$4"
VF_REFERENCE_FASTA="$5"

OUTPUT_DIR="./vfdb_mapping_out"
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

DNA_PREFIX="${DNA_OUT}/${SAMPLE_NAME}.${PATHOGEN_NAME}"
RNA_PREFIX="${RNA_OUT}/${SAMPLE_NAME}.${PATHOGEN_NAME}"

log "Sample:         $SAMPLE_NAME"
log "Pathogen:       $PATHOGEN_NAME"
log "VF reference:   $VF_REFERENCE_FASTA"
log "DNA reads:      $DNA_TRIMMED_READS"
log "RNA reads:      $RNA_TRIMMED_READS"
log "Threads:        $THREADS"

# ─────────────────────────────────────────────
# Validate inputs
# ─────────────────────────────────────────────
[[ ! -f "$VF_REFERENCE_FASTA" ]]  && err "VF reference FASTA not found: $VF_REFERENCE_FASTA"
[[ ! -f "$DNA_TRIMMED_READS" ]]   && err "DNA reads not found: $DNA_TRIMMED_READS"
[[ ! -f "$RNA_TRIMMED_READS" ]]   && err "RNA reads not found: $RNA_TRIMMED_READS"

# ─────────────────────────────────────────────
# Step 0a — Cluster sequences with CD-HIT
# ─────────────────────────────────────────────
CLUSTERED_FASTA="${VF_REFERENCE_FASTA%.*}.clustered.fa"

if [[ ! -f "$CLUSTERED_FASTA" ]]; then
    log "Clustering sequences with CD-HIT (c=0.9, s=0.9): $CLUSTERED_FASTA"
    cd-hit \
        -i "$VF_REFERENCE_FASTA" \
        -o "$CLUSTERED_FASTA" \
        -c 0.9 \
        -s 0.9 \
        -T "$THREADS" \
        -M 0 \
        >> "$LOG_FILE" 2>&1 \
        || err "CD-HIT clustering failed"
    log "Clustering complete: $CLUSTERED_FASTA"
else
    log "Clustered FASTA already exists, skipping: $CLUSTERED_FASTA"
fi

# Use clustered FASTA as input to header cleaning
VF_REFERENCE_FASTA="$CLUSTERED_FASTA"

# ─────────────────────────────────────────────
# Step 0 — Sanitize FASTA headers
# ─────────────────────────────────────────────
CLEAN_FASTA="${VF_REFERENCE_FASTA%.*}.clean.fa"

if [[ ! -f "$CLEAN_FASTA" ]]; then
    log "Creating clean FASTA with sanitized headers: $CLEAN_FASTA"
    awk '
        /^>/ {
            header = substr($0, 2)

            # Extract VFG ID (first token)
            vfg_id = header
            sub(/[^A-Za-z0-9].*/, "", vfg_id)

            # Extract gene name from second parenthetical group
            gene = ""
            tmp = header
            if (match(tmp, /\(gb\|[^)]*\)/)) {
                sub(/\(gb\|[^)]*\)/, "", tmp)
            }
            if (match(tmp, /\(([^)]+)\)/, arr)) {
                gene = arr[1]
            }

            # Extract strain from the LAST [...] group
            strain = ""
            tmp2 = header
            while (match(tmp2, /\[([^\]]+)\]/, arr)) {
                strain = arr[1]
                tmp2 = substr(tmp2, RSTART + RLENGTH)
            }

            # Build new header: gene_VFGid_strain
            if (gene != "")
                new_header = gene "_" vfg_id
            else
                new_header = vfg_id

            if (strain != "") new_header = new_header "_" strain

            # Sanitize: replace spaces and special chars with _
            gsub(/[^0-9A-Za-z_]/, "_", new_header)
            # Collapse multiple underscores
            gsub(/_+/, "_", new_header)
            # Strip trailing underscore
            sub(/_$/, "", new_header)

            print ">" new_header
            next
        }
        { print }
    ' "$VF_REFERENCE_FASTA" > "$CLEAN_FASTA" \
        || err "Failed to create clean FASTA"
    log "Clean FASTA ready: $CLEAN_FASTA"
else
    log "Clean FASTA already exists, skipping: $CLEAN_FASTA"
fi

# Use clean FASTA for all downstream steps
VF_REFERENCE_FASTA="$CLEAN_FASTA"
BT2_INDEX="${CLEAN_FASTA%.*}"


# ─────────────────────────────────────────────
# Step 1 — Index VF reference for BWA (DNA)
# ─────────────────────────────────────────────
if [[ ! -f "${VF_REFERENCE_FASTA}.bwt" ]]; then
    log "Indexing VF reference for BWA"
    bwa index "$VF_REFERENCE_FASTA" 2>> "$LOG_FILE" || err "BWA indexing failed"
    log "BWA indexing complete"
else
    log "BWA index already exists, skipping"
fi

# ─────────────────────────────────────────────
# Step 2 — Index VF reference for Bowtie2 (RNA)
# ─────────────────────────────────────────────


if [[ ! -f "${BT2_INDEX}.1.bt2" && ! -f "${BT2_INDEX}.1.bt2l" ]]; then
    log "Indexing VF reference for Bowtie2"
    bowtie2-build \
        --threads "$THREADS" \
        "$VF_REFERENCE_FASTA" \
        "$BT2_INDEX" \
        2>> "$LOG_FILE" || err "Bowtie2 indexing failed"
    log "Bowtie2 indexing complete"
else
    log "Bowtie2 index already exists, skipping"
fi

# ─────────────────────────────────────────────
# Step 3 — Map DNA reads (BWA aln/samse, ancient DNA params)
# ─────────────────────────────────────────────
DNA_BAM="${DNA_PREFIX}.bam"
DNA_DEDUP_BAM="${DNA_PREFIX}.dedup.bam"
DNA_DEDUP_METRICS="${DNA_PREFIX}.dedup_metrics.txt"

log "Mapping DNA reads with bwa aln/samse"
bwa aln \
    -n 0.1 \
    -l 1000 \
    -t "$THREADS" \
    "$VF_REFERENCE_FASTA" \
    "$DNA_TRIMMED_READS" \
    2>> "$LOG_FILE" \
| bwa samse \
    -r "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}" \
    "$VF_REFERENCE_FASTA" \
    - \
    "$DNA_TRIMMED_READS" \
    2>> "$LOG_FILE" \
| samtools view -F 4 -q 20 -bSh \
| samtools sort -@ "$THREADS" -o "$DNA_BAM" \
    2>> "$LOG_FILE" \
|| err "DNA mapping failed"

samtools index "$DNA_BAM" 2>> "$LOG_FILE"
log "DNA BAM ready: $DNA_BAM"

# ─────────────────────────────────────────────
# Step 4 — Dedup DNA BAM
# ─────────────────────────────────────────────
log "Removing duplicates from DNA BAM with Picard"
java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
    -I "$DNA_BAM" \
    -O "$DNA_DEDUP_BAM" \
    -M "$DNA_DEDUP_METRICS" \
    --REMOVE_DUPLICATES true \
    2>> "$LOG_FILE" \
|| err "Picard dedup failed on DNA BAM"

samtools index "$DNA_DEDUP_BAM" 2>> "$LOG_FILE"
log "DNA dedup BAM ready: $DNA_DEDUP_BAM"

# ─────────────────────────────────────────────
# Step 5 — Map RNA reads (Bowtie2, single-end)
# ─────────────────────────────────────────────
RNA_BAM="${RNA_PREFIX}.bam"
RNA_DEDUP_BAM="${RNA_PREFIX}.dedup.bam"
RNA_DEDUP_METRICS="${RNA_PREFIX}.dedup_metrics.txt"
RNA_SAM="${RNA_PREFIX}.sam"
RNA_RG_SAM="${RNA_PREFIX}.rg.sam"

log "Mapping RNA reads with Bowtie2"
bowtie2 \
    -x "$BT2_INDEX" \
    -U "$RNA_TRIMMED_READS" \
    --end-to-end --sensitive \
    --no-unal \
    --threads "$THREADS" \
    -S "$RNA_SAM" \
    2>> "$LOG_FILE" \
    || err "RNA mapping failed"

log "Adding read group to RNA SAM"
samtools addreplacerg \
    -r "ID:${SAMPLE_NAME}" \
    -r "SM:${SAMPLE_NAME}" \
    -o "$RNA_RG_SAM" \
    "$RNA_SAM" \
    2>> "$LOG_FILE" \
    || err "samtools addreplacerg failed"

log "Filtering and sorting RNA BAM"
samtools view -bh -F 4 -q 20 "$RNA_RG_SAM" \
    | samtools sort -@ "$THREADS" -o "$RNA_BAM" \
    2>> "$LOG_FILE" \
    || err "samtools view/sort failed on RNA"

rm -f "$RNA_SAM" "$RNA_RG_SAM"
log "Removed intermediate SAM files"

samtools index "$RNA_BAM" 2>> "$LOG_FILE"
log "RNA BAM ready: $RNA_BAM"

# ─────────────────────────────────────────────
# Step 6 — Dedup RNA BAM
# ─────────────────────────────────────────────
log "Removing duplicates from RNA BAM with Picard"
java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
    -I "$RNA_BAM" \
    -O "$RNA_DEDUP_BAM" \
    -M "$RNA_DEDUP_METRICS" \
    --REMOVE_DUPLICATES true \
    2>> "$LOG_FILE" \
|| err "Picard dedup failed on RNA BAM"

samtools index "$RNA_DEDUP_BAM" 2>> "$LOG_FILE"
log "RNA dedup BAM ready: $RNA_DEDUP_BAM"

log "Pipeline complete for $SAMPLE_NAME / $PATHOGEN_NAME"
log "  DNA dedup BAM: $DNA_DEDUP_BAM"
log "  RNA dedup BAM: $RNA_DEDUP_BAM"


# ─────────────────────────────────────────────
# Step 7 — Coverage summary table (DNA + RNA)
# ─────────────────────────────────────────────
COVERAGE_TABLE="${SAMPLE_OUT}/${SAMPLE_NAME}.${PATHOGEN_NAME}.coverage_summary.tsv"

log "Running samtools coverage on DNA dedup BAM"

# Get sequence IDs with at least one read mapped in DNA (numreads col 4 > 0)
DNA_COVERED_SEQS=$(samtools coverage "$DNA_DEDUP_BAM" 2>> "$LOG_FILE" \
    | awk -F'\t' 'NR>1 && $4 > 0 { print $1 }')

if [[ -z "$DNA_COVERED_SEQS" ]]; then
    log "No sequences with DNA coverage found — skipping coverage table"
else
    log "Building coverage table for sequences with DNA reads mapped"

    # Write header
    printf "seq_id\tfull_header\tdna_coverage\tdna_meandepth\trna_coverage\trna_meandepth\n" \
        > "$COVERAGE_TABLE"

    # Pull full headers from reference FASTA into an associative array
    declare -A FASTA_HEADERS
    while IFS= read -r line; do
        if [[ "$line" == ">"* ]]; then
            seq_id="${line:1}"
            seq_id="${seq_id%% *}"
            full_header="${line:1}"
            FASTA_HEADERS["$seq_id"]="$full_header"
        fi
    done < "$VF_REFERENCE_FASTA"

    # Get RNA coverage stats into associative arrays keyed by seq_id
    # cols: #rname(1) startpos(2) endpos(3) numreads(4) covbases(5) coverage(6) meandepth(7)
    declare -A RNA_COV RNA_DEPTH
    while IFS=$'\t' read -r rseq _ _ _ _ rcov rdepth _; do
        RNA_COV["$rseq"]="$rcov"
        RNA_DEPTH["$rseq"]="$rdepth"
    done < <(samtools coverage "$RNA_DEDUP_BAM" 2>> "$LOG_FILE" | tail -n +2)

    # Get DNA coverage stats for covered sequences and build table
    # Filter to numreads(4) > 0, extract coverage(6) and meandepth(7)
    while IFS=$'\t' read -r seq _ _ _ _ dcov ddepth _; do
        full_header="${FASTA_HEADERS[$seq]:-$seq}"
        rcov="${RNA_COV[$seq]:-0}"
        rdepth="${RNA_DEPTH[$seq]:-0}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$seq" "$full_header" "$dcov" "$ddepth" "$rcov" "$rdepth"
    done < <(samtools coverage "$DNA_DEDUP_BAM" 2>> "$LOG_FILE" \
        | awk -F'\t' 'NR>1 && $4 > 0') \
        >> "$COVERAGE_TABLE"

    log "Coverage table written: $COVERAGE_TABLE"
fi