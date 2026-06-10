#!/usr/bin/env bash
# =============================================================================
# Genome Assembly Pipeline: FASTP → bwa-mem2 → Shovill/SPAdes (de novo)
#                          → RagTag scaffold → Pilon → Coverage map → Report
# Platform: Docker (mambaorg/micromamba base)
# Checkpoint: each step skipped if output already exists (re-run safe)
# Usage: docker run --rm \
#          -v /your/data:/data -v /your/results:/results -v /your/refs:/ref \
#          seqforge \
#          -1 /data/R1.fastq.gz -2 /data/R2.fastq.gz \
#          -r /ref/reference.fasta \
#          -o /results [-t threads] [-p prefix] [-g genome_size]
# =============================================================================

set -euo pipefail

THREADS=6
PREFIX="sample"
MIN_MAP_QUAL=10
GENOME_SIZE=""  # required; must provide with -g
SHOVILL_DEPTH=150   # subsample to this depth for assembly (0 = use all reads)
REFERENCE=""

log_step() {
    echo ""
    echo "╔══════════════════════════════════════════════════╗"
    printf  "║  %-48s║\n" "$1"
    echo "╚══════════════════════════════════════════════════╝"
}
log_info() { echo "    [INFO]  $1"; }
log_ok()   { echo "    [DONE]  $1"; }
log_file() { echo "    [FILE]  $1"; }
log_skip() { echo "    [SKIP]  $1 — output exists, skipping"; }

usage() {
    echo ""
    echo "Usage: $0 -1 READ1 -2 READ2 -r REFERENCE -g GENOME_SIZE -o OUTDIR [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  -1  PATH    R1 FASTQ file (gzipped)"
    echo "  -2  PATH    R2 FASTQ file (gzipped)"
    echo "  -r  PATH    Reference genome FASTA"
    echo "  -g  INT     Genome size in bp (estimated size of your target organism)"
    echo "  -o  PATH    Output directory"
    echo ""
    echo "Optional:"
    echo "  -t  INT     CPU threads        (default: 6)"
    echo "  -p  NAME    Sample prefix      (default: sample)"
    echo "  -d  INT     Assembly depth x   (default: 150)"
    echo "  -h          Show this help"
    echo ""
    echo "Genome size examples:"
    echo "  SARS-CoV-2: 29903"
    echo "  E. coli:    4600000"
    echo "  Arabidopsis: 135000000"
    echo "  Human:      3200000000"
    echo "  Wheat:      17000000000"
    echo ""
    exit 1
}

READ1=""
READ2=""
OUTDIR=""

while getopts "1:2:r:o:t:p:g:d:h" opt; do
    case $opt in
        1) READ1="$OPTARG" ;;
        2) READ2="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        p) PREFIX="$OPTARG" ;;
        g) GENOME_SIZE="$OPTARG" ;;
        d) SHOVILL_DEPTH="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

# Validate required arguments
for var in READ1 READ2 REFERENCE GENOME_SIZE OUTDIR; do
    if [ -z "${!var:-}" ]; then
        echo "ERROR: Missing required argument for $var (use -h for help)"
        usage
    fi
done

# Validate reference exists
if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference file not found: $REFERENCE"
    exit 1
fi

# ---------------------------------------------------------------------------
# Directory structure
# ---------------------------------------------------------------------------
FASTQC_RAW_DIR="${OUTDIR}/00_fastqc_raw"
TRIM_DIR="${OUTDIR}/01_trimmed"
MAP_DIR="${OUTDIR}/02_mapped"
MAPPED_FQ_DIR="${OUTDIR}/03_mapped_fastq"
DENOVO_DIR="${OUTDIR}/04_shovill_denovo"
RAGTAG_DIR="${OUTDIR}/05_ragtag_scaffold"
FINAL_DIR="${OUTDIR}/06_final_consensus"
COVERAGE_DIR="${OUTDIR}/07_coverage_map"
FASTQC_FINAL_DIR="${OUTDIR}/08_fastqc_final"
LOG_DIR="${OUTDIR}/logs"

mkdir -p "$FASTQC_RAW_DIR" "$TRIM_DIR" "$MAP_DIR" "$MAPPED_FQ_DIR" "$DENOVO_DIR" \
         "$RAGTAG_DIR" "$FINAL_DIR" "$COVERAGE_DIR" \
         "$FASTQC_FINAL_DIR" "$LOG_DIR"

# ---------------------------------------------------------------------------
# Cleanup: remove old folders from previous pipeline versions
# ---------------------------------------------------------------------------
for OLD_DIR in \
    "${OUTDIR}/07_annotation" \
    "${OUTDIR}/08_coverage_map" \
; do
    if [ -d "$OLD_DIR" ] && [ "$OLD_DIR" != "$COVERAGE_DIR" ]; then
        echo "    [CLEAN] Removing old folder: $OLD_DIR"
        rm -rf "$OLD_DIR"
    fi
done

# Checkpoint file paths
TRIM_R1="${TRIM_DIR}/${PREFIX}_R1_trimmed.fastq.gz"
TRIM_R2="${TRIM_DIR}/${PREFIX}_R2_trimmed.fastq.gz"
SORTED_BAM="${MAP_DIR}/${PREFIX}_sorted.bam"
MAPPED_R1="${MAPPED_FQ_DIR}/${PREFIX}_mapped_R1.fastq.gz"
MAPPED_R2="${MAPPED_FQ_DIR}/${PREFIX}_mapped_R2.fastq.gz"
DENOVO_CONTIGS="${DENOVO_DIR}/${PREFIX}_contigs.fasta"
SCAFFOLD="${RAGTAG_DIR}/ragtag.scaffold.fasta"
SCAFFOLD_BAM="${RAGTAG_DIR}/${PREFIX}_scaffold_sorted.bam"
FINAL_CONSENSUS="${FINAL_DIR}/${PREFIX}_final_consensus.fasta"
FINAL_BAM="${COVERAGE_DIR}/${PREFIX}_final_sorted.bam"
DEPTH_TXT="${COVERAGE_DIR}/${PREFIX}_depth.txt"
COVERAGE_PLOT="${COVERAGE_DIR}/${PREFIX}_coverage"
REPORT_FILE="${OUTDIR}/${PREFIX}_assembly_report.html"

echo ""
echo "╔══════════════════════════════════════════════════╗"
echo "║         Genome Assembly Pipeline (Docker)        ║"
echo "╠══════════════════════════════════════════════════╣"
printf "║  Sample      : %-33s║\n" "$PREFIX"
printf "║  Threads     : %-33s║\n" "$THREADS"
printf "║  Genome size : %-33s║\n" "${GENOME_SIZE} bp"
printf "║  Assembly depth: %-31s║\n" "${SHOVILL_DEPTH}×"
echo "╠══════════════════════════════════════════════════╣"
printf "║  Reference   : %-33s║\n" "$(basename "$REFERENCE")"
printf "║  OUTPUT DIR  : %-33s║\n" "$OUTDIR"
echo "╠══════════════════════════════════════════════════╣"
printf "║  00_fastqc_raw      : %-26s║\n" "$FASTQC_RAW_DIR"
printf "║  01_trimmed         : %-26s║\n" "$TRIM_DIR"
printf "║  02_mapped          : %-26s║\n" "$MAP_DIR"
printf "║  03_mapped_fastq    : %-26s║\n" "$MAPPED_FQ_DIR"
printf "║  04_shovill_denovo  : %-26s║\n" "$DENOVO_DIR"
printf "║  05_ragtag_scaffold : %-26s║\n" "$RAGTAG_DIR"
printf "║  06_final_consensus : %-26s║\n" "$FINAL_DIR"
printf "║  07_coverage_map    : %-26s║\n" "$COVERAGE_DIR"
printf "║  08_fastqc_final    : %-26s║\n" "$FASTQC_FINAL_DIR"
printf "║  logs               : %-26s║\n" "$LOG_DIR"
echo "╚══════════════════════════════════════════════════╝"
echo ""

START_TIME=$(date +%s)

# ---------------------------------------------------------------------------
# STEP 0: FastQC on raw reads — quality before any processing
# ---------------------------------------------------------------------------
log_step "STEP 0/9 · FastQC on raw reads"

FASTQC_RAW_LOG="${LOG_DIR}/${PREFIX}_fastqc_raw.log"

if ls "${FASTQC_RAW_DIR}"/*_fastqc.html 2>/dev/null | grep -q .; then
    log_skip "FastQC raw reads"
else
    log_info "Input R1 : $READ1"
    log_info "Input R2 : $READ2"
    log_info "Output   : $FASTQC_RAW_DIR"

    fastqc \
        "$READ1" "$READ2" \
        --outdir "$FASTQC_RAW_DIR" \
        --threads "$THREADS" \
        --extract \
        2> "$FASTQC_RAW_LOG"

    log_ok "FastQC raw reads complete"
fi
log_file "${FASTQC_RAW_DIR}"

# ---------------------------------------------------------------------------
# STEP 1: FASTP trimming
# ---------------------------------------------------------------------------
log_step "STEP 1/9 · FASTP trimming"

if [ -f "$TRIM_R1" ] && [ -f "$TRIM_R2" ]; then
    log_skip "Trimmed reads"
else
    log_info "Input R1 : $READ1"
    log_info "Input R2 : $READ2"
    log_info "Log      : ${LOG_DIR}/${PREFIX}_fastp.log"

    fastp \
        --in1 "$READ1" --in2 "$READ2" \
        --out1 "$TRIM_R1" --out2 "$TRIM_R2" \
        --html "${TRIM_DIR}/${PREFIX}_fastp.html" \
        --json "${TRIM_DIR}/${PREFIX}_fastp.json" \
        --thread 4 \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --qualified_quality_phred 20 \
        --length_required 50 \
        2> "${LOG_DIR}/${PREFIX}_fastp.log"

    log_ok "Trimming complete"
fi
log_file "$TRIM_R1"
log_file "$TRIM_R2"

# ---------------------------------------------------------------------------
# STEP 2: bwa-mem2 mapping to reference (-F 12: proper pairs only)
# ---------------------------------------------------------------------------
log_step "STEP 2/9 · Map to reference genome"

if [ -f "$SORTED_BAM" ]; then
    log_skip "Reference mapping"
else
    log_info "Indexing reference..."
    bwa-mem2 index "$REFERENCE" \
        2> "${LOG_DIR}/${PREFIX}_bwa_index.log"

    log_info "Mapping trimmed reads..."
    log_info "Filter    : -F 12 (proper pairs only), -q 10 (MAPQ)"
    log_info "Log       : ${LOG_DIR}/${PREFIX}_bwamem2.log"

    bwa-mem2 mem \
        -t "$THREADS" \
        "$REFERENCE" \
        "$TRIM_R1" "$TRIM_R2" \
        2> "${LOG_DIR}/${PREFIX}_bwamem2.log" \
    | samtools view -bS -F 12 -q "$MIN_MAP_QUAL" - \
    | samtools sort -@ 4 -o "$SORTED_BAM" -

    samtools index "$SORTED_BAM"
    samtools flagstat "$SORTED_BAM" > "${MAP_DIR}/${PREFIX}_ref_flagstat.txt"

    log_ok "Reference mapping complete"
fi
log_file "$SORTED_BAM"
log_file "${MAP_DIR}/${PREFIX}_ref_flagstat.txt"

# ---------------------------------------------------------------------------
# STEP 3: Extract mapped reads as FASTQ pairs (avoids name-sort)
# ---------------------------------------------------------------------------
log_step "STEP 3/9 · Extract mapped reads"

if [ -f "$MAPPED_R1" ] && [ -f "$MAPPED_R2" ]; then
    log_skip "Mapped read extraction"
else
    log_info "Extracting paired reads from BAM..."
    log_info "Method   : samtools collate + samtools fastq"

    samtools collate -u -O "$SORTED_BAM" \
    | samtools fastq \
        -@ 4 \
        -1 "$MAPPED_R1" \
        -2 "$MAPPED_R2" \
        -0 /dev/null -s /dev/null -n - \
        2> "${LOG_DIR}/${PREFIX}_extract.log"

    log_ok "Mapped read extraction complete"
fi
log_file "$MAPPED_R1"
log_file "$MAPPED_R2"

# ---------------------------------------------------------------------------
# STEP 4: Shovill de novo assembly
# ---------------------------------------------------------------------------
log_step "STEP 4/9 · De novo assembly (Shovill/SPAdes)"

if [ -f "$DENOVO_CONTIGS" ]; then
    log_skip "De novo assembly"
else
    log_info "Assembling mapped reads..."
    log_info "Assembler: SPAdes (via Shovill)"
    log_info "K-mers   : 71, 91, 111"
    log_info "Target depth: ${SHOVILL_DEPTH}×"
    log_info "Watch    : tail -f ${LOG_DIR}/${PREFIX}_shovill.log"

    shovill \
        --R1 "$MAPPED_R1" --R2 "$MAPPED_R2" \
        --outdir shovill_out \
        --assembler spades \
        --gsize "$GENOME_SIZE" \
        --depth "$SHOVILL_DEPTH" \
        --kmers 71,91,111 \
        --cpus "$THREADS" \
        --ram 8 \
        --minlen 200 \
        2> "${LOG_DIR}/${PREFIX}_shovill.log"

    cp shovill_out/contigs.fa "$DENOVO_CONTIGS"

    log_ok "De novo assembly complete"
fi

CONTIG_COUNT=$(grep -c ">" "$DENOVO_CONTIGS" || echo "0")
log_info "Contigs: $CONTIG_COUNT"
log_file "$DENOVO_CONTIGS"

# ---------------------------------------------------------------------------
# STEP 5: RagTag scaffolding against reference
# ---------------------------------------------------------------------------
log_step "STEP 5/9 · Scaffold contigs"

if [ -f "$SCAFFOLD" ]; then
    log_skip "RagTag scaffolding"
else
    log_info "Running RagTag..."
    log_info "Watch    : tail -f ${LOG_DIR}/${PREFIX}_ragtag.log"

    ragtag.py scaffold \
        "$REFERENCE" \
        "$DENOVO_CONTIGS" \
        -o "$RAGTAG_DIR" \
        -t "$THREADS" \
        -u \
        2> "${LOG_DIR}/${PREFIX}_ragtag.log"

    log_ok "RagTag scaffolding complete"
fi
log_file "$SCAFFOLD"

# ---------------------------------------------------------------------------
# STEP 6: Remap to scaffold + Pilon polishing
# ---------------------------------------------------------------------------
log_step "STEP 6/9 · Pilon consensus polishing"

if [ -f "$FINAL_CONSENSUS" ]; then
    log_skip "Pilon final consensus"
else
    log_info "Indexing scaffold..."
    bwa-mem2 index "$SCAFFOLD" \
        2> "${LOG_DIR}/${PREFIX}_scaffold_index.log"

    log_info "Mapping trimmed reads to scaffold..."
    log_info "Filter    : -F 12 (proper pairs only)"

    bwa-mem2 mem \
        -t "$THREADS" \
        "$SCAFFOLD" \
        "$TRIM_R1" "$TRIM_R2" \
        2> "${LOG_DIR}/${PREFIX}_scaffold_bwa.log" \
    | samtools view -bS -F 12 -q "$MIN_MAP_QUAL" - \
    | samtools sort -@ 4 -o "$SCAFFOLD_BAM" -

    samtools index "$SCAFFOLD_BAM"
    samtools flagstat "$SCAFFOLD_BAM" > "${RAGTAG_DIR}/${PREFIX}_scaffold_flagstat.txt"

    log_info "Running Pilon..."
    log_info "Watch live: tail -f ${LOG_DIR}/${PREFIX}_pilon.log"

    pilon \
        --genome "$SCAFFOLD" \
        --frags "$SCAFFOLD_BAM" \
        --output "${PREFIX}_final_consensus" \
        --outdir "$FINAL_DIR" \
        --changes \
        --vcf \
        --fix all \
        --threads "$THREADS" \
        2> "${LOG_DIR}/${PREFIX}_pilon.log"

    log_ok "Pilon complete"
fi

N_COUNT=$(grep -v ">" "$FINAL_CONSENSUS" | tr -cd 'N' | wc -c || true)
log_info "Remaining Ns: ${N_COUNT}"
log_file "$FINAL_CONSENSUS"
log_file "${FINAL_DIR}/${PREFIX}_final_consensus.vcf"

# ---------------------------------------------------------------------------
# STEP 7: Coverage map — fresh mapping to final consensus
# ---------------------------------------------------------------------------
log_step "STEP 7/9 · Coverage map"

if [ -f "${COVERAGE_PLOT}.jpg" ]; then
    log_skip "Coverage map"
else
    log_info "Indexing final consensus..."
    bwa-mem2 index "$FINAL_CONSENSUS" \
        2> "${LOG_DIR}/${PREFIX}_final_index.log"

    log_info "Mapping trimmed reads to final consensus..."
    log_info "Filter    : -F 12 (proper pairs only)"
    log_info "Log       : ${LOG_DIR}/${PREFIX}_final_bwa.log"

    bwa-mem2 mem \
        -t "$THREADS" \
        "$FINAL_CONSENSUS" \
        "$TRIM_R1" "$TRIM_R2" \
        2> "${LOG_DIR}/${PREFIX}_final_bwa.log" \
    | samtools view -bS -F 12 -q "$MIN_MAP_QUAL" - \
    | samtools sort -@ "$THREADS" -o "$FINAL_BAM" -

    samtools index "$FINAL_BAM"
    samtools flagstat "$FINAL_BAM" > "${COVERAGE_DIR}/${PREFIX}_final_flagstat.txt"
    samtools stats "$FINAL_BAM" > "${COVERAGE_DIR}/${PREFIX}_final_stats.txt"

    log_info "Calculating depth..."
    samtools depth "$FINAL_BAM" > "$DEPTH_TXT"

    log_info "Generating coverage plot..."
    python3 /usr/local/bin/Draw_SequencingDepth.py \
        "$DEPTH_TXT" \
        "$COVERAGE_PLOT"

    log_ok "Coverage map complete"
fi

FINAL_MAPPED=$(samtools view -c "$FINAL_BAM" 2>/dev/null || echo "0")
log_file "$FINAL_BAM"
log_file "$DEPTH_TXT"
log_file "${COVERAGE_PLOT}.jpg"

# ---------------------------------------------------------------------------
# STEP 8: FastQC on final mapped reads (quality of reads on consensus)
# ---------------------------------------------------------------------------
log_step "STEP 8/9 · FastQC on trimmed reads"

if ls "${FASTQC_FINAL_DIR}"/*_fastqc.html 2>/dev/null | grep -q .; then
    log_skip "FastQC final trimmed reads"
else
    log_info "Running FastQC on trimmed reads (post-trimming quality)..."
    log_info "Output : $FASTQC_FINAL_DIR"

    fastqc \
        "$TRIM_R1" "$TRIM_R2" \
        --outdir "$FASTQC_FINAL_DIR" \
        --threads "$THREADS" \
        --extract \
        2> "${LOG_DIR}/${PREFIX}_fastqc_final.log"

    log_ok "FastQC final reads complete"
fi

FINAL_FASTQC_REPORT=$(ls "${FASTQC_FINAL_DIR}"/*_fastqc.html 2>/dev/null | head -1 || echo "")
[ -n "$FINAL_FASTQC_REPORT" ] && log_file "$FINAL_FASTQC_REPORT"

# ---------------------------------------------------------------------------
# Generate comprehensive assembly report
# ---------------------------------------------------------------------------
CONSENSUS_LEN=$(grep -v ">" "$FINAL_CONSENSUS" | tr -d '\n' | wc -c || echo "0")
CONSENSUS_N=$(grep -v ">" "$FINAL_CONSENSUS" | tr -cd 'N' | wc -c || echo "0")
CONSENSUS_GC=$(grep -v ">" "$FINAL_CONSENSUS" | tr -d '\n' | \
    awk '{g=gsub(/[GCgc]/,""); t=length($0); printf "%.2f", (g/t)*100}' || echo "0")
AVG_DEPTH=$(awk '{sum+=$3; n++} END {printf "%.1f", sum/n}' "$DEPTH_TXT" 2>/dev/null || echo "0")
MIN_DEPTH=$(awk 'NR==1{m=$3} $3<m{m=$3} END{print m}' "$DEPTH_TXT" 2>/dev/null || echo "0")
MAX_DEPTH=$(awk 'BEGIN{m=0} $3>m{m=$3} END{print m}' "$DEPTH_TXT" 2>/dev/null || echo "0")
COVERED_PCT=$(awk -v len="$CONSENSUS_LEN" '$3>0{n++} END{printf "%.2f", (n/len)*100}' "$DEPTH_TXT" 2>/dev/null || echo "0")
LOW_DEPTH_BASES=$(awk '$3<10{n++} END{print n+0}' "$DEPTH_TXT" 2>/dev/null || echo "0")

SNPS_CORRECTED=$(grep -c "^" "${FINAL_DIR}/${PREFIX}_final_consensus.changes" 2>/dev/null || echo "0")
GAPS_FILLED=$(grep -c "ClosedGap" "${LOG_DIR}/${PREFIX}_pilon.log" 2>/dev/null || echo "0")

# Calculate elapsed time for report
REPORT_END=$(date +%s)
REPORT_ELAPSED=$(( REPORT_END - START_TIME ))
REPORT_MIN=$(( REPORT_ELAPSED / 60 ))
REPORT_SEC=$(( REPORT_ELAPSED % 60 ))

log_step "STEP 9/9 · Generate HTML report"
log_info "Generating assembly report..."

python3 /usr/local/bin/generate_report.py \
    --prefix          "$PREFIX" \
    --outdir          "$OUTDIR" \
    --depth           "$DEPTH_TXT" \
    --consensus_len   "$CONSENSUS_LEN" \
    --consensus_gc    "$CONSENSUS_GC" \
    --consensus_n     "$CONSENSUS_N" \
    --avg_depth       "$AVG_DEPTH" \
    --min_depth       "$MIN_DEPTH" \
    --max_depth       "$MAX_DEPTH" \
    --covered_pct     "$COVERED_PCT" \
    --low_depth_bases "$LOW_DEPTH_BASES" \
    --n_count         "$N_COUNT" \
    --contig_count    "$CONTIG_COUNT" \
    --snps_corrected  "$SNPS_CORRECTED" \
    --gaps_filled     "$GAPS_FILLED" \
    --final_mapped    "$FINAL_MAPPED" \
    --elapsed_min     "$REPORT_MIN" \
    --elapsed_sec     "$REPORT_SEC" \
    --threads         "$THREADS" \
    --denovo_reads    "0" \
    --fastqc_raw_dir  "$FASTQC_RAW_DIR" \
    --fastqc_final_dir "$FASTQC_FINAL_DIR"

REPORT_FILE="${OUTDIR}/${PREFIX}_assembly_report.html"
log_ok "HTML assembly report written"
log_file "$REPORT_FILE"

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
END_TIME=$(date +%s)
ELAPSED=$(( END_TIME - START_TIME ))
ELAPSED_MIN=$(( ELAPSED / 60 ))
ELAPSED_SEC=$(( ELAPSED % 60 ))

echo ""
echo "╔══════════════════════════════════════════════════╗"
echo "║              Pipeline Complete                   ║"
echo "╠══════════════════════════════════════════════════╣"
printf "║  Sample        : %-31s║\n" "$PREFIX"
printf "║  Runtime       : %-31s║\n" "${ELAPSED_MIN}m ${ELAPSED_SEC}s"
echo "╠══════════════════════════════════════════════════╣"
echo "║  CONSENSUS                                       ║"
echo "╠══════════════════════════════════════════════════╣"
printf "║  Length        : %-31s║\n" "${CONSENSUS_LEN} bp"
printf "║  GC content    : %-31s║\n" "${CONSENSUS_GC} %"
printf "║  Remaining Ns  : %-31s║\n" "$N_COUNT"
printf "║  Contigs       : %-31s║\n" "$CONTIG_COUNT"
echo "╠══════════════════════════════════════════════════╣"
echo "║  COVERAGE                                        ║"
echo "╠══════════════════════════════════════════════════╣"
printf "║  Average depth : %-31s║\n" "${AVG_DEPTH} x"
printf "║  Min depth     : %-31s║\n" "${MIN_DEPTH} x"
printf "║  Max depth     : %-31s║\n" "${MAX_DEPTH} x"
printf "║  Breadth       : %-31s║\n" "${COVERED_PCT} %"
printf "║  Bases < 10x   : %-31s║\n" "${LOW_DEPTH_BASES} bp"
echo "╠══════════════════════════════════════════════════╣"
echo "║  OUTPUT FILES                                    ║"
echo "╠══════════════════════════════════════════════════╣"
printf "║  Assembly report: %-30s║\n" "${OUTDIR}/${PREFIX}_assembly_report.html"
printf "║  Coverage plot  : %-30s║\n" "${COVERAGE_PLOT}.jpg"
printf "║  Variant calls  : %-30s║\n" "${FINAL_DIR}/${PREFIX}_final_consensus.vcf"
echo "╠══════════════════════════════════════════════════╣"
printf "║  All logs in    : %-30s║\n" "$LOG_DIR"
echo "╚══════════════════════════════════════════════════╝"
echo ""
