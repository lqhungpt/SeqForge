#!/usr/bin/env nextflow
// =============================================================================
// SeqForge — Short-read Genome Assembly Pipeline — Nextflow
// Steps: FastQC → FASTP → bwa-mem2 → Shovill/SPAdes → RagTag → Pilon
//        → FastQC → Coverage → HTML Report
// Author: hungluong
// =============================================================================

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------
params.reads        = null          // glob: "/data/*_{R1,R2}.fastq.gz"
params.reference    = null          // FASTA reference (required)
params.genome_size  = null          // genome size in bp (required)
params.outdir       = "results"
params.depth        = 150
params.threads      = 6
params.help         = false

def helpMsg() {
    log.info """
    ╔════════════════════════════════════════════════════════╗
    ║   SeqForge — Short-read Genome Assembly (Nextflow)    ║
    ╠════════════════════════════════════════════════════════╣
    ║  Usage:                                                ║
    ║    nextflow run main.nf \\                             ║
    ║      --reads '/data/*_{R1,R2}.fastq.gz' \\             ║
    ║      --reference /data/reference.fasta \\              ║
    ║      --genome_size 3200000000                          ║
    ║                                                        ║
    ║  Required:                                             ║
    ║    --reads           FASTQ glob pattern                ║
    ║    --reference       Reference genome FASTA            ║
    ║    --genome_size     Estimated genome size (bp)        ║
    ║                                                        ║
    ║  Optional:                                             ║
    ║    --outdir          Output directory (default: results)║
    ║    --depth           Assembly depth x (default: 150)   ║
    ║    --threads         CPU threads      (default: 6)     ║
    ║                                                        ║
    ║  Genome size examples:                                 ║
    ║    SARS-CoV-2: 29903                                   ║
    ║    E. coli: 4600000                                    ║
    ║    Arabidopsis: 135000000                              ║
    ║    Human: 3200000000                                   ║
    ║    Wheat: 17000000000                                  ║
    ╚════════════════════════════════════════════════════════╝
    """.stripIndent()
}

if (params.help) { helpMsg(); exit 0 }
if (!params.reads) {
    log.error "ERROR: --reads is required. Use --help for usage."
    exit 1
}
if (!params.reference) {
    log.error "ERROR: --reference is required. Use --help for usage."
    exit 1
}
if (!params.genome_size) {
    log.error "ERROR: --genome_size is required. Use --help for usage."
    exit 1
}

// ---------------------------------------------------------------------------
// Input channels
// ---------------------------------------------------------------------------
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .set { reads_ch }

Channel
    .fromPath(params.reference, checkIfExists: true)
    .set { ref_ch }

// ---------------------------------------------------------------------------
// PROCESS 0: FastQC on raw reads
// ---------------------------------------------------------------------------
process FASTQC_RAW {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/00_fastqc_raw", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.{html,zip}"), emit: reports
    path "*.log",                                       emit: log

    script:
    """
    fastqc ${reads[0]} ${reads[1]} \\
        --outdir . \\
        --threads ${params.threads} \\
        --extract \\
        2> ${sample_id}_fastqc_raw.log
    """
}

// ---------------------------------------------------------------------------
// PROCESS 1: FASTP trimming
// ---------------------------------------------------------------------------
process FASTP {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/01_trimmed", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fastq.gz"),
                          path("${sample_id}_R2_trimmed.fastq.gz"), emit: trimmed
    path "${sample_id}_fastp.{html,json}",                           emit: qc
    path "*.log",                                                    emit: log

    script:
    """
    fastp \\
        --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${sample_id}_R1_trimmed.fastq.gz \\
        --out2 ${sample_id}_R2_trimmed.fastq.gz \\
        --html ${sample_id}_fastp.html \\
        --json ${sample_id}_fastp.json \\
        --thread 4 \\
        --detect_adapter_for_pe \\
        --trim_poly_g \\
        --qualified_quality_phred 20 \\
        --length_required 50 \\
        2> ${sample_id}_fastp.log
    """
}

// ---------------------------------------------------------------------------
// PROCESS 2: bwa-mem2 mapping to reference
// ---------------------------------------------------------------------------
process BWA_MAP {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/02_mapped", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"),
                          path("${sample_id}_sorted.bam.bai"), emit: bam
    path "${sample_id}_ref_flagstat.txt",                      emit: flagstat
    path "*.log",                                              emit: log

    script:
    """
    bwa-mem2 index ${reference} 2> ${sample_id}_bwa_index.log

    bwa-mem2 mem \\
        -t ${params.threads} \\
        ${reference} \\
        ${r1} ${r2} \\
        2> ${sample_id}_bwamem2.log \\
    | samtools view -bS -F 12 -q 10 - \\
    | samtools sort -@ 4 -o ${sample_id}_sorted.bam -

    samtools index ${sample_id}_sorted.bam
    samtools flagstat ${sample_id}_sorted.bam > ${sample_id}_ref_flagstat.txt
    """
}

// ---------------------------------------------------------------------------
// PROCESS 3: Extract mapped reads as FASTQ
// ---------------------------------------------------------------------------
process EXTRACT_MAPPED {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/03_mapped_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}_mapped_R1.fastq.gz"),
                          path("${sample_id}_mapped_R2.fastq.gz"), emit: fastq
    path "*.log",                                                   emit: log

    script:
    """
    samtools collate -u -O ${bam} \\
    | samtools fastq \\
        -@ 4 \\
        -1 ${sample_id}_mapped_R1.fastq.gz \\
        -2 ${sample_id}_mapped_R2.fastq.gz \\
        -0 /dev/null -s /dev/null -n - \\
        2> ${sample_id}_extract.log
    """
}

// ---------------------------------------------------------------------------
// PROCESS 4: Shovill de novo assembly
// ---------------------------------------------------------------------------
process SHOVILL {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/04_shovill_denovo", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta"), emit: contigs
    path "*.log",                                              emit: log

    script:
    """
    shovill \\
        --R1 ${r1} --R2 ${r2} \\
        --outdir shovill_out \\
        --assembler spades \\
        --gsize ${params.genome_size} \\
        --depth ${params.depth} \\
        --kmers 71,91,111 \\
        --cpus ${params.threads} \\
        --ram 8 \\
        --minlen 200 \\
        2> ${sample_id}_shovill.log

    cp shovill_out/contigs.fa ${sample_id}_contigs.fasta
    """
}

// ---------------------------------------------------------------------------
// PROCESS 5: RagTag scaffolding
// ---------------------------------------------------------------------------
process RAGTAG {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/05_ragtag_scaffold", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)
    path reference

    output:
    tuple val(sample_id), path("ragtag.scaffold.fasta"), emit: scaffold
    path "ragtag.scaffold.stats",                        emit: stats
    path "*.log",                                        emit: log

    script:
    """
    ragtag.py scaffold \\
        ${reference} \\
        ${contigs} \\
        -o . \\
        -t ${params.threads} \\
        -u \\
        2> ${sample_id}_ragtag.log
    """
}

// ---------------------------------------------------------------------------
// PROCESS 6: Remap to scaffold + Pilon polishing
// ---------------------------------------------------------------------------
process PILON {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/06_final_consensus", mode: 'copy'

    input:
    tuple val(sample_id), path(scaffold)
    tuple val(sample_id2), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_final_consensus.fasta"), emit: consensus
    path "${sample_id}_final_consensus.vcf",                          emit: vcf
    path "${sample_id}_final_consensus.changes",                      emit: changes
    path "${sample_id}_scaffold_flagstat.txt",                        emit: flagstat
    path "*.log",                                                     emit: log

    script:
    """
    bwa-mem2 index ${scaffold} 2> ${sample_id}_scaffold_index.log

    bwa-mem2 mem \\
        -t ${params.threads} \\
        ${scaffold} \\
        ${r1} ${r2} \\
        2> ${sample_id}_scaffold_bwa.log \\
    | samtools view -bS -F 12 -q 10 - \\
    | samtools sort -@ 4 -o ${sample_id}_scaffold_sorted.bam -

    samtools index ${sample_id}_scaffold_sorted.bam
    samtools flagstat ${sample_id}_scaffold_sorted.bam > ${sample_id}_scaffold_flagstat.txt

    pilon \\
        --genome ${scaffold} \\
        --frags ${sample_id}_scaffold_sorted.bam \\
        --output ${sample_id}_final_consensus \\
        --outdir . \\
        --changes \\
        --vcf \\
        --fix all \\
        --threads ${params.threads} \\
        2> ${sample_id}_pilon.log
    """
}

// ---------------------------------------------------------------------------
// PROCESS 7: Coverage map
// ---------------------------------------------------------------------------
process COVERAGE {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}/${sample_id}/07_coverage_map", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus)
    tuple val(sample_id2), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_final_sorted.bam"),
                          path("${sample_id}_final_sorted.bam.bai"), emit: bam
    path "${sample_id}_depth.txt",                                   emit: depth
    path "${sample_id}_coverage.jpg",                                emit: plot
    path "${sample_id}_final_flagstat.txt",                          emit: flagstat
    path "${sample_id}_final_stats.txt",                             emit: stats
    path "*.log",                                                    emit: log

    script:
    """
    bwa-mem2 index ${consensus} 2> ${sample_id}_final_index.log

    bwa-mem2 mem \\
        -t ${params.threads} \\
        ${consensus} \\
        ${r1} ${r2} \\
        2> ${sample_id}_final_bwa.log \\
    | samtools view -bS -F 12 -q 10 - \\
    | samtools sort -@ ${params.threads} -o ${sample_id}_final_sorted.bam -

    samtools index ${sample_id}_final_sorted.bam
    samtools flagstat ${sample_id}_final_sorted.bam > ${sample_id}_final_flagstat.txt
    samtools stats   ${sample_id}_final_sorted.bam > ${sample_id}_final_stats.txt
    samtools depth   ${sample_id}_final_sorted.bam > ${sample_id}_depth.txt

    python3 /usr/local/bin/Draw_SequencingDepth.py \\
        ${sample_id}_depth.txt \\
        ${sample_id}_coverage
    """
}

// ---------------------------------------------------------------------------
// PROCESS 8: FastQC on trimmed reads
// ---------------------------------------------------------------------------
process FASTQC_TRIMMED {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/08_fastqc_final", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("*_fastqc.{html,zip}"), emit: reports
    path "*.log",                                       emit: log

    script:
    """
    fastqc ${r1} ${r2} \\
        --outdir . \\
        --threads ${params.threads} \\
        --extract \\
        2> ${sample_id}_fastqc_trimmed.log
    """
}

// ---------------------------------------------------------------------------
// PROCESS 9: Generate HTML report - SIMPLIFIED VERSION
// ---------------------------------------------------------------------------
process REPORT {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus)
    path depth

    output:
    path "${sample_id}_assembly_report.html", emit: report

    script:
    """
    CONSENSUS_LEN=\$(grep -v ">" ${consensus} | tr -d '\\n' | wc -c)
    CONSENSUS_GC=\$(grep -v ">" ${consensus} | tr -d '\\n' | awk '{g=gsub(/[GCgc]/,""); t=length(\$0); printf "%.1f", (g/t)*100}')
    CONSENSUS_N=\$(grep -v ">" ${consensus} | tr -cd 'N' | wc -c)
    AVG_DEPTH=\$(awk '{sum+=\$3; n++} END {printf "%.1f", sum/n}' ${depth})
    MIN_DEPTH=\$(awk 'NR==1{m=\$3} \$3<m{m=\$3} END{print m}' ${depth})
    MAX_DEPTH=\$(awk 'BEGIN{m=0} \$3>m{m=\$3} END{print m}' ${depth})
    COVERED_PCT=\$(awk -v len="\$CONSENSUS_LEN" '\$3>0{n++} END{printf "%.1f", (n/len)*100}' ${depth})

    cat > ${sample_id}_assembly_report.html << 'EOFHTML'
<!DOCTYPE html>
<html>
<head><title>Genome Assembly Report</title>
<style>
body { font-family: Arial; margin: 20px; background: #f0f4f8; color: #2d3748; }
h1 { color: #2b6cb0; border-bottom: 2px solid #2b6cb0; padding-bottom: 10px; }
.card { background: white; padding: 20px; margin: 15px 0; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,.1); }
h2 { color: #718096; font-size: 16px; margin-top: 0; text-transform: uppercase; }
table { width: 100%; border-collapse: collapse; }
th, td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }
th { background: #ebf4ff; color: #2b6cb0; font-weight: 600; }
tr:hover { background: #f7fafc; }
.stat-row { font-weight: bold; }
footer { margin-top: 30px; padding: 20px; text-align: center; color: #a0aec0; font-size: 12px; }
</style></head>
<body>
<h1>🧬 Genome Assembly Report - ${sample_id}</h1>
<div class="card">
<h2>Assembly Metrics</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr class="stat-row"><td>Consensus Length</td><td>\$CONSENSUS_LEN bp</td></tr>
<tr><td>GC Content</td><td>\$CONSENSUS_GC %</td></tr>
<tr><td>Remaining Ns</td><td>\$CONSENSUS_N</td></tr>
<tr class="stat-row"><td>Average Depth</td><td>\$AVG_DEPTH ×</td></tr>
<tr><td>Min Depth</td><td>\$MIN_DEPTH ×</td></tr>
<tr><td>Max Depth</td><td>\$MAX_DEPTH ×</td></tr>
<tr class="stat-row"><td>Breadth of Coverage</td><td>\$COVERED_PCT %</td></tr>
</table>
</div>
<footer>
<p>SeqForge v2.0 - Portable Genome Assembly Pipeline</p>
<p>For support: https://github.com/lqhungpt/SeqForge</p>
</footer>
</body>
</html>
EOFHTML
    """
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    FASTQC_RAW(reads_ch)
    FASTP(reads_ch)
    BWA_MAP(FASTP.out.trimmed, ref_ch.first())
    EXTRACT_MAPPED(BWA_MAP.out.bam)
    SHOVILL(EXTRACT_MAPPED.out.fastq)
    RAGTAG(SHOVILL.out.contigs, ref_ch.first())
    PILON(RAGTAG.out.scaffold, FASTP.out.trimmed)
    COVERAGE(PILON.out.consensus, FASTP.out.trimmed)
    FASTQC_TRIMMED(FASTP.out.trimmed)
    REPORT(PILON.out.consensus, COVERAGE.out.depth)
}
