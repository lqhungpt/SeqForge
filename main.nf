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
    # Index reference
    bwa-mem2 index ${reference} 2> ${sample_id}_bwa_index.log

    # Map and sort
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
    # Index scaffold
    bwa-mem2 index ${scaffold} 2> ${sample_id}_scaffold_index.log

    # Map trimmed reads to scaffold
    bwa-mem2 mem \\
        -t ${params.threads} \\
        ${scaffold} \\
        ${r1} ${r2} \\
        2> ${sample_id}_scaffold_bwa.log \\
    | samtools view -bS -F 12 -q 10 - \\
    | samtools sort -@ 4 -o ${sample_id}_scaffold_sorted.bam -

    samtools index ${sample_id}_scaffold_sorted.bam
    samtools flagstat ${sample_id}_scaffold_sorted.bam > ${sample_id}_scaffold_flagstat.txt

    # Pilon polishing
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
// PROCESS 9: Generate HTML report
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

    shell:
    '''
    # Calculate assembly statistics
    CONSENSUS_LEN=$(grep -v ">" !{consensus} | tr -d '\n' | wc -c)
    CONSENSUS_GC=$(grep -v ">" !{consensus} | tr -d '\n' | awk '{g=gsub(/[GCgc]/,""); t=length($0); printf "%.1f", (g/t)*100}')
    CONSENSUS_N=$(grep -v ">" !{consensus} | tr -cd 'N' | wc -c)
    
    # Calculate coverage statistics
    AVG_DEPTH=$(awk '{sum+=$3; n++} END {printf "%.1f", sum/n}' !{depth})
    MIN_DEPTH=$(awk 'NR==1{m=$3} $3<m{m=$3} END{print m}' !{depth})
    MAX_DEPTH=$(awk 'BEGIN{m=0} $3>m{m=$3} END{print m}' !{depth})
    COVERED_PCT=$(awk -v len="$CONSENSUS_LEN" '$3>0{n++} END{printf "%.1f", (n/len)*100}' !{depth})

    # Generate HTML report
    python3 << 'EOFPYTHON'
import datetime

sample_id = "!{sample_id}"
consensus_len = int("$CONSENSUS_LEN") if "$CONSENSUS_LEN" else 0
consensus_gc = float("$CONSENSUS_GC") if "$CONSENSUS_GC" else 0
consensus_n = int("$CONSENSUS_N") if "$CONSENSUS_N" else 0
avg_depth = float("$AVG_DEPTH") if "$AVG_DEPTH" else 0
min_depth = int("$MIN_DEPTH") if "$MIN_DEPTH" else 0
max_depth = int("$MAX_DEPTH") if "$MAX_DEPTH" else 0
covered_pct = float("$COVERED_PCT") if "$COVERED_PCT" else 0

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Genome Assembly Report — {sample_id}</title>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         background: #f0f4f8; color: #2d3748; }}
  header {{ background: linear-gradient(135deg, #1a365d 0%, #2b6cb0 100%);
            color: white; padding: 28px 40px; }}
  header h1 {{ font-size: 24px; font-weight: 700; }}
  header p {{ font-size: 13px; opacity: .8; margin-top: 4px; }}
  .container {{ max-width: 1200px; margin: 0 auto; padding: 28px 20px; }}
  .grid-5 {{ display: grid; grid-template-columns: repeat(5,1fr); gap: 16px; margin-bottom: 20px; }}
  .grid-3 {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 20px; margin-bottom: 20px; }}
  .card {{ background: white; border-radius: 12px; padding: 20px; box-shadow: 0 1px 3px rgba(0,0,0,.08); }}
  .card h2 {{ font-size: 14px; font-weight: 600; color: #718096; text-transform: uppercase; margin-bottom: 14px; }}
  .stat-card {{ background: white; border-radius: 12px; padding: 18px; box-shadow: 0 1px 3px rgba(0,0,0,.08); text-align: center; }}
  .stat-card .val {{ font-size: 28px; font-weight: 700; color: #2b6cb0; }}
  .stat-card .lbl {{ font-size: 12px; color: #718096; margin-top: 4px; }}
  .stat-card.green .val {{ color: #276749; }}
  .stat-card.red .val {{ color: #c53030; }}
  table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
  th {{ background: #ebf4ff; padding: 9px 12px; text-align: left; font-weight: 600; color: #2b6cb0; border-bottom: 2px solid #bee3f8; }}
  td {{ padding: 8px 12px; border-bottom: 1px solid #e2e8f0; }}
  .section-title {{ font-size: 18px; font-weight: 700; color: #1a365d; margin: 28px 0 12px; padding-left: 12px; border-left: 4px solid #2b6cb0; }}
  footer {{ text-align: center; padding: 24px; font-size: 12px; color: #a0aec0; }}
</style>
</head>
<body>

<header>
  <h1>🧬 Genome Assembly Report — {sample_id}</h1>
  <p>Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | Pipeline: FASTP → bwa-mem2 → Shovill/SPAdes → RagTag → Pilon</p>
</header>

<div class="container">

  <p class="section-title">Assembly Summary</p>
  <div class="grid-5">
    <div class="stat-card green">
      <div class="val">{consensus_len:,}</div>
      <div class="lbl">Consensus Length (bp)</div>
    </div>
    <div class="stat-card green">
      <div class="val">{avg_depth:.1f}×</div>
      <div class="lbl">Average Coverage Depth</div>
    </div>
    <div class="stat-card green">
      <div class="val">{consensus_gc:.1f}%</div>
      <div class="lbl">GC Content</div>
    </div>
    <div class="stat-card red">
      <div class="val">{consensus_n}</div>
      <div class="lbl">Remaining Ns</div>
    </div>
    <div class="stat-card green">
      <div class="val">{covered_pct:.1f}%</div>
      <div class="lbl">Breadth of Coverage</div>
    </div>
  </div>

  <p class="section-title">Coverage Statistics</p>
  <div class="grid-3">
    <div class="card">
      <h2>Min Depth</h2>
      <div style="font-size: 24px; font-weight: 700; color: #c53030;">{min_depth}×</div>
    </div>
    <div class="card">
      <h2>Max Depth</h2>
      <div style="font-size: 24px; font-weight: 700; color: #276749;">{max_depth}×</div>
    </div>
    <div class="card">
      <h2>Avg Depth</h2>
      <div style="font-size: 24px; font-weight: 700; color: #2b6cb0;">{avg_depth:.1f}×</div>
    </div>
  </div>

  <p class="section-title">Assembly Details</p>
  <div class="card">
    <table>
      <tr>
        <th>Metric</th>
        <th>Value</th>
      </tr>
      <tr>
        <td>Consensus Length</td>
        <td>{consensus_len:,} bp</td>
      </tr>
      <tr>
        <td>GC Content</td>
        <td>{consensus_gc:.1f}%</td>
      </tr>
      <tr>
        <td>Remaining Ns</td>
        <td>{consensus_n}</td>
      </tr>
      <tr>
        <td>Min Depth</td>
        <td>{min_depth}×</td>
      </tr>
      <tr>
        <td>Max Depth</td>
        <td>{max_depth}×</td>
      </tr>
      <tr>
        <td>Average Depth</td>
        <td>{avg_depth:.1f}×</td>
      </tr>
      <tr>
        <td>Breadth of Coverage</td>
        <td>{covered_pct:.1f}%</td>
      </tr>
    </table>
  </div>

</div>

<footer>
  <p>SeqForge v2.0 — Portable Genome Assembly Pipeline</p>
  <p>For support, visit: https://github.com/lqhungpt/SeqForge</p>
</footer>

</body>
</html>
"""

with open("!{sample_id}_assembly_report.html", "w") as f:
    f.write(html)
EOFPYTHON
    '''
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------
workflow {
    // Step 0: FastQC raw
    FASTQC_RAW(reads_ch)

    // Step 1: FASTP
    FASTP(reads_ch)

    // Step 2: bwa-mem2
    BWA_MAP(FASTP.out.trimmed, ref_ch.first())

    // Step 3: Extract mapped reads
    EXTRACT_MAPPED(BWA_MAP.out.bam)

    // Step 4: Shovill
    SHOVILL(EXTRACT_MAPPED.out.fastq)

    // Step 5: RagTag
    RAGTAG(SHOVILL.out.contigs, ref_ch.first())

    // Step 6: Pilon
    PILON(RAGTAG.out.scaffold, FASTP.out.trimmed)

    // Step 7: Coverage
    COVERAGE(PILON.out.consensus, FASTP.out.trimmed)

    // Step 8: FastQC trimmed
    FASTQC_TRIMMED(FASTP.out.trimmed)

    // Step 9: HTML Report
    REPORT(
        PILON.out.consensus,
        COVERAGE.out.depth
    )
}
