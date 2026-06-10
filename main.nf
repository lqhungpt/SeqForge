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

    script:
    """
    # Calculate stats
    CONSENSUS_LEN=\$(grep -v ">" ${consensus} | tr -d '\\n' | wc -c || echo 0)
    CONSENSUS_GC=\$(grep -v ">" ${consensus} | tr -d '\\n' | awk '{g=gsub(/[GCgc]/,""); t=length(\$0); printf "%.2f", (g/t)*100}')
    CONSENSUS_N=\$(grep -v ">" ${consensus} | tr -cd 'N' | wc -c || echo 0)
    AVG_DEPTH=\$(awk '{sum+=\$3; n++} END {printf "%.1f", sum/n}' ${depth} 2>/dev/null || echo 0)
    MIN_DEPTH=\$(awk 'NR==1{m=\$3} \$3<m{m=\$3} END{print m}' ${depth} 2>/dev/null || echo 0)
    MAX_DEPTH=\$(awk 'BEGIN{m=0} \$3>m{m=\$3} END{print m}' ${depth} 2>/dev/null || echo 0)
    COVERED_PCT=\$(awk -v len="\$CONSENSUS_LEN" '\$3>0{n++} END{printf "%.2f", (n/len)*100}' ${depth} 2>/dev/null || echo 0)
    LOW_DEPTH=\$(awk '\$3<10{n++} END{print n+0}' ${depth} 2>/dev/null || echo 0)

    python3 /usr/local/bin/generate_report.py \\
        --prefix          ${sample_id} \\
        --outdir          . \\
        --depth           ${depth} \\
        --consensus_len   \$CONSENSUS_LEN \\
        --consensus_gc    \$CONSENSUS_GC \\
        --consensus_n     \$CONSENSUS_N \\
        --avg_depth       \$AVG_DEPTH \\
        --min_depth       \$MIN_DEPTH \\
        --max_depth       \$MAX_DEPTH \\
        --covered_pct     \$COVERED_PCT \\
        --low_depth_bases \$LOW_DEPTH \\
        --n_count         \$CONSENSUS_N \\
        --contig_count    0 \\
        --snps_corrected  0 \\
        --gaps_filled     0 \\
        --final_mapped    0 \\
        --elapsed_min     0 \\
        --elapsed_sec     0 \\
        --threads         ${params.threads} \\
        --fastqc_raw_dir  . \\
        --fastqc_final_dir .
    """
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
