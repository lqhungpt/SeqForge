# SeqForge v2.0 — Local Testing Guide

Quick steps to test the updated pipeline on your local machine.

---

## Prerequisites

- Docker installed and running
- Nextflow installed (optional, for Nextflow tests)
- ~30 GB free disk space for build + test run
- 15+ GB RAM available

---

## Step 1: Build Docker Image

```bash
cd /path/to/seqforge_v2.0
docker build -t seqforge:v2.0 .
```

Expected output:
```
Successfully built [hash]
Successfully tagged seqforge:v2.0
```

**Troubleshooting:**
- If build fails on conda installs, network issue likely — retry
- If `shovill` download fails, check internet connection
- If memory error: reduce background processes

Build takes **15–30 minutes** depending on internet speed and CPU.

---

## Step 2: Prepare Test Data

For testing, you need:
- 2 FASTQ files (R1 and R2) — can be small subset
- 1 reference FASTA file

### Option A: Use Public Test Data (SARS-CoV-2)

```bash
mkdir -p ~/seqforge_test/data ~/seqforge_test/refs

# Download small SARS-CoV-2 reads (10K reads, ~600 MB)
cd ~/seqforge_test/data
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra30/SRR/SRR035/SRR035417/SRR035417_1.fastq.gz
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra30/SRR/SRR035/SRR035417/SRR035417_2.fastq.gz

# Download SARS-CoV-2 reference
cd ~/seqforge_test/refs
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/805/GCF_009858805.2_ASM985880v2/GCF_009858805.2_ASM985880v2_genomic.fna.gz
gunzip GCF_009858805.2_ASM985880v2_genomic.fna.gz
mv GCF_009858805.2_ASM985880v2_genomic.fna reference.fasta
```

### Option B: Use Your Own Data

```bash
# Just copy/symlink your FASTQ and reference
mkdir -p ~/seqforge_test/data ~/seqforge_test/refs
cp /your/data/*_R1.fastq.gz ~/seqforge_test/data/
cp /your/data/*_R2.fastq.gz ~/seqforge_test/data/
cp /your/reference.fasta ~/seqforge_test/refs/
```

---

## Step 3: Test with Docker (Direct)

### Single sample test

```bash
cd ~/seqforge_test

docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/refs:/ref \
  -v $(pwd)/results:/results \
  seqforge:v2.0 \
  -1 /data/SRR035417_1.fastq.gz \
  -2 /data/SRR035417_2.fastq.gz \
  -r /ref/reference.fasta \
  -g 29903 \
  -o /results/test_sample \
  -p test_sample \
  -t 4 \
  -d 100
```

**Expected output:**
```
╔══════════════════════════════════════════════════╗
║         Genome Assembly Pipeline (Docker)        ║
╠══════════════════════════════════════════════════╣
║  Sample      : test_sample                       ║
║  Threads     : 4                                 ║
║  Genome size : 29903 bp                          ║
║  Assembly depth: 100×                            ║
...
╚══════════════════════════════════════════════════╝

STEP 0/9 · FastQC on raw reads
    [DONE]  FastQC raw reads complete
...
STEP 9/9 · Generate HTML report
    [DONE]  HTML assembly report written

Pipeline complete. Results in: /results/test_sample
```

**Runtime:** ~5–15 minutes depending on data size & CPU.

### Check results

```bash
ls -lh ~/seqforge_test/results/test_sample/

# Expected output folders:
# 00_fastqc_raw/
# 01_trimmed/
# 02_mapped/
# 03_mapped_fastq/
# 04_shovill_denovo/
# 05_ragtag_scaffold/
# 06_final_consensus/
# 07_coverage_map/
# 08_fastqc_final/
# logs/
# test_sample_assembly_report.html

# Open HTML report
open ~/seqforge_test/results/test_sample/test_sample_assembly_report.html
```

---

## Step 4: Test with Nextflow (Optional)

Create test config:

```bash
cd ~/seqforge_test

# Copy Nextflow files
cp /path/to/seqforge_v2.0/main.nf .
cp /path/to/seqforge_v2.0/nextflow.config .

# Run Nextflow test
nextflow run main.nf \
  --reads 'data/*_{1,2}.fastq.gz' \
  --reference refs/reference.fasta \
  --genome_size 29903 \
  --outdir results_nextflow \
  -profile docker \
  -resume
```

**Note:** Replace `{1,2}` with actual R1/R2 pattern from your files.

---

## Step 5: Verify Required Parameters

Test that the pipeline **rejects** missing parameters:

### Missing reference (should fail)
```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  seqforge:v2.0 \
  -1 /data/SRR035417_1.fastq.gz \
  -2 /data/SRR035417_2.fastq.gz \
  -o /results/test \
  -g 29903

# Expected error:
# ERROR: Missing required argument for REFERENCE (use -h for help)
```

### Missing genome size (should fail)
```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/refs:/ref \
  -v $(pwd)/results:/results \
  seqforge:v2.0 \
  -1 /data/SRR035417_1.fastq.gz \
  -2 /data/SRR035417_2.fastq.gz \
  -r /ref/reference.fasta \
  -o /results/test

# Expected error:
# ERROR: Missing required argument for GENOME_SIZE (use -h for help)
```

✅ **If both errors appear → v2.0 is correctly enforcing required parameters**

---

## Step 6: Check Help Output

```bash
# Docker help
docker run --rm seqforge:v2.0 -h

# Nextflow help
nextflow run main.nf --help

# Both should show genome size examples:
# SARS-CoV-2: 29903
# E. coli: 4600000
# Human: 3200000000
# Wheat: 17000000000
```

✅ **If examples appear → v2.0 documentation is updated**

---

## Common Issues & Fixes

### Docker image not found
```bash
# Rebuild with correct tag
docker build -t seqforge:v2.0 .

# List images to verify
docker images | grep seqforge
```

### Permission denied on output folder
```bash
# Fix ownership (Linux/macOS)
sudo chown -R $(whoami):$(whoami) ~/seqforge_test/results/
```

### Out of disk space
```bash
# Clean old Docker data
docker system prune -a

# Check disk usage
df -h
```

### Pipeline hangs at assembly step
- Increase `-d` (assembly depth) if small dataset
- Or lower it if you're on limited RAM
- Check `logs/*shovill.log` for details

### Coverage plot not generated
- Ensure final consensus FASTA has sequence data
- Check `logs/*_pilon.log` for errors
- Manually verify: `samtools stats results/*/07_coverage_map/*final*.bam`

---

## Sample Test Results (Expected)

For **SARS-CoV-2** (29,903 bp) with 100× depth target:
- Consensus length: ~29,800–29,900 bp (should match reference)
- GC content: ~37–38%
- Remaining Ns: 0–100 (depending on coverage gaps)
- Average depth: ~100–150×
- Contigs (de novo): 1–5 (depending on breaks)
- Runtime: 5–15 minutes on 4 CPU cores

For **E. coli** (4.6 Mb):
- Consensus length: ~4,600,000 bp
- GC content: ~50–51%
- Runtime: 30–60 minutes on 6 CPU cores

---

## Next Steps After Testing

✅ If all tests pass:
1. Update documentation with any local findings
2. Test on your actual data (not public datasets)
3. Commit to GitHub: `git add -A && git commit -m "Release v2.0"`
4. Tag release: `git tag v2.0 && git push --tags`
5. Push Docker image: `docker tag seqforge:v2.0 yourusername/seqforge:v2.0 && docker push yourusername/seqforge:v2.0`

❌ If tests fail:
1. Check logs in `results/*/logs/`
2. Report errors with:
   - Docker version: `docker --version`
   - OS/CPU info: `uname -a`
   - Nextflow version (if testing Nextflow): `nextflow -version`
3. Run with verbose mode: add `-v` flag to bash script

---

## Quick Ref: Test Commands

```bash
# Build
docker build -t seqforge:v2.0 .

# Single sample (SARS-CoV-2)
docker run --rm -v $(pwd)/data:/data -v $(pwd)/refs:/ref -v $(pwd)/results:/results seqforge:v2.0 \
  -1 /data/*_1.fastq.gz -2 /data/*_2.fastq.gz -r /ref/reference.fasta -g 29903 -o /results/test -t 4

# Nextflow single
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --reference refs/reference.fasta --genome_size 29903 -profile docker

# Test help
docker run --rm seqforge:v2.0 -h
nextflow run main.nf --help

# View results
open results/*/test_assembly_report.html
```
