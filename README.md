[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
# Genome Assembly Pipeline

FASTP → bwa-mem2 → Shovill/SPAdes → RagTag → Pilon → FastQC → Coverage → HTML Report

**You provide raw FASTQ files and a reference genome — everything else is automated.**

Works for any genome size: viruses (~190 kb) to large genomes (gigabases, e.g., plants, vertebrates).

---

## Requirements

| Requirement | Notes |
|-------------|-------|
| **Docker** | Only software needed on your machine |
| **Nextflow** | Optional — for scalable/HPC runs |
| **Reference genome** | FASTA file (required) |
| RAM | ≥ 16 GB recommended; scale with genome size |
| CPU | ≥ 6 cores recommended |
| Disk | ~20 GB per sample (adjust for your genome size) |

### Install Docker

| OS | Link |
|----|------|
| macOS | https://docs.docker.com/desktop/mac/install/ |
| Windows | https://docs.docker.com/desktop/windows/install/ |
| Linux (Ubuntu) | `sudo apt install docker.io && sudo systemctl start docker` |

### Install Nextflow (optional)

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
```

---

## Quick Start

### Step 1 — Build the Docker image (once only)

```bash
git clone https://github.com/lqhungpt/SeqForge.git
cd SeqForge
docker build -t seqforge .
```

> Build takes 15–30 minutes. Only needed once, or when Dockerfile changes.

---

## Option A — Run with Nextflow (recommended)

Nextflow handles parallelism, checkpointing, and HPC submission automatically.

### Single sample

```bash
nextflow run lqhungpt/SeqForge \
    -r main \
    --reads '/path/to/FASTQ_{R1,R2}.fastq.gz' \
    --reference /path/to/reference.fasta \
    --genome_size 3200000000 \
    --outdir /path/to/results \
    -profile docker \
    -resume
```

### Multiple samples at once

```bash
nextflow run lqhungpt/SeqForge \
    -r main \
    --reads '/path/to/data/*_{R1,R2}.fastq.gz' \
    --reference /path/to/reference.fasta \
    --genome_size 3200000000 \
    --outdir /path/to/results \
    -profile docker \
    -resume
```

### On HPC (SLURM)

```bash
nextflow run lqhungpt/SeqForge \
    -r main \
    --reads '/path/to/data/*_{R1,R2}.fastq.gz' \
    --reference /path/to/reference.fasta \
    --genome_size 3200000000 \
    --outdir /path/to/results \
    -profile slurm \
    -resume
```

### Nextflow parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--reads` | FASTQ glob pattern (single quotes required) | `/data/*_{R1,R2}.fastq.gz` |
| `--reference` | Reference genome FASTA (required) | `/data/reference.fasta` |
| `--genome_size` | Estimated genome size in bp (required) | `29903` (SARS-CoV-2), `3200000000` (human) |
| `--outdir` | Output directory | `/data/results` |
| `--depth` | Assembly depth (×) | `50`, `100`, `200` |
| `--threads` | CPU threads per task | `4`, `8`, `16` |

### Nextflow profiles

| Profile | Use case |
|---------|----------|
| `docker` | Local machine (macOS, Linux, Windows) |
| `singularity` | HPC without Docker |
| `slurm` | SLURM cluster |
| `pbs` | PBS cluster |
| `test` | Quick test run (small dataset) |

---

## Option B — Run with launcher scripts (macOS / Linux / Windows)

### macOS / Linux

```bash
bash run_pipeline.sh \
    -1 /path/to/sample_R1.fastq.gz \
    -2 /path/to/sample_R2.fastq.gz \
    -r /path/to/reference.fasta \
    -g 3200000000 \
    -o /path/to/results/sample \
    -p sample \
    -t 6
```

### Windows PowerShell

```powershell
.\run_pipeline.ps1 `
    -Read1 C:\data\sample_R1.fastq.gz `
    -Read2 C:\data\sample_R2.fastq.gz `
    -Reference C:\refs\reference.fasta `
    -GenomeSize 3200000000 `
    -OutDir C:\results\sample `
    -Prefix sample `
    -Threads 6
```

### Windows Command Prompt

```bat
run_pipeline.bat -1 C:\data\R1.fastq.gz -2 C:\data\R2.fastq.gz -r C:\refs\reference.fasta -g 3200000000 -o C:\results\sample -p sample -t 6
```

### Launcher script parameters

| Flag | Description | Required | Example |
|------|-------------|----------|---------|
| `-1` | R1 FASTQ file (gzipped) | **YES** | `/data/R1.fastq.gz` |
| `-2` | R2 FASTQ file (gzipped) | **YES** | `/data/R2.fastq.gz` |
| `-r` | Reference genome FASTA | **YES** | `/data/reference.fasta` |
| `-g` | Genome size (bp) | **YES** | `29903`, `4600000`, `3200000000` |
| `-o` | Output directory | **YES** | `/results/sample` |
| `-p` | Sample prefix/name | NO | `sample` |
| `-t` | CPU threads | NO | `6` |
| `-d` | Assembly depth (×) | NO | `150` |

---

## Pipeline steps

| Step | Tool | Description |
|------|------|-------------|
| 0 | FastQC | Raw read quality assessment |
| 1 | FASTP | Adapter trimming, quality filtering |
| 2 | bwa-mem2 | Map trimmed reads to reference |
| 3 | samtools | Extract mapped reads |
| 4 | Shovill/SPAdes | De novo assembly (k=71,91,111) |
| 5 | RagTag | Scaffold contigs against reference |
| 6 | Pilon | Polish consensus, fill gaps |
| 7 | bwa-mem2 + samtools | Coverage depth map |
| 8 | FastQC | Post-trimming read quality |

---

## Output structure

```
results/SAMPLE/
├── 00_fastqc_raw/               # FastQC raw reads
├── 01_trimmed/                  # FASTP trimmed reads + QC report
├── 02_mapped/                   # BWA-mem2 BAM (reference mapping)
├── 03_mapped_fastq/             # Extracted mapped reads
├── 04_shovill_denovo/           # SPAdes contigs
├── 05_ragtag_scaffold/          # RagTag scaffold
├── 06_final_consensus/          # Pilon FASTA + VCF + changes
├── 07_coverage_map/             # Coverage depth + plot
├── 08_fastqc_final/             # FastQC trimmed reads
├── logs/                        # Step-by-step logs
└── SAMPLE_assembly_report.html  # Interactive HTML report
```

---

## Checkpoint / resume

Each step is skipped if output already exists.

**With Nextflow** — use `-resume` flag (automatic).

**With launcher scripts** — delete the folder to rerun a step:

```bash
# Rerun from de novo assembly onwards
rm -rf results/sample/04_shovill_denovo/
rm -rf results/sample/05_ragtag_scaffold/
rm -rf results/sample/06_final_consensus/
rm -f  results/sample/sample_assembly_report.html
```

---

## Typical genome sizes for reference

| Organism | Genome size |
|----------|-------------|
| ASFV, influenza | 190 kb, 13 kb |
| SARS-CoV-2 | 29.9 kb |
| *E. coli* | 4.6 Mb |
| *Arabidopsis* | 135 Mb |
| *Drosophila* | 140 Mb |
| Human | 3.2 Gb |
| Wheat | 17 Gb |

Pass `-g` / `--genome_size` matching your organism. This helps Shovill estimate coverage and optimize assembly parameters.

---

## Citation

If you use this pipeline please cite the tools:

- **SPAdes**: Prjibelski et al. (2020) *Current Protocols in Bioinformatics* doi:10.1002/cpbi.102
- **Pilon**: Walker et al. (2014) *PLOS ONE* doi:10.1371/journal.pone.0112963
- **bwa-mem2**: Vasimuddin et al. (2019) *IPDPS* doi:10.1109/IPDPS.2019.00041
- **RagTag**: Alonge et al. (2022) *Genome Biology* doi:10.1186/s13059-022-02823-7
- **FASTP**: Chen et al. (2018) *Bioinformatics* doi:10.1093/bioinformatics/bty560
- **FastQC**: Andrews (2010) https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

---

## Author

Hung Luong — University of Nebraska-Lincoln
GitHub: [@lqhungpt](https://github.com/lqhungpt)
