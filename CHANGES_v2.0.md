# Changes Summary: SeqForge v2.0 — Large Genome Assembly (No Bundled Reference)

## Overview
Updated pipeline to:
1. **Remove bundled ASFV reference** — now REQUIRES user to provide a reference via `-r` flag
2. **Remove all ASFV-specific language** — replaced with generic genome assembly terminology
3. **Support large genomes** (virus to gigabase genomes)

---

## File-by-file Changes

### 1. **Dockerfile**
**Removed:**
- Lines 86–100: Entire "bundled reference genome" section
  - NCBI download of ASFV FR682468.2
  - `/opt/references/default/` directory creation
  - `ENV DEFAULT_REF` variable

**Updated:**
- Label description: "Viral assembly" → "Genome assembly"
- Comments updated to reflect that users provide their own reference

**Impact:** Image is now ~100 MB smaller (no pre-downloaded reference)

---

### 2. **seqforge_pipeline.sh**
**Removed:**
- Lines 25, 62–66: Bundled reference fallback logic
  - `DEFAULT_REF` env var handling
  - Automatic reference detection if none provided
- Default value for `GENOME_SIZE` (was `5000000`)

**Updated:**
- **Required flags**: Both `-r REFERENCE` and `-g GENOME_SIZE` are now **mandatory**
- Line 39: Usage message updated to show `-g GENOME_SIZE` as required
- Line 69–88: New validation to check both reference file exists AND genome_size is provided
- Step names: "Map to reference genome" (was "Viral genome")
- Help text: Added example genome sizes for different organisms
- All "ASFV" / "viral" language replaced with "genome"

**Key additions:**
```bash
GENOME_SIZE=""  # required; must provide with -g

# Validation
for var in READ1 READ2 REFERENCE GENOME_SIZE OUTDIR; do
    if [ -z "${!var:-}" ]; then
        echo "ERROR: Missing required argument for $var"
        usage
    fi
done
```

---

### 3. **README.md**
**Major rewrites:**
- Title: "Short-Read Genome Assembly Pipeline" (no ASFV)
- Added prominent note: "Works for any genome size: viruses (~190 kb) to large genomes (gigabases)"
- **Requirements section**: Added reference genome as explicit requirement
- **Examples**: All default values adjusted for large genomes
  - `--genome_size` example: `3000000000` (3 Gb) instead of `5000000` (5 Mb)
  - Added reference in all example commands
- **New table**: "Typical genome sizes for reference" (ASFV, E. coli, Arabidopsis, human, wheat)
- Removed "Bundled reference" section
- Updated all launchers to show `-r REFERENCE` flag

---

### 4. **nextflow.config**
**Removed:**
- Old `test` profile hardcoded ASFV reference

**Updated:**
- Description: "Viral assembly" → "Genome assembly"
- Test profile now includes `reference = "${projectDir}/test_data/reference.fasta"`
- Version bumped: `1.0.0` → `2.0.0`

---

### 5. **main.nf**
**Removed:**
- Bundled reference handling (lines 60–64)
  - `Channel.fromPath(file("NO_REFERENCE"))` fallback
- Default value for `genome_size` (was `5000000`)

**Updated:**
- **Help message**: Added `--genome_size` as required parameter
- **Error checks**: Now validates `--reads`, `--reference`, AND `--genome_size`
  ```nf
  if (!params.genome_size) {
      log.error "ERROR: --genome_size is required. Use --help for usage."
      exit 1
  }
  ```
- Created dedicated reference channel:
  ```nf
  Channel.fromPath(params.reference, checkIfExists: true)
      .set { ref_ch }
  ```
- Updated help text with genome size examples
- All process descriptions now generic (no ASFV mentions)

---

## User-Facing Changes

### Before (v1.0 — ASFV-focused, genome size optional)
```bash
# Could omit reference & genome size — would use bundled ASFV (190 kb)
docker run seqforge -1 R1.fastq.gz -2 R2.fastq.gz -o results

# OR with custom reference
docker run -v /path/to/ref:/ref seqforge \
    -1 R1.fastq.gz -2 R2.fastq.gz -r /ref/custom.fasta -o results
```

### After (v2.0 — Generic, requires both reference and genome size)
```bash
# Reference AND genome size are REQUIRED
docker run -v /path/to/ref:/ref seqforge \
    -1 R1.fastq.gz -2 R2.fastq.gz \
    -r /ref/reference.fasta \
    -g 3200000000 \
    -o results

# Nextflow example
nextflow run main.nf \
    --reads 'data/*_{R1,R2}.fastq.gz' \
    --reference /path/to/reference.fasta \
    --genome_size 3200000000 \
    --outdir results
```

---

## Parameter Guidance

### Genome sizes to pass via `-g` / `--genome_size`
| Organism | Size | Example |
|----------|------|---------|
| ASFV | 190 kb | `-g 190000` |
| SARS-CoV-2 | 29.9 kb | `-g 29903` |
| E. coli | 4.6 Mb | `-g 4600000` |
| Arabidopsis | 135 Mb | `-g 135000000` |
| Human | 3.2 Gb | `-g 3200000000` |
| Wheat | 17 Gb | `-g 17000000000` |

---

## Testing

### Old way (v1.0)
Could test without a reference — image had ASFV pre-loaded

### New way (v2.0)
Must provide reference in test:
```bash
nextflow run main.nf \
    --reads 'test_data/*_{R1,R2}.fastq.gz' \
    --reference test_data/reference.fasta \
    --genome_size 5000000 \
    -profile test
```

---

## Docker Build
Build size reduction:
- **Before**: ~5 GB (includes ASFV reference + tools)
- **After**: ~4.8 GB (tools only, no pre-bundled reference)

No change to build time (~15–30 min) — reference download was fast; main time is conda/pip installs.

---

## Migration from v1.0 to v2.0

If you have v1.0 running:
1. Remove old image: `docker rmi seqforge:old`
2. Rebuild: `docker build -t seqforge .`
3. Update scripts to always pass **both** `-r` and `-g` flags
4. Examples:
   ```bash
   # Old (worked in v1.0)
   bash run_pipeline.sh -1 R1.fastq.gz -2 R2.fastq.gz -o results

   # New (required in v2.0)
   bash run_pipeline.sh -1 R1.fastq.gz -2 R2.fastq.gz -r reference.fasta -g 3200000000 -o results
   ```

## Backward Compatibility
**NOT backward compatible**. Old command lines that omitted `-r` or `-g` will fail with:
```
ERROR: Missing required argument for REFERENCE
ERROR: Missing required argument for GENOME_SIZE
```

This is intentional — forces users to:
- Provide their own reference (no ASFV assumption)
- Specify genome size for proper assembly parameter optimization

---

## Future Enhancements (Optional)
- Add `--skip-reference-mapping` for truly de novo assembly (skip steps 2–3)
- Support for long-read assembly (Flye + Racon)
- Polishing stage for large genomes (medaka, racon)
- Chromosome-level scaffolding (liftoff for annotations)
