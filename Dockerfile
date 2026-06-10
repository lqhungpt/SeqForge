# =============================================================================
# Genome Assembly Pipeline - Docker Image
# Tools: fastp, bwa-mem2, samtools>=1.17, shovill, spades, pilon,
#        ragtag, minimap2, liftoff, multiqc, Draw_SequencingDepth
# =============================================================================
FROM mambaorg/micromamba:1.5.8

LABEL maintainer="hungluong"
LABEL description="Genome assembly: FASTP -> bwa-mem2 -> Shovill/SPAdes -> RagTag -> Pilon -> Liftoff -> MultiQC"

USER root

# System deps
RUN apt-get update && apt-get install -y \
    curl wget procps default-jre \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Install all tools in one environment — no prokka conflict
# ---------------------------------------------------------------------------
RUN micromamba install -y -n base \
    -c bioconda -c conda-forge \
    fastp \
    bwa-mem2 \
    bwa \
    "samtools>=1.17" \
    spades \
    pilon \
    flash2 \
    lighter \
    samclip \
    seqtk \
    pigz \
    trimmomatic \
    csvtk \
    megahit \
    minimap2 \
    mummer4 \
    ragtag \
    liftoff \
    fastqc \
    multiqc \
    perl \
    && micromamba clean --all --yes

# Ensure conda bins are in PATH
ENV PATH="/usr/local/bin:/opt/conda/bin:$PATH"

# Verify key tools and samtools version
RUN which fastp && which bwa-mem2 && which samtools && which spades.py && \
    which pilon && which ragtag.py && which minimap2 && which liftoff && which multiqc
RUN samtools version | head -1

# ---------------------------------------------------------------------------
# Install Python plotting dependencies
# ---------------------------------------------------------------------------
RUN pip install --no-cache-dir --root-user-action=ignore \
    numpy matplotlib pandas

# ---------------------------------------------------------------------------
# Install shovill 1.1.0
# ---------------------------------------------------------------------------
RUN curl -fsSL https://raw.githubusercontent.com/tseemann/shovill/v1.1.0/bin/shovill \
    -o /usr/local/bin/shovill \
    && chmod +x /usr/local/bin/shovill

# Shovill adapter DB
RUN mkdir -p /usr/local/db \
    && curl -fsSL https://raw.githubusercontent.com/tseemann/shovill/v1.1.0/db/trimmomatic.fa \
    -o /usr/local/db/trimmomatic.fa

# Symlink flash2 -> flash for shovill
RUN ln -sf $(which flash2) /usr/local/bin/flash

# Stubs for unused assemblers (shovill dependency check)
RUN printf '#!/usr/bin/env bash\necho "KMC ver. 3.2.1"\n'  > /usr/local/bin/kmc             && chmod +x /usr/local/bin/kmc
RUN printf '#!/usr/bin/env bash\necho "1.2.9"\n'            > /usr/local/bin/megahit_toolkit  && chmod +x /usr/local/bin/megahit_toolkit
RUN printf '#!/usr/bin/env bash\necho "SKESA 2.5.1"\n'      > /usr/local/bin/skesa            && chmod +x /usr/local/bin/skesa
RUN printf '#!/usr/bin/env bash\necho "Version 1.2.10"\n'   > /usr/local/bin/velvetg          && chmod +x /usr/local/bin/velvetg
RUN printf '#!/usr/bin/env bash\necho "Version 1.2.10"\n'   > /usr/local/bin/velveth          && chmod +x /usr/local/bin/velveth

# Patch shovill db path
RUN sed -i 's|$FindBin::RealBin/../db/|/usr/local/db/|g' /usr/local/bin/shovill

# ---------------------------------------------------------------------------
# Copy pipeline scripts
# ---------------------------------------------------------------------------
COPY seqforge_pipeline.sh /usr/local/bin/seqforge_pipeline.sh
COPY Draw_SequencingDepth.py /usr/local/bin/Draw_SequencingDepth.py
COPY generate_report.py /usr/local/bin/generate_report.py
RUN chmod +x /usr/local/bin/seqforge_pipeline.sh \
    && chmod +x /usr/local/bin/Draw_SequencingDepth.py \
    && chmod +x /usr/local/bin/generate_report.py

WORKDIR /data

# For direct Docker use: docker run seqforge -1 R1 -2 R2 -r reference.fasta ...
# Default entrypoint: pipeline script
ENTRYPOINT ["/usr/local/bin/seqforge_pipeline.sh"]
