FROM mambaorg/micromamba:2.0.5

LABEL org.opencontainers.image.title="Viromics" \
      org.opencontainers.image.version="2.0.0" \
      org.opencontainers.image.description="Amplicon-based viral whole-genome sequencing and genomic epidemiology pipeline" \
      org.opencontainers.image.licenses="MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER pipeline/environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes && \
    rm -f /tmp/environment.yml

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH="/opt/conda/bin:$PATH"

RUN fastp --version && \
    fastplong --version && \
    minimap2 --version && \
    samtools --version | head -n 1 && \
    bcftools --version | head -n 1 && \
    tabix --version | head -n 1 && \
    python3 -c "import jsonschema; print(jsonschema.__version__)"
