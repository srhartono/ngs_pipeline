# XRN2 Multiomics Pipeline Docker Image
# Multi-stage build for optimized image size

# Stage 1: Base environment with conda
FROM mambaorg/micromamba:1.5.1 as base

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    pkg-config \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Set up micromamba environment
COPY env/environment.yml /tmp/environment.yml

USER $MAMBA_USER

# Create conda environment
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Activate conda environment by default
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Stage 2: Add pipeline code
FROM base as pipeline

USER root

# Set working directory
WORKDIR /opt/ngs_pipeline

# Copy pipeline code
COPY . .

# Install Python package
RUN pip install -e .

# Create directories
RUN mkdir -p /workspace /data /results && \
    chown -R $MAMBA_USER:$MAMBA_USER /opt/ngs_pipeline /workspace /data /results

# Switch to non-root user
USER $MAMBA_USER

# Set up environment variables
ENV PATH="/opt/conda/bin:$PATH"
ENV PYTHONPATH="/opt/ngs_pipeline:$PYTHONPATH"
ENV R_LIBS_USER="/opt/conda/lib/R/library"

# Expose common ports
EXPOSE 8080 8888

# Set entrypoints
COPY docker/entrypoint.sh /opt/entrypoint.sh
USER root
RUN chmod +x /opt/entrypoint.sh
USER $MAMBA_USER

ENTRYPOINT ["/opt/entrypoint.sh"]
CMD ["bash"]

# Metadata
LABEL maintainer="XRN2 Team <xrn2-team@example.com>"
LABEL description="XRN2 Multiomics Pipeline - Comprehensive analysis of multiomics data"
LABEL version="1.0.0"
LABEL org.opencontainers.image.source="https://github.com/xrn2-team/ngs_pipeline"