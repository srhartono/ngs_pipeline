#!/bin/bash
# Docker entrypoint script for ngs Pipeline

set -e

# Activate conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate base

# Set up R library path
export R_LIBS_USER="/opt/conda/lib/R/library"

# Function to display help
show_help() {
    cat << EOF
ngs Pipeline Docker Container

Available commands:
  nextflow    - Run Nextflow pipeline
  python      - Run Python CLI
  R           - Start R session
  Rscript     - Run R script
  bash        - Start bash shell
  jupyter     - Start Jupyter Lab (use with -p 8888:8888)
  test        - Run pipeline tests
  help        - Show this help

Examples:
  # Run full pipeline
  docker run -v /path/to/data:/data -v /path/to/results:/results \\
    ngs_pipeline nextflow run /opt/xrn2_multiomics/main.nf \\
    --samples /data/samplesheet.tsv --outdir /results

  # Run Python CLI
  docker run -v /path/to/data:/data \\
    ngs_pipeline python -m xrn2_multiomics qc /data

  # Interactive session
  docker run -it -v /path/to/data:/data \\
    ngs_pipeline bash

  # Jupyter Lab
  docker run -p 8888:8888 -v /path/to/data:/data \\
    ngs_pipeline jupyter

EOF
}

# Handle command
case "$1" in
    nextflow)
        shift
        exec nextflow "$@"
        ;;
    python)
        shift
        exec python "$@"
        ;;
    R)
        shift
        exec R "$@"
        ;;
    Rscript)
        shift
        exec Rscript "$@"
        ;;
    jupyter)
        shift
        cd /workspace
        exec jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --allow-root "$@"
        ;;
    test)
        cd /opt/ngs_pipeline
        exec make test
        ;;
    help|--help|-h)
        show_help
        ;;
    bash|shell)
        exec /bin/bash
        ;;
    "")
        # No command provided, show help and start bash
        show_help
        echo ""
        echo "Starting interactive bash session..."
        exec /bin/bash
        ;;
    *)
        # Unknown command, try to execute it
        exec "$@"
        ;;
esac