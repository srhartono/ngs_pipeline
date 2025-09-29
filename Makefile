# XRN2 Multiomics Pipeline Makefile
# Convenience targets for setup, testing, and execution

.PHONY: help setup clean test demo install-deps build-refs run-local run-slurm check-env

# Default target
help:
	@echo "XRN2 Multiomics Pipeline"
	@echo "========================"
	@echo ""
	@echo "Available targets:"
	@echo "  setup        - Install dependencies and setup environment"
	@echo "  install-deps - Install conda/mamba dependencies"
	@echo "  build-refs   - Download and build reference files"
	@echo "  test         - Run quick test with toy data"
	@echo "  demo         - Run full demo pipeline"
	@echo "  run-local    - Run pipeline locally"
	@echo "  run-slurm    - Run pipeline on SLURM cluster"
	@echo "  check-env    - Check environment and dependencies"
	@echo "  clean        - Clean up temporary files"
	@echo "  clean-all    - Clean everything including results"

# Setup everything
setup: install-deps build-refs
	@echo "Setup complete! Ready to run pipeline."

# Install dependencies
install-deps:
	@echo "Installing conda environment..."
	@if command -v mamba >/dev/null 2>&1; then \
		mamba env create -f env/environment.yml -n xrn2_multiomics || mamba env update -f env/environment.yml -n xrn2_multiomics; \
	elif command -v conda >/dev/null 2>&1; then \
		conda env create -f env/environment.yml -n xrn2_multiomics || conda env update -f env/environment.yml -n xrn2_multiomics; \
	else \
		echo "Error: Neither mamba nor conda found. Please install Anaconda/Miniconda first."; \
		exit 1; \
	fi
	@echo "Installing Python package in development mode..."
	@bash -c "source $$(conda info --base)/etc/profile.d/conda.sh && conda activate xrn2_multiomics && pip install -e ."

# Build reference files
build-refs:
	@echo "Downloading and building reference files..."
	@bash references/get_references.sh --hs1 --build-index star bowtie2 bismark --dbsnp
	@echo "Reference files ready."

# Quick test with toy data
test: check-env
	@echo "Running quick test with toy data..."
	@bash -c "source $$(conda info --base)/etc/profile.d/conda.sh && conda activate xrn2_multiomics && \
		nextflow run main.nf -profile test -with-report test_report.html -with-timeline test_timeline.html"
	@echo "Test completed successfully!"

# Full demo pipeline
demo: check-env
	@echo "Running full demo pipeline..."
	@bash -c "source $$(conda info --base)/etc/profile.d/conda.sh && conda activate xrn2_multiomics && \
		nextflow run main.nf -profile local --samples config/samplesheet.tsv --enable_tcga false --enable_string false \
		-with-report demo_report.html -with-timeline demo_timeline.html"
	@echo "Demo completed successfully!"

# Run pipeline locally
run-local: check-env
	@echo "Running pipeline locally..."
	@bash -c "source $$(conda info --base)/etc/profile.d/conda.sh && conda activate xrn2_multiomics && \
		nextflow run main.nf -profile local --samples config/samplesheet.tsv"

# Run pipeline on SLURM
run-slurm: check-env
	@echo "Running pipeline on SLURM cluster..."
	@bash -c "source $$(conda info --base)/etc/profile.d/conda.sh && conda activate xrn2_multiomics && \
		nextflow run main.nf -profile slurm --samples config/samplesheet.tsv"

# Check environment
check-env:
	@echo "Checking environment and dependencies..."
	@bash -c "source $$(conda info --base)/etc/profile.d/conda.sh && conda activate xrn2_multiomics && \
		python -c 'import xrn2_multiomics; print(\"Python package: OK\")' && \
		nextflow -version && \
		echo 'Environment check: PASSED'"

# Clean temporary files
clean:
	@echo "Cleaning temporary files..."
	@rm -rf work/
	@rm -rf .nextflow*
	@rm -f *.html
	@rm -f *.svg
	@rm -f trace.txt
	@rm -f timeline.html
	@rm -f report.html
	@find . -name "*.pyc" -delete
	@find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	@echo "Cleanup complete."

# Clean everything including results
clean-all: clean
	@echo "Cleaning all results..."
	@rm -rf results/
	@rm -rf references/hs1/
	@rm -rf references/indices/
	@echo "Full cleanup complete."

# Development targets
lint:
	@echo "Running code linting..."
	@bash -c "source $$(conda info --base)/etc/profile.d/conda.sh && conda activate xrn2_multiomics && \
		python -m flake8 xrn2_multiomics/ tests/ && \
		python -m black --check xrn2_multiomics/ tests/"

format:
	@echo "Formatting code..."
	@bash -c "source $$(conda info --base)/etc/profile.d/conda.sh && conda activate xrn2_multiomics && \
		python -m black xrn2_multiomics/ tests/"

# Docker targets
docker-build:
	@echo "Building Docker image..."
	@docker build -t xrn2_multiomics:latest docker/

docker-run:
	@echo "Running pipeline in Docker..."
	@docker run --rm -v $$(pwd):/workspace -w /workspace xrn2_multiomics:latest \
		nextflow run main.nf -profile docker --samples config/samplesheet.tsv

# Utility targets
show-config:
	@echo "Current configuration:"
	@bash -c "source $$(conda info --base)/etc/profile.d/conda.sh && conda activate xrn2_multiomics && \
		nextflow config main.nf -profile local"

show-dag:
	@echo "Generating workflow DAG..."
	@bash -c "source $$(conda info --base)/etc/profile.d/conda.sh && conda activate xrn2_multiomics && \
		nextflow run main.nf --samples config/samplesheet.tsv -with-dag flowchart.html"
	@echo "DAG saved to flowchart.html"

# Version information
version:
	@echo "XRN2 Multiomics Pipeline v1.0.0"
	@bash -c "source $$(conda info --base)/etc/profile.d/conda.sh && conda activate xrn2_multiomics && \
		nextflow -version && \
		python --version"