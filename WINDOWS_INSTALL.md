# NGS Pipeline - Windows Installation Guide

## Overview

The NGS Pipeline is primarily designed for Linux/macOS systems, but can be run on Windows using several approaches. This guide will help you set up the pipeline on Windows.

## Option 1: Python Components Only (Recommended for Windows)

This approach installs the Python analysis components and uses Docker or WSL for bioinformatics tools.

### Prerequisites

1. **Anaconda or Miniconda**
   - Download from: https://docs.conda.io/en/latest/miniconda.html
   - Install and restart your terminal

2. **Git** (if not already installed)
   - Download from: https://git-scm.com/download/win
   - Or install via conda: `conda install git`

### Installation Steps

1. **Install Python Environment**
   ```powershell
   .\setup.ps1 install-deps
   ```

2. **Activate Environment**
   ```powershell
   conda activate ngs_pipeline
   ```

3. **Test Python Components**
   ```powershell
   python -c "import ngs_pipeline; print('NGS package installed successfully!')"
   ```

4. **Start Jupyter for Interactive Analysis**
   ```powershell
   jupyter lab notebooks/NGS_analysis_template.ipynb
   ```

## Option 2: Full Pipeline with Docker (Advanced)

For the complete bioinformatics pipeline including Nextflow workflows.

### Prerequisites

1. **Docker Desktop**
   - Download from: https://www.docker.com/products/docker-desktop
   - Install and start Docker Desktop

2. **Nextflow**
   - Download from: https://www.nextflow.io/docs/latest/getstarted.html#installation
   - Add to your PATH

### Installation Steps

1. **Build Docker Image**
   ```powershell
   docker build -t ngs_pipeline:latest docker/
   ```

2. **Run Pipeline in Docker**
   ```powershell
   docker run -it --rm -v ${PWD}:/workspace -w /workspace ngs_pipeline:latest nextflow run main.nf -profile docker
   ```

## Option 3: WSL2 with Ubuntu (Linux Compatibility)

For full Linux compatibility on Windows.

### Prerequisites

1. **Install WSL2**
   ```powershell
   wsl --install Ubuntu-22.04
   ```

2. **Restart and Set Up Ubuntu**
   - Follow Ubuntu setup prompts
   - Install conda in WSL2

### Installation Steps

1. **Enter WSL2**
   ```powershell
   wsl
   ```

2. **Run Linux Installation**
   ```bash
   make setup
   ```

## What's Available in Each Option

### Python Components Only ✅
- Data quality control analysis
- RNA-seq differential expression (DESeq2)
- Statistical analysis and visualization
- Interactive Jupyter notebooks
- HTML report generation
- Gene ontology analysis

### Missing on Windows (requires Docker/WSL) ❌
- Raw FASTQ processing (FastQC, trimming)
- Read alignment (STAR, Bowtie2)
- Peak calling (MACS2)
- ENDseq analysis (Bowtie2)
- Complete Nextflow pipeline execution

## Recommended Workflow for Windows Users

1. **Start with Python components** to analyze pre-processed data
2. **Use the Jupyter notebook** for interactive exploration
3. **Process raw data** using:
   - Galaxy Project (https://usegalaxy.org/) for web-based analysis
   - Cloud platforms (AWS, Google Cloud, Azure) for compute-intensive steps
   - Collaborate with Linux/macOS users for raw data processing

## File Structure After Installation

```
deseq2_claude/
├── ngs_pipeline/        # Python package
├── notebooks/              # Jupyter analysis templates
├── examples/               # Sample configurations
├── test/                   # Test data generators
├── templates/              # HTML report templates
└── setup.ps1              # Windows setup script
```

## Testing Your Installation

1. **Test Python Package**
   ```powershell
   .\setup.ps1 check-env
   ```

2. **Generate Test Data**
   ```powershell
   conda activate ngs_pipeline
   python test/generate_test_data.py
   ```

3. **Run Analysis Notebook**
   ```powershell
   jupyter lab notebooks/NGS_analysis_template.ipynb
   ```

## Troubleshooting

### Common Issues

1. **Conda environment creation fails**
   - Ensure you're using `environment_windows.yml`
   - Try: `conda clean --all` then retry

2. **Package import errors**
   - Activate environment: `conda activate ngs_pipeline`
   - Reinstall package: `pip install -e .`

3. **R packages missing**
   - Install manually: `conda install r-essentials`
   - Or use: `install.packages("package_name")` in R

### Getting Help

1. Check the full documentation in `README.md`
2. Review error logs in the terminal
3. Ensure all prerequisites are installed
4. Try the Docker option for full functionality

## Next Steps

Once installed, proceed with:
1. **Data preparation** - format your sample sheets
2. **Configuration** - adjust parameters in `config/`
3. **Analysis** - start with the Jupyter notebook template
4. **Results** - export publication-ready figures and reports

For questions or issues, please refer to the main documentation or consider using the Docker/WSL2 options for full pipeline functionality.