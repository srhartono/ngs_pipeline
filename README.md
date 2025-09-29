# NGS Pipeline

A comprehensive bioinformatics pipeline for analyzing RNA-seq, EU-seq, sDRIP-seq, and ENDseq data for multi-condition NGS experiments.

## Overview

This pipeline processes multiomics NGS data with flexible experimental designs supporting:
- Multiple treatment conditions (Treat1, Treat2, etc.)
- Multiple control groups (Control1, Control2, etc.)
- Batch effects and replicates
- Project-based organization
- Configurable comparisons

## Features

- **Multi-assay support**: RNA-seq, EU-seq, sDRIP-seq, ENDseq
- **Flexible experimental design**: Support for multiple treatments, controls, batches, and replicates
- **Project organization**: Structured data/results/notebooks directories by project
- **Dual workflows**: Both Nextflow and Snakemake implementations
- **Comprehensive QC**: FastQC, MultiQC, adapter trimming, contamination checks
- **Assay-specific mapping**: STAR for RNA/EU-seq, Bowtie2 for sDRIP/ENDseq
- **Peak calling**: MACS2 with assay-appropriate parameters
- **Differential analysis**: DESeq2 with flexible comparison specifications
- **Integrative analysis**: PCA, clustering, pathway analysis across assays
- **Visualization**: Interactive plots, bigWig tracks, IGV sessions
- **Reporting**: Comprehensive analysis reports and dashboards

## Quick Start

### 1. Setup Environment

```bash
# Clone and setup
git clone <repository>
cd ngs_pipeline
pip install -e .

# Download references (requires ~50GB disk space)
bash references/get_references.sh --hs1 --build-index star bowtie2 --dbsnp
```

### 2. Configure Project

Create or edit `config/project_config.yml`:

```yaml
project:
  name: "project1"
  
samples:
  treatments: ["Treat1", "Treat2"]
  controls: ["Control1", "Control2"]
  batches: ["batch1", "batch2"]
  replicates: [1, 2]

assays: ["rnaseq", "euseq", "sdripseq", "endseq"]

comparisons:
  primary:
    - name: "Treat1_vs_Control1"
      case: "Treat1"
      control: "Control1"
      description: "Primary treatment 1 vs control 1"
```

### 3. Prepare Sample Data

Organize your data with the new naming scheme:
- Sample format: `{condition}_{batch}_rep{replicate}_{assay}.fastq.gz`
- Examples: `Treat1_batch1_rep1_rnaseq.fastq.gz`, `Control1_batch2_rep2_euseq.fastq.gz`

Place FASTQ files in: `data/{project_name}/`

### 4. Run Workflows

#### Nextflow (Original)
```bash
# Local execution
nextflow run main.nf -profile local --samples config/samplesheet.tsv

# Cluster execution
nextflow run main.nf -profile slurm --samples config/samplesheet.tsv
```

#### Snakemake (New Parallel Implementation)
```bash
# Local execution
snakemake -s workflow/Snakefile --configfile config/project_config.yml --cores 8

# Cluster execution
snakemake -s workflow/Snakefile --configfile config/project_config.yml --profile slurm
```

### 5. Analysis Notebooks

Interactive analysis with Jupyter notebooks:
```bash
# Start analysis environment
jupyter lab notebooks/{project_name}/

# Use template
cp notebooks/project1/ngs_analysis_template.ipynb notebooks/{your_project}/
```

### 6. Test Run

```bash
# Generate test data
python -m ngs_pipeline.test.generate_test_data data/test_project/

# Quick test with toy data
ngs_pipeline test --project test_project

# Run pipeline on test data
snakemake -s workflow/Snakefile --configfile config/project_config.yml --cores 4
```

## Input Requirements

### Required Files
- FASTQ files (gzipped) with proper naming convention
- Project configuration YAML file

### Reference Files (auto-downloaded)
- hs1 (T2T/CHM13) genome FASTA
- hs1 annotation GTF
- dbSNP VCF for variant annotation
- Pre-built indices for STAR, Bowtie2, Bismark

## Output Structure

```
results/
├── qc/
│   ├── fastqc/
│   ├── multiqc/
│   └── trimming/
├── mapping/
│   ├── star/
│   ├── bowtie2/
│   └── bismark/
├── tracks/
│   ├── bigwig/
│   └── igv_sessions/
├── peaks/
│   ├── macs2/
│   └── frip_scores/
├── quantification/
│   ├── counts/
## Project Structure

The pipeline organizes data and results by project:

```
ngs_pipeline/
├── config/
│   └── project_config.yml      # Project configuration
├── data/
│   ├── project1/               # Raw FASTQ files
│   │   ├── Treat1_batch1_rep1_rnaseq.fastq.gz
│   │   ├── Treat1_batch1_rep1_euseq.fastq.gz
│   │   └── ...
│   └── project2/
├── results/
│   ├── project1/               # Analysis results
│   │   ├── qc/                 # Quality control
│   │   ├── mapping/            # Alignment results
│   │   ├── peaks/              # Peak calling
│   │   ├── quantification/     # Gene/peak counts
│   │   ├── differential/       # DE analysis
│   │   └── reports/            # Final reports
│   └── project2/
├── notebooks/
│   ├── project1/               # Interactive analysis
│   │   └── ngs_analysis_template.ipynb
│   └── project2/
├── workflow/                   # Snakemake workflow
│   ├── Snakefile
│   ├── envs/                   # Conda environments
│   └── scripts/                # Analysis scripts
└── main.nf                     # Nextflow workflow
```

## Configuration

### Project Configuration (`config/project_config.yml`)

```yaml
project:
  name: "my_experiment"

paths:
  data_dir: "data/my_experiment"
  results_dir: "results/my_experiment"
  notebooks_dir: "notebooks/my_experiment"

samples:
  treatments: ["Treat1", "Treat2"]
  controls: ["Control1", "Control2"]
  batches: ["batch1", "batch2"]
  replicates: [1, 2]

assays: ["rnaseq", "euseq", "sdripseq", "endseq"]

comparisons:
  primary:
    - name: "Treat1_vs_Control1"
      case: "Treat1"
      control: "Control1"
      description: "Treatment 1 vs Control 1"
  
  batch_specific:
    - name: "Treat1_batch1_vs_Control1_batch1"
      case: "Treat1_batch1"
      control: "Control1_batch1"
      description: "Treatment 1 vs Control 1 in batch 1"
```

### Workflow Parameters

**Nextflow:**
```bash
--samples          # Path to samplesheet TSV
--assays           # Comma-separated assay list
--genome_fasta     # hs1 genome FASTA
--genome_gtf       # hs1 annotation GTF
--outdir           # Output directory
```

**Snakemake:**
```bash
--configfile       # Project configuration YAML
--cores            # Number of CPU cores
--use-conda        # Use conda environments
```
    --star_index /path/to/star_index
```

### Assay-Specific Parameters

```bash
# Custom MACS2 parameters for peak calling
nextflow run main.nf \
    --macs2_endseq_params "--nomodel --shift -100 --extsize 200" \
    --macs2_sdripseq_params "--broad --broad-cutoff 0.1"
```

### Resource Allocation

```bash
# Custom resource limits
nextflow run main.nf -profile slurm \
    --max_cpus 48 \
    --max_memory 192.GB \
    --max_time 24.h
```

## Biological Insights

The pipeline automatically generates "Top 5 Biology Insights" by integrating:

1. **Rescue Analysis**: Genes with significant DE in Treated vs Untreated that are rescued by Treated2 vs. Untreated2
2. **R-loop Integration**: sDRIP-seq peaks near transcription start/end sites
3. **Termination Defects**: END-seq breaks concordant with RNA changes
4. **Methylation Changes**: CpG methylation at promoters/gene bodies
5. **Pathway Enrichment**: GO/KEGG terms for RNA processing, R-loops, DNA damage

## Extending the Pipeline

### Adding New Assays

1. Add assay-specific mapping process in `main.nf`
2. Create assay module in `ngs_pipeline/`
3. Update samplesheet validation
4. Add visualization components

Example for ChIP-seq:
```nextflow
process MAP_CHIPSEQ {
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.bam"), path("*.bai")
    
    script:
    """
    bowtie2 -x $genome_index -U $reads | samtools sort > ${meta.id}.bam
    samtools index ${meta.id}.bam
    """
}
```

### Custom Analyses

Add custom analysis scripts to `bin/` directory and call from processes:

```nextflow
process CUSTOM_ANALYSIS {
    script:
    """
    custom_script.py --input $input --output $output
    """
}
```

## Troubleshooting

### Common Issues

1. **Out of Memory**: Increase `--max_memory` or use `--slurm_queue highmem`
2. **STAR Index**: Ensure sufficient RAM (>30GB) for hs1 indexing
3. **Network Issues**: Use `--enable_tcga false` if API access fails
4. **Disk Space**: Pipeline requires ~100GB for full hs1 analysis

### Resource Requirements

- **Minimum**: 16 cores, 64GB RAM, 200GB disk
- **Recommended**: 32 cores, 128GB RAM, 500GB disk
- **STAR Indexing**: 48GB RAM for hs1 genome

### Contact & Support

- Issues: GitHub Issues
- Documentation: See `docs/` directory
- Examples: See `examples/` directory

## Citation

If you use this pipeline, please cite:
- Nextflow: doi.org/10.1038/nbt.3820
- STAR: doi.org/10.1093/bioinformatics/bts635
- DESeq2: doi.org/10.1186/s13059-014-0550-8
- MACS2: doi.org/10.1186/gb-2008-9-9-r137

## License

MIT License - see LICENSE file for details.