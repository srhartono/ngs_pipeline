# NGS Pipeline Examples

This directory contains example configurations and datasets for the NGS pipeline.

## Overview

The NGS pipeline supports multi-condition experiments with flexible naming schemes:

- **Conditions**: Multiple treatments (Treat1, Treat2) and controls (Control1, Control2)
- **Batches**: Support for batch effects (batch1, batch2)
- **Replicates**: Multiple biological replicates per condition
- **Assays**: RNA-seq, EU-seq, sDRIP-seq, ENDseq

## Configuration Examples

### `project_config.yml`
Complete project configuration demonstrating:
- **Treatments**: Treat1, Treat2
- **Controls**: Control1, Control2
- **Batches**: batch1, batch2
- **Replicates**: 2 biological replicates per condition
- **Assays**: RNA-seq, EU-seq, sDRIP-seq, ENDseq
- **Comparisons**: Primary, batch-specific, and cross-treatment

### `basic_config.yml`
Minimal configuration for testing:
- Only RNA-seq data
- 1 treatment vs 1 control
- Single batch, 2 replicates each
- **Total samples**: 4

### `multiomics_config.yml`
Full multiomics configuration:
- All 4 assay types
- Multiple treatments and controls
- Batch effect modeling
- Comprehensive comparison matrix

## Parameter Files

### `params_quick.json`
Fast execution parameters for testing:
```json
{
  "skip_multiqc": false,
  "skip_deseq2": false,
  "trimgalore_args": "--length 50",
  "star_index": "references/hs1/indices/star",
  "macs2_genome_size": "hs",
  "max_cpus": 4,
  "max_memory": "8.GB",
  "max_time": "2.h"
}
```

### `params_production.json`
Production-ready parameters:
```json
{
  "skip_multiqc": false,
  "skip_deseq2": false,
  "publish_dir_mode": "copy",
  "star_index": "references/hs1/indices/star",
  "bowtie2_index": "references/hs1/indices/bowtie2/hs1",
  "bismark_index": "references/hs1/indices/bismark",
  "genome_fasta": "references/hs1/hs1.fasta",
  "genome_gtf": "references/hs1/hs1.gtf",
  "dbsnp_vcf": "references/hs1/dbsnp.vcf.gz",
  "macs2_genome_size": "hs",
  "deseq2_fdr": 0.05,
  "deseq2_lfc_threshold": 1.0,
  "max_cpus": 16,
  "max_memory": "64.GB",
  "max_time": "24.h"
}
```

### `params_cluster.json`
HPC cluster optimized parameters:
```json
{
  "skip_multiqc": false,
  "skip_deseq2": false,
  "publish_dir_mode": "copy",
  "max_cpus": 32,
  "max_memory": "128.GB",
  "max_time": "48.h",
  "star_args": "--runThreadN 16 --limitBAMsortRAM 32000000000",
  "macs2_args": "--broad --broad-cutoff 0.1",
  "deseq2_shrinkage_type": "apeglm"
}
```

## Execution Examples

### 1. Quick Test Run
```bash
# Generate test data
python test/generate_test_data.py --output-dir test_data --num-reads 10000

# Run with test parameters
nextflow run main.nf \
  --input examples/samplesheet_minimal.tsv \
  --outdir results_test \
  -params-file examples/params_quick.json \
  -profile test
```

### 2. Full Production Run
```bash
# Download and index reference genome
bash references/get_references.sh --hs1 --build-index all

# Run full pipeline
nextflow run main.nf \
  --input examples/samplesheet.tsv \
  --outdir results_production \
  -params-file examples/params_production.json \
  -profile standard
```

### 3. HPC Cluster Execution
```bash
# SLURM cluster
nextflow run main.nf \
  --input examples/samplesheet.tsv \
  --outdir results_cluster \
  -params-file examples/params_cluster.json \
  -profile slurm

# PBS/Torque cluster
nextflow run main.nf \
  --input examples/samplesheet.tsv \
  --outdir results_cluster \
  -params-file examples/params_cluster.json \
  -profile pbs
```

### 4. Cloud Execution (AWS Batch)
```bash
# Set up AWS credentials and configuration
export AWS_DEFAULT_REGION=us-east-1

# Run on AWS Batch
nextflow run main.nf \
  --input examples/samplesheet.tsv \
  --outdir s3://my-bucket/results \
  -params-file examples/params_production.json \
  -profile aws_batch \
  -w s3://my-bucket/work
```

### 5. Docker/Singularity Execution
```bash
# Docker
nextflow run main.nf \
  --input examples/samplesheet.tsv \
  --outdir results_docker \
  -profile docker

# Singularity
nextflow run main.nf \
  --input examples/samplesheet.tsv \
  --outdir results_singularity \
  -profile singularity
```

### 6. Specific Analysis Subsets
```bash
# RNA-seq only (skip other assays)
nextflow run main.nf \
  --input examples/samplesheet_rnaseq_only.tsv \
  --outdir results_rnaseq \
  --skip_sdripseq \
  --skip_bsseq \
  --skip_endseq

# Skip differential expression
nextflow run main.nf \
  --input examples/samplesheet.tsv \
  --outdir results_no_deseq2 \
  --skip_deseq2

# QC and mapping only
nextflow run main.nf \
  --input examples/samplesheet.tsv \
  --outdir results_qc_mapping \
  --skip_peaks \
  --skip_variants \
  --skip_deseq2
```

## Resource Scaling Guidelines

### Small Dataset (< 1M reads per sample)
- CPUs: 2-4
- Memory: 4-8 GB
- Time: 1-2 hours
- Use: `params_quick.json`

### Medium Dataset (1-10M reads per sample)
- CPUs: 8-16
- Memory: 16-32 GB
- Time: 4-8 hours
- Use: `params_production.json`

### Large Dataset (> 10M reads per sample)
- CPUs: 16-32
- Memory: 32-128 GB
- Time: 8-24 hours
- Use: `params_cluster.json`

## Customization Tips

### 1. Modify Trimming Parameters
```json
{
  "trimgalore_args": "--quality 30 --length 100 --stringency 3"
}
```

### 2. Adjust Mapping Sensitivity
```json
{
  "star_args": "--outFilterMultimapNmax 1 --outFilterMismatchNmax 2",
  "bowtie2_args": "--very-sensitive --no-mixed --no-discordant"
}
```

### 3. Fine-tune Peak Calling
```json
{
  "macs2_args": "--broad --broad-cutoff 0.05 --keep-dup auto",
  "macs2_pvalue": "1e-5"
}
```

### 4. Customize Differential Expression
```json
{
  "deseq2_fdr": 0.01,
  "deseq2_lfc_threshold": 1.5,
  "deseq2_shrinkage_type": "ashr"
}
```

### 5. Add Custom Contrasts
Edit the DESeq2 script to include additional comparisons:
```r
# Custom contrasts
contrasts <- list(
  # Add your custom contrast here
  c("condition", "your_condition", "your_reference")
)
```

## Troubleshooting Common Issues

### 1. Memory Issues
- Reduce `max_memory` parameter
- Use `--star_args "--limitBAMsortRAM 10000000000"`
- Split large samples into smaller chunks

### 2. Time Limits
- Increase `max_time` parameter
- Use faster alignment parameters
- Reduce dataset size for testing

### 3. Reference Genome Issues
- Ensure indices are built correctly
- Check file paths in parameter files
- Verify genome version compatibility

### 4. Sample Sheet Validation
```bash
# Validate sample sheet format
python -c "
import pandas as pd
df = pd.read_csv('examples/samplesheet.tsv', sep='\t')
print('Columns:', df.columns.tolist())
print('Samples:', len(df))
print('Assays:', df['assay'].unique())
print('Conditions:', df['condition'].unique())
"
```

For more detailed documentation, see the main [README.md](../README.md) file.