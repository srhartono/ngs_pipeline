"""
Quality Control module for NGS Pipeline.

This module provides functions for analyzing FASTQ quality, adapter content,
contamination screening, and generating QC reports.
"""

import logging
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import json

from .utils import (
    validate_file_exists, validate_directory_exists, 
    save_metrics_json, create_output_dirs, format_number,
    count_fastq_reads, is_gzipped, ProgressReporter
)

logger = logging.getLogger(__name__)

def run_fastqc(
    fastq_files: List[Path], 
    output_dir: Path, 
    threads: int = 4
) -> List[Path]:
    """
    Run FastQC on FASTQ files.
    
    Args:
        fastq_files: List of FASTQ file paths
        output_dir: Output directory for FastQC results
        threads: Number of threads to use
        
    Returns:
        List of FastQC output HTML files
    """
    logger.info(f"Running FastQC on {len(fastq_files)} files")
    
    output_dir = validate_directory_exists(output_dir, create=True)
    html_files = []
    
    for fastq_file in fastq_files:
        validate_file_exists(fastq_file)
        
        cmd = [
            'fastqc', 
            str(fastq_file), 
            '--outdir', str(output_dir),
            '--threads', str(threads),
            '--quiet'
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.debug(f"FastQC completed for {fastq_file}")
            
            # Determine output HTML file
            base_name = fastq_file.stem
            if base_name.endswith('.fastq'):
                base_name = base_name[:-6]
            elif base_name.endswith('.fq'):
                base_name = base_name[:-3]
                
            html_file = output_dir / f"{base_name}_fastqc.html"
            if html_file.exists():
                html_files.append(html_file)
                
        except subprocess.CalledProcessError as e:
            logger.error(f"FastQC failed for {fastq_file}: {e}")
            raise
    
    return html_files

def run_trimming(
    fastq_files: List[Path],
    output_dir: Path,
    adapter_seq: Optional[str] = None,
    quality_cutoff: int = 20,
    min_length: int = 20,
    threads: int = 4
) -> Tuple[List[Path], Dict[str, Any]]:
    """
    Run adapter trimming and quality filtering.
    
    Args:
        fastq_files: List of FASTQ files (1 or 2 for SE/PE)
        output_dir: Output directory
        adapter_seq: Adapter sequence to trim
        quality_cutoff: Quality score cutoff
        min_length: Minimum read length after trimming
        threads: Number of threads
        
    Returns:
        Tuple of (trimmed_files, trimming_stats)
    """
    logger.info("Running adapter trimming and quality filtering")
    
    output_dir = validate_directory_exists(output_dir, create=True)
    
    if len(fastq_files) == 1:
        # Single-end
        return _run_trimming_se(
            fastq_files[0], output_dir, adapter_seq, 
            quality_cutoff, min_length, threads
        )
    elif len(fastq_files) == 2:
        # Paired-end
        return _run_trimming_pe(
            fastq_files[0], fastq_files[1], output_dir, 
            adapter_seq, quality_cutoff, min_length, threads
        )
    else:
        raise ValueError("Expected 1 (SE) or 2 (PE) FASTQ files")

def _run_trimming_se(
    fastq_file: Path, 
    output_dir: Path, 
    adapter_seq: Optional[str],
    quality_cutoff: int, 
    min_length: int, 
    threads: int
) -> Tuple[List[Path], Dict[str, Any]]:
    """Run trimming for single-end reads."""
    
    base_name = fastq_file.stem
    if base_name.endswith('.fastq'):
        base_name = base_name[:-6]
    elif base_name.endswith('.fq'):
        base_name = base_name[:-3]
    
    trimmed_file = output_dir / f"{base_name}_trimmed.fq.gz"
    report_file = output_dir / f"{base_name}_trimming_report.txt"
    
    cmd = [
        'trim_galore',
        '--quality', str(quality_cutoff),
        '--length', str(min_length),
        '--output_dir', str(output_dir),
        '--fastqc',
        '--gzip'
    ]
    
    if adapter_seq:
        cmd.extend(['--adapter', adapter_seq])
    
    cmd.append(str(fastq_file))
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.debug(f"Trimming completed for {fastq_file}")
        
        # Parse trimming statistics
        stats = _parse_trim_galore_report(report_file)
        
        return [trimmed_file], stats
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Trimming failed for {fastq_file}: {e}")
        raise

def _run_trimming_pe(
    fastq_r1: Path, 
    fastq_r2: Path, 
    output_dir: Path,
    adapter_seq: Optional[str], 
    quality_cutoff: int, 
    min_length: int, 
    threads: int
) -> Tuple[List[Path], Dict[str, Any]]:
    """Run trimming for paired-end reads."""
    
    base_name = fastq_r1.stem
    if base_name.endswith('.fastq'):
        base_name = base_name[:-6]
    elif base_name.endswith('.fq'):
        base_name = base_name[:-3]
    
    # Remove R1/R2 suffix if present
    if base_name.endswith('_R1') or base_name.endswith('_1'):
        base_name = base_name[:-3]
    elif base_name.endswith('_R2') or base_name.endswith('_2'):
        base_name = base_name[:-3]
    
    trimmed_r1 = output_dir / f"{base_name}_R1_trimmed.fq.gz"
    trimmed_r2 = output_dir / f"{base_name}_R2_trimmed.fq.gz"
    
    cmd = [
        'trim_galore',
        '--paired',
        '--quality', str(quality_cutoff),
        '--length', str(min_length),
        '--output_dir', str(output_dir),
        '--fastqc',
        '--gzip'
    ]
    
    if adapter_seq:
        cmd.extend(['--adapter', adapter_seq])
    
    cmd.extend([str(fastq_r1), str(fastq_r2)])
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.debug(f"Paired-end trimming completed for {fastq_r1}, {fastq_r2}")
        
        # Parse trimming statistics
        report_r1 = output_dir / f"{fastq_r1.stem}_trimming_report.txt"
        report_r2 = output_dir / f"{fastq_r2.stem}_trimming_report.txt"
        
        stats_r1 = _parse_trim_galore_report(report_r1)
        stats_r2 = _parse_trim_galore_report(report_r2)
        
        # Combine stats
        combined_stats = {
            'R1': stats_r1,
            'R2': stats_r2,
            'total_pairs_processed': stats_r1.get('sequences_processed', 0),
            'pairs_surviving': min(
                stats_r1.get('sequences_surviving', 0),
                stats_r2.get('sequences_surviving', 0)
            )
        }
        
        return [trimmed_r1, trimmed_r2], combined_stats
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Paired-end trimming failed: {e}")
        raise

def _parse_trim_galore_report(report_file: Path) -> Dict[str, Any]:
    """Parse Trim Galore report file for statistics."""
    
    stats = {}
    
    if not report_file.exists():
        logger.warning(f"Trim Galore report not found: {report_file}")
        return stats
    
    try:
        with open(report_file, 'r') as f:
            content = f.read()
            
        # Extract key statistics using string parsing
        lines = content.split('\n')
        for line in lines:
            if 'Total reads processed:' in line:
                stats['sequences_processed'] = int(line.split(':')[1].strip().replace(',', ''))
            elif 'Reads with adapters:' in line:
                stats['sequences_with_adapters'] = int(line.split(':')[1].strip().split()[0].replace(',', ''))
            elif 'Reads written (passing filters):' in line:
                stats['sequences_surviving'] = int(line.split(':')[1].strip().split()[0].replace(',', ''))
                
    except Exception as e:
        logger.warning(f"Could not parse Trim Galore report {report_file}: {e}")
    
    return stats

def check_contamination(
    fastq_files: List[Path],
    output_dir: Path,
    reference_db: Optional[Path] = None,
    threads: int = 4
) -> Dict[str, Any]:
    """
    Check for contamination using FastQ Screen.
    
    Args:
        fastq_files: List of FASTQ files
        output_dir: Output directory
        reference_db: Optional custom reference database
        threads: Number of threads
        
    Returns:
        Contamination screening results
    """
    logger.info("Running contamination screening")
    
    output_dir = validate_directory_exists(output_dir, create=True)
    
    contamination_results = {}
    
    for fastq_file in fastq_files:
        validate_file_exists(fastq_file)
        
        cmd = [
            'fastq_screen',
            '--threads', str(threads),
            '--outdir', str(output_dir),
            '--subset', '100000',  # Screen subset of reads for speed
            str(fastq_file)
        ]
        
        if reference_db:
            cmd.extend(['--conf', str(reference_db)])
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            # Parse results
            base_name = fastq_file.stem
            if base_name.endswith('.fastq'):
                base_name = base_name[:-6]
            elif base_name.endswith('.fq'):
                base_name = base_name[:-3]
            
            screen_file = output_dir / f"{base_name}_screen.txt"
            if screen_file.exists():
                contamination_results[str(fastq_file)] = _parse_fastq_screen_results(screen_file)
                
        except subprocess.CalledProcessError as e:
            logger.warning(f"FastQ Screen failed for {fastq_file}: {e}")
            contamination_results[str(fastq_file)] = {'error': str(e)}
    
    return contamination_results

def _parse_fastq_screen_results(screen_file: Path) -> Dict[str, Any]:
    """Parse FastQ Screen results file."""
    
    results = {}
    
    try:
        with open(screen_file, 'r') as f:
            lines = f.readlines()
        
        # Skip header lines and parse data
        data_lines = [line.strip() for line in lines if not line.startswith('#') and line.strip()]
        
        for line in data_lines[1:]:  # Skip column headers
            parts = line.split('\t')
            if len(parts) >= 5:
                genome = parts[0]
                reads_processed = int(parts[1])
                unmapped = float(parts[2])
                one_hit_one_genome = float(parts[3])
                multiple_hits_one_genome = float(parts[4])
                
                results[genome] = {
                    'reads_processed': reads_processed,
                    'unmapped_percent': unmapped,
                    'one_hit_one_genome_percent': one_hit_one_genome,
                    'multiple_hits_one_genome_percent': multiple_hits_one_genome
                }
                
    except Exception as e:
        logger.warning(f"Could not parse FastQ Screen results {screen_file}: {e}")
    
    return results

def infer_strandedness(
    bam_file: Path,
    bed_file: Path,
    output_dir: Path
) -> Dict[str, Any]:
    """
    Infer library strandedness using RSeQC infer_experiment.py.
    
    Args:
        bam_file: Aligned BAM file
        bed_file: Gene model BED file
        output_dir: Output directory
        
    Returns:
        Strandedness inference results
    """
    logger.info("Inferring library strandedness")
    
    output_dir = validate_directory_exists(output_dir, create=True)
    output_file = output_dir / f"{bam_file.stem}_strandedness.txt"
    
    cmd = [
        'infer_experiment.py',
        '-r', str(bed_file),
        '-i', str(bam_file)
    ]
    
    try:
        with open(output_file, 'w') as out_f:
            result = subprocess.run(cmd, check=True, stdout=out_f, stderr=subprocess.PIPE, text=True)
        
        # Parse results
        strandedness_results = _parse_strandedness_results(output_file)
        return strandedness_results
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Strandedness inference failed for {bam_file}: {e}")
        return {'error': str(e)}

def _parse_strandedness_results(output_file: Path) -> Dict[str, Any]:
    """Parse RSeQC infer_experiment.py output."""
    
    results = {}
    
    try:
        with open(output_file, 'r') as f:
            content = f.read()
        
        lines = content.strip().split('\n')
        
        for line in lines:
            if 'Fraction of reads failed to determine:' in line:
                results['failed_to_determine'] = float(line.split(':')[1].strip())
            elif 'Fraction of reads explained by "++,--":' in line:
                results['forward_strand'] = float(line.split(':')[1].strip())
            elif 'Fraction of reads explained by "+-,-+":' in line:
                results['reverse_strand'] = float(line.split(':')[1].strip())
        
        # Determine strandedness
        if 'forward_strand' in results and 'reverse_strand' in results:
            forward = results['forward_strand']
            reverse = results['reverse_strand']
            
            if forward > 0.8:
                results['library_type'] = 'fr-secondstrand'  # Same as fr
                results['strandedness'] = 'forward'
            elif reverse > 0.8:
                results['library_type'] = 'fr-firststrand'   # Same as rf
                results['strandedness'] = 'reverse'
            else:
                results['library_type'] = 'unstranded'
                results['strandedness'] = 'unstranded'
                
    except Exception as e:
        logger.warning(f"Could not parse strandedness results {output_file}: {e}")
    
    return results

def generate_qc_summary(
    fastqc_results: List[Path],
    trimming_stats: Dict[str, Any],
    contamination_results: Dict[str, Any],
    output_dir: Path
) -> Path:
    """
    Generate QC summary report.
    
    Args:
        fastqc_results: List of FastQC HTML files
        trimming_stats: Trimming statistics
        contamination_results: Contamination screening results
        output_dir: Output directory
        
    Returns:
        Path to summary report file
    """
    logger.info("Generating QC summary report")
    
    output_dir = validate_directory_exists(output_dir, create=True)
    summary_file = output_dir / "qc_summary.json"
    
    summary = {
        'fastqc_files': [str(f) for f in fastqc_results],
        'trimming_statistics': trimming_stats,
        'contamination_screening': contamination_results,
        'total_files_processed': len(fastqc_results)
    }
    
    save_metrics_json(summary, summary_file)
    
    # Also create a human-readable text summary
    text_summary_file = output_dir / "qc_summary.txt"
    with open(text_summary_file, 'w') as f:
        f.write("NGS Pipeline - QC Summary\n")
        f.write("=" * 40 + "\n\n")
        
        f.write(f"Total files processed: {len(fastqc_results)}\n\n")
        
        if trimming_stats:
            f.write("Trimming Statistics:\n")
            f.write("-" * 20 + "\n")
            for key, value in trimming_stats.items():
                f.write(f"{key}: {value}\n")
            f.write("\n")
        
        if contamination_results:
            f.write("Contamination Screening:\n")
            f.write("-" * 25 + "\n")
            for file_path, results in contamination_results.items():
                f.write(f"File: {Path(file_path).name}\n")
                if 'error' in results:
                    f.write(f"  Error: {results['error']}\n")
                else:
                    for genome, stats in results.items():
                        f.write(f"  {genome}: {stats.get('one_hit_one_genome_percent', 0):.1f}% mapped\n")
                f.write("\n")
    
    logger.info(f"QC summary saved to {summary_file}")
    return summary_file

def run_qc_analysis(
    input_dir: Path,
    output_dir: Path,
    assay: str = "RNAseq",
    threads: int = 4,
    adapter_seq: Optional[str] = None,
    quality_cutoff: int = 20,
    min_length: int = 20
) -> Dict[str, Any]:
    """
    Run complete QC analysis pipeline.
    
    Args:
        input_dir: Directory containing FASTQ files
        output_dir: Output directory for QC results
        assay: Assay type
        threads: Number of threads to use
        adapter_seq: Adapter sequence for trimming
        quality_cutoff: Quality score cutoff for trimming
        min_length: Minimum read length after trimming
        
    Returns:
        Dictionary with QC results
    """
    logger.info(f"Starting QC analysis for {assay}")
    
    # Create output subdirectories
    output_dirs = create_output_dirs(output_dir, [
        'fastqc', 'trimming', 'contamination', 'strandedness', 'summary'
    ])
    
    # Find FASTQ files
    fastq_files = []
    for pattern in ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']:
        fastq_files.extend(input_dir.glob(pattern))
    
    if not fastq_files:
        raise ValueError(f"No FASTQ files found in {input_dir}")
    
    logger.info(f"Found {len(fastq_files)} FASTQ files")
    
    # Run FastQC
    fastqc_results = run_fastqc(fastq_files, output_dirs['fastqc'], threads)
    
    # Run trimming on first file as example
    if fastq_files:
        trimmed_files, trimming_stats = run_trimming(
            [fastq_files[0]], 
            output_dirs['trimming'], 
            adapter_seq, 
            quality_cutoff, 
            min_length, 
            threads
        )
    else:
        trimming_stats = {}
    
    # Run contamination screening (optional)
    contamination_results = {}
    try:
        contamination_results = check_contamination(
            fastq_files[:2],  # Screen first 2 files only for speed
            output_dirs['contamination'], 
            threads=threads
        )
    except Exception as e:
        logger.warning(f"Contamination screening failed: {e}")
    
    # Generate summary
    summary_file = generate_qc_summary(
        fastqc_results,
        trimming_stats,
        contamination_results,
        output_dirs['summary']
    )
    
    results = {
        'fastqc_results': fastqc_results,
        'trimming_stats': trimming_stats,
        'contamination_results': contamination_results,
        'summary_file': summary_file,
        'output_directories': output_dirs
    }
    
    logger.info("QC analysis completed successfully")
    return results