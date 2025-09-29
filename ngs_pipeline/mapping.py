"""
Mapping quality control module for NGS Pipeline.

This module provides functions for analyzing mapping quality, calculating
alignment statistics, and generating mapping-related QC metrics.
"""

import logging
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from .utils import (
    validate_file_exists, validate_directory_exists, 
    save_metrics_json, create_output_dirs, format_number
)

logger = logging.getLogger(__name__)

def get_bam_stats(bam_file: Path, output_dir: Path) -> Dict[str, Any]:
    """
    Get basic BAM file statistics using samtools stats.
    
    Args:
        bam_file: Path to BAM file
        output_dir: Output directory for stats files
        
    Returns:
        Dictionary with alignment statistics
    """
    logger.info(f"Getting BAM statistics for {bam_file}")
    
    validate_file_exists(bam_file)
    output_dir = validate_directory_exists(output_dir, create=True)
    
    stats_file = output_dir / f"{bam_file.stem}_stats.txt"
    
    # Run samtools stats
    cmd = ['samtools', 'stats', str(bam_file)]
    
    try:
        with open(stats_file, 'w') as f:
            result = subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
        
        # Parse stats file
        stats = _parse_samtools_stats(stats_file)
        return stats
        
    except subprocess.CalledProcessError as e:
        logger.error(f"samtools stats failed for {bam_file}: {e}")
        raise

def _parse_samtools_stats(stats_file: Path) -> Dict[str, Any]:
    """Parse samtools stats output file."""
    
    stats = {}
    
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                if line.startswith('SN\t'):
                    # Summary statistics
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        key = parts[1].rstrip(':')
                        value = parts[2]
                        
                        # Convert to appropriate type
                        try:
                            if '.' in value:
                                stats[key] = float(value)
                            else:
                                stats[key] = int(value)
                        except ValueError:
                            stats[key] = value
                            
    except Exception as e:
        logger.warning(f"Could not parse samtools stats {stats_file}: {e}")
    
    return stats

def get_bam_flagstat(bam_file: Path, output_dir: Path) -> Dict[str, Any]:
    """
    Get BAM flagstat statistics.
    
    Args:
        bam_file: Path to BAM file
        output_dir: Output directory
        
    Returns:
        Dictionary with flagstat results
    """
    logger.info(f"Getting flagstat for {bam_file}")
    
    validate_file_exists(bam_file)
    output_dir = validate_directory_exists(output_dir, create=True)
    
    flagstat_file = output_dir / f"{bam_file.stem}_flagstat.txt"
    
    # Run samtools flagstat
    cmd = ['samtools', 'flagstat', str(bam_file)]
    
    try:
        with open(flagstat_file, 'w') as f:
            result = subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
        
        # Parse flagstat file
        flagstat = _parse_samtools_flagstat(flagstat_file)
        return flagstat
        
    except subprocess.CalledProcessError as e:
        logger.error(f"samtools flagstat failed for {bam_file}: {e}")
        raise

def _parse_samtools_flagstat(flagstat_file: Path) -> Dict[str, Any]:
    """Parse samtools flagstat output."""
    
    flagstat = {}
    
    try:
        with open(flagstat_file, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            line = line.strip()
            if 'in total' in line:
                flagstat['total_reads'] = int(line.split()[0])
            elif 'secondary' in line:
                flagstat['secondary_alignments'] = int(line.split()[0])
            elif 'supplementary' in line:
                flagstat['supplementary_alignments'] = int(line.split()[0])
            elif 'duplicates' in line:
                flagstat['duplicate_reads'] = int(line.split()[0])
            elif 'mapped (' in line:
                mapped = int(line.split()[0])
                percent = float(line.split('(')[1].split('%')[0])
                flagstat['mapped_reads'] = mapped
                flagstat['mapping_rate'] = percent
            elif 'paired in sequencing' in line:
                flagstat['paired_reads'] = int(line.split()[0])
            elif 'properly paired' in line:
                paired = int(line.split()[0])
                percent = float(line.split('(')[1].split('%')[0])
                flagstat['properly_paired'] = paired
                flagstat['properly_paired_rate'] = percent
                
    except Exception as e:
        logger.warning(f"Could not parse flagstat {flagstat_file}: {e}")
    
    return flagstat

def get_insert_size_metrics(bam_file: Path, output_dir: Path) -> Dict[str, Any]:
    """
    Calculate insert size metrics for paired-end data using Picard.
    
    Args:
        bam_file: Path to BAM file
        output_dir: Output directory
        
    Returns:
        Dictionary with insert size metrics
    """
    logger.info(f"Calculating insert size metrics for {bam_file}")
    
    validate_file_exists(bam_file)
    output_dir = validate_directory_exists(output_dir, create=True)
    
    metrics_file = output_dir / f"{bam_file.stem}_insert_size_metrics.txt"
    histogram_file = output_dir / f"{bam_file.stem}_insert_size_histogram.pdf"
    
    cmd = [
        'picard', 'CollectInsertSizeMetrics',
        f'INPUT={bam_file}',
        f'OUTPUT={metrics_file}',
        f'HISTOGRAM_FILE={histogram_file}',
        'VALIDATION_STRINGENCY=LENIENT'
    ]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Parse metrics file
        metrics = _parse_picard_insert_size_metrics(metrics_file)
        return metrics
        
    except subprocess.CalledProcessError as e:
        logger.warning(f"Picard CollectInsertSizeMetrics failed for {bam_file}: {e}")
        return {}

def _parse_picard_insert_size_metrics(metrics_file: Path) -> Dict[str, Any]:
    """Parse Picard insert size metrics file."""
    
    metrics = {}
    
    try:
        with open(metrics_file, 'r') as f:
            lines = f.readlines()
        
        # Find metrics section
        metrics_start = False
        for i, line in enumerate(lines):
            if line.startswith('## METRICS CLASS'):
                metrics_start = True
                continue
            elif metrics_start and line.startswith('## HISTOGRAM'):
                break
            elif metrics_start and not line.startswith('#'):
                if 'MEDIAN_INSERT_SIZE' in lines[i-1]:
                    # Header line
                    headers = lines[i-1].strip().split('\t')
                    values = line.strip().split('\t')
                    
                    for header, value in zip(headers, values):
                        try:
                            if '.' in value:
                                metrics[header] = float(value)
                            else:
                                metrics[header] = int(value)
                        except ValueError:
                            metrics[header] = value
                    break
                    
    except Exception as e:
        logger.warning(f"Could not parse insert size metrics {metrics_file}: {e}")
    
    return metrics

def calculate_coverage_stats(
    bam_file: Path, 
    bed_file: Optional[Path],
    output_dir: Path
) -> Dict[str, Any]:
    """
    Calculate coverage statistics using bedtools.
    
    Args:
        bam_file: Path to BAM file
        bed_file: Optional BED file for regions of interest
        output_dir: Output directory
        
    Returns:
        Dictionary with coverage statistics
    """
    logger.info(f"Calculating coverage statistics for {bam_file}")
    
    validate_file_exists(bam_file)
    output_dir = validate_directory_exists(output_dir, create=True)
    
    coverage_file = output_dir / f"{bam_file.stem}_coverage.txt"
    
    if bed_file:
        validate_file_exists(bed_file)
        cmd = ['bedtools', 'coverage', '-a', str(bed_file), '-b', str(bam_file)]
    else:
        cmd = ['bedtools', 'genomecov', '-ibam', str(bam_file), '-bg']
    
    try:
        with open(coverage_file, 'w') as f:
            result = subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
        
        # Calculate basic coverage statistics
        coverage_stats = _calculate_coverage_summary(coverage_file, bed_file is not None)
        return coverage_stats
        
    except subprocess.CalledProcessError as e:
        logger.warning(f"Coverage calculation failed for {bam_file}: {e}")
        return {}

def _calculate_coverage_summary(coverage_file: Path, is_region_based: bool) -> Dict[str, Any]:
    """Calculate summary statistics from coverage file."""
    
    stats = {}
    
    try:
        df = pd.read_csv(coverage_file, sep='\t', header=None)
        
        if is_region_based:
            # bedtools coverage format
            if len(df.columns) >= 7:
                coverage_col = df.iloc[:, -1]  # Last column is usually coverage
                stats['mean_coverage'] = coverage_col.mean()
                stats['median_coverage'] = coverage_col.median()
                stats['std_coverage'] = coverage_col.std()
                stats['regions_with_coverage'] = (coverage_col > 0).sum()
                stats['total_regions'] = len(coverage_col)
        else:
            # bedtools genomecov format
            if len(df.columns) >= 4:
                coverage_col = df.iloc[:, 3]  # 4th column is coverage
                stats['mean_coverage'] = coverage_col.mean()
                stats['median_coverage'] = coverage_col.median()
                stats['std_coverage'] = coverage_col.std()
                
    except Exception as e:
        logger.warning(f"Could not calculate coverage summary from {coverage_file}: {e}")
    
    return stats

def analyze_mapping_quality(
    bam_file: Path,
    output_dir: Path,
    min_mapq: int = 10
) -> Dict[str, Any]:
    """
    Analyze mapping quality distribution.
    
    Args:
        bam_file: Path to BAM file
        output_dir: Output directory
        min_mapq: Minimum mapping quality threshold
        
    Returns:
        Dictionary with mapping quality statistics
    """
    logger.info(f"Analyzing mapping quality for {bam_file}")
    
    validate_file_exists(bam_file)
    output_dir = validate_directory_exists(output_dir, create=True)
    
    # Extract mapping qualities using samtools
    cmd = ['samtools', 'view', str(bam_file)]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        mapq_values = []
        for line in result.stdout.split('\n'):
            if line.strip():
                fields = line.split('\t')
                if len(fields) > 4:
                    try:
                        mapq = int(fields[4])
                        mapq_values.append(mapq)
                    except (ValueError, IndexError):
                        continue
        
        if mapq_values:
            mapq_array = np.array(mapq_values)
            
            stats = {
                'total_reads': len(mapq_values),
                'mean_mapq': float(np.mean(mapq_array)),
                'median_mapq': float(np.median(mapq_array)),
                'std_mapq': float(np.std(mapq_array)),
                'reads_above_threshold': int(np.sum(mapq_array >= min_mapq)),
                'percent_above_threshold': float(np.sum(mapq_array >= min_mapq) / len(mapq_values) * 100)
            }
            
            # Create mapping quality histogram
            _plot_mapq_histogram(mapq_array, output_dir / f"{bam_file.stem}_mapq_histogram.png")
            
            return stats
        else:
            return {'error': 'No mapping quality values found'}
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Mapping quality analysis failed for {bam_file}: {e}")
        return {'error': str(e)}

def _plot_mapq_histogram(mapq_values: np.ndarray, output_file: Path) -> None:
    """Plot mapping quality histogram."""
    
    plt.figure(figsize=(10, 6))
    plt.hist(mapq_values, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('Mapping Quality (MAPQ)')
    plt.ylabel('Number of Reads')
    plt.title('Mapping Quality Distribution')
    plt.grid(True, alpha=0.3)
    
    # Add statistics text
    stats_text = f'Mean: {np.mean(mapq_values):.2f}\nMedian: {np.median(mapq_values):.2f}'
    plt.text(0.7, 0.8, stats_text, transform=plt.gca().transAxes, 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def generate_mapping_qc_report(
    bam_files: List[Path],
    output_dir: Path,
    genome_gtf: Optional[Path] = None
) -> Dict[str, Any]:
    """
    Generate comprehensive mapping QC report.
    
    Args:
        bam_files: List of BAM files to analyze
        output_dir: Output directory
        genome_gtf: Optional GTF file for coverage analysis
        
    Returns:
        Dictionary with all mapping QC results
    """
    logger.info(f"Generating mapping QC report for {len(bam_files)} BAM files")
    
    output_dir = validate_directory_exists(output_dir, create=True)
    
    # Create subdirectories
    output_dirs = create_output_dirs(output_dir, [
        'stats', 'flagstat', 'insert_size', 'coverage', 'mapq', 'summary'
    ])
    
    all_results = {}
    
    for bam_file in bam_files:
        sample_name = bam_file.stem
        logger.info(f"Processing {sample_name}")
        
        sample_results = {
            'bam_file': str(bam_file),
            'sample_name': sample_name
        }
        
        # Basic BAM statistics
        try:
            bam_stats = get_bam_stats(bam_file, output_dirs['stats'])
            sample_results['bam_stats'] = bam_stats
        except Exception as e:
            logger.warning(f"Failed to get BAM stats for {sample_name}: {e}")
            sample_results['bam_stats'] = {'error': str(e)}
        
        # Flagstat
        try:
            flagstat = get_bam_flagstat(bam_file, output_dirs['flagstat'])
            sample_results['flagstat'] = flagstat
        except Exception as e:
            logger.warning(f"Failed to get flagstat for {sample_name}: {e}")
            sample_results['flagstat'] = {'error': str(e)}
        
        # Insert size metrics (for paired-end data)
        try:
            insert_metrics = get_insert_size_metrics(bam_file, output_dirs['insert_size'])
            sample_results['insert_size_metrics'] = insert_metrics
        except Exception as e:
            logger.warning(f"Failed to get insert size metrics for {sample_name}: {e}")
            sample_results['insert_size_metrics'] = {'error': str(e)}
        
        # Coverage statistics
        try:
            coverage_stats = calculate_coverage_stats(bam_file, genome_gtf, output_dirs['coverage'])
            sample_results['coverage_stats'] = coverage_stats
        except Exception as e:
            logger.warning(f"Failed to calculate coverage for {sample_name}: {e}")
            sample_results['coverage_stats'] = {'error': str(e)}
        
        # Mapping quality analysis
        try:
            mapq_stats = analyze_mapping_quality(bam_file, output_dirs['mapq'])
            sample_results['mapq_stats'] = mapq_stats
        except Exception as e:
            logger.warning(f"Failed to analyze mapping quality for {sample_name}: {e}")
            sample_results['mapq_stats'] = {'error': str(e)}
        
        all_results[sample_name] = sample_results
    
    # Create summary report
    summary_file = output_dirs['summary'] / 'mapping_qc_summary.json'
    save_metrics_json(all_results, summary_file)
    
    # Create summary plot
    _create_mapping_summary_plots(all_results, output_dirs['summary'])
    
    logger.info(f"Mapping QC report saved to {summary_file}")
    return all_results

def _create_mapping_summary_plots(results: Dict[str, Any], output_dir: Path) -> None:
    """Create summary plots for mapping QC."""
    
    # Extract mapping rates for all samples
    sample_names = []
    mapping_rates = []
    
    for sample_name, sample_data in results.items():
        flagstat = sample_data.get('flagstat', {})
        if 'mapping_rate' in flagstat and 'error' not in flagstat:
            sample_names.append(sample_name)
            mapping_rates.append(flagstat['mapping_rate'])
    
    if mapping_rates:
        # Mapping rate bar plot
        plt.figure(figsize=(12, 6))
        bars = plt.bar(range(len(sample_names)), mapping_rates)
        plt.xlabel('Samples')
        plt.ylabel('Mapping Rate (%)')
        plt.title('Mapping Rates Across Samples')
        plt.xticks(range(len(sample_names)), sample_names, rotation=45, ha='right')
        plt.ylim(0, 100)
        
        # Color bars based on mapping rate
        for bar, rate in zip(bars, mapping_rates):
            if rate >= 90:
                bar.set_color('green')
            elif rate >= 70:
                bar.set_color('orange')
            else:
                bar.set_color('red')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'mapping_rates_summary.png', dpi=300, bbox_inches='tight')
        plt.close()

def run_mapping_qc(
    bam_files: List[Path],
    output_dir: Path,
    genome_gtf: Optional[Path] = None
) -> Dict[str, Any]:
    """
    Main function to run mapping QC analysis.
    
    Args:
        bam_files: List of BAM files to analyze
        output_dir: Output directory
        genome_gtf: Optional GTF file for coverage analysis
        
    Returns:
        Dictionary with mapping QC results
    """
    logger.info("Starting mapping QC analysis")
    
    results = generate_mapping_qc_report(bam_files, output_dir, genome_gtf)
    
    logger.info("Mapping QC analysis completed successfully")
    return results