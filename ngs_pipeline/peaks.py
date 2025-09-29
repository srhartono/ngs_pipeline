"""
Peak calling analysis module for NGS Pipeline.

This module provides functions for analyzing peak calling results,
calculating peak statistics, and generating peak-related visualizations.
"""

import logging
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict

from .utils import (
    validate_file_exists, validate_directory_exists, 
    save_metrics_json, create_output_dirs, format_number
)

logger = logging.getLogger(__name__)

def parse_peak_file(peak_file: Path, file_format: str = 'auto') -> pd.DataFrame:
    """
    Parse peak file into DataFrame.
    
    Args:
        peak_file: Path to peak file
        file_format: File format ('narrowPeak', 'broadPeak', 'bed', 'auto')
        
    Returns:
        DataFrame with peak information
    """
    validate_file_exists(peak_file)
    
    if file_format == 'auto':
        if peak_file.suffix in ['.narrowPeak', '.broadPeak']:
            file_format = peak_file.suffix[1:]
        else:
            file_format = 'bed'
    
    if file_format == 'narrowPeak':
        columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 
                  'signal_value', 'pvalue', 'qvalue', 'peak_summit']
    elif file_format == 'broadPeak':
        columns = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                  'signal_value', 'pvalue', 'qvalue']
    else:
        # Generic BED format
        columns = ['chrom', 'start', 'end']
    
    try:
        df = pd.read_csv(peak_file, sep='\t', header=None, comment='#')
        
        # Assign column names based on available columns
        df.columns = columns[:len(df.columns)]
        
        # Ensure basic columns are present
        if 'start' in df.columns and 'end' in df.columns:
            df['length'] = df['end'] - df['start']
        
        return df
        
    except Exception as e:
        logger.error(f"Failed to parse peak file {peak_file}: {e}")
        raise

def calculate_peak_statistics(peaks_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Calculate basic peak statistics.
    
    Args:
        peaks_df: DataFrame with peak information
        
    Returns:
        Dictionary with peak statistics
    """
    stats = {
        'total_peaks': len(peaks_df),
        'total_bp_covered': int(peaks_df['length'].sum()) if 'length' in peaks_df.columns else 0
    }
    
    if 'length' in peaks_df.columns:
        stats.update({
            'mean_peak_length': float(peaks_df['length'].mean()),
            'median_peak_length': float(peaks_df['length'].median()),
            'std_peak_length': float(peaks_df['length'].std()),
            'min_peak_length': int(peaks_df['length'].min()),
            'max_peak_length': int(peaks_df['length'].max())
        })
    
    if 'score' in peaks_df.columns:
        stats.update({
            'mean_peak_score': float(peaks_df['score'].mean()),
            'median_peak_score': float(peaks_df['score'].median()),
            'std_peak_score': float(peaks_df['score'].std())
        })
    
    if 'signal_value' in peaks_df.columns:
        stats.update({
            'mean_signal_value': float(peaks_df['signal_value'].mean()),
            'median_signal_value': float(peaks_df['signal_value'].median()),
            'std_signal_value': float(peaks_df['signal_value'].std())
        })
    
    if 'qvalue' in peaks_df.columns:
        # Count peaks by significance thresholds
        stats.update({
            'peaks_qval_0.01': int((peaks_df['qvalue'] >= -np.log10(0.01)).sum()),
            'peaks_qval_0.05': int((peaks_df['qvalue'] >= -np.log10(0.05)).sum()),
            'peaks_qval_0.1': int((peaks_df['qvalue'] >= -np.log10(0.1)).sum())
        })
    
    # Chromosome distribution
    if 'chrom' in peaks_df.columns:
        chrom_counts = peaks_df['chrom'].value_counts().to_dict()
        stats['peaks_per_chromosome'] = chrom_counts
    
    return stats

def calculate_frip_score(
    bam_file: Path,
    peak_file: Path,
    output_dir: Path
) -> Dict[str, Any]:
    """
    Calculate Fraction of Reads in Peaks (FRiP) score.
    
    Args:
        bam_file: Path to BAM file
        peak_file: Path to peak file
        output_dir: Output directory
        
    Returns:
        Dictionary with FRiP statistics
    """
    logger.info(f"Calculating FRiP score for {bam_file} and {peak_file}")
    
    validate_file_exists(bam_file)
    validate_file_exists(peak_file)
    output_dir = validate_directory_exists(output_dir, create=True)
    
    # Count total reads
    cmd_total = ['samtools', 'view', '-c', str(bam_file)]
    
    try:
        result_total = subprocess.run(cmd_total, check=True, capture_output=True, text=True)
        total_reads = int(result_total.stdout.strip())
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to count total reads: {e}")
        return {'error': str(e)}
    
    # Count reads in peaks using bedtools intersect
    cmd_intersect = [
        'bedtools', 'intersect', 
        '-a', str(bam_file), 
        '-b', str(peak_file),
        '-c'
    ]
    
    try:
        result_intersect = subprocess.run(cmd_intersect, check=True, capture_output=True, text=True)
        
        # Count reads that intersect with peaks
        reads_in_peaks = 0
        for line in result_intersect.stdout.strip().split('\n'):
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) > 0:
                    try:
                        count = int(parts[-1])  # Last column should be overlap count
                        if count > 0:
                            reads_in_peaks += 1
                    except (ValueError, IndexError):
                        continue
        
        frip_score = reads_in_peaks / total_reads if total_reads > 0 else 0
        
        frip_stats = {
            'total_reads': total_reads,
            'reads_in_peaks': reads_in_peaks,
            'frip_score': frip_score,
            'frip_percentage': frip_score * 100
        }
        
        # Save FRiP results
        frip_file = output_dir / f"{bam_file.stem}_{peak_file.stem}_frip.json"
        save_metrics_json(frip_stats, frip_file)
        
        return frip_stats
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to calculate FRiP score: {e}")
        return {'error': str(e)}

def annotate_peaks_with_genes(
    peak_file: Path,
    gtf_file: Path,
    output_dir: Path,
    promoter_upstream: int = 2000,
    promoter_downstream: int = 500
) -> pd.DataFrame:
    """
    Annotate peaks with nearby genes.
    
    Args:
        peak_file: Path to peak file
        gtf_file: Path to GTF file
        output_dir: Output directory
        promoter_upstream: Distance upstream of TSS for promoter region
        promoter_downstream: Distance downstream of TSS for promoter region
        
    Returns:
        DataFrame with annotated peaks
    """
    logger.info(f"Annotating peaks with genes")
    
    validate_file_exists(peak_file)
    validate_file_exists(gtf_file)
    output_dir = validate_directory_exists(output_dir, create=True)
    
    # Use bedtools closest to find nearest genes
    annotated_file = output_dir / f"{peak_file.stem}_annotated.txt"
    
    cmd = [
        'bedtools', 'closest',
        '-a', str(peak_file),
        '-b', str(gtf_file),
        '-d'  # Add distance to output
    ]
    
    try:
        with open(annotated_file, 'w') as f:
            result = subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
        
        # Parse annotated file
        annotated_df = pd.read_csv(annotated_file, sep='\t', header=None)
        
        # Add column names based on input format
        peak_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
        gtf_cols = ['gene_chrom', 'source', 'feature', 'gene_start', 'gene_end', 
                   'gene_score', 'gene_strand', 'frame', 'attributes']
        
        n_peak_cols = min(len(peak_cols), annotated_df.shape[1] - len(gtf_cols) - 1)
        n_gtf_cols = min(len(gtf_cols), annotated_df.shape[1] - n_peak_cols - 1)
        
        column_names = peak_cols[:n_peak_cols] + gtf_cols[:n_gtf_cols] + ['distance']
        annotated_df.columns = column_names[:annotated_df.shape[1]]
        
        # Extract gene names from attributes if present
        if 'attributes' in annotated_df.columns:
            annotated_df['gene_name'] = annotated_df['attributes'].str.extract(r'gene_name "([^"]+)"')
            annotated_df['gene_id'] = annotated_df['attributes'].str.extract(r'gene_id "([^"]+)"')
        
        # Classify peak locations
        if 'distance' in annotated_df.columns:
            conditions = [
                annotated_df['distance'] == 0,
                annotated_df['distance'] <= promoter_upstream,
                annotated_df['distance'] <= 10000,
                annotated_df['distance'] <= 50000
            ]
            choices = ['Overlapping', 'Promoter', 'Near gene (<10kb)', 'Distal (<50kb)']
            annotated_df['peak_annotation'] = np.select(conditions, choices, default='Intergenic')
        
        return annotated_df
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Peak annotation failed: {e}")
        return pd.DataFrame()

def compare_peak_sets(
    peak_files: List[Path],
    sample_names: List[str],
    output_dir: Path
) -> Dict[str, Any]:
    """
    Compare multiple peak sets and find overlaps.
    
    Args:
        peak_files: List of peak files
        sample_names: List of sample names
        output_dir: Output directory
        
    Returns:
        Dictionary with comparison results
    """
    logger.info(f"Comparing {len(peak_files)} peak sets")
    
    output_dir = validate_directory_exists(output_dir, create=True)
    
    if len(peak_files) != len(sample_names):
        raise ValueError("Number of peak files must match number of sample names")
    
    # Load all peak sets
    peak_sets = {}
    for peak_file, sample_name in zip(peak_files, sample_names):
        peaks_df = parse_peak_file(peak_file)
        peak_sets[sample_name] = peaks_df
    
    # Calculate pairwise overlaps using bedtools
    overlap_matrix = pd.DataFrame(
        index=sample_names, 
        columns=sample_names, 
        dtype=float
    )
    
    for i, sample1 in enumerate(sample_names):
        for j, sample2 in enumerate(sample_names):
            if i == j:
                overlap_matrix.loc[sample1, sample2] = 1.0
            else:
                overlap_pct = _calculate_peak_overlap_percentage(
                    peak_files[i], peak_files[j], output_dir
                )
                overlap_matrix.loc[sample1, sample2] = overlap_pct
    
    # Create overlap heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(overlap_matrix.astype(float), annot=True, cmap='Blues', 
                square=True, cbar_kws={'label': 'Overlap Percentage'})
    plt.title('Peak Set Overlaps')
    plt.tight_layout()
    plt.savefig(output_dir / 'peak_overlap_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Find consensus peaks (present in multiple samples)
    consensus_results = _find_consensus_peaks(peak_files, sample_names, output_dir)
    
    comparison_results = {
        'overlap_matrix': overlap_matrix.to_dict(),
        'consensus_peaks': consensus_results,
        'individual_peak_counts': {name: len(peak_sets[name]) for name in sample_names}
    }
    
    # Save results
    results_file = output_dir / 'peak_comparison_results.json'
    save_metrics_json(comparison_results, results_file)
    
    return comparison_results

def _calculate_peak_overlap_percentage(
    peak_file1: Path, 
    peak_file2: Path, 
    output_dir: Path
) -> float:
    """Calculate percentage overlap between two peak sets."""
    
    temp_overlap_file = output_dir / f"temp_overlap_{peak_file1.stem}_{peak_file2.stem}.txt"
    
    cmd = [
        'bedtools', 'intersect',
        '-a', str(peak_file1),
        '-b', str(peak_file2),
        '-wo'  # Write original entry plus overlap length
    ]
    
    try:
        with open(temp_overlap_file, 'w') as f:
            result = subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
        
        # Count overlapping peaks
        total_peaks1 = sum(1 for _ in open(peak_file1))
        overlapping_peaks = sum(1 for _ in open(temp_overlap_file))
        
        overlap_percentage = (overlapping_peaks / total_peaks1) * 100 if total_peaks1 > 0 else 0
        
        # Clean up temp file
        temp_overlap_file.unlink()
        
        return overlap_percentage
        
    except subprocess.CalledProcessError as e:
        logger.warning(f"Failed to calculate overlap between {peak_file1} and {peak_file2}: {e}")
        return 0.0

def _find_consensus_peaks(
    peak_files: List[Path], 
    sample_names: List[str], 
    output_dir: Path
) -> Dict[str, Any]:
    """Find consensus peaks present in multiple samples."""
    
    if len(peak_files) < 2:
        return {'message': 'Need at least 2 samples for consensus analysis'}
    
    consensus_results = {}
    
    # Merge all peaks and find those present in multiple samples
    merged_file = output_dir / 'all_peaks_merged.bed'
    
    # First, concatenate all peak files
    temp_concat_file = output_dir / 'temp_all_peaks.bed'
    
    with open(temp_concat_file, 'w') as outf:
        for i, peak_file in enumerate(peak_files):
            with open(peak_file, 'r') as inf:
                for line in inf:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            # Add sample info to peak name
                            if len(parts) > 3:
                                parts[3] = f"{sample_names[i]}_{parts[3]}"
                            else:
                                parts.append(f"{sample_names[i]}_peak_{i}")
                            outf.write('\t'.join(parts) + '\n')
    
    # Merge overlapping peaks
    cmd_merge = ['bedtools', 'merge', '-i', str(temp_concat_file), '-c', '4', '-o', 'collapse']
    
    try:
        with open(merged_file, 'w') as f:
            result = subprocess.run(cmd_merge, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
        
        # Analyze merged peaks to find consensus
        consensus_peaks = []
        with open(merged_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    sample_list = parts[3].split(',')
                    unique_samples = set()
                    for sample_peak in sample_list:
                        sample_name = sample_peak.split('_')[0]
                        unique_samples.add(sample_name)
                    
                    if len(unique_samples) >= 2:  # Present in at least 2 samples
                        consensus_peaks.append({
                            'chrom': parts[0],
                            'start': int(parts[1]),
                            'end': int(parts[2]),
                            'samples': list(unique_samples),
                            'n_samples': len(unique_samples)
                        })
        
        consensus_results = {
            'total_consensus_peaks': len(consensus_peaks),
            'consensus_peaks': consensus_peaks[:100],  # Limit output size
            'peaks_in_all_samples': len([p for p in consensus_peaks if p['n_samples'] == len(sample_names)])
        }
        
        # Clean up temp files
        temp_concat_file.unlink()
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to find consensus peaks: {e}")
        consensus_results = {'error': str(e)}
    
    return consensus_results

def create_peak_visualizations(
    peak_files: List[Path],
    sample_names: List[str],
    output_dir: Path
) -> None:
    """Create various peak-related visualizations."""
    
    logger.info("Creating peak visualizations")
    
    output_dir = validate_directory_exists(output_dir, create=True)
    
    all_stats = []
    
    for peak_file, sample_name in zip(peak_files, sample_names):
        peaks_df = parse_peak_file(peak_file)
        stats = calculate_peak_statistics(peaks_df)
        stats['sample'] = sample_name
        all_stats.append(stats)
    
    stats_df = pd.DataFrame(all_stats)
    
    # Peak count comparison
    if 'total_peaks' in stats_df.columns:
        plt.figure(figsize=(10, 6))
        bars = plt.bar(stats_df['sample'], stats_df['total_peaks'])
        plt.xlabel('Sample')
        plt.ylabel('Number of Peaks')
        plt.title('Peak Counts Across Samples')
        plt.xticks(rotation=45, ha='right')
        
        # Color bars by peak count
        max_peaks = stats_df['total_peaks'].max()
        for bar, count in zip(bars, stats_df['total_peaks']):
            color_intensity = count / max_peaks
            bar.set_color(plt.cm.Blues(color_intensity))
        
        plt.tight_layout()
        plt.savefig(output_dir / 'peak_counts_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # Peak length distributions
    if 'mean_peak_length' in stats_df.columns:
        plt.figure(figsize=(10, 6))
        plt.bar(stats_df['sample'], stats_df['mean_peak_length'])
        plt.xlabel('Sample')
        plt.ylabel('Mean Peak Length (bp)')
        plt.title('Mean Peak Lengths Across Samples')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(output_dir / 'peak_lengths_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()

def analyze_peaks(
    peak_files: List[Path],
    output_dir: Path,
    assay: str = "sDRIPseq",
    genome_gtf: Optional[Path] = None
) -> Dict[str, Any]:
    """
    Main function to analyze peak calling results.
    
    Args:
        peak_files: List of peak files to analyze
        output_dir: Output directory
        assay: Assay type
        genome_gtf: Optional GTF file for gene annotation
        
    Returns:
        Dictionary with peak analysis results
    """
    logger.info(f"Starting peak analysis for {assay}")
    
    output_dir = validate_directory_exists(output_dir, create=True)
    
    # Create subdirectories
    output_dirs = create_output_dirs(output_dir, [
        'statistics', 'annotation', 'comparison', 'frip', 'visualizations'
    ])
    
    sample_names = [f.stem for f in peak_files]
    
    # Calculate statistics for each peak file
    all_stats = {}
    for peak_file in peak_files:
        sample_name = peak_file.stem
        peaks_df = parse_peak_file(peak_file)
        stats = calculate_peak_statistics(peaks_df)
        all_stats[sample_name] = stats
        
        # Save individual statistics
        stats_file = output_dirs['statistics'] / f"{sample_name}_peak_stats.json"
        save_metrics_json(stats, stats_file)
    
    # Gene annotation if GTF provided
    annotated_results = {}
    if genome_gtf:
        for peak_file in peak_files:
            sample_name = peak_file.stem
            try:
                annotated_df = annotate_peaks_with_genes(
                    peak_file, genome_gtf, output_dirs['annotation']
                )
                annotated_file = output_dirs['annotation'] / f"{sample_name}_annotated.tsv"
                annotated_df.to_csv(annotated_file, sep='\t', index=False)
                annotated_results[sample_name] = str(annotated_file)
            except Exception as e:
                logger.warning(f"Peak annotation failed for {sample_name}: {e}")
    
    # Compare peak sets if multiple files
    comparison_results = {}
    if len(peak_files) > 1:
        try:
            comparison_results = compare_peak_sets(
                peak_files, sample_names, output_dirs['comparison']
            )
        except Exception as e:
            logger.warning(f"Peak comparison failed: {e}")
    
    # Create visualizations
    try:
        create_peak_visualizations(peak_files, sample_names, output_dirs['visualizations'])
    except Exception as e:
        logger.warning(f"Peak visualization creation failed: {e}")
    
    # Compile final results
    results = {
        'peak_statistics': all_stats,
        'annotated_peaks': annotated_results,
        'comparison_results': comparison_results,
        'output_directories': {k: str(v) for k, v in output_dirs.items()}
    }
    
    # Save final results
    final_results_file = output_dir / 'peak_analysis_results.json'
    save_metrics_json(results, final_results_file)
    
    logger.info("Peak analysis completed successfully")
    return results