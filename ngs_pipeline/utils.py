"""
Utility functions for the NGS Pipeline.

This module provides common utility functions used across the pipeline,
including logging setup, file validation, and data processing helpers.
"""

import logging
import sys
from pathlib import Path
from typing import Optional, Dict, Any, List, Union
import pandas as pd
import json
import hashlib
import gzip
from rich.logging import RichHandler
from rich.console import Console

console = Console()

def setup_logging(level: int = logging.INFO) -> None:
    """
    Set up logging with Rich handler for colored output.
    
    Args:
        level: Logging level (default: INFO)
    """
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True)]
    )

def validate_file_exists(file_path: Union[str, Path]) -> Path:
    """
    Validate that a file exists and return Path object.
    
    Args:
        file_path: Path to file
        
    Returns:
        Path object if file exists
        
    Raises:
        FileNotFoundError: If file doesn't exist
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    return path

def validate_directory_exists(dir_path: Union[str, Path], create: bool = False) -> Path:
    """
    Validate that a directory exists, optionally create it.
    
    Args:
        dir_path: Path to directory
        create: Whether to create directory if it doesn't exist
        
    Returns:
        Path object
        
    Raises:
        FileNotFoundError: If directory doesn't exist and create=False
    """
    path = Path(dir_path)
    if not path.exists():
        if create:
            path.mkdir(parents=True, exist_ok=True)
        else:
            raise FileNotFoundError(f"Directory not found: {path}")
    return path

def get_file_md5(file_path: Union[str, Path]) -> str:
    """
    Calculate MD5 hash of a file.
    
    Args:
        file_path: Path to file
        
    Returns:
        MD5 hash string
    """
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def is_gzipped(file_path: Union[str, Path]) -> bool:
    """
    Check if a file is gzipped.
    
    Args:
        file_path: Path to file
        
    Returns:
        True if file is gzipped
    """
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def count_fastq_reads(file_path: Union[str, Path]) -> int:
    """
    Count the number of reads in a FASTQ file.
    
    Args:
        file_path: Path to FASTQ file
        
    Returns:
        Number of reads
    """
    if is_gzipped(file_path):
        opener = gzip.open
    else:
        opener = open
        
    with opener(file_path, 'rt') as f:
        lines = sum(1 for _ in f)
    
    return lines // 4

def validate_samplesheet(samplesheet_path: Union[str, Path]) -> pd.DataFrame:
    """
    Validate samplesheet format and content.
    
    Args:
        samplesheet_path: Path to samplesheet TSV file
        
    Returns:
        Validated DataFrame
        
    Raises:
        ValueError: If samplesheet format is invalid
    """
    required_columns = [
        'sample_id', 'condition', 'assay', 'layout', 'fastq_1', 'replicate'
    ]
    
    optional_columns = [
        'fastq_2', 'adapter', 'read_length', 'library_strandedness', 'batch'
    ]
    
    valid_assays = ['RNAseq', 'EUseq', 'sDRIPseq', 'BSseq', 'ENDseq']
    valid_layouts = ['SE', 'PE']
    valid_conditions = ['Untreated', 'siXRN2', 'siXRN2_XRN2WT_OE', 'siXRN2_XRN2MT_OE']
    valid_strandedness = ['auto', 'fr', 'rf', 'unstranded']
    
    # Read samplesheet
    try:
        df = pd.read_csv(samplesheet_path, sep='\t')
    except Exception as e:
        raise ValueError(f"Could not read samplesheet: {e}")
    
    # Check required columns
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Validate values
    invalid_assays = df[~df['assay'].isin(valid_assays)]['assay'].unique()
    if len(invalid_assays) > 0:
        raise ValueError(f"Invalid assay types: {invalid_assays}")
    
    invalid_layouts = df[~df['layout'].isin(valid_layouts)]['layout'].unique()
    if len(invalid_layouts) > 0:
        raise ValueError(f"Invalid layout types: {invalid_layouts}")
    
    invalid_conditions = df[~df['condition'].isin(valid_conditions)]['condition'].unique()
    if len(invalid_conditions) > 0:
        raise ValueError(f"Invalid conditions: {invalid_conditions}")
    
    # Check paired-end samples have fastq_2
    pe_samples = df[df['layout'] == 'PE']
    if len(pe_samples) > 0:
        missing_r2 = pe_samples[pe_samples['fastq_2'].isna()]
        if len(missing_r2) > 0:
            raise ValueError(f"Paired-end samples missing fastq_2: {missing_r2['sample_id'].tolist()}")
    
    # Validate file paths exist
    for _, row in df.iterrows():
        try:
            validate_file_exists(row['fastq_1'])
        except FileNotFoundError:
            raise ValueError(f"FASTQ file not found for sample {row['sample_id']}: {row['fastq_1']}")
        
        if row['layout'] == 'PE' and pd.notna(row['fastq_2']):
            try:
                validate_file_exists(row['fastq_2'])
            except FileNotFoundError:
                raise ValueError(f"FASTQ R2 file not found for sample {row['sample_id']}: {row['fastq_2']}")
    
    # Fill in defaults
    if 'library_strandedness' not in df.columns:
        df['library_strandedness'] = 'auto'
    else:
        df['library_strandedness'] = df['library_strandedness'].fillna('auto')
        invalid_strand = df[~df['library_strandedness'].isin(valid_strandedness)]['library_strandedness'].unique()
        if len(invalid_strand) > 0:
            raise ValueError(f"Invalid strandedness values: {invalid_strand}")
    
    return df

def save_metrics_json(metrics: Dict[str, Any], output_file: Union[str, Path]) -> None:
    """
    Save metrics dictionary to JSON file.
    
    Args:
        metrics: Dictionary of metrics
        output_file: Output JSON file path
    """
    with open(output_file, 'w') as f:
        json.dump(metrics, f, indent=2)

def load_metrics_json(json_file: Union[str, Path]) -> Dict[str, Any]:
    """
    Load metrics from JSON file.
    
    Args:
        json_file: Path to JSON file
        
    Returns:
        Dictionary of metrics
    """
    with open(json_file, 'r') as f:
        return json.load(f)

def merge_count_files(count_files: List[Path], output_file: Path) -> pd.DataFrame:
    """
    Merge multiple count files into a single matrix.
    
    Args:
        count_files: List of count file paths
        output_file: Output merged count file
        
    Returns:
        Merged count DataFrame
    """
    count_dfs = []
    
    for count_file in count_files:
        # Assume featureCounts format with first column as gene_id
        df = pd.read_csv(count_file, sep='\t', comment='#')
        
        # Extract sample name from filename
        sample_name = Path(count_file).stem
        
        # Keep only gene_id and count columns
        gene_col = df.columns[0]  # Usually 'Geneid'
        count_col = df.columns[-1]  # Usually the last column
        
        sample_df = df[[gene_col, count_col]].copy()
        sample_df.columns = ['gene_id', sample_name]
        sample_df.set_index('gene_id', inplace=True)
        
        count_dfs.append(sample_df)
    
    # Merge all DataFrames
    merged_df = pd.concat(count_dfs, axis=1, join='outer')
    merged_df.fillna(0, inplace=True)
    
    # Save merged counts
    merged_df.to_csv(output_file, sep='\t')
    
    return merged_df

def format_number(num: Union[int, float], precision: int = 2) -> str:
    """
    Format number with appropriate precision and units.
    
    Args:
        num: Number to format
        precision: Decimal precision
        
    Returns:
        Formatted number string
    """
    if num >= 1e9:
        return f"{num/1e9:.{precision}f}B"
    elif num >= 1e6:
        return f"{num/1e6:.{precision}f}M"
    elif num >= 1e3:
        return f"{num/1e3:.{precision}f}K"
    else:
        return f"{num:.{precision}f}"

def get_file_size_mb(file_path: Union[str, Path]) -> float:
    """
    Get file size in megabytes.
    
    Args:
        file_path: Path to file
        
    Returns:
        File size in MB
    """
    return Path(file_path).stat().st_size / (1024 * 1024)

class ProgressReporter:
    """Simple progress reporter for long-running operations."""
    
    def __init__(self, total: int, description: str = "Processing"):
        self.total = total
        self.current = 0
        self.description = description
    
    def update(self, increment: int = 1) -> None:
        """Update progress by increment."""
        self.current += increment
        percent = (self.current / self.total) * 100
        console.print(f"{self.description}: {self.current}/{self.total} ({percent:.1f}%)")
    
    def finish(self) -> None:
        """Mark as complete."""
        console.print(f"[bold green]{self.description}: Complete![/bold green]")

def create_output_dirs(base_dir: Path, subdirs: List[str]) -> Dict[str, Path]:
    """
    Create output directory structure.
    
    Args:
        base_dir: Base output directory
        subdirs: List of subdirectory names
        
    Returns:
        Dictionary mapping subdir names to Path objects
    """
    dirs = {}
    
    for subdir in subdirs:
        dir_path = base_dir / subdir
        dir_path.mkdir(parents=True, exist_ok=True)
        dirs[subdir] = dir_path
    
    return dirs