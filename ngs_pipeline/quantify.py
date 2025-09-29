"""
Quantification module for gene expression analysis.

This module provides functions for processing count matrices and 
preparing data for differential expression analysis.
"""

import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
import pandas as pd
import numpy as np
from .utils import validate_file_exists, save_metrics_json, merge_count_files

logger = logging.getLogger(__name__)

def run_quantification(
    count_files: List[Path],
    samplesheet: Path,
    output_dir: Path,
    design_formula: str = "~ batch + condition"
) -> Dict[str, Any]:
    """Run gene quantification analysis."""
    logger.info("Running quantification analysis")
    
    # Merge count files
    merged_counts_file = output_dir / "merged_counts.tsv"
    merged_df = merge_count_files(count_files, merged_counts_file)
    
    return {"merged_counts": str(merged_counts_file), "n_genes": len(merged_df)}