"""
Visualization module.
"""

import logging
from pathlib import Path
from typing import Dict, Any

logger = logging.getLogger(__name__)

def create_visualizations(
    data_dir: Path,
    output_dir: Path,
    genome_gtf: Path,
    create_igv: bool = True,
    create_ucsc: bool = True
) -> Dict[str, Any]:
    """Create visualizations."""
    logger.info("Creating visualizations")
    return {"igv": create_igv, "ucsc": create_ucsc}