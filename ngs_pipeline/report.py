"""
Reporting module.
"""

import logging
from pathlib import Path
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)

def generate_final_report(
    results_dir: Path,
    output_file: Path,
    title: str = "NGS Analysis Report",
    samplesheet: Optional[Path] = None
) -> Path:
    """Generate final HTML report."""
    logger.info("Generating final report")
    
    # Create basic HTML content
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head><title>{title}</title></head>
    <body>
    <h1>{title}</h1>
    <p>Analysis completed successfully.</p>
    <p>Results directory: {results_dir}</p>
    </body>
    </html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    return output_file