#!/usr/bin/env python3
"""
NGS Pipeline CLI

Command-line interface for the NGS analysis pipeline.
Supports multi-condition experiments with RNA-seq, EU-seq, ENDseq, and sDRIP-seq.
Provides utilities for QC, analysis, visualization, and reporting.
"""

import typer
import sys
from pathlib import Path
from typing import Optional, List
from rich.console import Console
from rich.logging import RichHandler
import logging

from . import __version__
from .utils import setup_logging, validate_file_exists
from .qc import run_qc_analysis
from .mapping import run_mapping_qc
from .peaks import analyze_peaks
from .quantify import run_quantification
from .variants import analyze_variants
from .methylation import analyze_methylation
from .ml import run_ml_analysis
from .go_string import run_pathway_analysis
from .report import generate_final_report
from .viz import create_visualizations

app = typer.Typer(
    name="ngs_pipeline",
    help="NGS Pipeline - Comprehensive analysis of NGS data with multiple conditions",
    add_completion=False,
)

console = Console()

# Global options
def version_callback(value: bool):
    if value:
        console.print(f"NGS Pipeline v{__version__}")
        raise typer.Exit()

def verbose_callback(value: bool):
    if value:
        setup_logging(level=logging.DEBUG)
    else:
        setup_logging(level=logging.INFO)

@app.callback()
def main(
    version: Optional[bool] = typer.Option(
        None, "--version", "-v", 
        callback=version_callback, 
        is_eager=True,
        help="Show version and exit"
    ),
    verbose: bool = typer.Option(
        False, "--verbose", 
        callback=verbose_callback,
        help="Enable verbose logging"
    ),
):
    """NGS Pipeline CLI"""
    pass

@app.command()
def qc(
    input_dir: Path = typer.Argument(..., help="Input directory with FASTQ files"),
    output_dir: Path = typer.Option("./qc_results", help="Output directory"),
    assay: str = typer.Option("RNAseq", help="Assay type"),
    threads: int = typer.Option(4, help="Number of threads"),
):
    """Run quality control analysis on FASTQ files."""
    console.print(f"[bold blue]Running QC analysis for {assay}[/bold blue]")
    
    try:
        results = run_qc_analysis(
            input_dir=input_dir,
            output_dir=output_dir,
            assay=assay,
            threads=threads
        )
        console.print(f"[bold green]QC analysis completed successfully![/bold green]")
        console.print(f"Results saved to: {output_dir}")
        
    except Exception as e:
        console.print(f"[bold red]Error in QC analysis: {e}[/bold red]")
        sys.exit(1)

@app.command()
def mapping_qc(
    bam_files: List[Path] = typer.Argument(..., help="BAM files to analyze"),
    output_dir: Path = typer.Option("./mapping_qc", help="Output directory"),
    genome_gtf: Optional[Path] = typer.Option(None, help="GTF file for feature annotation"),
):
    """Analyze mapping quality and generate statistics."""
    console.print("[bold blue]Running mapping QC analysis[/bold blue]")
    
    try:
        results = run_mapping_qc(
            bam_files=bam_files,
            output_dir=output_dir,
            genome_gtf=genome_gtf
        )
        console.print("[bold green]Mapping QC completed successfully![/bold green]")
        
    except Exception as e:
        console.print(f"[bold red]Error in mapping QC: {e}[/bold red]")
        sys.exit(1)

@app.command()
def peaks(
    peak_files: List[Path] = typer.Argument(..., help="Peak files (narrowPeak/broadPeak)"),
    output_dir: Path = typer.Option("./peak_analysis", help="Output directory"),
    assay: str = typer.Option("sDRIPseq", help="Assay type"),
    genome_gtf: Optional[Path] = typer.Option(None, help="GTF file for annotation"),
):
    """Analyze peak calling results and generate summaries."""
    console.print(f"[bold blue]Analyzing {assay} peaks[/bold blue]")
    
    try:
        results = analyze_peaks(
            peak_files=peak_files,
            output_dir=output_dir,
            assay=assay,
            genome_gtf=genome_gtf
        )
        console.print("[bold green]Peak analysis completed successfully![/bold green]")
        
    except Exception as e:
        console.print(f"[bold red]Error in peak analysis: {e}[/bold red]")
        sys.exit(1)

@app.command()
def quantify(
    count_files: List[Path] = typer.Argument(..., help="Count matrix files"),
    samplesheet: Path = typer.Argument(..., help="Sample metadata file"),
    output_dir: Path = typer.Option("./quantification", help="Output directory"),
    design_formula: str = typer.Option("~ batch + condition", help="DESeq2 design formula"),
):
    """Run gene quantification and differential expression analysis."""
    console.print("[bold blue]Running quantification analysis[/bold blue]")
    
    try:
        results = run_quantification(
            count_files=count_files,
            samplesheet=samplesheet,
            output_dir=output_dir,
            design_formula=design_formula
        )
        console.print("[bold green]Quantification analysis completed![/bold green]")
        
    except Exception as e:
        console.print(f"[bold red]Error in quantification: {e}[/bold red]")
        sys.exit(1)

@app.command()
def variants(
    vcf_files: List[Path] = typer.Argument(..., help="VCF files to analyze"),
    output_dir: Path = typer.Option("./variant_analysis", help="Output directory"),
    dbsnp_vcf: Optional[Path] = typer.Option(None, help="dbSNP VCF for annotation"),
):
    """Analyze variant calling results."""
    console.print("[bold blue]Running variant analysis[/bold blue]")
    
    try:
        results = analyze_variants(
            vcf_files=vcf_files,
            output_dir=output_dir,
            dbsnp_vcf=dbsnp_vcf
        )
        console.print("[bold green]Variant analysis completed![/bold green]")
        
    except Exception as e:
        console.print(f"[bold red]Error in variant analysis: {e}[/bold red]")
        sys.exit(1)

@app.command()
def methylation(
    methylation_files: List[Path] = typer.Argument(..., help="Methylation context files"),
    output_dir: Path = typer.Option("./methylation_analysis", help="Output directory"),
    genome_gtf: Path = typer.Option(..., help="GTF file for gene annotation"),
):
    """Analyze methylation patterns and generate metagene plots."""
    console.print("[bold blue]Running methylation analysis[/bold blue]")
    
    try:
        results = analyze_methylation(
            methylation_files=methylation_files,
            output_dir=output_dir,
            genome_gtf=genome_gtf
        )
        console.print("[bold green]Methylation analysis completed![/bold green]")
        
    except Exception as e:
        console.print(f"[bold red]Error in methylation analysis: {e}[/bold red]")
        sys.exit(1)

@app.command()
def ml_analysis(
    data_files: List[Path] = typer.Argument(..., help="Data files for ML analysis"),
    samplesheet: Path = typer.Argument(..., help="Sample metadata file"),
    output_dir: Path = typer.Option("./ml_analysis", help="Output directory"),
    methods: List[str] = typer.Option(["pca", "lda", "kmeans"], help="ML methods to run"),
):
    """Run machine learning integration analysis."""
    console.print("[bold blue]Running ML integration analysis[/bold blue]")
    
    try:
        results = run_ml_analysis(
            data_files=data_files,
            samplesheet=samplesheet,
            output_dir=output_dir,
            methods=methods
        )
        console.print("[bold green]ML analysis completed![/bold green]")
        
    except Exception as e:
        console.print(f"[bold red]Error in ML analysis: {e}[/bold red]")
        sys.exit(1)

@app.command()
def pathways(
    de_results: List[Path] = typer.Argument(..., help="Differential expression results"),
    output_dir: Path = typer.Option("./pathway_analysis", help="Output directory"),
    enable_string: bool = typer.Option(True, help="Enable STRING network analysis"),
    enable_tcga: bool = typer.Option(False, help="Enable TCGA cross-reference"),
    string_species: int = typer.Option(9606, help="STRING species ID (9606=human)"),
):
    """Run pathway and network analysis."""
    console.print("[bold blue]Running pathway analysis[/bold blue]")
    
    try:
        results = run_pathway_analysis(
            de_results=de_results,
            output_dir=output_dir,
            enable_string=enable_string,
            enable_tcga=enable_tcga,
            string_species=string_species
        )
        console.print("[bold green]Pathway analysis completed![/bold green]")
        
    except Exception as e:
        console.print(f"[bold red]Error in pathway analysis: {e}[/bold red]")
        sys.exit(1)

@app.command()
def visualize(
    data_dir: Path = typer.Argument(..., help="Directory with analysis results"),
    output_dir: Path = typer.Option("./visualizations", help="Output directory"),
    genome_gtf: Path = typer.Option(..., help="GTF file for annotations"),
    create_igv: bool = typer.Option(True, help="Create IGV session files"),
    create_ucsc: bool = typer.Option(True, help="Create UCSC track hub"),
):
    """Create visualizations and browser sessions."""
    console.print("[bold blue]Creating visualizations[/bold blue]")
    
    try:
        results = create_visualizations(
            data_dir=data_dir,
            output_dir=output_dir,
            genome_gtf=genome_gtf,
            create_igv=create_igv,
            create_ucsc=create_ucsc
        )
        console.print("[bold green]Visualizations created successfully![/bold green]")
        
    except Exception as e:
        console.print(f"[bold red]Error creating visualizations: {e}[/bold red]")
        sys.exit(1)

@app.command()
def report(
    results_dir: Path = typer.Argument(..., help="Directory with all analysis results"),
    output_file: Path = typer.Option("final_report.html", help="Output HTML report file"),
    title: str = typer.Option("NGS Analysis Report", help="Report title"),
    samplesheet: Optional[Path] = typer.Option(None, help="Sample metadata file"),
):
    """Generate final consolidated HTML report."""
    console.print("[bold blue]Generating final report[/bold blue]")
    
    try:
        report_file = generate_final_report(
            results_dir=results_dir,
            output_file=output_file,
            title=title,
            samplesheet=samplesheet
        )
        console.print(f"[bold green]Report generated successfully![/bold green]")
        console.print(f"Report saved to: {report_file}")
        
    except Exception as e:
        console.print(f"[bold red]Error generating report: {e}[/bold red]")
        sys.exit(1)

@app.command()
def validate_samplesheet(
    samplesheet: Path = typer.Argument(..., help="Samplesheet file to validate"),
    output_file: Optional[Path] = typer.Option(None, help="Output validated samplesheet"),
):
    """Validate and optionally reformat samplesheet."""
    console.print("[bold blue]Validating samplesheet[/bold blue]")
    
    try:
        from .utils import validate_samplesheet
        
        valid_samples = validate_samplesheet(samplesheet)
        console.print(f"[bold green]Samplesheet is valid![/bold green]")
        console.print(f"Found {len(valid_samples)} valid samples")
        
        if output_file:
            valid_samples.to_csv(output_file, sep='\t', index=False)
            console.print(f"Validated samplesheet saved to: {output_file}")
            
    except Exception as e:
        console.print(f"[bold red]Samplesheet validation failed: {e}[/bold red]")
        sys.exit(1)

if __name__ == "__main__":
    app()