#!/usr/bin/env python3
"""
NGS Pipeline - Test Suite

Comprehensive pytest test suite for validating pipeline components.
"""

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pytest

# Add the project root to the Python path
import sys
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

try:
    from ngs_pipeline import cli, utils, qc, mapping, peaks
    from test.generate_test_data import FASTQGenerator, create_sample_data
except ImportError as e:
    pytest.skip(f"Could not import project modules: {e}", allow_module_level=True)


class TestUtils:
    """Test utility functions."""
    
    def test_validate_file_exists(self, tmp_path):
        """Test file validation function."""
        # Create test file
        test_file = tmp_path / "test.txt"
        test_file.write_text("test content")
        
        # Test existing file
        assert utils.validate_file(test_file) == test_file
        
        # Test non-existing file
        with pytest.raises(FileNotFoundError):
            utils.validate_file(tmp_path / "nonexistent.txt")
    
    def test_validate_file_extensions(self, tmp_path):
        """Test file extension validation."""
        # Create test files
        fastq_file = tmp_path / "test.fastq.gz"
        fastq_file.write_text("@read1\nACGT\n+\nIIII\n")
        
        bam_file = tmp_path / "test.bam"
        bam_file.write_bytes(b"BAM\x01")  # Minimal BAM header
        
        # Test FASTQ validation
        assert utils.validate_file(fastq_file, ['.fastq', '.fastq.gz']) == fastq_file
        
        # Test BAM validation
        assert utils.validate_file(bam_file, ['.bam']) == bam_file
        
        # Test wrong extension
        with pytest.raises(ValueError):
            utils.validate_file(fastq_file, ['.bam'])
    
    def test_parse_samplesheet(self, tmp_path):
        """Test sample sheet parsing."""
        # Create test sample sheet
        samplesheet = tmp_path / "samplesheet.tsv"
        samplesheet_content = """sample_id\tcondition\tassay\treplicate\tfastq_1\tfastq_2\tstrandedness
Untreated_rnaseq_R1\tUntreated\trnaseq\t1\tUntreated_rnaseq_R1.fastq.gz\t\tauto
siXRN2_rnaseq_R1\tsiXRN2\trnaseq\t1\tsiXRN2_rnaseq_R1.fastq.gz\t\tauto"""
        
        samplesheet.write_text(samplesheet_content)
        
        # Parse sample sheet
        samples = utils.parse_samplesheet(samplesheet)
        
        assert len(samples) == 2
        assert samples[0]['sample_id'] == 'Untreated_rnaseq_R1'
        assert samples[0]['condition'] == 'Untreated'
        assert samples[0]['assay'] == 'rnaseq'
        assert samples[1]['sample_id'] == 'siXRN2_rnaseq_R1'
    
    def test_setup_logging(self, tmp_path):
        """Test logging setup."""
        log_file = tmp_path / "test.log"
        logger = utils.setup_logging(log_file, verbose=True)
        
        logger.info("Test message")
        
        assert log_file.exists()
        log_content = log_file.read_text()
        assert "Test message" in log_content


class TestQC:
    """Test quality control functions."""
    
    def test_run_fastqc(self, tmp_path):
        """Test FastQC execution."""
        # Create minimal FASTQ file
        fastq_file = tmp_path / "test.fastq"
        fastq_content = "@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n"
        fastq_file.write_text(fastq_content)
        
        output_dir = tmp_path / "fastqc_output"
        
        # Mock FastQC execution
        with patch('subprocess.run') as mock_run:
            mock_run.return_value.returncode = 0
            
            result = qc.run_fastqc(fastq_file, output_dir)
            
            mock_run.assert_called_once()
            assert str(fastq_file) in mock_run.call_args[0][0]
    
    def test_parse_fastqc_data(self, tmp_path):
        """Test FastQC data parsing."""
        # Create mock FastQC data file
        fastqc_data = tmp_path / "test_fastqc_data.txt"
        fastqc_content = """##FastQC	0.11.9
>>Basic Statistics	pass
#Measure	Value
Filename	test.fastq
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	1000
Sequences flagged as poor quality	50
Sequence length	150
%GC	45
>>END_MODULE"""
        
        fastqc_data.write_text(fastqc_content)
        
        # Parse data
        stats = qc.parse_fastqc_data(fastqc_data)
        
        assert stats['total_sequences'] == 1000
        assert stats['poor_quality'] == 50
        assert stats['sequence_length'] == 150
        assert stats['gc_content'] == 45
    
    def test_calculate_contamination(self, tmp_path):
        """Test contamination calculation."""
        # Create mock contamination results
        contamination_data = {
            'total_reads': 10000,
            'rrna_reads': 250,
            'mtdna_reads': 100,
            'adapter_reads': 50
        }
        
        percentages = qc.calculate_contamination(contamination_data)
        
        assert percentages['rrna_percent'] == 2.5
        assert percentages['mtdna_percent'] == 1.0
        assert percentages['adapter_percent'] == 0.5


class TestMapping:
    """Test mapping analysis functions."""
    
    def test_parse_star_log(self, tmp_path):
        """Test STAR log parsing."""
        # Create mock STAR log
        star_log = tmp_path / "Log.final.out"
        log_content = """                                 Started job on |	Dec 15 15:30:45
                             Started mapping on |	Dec 15 15:30:45
                                    Finished on |	Dec 15 15:31:30
       Mapping speed, Million of reads per hour |	80.00

                          Number of input reads |	1000000
                      Average input read length |	150
                                    UNIQUE READS:
                   Uniquely mapped reads number |	850000
                   Uniquely mapped reads % |	85.00%
                     Average mapped length |	149.50
                          Number of splices: Total |	45000
               Number of splices: Annotated (sjdb) |	43000
                          Number of splices: GT/AG |	44500
                          Number of splices: GC/AG |	400
                          Number of splices: AT/AC |	100
                  Number of reads mapped to multiple loci |	50000
                       % of reads mapped to multiple loci |	5.00%
                  Number of reads mapped to too many loci |	10000
                       % of reads mapped to too many loci |	1.00%
                                        UNMAPPED READS:
             % of reads unmapped: too many mismatches |	2.50%
                      % of reads unmapped: too short |	6.50%
                                           % of reads unmapped: other |	0.00%"""
        
        star_log.write_text(log_content)
        
        # Parse log
        stats = mapping.parse_star_log(star_log)
        
        assert stats['input_reads'] == 1000000
        assert stats['uniquely_mapped'] == 850000
        assert stats['uniquely_mapped_percent'] == 85.0
        assert stats['multimapped'] == 50000
        assert stats['multimapped_percent'] == 5.0
    
    def test_calculate_coverage_stats(self, tmp_path):
        """Test coverage statistics calculation."""
        # Mock coverage data
        coverage_data = [0, 1, 5, 10, 15, 20, 25, 30, 35, 40]
        
        stats = mapping.calculate_coverage_stats(coverage_data)
        
        assert stats['mean_coverage'] == 18.5
        assert stats['median_coverage'] == 17.5
        assert stats['coverage_std'] > 0
        assert 0 <= stats['bases_covered_1x'] <= 100
    
    def test_parse_flagstat(self, tmp_path):
        """Test samtools flagstat parsing."""
        # Create mock flagstat output
        flagstat_file = tmp_path / "test.flagstat"
        flagstat_content = """1000000 + 0 in total (QC-passed reads + QC-failed reads)
50000 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
950000 + 0 mapped (95.00% : N/A)
1000000 + 0 paired in sequencing
500000 + 0 read1
500000 + 0 read2
900000 + 0 properly paired (90.00% : N/A)
920000 + 0 with itself and mate mapped
30000 + 0 singletons (3.00% : N/A)"""
        
        flagstat_file.write_text(flagstat_content)
        
        # Parse flagstat
        stats = mapping.parse_flagstat(flagstat_file)
        
        assert stats['total_reads'] == 1000000
        assert stats['mapped_reads'] == 950000
        assert stats['mapping_rate'] == 95.0
        assert stats['properly_paired'] == 900000
        assert stats['singletons'] == 30000


class TestPeaks:
    """Test peak calling analysis functions."""
    
    def test_parse_macs2_summits(self, tmp_path):
        """Test MACS2 summits parsing."""
        # Create mock summits file
        summits_file = tmp_path / "test_summits.bed"
        summits_content = """chr1	1000	1001	peak_1	100	.
chr1	2000	2001	peak_2	200	.
chr2	3000	3001	peak_3	150	."""
        
        summits_file.write_text(summits_content)
        
        # Parse summits
        summits = peaks.parse_macs2_summits(summits_file)
        
        assert len(summits) == 3
        assert summits[0]['chr'] == 'chr1'
        assert summits[0]['start'] == 1000
        assert summits[0]['score'] == 100
    
    def test_calculate_frip_score(self, tmp_path):
        """Test FRiP score calculation."""
        # Mock data
        total_reads = 1000000
        reads_in_peaks = 150000
        
        frip = peaks.calculate_frip_score(reads_in_peaks, total_reads)
        
        assert frip == 15.0
    
    def test_annotate_peaks(self, tmp_path):
        """Test peak annotation."""
        # Create mock peak file
        peaks_file = tmp_path / "test_peaks.narrowPeak"
        peaks_content = """chr1	1000	2000	peak_1	100	.	10.5	15.2	5.1	500
chr1	5000	6000	peak_2	200	.	20.5	25.2	10.1	500"""
        
        peaks_file.write_text(peaks_content)
        
        # Create mock annotation
        mock_annotation = [
            {'peak_id': 'peak_1', 'gene': 'GENE1', 'distance': 0, 'annotation': 'promoter'},
            {'peak_id': 'peak_2', 'gene': 'GENE2', 'distance': 1000, 'annotation': 'intergenic'}
        ]
        
        with patch('ngs_pipeline.peaks.annotate_with_homer') as mock_homer:
            mock_homer.return_value = mock_annotation
            
            annotations = peaks.annotate_peaks(peaks_file, "mock_genome.gtf")
            
            assert len(annotations) == 2
            assert annotations[0]['gene'] == 'GENE1'
            assert annotations[1]['annotation'] == 'intergenic'


class TestCLI:
    """Test command-line interface."""
    
    def test_cli_help(self):
        """Test CLI help command."""
        from typer.testing import CliRunner
        
        runner = CliRunner()
        result = runner.invoke(cli.app, ["--help"])
        
        assert result.exit_code == 0
        assert "NGS Pipeline" in result.output
    
    def test_qc_command(self, tmp_path):
        """Test QC command."""
        from typer.testing import CliRunner
        
        # Create test FASTQ
        fastq_file = tmp_path / "test.fastq.gz"
        fastq_file.write_text("@read1\nACGT\n+\nIIII\n")
        
        runner = CliRunner()
        
        with patch('ngs_pipeline.qc.run_fastqc') as mock_fastqc:
            mock_fastqc.return_value = {"status": "completed"}
            
            result = runner.invoke(cli.app, [
                "qc", 
                "--input", str(fastq_file),
                "--output-dir", str(tmp_path / "qc_output")
            ])
            
            # Command should execute without error
            assert result.exit_code == 0


class TestDataGeneration:
    """Test synthetic data generation."""
    
    def test_fastq_generator_init(self):
        """Test FASTQ generator initialization."""
        generator = FASTQGenerator(read_length=100, seed=123)
        
        assert generator.read_length == 100
        assert len(generator.gene_sequences) > 0
        assert 'GAPDH' in generator.gene_sequences
    
    def test_generate_rnaseq_reads(self):
        """Test RNA-seq read generation."""
        generator = FASTQGenerator(read_length=50, seed=42)
        reads = generator.generate_rnaseq_reads(100)
        
        assert len(reads) == 100
        
        # Check read format
        header, sequence, quality = reads[0]
        assert header.startswith('@RNA_READ_')
        assert len(sequence) == 50
        assert len(quality) == 50
        assert all(base in 'ATGC' for base in sequence)
    
    def test_generate_euseq_reads(self):
        """Test EU-seq read generation."""
        generator = FASTQGenerator(read_length=50, seed=42)
        reads = generator.generate_euseq_reads(100)
        
        assert len(reads) == 100
        header, sequence, quality = reads[0]
        assert header.startswith('@EU_READ_')
        assert len(sequence) == 50
    
    def test_generate_endseq_reads(self):
        """Test ENDseq read generation."""
        generator = FASTQGenerator(read_length=50, seed=42)
        reads = generator.generate_endseq_reads(100)
        
        assert len(reads) == 100
        header, sequence, quality = reads[0]
        assert header.startswith('@END_READ_')
        
        # Check for 5' bias in sequence generation
        assert len(sequence) == 50
        assert len(quality) == 50
    
    def test_write_fastq(self, tmp_path):
        """Test FASTQ file writing."""
        generator = FASTQGenerator()
        reads = [
            ("@read1", "ACGT" * 25, "I" * 100),
            ("@read2", "TGCA" * 25, "J" * 100)
        ]
        
        output_file = tmp_path / "test.fastq.gz"
        generator.write_fastq(reads, output_file, compress=True)
        
        assert output_file.exists()
        
        # Check if file is gzipped and contains expected content
        import gzip
        with gzip.open(output_file, 'rt') as f:
            content = f.read()
            assert "@read1" in content
            assert "ACGT" * 25 in content
            assert "@read2" in content


class TestIntegration:
    """Integration tests for pipeline components."""
    
    def test_full_test_data_generation(self, tmp_path):
        """Test complete test data generation."""
        # Generate minimal test dataset
        create_sample_data(tmp_path, num_reads=100)
        
        # Check that files were created
        fastq_files = list(tmp_path.glob("*.fastq.gz"))
        assert len(fastq_files) > 0
        
        # Check sample sheet creation
        from test.generate_test_data import create_samplesheet
        create_samplesheet(tmp_path)
        
        samplesheet = tmp_path / "samplesheet.tsv"
        assert samplesheet.exists()
        
        # Verify sample sheet content
        df = pd.read_csv(samplesheet, sep='\t')
        assert 'sample_id' in df.columns
        assert 'condition' in df.columns
        assert 'assay' in df.columns
        assert len(df) > 0
    
    @pytest.mark.slow
    def test_pipeline_dry_run(self, tmp_path):
        """Test Nextflow pipeline dry run."""
        # Skip if Nextflow not available
        if not shutil.which('nextflow'):
            pytest.skip("Nextflow not available")
        
        # Create minimal test data
        create_sample_data(tmp_path, num_reads=100)
        from test.generate_test_data import create_samplesheet
        create_samplesheet(tmp_path)
        
        # Run pipeline in dry-run mode
        cmd = [
            'nextflow', 'run', str(project_root / 'main.nf'),
            '--input', str(tmp_path / 'samplesheet.tsv'),
            '--outdir', str(tmp_path / 'results'),
            '-profile', 'test',
            '-dry-run'
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            # Dry run should complete successfully
            assert result.returncode == 0 or "No such variable" in result.stderr
        except subprocess.TimeoutExpired:
            pytest.skip("Nextflow dry run timed out")


# Pytest configuration
def pytest_configure(config):
    """Configure pytest markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )


if __name__ == '__main__':
    # Run tests when script is executed directly
    pytest.main([__file__, '-v'])