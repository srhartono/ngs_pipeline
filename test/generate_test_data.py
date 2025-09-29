#!/usr/bin/env python3
"""
NGS Pipeline - Test Data Generator

Generates minimal synthetic FASTQ files for testing the pipeline.
Creates realistic-looking reads that will pass through the analysis workflow.
"""

import gzip
import random
import sys
from pathlib import Path
from typing import List, Tuple, Dict
import argparse


class FASTQGenerator:
    """Generate synthetic FASTQ files for testing."""
    
    def __init__(self, read_length: int = 150, seed: int = 42):
        """
        Initialize FASTQ generator.
        
        Args:
            read_length: Length of generated reads
            seed: Random seed for reproducibility
        """
        self.read_length = read_length
        random.seed(seed)
        
        # Define nucleotide probabilities (slightly biased toward A/T for realism)
        self.nucleotides = ['A', 'T', 'G', 'C']
        self.weights = [0.28, 0.28, 0.22, 0.22]
        
        # Quality score characters (Phred+33)
        self.quality_chars = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
        
        # Simulated gene sequences for realistic mapping
        self.gene_sequences = self._generate_gene_sequences()
    
    def _generate_gene_sequences(self) -> Dict[str, str]:
        """Generate some realistic gene sequences."""
        genes = {}
        
        # Add some common housekeeping genes
        housekeeping_genes = ['GAPDH', 'ACTB', 'TUBB', 'RPL13A', 'HPRT1']
        for gene in housekeeping_genes:
            genes[gene] = self._random_sequence(random.randint(1000, 3000))
        
        # Add some random genes
        for i in range(50):
            gene_name = f'GENE{i:03d}'
            genes[gene_name] = self._random_sequence(random.randint(500, 4000))
        
        return genes
    
    def _random_sequence(self, length: int) -> str:
        """Generate random DNA sequence."""
        return ''.join(random.choices(self.nucleotides, weights=self.weights, k=length))
    
    def _generate_quality_scores(self, length: int, mean_quality: int = 35) -> str:
        """Generate realistic quality scores."""
        qualities = []
        for i in range(length):
            # Quality decreases slightly toward the end of reads
            pos_effect = max(0, mean_quality - (i / length) * 5)
            quality = max(2, int(random.normalvariate(pos_effect, 8)))
            quality = min(quality, 41)  # Cap at Phred 41
            qualities.append(self.quality_chars[quality])
        return ''.join(qualities)
    
    def _introduce_errors(self, sequence: str, error_rate: float = 0.01) -> str:
        """Introduce sequencing errors into a sequence."""
        sequence = list(sequence)
        num_errors = int(len(sequence) * error_rate)
        
        for _ in range(num_errors):
            pos = random.randint(0, len(sequence) - 1)
            # 80% substitutions, 10% insertions, 10% deletions
            error_type = random.choices(['sub', 'ins', 'del'], weights=[0.8, 0.1, 0.1])[0]
            
            if error_type == 'sub':
                original = sequence[pos]
                new_base = random.choice([n for n in self.nucleotides if n != original])
                sequence[pos] = new_base
            elif error_type == 'ins':
                sequence.insert(pos, random.choice(self.nucleotides))
            elif error_type == 'del' and len(sequence) > 100:
                sequence.pop(pos)
        
        return ''.join(sequence)
    
    def _simulate_read_from_gene(self, gene_seq: str, strand: str = '+') -> str:
        """Simulate a read from a gene sequence."""
        if len(gene_seq) < self.read_length:
            # For short genes, just use the whole sequence padded
            read_seq = gene_seq + self._random_sequence(self.read_length - len(gene_seq))
        else:
            # Random position within the gene
            start_pos = random.randint(0, len(gene_seq) - self.read_length)
            read_seq = gene_seq[start_pos:start_pos + self.read_length]
        
        # Reverse complement for minus strand
        if strand == '-':
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
            read_seq = ''.join(complement.get(base, base) for base in reversed(read_seq))
        
        return read_seq
    
    def generate_rnaseq_reads(self, num_reads: int) -> List[Tuple[str, str, str]]:
        """Generate RNA-seq reads with realistic gene expression patterns."""
        reads = []
        
        # Expression levels for different genes (log-normal distribution)
        gene_names = list(self.gene_sequences.keys())
        expression_levels = {gene: max(1, random.lognormvariate(2, 1)) for gene in gene_names}
        
        # Normalize to sum to num_reads
        total_expression = sum(expression_levels.values())
        for gene in expression_levels:
            expression_levels[gene] = int((expression_levels[gene] / total_expression) * num_reads)
        
        read_id = 1
        for gene, count in expression_levels.items():
            gene_seq = self.gene_sequences[gene]
            
            for _ in range(count):
                # Random strand
                strand = random.choice(['+', '-'])
                read_seq = self._simulate_read_from_gene(gene_seq, strand)
                read_seq = self._introduce_errors(read_seq, 0.005)  # Low error rate
                
                # Trim to exact read length
                read_seq = read_seq[:self.read_length]
                if len(read_seq) < self.read_length:
                    read_seq += self._random_sequence(self.read_length - len(read_seq))
                
                quality = self._generate_quality_scores(len(read_seq), 37)
                header = f"@RNA_READ_{read_id}_{gene}_{strand}"
                
                reads.append((header, read_seq, quality))
                read_id += 1
                
                if len(reads) >= num_reads:
                    break
            
            if len(reads) >= num_reads:
                break
        
        return reads[:num_reads]
    
    def generate_euseq_reads(self, num_reads: int) -> List[Tuple[str, str, str]]:
        """Generate EU-seq reads (nascent RNA, more biased toward active genes)."""
        reads = []
        
        # EU-seq should be enriched for highly expressed genes
        active_genes = ['GAPDH', 'ACTB', 'TUBB'] + [f'GENE{i:03d}' for i in range(10)]
        
        read_id = 1
        for _ in range(num_reads):
            # Bias toward active genes
            gene = random.choices(active_genes, weights=[3, 2, 2] + [1] * 10)[0]
            gene_seq = self.gene_sequences[gene]
            
            strand = random.choice(['+', '-'])
            read_seq = self._simulate_read_from_gene(gene_seq, strand)
            read_seq = self._introduce_errors(read_seq, 0.008)  # Slightly higher error
            
            read_seq = read_seq[:self.read_length]
            if len(read_seq) < self.read_length:
                read_seq += self._random_sequence(self.read_length - len(read_seq))
            
            quality = self._generate_quality_scores(len(read_seq), 35)
            header = f"@EU_READ_{read_id}_{gene}_{strand}"
            
            reads.append((header, read_seq, quality))
            read_id += 1
        
        return reads
    
    def generate_sdripseq_reads(self, num_reads: int) -> List[Tuple[str, str, str]]:
        """Generate sDRIP-seq reads (enriched for R-loops)."""
        reads = []
        
        # sDRIP-seq enriched regions (simulate R-loop prone sequences)
        # Add G-rich sequences and transcription start sites
        read_id = 1
        for _ in range(num_reads):
            if random.random() < 0.3:
                # G-rich R-loop sequence
                read_seq = self._generate_g_rich_sequence()
            else:
                # Regular genomic sequence
                gene = random.choice(list(self.gene_sequences.keys()))
                gene_seq = self.gene_sequences[gene]
                read_seq = self._simulate_read_from_gene(gene_seq)
            
            read_seq = self._introduce_errors(read_seq, 0.01)
            read_seq = read_seq[:self.read_length]
            if len(read_seq) < self.read_length:
                read_seq += self._random_sequence(self.read_length - len(read_seq))
            
            quality = self._generate_quality_scores(len(read_seq), 33)
            header = f"@SDRIP_READ_{read_id}"
            
            reads.append((header, read_seq, quality))
            read_id += 1
        
        return reads
    
    def _generate_g_rich_sequence(self) -> str:
        """Generate G-rich sequence typical of R-loops."""
        # Create sequence with elevated G content
        nucleotides = ['A', 'T', 'G', 'C']
        weights = [0.15, 0.15, 0.45, 0.25]  # G-rich
        return ''.join(random.choices(nucleotides, weights=weights, k=self.read_length))
    
    def generate_endseq_reads(self, num_reads: int) -> List[Tuple[str, str, str]]:
        """Generate ENDseq reads (5' end mapping with specific patterns)."""
        reads = []
        
        read_id = 1
        for _ in range(num_reads):
            gene = random.choice(list(self.gene_sequences.keys()))
            gene_seq = self.gene_sequences[gene]
            strand = random.choice(['+', '-'])
            
            # ENDseq focuses on 5' ends, so simulate reads from transcript start sites
            if strand == '+':
                # Take reads from 5' end of gene (start of sequence)
                start_pos = random.randint(0, min(500, len(gene_seq) - self.read_length))
            else:
                # Take reads from 3' end of gene (which represents 5' of antisense)
                start_pos = random.randint(max(0, len(gene_seq) - 500), len(gene_seq) - self.read_length)
            
            read_seq = gene_seq[start_pos:start_pos + self.read_length]
            
            # Add some noise
            read_seq = self._introduce_errors(read_seq, error_rate=0.005)
            quality = self._generate_quality_scores(len(read_seq), mean_quality=38)
            
            header = f"@END_READ_{read_id}_{gene}_{strand}"
            reads.append((header, read_seq, quality))
            read_id += 1
        
        return reads
    
    def generate_endseq_reads_legacy(self, num_reads: int) -> List[Tuple[str, str, str]]:
        """Generate END-seq reads (3' end enriched)."""
        reads = []
        
        read_id = 1
        for _ in range(num_reads):
            gene = random.choice(list(self.gene_sequences.keys()))
            gene_seq = self.gene_sequences[gene]
            
            # Bias toward 3' end of genes
            if len(gene_seq) > self.read_length * 2:
                # Take from last 1/3 of gene
                start_region = len(gene_seq) * 2 // 3
                start_pos = random.randint(start_region, len(gene_seq) - self.read_length)
                read_seq = gene_seq[start_pos:start_pos + self.read_length]
            else:
                read_seq = self._simulate_read_from_gene(gene_seq)
            
            read_seq = self._introduce_errors(read_seq, 0.007)
            read_seq = read_seq[:self.read_length]
            if len(read_seq) < self.read_length:
                read_seq += self._random_sequence(self.read_length - len(read_seq))
            
            quality = self._generate_quality_scores(len(read_seq), 36)
            header = f"@END_READ_{read_id}"
            
            reads.append((header, read_seq, quality))
            read_id += 1
        
        return reads
    
    def write_fastq(self, reads: List[Tuple[str, str, str]], output_file: Path, compress: bool = True):
        """
        Write reads to FASTQ file.
        
        Args:
            reads: List of (header, sequence, quality) tuples
            output_file: Output file path
            compress: Whether to gzip compress the output
        """
        output_file = Path(output_file)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        if compress and not str(output_file).endswith('.gz'):
            output_file = Path(str(output_file) + '.gz')
        
        opener = gzip.open if compress else open
        mode = 'wt' if compress else 'w'
        
        with opener(output_file, mode) as f:
            for header, sequence, quality in reads:
                f.write(f"{header}\n")
                f.write(f"{sequence}\n")
                f.write("+\n")
                f.write(f"{quality}\n")
        
        print(f"Generated {len(reads)} reads in {output_file}")


def create_sample_data(output_dir: Path, num_reads: int = 10000):
    """
    Create complete test dataset with all assays and conditions.
    
    Args:
        output_dir: Directory to create test data
        num_reads: Number of reads per FASTQ file
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    generator = FASTQGenerator()
    
    # Define conditions and samples
    # Updated naming scheme: Treat1/Treat2/Control1/Control2 with batch and replicate structure
    treatments = ['Treat1', 'Treat2']
    controls = ['Control1', 'Control2']
    batches = ['batch1', 'batch2']
    replicates = [1, 2]
    assays = ['rnaseq', 'euseq', 'sdripseq', 'endseq']  # Removed bsseq
    
    # Generate all condition combinations
    all_conditions = []
    for treatment in treatments + controls:
        for batch in batches:
            for rep in replicates:
                all_conditions.append(f"{treatment}_{batch}_rep{rep}")
    
    # Generate data for each condition and assay
    for condition in all_conditions:
        for assay in assays:
            sample_id = f"{condition}_{assay}"
            
            # Generate reads based on assay type
            if assay == 'rnaseq':
                reads = generator.generate_rnaseq_reads(num_reads)
            elif assay == 'euseq':
                reads = generator.generate_euseq_reads(num_reads)
            elif assay == 'sdripseq':
                reads = generator.generate_sdripseq_reads(num_reads)
            elif assay == 'endseq':
                reads = generator.generate_endseq_reads(num_reads)
            
            # Write FASTQ file
            fastq_file = output_dir / f"{sample_id}.fastq.gz"
            generator.write_fastq(reads, fastq_file)
    
    print(f"\nTest data generated in: {output_dir}")
    print(f"Total files: {len(all_conditions) * len(assays)}")
    print(f"Conditions: {len(all_conditions)} ({len(treatments)} treatments + {len(controls)} controls)")
    print(f"Assays: {len(assays)} (RNA-seq, EU-seq, sDRIP-seq, ENDseq)")
    print(f"Structure: Treatment/Control_BatchX_repY_assay.fastq.gz")
    print(f"Reads per file: {num_reads}")


def create_samplesheet(output_dir: Path):
    """Create sample sheet for test data."""
    output_dir = Path(output_dir)
    
    # Updated naming scheme
    treatments = ['Treat1', 'Treat2']
    controls = ['Control1', 'Control2']
    batches = ['batch1', 'batch2']
    replicates = [1, 2]
    assays = ['rnaseq', 'euseq', 'sdripseq', 'endseq']
    
    samplesheet_file = output_dir / 'samplesheet.tsv'
    
    with open(samplesheet_file, 'w') as f:
        # Header
        f.write("sample_id\tcondition\tbatch\treplicate\tassay\tfastq_1\tfastq_2\tstrandedness\n")
        
        # Data rows
        for treatment in treatments + controls:
            for batch in batches:
                for rep in replicates:
                    for assay in assays:
                        sample_id = f"{treatment}_{batch}_rep{rep}_{assay}"
                        fastq_file = f"{sample_id}.fastq.gz"
                        
                        # Parse condition type (treatment vs control)
                        condition_type = "treatment" if treatment.startswith("Treat") else "control"
                        
                        f.write(f"{sample_id}\t{treatment}\t{batch}\t{rep}\t{assay}\t{fastq_file}\t\tauto\n")
    
    print(f"Sample sheet created: {samplesheet_file}")


def main():
    """Command-line interface for test data generation."""
    parser = argparse.ArgumentParser(
        description='Generate synthetic test data for NGS Pipeline'
    )
    parser.add_argument(
        '--output-dir', '-o',
        type=Path,
        default=Path.cwd() / 'test_data',
        help='Output directory for test data'
    )
    parser.add_argument(
        '--num-reads', '-n',
        type=int,
        default=10000,
        help='Number of reads per FASTQ file'
    )
    parser.add_argument(
        '--read-length', '-l',
        type=int,
        default=150,
        help='Length of generated reads'
    )
    parser.add_argument(
        '--seed', '-s',
        type=int,
        default=42,
        help='Random seed for reproducibility'
    )
    parser.add_argument(
        '--samplesheet-only',
        action='store_true',
        help='Generate only the sample sheet'
    )
    
    args = parser.parse_args()
    
    if args.samplesheet_only:
        create_samplesheet(args.output_dir)
    else:
        print(f"Generating test data with {args.num_reads} reads per file...")
        print(f"Read length: {args.read_length}")
        print(f"Random seed: {args.seed}")
        print(f"Output directory: {args.output_dir}")
        print()
        
        # Set up generator
        global generator
        generator = FASTQGenerator(args.read_length, args.seed)
        
        # Generate test data
        create_sample_data(args.output_dir, args.num_reads)
        
        # Create sample sheet
        create_samplesheet(args.output_dir)
        
        print("\nTest data generation complete!")
        print(f"To run pipeline with test data:")
        print(f"  nextflow run main.nf --input {args.output_dir}/samplesheet.tsv --outdir results")


if __name__ == '__main__':
    main()