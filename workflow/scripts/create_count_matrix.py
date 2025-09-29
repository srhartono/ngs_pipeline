#!/usr/bin/env python3
"""
Create count matrix from individual featureCounts output files
"""

import pandas as pd
import os
from pathlib import Path

def main():
    # Get input files from snakemake
    count_files = snakemake.input.counts
    output_matrix = snakemake.output.matrix
    output_metadata = snakemake.output.metadata
    
    # Parse count files and create matrix
    count_data = {}
    sample_metadata = []
    
    for count_file in count_files:
        # Extract sample name from path
        sample_name = Path(count_file).stem.replace('_counts', '')
        
        # Parse sample information
        parts = sample_name.split('_')
        condition = parts[0]
        batch = parts[1]
        replicate = parts[2].replace('rep', '')
        assay = parts[3]
        
        # Read featureCounts output (skip first line which is comment)
        df = pd.read_csv(count_file, sep='\t', skiprows=1)
        
        # Extract gene counts (last column)
        counts = df.iloc[:, -1]
        gene_ids = df.iloc[:, 0]
        
        count_data[sample_name] = counts
        
        # Add to metadata
        sample_metadata.append({
            'sample_id': sample_name,
            'condition': condition,
            'batch': batch,
            'replicate': int(replicate),
            'assay': assay,
            'treatment_type': 'treatment' if condition.startswith('Treat') else 'control'
        })
    
    # Create count matrix
    count_matrix = pd.DataFrame(count_data, index=gene_ids)
    count_matrix.index.name = 'gene_id'
    
    # Create metadata dataframe
    metadata_df = pd.DataFrame(sample_metadata)
    
    # Save outputs
    count_matrix.to_csv(output_matrix, sep='\t')
    metadata_df.to_csv(output_metadata, sep='\t', index=False)
    
    print(f"Count matrix saved to: {output_matrix}")
    print(f"Sample metadata saved to: {output_metadata}")
    print(f"Matrix shape: {count_matrix.shape}")
    print(f"Samples: {len(metadata_df)}")

if __name__ == "__main__":
    main()