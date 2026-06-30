#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comprehensive visualization for both MoC and AMI similarity matrices.
Enhanced version with detailed plots for both metrics.
Filtered for cell types: ASC, OGC, OPC
"""

import numpy as np
import pandas as pd
import os
from pathlib import Path
import sklearn.manifold
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import warnings
warnings.filterwarnings('ignore')

# Input and output directories
input_dir = "/Users/mozhganoroujlu/Desktop/MOZHGUN/cell_fate/hi_c/codes_figures/folders/normalized_contacts/bandnorm/TADs_non_neuron/"
output_dir = "/Users/mozhganoroujlu/Desktop/MOZHGUN/cell_fate/hi_c/codes_figures/folders/normalized_contacts/similarity_results"

# Create output directory
Path(output_dir).mkdir(parents=True, exist_ok=True)

# Load cell type mapping file
celltype_file = "/Users/mozhganoroujlu/Desktop/MOZHGUN/cell_fate/hi_c/codes_figures/folders/non_neuron.tsv"
print("Loading cell type mapping...")

# Check if cell type file exists
if not os.path.exists(celltype_file):
    print(f"ERROR: Cell type file not found at: {celltype_file}")
    exit(1)

# Load the cell type file
celltype_df = pd.read_csv(celltype_file, sep='\t', header=0, names=['barcode', 'celltype'])
print(f"Loaded {len(celltype_df)} total cells from cell type file")

# Filter for desired cell types
desired_celltypes = ['ASC', 'OGC', 'OPC']
filtered_cells = celltype_df[celltype_df['celltype'].isin(desired_celltypes)]
print(f"Found {len(filtered_cells)} cells with cell types {desired_celltypes}")

# Check if similarity matrix files exist
moc_file = os.path.join(output_dir, "simi_moc.xlsx")
ami_file = os.path.join(output_dir, "simi_ami.xlsx")

if not os.path.exists(moc_file):
    print(f"ERROR: MoC similarity matrix file not found at: {moc_file}")
    print("Please make sure the similarity matrices have been generated first.")
    exit(1)

if not os.path.exists(ami_file):
    print(f"ERROR: AMI similarity matrix file not found at: {ami_file}")
    print("Please make sure the similarity matrices have been generated first.")
    exit(1)

# Load similarity matrices
print("Loading similarity matrices...")
moc_df = pd.read_excel(moc_file, index_col=0)
ami_df = pd.read_excel(ami_file, index_col=0)

# Filter barcodes to only include those in filtered_cells
# First, clean the barcodes from the matrices (remove '_tads' suffix if present)
matrix_barcodes_clean = [barcode.replace('_tads', '') for barcode in moc_df.index.tolist()]

# Create a mapping to find which indices to keep
keep_indices = []
keep_barcodes = []
missing_barcodes = []

for i, clean_barcode in enumerate(matrix_barcodes_clean):
    if clean_barcode in filtered_cells['barcode'].values:
        keep_indices.append(i)
        keep_barcodes.append(moc_df.index.tolist()[i])  # Keep original barcode with _tads
    else:
        missing_barcodes.append(clean_barcode)

print(f"Keeping {len(keep_indices)} cells after filtering by cell type")
if len(missing_barcodes) > 0:
    print(f"Note: {len(missing_barcodes)} barcodes from similarity matrices not found in cell type file")
    print(f"First few missing barcodes: {missing_barcodes[:5]}")

if len(keep_indices) == 0:
    print("ERROR: No matching cells found after filtering. Please check barcode formats.")
    print("Sample barcodes from cell type file:", filtered_cells['barcode'].head(10).tolist())
    print("Sample barcodes from similarity matrix:", matrix_barcodes_clean[:10])
    exit(1)

# Filter the matrices
simi_moc = moc_df.iloc[keep_indices, keep_indices].values
simi_ami = ami_df.iloc[keep_indices, keep_indices].values
barcodes = keep_barcodes  # Keep original barcodes with _tads suffix

n_cells = simi_moc.shape[0]
print(f"Filtered matrices for {n_cells} cells")

def create_comprehensive_figure(similarity_matrix, matrix_name, color_scheme):
    """
    Create comprehensive figure for a similarity matrix
    """
    fig = plt.figure(figsize=(15, 10))
    
    # Flattened matrix (off-diagonal elements)
    matrix_flat = similarity_matrix[np.triu_indices_from(similarity_matrix, k=1)]
    
    # 1. Similarity Matrix Heatmap
    plt.subplot(2, 3, 1)
    im = plt.imshow(similarity_matrix, cmap='viridis', aspect='auto', interpolation='nearest')
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.title(f'{matrix_name} Similarity Matrix\n{n_cells} cells (ASC, OGC, OPC)', fontsize=12, fontweight='bold')
    plt.xlabel('Cell Index')
    plt.ylabel('Cell Index')
    
    # 2. Distribution
    plt.subplot(2, 3, 2)
    plt.hist(matrix_flat, bins=50, color=color_scheme['hist'], edgecolor='black', alpha=0.7)
    plt.title(f'{matrix_name} Distribution (off-diagonal)')
    plt.xlabel(f'{matrix_name} Value')
    plt.ylabel('Frequency')
    plt.grid(True, alpha=0.3)
    
    # Add distribution statistics
    dist_stats = f'Mean: {np.mean(matrix_flat):.3f}\nStd: {np.std(matrix_flat):.3f}'
    plt.text(0.95, 0.95, dist_stats, transform=plt.gca().transAxes, 
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # 3. MDS - Main plot
    plt.subplot(2, 3, (3, 6))
    
    try:
        mds = sklearn.manifold.MDS(
            n_components=2, 
            dissimilarity='precomputed', 
            random_state=42,
            max_iter=1000,
            eps=1e-9,
            normalized_stress='auto'
        )
        embedding = mds.fit_transform(1 - similarity_matrix)
        
        # Calculate optimal point size
        point_size = max(3, 2000 / n_cells)
        
        # Create scatter plot with color based on density
        if len(embedding) > 10:
            xy = embedding.T
            z = gaussian_kde(xy)(xy)
            scatter = plt.scatter(embedding[:, 0], embedding[:, 1], 
                                c=z, cmap=color_scheme['cmap'], s=point_size, alpha=0.7)
            plt.colorbar(scatter, label='Point Density')
        else:
            plt.scatter(embedding[:, 0], embedding[:, 1], 
                       c=color_scheme['points'], s=point_size*2, alpha=0.8)
        
        plt.title(f'{matrix_name} + MDS\n{n_cells} cells (ASC, OGC, OPC) | Stress: {mds.stress_:.4f}', 
                  fontsize=14, fontweight='bold')
        plt.xlabel('MDS Component 1')
        plt.ylabel('MDS Component 2')
        plt.grid(True, alpha=0.3)
        
        # Add statistics text box
        stats_text = f'{matrix_name} Statistics:\nMean: {np.mean(matrix_flat):.3f}\nStd: {np.std(matrix_flat):.3f}\nMin: {np.min(matrix_flat):.3f}\nMax: {np.max(matrix_flat):.3f}'
        plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
    except Exception as e:
        plt.text(0.5, 0.5, f'MDS failed:\n{str(e)}', ha='center', va='center', 
                 transform=plt.gca().transAxes, bbox=dict(facecolor='red', alpha=0.3))
        plt.title(f'{matrix_name} MDS - Failed')
    
    plt.tight_layout()
    return fig

# Create MoC comprehensive figure
print("Creating MoC comprehensive figure...")
moc_colors = {
    'hist': 'lightblue',
    'cmap': 'Blues',
    'points': 'blue'
}
moc_fig = create_comprehensive_figure(simi_moc, "MoC", moc_colors)
moc_fig.savefig(os.path.join(output_dir, "moc_comprehensive_analysis_ASC_OGC_OPC.png"), 
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close(moc_fig)

# Create AMI comprehensive figure
print("Creating AMI comprehensive figure...")
ami_colors = {
    'hist': 'lightcoral',
    'cmap': 'Reds',
    'points': 'red'
}
ami_fig = create_comprehensive_figure(simi_ami, "AMI", ami_colors)
ami_fig.savefig(os.path.join(output_dir, "ami_comprehensive_analysis_ASC_OGC_OPC.png"), 
                dpi=300, bbox_inches='tight', facecolor='white')
plt.close(ami_fig)

# Create combined comparison figure
print("Creating combined comparison figure...")
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Compute MDS for both matrices
try:
    mds = sklearn.manifold.MDS(
        n_components=2, 
        dissimilarity='precomputed', 
        random_state=42,
        max_iter=1000,
        eps=1e-9
    )
    
    embedding_moc = mds.fit_transform(1 - simi_moc)
    embedding_ami = mds.fit_transform(1 - simi_ami)
    
    point_size = max(3, 1500 / n_cells)
    
    # Row 1: MoC visualizations
    # MoC Heatmap
    im1 = axes[0,0].imshow(simi_moc, cmap='Blues', aspect='auto', interpolation='nearest')
    plt.colorbar(im1, ax=axes[0,0])
    axes[0,0].set_title(f'MoC Similarity Matrix\n({n_cells} cells, ASC/OGC/OPC)')
    axes[0,0].set_xlabel('Cell Index')
    axes[0,0].set_ylabel('Cell Index')
    
    # MoC Distribution
    moc_flat = simi_moc[np.triu_indices_from(simi_moc, k=1)]
    axes[0,1].hist(moc_flat, bins=50, color='lightblue', edgecolor='black', alpha=0.7)
    axes[0,1].set_title('MoC Distribution')
    axes[0,1].set_xlabel('MoC Value')
    axes[0,1].set_ylabel('Frequency')
    axes[0,1].grid(True, alpha=0.3)
    
    # MoC MDS with density coloring
    if len(embedding_moc) > 10:
        xy_moc = embedding_moc.T
        z_moc = gaussian_kde(xy_moc)(xy_moc)
        scatter1 = axes[0,2].scatter(embedding_moc[:, 0], embedding_moc[:, 1], 
                                    c=z_moc, cmap='Blues', s=point_size, alpha=0.7)
        plt.colorbar(scatter1, ax=axes[0,2], label='Point Density')
    else:
        scatter1 = axes[0,2].scatter(embedding_moc[:, 0], embedding_moc[:, 1], 
                                    c='blue', s=point_size, alpha=0.7)
    axes[0,2].set_title(f'MoC MDS ({n_cells} cells, ASC/OGC/OPC)')
    axes[0,2].set_xlabel('Component 1')
    axes[0,2].set_ylabel('Component 2')
    axes[0,2].grid(True, alpha=0.3)
    
    # Row 2: AMI visualizations
    # AMI Heatmap
    im2 = axes[1,0].imshow(simi_ami, cmap='Reds', aspect='auto', interpolation='nearest')
    plt.colorbar(im2, ax=axes[1,0])
    axes[1,0].set_title(f'AMI Similarity Matrix\n({n_cells} cells, ASC/OGC/OPC)')
    axes[1,0].set_xlabel('Cell Index')
    axes[1,0].set_ylabel('Cell Index')
    
    # AMI Distribution
    ami_flat = simi_ami[np.triu_indices_from(simi_ami, k=1)]
    axes[1,1].hist(ami_flat, bins=50, color='lightcoral', edgecolor='black', alpha=0.7)
    axes[1,1].set_title('AMI Distribution')
    axes[1,1].set_xlabel('AMI Value')
    axes[1,1].set_ylabel('Frequency')
    axes[1,1].grid(True, alpha=0.3)
    
    # AMI MDS with density coloring
    if len(embedding_ami) > 10:
        xy_ami = embedding_ami.T
        z_ami = gaussian_kde(xy_ami)(xy_ami)
        scatter2 = axes[1,2].scatter(embedding_ami[:, 0], embedding_ami[:, 1], 
                                    c=z_ami, cmap='Reds', s=point_size, alpha=0.7)
        plt.colorbar(scatter2, ax=axes[1,2], label='Point Density')
    else:
        scatter2 = axes[1,2].scatter(embedding_ami[:, 0], embedding_ami[:, 1], 
                                    c='red', s=point_size, alpha=0.7)
    axes[1,2].set_title(f'AMI MDS ({n_cells} cells, ASC/OGC/OPC)')
    axes[1,2].set_xlabel('Component 1')
    axes[1,2].set_ylabel('Component 2')
    axes[1,2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "combined_moc_ami_comparison_ASC_OGC_OPC.png"), 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
except Exception as e:
    print(f"Combined figure MDS failed: {e}")

# Create summary statistics
print("Creating summary statistics...")
summary_stats = {
    'Metric': ['MoC', 'AMI'],
    'Mean': [np.mean(simi_moc[np.triu_indices_from(simi_moc, k=1)]), 
             np.mean(simi_ami[np.triu_indices_from(simi_ami, k=1)])],
    'Std': [np.std(simi_moc[np.triu_indices_from(simi_moc, k=1)]), 
            np.std(simi_ami[np.triu_indices_from(simi_ami, k=1)])],
    'Min': [np.min(simi_moc[np.triu_indices_from(simi_moc, k=1)]), 
            np.min(simi_ami[np.triu_indices_from(simi_ami, k=1)])],
    'Max': [np.max(simi_moc[np.triu_indices_from(simi_moc, k=1)]), 
            np.max(simi_ami[np.triu_indices_from(simi_ami, k=1)])],
    'Median': [np.median(simi_moc[np.triu_indices_from(simi_moc, k=1)]), 
               np.median(simi_ami[np.triu_indices_from(simi_ami, k=1)])]
}

summary_df = pd.DataFrame(summary_stats)
summary_df.to_excel(os.path.join(output_dir, "similarity_summary_statistics_ASC_OGC_OPC.xlsx"), index=False)

# Export AMI MDS coordinates to TSV file
print("Exporting AMI MDS coordinates to TSV...")
try:
    # Compute AMI MDS with the same parameters as used in the figures
    mds_ami = sklearn.manifold.MDS(
        n_components=2, 
        dissimilarity='precomputed', 
        random_state=42,
        max_iter=1000,
        eps=1e-9
    )
    embedding_ami = mds_ami.fit_transform(1 - simi_ami)
    
    # Keep original barcodes with _tads suffix (don't clean them)
    # Create DataFrame with required columns
    # Remove '_tads' suffix from barcodes
    cleaned_barcodes = [barcode.replace('_tads', '') for barcode in barcodes]

# Create DataFrame
    ami_mds_df = pd.DataFrame({
        'barcode': cleaned_barcodes,
        'mds_component_1': embedding_ami[:, 0],
        'mds_component_2': embedding_ami[:, 1]
        })

# Add cell type information
    celltype_map = dict(zip(filtered_cells['barcode'], filtered_cells['celltype']))
    ami_mds_df['celltype'] = [
        celltype_map.get(bc, 'Unknown') for bc in cleaned_barcodes
]
    
    # Save as TSV
    ami_mds_df.to_csv(os.path.join(output_dir, "ami_mds_coordinates_ASC_OGC_OPC.tsv"), 
                     sep='\t', index=False)
    
    print(f"AMI MDS coordinates exported for {len(ami_mds_df)} cells with cell types")
    
except Exception as e:
    print(f"Error exporting AMI MDS coordinates: {e}")

print("\n" + "="*60)
print("VISUALIZATION COMPLETE!")
print("="*60)
print(f"Processed: {n_cells} cells (cell types: ASC, OGC, OPC)")
print(f"Files created:")
print(f"  - moc_comprehensive_analysis_ASC_OGC_OPC.png: MoC detailed analysis")
print(f"  - ami_comprehensive_analysis_ASC_OGC_OPC.png: AMI detailed analysis")
print(f"  - combined_moc_ami_comparison_ASC_OGC_OPC.png: Side-by-side comparison with density coloring")
print(f"  - similarity_summary_statistics_ASC_OGC_OPC.xlsx: Summary statistics")
print(f"  - ami_mds_coordinates_ASC_OGC_OPC.tsv: AMI MDS coordinates with cell types")
print(f"Results saved to: {output_dir}")
print("="*60)