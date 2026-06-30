#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 15:27:31 2026

@author: mozhganoroujlu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compartment Lattice Plots for Cell Types (ASC, OGC, OPC)
Creates sign-based and average compartment lattice plots for bulk Hi-C data
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import quantile_transform
from matplotlib.patches import Rectangle

EPS = 1e-5
N_COMPONENTS = 10
EPS_LOG = 1e-3

# -------------------- USER PATHS --------------------
txt_dir = "/bandnorm/bandnorm_txt_non_neuron/" #follow the pipeline for normalization in /scripts/scHiC normalization/ to get txt normalized contacts
cpg_txt = "/chrom1_cpg_ratios.txt" #run higashi_cpg.py to get this
celltype_file = "/non_neuron.tsv" # download it from /data directory
out_dir = "/output_compartments/lattice_plots/new/"
chrom = "chr1"
# -----------------------------------------------------

os.makedirs(out_dir, exist_ok=True)

# -------------------- FUNCTIONS (from your algorithm) --------------------
def load_cpg_txt(cpg_txt):
    """Load CpG ratios from text file with columns: bin_id, cpg_ratio"""
    df = pd.read_csv(cpg_txt, sep='\t')
    if 'bin_id' not in df.columns or 'cpg_ratio' not in df.columns:
        raise ValueError("CpG text file must contain 'bin_id' and 'cpg_ratio' columns")
    return df

def read_txt_contacts(path):
    try:
        df = pd.read_csv(path, sep='\t', header=None, names=['bin1','bin2','c'], engine='c', dtype={'bin1':int,'bin2':int,'c':float})
    except Exception:
        df = pd.read_csv(path, sep=r'\s+', header=None, names=['bin1','bin2','c'], engine='c', dtype={'bin1':int,'bin2':int,'c':float})
    return df

def build_dense_Q(contact_df, num_bins):
    Q = np.zeros((num_bins, num_bins), dtype=float)
    for i,j,v in contact_df[['bin1','bin2','c']].itertuples(index=False):
        if i < 0 or j < 0 or i >= num_bins or j >= num_bins:
            continue
        Q[i,j] += v
        Q[j,i] += v
    np.fill_diagonal(Q, 0.0)
    return Q

def compute_decay(Q):
    n = Q.shape[0]
    decay = np.zeros(n, dtype=float)
    for d in range(n):
        diag = np.diag(Q, d)
        decay[d] = np.nanmean(diag) if diag.size > 0 else 0.0
    decay = np.where(decay <= 0, EPS, decay)
    return decay

def build_E_from_Q(Q, decay):
    n = Q.shape[0]
    E = np.zeros_like(Q, dtype=float)
    for d in range(1, n):
        if decay[d] == 0: continue
        rows = np.arange(0, n-d)
        cols = rows + d
        E[rows, cols] = (Q[rows, cols] + EPS) / (decay[d] + EPS)
        E[cols, rows] = E[rows, cols]
    np.fill_diagonal(E, 1.0)
    return E

def safe_corr_from_E(E, eps_log):
    M = np.log2(E + eps_log)
    row_std = np.nanstd(M, axis=1)
    const_rows = np.where(row_std == 0)[0]
    if const_rows.size > 0:
        M[const_rows, :] += np.random.normal(scale=1e-6, size=(len(const_rows), M.shape[1]))
    C = np.corrcoef(M)
    return np.nan_to_num(C, nan=0.0, posinf=0.0, neginf=0.0)

def compute_bulk_Q(txt_paths, num_bins):
    """Compute bulk contact matrix by aggregating all single-cell matrices"""
    bulk_Q = np.zeros((num_bins, num_bins), dtype=float)
    used_barcodes = []
    for p in txt_paths:
        try:
            df = read_txt_contacts(p)
            Q = build_dense_Q(df, num_bins)
            bulk_Q += Q
            used_barcodes.append(p)
        except Exception as e:
            print(f"Error processing {p}: {e}")
            continue
    if len(used_barcodes) == 0:
        raise RuntimeError("No barcodes could be processed")
    return bulk_Q, used_barcodes

def fit_bulk_pca(bulk_Q, n_components, eps_log):
    decay = compute_decay(bulk_Q)
    E_bulk = build_E_from_Q(bulk_Q, decay)
    C_bulk = safe_corr_from_E(E_bulk, eps_log)
    pca = PCA(n_components=n_components)
    pcs_bulk = pca.fit_transform(C_bulk)
    return pca, pcs_bulk

def orient_bulk_by_cpg_extremes(pcs_bulk, cpg_values):
    """
    Use Top 10% vs Bottom 10% means for orientation
    Returns: flip = True if compartments need to be reversed
    """
    pc1 = pcs_bulk[:, 0]
    
    # Get top 10% and bottom 10% of PC1 values
    top_10_threshold = np.quantile(pc1, 0.9)
    bottom_10_threshold = np.quantile(pc1, 0.1)
    
    top_10_mask = pc1 > top_10_threshold
    bottom_10_mask = pc1 < bottom_10_threshold
    
    # Calculate mean CpG ratios for extremes
    top_10_cpg_mean = np.nanmean(cpg_values[top_10_mask])
    bottom_10_cpg_mean = np.nanmean(cpg_values[bottom_10_mask])
    
    print(f"Top 10% CpG mean: {top_10_cpg_mean:.4f}, Bottom 10% CpG mean: {bottom_10_cpg_mean:.4f}")
    
    # If top compartments (putative A) have LOWER CpG than bottom compartments (putative B), flip
    flip = top_10_cpg_mean < bottom_10_cpg_mean
    
    if flip:
        print("  → FLIPPING compartments (A compartments should have higher CpG)")
    else:
        print("  → Keeping original orientation")
    
    return flip

def quantile_normalize_scores(scores, valid_idx, num_bins):
    """
    Quantile normalize scores using maximum number of quantiles as 20% of number of bins
    and map to range [-2.5, 2.5]
    """
    n_quantiles = max(2, int(0.2 * num_bins))  # 20% of bins, minimum 2 quantiles
    
    if len(valid_idx) == 0:
        return scores
    
    # Extract valid scores for normalization
    valid_scores = scores[valid_idx]
    
    # Perform quantile normalization
    try:
        # Use uniform distribution as base and then map to desired range
        q_vals = quantile_transform(
            valid_scores.reshape(1, -1), 
            axis=1, 
            output_distribution='uniform',
            n_quantiles=min(n_quantiles, len(valid_scores))
        )[0]
        
        # Map from [0,1] to [-2.5, 2.5]
        # 0 -> -2.5, 0.5 -> 0, 1 -> 2.5
        normalized_scores = q_vals * 5 - 2.5
        
        # Apply normalized scores back to valid positions
        scores_normalized = scores.copy()
        scores_normalized[valid_idx] = normalized_scores
        
        return scores_normalized
        
    except Exception as e:
        print(f"Warning: Quantile normalization failed, using original scores: {e}")
        return scores

# -------------------- NEW PLOTTING FUNCTIONS (your requests) --------------------
def plot_compartment_sign_lattice(celltype, compartment_scores, out_dir):
    """
    Create a lattice plot showing sign matching between compartment scores.
    Blue: Both positive, Red: Both negative, White: Different signs
    """
    print(f"Creating sign-based lattice plot for {celltype}")
    
    n_bins = len(compartment_scores)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create an RGB image array (height x width x 3)
    # Initialize with white (empty cells)
    rgb_image = np.ones((n_bins, n_bins, 3))
    
    # Get sign of each bin (True for positive, False for negative)
    signs = compartment_scores > 0
    
    # Fill the lattice based on sign matching
    for i in range(n_bins):
        for j in range(n_bins):
            if signs[i] and signs[j]:  # Both positive
                rgb_image[i, j] = [0, 0, 1]  # Blue
            elif not signs[i] and not signs[j]:  # Both negative
                rgb_image[i, j] = [1, 0, 0]  # Red
            # Different signs remain white (already set)
    
    # Display the image
    ax.imshow(rgb_image, interpolation='none', aspect='auto')
    
    # Customize the plot
    ax.set_title(f'Compartment Sign Matching - {celltype}\n(Blue: Both +, Red: Both -, White: Different signs)', 
                fontsize=14, fontweight='bold')
    ax.set_xlabel('Bin ID', fontsize=12)
    ax.set_ylabel('Bin ID', fontsize=12)
    
    # Add color legend
    legend_elements = [
        Rectangle((0, 0), 1, 1, facecolor='blue', edgecolor='black', label='Both Positive'),
        Rectangle((0, 0), 1, 1, facecolor='red', edgecolor='black', label='Both Negative'),
        Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='black', label='Different Signs')
    ]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1))
    
    # Add statistics
    n_total = n_bins * n_bins
    n_both_pos = np.sum(signs[:, None] & signs[None, :])
    n_both_neg = np.sum((~signs[:, None]) & (~signs[None, :]))
    n_diff = n_total - n_both_pos - n_both_neg
    
    stats_text = f'Total pairs: {n_total:,}\nBoth +: {n_both_pos:,} ({n_both_pos/n_total*100:.1f}%)\nBoth -: {n_both_neg:,} ({n_both_neg/n_total*100:.1f}%)\nDifferent: {n_diff:,} ({n_diff/n_total*100:.1f}%)'
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'sign_lattice_{celltype}.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()


def plot_average_compartment_lattice(celltype, compartment_scores, out_dir):
    """
    Create a lattice plot showing average compartment scores for bins with matching signs.
    More negative = darker red, More positive = darker blue
    """
    print(f"Creating average compartment lattice plot for {celltype}")
    
    n_bins = len(compartment_scores)
    
    # Create matrix for average scores (initialize with NaN for different signs)
    avg_scores = np.full((n_bins, n_bins), np.nan)
    
    # Get sign of each bin
    signs = compartment_scores > 0
    
    # Fill the lattice with average scores where signs match
    for i in range(n_bins):
        for j in range(n_bins):
            if signs[i] == signs[j]:  # Same sign
                avg_scores[i, j] = (compartment_scores[i] + compartment_scores[j]) / 2
    
    # Create figure (single heatmap only)
    fig, ax1 = plt.subplots(1, 1, figsize=(10, 9))
    
    # Heatmap of average scores
    masked_avg = np.ma.masked_where(np.isnan(avg_scores), avg_scores)
    
    cmap = plt.cm.RdBu_r
    cmap.set_bad('white')  # Set NaN values to white
    
    im = ax1.imshow(masked_avg, cmap=cmap, vmin=-2.5, vmax=2.5, 
                   interpolation='none', aspect='auto')
    
    ax1.set_title(f'Average Compartment Score (Same Sign Only) - {celltype}\n(More negative = darker blue, More positive = darker red)', 
                 fontsize=14, fontweight='bold')
    ax1.set_xlabel('Bin ID', fontsize=12)
    ax1.set_ylabel('Bin ID', fontsize=12)
    
    cbar = plt.colorbar(im, ax=ax1, shrink=0.8)
    cbar.set_label('Average Compartment Score', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'avg_compartment_lattice_{celltype}.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()


# -------------------- MAIN --------------------
def main():
    print("="*60)
    print("COMPARTMENT LATTICE PLOTS FOR CELL TYPES")
    print("="*60)
    
    # Load CpG bins
    print("\nLoading CpG bins...")
    cpg_bins = load_cpg_txt(cpg_txt)
    num_bins = cpg_bins.shape[0]
    print(f"Number of bins: {num_bins}")
    
    # Load cell type info
    print("\nLoading cell type information...")
    celltype_df = pd.read_csv(celltype_file, sep='\t')
    barcode_to_celltype = dict(zip(celltype_df['barcode'], celltype_df['celltype']))
    
    # Collect txt files and group by cell type
    txt_files = [f for f in os.listdir(txt_dir) if f.endswith('.txt')]
    celltype_groups = {'ASC': [], 'OGC': [], 'OPC': []}
    
    for f in txt_files:
        barcode = os.path.splitext(f)[0]
        celltype = barcode_to_celltype.get(barcode)
        if celltype in celltype_groups:
            celltype_groups[celltype].append(os.path.join(txt_dir, f))
    
    print(f"\nFound files per cell type:")
    for celltype, files in celltype_groups.items():
        print(f"  {celltype}: {len(files)} files")
    
    # Dictionary to store bulk compartment scores
    bulk_scores = {}
    
    # Process each cell type
    for celltype in ['ASC', 'OGC', 'OPC']:
        print(f"\n{'-'*40}")
        print(f"Processing {celltype}")
        print(f"{'-'*40}")
        
        files = celltype_groups[celltype]
        if len(files) == 0:
            print(f"Warning: No files found for {celltype}")
            continue
        
        # Step 1: Compute bulk contact matrix by aggregating all single-cell matrices
        print(f"Aggregating {len(files)} single-cell matrices...")
        bulk_Q, used_files = compute_bulk_Q(files, num_bins)
        print(f"Successfully aggregated {len(used_files)} files")
        
        # Step 2: Fit PCA on bulk data
        print("Fitting PCA on bulk data...")
        pca, pcs_bulk = fit_bulk_pca(bulk_Q, N_COMPONENTS, EPS_LOG)
        
        # Step 3: Orient compartments using CpG extremes
        print("Orienting compartments using CpG extremes...")
        flip = orient_bulk_by_cpg_extremes(pcs_bulk, cpg_bins['cpg_ratio'].values)
        
        # Step 4: Get oriented PC1
        pc1 = pcs_bulk[:, 0]
        if flip:
            pc1 = -pc1
        
        # Step 5: Apply quantile normalization (using your algorithm)
        print("Applying quantile normalization...")
        binfilter = np.any(bulk_Q != 0, axis=1)
        valid_idx = np.where(binfilter)[0]
        pc1_normalized = quantile_normalize_scores(pc1, valid_idx, num_bins)
        
        print(f"PC1 range after normalization: [{pc1_normalized.min():.3f}, {pc1_normalized.max():.3f}]")
        
        # Store the scores
        bulk_scores[celltype] = pc1_normalized
        
        # Step 6: Create your requested plots
        print(f"Creating lattice plots for {celltype}...")
        plot_compartment_sign_lattice(celltype, pc1_normalized, out_dir)
        plot_average_compartment_lattice(celltype, pc1_normalized, out_dir)
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Output directory: {out_dir}")
    print("\nGenerated files:")
    for celltype in bulk_scores.keys():
        print(f"  {celltype}:")
        print(f"    - sign_lattice_{celltype}.png")
        print(f"    - avg_compartment_lattice_{celltype}.png")
    print("\nDone!")


if __name__ == "__main__":
    main()
