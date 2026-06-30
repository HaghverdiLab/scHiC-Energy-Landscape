import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import os

# Set random seed for reproducibility
np.random.seed(42)

def generate_synthetic_barcodes(input_file, output_file, 
                                asc_multiplier=4, other_multiplier=2, 
                                noise_scale=0.1, energy_noise_scale=0.15):
    """
    Generate synthetic barcodes to balance cell type counts and create visualization
    
    Parameters:
    - input_file: path to input TSV file
    - output_file: path to output TSV file with synthetic data
    - asc_multiplier: number of synthetic points per original ASC barcode
    - other_multiplier: number of synthetic points per original OGC/OPC barcode
    - noise_scale: scale of Gaussian noise for MDS coordinates
    - energy_noise_scale: scale of Gaussian noise for energy values
    """
    
    # Read the data
    df = pd.read_csv(input_file, sep='\t')
    
    # Get counts per cell type
    celltype_counts = df['celltype'].value_counts()
    print(f"Original cell type counts:")
    for celltype, count in celltype_counts.items():
        print(f"  {celltype}: {count}")
    
    # Create synthetic data
    synthetic_rows = []
    
    for idx, row in df.iterrows():
        celltype = row['celltype']
        
        # Determine number of synthetic points based on cell type
        if celltype == 'ASC':
            n_synthetic = asc_multiplier
        elif celltype in ['OGC', 'OPC']:
            n_synthetic = other_multiplier
        else:
            n_synthetic = 1  # Default for any other cell types
        
        # Generate synthetic barcodes
        for i in range(n_synthetic):
            # Create a new barcode name
            synthetic_barcode = f"{row['barcode']}_syn_{i+1}"
            
            # Add Gaussian noise to MDS components
            synthetic_mds1 = row['mds_component_1'] + np.random.normal(0, noise_scale)
            synthetic_mds2 = row['mds_component_2'] + np.random.normal(0, noise_scale)
            
            # Add Gaussian noise to energy
            synthetic_energy = row['energy'] + np.random.normal(0, energy_noise_scale)
            
            # Create synthetic row
            synthetic_row = {
                'barcode': synthetic_barcode,
                'mds_component_1': synthetic_mds1,
                'mds_component_2': synthetic_mds2,
                'celltype': celltype,
                'energy': synthetic_energy,
                'category': row['category'],
                'state': 'generated'
            }
            synthetic_rows.append(synthetic_row)
    
    # Create synthetic dataframe
    df_synthetic = pd.DataFrame(synthetic_rows)
    
    # Add state column to original data
    df_original = df.copy()
    df_original['state'] = 'original'
    
    # Combine original and synthetic data
    df_combined = pd.concat([df_original, df_synthetic], ignore_index=True)
    
    # Save to TSV
    df_combined.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved combined data to {output_file}")
    
    # Print new counts
    print(f"\nNew cell type counts after generation:")
    for celltype in ['ASC', 'OGC', 'OPC']:
        count = len(df_combined[df_combined['celltype'] == celltype])
        print(f"  {celltype}: {count}")
    
    return df_combined, df_original

def plot_original_mds(df_original, plot_file):
    """
    Create MDS plot with only original data points (all circles)
    Highlight the 5 lowest energy points with stars
    """
    # Set up the plot
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Define more distinct colors: yellow, green, red
    celltype_colors = {
        'ASC': '#FFD700',      # Bright yellow/gold
        'OGC': '#2ECC71',      # Emerald green
        'OPC': '#E74C3C',      # Bright red
        'Other': '#95A5A6'     # Gray for any other types
    }
    
    # Plot original data (all circles)
    for celltype in df_original['celltype'].unique():
        subset = df_original[df_original['celltype'] == celltype]
        color = celltype_colors.get(celltype, celltype_colors['Other'])
        ax.scatter(subset['mds_component_1'], subset['mds_component_2'],
                  c=color, label=f'{celltype}', marker='o',
                  s=70, alpha=0.7, edgecolors='black', linewidth=0.5)
    
    # Find the 5 lowest energy points in original data
    lowest_energy_points = df_original.nsmallest(5, 'energy')
    
    # Highlight these points with stars
    star_scatter = ax.scatter(lowest_energy_points['mds_component_1'], 
                             lowest_energy_points['mds_component_2'],
                             marker='*', s=300, c='black', 
                             edgecolors='white', linewidth=1.5,
                             zorder=10, label='5 Lowest Energy Points')
    
    # Create a custom legend for the lowest energy points
    energy_info = []
    for idx, row in lowest_energy_points.iterrows():
        info = f"Energy: {row['energy']:.3f}, {row['celltype']}"
        energy_info.append(info)
    
    # Add text box with energy information
    textstr = '\n'.join(energy_info)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=props)
    
    # Customize the plot
    ax.set_xlabel('MDS Component 1', fontsize=12)
    ax.set_ylabel('MDS Component 2', fontsize=12)
    ax.set_title('MDS Plot: Original Barcodes Only\nStars indicate 5 lowest energy values', 
                 fontsize=14, fontweight='bold')
    
    # Add legend
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), bbox_to_anchor=(1.05, 1), 
              loc='upper left', fontsize=10)
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Adjust layout to prevent legend cutoff
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Saved original data MDS plot to {plot_file}")
    
    # Print the lowest energy points to console
    print("\n" + "="*50)
    print("5 LOWEST ENERGY POINTS (ORIGINAL DATA):")
    print("="*50)
    for idx, row in lowest_energy_points.iterrows():
        print(f"  Energy: {row['energy']:.6f} | Celltype: {row['celltype']} | Barcode: {row['barcode']}")
    
    # Show the plot (if running interactively)
    plt.show()

def plot_combined_mds(df_combined, plot_file):
    """
    Create MDS plot with original (circles) and generated (triangles) points
    Highlight the 5 lowest energy points with stars
    """
    # Set up the plot
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Define more distinct colors: yellow, green, red
    celltype_colors = {
        'ASC': '#FFD700',      # Bright yellow/gold
        'OGC': '#2ECC71',      # Emerald green
        'OPC': '#E74C3C',      # Bright red
        'Other': '#95A5A6'     # Gray for any other types
    }
    
    # Plot original data (circles)
    for celltype in df_combined['celltype'].unique():
        subset = df_combined[(df_combined['celltype'] == celltype) & 
                            (df_combined['state'] == 'original')]
        color = celltype_colors.get(celltype, celltype_colors['Other'])
        ax.scatter(subset['mds_component_1'], subset['mds_component_2'],
                  c=color, label=f'{celltype} (original)', marker='o',
                  s=70, alpha=0.7, edgecolors='black', linewidth=0.5)
    
    # Plot generated data (triangles)
    for celltype in df_combined['celltype'].unique():
        subset = df_combined[(df_combined['celltype'] == celltype) & 
                            (df_combined['state'] == 'generated')]
        color = celltype_colors.get(celltype, celltype_colors['Other'])
        ax.scatter(subset['mds_component_1'], subset['mds_component_2'],
                  c=color, label=f'{celltype} (generated)', marker='^',
                  s=50, alpha=0.5, edgecolors='black', linewidth=0.3)
    
    # Find the 5 lowest energy points in combined data
    lowest_energy_points = df_combined.nsmallest(5, 'energy')
    
    # Highlight these points with stars
    star_scatter = ax.scatter(lowest_energy_points['mds_component_1'], 
                             lowest_energy_points['mds_component_2'],
                             marker='*', s=300, c='black', 
                             edgecolors='white', linewidth=1.5,
                             zorder=10, label='5 Lowest Energy Points')
    
    # Create a custom legend for the lowest energy points
    energy_info = []
    for idx, row in lowest_energy_points.iterrows():
        info = f"Energy: {row['energy']:.3f}, {row['celltype']} ({row['state']})"
        energy_info.append(info)
    
    # Add text box with energy information
    textstr = '\n'.join(energy_info)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=props)
    
    # Customize the plot
    ax.set_xlabel('MDS Component 1', fontsize=12)
    ax.set_ylabel('MDS Component 2', fontsize=12)
    ax.set_title('MDS Plot: Original (circles) vs Generated (triangles) Barcodes\nStars indicate 5 lowest energy values', 
                 fontsize=14, fontweight='bold')
    
    # Add legend
    handles, labels = ax.get_legend_handles_labels()
    # Remove duplicate entries
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), bbox_to_anchor=(1.05, 1), 
              loc='upper left', fontsize=10)
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Adjust layout to prevent legend cutoff
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Saved combined data MDS plot to {plot_file}")
    
    # Print the lowest energy points to console
    print("\n" + "="*50)
    print("5 LOWEST ENERGY POINTS (COMBINED DATA):")
    print("="*50)
    for idx, row in lowest_energy_points.iterrows():
        print(f"  Energy: {row['energy']:.6f} | Celltype: {row['celltype']} | State: {row['state']} | Barcode: {row['barcode']}")
    
    # Show the plot (if running interactively)
    plt.show()

def main():
    # File paths
    input_file = '/Users/mozhganoroujlu/Desktop/MOZHGUN/cell_fate/hi_c/codes_figures/folders/normalized_contacts/non_neuron_filtered_celltypes.tsv'  # Replace with your actual input file
    output_file = '/Users/mozhganoroujlu/Desktop/MOZHGUN/cell_fate/hi_c/codes_figures/folders/normalized_contacts/generated_barcodes.tsv'
    plot_original_file = '/Users/mozhganoroujlu/Desktop/MOZHGUN/cell_fate/hi_c/codes_figures/folders/normalized_contacts/mds_plot.png'
    plot_combined_file = '/Users/mozhganoroujlu/Desktop/MOZHGUN/cell_fate/hi_c/codes_figures/folders/normalized_contacts/mds_plot_combined.png'
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        print("Please update the 'input_file' variable with the correct path to your data.")
        return
    
    # Generate synthetic data
    # Adjust these parameters as needed:
    # - asc_multiplier: number of synthetic points per ASC barcode (default: 4)
    # - other_multiplier: number of synthetic points per OGC/OPC barcode (default: 2)
    # - noise_scale: amount of spatial noise (default: 0.1)
    # - energy_noise_scale: amount of energy noise (default: 0.15)
    df_combined, df_original = generate_synthetic_barcodes(
        input_file=input_file,
        output_file=output_file,
        asc_multiplier=5,
        other_multiplier=2,
        noise_scale=0.1,
        energy_noise_scale=0.15
    )
    
    # Create the first MDS plot (original data only)
    plot_original_mds(df_original, plot_original_file)
    
    # Create the second MDS plot (combined original + generated)
    plot_combined_mds(df_combined, plot_combined_file)
    
    # Display statistics
    print("\n" + "="*50)
    print("SUMMARY STATISTICS")
    print("="*50)
    print(f"Total original barcodes: {len(df_combined[df_combined['state'] == 'original'])}")
    print(f"Total generated barcodes: {len(df_combined[df_combined['state'] == 'generated'])}")
    print(f"Total combined barcodes: {len(df_combined)}")
    print("\nCell type distribution after generation:")
    for celltype in sorted(df_combined['celltype'].unique()):
        original = len(df_combined[(df_combined['celltype'] == celltype) & 
                                   (df_combined['state'] == 'original')])
        generated = len(df_combined[(df_combined['celltype'] == celltype) & 
                                    (df_combined['state'] == 'generated')])
        total = original + generated
        print(f"  {celltype}: {total} total ({original} original + {generated} generated)")

if __name__ == "__main__":
    main()
