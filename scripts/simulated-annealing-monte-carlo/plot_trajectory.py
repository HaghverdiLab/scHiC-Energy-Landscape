#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spatial Simulated Annealing - Single Run with Smooth Trajectory Plot
Optimized for Apple Silicon (M1/M2/M3/M4)
"""
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from scipy.interpolate import CubicSpline
import random
import math

# Try to import tqdm for progress bar, but provide fallback
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    print("tqdm not installed. Progress bar will be disabled.")
    # Create a simple fallback
    class tqdm:
        def __init__(self, *args, **kwargs):
            self.total = kwargs.get('total', 0)
            self.desc = kwargs.get('desc', '')
            self.n = 0
            
        def update(self, n=1):
            self.n += n
            if self.n % 100 == 0:  # Print every 100 iterations
                print(f"{self.desc}: {self.n}/{self.total}", end='\r')
                
        def close(self):
            print()  # New line at end

# ────────────────────────────────────────────────
# Helper: ensure output directory exists
# ────────────────────────────────────────────────
def ensure_dir(path: str):
    Path(os.path.dirname(path)).mkdir(parents=True, exist_ok=True)


# ────────────────────────────────────────────────
# Vectorized neighbor scoring (hot path)
# ────────────────────────────────────────────────
def get_best_neighbor_probs(
    current_idx: int,
    neighbor_indices: np.ndarray,
    coords: np.ndarray,
    energies: np.ndarray,
    spatial_weight: float = 0.7,
    energy_weight: float = 0.3
) -> tuple[np.ndarray, float]:
    """
    Returns selection probabilities and total score sum (vectorized)
    """
    current_coord = coords[current_idx]
    current_energy = energies[current_idx]

    # Vectorized spatial distance
    neigh_coords = coords[neighbor_indices]
    spatial_dists = np.sqrt(np.sum((neigh_coords - current_coord) ** 2, axis=1))
    spatial_scores = 1.0 / (1.0 + spatial_dists)

    # Vectorized energy similarity
    neigh_energies = energies[neighbor_indices]
    energy_diffs = np.abs(neigh_energies - current_energy)
    energy_scores = 1.0 / (1.0 + energy_diffs)

    # Combined score
    total_scores = spatial_weight * spatial_scores + energy_weight * energy_scores
    score_sum = total_scores.sum()

    if score_sum > 0:
        probs = total_scores / score_sum
        return probs, score_sum
    else:
        # uniform fallback
        probs = np.full(len(neighbor_indices), 1.0 / len(neighbor_indices))
        return probs, 1.0


def select_neighbor(
    current_idx: int,
    neighborhoods: dict,
    coords: np.ndarray,
    energies: np.ndarray,
    spatial_weight: float,
    energy_weight: float
) -> int | None:
    """Selects next cell index or returns None"""
    neigh_info = neighborhoods.get(current_idx)
    if not neigh_info or len(neigh_info['indices']) == 0:
        return None

    indices = neigh_info['indices']
    probs, _ = get_best_neighbor_probs(
        current_idx, indices, coords, energies, spatial_weight, energy_weight
    )

    # Choose one neighbor
    chosen_pos = np.random.choice(len(indices), p=probs)
    return indices[chosen_pos]


# ────────────────────────────────────────────────
# Smooth trajectory plotting with energy gradient using average energy per segment
# ────────────────────────────────────────────────
def plot_smooth_trajectory_with_energy(
    original_df: pd.DataFrame, 
    trajectory_data: list, 
    energies_at_states: list,
    skip_points: int = 5  # Connect every Nth point and average energies in between
):
    """
    Plot a smooth trajectory line colored by average energy of skipped segments
    
    Parameters:
    - original_df: DataFrame with all cell data
    - trajectory_data: List of [x, y, barcode] coordinates for accepted states
    - energies_at_states: List of energy values corresponding to each point
    - skip_points: Number of points to skip for smoother line (connect every Nth point)
                 and average energies of the skipped points for each segment
    """
    plt.figure(figsize=(14, 10))

    # Plot all cells in background
    cell_types = original_df['celltype'].unique()
    colors = plt.cm.Set3(np.linspace(0, 1, len(cell_types)))

    for i, ct in enumerate(cell_types):
        mask = original_df['celltype'] == ct
        plt.scatter(
            original_df.loc[mask, 'mds_component_1'],
            original_df.loc[mask, 'mds_component_2'],
            c=[colors[i]], 
            label=ct, 
            alpha=0.2, 
            s=10,
            edgecolors='none'
        )

    # Convert trajectory data to numpy arrays
    traj_coords = np.array([[x, y] for x, y, b in trajectory_data])
    traj_energies = np.array(energies_at_states)
    
    if len(traj_coords) < 2:
        print("Not enough trajectory points to plot")
        return
    
    # Create segments based on skip_points
    # Each segment will connect point i to point i+skip_points
    # and its color will be based on the average energy of points between them
    
    segment_starts = []
    segment_ends = []
    segment_avg_energies = []
    
    for i in range(0, len(traj_coords) - skip_points, skip_points):
        start_idx = i
        end_idx = i + skip_points
        
        # Get coordinates for start and end
        start_coord = traj_coords[start_idx]
        end_coord = traj_coords[end_idx]
        
        # Calculate average energy of all points between start and end (inclusive)
        energies_in_segment = traj_energies[start_idx:end_idx + 1]
        avg_energy = np.mean(energies_in_segment)
        
        segment_starts.append(start_coord)
        segment_ends.append(end_coord)
        segment_avg_energies.append(avg_energy)
    
    # Add the final segment if there are remaining points
    if len(traj_coords) > skip_points and (len(traj_coords) - 1) % skip_points != 0:
        last_start = (len(traj_coords) - 1) - ((len(traj_coords) - 1) % skip_points)
        if last_start < len(traj_coords) - 1:
            start_coord = traj_coords[last_start]
            end_coord = traj_coords[-1]
            energies_in_segment = traj_energies[last_start:]
            avg_energy = np.mean(energies_in_segment)
            
            segment_starts.append(start_coord)
            segment_ends.append(end_coord)
            segment_avg_energies.append(avg_energy)
    
    # Convert to arrays for plotting
    segment_starts = np.array(segment_starts)
    segment_ends = np.array(segment_ends)
    segment_avg_energies = np.array(segment_avg_energies)
    
    # Create smooth interpolation for each segment
    from matplotlib.collections import LineCollection
    
    all_segments = []
    segment_colors = []
    
    for i in range(len(segment_starts)):
        start = segment_starts[i]
        end = segment_ends[i]
        avg_energy = segment_avg_energies[i]
        
        # Create a few intermediate points for smooth curve
        # Use linear interpolation for simplicity
        t = np.linspace(0, 1, 10)
        x = start[0] + t * (end[0] - start[0])
        y = start[1] + t * (end[1] - start[1])
        
        # Create line segments for this portion
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
        for segment in segments:
            all_segments.append(segment)
            segment_colors.append(avg_energy)
    
    # Create LineCollection
    lc = LineCollection(all_segments, cmap='viridis_r', linewidth=1, alpha=0.8)
    lc.set_array(np.array(segment_colors))
    plt.gca().add_collection(lc)
    
    # Add colorbar
    cbar = plt.colorbar(lc, label='Average Energy (per segment)', shrink=0.8)
    cbar.ax.tick_params(labelsize=10)
    
    # Plot start and end points
    plt.scatter(traj_coords[0, 0], traj_coords[0, 1], 
               marker='o', s=200, color='green', edgecolors='white', 
               linewidth=2, zorder=5, label='Start')
    
    plt.scatter(traj_coords[-1, 0], traj_coords[-1, 1], 
               marker='s', s=200, color='red', edgecolors='white', 
               linewidth=2, zorder=5, label='End')
    
    # Add energy text for start and end
    plt.annotate(f'Energy: {traj_energies[0]:.3f}', 
                (traj_coords[0, 0], traj_coords[0, 1]),
                xytext=(10, 10), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
    
    plt.annotate(f'Energy: {traj_energies[-1]:.3f}', 
                (traj_coords[-1, 0], traj_coords[-1, 1]),
                xytext=(10, -10), textcoords='offset points',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))

    plt.legend(loc='upper right', fontsize=12, framealpha=0.9)
    
    # Get cell types for start and end
    if len(trajectory_data) > 0:
        start_barcode = trajectory_data[0][2]
        end_barcode = trajectory_data[-1][2]
        start_type = original_df[original_df['barcode'] == start_barcode]['celltype'].iloc[0] if not original_df[original_df['barcode'] == start_barcode].empty else '?'
        end_type = original_df[original_df['barcode'] == end_barcode]['celltype'].iloc[0] if not original_df[original_df['barcode'] == end_barcode].empty else '?'
    else:
        start_type = end_type = '?'
    
    plt.title(f'Spatial Trajectory - Smooth Path with Energy Averages\nStart: {start_type} → End: {end_type}\n(Averaging over {skip_points} points per segment)', 
              fontsize=14, fontweight='bold')
    plt.xlabel('MDS Component 1', fontsize=12)
    plt.ylabel('MDS Component 2', fontsize=12)
    plt.tight_layout()
    
    # Save the figure
    output_dir = "/Users/mozhganoroujlu/Desktop/SA/score/results/modif/generated/"
    ensure_dir(output_dir)
    plt.savefig(os.path.join(output_dir, f'smooth_trajectory_energy_gradient_avg{skip_points}.png'), 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print some statistics about the segments
    print(f"\nSegment Statistics (averaging over {skip_points} points per segment):")
    print(f"Number of segments: {len(segment_avg_energies)}")
    print(f"Average energy range: {segment_avg_energies.min():.4f} to {segment_avg_energies.max():.4f}")
    print(f"Overall energy change: {traj_energies[-1] - traj_energies[0]:.4f}")


# ────────────────────────────────────────────────
# Main SA function - Single Run Only
# ────────────────────────────────────────────────
def run_single_spatial_simulated_annealing(
    tsv_file: str,
    k_neighbors: int = 150,
    spatial_weight: float = 0.6,
    energy_weight: float = 0.4,
    T_start: float = 0.2,
    T_min: float = 1e-6,
    alpha: float = 0.98,
    outer_iterations: int = 200,
    inner_iterations: int = 50,
    skip_points_for_plot: int = 5,  # Connect every 5 accepted states and average their energies
):
    """
    Runs spatial neighborhood-based simulated annealing once and creates
    a smooth trajectory plot colored by average energy of skipped segments.
    """
    print("Loading data...")
    df = pd.read_csv(tsv_file, sep='\t')

    # Filter Non-neuron
    glut_df = df[df['category'] == 'Non-neuron'].reset_index(drop=True)
    print(f"Working with {len(glut_df):,} Non-neuron cells")

    if 'energy' not in glut_df.columns:
        raise ValueError("Column 'energy' not found in data")

    # ── Precompute fast numpy arrays ──
    coords = glut_df[['mds_component_1', 'mds_component_2']].to_numpy()
    energies = glut_df['energy'].to_numpy()
    barcodes = glut_df['barcode'].to_numpy()
    celltypes = glut_df['celltype'].to_numpy()

    # ── Build spatial k-NN graph ──
    print("Building spatial k-NN graph...")
    start = time.time()
    nbrs = NearestNeighbors(n_neighbors=k_neighbors + 1, algorithm='ball_tree', n_jobs=-1)
    nbrs.fit(coords)
    distances, indices = nbrs.kneighbors(coords)
    print(f"k-NN built in {time.time() - start:.2f} seconds")

    # Build neighborhoods (exclude self)
    neighborhoods = {}
    for i in range(len(glut_df)):
        neighborhoods[i] = {
            'indices': indices[i, 1:],
            'distances': distances[i, 1:],
            'barcode': barcodes[i]
        }

    # ── Find OPC starting point ──
    opc_mask = celltypes == 'OPC'
    opc_indices = np.where(opc_mask)[0]
    if len(opc_indices) == 0:
        raise ValueError("No OPC cells found in data")

    # Randomly select one OPC as starting point
    start_idx = np.random.choice(opc_indices)
    
    # Initialize trajectory storage
    trajectory_coords = []  # Store [x, y, barcode]
    trajectory_energies = []

    current_idx = int(start_idx)
    current_energy = energies[current_idx]
    current_barcode = barcodes[current_idx]
    current_celltype = celltypes[current_idx]
    
    # Store initial state
    trajectory_coords.append([coords[current_idx, 0], coords[current_idx, 1], current_barcode])
    trajectory_energies.append(current_energy)
    
    print(f"Starting from OPC cell: {current_barcode}")
    print(f"Initial energy: {current_energy:.4f}")

    T = T_start
    accept_count_total = 0
    
    # Track accepted states for trajectory
    accepted_states = [(coords[current_idx, 0], coords[current_idx, 1], current_barcode, current_energy)]

    # Progress bar
    total_iterations = outer_iterations * inner_iterations
    pbar = tqdm(total=total_iterations, desc="SA iterations")

    for outer in range(outer_iterations):
        for inner in range(1, inner_iterations + 1):
            next_idx = select_neighbor(
                current_idx,
                neighborhoods,
                coords,
                energies,
                spatial_weight,
                energy_weight
            )

            if next_idx is None:
                pbar.update(1)
                continue

            next_energy = energies[next_idx]
            delta_E = next_energy - current_energy

            accepted = False
            if delta_E < 0:
                accepted = True
            elif T > 0:
                prob = math.exp(-delta_E / T)
                if random.random() < prob:
                    accepted = True

            if accepted:
                current_idx = next_idx
                current_energy = next_energy
                current_barcode = barcodes[next_idx]
                current_celltype = celltypes[next_idx]
                accept_count_total += 1
                
                # Store accepted state
                accepted_states.append((
                    coords[current_idx, 0], 
                    coords[current_idx, 1], 
                    current_barcode, 
                    current_energy
                ))
                
                # Store for trajectory
                trajectory_coords.append([coords[current_idx, 0], coords[current_idx, 1], current_barcode])
                trajectory_energies.append(current_energy)

            pbar.update(1)

        # Cool down
        T *= alpha
        if T <= T_min:
            break

    pbar.close()

    print(f"\nSimulation complete!")
    print(f"Total accepted states: {len(accepted_states)}")
    print(f"Final cell type: {current_celltype}")
    print(f"Final energy: {current_energy:.4f}")
    print(f"Acceptance rate: {accept_count_total/(outer_iterations*inner_iterations)*100:.2f}%")

    # Create smooth trajectory plot with energy gradient using average energies
    print(f"\nCreating smooth trajectory plot with energy gradient (averaging over {skip_points_for_plot} points per segment)...")
    
    # Extract coordinates and energies for plotting
    plot_coords = [[x, y, b] for x, y, b in trajectory_coords]
    plot_energies = trajectory_energies
    
    plot_smooth_trajectory_with_energy(
        glut_df, 
        plot_coords, 
        plot_energies,
        skip_points=skip_points_for_plot
    )

    return accepted_states, trajectory_coords, trajectory_energies


if __name__ == "__main__":
    # Run a single simulation and create smooth trajectory plot
    accepted_states, trajectory_coords, trajectory_energies = run_single_spatial_simulated_annealing(
        tsv_file='/Users/mozhganoroujlu/Desktop/MOZHGUN/cell_fate/hi_c/codes_figures/folders/normalized_contacts/generated_barcodes.tsv',
        k_neighbors=150,
        spatial_weight=0.6,
        energy_weight=0.4,
        outer_iterations=200,
        inner_iterations=50,
        skip_points_for_plot=5  # Connect every 5 accepted states and average their energies
    )