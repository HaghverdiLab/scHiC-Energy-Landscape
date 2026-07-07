#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 16:05:02 2026

@author: mozhganoroujlu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Improved & Accelerated Spatial Simulated Annealing
Optimized for Apple Silicon (M1/M2/M3/M4)
"""
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm
import random
import math
from openpyxl import load_workbook
# Optional: uncomment for extra speedup (needs numba installed: pip install numba)
# from numba import njit, prange

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
# Function to add detailed sheet to final summary
# ────────────────────────────────────────────────
def add_detailed_sheet_to_summary(final_summary_file: str, run_data: list, original_df: pd.DataFrame):
    """
    Add a detailed sheet to the final summary Excel file with columns:
    barcode, energy, cell type, state, total number counted as final state
    """
    # Create detailed dataframe
    detailed_data = []
    
    # Group by barcode to count occurrences as final state
    barcode_counts = {}
    for run in run_data:
        barcode = run['final_barcode']
        if barcode in barcode_counts:
            barcode_counts[barcode] += 1
        else:
            barcode_counts[barcode] = 1
    
    # Get unique final barcodes
    unique_barcodes = set([run['final_barcode'] for run in run_data])
    
    for barcode in unique_barcodes:
        # Find the run data for this barcode (use the first occurrence)
        run_info = next(run for run in run_data if run['final_barcode'] == barcode)
        
        # Get cell info from original dataframe
        cell_info = original_df[original_df['barcode'] == barcode].iloc[0] if not original_df[original_df['barcode'] == barcode].empty else None
        
        if cell_info is not None:
            row = {
                'barcode': barcode,
                'energy': cell_info['energy'],
                'cell type': cell_info['celltype'],
                'total number counted as final state': barcode_counts[barcode]
            }
            
            # Add state column if it exists in original data
            if 'state' in original_df.columns:
                row['state'] = cell_info['state']
            
            detailed_data.append(row)
    
    # Create dataframe
    detailed_df = pd.DataFrame(detailed_data)
    
    # Load existing workbook or create new one
    if os.path.exists(final_summary_file):
        with pd.ExcelWriter(final_summary_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
            detailed_df.to_excel(writer, sheet_name='Detailed_Final_States', index=False)
    else:
        with pd.ExcelWriter(final_summary_file, engine='openpyxl') as writer:
            detailed_df.to_excel(writer, sheet_name='Detailed_Final_States', index=False)


# ────────────────────────────────────────────────
# Function to create run statistics report
# ────────────────────────────────────────────────
def create_run_statistics_report(run_data: list, output_file: str, original_df: pd.DataFrame):
    """
    Create a separate Excel file with run statistics:
    initial state, final state, initial state energy, initial state cell type,
    initial state state, final state energy, final state cell type, final state state,
    total barcode visited, total barcode accepted, average acceptance rate
    """
    statistics_data = []
    
    for run in run_data:
        # Get initial cell info
        initial_cell_info = original_df[original_df['barcode'] == run['initial_barcode']].iloc[0] if not original_df[original_df['barcode'] == run['initial_barcode']].empty else None
        
        # Get final cell info
        final_cell_info = original_df[original_df['barcode'] == run['final_barcode']].iloc[0] if not original_df[original_df['barcode'] == run['final_barcode']].empty else None
        
        if initial_cell_info is not None and final_cell_info is not None:
            row = {
                'run_number': run['run_number'],
                'initial_state': run['initial_barcode'],
                'final_state': run['final_barcode'],
                'initial_state_energy': initial_cell_info['energy'],
                'initial_state_cell_type': initial_cell_info['celltype'],
                'final_state_energy': final_cell_info['energy'],
                'final_state_cell_type': final_cell_info['celltype'],
                'total_barcode_visited': run['total_barcodes_visited'],
                'total_barcode_accepted': run['total_barcodes_accepted'],
                'average_acceptance_rate': run['average_acceptance_rate']
            }
            
            # Add state columns if they exist in original data
            if 'state' in original_df.columns:
                row['initial_state_state'] = initial_cell_info['state']
                row['final_state_state'] = final_cell_info['state']
            
            statistics_data.append(row)
    
    # Create dataframe
    statistics_df = pd.DataFrame(statistics_data)
    
    # Save to Excel
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        statistics_df.to_excel(writer, sheet_name='Run_Statistics', index=False)
        
        # Add summary statistics
        summary_data = {
            'Metric': ['Total Runs', 'Average Acceptance Rate', 'Unique Final Barcodes'],
            'Value': [
                len(run_data),
                statistics_df['average_acceptance_rate'].mean(),
                statistics_df['final_state'].nunique()
            ]
        }
        pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)


# ────────────────────────────────────────────────
# Main SA function
# ────────────────────────────────────────────────
def run_spatial_simulated_annealing(
    tsv_file: str,
    output_excel: str = "/Users/mozhganoroujlu/Desktop/SA/score/results/modif/generated/all_opcs/150n_e4_all_runs.xlsx",
    output_rates: str = "/Users/mozhganoroujlu/Desktop/SA/score/results/modif/generated/all_opcs/150n_e4_acceptance_rates.xlsx",
    output_final: str = "/Users/mozhganoroujlu/Desktop/SA/score/results/modif/generated/all_opcs/150n_e4_final_summary.xlsx",
    output_run_stats: str = "/Users/mozhganoroujlu/Desktop/SA/score/results/modif/generated/all_opcs/150n_e4_run_statistics.xlsx",
    num_runs: int = 5,
    k_neighbors: int = 100,
    spatial_weight: float = 0.7,
    energy_weight: float = 0.3,
    T_start: float = 0.4,
    T_min: float = 1e-6,
    alpha: float = 0.98,
    outer_iterations: int = 300,
    inner_iterations: int = 70,
    plot_first_only: bool = True,
) -> tuple[pd.DataFrame, dict, np.ndarray]:
    """
    Runs spatial neighborhood-based simulated annealing from OPC cells.
    Optimized version for speed on Apple Silicon.
    """
    # ── Prepare output paths ──
    ensure_dir(output_excel)
    ensure_dir(output_rates)
    ensure_dir(output_final)
    ensure_dir(output_run_stats)

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

    # ── Find OPC starting points ──
    opc_mask = celltypes == 'OPC'
    opc_indices = np.where(opc_mask)[0]
    if len(opc_indices) < num_runs:
        raise ValueError(f"Not enough OPC cells ({len(opc_indices)}) for {num_runs} runs")

    # Sample unique starting points
    start_indices = np.random.choice(opc_indices, num_runs, replace=False)

    # ── Cooling schedule (shared) ──
    Ts = []
    T = T_start
    for _ in range(outer_iterations):
        Ts.append(T)
        T *= alpha
        if T <= T_min:
            break

    # ── Prepare Excel writers ──
    all_acceptance_rates = []
    final_results = []
    
    # Store detailed run data for the new output formats
    detailed_run_data = []

    with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
        for run_idx in tqdm(range(num_runs), desc="SA runs"):
            current_idx = int(start_indices[run_idx])
            current_energy = energies[current_idx]
            current_barcode = barcodes[current_idx]
            current_celltype = celltypes[current_idx]
            
            # Store initial state for run statistics
            initial_barcode = current_barcode
            initial_celltype = current_celltype
            initial_energy = current_energy

            results = []
            results.append({
                'T': None,
                'iteration': 'initial',
                'final_cell_energy': current_energy,
                'final_cell_type': current_celltype,
                'acceptance_state': None,
                'barcode': current_barcode
            })

            T = T_start
            accept_count_total = 0
            acceptance_rates_this_run = []
            
            # Track unique barcodes visited
            visited_barcodes = set([current_barcode])

            for outer in range(outer_iterations):
                accept_count = 0

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
                        # stay
                        results.append({
                            'T': T,
                            'iteration': inner,
                            'final_cell_energy': current_energy,
                            'final_cell_type': current_celltype,
                            'acceptance_state': 0,
                            'barcode': current_barcode
                        })
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

                    acceptance_state = 1 if accepted else 0
                    if accepted:
                        current_idx = next_idx
                        current_energy = next_energy
                        current_barcode = barcodes[next_idx]
                        current_celltype = celltypes[next_idx]
                        accept_count += 1
                        accept_count_total += 1
                        
                        # Add to visited barcodes
                        visited_barcodes.add(current_barcode)

                    results.append({
                        'T': T,
                        'iteration': inner,
                        'final_cell_energy': current_energy,
                        'final_cell_type': current_celltype,
                        'acceptance_state': acceptance_state,
                        'barcode': current_barcode
                    })

                # Record acceptance rate per temperature
                acceptance_rate = accept_count / inner_iterations
                acceptance_rates_this_run.append(acceptance_rate)

                # Cool down
                T *= alpha
                if T <= T_min:
                    break

            # Save this run
            run_df = pd.DataFrame(results)
            run_df.to_excel(writer, sheet_name=f'Run_{run_idx+1}', index=False)

            all_acceptance_rates.append(acceptance_rates_this_run)

            avg_accept = np.mean(acceptance_rates_this_run) if acceptance_rates_this_run else 0

            final_results.append({
                'run_number': run_idx + 1,
                'final_celltype': current_celltype,
                'energy': current_energy,
                'average_acceptance_rate': avg_accept,
                'barcode': current_barcode
            })
            
            # Store detailed run data for new output formats
            detailed_run_data.append({
                'run_number': run_idx + 1,
                'initial_barcode': initial_barcode,
                'initial_celltype': initial_celltype,
                'initial_energy': initial_energy,
                'final_barcode': current_barcode,
                'final_celltype': current_celltype,
                'final_energy': current_energy,
                'total_barcodes_visited': len(visited_barcodes),
                'total_barcodes_accepted': accept_count_total,
                'average_acceptance_rate': avg_accept
            })

            # Plot only first run if requested
            if plot_first_only and run_idx == 0:
                print("Plotting trajectory for run 1...")
                plot_spatial_trajectory(glut_df, run_df, run_number=1)

    # ── Save acceptance rates ──
    print("Saving acceptance rates...")
    rates_data = {'T': Ts}
    max_len = len(Ts)
    for i, rates in enumerate(all_acceptance_rates, 1):
        padded = rates + [np.nan] * (max_len - len(rates)) if len(rates) < max_len else rates[:max_len]
        rates_data[f'Run_{i}'] = padded

    pd.DataFrame(rates_data).to_excel(output_rates, index=False)

    # ── Save final summary ──
    print("Saving final summary...")
    df_final = pd.DataFrame(final_results)
    with pd.ExcelWriter(output_final, engine='openpyxl') as writer:
        df_final.to_excel(writer, sheet_name='Runs', index=False)
        df_final['final_celltype'].value_counts().reset_index().to_excel(
            writer, sheet_name='Summary', index=False
        )
    
    # ── Add detailed sheet to final summary ──
    print("Adding detailed sheet to final summary...")
    add_detailed_sheet_to_summary(output_final, detailed_run_data, glut_df)
    
    # ── Create run statistics report ──
    print("Creating run statistics report...")
    create_run_statistics_report(detailed_run_data, output_run_stats, glut_df)

    print(f"Done. Results in:\n  {output_excel}\n  {output_rates}\n  {output_final}\n  {output_run_stats}")

    return df_final, neighborhoods, coords


# ────────────────────────────────────────────────
# Plotting (unchanged logic, minor cleanup)
# ────────────────────────────────────────────────
def plot_spatial_trajectory(original_df: pd.DataFrame, trajectory_df: pd.DataFrame, run_number: int = 1):
    plt.figure(figsize=(14, 10))

    cell_types = original_df['celltype'].unique()
    colors = plt.cm.Set3(np.linspace(0, 1, len(cell_types)))

    for i, ct in enumerate(cell_types):
        mask = original_df['celltype'] == ct
        plt.scatter(
            original_df.loc[mask, 'mds_component_1'],
            original_df.loc[mask, 'mds_component_2'],
            c=[colors[i]], label=ct, alpha=0.3, s=20
        )

    # Trajectory path
    traj = trajectory_df[trajectory_df['iteration'] != 'initial']
    traj_coords = original_df.set_index('barcode').loc[traj['barcode'], ['mds_component_1', 'mds_component_2']].to_numpy()

    plt.plot(traj_coords[:, 0], traj_coords[:, 1], 'k-', alpha=0.5, linewidth=1, label='Path')

    # Points along trajectory
    for i, (x, y) in enumerate(traj_coords):
        ct = traj.iloc[i]['final_cell_type']
        color_idx = list(cell_types).index(ct)
        plt.scatter(x, y, c=[colors[color_idx]], s=30, alpha=0.7)

    # Start & End
    if len(traj_coords) > 0:
        plt.scatter(*traj_coords[0], marker='*', s=200, color='green', label='Start')
        plt.scatter(*traj_coords[-1], marker='*', s=200, color='red', label='End')

    plt.legend()
    end_type = traj.iloc[-1]['final_cell_type'] if not traj.empty else '?'
    plt.title(f'Spatial Trajectory - Run {run_number} - End: {end_type}')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Customize paths & parameters here
    run_spatial_simulated_annealing(
        tsv_file='/Users/mozhganoroujlu/Desktop/MOZHGUN/cell_fate/hi_c/codes_figures/folders/normalized_contacts/generated_barcodes.tsv',
        num_runs=1068,
        k_neighbors=150,           # ← lowered default (increase if needed)
        spatial_weight=0.6,
        energy_weight=0.4,
        outer_iterations=200,
        inner_iterations=50,
    )
