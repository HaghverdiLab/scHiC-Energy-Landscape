#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updates:
    computes Ising energy
- 2-column compartment score format: bin_id  score
- Energy distribution figure colored by cell type (from TSV metadata)
- Violin plot figure: OPC, ASC, OGC side by side

"""
import os
import csv
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde


def save_fig(fig, base_path):
    """Save figure as both PNG (screen) and PDF (paper)."""
    for ext in ['.png', '.pdf']:
        path = base_path + ext
        fig.savefig(path, dpi=150, bbox_inches='tight')
        print(f"Saved: {path}")


def plot_energy_by_celltype(results, barcode_to_celltype, output_base):
    """
    Overlapping histograms + KDE curves per cell type, with rug plot.
    No grid. Tighter KDE bandwidth to resolve multiple peaks.
    """
    typed = []
    for barcode, energy in results:
        ct = barcode_to_celltype.get(barcode, "Unknown")
        typed.append((barcode, energy, ct))

    # Draw order controls overlap stacking — items drawn later sit on top.
    # OGC and OPC are drawn first (pushed backward), ASC drawn last (brought forward).
    draw_order = ['OGC', 'OPC', 'ASC']
    present = set(t[2] for t in typed)
    celltypes = [ct for ct in draw_order if ct in present]
    celltypes += sorted(ct for ct in present if ct not in draw_order)  # any extras, just in case

    color_map = {
        'OPC': '#4C72B0',
        'ASC': '#55A868',
        'OGC': '#C44E52',
    }

    zorder_map = {'OGC': 1, 'OPC': 2, 'ASC': 3}

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_facecolor('white')
    ax.grid(False)                          # no grid

    all_energies = [e for _, e, _ in typed]
    global_min = min(all_energies)
    global_max = max(all_energies)
    bins = np.linspace(global_min, global_max, 35)

    lowest = min(typed, key=lambda x: x[1])
    lowest_barcode, lowest_energy, lowest_ct = lowest

    for ct in celltypes:
        ct_energies = np.array([e for _, e, c in typed if c == ct])
        color = color_map[ct]
        z = zorder_map.get(ct, 1)

        ax.hist(ct_energies, bins=bins, alpha=0.45, color=color,
                edgecolor='white', linewidth=0.5, label=ct, zorder=z)

        if len(ct_energies) >= 2:
            kde = gaussian_kde(ct_energies, bw_method=0.25)
            x_range = np.linspace(global_min, global_max, 500)
            bin_width = bins[1] - bins[0]
            ax.plot(x_range, kde(x_range) * len(ct_energies) * bin_width,
                    color=color, linewidth=2.0, zorder=z + 0.5)

        ax.plot(ct_energies, np.full_like(ct_energies, -2.5),
                '|', color=color, alpha=0.6, markersize=6, zorder=z)

    short_bc = lowest_barcode[:16] + ('...' if len(lowest_barcode) > 16 else '')
    annotation = (
        f"Lowest energy cell:\n"
        f"Barcode: {short_bc}\n"
        f"Type: {lowest_ct}\n"
        f"Energy: {lowest_energy:.2f}"
    )
    ax.text(
        0.02, 0.97, annotation,
        transform=ax.transAxes,
        va='top', ha='left', fontsize=9,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFF3CD',
                  edgecolor='#CCAA44', alpha=0.9)
    )

    ax.set_xlabel('Energy', fontsize=13)
    ax.set_ylabel('Number of Cells', fontsize=13)
    ax.set_title('Ising Energy Distribution by Cell Type',
                 fontsize=14, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(bottom=-5)
    ax.legend(title='Cell Type', frameon=True, fontsize=10,
              title_fontsize=10, loc='upper right')

    plt.tight_layout()
    save_fig(fig, output_base)
    plt.close()


def plot_violin_by_celltype(results, barcode_to_celltype, output_base):
    """
    Three-panel violin plot: OPC → ASC → OGC.
    No grid. Tighter KDE bandwidth on violin bodies.
    """
    CELLTYPE_ORDER = ['OPC', 'ASC', 'OGC']
    COLOR_MAP = {
        'OPC': '#4C72B0',
        'ASC': '#55A868',
        'OGC': '#C44E52',
    }

    groups = {ct: [] for ct in CELLTYPE_ORDER}
    for barcode, energy in results:
        ct = barcode_to_celltype.get(barcode, None)
        if ct in groups:
            groups[ct].append(energy)

    fig, axes = plt.subplots(1, 3, figsize=(12, 6), sharey=True)
    fig.suptitle('Ising Energy Distribution per Cell Type',
                 fontsize=15, fontweight='bold', y=1.01)

    for ax, ct in zip(axes, CELLTYPE_ORDER):
        ax.set_facecolor('white')
        ax.grid(False)                      # no grid

        data = np.array(groups[ct])
        color = COLOR_MAP[ct]

        if len(data) >= 2:
            # Standard violin — no KDE, uses raw data distribution
            parts = ax.violinplot(data, positions=[0], showmedians=False,
                                  showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(color)
                pc.set_edgecolor(color)
                pc.set_alpha(0.65)
                pc.set_linewidth(1.2)

            # Box plot elements
            q1, median, q3 = np.percentile(data, [25, 50, 75])
            iqr = q3 - q1
            whisker_lo = max(data.min(), q1 - 1.5 * iqr)
            whisker_hi = min(data.max(), q3 + 1.5 * iqr)

            ax.add_patch(plt.Rectangle(
                (-0.08, q1), 0.16, iqr,
                facecolor='none', edgecolor=color, linewidth=1.5, zorder=3
            ))
            ax.hlines(median, -0.08, 0.08, color=color, linewidth=2.5, zorder=4)
            ax.vlines(0, whisker_lo, q1, color=color, linewidth=1.5, zorder=3)
            ax.vlines(0, q3, whisker_hi, color=color, linewidth=1.5, zorder=3)
            ax.hlines(whisker_lo, -0.04, 0.04, color=color, linewidth=1.5, zorder=3)
            ax.hlines(whisker_hi, -0.04, 0.04, color=color, linewidth=1.5, zorder=3)

            # Jittered points
            jitter = np.random.default_rng(42).uniform(-0.06, 0.06, len(data))
            ax.scatter(jitter, data, color=color, alpha=0.35,
                       s=8, zorder=2, edgecolors='none')

            ax.text(0.97, 0.02,
                    f'n = {len(data)}\nmedian = {median:.2f}',
                    transform=ax.transAxes, ha='right', va='bottom',
                    fontsize=9, color='#444444')
        else:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=11, color='gray')

        ax.set_title(ct, fontsize=13, fontweight='bold', color=color)
        ax.set_xticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

    # Only the leftmost panel keeps the Y-axis spine and label
    axes[0].spines['left'].set_visible(True)
    axes[0].set_ylabel('Energy', fontsize=12)
    plt.tight_layout()
    save_fig(fig, output_base)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Compute compartment score sums from paired files.")
    parser.add_argument(
        "contacts_folder",
        default="/filtered_tads_contacts/",
        nargs="?",
        help="Path to the folder containing contact files (*_filtered.txt)"
    )
    parser.add_argument(
        "scores_folder",
        default="/output_compartments/new/new/quantile_scores/",
        nargs="?",
        help="Path to the folder containing score files (e.g., *.txt)"
    )
    parser.add_argument(
        "output_csv",
        default="/output_compartments/new/new/new_energy.csv",
        nargs="?",
        help="Path to the output CSV file (e.g., output.csv)"
    )
    parser.add_argument(
        "metadata_tsv",
        default="/non_neuron.tsv", # data directory
        nargs="?",
        help="(Optional) Path to TSV file with columns: barcode  celltype"
    )
    args = parser.parse_args()

    output_dir = os.path.dirname(args.output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    barcode_to_celltype = {}
    if args.metadata_tsv and os.path.exists(args.metadata_tsv):
        with open(args.metadata_tsv, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader, None)
            for row in reader:
                if len(row) >= 2:
                    barcode_to_celltype[row[0].strip()] = row[1].strip()
        print(f"Loaded {len(barcode_to_celltype)} barcode-celltype mappings.")
    else:
        print("No metadata TSV provided — all cells labeled 'Unknown'.")

    results = []
    missing_scores = []

    for filename in os.listdir(args.contacts_folder):
        if filename.endswith("_filtered.txt"):
            barcode = filename.replace("_filtered.txt", "")
            contact_path = os.path.join(args.contacts_folder, filename)
            score_path = os.path.join(args.scores_folder, barcode + ".txt")

            if not os.path.exists(score_path):
                missing_scores.append(barcode)
                continue

            scores = {}
            try:
                with open(score_path, 'r') as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) == 2:
                            bin_id, score = parts
                            scores[int(bin_id)] = float(score)
            except Exception as e:
                print(f"Error reading score file {score_path}: {e}")
                continue

            sum_prod = 0.0
            try:
                with open(contact_path, 'r') as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) == 3:
                            i, j, c = int(parts[0]), int(parts[1]), float(parts[2])
                            if i in scores and j in scores:
                                sum_prod += scores[i] * scores[j] * c
            except Exception as e:
                print(f"Error reading contact file {contact_path}: {e}")
                continue

            results.append((barcode, -sum_prod))

    try:
        with open(args.output_csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['barcode', 'Energy'])
            for row in results:
                writer.writerow(row)
        print(f"Output written to {args.output_csv}")
    except Exception as e:
        print(f"Error writing to {args.output_csv}: {e}")

    if len(results) >= 2:
        base = os.path.splitext(args.output_csv)[0]
        try:
            plot_energy_by_celltype(
                results, barcode_to_celltype,
                base + "_energy_by_celltype"      # no extension — save_fig adds .png/.pdf
            )
        except Exception as e:
            print(f"Warning: histogram figure failed: {e}")

        try:
            plot_violin_by_celltype(
                results, barcode_to_celltype,
                base + "_energy_violin"
            )
        except Exception as e:
            print(f"Warning: violin figure failed: {e}")
    else:
        print("Not enough results to plot (need at least 2 cells).")

    if missing_scores:
        print(f"Missing score files for {len(missing_scores)} barcodes: {', '.join(missing_scores)}")


if __name__ == "__main__":
    main()
