#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 01:01:23 2025

@author: mozhganoroujlu
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np
import sklearn
import os
from scKTLD import edge2adj, callTLD, displayTLD

# Set a macOS-compatible font to avoid findfont warnings
matplotlib.rcParams['font.family'] = 'Arial'

# Pre-check and register colormaps only if not already registered
try:
    if 'hicheat' not in matplotlib.colormaps:
        hiccolors = ["lightyellow", "red"]
        my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('hicheat', hiccolors)
        matplotlib.colormaps.register(my_cmap)
    if 'featureheat' not in matplotlib.colormaps:
        featurecolors = ["blue", "white", "red"]
        my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('featureheat', featurecolors)
        matplotlib.colormaps.register(my_cmap)
    print("Colormaps checked/registered successfully")
except Exception as e:
    print(f"Error registering colormaps: {e}")
    pass

# Path to your edge list file
path_input = "/Users/mozhganoroujlu/Desktop/MOZHGUN/cell_fate/hi_c/codes_figures/folders/normalized_contacts/bandnorm/bandnorm_txt/AAACGAAAGACCGCAA.txt"

# Verify file exists
if not os.path.exists(path_input):
    raise FileNotFoundError(f"File not found: {path_input}. Ensure it's unzipped and check path.")

# Load edge list
try:
    graph_edge = np.loadtxt(path_input, dtype=np.float32)  # float32 to save memory
    print(f"Loaded edge list shape: {graph_edge.shape}")  # Expect (45301, 3)
    if graph_edge.shape[1] != 3:
        raise ValueError("Edge list must have 3 columns (source, target, weight).")
except Exception as e:
    print(f"Error loading file: {e}")
    raise

# Adjust bin indices (convert 1-based to 0-based if needed)
try:
    if graph_edge[:, :2].min() > 0:  # Check if 1-based
        print("Converting 1-based to 0-based indices")
        graph_edge[:, :2] -= 1
    print(f"Min index: {graph_edge[:, :2].min()}, Max index: {graph_edge[:, :2].max()}")
except Exception as e:
    print(f"Error adjusting indices: {e}")
    raise

# Convert edge list to adjacency matrix
try:
    graph_adj = edge2adj(graph_edge, chr='chr1', resolution=500000, reference='mm10')
    print(f"Converted to adjacency matrix shape: {graph_adj.shape}")  # Expect ~391x391
except Exception as e:
    print(f"Error converting to adjacency matrix: {e}")
    raise

# Run scKTLD analysis
try:
    boundary_spec = callTLD(graph_adj)  # Compute boundaries (TADs)
    # Visualize entire chromosome (0 to ~195M bp, 391 bins)
    displayTLD(graph_adj, boundary_spec, start=0, stop=195554896, brecon=True)  # bp
    # Alternative: Visualize segment (8start= start bin, stop= end bin)
    displayTLD(graph_adj, boundary_spec, start=0, stop=200, brecon=True)  # kb
    plt.show()  # Ensure plot displays in Spyder
except Exception as e:
    print(f"Error in scKTLD functions: {e}")
    raise