Preprocessing Pipeline

This section describes the preprocessing workflow used to generate single-cell and cluster-level Hi-C contact maps from sci-Hi-C data downloaded from *GEO (GSE253407)*.

Two samples were used in this study: *GSM8020244 (LC462)* and *GSM8020245 (LC716)*. Each sample consists of `.sc.pairs.gz` contact files along with `valid_bc.txt.gz` files containing valid cell barcodes.

The preprocessing begins by filtering barcodes and restricting contacts to chromosome 1. For each sample, valid barcodes are extracted and used to filter the corresponding pairs file using the script `0.1_filter_valid_bc.py`. This step ensures that only high-quality barcodes and their associated chromosome 1 contacts are retained. The filtering process is applied separately to each dataset.

After filtering, the processed files from both samples are merged into a single combined pairs file using `0.2_merging_files.py`. This merged dataset contains contacts from both *LC462* and *LC716* and serves as the unified input for downstream processing.

Next, the merged contacts are split into single-cell and cluster-specific contact files. This is performed using `0.3.1_sc_splitting.py` and `0.3.1_cluster_splitting.py`, which organize the contact data based on individual barcodes and predefined cell clusters.

In the final preprocessing step, the resulting .pairs files are converted into .cool format contact matrices using `0.4.1_sc_pairs_to_cool.py` for single-cell data and `0.4.1_cluster_pairs_to_cool.py` for cluster-level data, use [non_neuron.tsv](https://github.com/HaghverdiLab/scHiC-Energy-Landscape/blob/main/data/non_neuron.tsv) for cell type embedding. These .cool files represent the final output of the preprocessing pipeline.

The final outputs consist of two main directories: one containing single-cell .cool contact maps and another containing cluster-level .cool contact maps. These datasets are subsequently used for all downstream analyses.


