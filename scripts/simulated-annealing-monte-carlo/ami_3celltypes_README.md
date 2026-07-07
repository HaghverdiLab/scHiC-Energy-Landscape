### Step 1: Generate MDS Coordinates from TAD Similarity

Run:

```text
scripts/ami_3celltypes.py
```

#### Prerequisites

Before running this script, you must first detect TAD boundaries for the non-neuronal cells (ASC, OGC, and OPC).

Run:

```text
scripts/TAD-detection/2.1_scKTLD_TADs_detection.py
```

This script generates the TAD boundary files required by `ami_3celltypes.py`.

Next, update the input directory in `ami_3celltypes.py` so that it points to the directory containing the generated TAD boundaries:

```python
tad_directory = "/path/to/TADs_non_neuron/"
```

The directory should contain the TAD boundary files for the non-neuronal barcodes (ASC, OGC, and OPC).

#### Similarity Matrices

Running `scripts/TAD-detection/2.1_scKTLD_TADs_detection.py` generates two Excel files:

* `simi_ami.xlsx`
* `simi_moc.xlsx`

Place both files in the directory specified by:

```python
output_dir = "/path/to/similarity_results/"
```

The script expects to find them using:

```python
moc_file = os.path.join(output_dir, "simi_moc.xlsx")
ami_file = os.path.join(output_dir, "simi_ami.xlsx")
```

Although both similarity matrices are produced, **only the AMI similarity matrix (****`simi_ami.xlsx`****) was used for the downstream analyses presented in the paper**. The MOC matrix is retained because the code supports both similarity measures.

#### Output

The script applies **Multidimensional Scaling (MDS)** to the AMI similarity matrix and saves the resulting two-dimensional barcode coordinates as:

```text
ami_mds_coordinates_ASC_OGC_OPC.tsv
```

This file is used as the input for the subsequent simulated annealing pipeline.
