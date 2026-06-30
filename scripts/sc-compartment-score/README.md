# Usage

## Step 1. Generate CpG ratios

Before running the compartment calling pipeline, you must first generate the CpG ratio for each genomic bin using the `higashi_cpg.py` script.

To do this, download the following files:

- **`mm10.fa.gz`** from the UCSC Genome Browser:
  https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/

- **`chrom_sizes.txt`** from the repository's `data/` directory.

After downloading these files, update the input paths in `higashi_cpg.py` and run:

```bash
def main():
    # Hardcoded file paths
    fasta_path = "mm10.fa"
    chrom_size_path = "chrom_sizes.txt"
    output_txt_path = "/Users/mozhganoroujlu/Desktop/MOZHGUN/cell_fate/hi_c/temp_files/chrom1_cpg_ratios.txt"
    resolution = 500000  # 500 kb resolution
```


```bash
python higashi_cpg.py
```

This will generate a tab-separated file containing the CpG ratio for every genomic bin, which is required for compartment orientation.

---

## Step 2. Prepare the cell type annotation

Download **`non_neuron.tsv`** from the repository's `data/` directory.

This is a two-column tab-separated file containing the barcode of each cell and its corresponding cell type:

```text
barcode    celltype
```

Example:

```text
AAACGAAAGACCGCAA    OPC
AAACGAAAGTTCTGTA    ASC
```

---

## Step 3. Run compartment calling

Update the input and output paths in `higashi_compartment.py`:

- directory containing the single-cell Hi-C contact matrices (`.txt` files),
- CpG ratio file generated in **Step 1**,
- `non_neuron.tsv` cell type annotation,
- output directories.

Then run:

```bash
python higashi_compartment_score.py
```

The script will:

- group cells by cell type,
- construct a bulk contact matrix for each cell type,
- compute bulk PCA following the Higashi compartment-calling strategy,
- orient PC1 using CpG enrichment,
- project each single cell into the bulk PCA space,
- quantile-normalize compartment scores to the range **[-2.5, 2.5]**,
- save compartment scores for every cell,
- generate quality-control figures and summary plots.
