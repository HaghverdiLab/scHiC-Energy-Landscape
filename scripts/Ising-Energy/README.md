### Computing Ising Energy

To calculate the Ising energy used in this study, follow the pipeline below.

### Step 1: Generate Normalized Contact Matrices

First, generate normalized single-cell contact matrices.

Follow the contact normalization pipeline provided in:

```text
    /scripts/scHiC normalization/
```

The output of this step will be the normalized contact matrices required for downstream analysis.

### Step 2: Detect TADs in Single Cells

Next, obtain TAD loci for each single cell.

Follow the scKTLD TAD detection pipeline provided in:

```text
    /scripts/TAD-detection/2.1_scKTLD_TADs_detection.py
```

This step generates the detected TAD regions for each individual cell.

### Step 3: Filter Intra-TAD Contacts

Filter contacts that occur within detected TAD regions by running:

```bash
    scripts/Ising-Energy/filter_inTADs_contacts.py
```

This step extracts intra-TAD contacts that are used for the Ising energy calculation.

### Step 4: Calculate Single-Cell Compartment Scores

Generate single-cell compartment scores by following the pipeline in:

```text
    /scripts/sc-compartment-score/
```

The resulting compartment scores are required as input for the Ising energy model.

### Step 5: Calculate Ising Energy

Run:

```bash
python compute_energy.py
```

Before running the script, update the input paths inside `compute_energy.py` to point to your own generated files from the previous steps.

The script will calculate the Ising energy values used in the study.

---

### Optional: Figure 3A — Intra-TAD Contact Proportion

To calculate and visualize the proportion of intra-TAD contacts shown in **Figure 3A**, run:

```bash
python intraTAD_contacts.py
```

This analysis requires the outputs from the previous steps:

1. Normalized contact matrices
2. Detected TAD loci
3. Filtered intra-TAD contacts

Follow Steps 1–3 before running this script.

