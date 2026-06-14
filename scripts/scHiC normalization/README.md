For normalization, we used the **BandNorm R package** (see the official GitHub repository for installation instructions and usage details: [BandNorm](https://github.com/sshen82/BandNorm)).

BandNorm requires input data in a long-format contact table, where each row represents an interaction between two genomic bins along with the corresponding contact count. The required format consists of five columns: `chr1, bin1, chr2, bin2, contacts` , with no header.

To generate this format, the .cool files (produced for each single cell using [0.4_sc_pairs_to_cool.py](https://github.com/Mozhgan0rujlu/scHiC-Energy-Landscape-/blob/main/scripts/preprocessing/0.4_sc_pairs_to_cool.py) were first converted into raw contact tables using the [1.0_cool_to_txt.py](https://github.com/HaghverdiLab/scHiC-Energy-Landscape/blob/main/scripts/scHiC%20normalization/1.0_cool_to_txt.py) script. This script extracts pairwise interactions from each .cool matrix and outputs files in the following format:

```text
chr1    6    chr1    6    258
chr1    6    chr1    7    10
chr1    6    chr1    8    6
chr1    6    chr1    9    2
```

Each file corresponds to a single cell and is named according to its barcode.

These raw contact files are then used as input for BandNorm. Normalization is performed by running `1.2_BandNorm.R`, which processes each cell independently and generates normalized contact matrices.

The final output is a directory containing normalized contact files for each single cell, stored as three-column text files in the following foramt:

```text
6	6	125.363345837984
6	7	3.91553584161911
6	8	2.40181410097148
6	9	0.687093425567686
```
Each file contains normalized interaction values between genomic bins and is used for downstream analyses.

