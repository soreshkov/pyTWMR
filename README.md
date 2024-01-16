# pyTWMR

Python implementation of [Transcriptome-Wide Mendelian Randomization(TWMR)](https://github.com/eleporcu/TWMR) and [revTWMR](https://github.com/eleporcu/revTWMR) methods by E. Porcu et al. with GPU computation support.

These two Mendelian Randomization methods are modified IVW MR, which account for Linkage Disequilibrium between QTLs(or, with account for correlation between instrumental variables, outside of genetics field).

![Data layout](data_layout.svg)

Methods estimate causal effect of QTLs and genes expression on outcome trait in case of TWMR and causal effect of trait on gene expression in case of revTWMR.

Methods require estimated effect sizes of QTLs on Gene expression(matrix) and on outcome(vector). Additionally, number of QTLs effecting gene expression and total number of samples in QTL-trait GWAS study required for accurate normalization.

- [Installation](#installation)
- [Usage](#usage)

## Installation

With pip:

1. `git clone` this repo into your local folder
2. `pip3 install .` from that folder

OR

`pip3 install git+https://github.com/soreshkov/pyTWMR.git`

## Usage

### In Python script

Core of methods are functions `TWMR` and `revTWMR`, implemented in corresponding files.

Import desired method:

```python
from pyRevTWMR import revTWMR
from pyTWMR import TWMR
```

#### TWMR

To use `TWMR`, just call function `TWMR`, which expects following parameters:

- `beta` - numpy 2D array, with effect sizes of QTLs on gene expression, with columns corresponding to genes and rows - to QTLs.
- `gamma` - numpy vector(1D array) with effect sizes of QTLs on outcome trait. Size of that vector is expected to be equal to number of rows in `beta`.
- `nEQTLs` - number of QTLs in study
- `nGWAS` - number of samples in QTL-to-trait GWAS
- `ldMatrix` - square numpy array with LD between QTLs. If not provided, identity matrix is used.
- `pseudoInverse` - True/False, True to use Moore-Penrose pseudo-inverse instead of inverse matrix
- `device` - If there is GPU with Cuda API available, you can pass name of GPU device(as enumerated by PyTorch framework) to utilize it instead of `cpu`(default option)

Output:

- `alpha` - array of causal effect values for each gene
- `se` - array of standard error estimations for each estimated causal effect

---

#### revTWMR

To use `revTWMR`, just call function `revTWMR`, which expects following parameters:

- `beta` - numpy vector(1D array) with effect sizes of QTLs on outcome trait. Size of that vector is expected to be equal to number of rows in `gamma`
- `gamma` - numpy 2D array, with effect sizes of QTLs on gene expression, with columns corresponding to genes and rows - to QTLs.
- `nGWAS` - number of samples in QTL-to-trait GWAS
- `ldMatrix` - square numpy array with LD between QTLs. If not provided, identity matrix is used.
- `pseudoInverse` - True/False, True to use Moore-Penrose pseudo-inverse instead of inverse matrix
- `device` - If there is GPU with Cuda API available, you can pass name of GPU device(as enumerated by PyTorch framework) to utilize it instead of `cpu`(default option)

Output:

- `alpha` - array of causal effect values for each gene
- `se` - array of standard error estimations for each estimated causal effect

See [demo.ipynb](https://github.com/soreshkov/pyTWMR/blob/master/demo.ipynb) for usage example.

### As console app

After installation it is possible to use package as console application.

For TWMR command would be:

`TWMR --beta BETA --gamma GAMMA --ld LD --nEQTLs NEQTLS --nGWAS NGWAS --output OUTPUT`

Where:

- `beta` - path to TSV file containing the beta matrix of effect sizes of SNPs on gene expression

- `gamma` - path to TSV file containing the gamma vector of SNP effect sizes on the trait

- `ld`- path to the file containing the LD matrix with correlation coefficients between SNPs
- `nEQTLs` -  number of samples in the eQTL study used to estimate SNP effects on gene expression
- `nGWAS` - number of samples in the GWAS used to estimate SNP effects on trait
- `output` Path to output file

For `beta` matrix we expect tab-separated table, with first column containing QTL labels, and rest of columns containing QTL effects on genes with gene name as column header.

For `gamma` GWAS vector we expect TSV files with first column containing QTLs and second column containing effect size of QTL on outcome trait.

For `ld` matrix we expect TSV file containing square matrix of correlations, without any column or row headers. Order of QTLs in matrix should match order of QTLs in `beta` and `gamma`.

Output TSV file contains of 3 columns, which contain list of genes, its estimated causal effect and standard error of this estimation.

---

Command for RevTWMR:

`RevTWMR --beta BETA --gamma GAMMA --ld LD --nGWAS NGWAS --output OUTPUT`

- `beta` - path to the file containing the beta vector
- `gamma` -    path to the file containing the gamma matrix
- `ld` -   path to the file containing the LD matrix with correlation coefficients between SNPs
- `nGWAS` -     number of samples in the GWAS used to estimate SNP effects on trait
- `output` -  path to output file

For `gamma` matrix we expect tab-separated table, with first column containing QTL labels, and rest of columns containing QTL effects on genes with gene name as column header.

For `beta` GWAS vector we expect TSV files with first column containing QTLs and second column containing effect size of QTL on outcome trait.

For `ld` matrix we expect TSV file containing square matrix of correlations, without any column or row headers. Order of QTLs in matrix should match order of QTLs in `beta` and `gamma`.

Output TSV file contains of 3 columns, which contain list of genes, its estimated causal effect and standard error of this estimation.