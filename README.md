# pyTWMR

Python implementation of [Transcriptome-Wide Mendelian Randomization(TWMR)](https://github.com/eleporcu/TWMR) and [revTWMR](https://github.com/eleporcu/revTWMR) methods by E. Porcu et al. with GPU computation support.

These two Mendelian Randomization methods are modified IVW MR, which account for Linkage Disequilibrium between QTLs(or, with account for correlation between instrumental variables, outside of genetics field).

![Data layout](data_layout.svg)

Methods estimate causal effect of QTLs and genes expression on outcome trait in case of TWMR and causal effect of trait on gene expression in case of revTWMR.

Methods require estimated effect sizes of QTLs on Gene expression(matrix) and on outcome(vector). Additionally, number of QTLs effecting gene expression and total number of samples in QTL-trait GWAS study required for accurate normalization.


- [Installation](#installation)
- [Usage](#usage)

## Installation

<details> <summary> As Python application/module </summary>

1. `git clone` this repo into your local folder
2. `pip3 install .` from that folder

</details>

## Usage

### In Python script.
Core of methods are functions `TWMR` and `revTWMR`, implemented in corresponding files. 

Import desired method:
```
from pyRevTWMR import revTWMR
from pyTWMR import TWMR
```
#### TWMR 
To use `TWMR`, just call function `TWMR`, which expects following parameters:
* `beta` - numpy 2D array, with effect sizes of QTLs on gene expression, with columns corresponding to genes and rows - to QTLs.
* `gamma` - numpy vector(1D array) with effect sizes of QTLs on outcome trait. Size of that vector is expected to be equal to number of rows in `beta`.
* `nEQTLs` - number of QTLs in study
* `nGWAS` - number of samples in QTL-to-trait GWAS
* `ldMatrix` - square numpy array with LD between QTLs. If not provided, identity matrix is used.
* `pseudoInverse` - True/False, True to use Moore-Penrose pseudo-inverse instead of inverse matrix
* `device` - If there is GPU with Cuda API available, you can pass name of GPU device(as enumerated by PyTorch framework) to utilize it instead of `cpu`(default option)

Output:
* `alpha` - array of causal effect values for each gene
* `se` - array of standard error estimations for each estimated causal effect

---

#### revTWMR
To use `revTWMR`, just call function `revTWMR`, which expects following parameters:
* `beta` - numpy vector(1D array) with effect sizes of QTLs on outcome trait. Size of that vector is expected to be equal to number of rows in `gamma`
* `gamma` - numpy 2D array, with effect sizes of QTLs on gene expression, with columns corresponding to genes and rows - to QTLs.
* `nGWAS` - number of samples in QTL-to-trait GWAS
* `ldMatrix` - square numpy array with LD between QTLs. If not provided, identity matrix is used.
* `pseudoInverse` - True/False, True to use Moore-Penrose pseudo-inverse instead of inverse matrix
* `device` - If there is GPU with Cuda API available, you can pass name of GPU device(as enumerated by PyTorch framework) to utilize it instead of `cpu`(default option)

Output:
* `alpha` - array of causal effect values for each gene
* `se` - array of standard error estimations for each estimated causal effect