import argparse

import numpy as np
import pandas as pd

from pyTWMR import TWMR

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--beta',required=True, help='Path to the file containing the beta matrix of effect sizes of SNPs on gene expression')
    parser.add_argument('--gamma',required=True, help='Path to the file containing the gamma vector of SNP effect sizes on the trait')
    parser.add_argument('--ld',required=True, help='Path to the file containing the LD matrix with correlation coefficients between SNPs')

    parser.add_argument('--nEQTLs',required=True,type=int, help='Number of samples in the eQTL study used to estimate SNP effects on gene expression')
    parser.add_argument('--nGWAS',required=True,type=int, help='Number of samples in the GWAS used to estimate SNP effects on trait')

    parser.add_argument('--output',required=True, help='Path to output file')
    return parser.parse_args()


def main():
    args = parse_args()
    beta = pd.read_csv(args.beta, sep='\t', index_col=0)
    gamma = pd.read_csv(args.gamma, sep='\t', squeeze=True, index_col=0)
    ld = np.loadtxt(args.ld, delimiter='\t')

    if beta.values.ndim != 2:
        raise ValueError("Beta should be a matrix")
    if gamma.values.ndim != 1:
        raise ValueError("Gamma should be a single column")
    if ld.ndim != 2 or ld.shape[0] != ld.shape[1]:
        raise ValueError("LD should be a square matrix")
    
    alpha, se = TWMR(beta.values, gamma.values, args.nEQTLs, args.nGWAS, ld)

    output = pd.DataFrame(alpha, index=beta.index, columns=['alpha'])
    output['se'] = se
    output.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()
