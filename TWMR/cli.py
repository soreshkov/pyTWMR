import argparse
import pathlib
import logging

import numpy as np
import pandas as pd

from TWMR import TWMR, QCorrectedTWMR

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--effect',required=True, type=str, help='Path to the file containing the beta matrix of standardized effect sizes of SNPs on gene expression')
    parser.add_argument('--ld',required=True, type=str, help='Path to the file containing the LD matrix with correlation coefficients between SNPs')

    parser.add_argument('--nEQTLs',required=True,type=int, help='Number of samples in the eQTL study used to estimate SNP effects on gene expression')
    parser.add_argument('--nGWAS',required=True,type=int, help='Number of samples in the GWAS used to estimate SNP effects on trait')
    
    parser.add_argument('--removeOutliers', action='store_true', help='Iteratively removes heterogenious outliers if this option is present.')
    
    parser.add_argument('--pseudoInverse',action='store_true', help='Uses Moore-Penrose pseudo-inverse instead of basic inverse matrix operation if present.')
    parser.add_argument('--device',required=False, default='cpu', help='Whether to use GPU or CPU for computations. Available options are `cpu` and `cuda`. Default is `cpu`.')
    parser.add_argument('--output', type=str,required=True, help='Path to output file')
    return parser.parse_args()



def main():
    args = parse_args()
    effect_path = args.effect
    ld_path = args.ld

    if not pathlib.Path(effect_path).exists():
        logging.error(f"File {effect_path} not found")
        exit(1)
    if ld_path and not pathlib.Path(ld_path).exists():
        logging.error(f"File {ld_path} not found")
        exit(1)

    effect_df = pd.read_csv(effect_path, sep='\t')
    ld_df = pd.read_csv(ld_path, sep='\t', header=None).to_numpy().astype(np.float32)


    if 'GWAS' not in effect_df.columns:
        logging.error(f"Column 'GWAS' not found in {effect_path}")
        exit(1)
    if 'SNPS' not in effect_df.columns:
        logging.error(f"Column 'SNPS' not found in {effect_path}")
        exit(1)
    else:
        effect_df.set_index('SNPS', inplace=True)

    if (len(ld_df.shape) != 2) and (ld_df.shape[0] != ld_df.shape[1]):
        logging.error(f"LD matrix is not square")
        exit(1)

    nEQTLs = args.nEQTLs
    nGWAS = args.nGWAS

    if args.removeOutliers:
        result,_ = QCorrectedTWMR(
            beta=effect_df.drop('GWAS', axis=1).to_numpy(), 
            gamma=effect_df['GWAS'].to_numpy(), 
            nEQTLs=nEQTLs, 
            NGwas=nGWAS, 
            rsnames= effect_df.index.to_list(),
            ldMatrix=ld_df, 
            pseudoInverse=args.pseudoInverse, 
            device=args.device
        )
    else:
        result = TWMR(
                beta=effect_df.drop('GWAS', axis=1).to_numpy(), 
                gamma=effect_df['GWAS'].to_numpy(), 
                nEQTLs=nEQTLs, 
                NGwas=nGWAS, 
                ldMatrix=ld_df, 
                pseudoInverse=args.pseudoInverse, 
                device=args.device)

#save results to output file
    
    result_df = pd.DataFrame()
    result_df['Gene'] = effect_df.columns.drop('GWAS')
    result_df['Alpha'] = result.Alpha
    result_df['Standard error'] = result.Se
    result_df['P-value'] = result.Pval
    result_df['Heterogenity P-value'] = result.HetP

    result_df.to_csv(args.output, index=False, sep='\t')


if __name__ == '__main__':
    main()

