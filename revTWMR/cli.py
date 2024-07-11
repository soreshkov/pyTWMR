import argparse, os

import pandas as pd

from revTWMR import revTWMR


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--effect',required=True, help='Path to the file containing the beta matrix of standardized effect sizes of SNPs on gene expression')    
    parser.add_argument('--ld',required=False, help='Path to the file containing the LD matrix with correlation coefficients between SNPs')
    parser.add_argument('--sampleSize', required=True, help='Sample size for each gene')
    parser.add_argument('--hetThreshold',required=False, type=float, default=0.05, help='Threshold to exclude heterogenic QTLs. Default is 0.05')
    parser.add_argument('--pseudoInverse',action='store_true', help='Uses Moore-Penrose pseudo-inverse instead of basic inverse matrix operation if present.')
    parser.add_argument('--device',required=False, default='cpu', help='Whether to use GPU or CPU for computations. Available options are `cpu` and `cuda`. Default is `cpu`.')
    parser.add_argument('--output',required=True, help='Path to output file')
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.effect):
        print(f"Error: File '{args.effect}' does not exist.")
        exit(1)

    if not os.path.exists(args.sampleSize):
        print(f"Error: File '{args.sampleSize}' does not exist.")
        exit(1)

    effect = pd.read_csv(args.effect, sep='\t', index_col=0)
    sample_size = pd.read_csv(args.sampleSize, sep='\t', index_col=0).N.to_dict()
    gwasEffect = effect.BETA_GWAS.values
    pseudoInverse = args.pseudoInverse
    device = args.device
    
    nGwas = effect.N.values
    
    results = {}
    for gene in effect.drop(['BETA_GWAS', 'SE', 'N'], axis=1, errors='ignore').columns:    
        effectTbl = effect.drop(['BETA_GWAS', 'SE', 'N'], axis=1, errors='ignore')[gene].values
        nQtls = sample_size[gene]
        
        result = revTWMR(
            effectTbl, 
            gwasEffect, 
            qtlLabels=effect.index.values, 
            gwasSizes=nGwas, 
            qtlExpSize=nQtls,
            pValIterativeThreshold=args.hetThreshold, 
            pseudoInverse=pseudoInverse, device=device)
        
        results[gene] = {
            'Alpha Original': result.Alpha,
            'SE Original': result.Se,
            'P Value Original': result.Pval,
            'N Original': result.N,
            'P heterogeneity Original': result.HetP,
            'Alpha': result.AlphaIterative,
            'SE': result.SeIterative,
            'P': result.PvalIterative,
            'P heterogeneity': result.HetPIterative,
            'N': len(result.rsname)
        }

    pd.DataFrame.from_dict(results, orient='index').to_csv(args.output, sep='\t', index_label='Gene')
    

if __name__ == '__main__':
    main()
