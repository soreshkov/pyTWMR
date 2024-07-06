from collections import namedtuple
import warnings

import numpy as np
from typing import Union, List
import scipy

import torch


revTWMRResult = namedtuple('revTWMRResult', ['Alpha', 'Se', 'Pval', 'N', 'HetP', 'AlphaIterative', 'SeIterative', 'PvalIterative', 'HetPIterative', 'rsname'])

def CohrainQ(
        alpha: float, 
        beta: Union[np.ndarray, torch.Tensor], 
        gamma: Union[np.ndarray, torch.Tensor], 
        se: float, nQTL: float, nGWAS: float, device='cpu'):
    """
    Calculate the Cohrain Q test statistic and p-value.

    Parameters
    ----------
    alpha : float
        The alpha value.
    beta : Union[np.ndarray, torch.Tensor]
        The beta value.
    gamma : Union[np.ndarray, torch.Tensor]
        The gamma value.
    se : float
        The standard error value.
    nQTL : float
        The number of QTL.
    nGWAS : float
        The number of GWAS.
    device : str
        Device to perform coputations on('cpu' or 'gpu')
    
    Returns
    -------
    phet : float
        The p-value.
    d : np.ndarray
        Test statistics.
    """
    if isinstance(beta, np.ndarray):
        beta = torch.from_numpy(beta).to(device)
    else:
        beta= beta.to(device)

    if isinstance(gamma, np.ndarray):
        gamma = torch.from_numpy(gamma).to(device)
    else:
        gamma = gamma.to(device)

    d = gamma - alpha * beta
    var_d = 1 / nQTL + se * se * beta * beta + alpha * alpha * 1 / nGWAS + (1 / nGWAS) * se * se
    var_d = torch.diag(var_d).to(device)
    z = d[:, None].T @ torch.linalg.inv(var_d) @ d
    N = len(d)
    
    d = torch.abs(d)
    if d.device.type != 'cpu':
        d = d.cpu()

    phet = 1 - scipy.stats.chi2.cdf(float(z), N - 1)

    return phet, d.numpy()


def Alpha(beta: torch.Tensor, gamma: torch.Tensor, D: torch.Tensor,pseudoInverse = False, device='cpu'):
    beta = beta.to(device)
    gamma = gamma.to(device)
    D = D.to(device)
    inverse = torch.linalg.pinv if pseudoInverse else torch.linalg.inv
    S = beta @ inverse(D) @ beta
    alpha = beta @ inverse(D) @ gamma / S
    return alpha.cpu()


def SE_iterative(beta: torch.Tensor, gamma: torch.Tensor, ngwas, nqtl : float, pseudoInverse = False, device='cpu'):
    
    inverse = torch.linalg.pinv if pseudoInverse else torch.linalg.inv

    D = torch.eye(len(beta), dtype=torch.float64).to(device)
    Dinv = inverse(D)
    GCGinv = 1/ (beta @  inverse(D) @ beta)
    df_dg = beta @ Dinv *GCGinv
    p1 = torch.kron(GCGinv, gamma @ Dinv @ (torch.eye(len(beta)).to(device) - (((beta *GCGinv)[:,None]) @ beta[:,None].T) @ Dinv))
    p2 = torch.kron(-gamma @ Dinv @ beta * GCGinv, GCGinv * beta @ Dinv)
    df_DG = p1 + p2
    alp = GCGinv * (beta @ Dinv @ gamma)
    Vgam = torch.maximum(torch.from_numpy(np.array(1/nqtl)), torch.var(gamma - beta *alp) )  * torch.ones_like(gamma)
    SEs = torch.hstack([torch.sqrt(1/ngwas ), torch.sqrt(Vgam)]).to(device)
    J = torch.hstack([df_DG, df_dg ]).to(device)
    R = torch.eye(2).to(device)
    Sigma = (SEs[:, None] @ SEs[:, None].T )   * torch.kron(D, R)
    return torch.sqrt(J @ Sigma @ J).cpu()


def revTWMR(
        qtlExposureEffects :Union[np.ndarray, torch.Tensor], 
        qtlTraitGWASEffects:Union[np.ndarray, torch.Tensor], 
        qtlLabels:Union[np.ndarray, List[str]], 
        gwasSizes:Union[np.ndarray, torch.Tensor], 
        qtlExpSize:float, 
        pValIterativeThreshold:float = 0.05,
        pseudoInverse = False, 
        device='cpu'):
    """
    revTWMR method for Mendelian Randomization anylysis. 
    Analyses the causal effect of trait(outcome) on gene expression(exposure) derived from GWAS and eQTL data.
    
    Parameters
    ----------
        qtlExposureEffect : Union[np.ndarray, torch.Tensor]
            Vector of standardized effect sizes trans-QTLs effect on gene expression, shape (QTLs,) 
        qtlTraitGWASEffects :  Union[np.ndarray, torch.Tensor]
            Vector of standardized effect sizes of independent SNPs on trait, shape (QTLs,)
        qtlLabels : Union[np.ndarray, List[str]]
            QTL name labels, shape (QTLs,)
        gwasSizes : Union[np.ndarray, torch.Tensor]
            GWAS sample size, shape (QTLs,)
        qtlExpSize : float
            QTL Exposure sample size, number
        pValIterativeThreshold : float, optional
            Threshold to exclude pleyotropic QTLs. Defaults to 0.05.
        pseudoInverse : bool, optional
            Forces Moore-Penrose pseudoinverse instead of default inverse matrix operation. Defaults to False.
        device : str, optional 
            Device to use for computations. Could be 'cpu' or 'cuda' if CUDA is available. Defaults to 'cpu'.

    Returns
    -------
        result: revTWMRResult
            Named tuple containing attributes:
            Alpha : float
                Causal effect estimated by revTWMR before applying heterogeneity test
            Se : float
                Standard error of Alpha
            Pval : float
                Pvalue calculated from Alpha and Se
            N : float
                number of SNPs used as instruments 
            HetP : float
                Original P-value for heterogeneity test
            AlphaIterative : float
                Causal effect estimated by revTWMR after removing the SNPs detected as outliers by the heterogenity test
            SeIterative : float
                standard error of AlphaIterative
            PvalIterative : float
                Pvalue calculated from AlphaIterative and SeIterative
            HetPIterative : float
                Pvalue of heterogenity test after outlier removal
            rsname : List[str]
                SNPs left after outlier removal
    """

    if not isinstance(qtlExposureEffects, np.ndarray) and not isinstance(qtlExposureEffects, torch.Tensor):
        raise ValueError('qtlExposureEffect must be either numpy array or torch tensor')
    
    if not isinstance(qtlTraitGWASEffects, np.ndarray) and not isinstance(qtlTraitGWASEffects, torch.Tensor):
        raise ValueError('qtlTraitGWASEffect must be either numpy array or torch tensor')
    
    if not isinstance(qtlLabels, np.ndarray) and not isinstance(qtlLabels, list):
        raise ValueError('qtlLabels must be either numpy array or list')
    
    if not isinstance(gwasSizes, np.ndarray) and not isinstance(gwasSizes, torch.Tensor):
        raise ValueError('gwasSize must be either numpy array or torch tensor')
    
    if not isinstance(pValIterativeThreshold, float):
        raise ValueError('pValIterativeThreshold must be a number')
    
    if not isinstance(qtlExpSize, float):
        raise ValueError('qtlExpSize must be a number')
    

    if torch.cuda.is_available() and device == 'cuda':
        device = 'cuda'
    elif not torch.cuda.is_available():
        device = 'cpu'
        warnings.warn('CUDA is not available, reverting to CPU')

    if len(qtlExposureEffects.shape) != 1:
        raise ValueError('qtlExposureEffects must be a 1D array/tensor')
    
    if len(qtlTraitGWASEffects.shape) != 1:
        raise ValueError('qtlTraitGWASEffects must be a 1D array/tensor')
    
    if len(qtlLabels.shape) != 1:
        raise ValueError('qtlLabels must be a 1D array/List')
    
    if len(gwasSizes.shape) != 1:
        raise ValueError('gwasSizes must be a 1D array/tensor')
        
    if len(np.unique([a.shape[0] for a in [qtlExposureEffects, qtlTraitGWASEffects, qtlLabels, gwasSizes]])) > 1:
        raise ValueError('All input arrays must contain same amount of QTLs')
    


    nonzeromask = (abs(qtlExposureEffects)  < abs(qtlTraitGWASEffects))
    if isinstance(qtlExposureEffects, np.ndarray):        
        qtlExposureEffects = torch.from_numpy(qtlExposureEffects)
    gamma = qtlExposureEffects[nonzeromask].to(device)

    if isinstance(qtlTraitGWASEffects, np.ndarray):
        qtlTraitGWASEffects = torch.from_numpy(qtlTraitGWASEffects)
    beta = qtlTraitGWASEffects[nonzeromask].to(device)

    if isinstance(gwasSizes, np.ndarray):
        gwasSizes =torch.from_numpy(gwasSizes)
    ngwas = gwasSizes[nonzeromask].to(device)

    rsname = qtlLabels[nonzeromask]

        

    D = torch.eye(len(beta), dtype=torch.float64)
        
    alpha = Alpha(beta,gamma, D, pseudoInverse, device)
    se = SE_iterative(beta, gamma, ngwas, qtlExpSize, pseudoInverse, device)

    pval = 2*(1 - scipy.stats.norm.cdf(float(abs(alpha/se))))

    pvalOriginal = pval
    seOriginal = se
    alphaOriginal = alpha
    Nstart = len(gamma)
    N = Nstart
    phetORIGINAL, d = CohrainQ(alpha, beta, gamma, se, qtlExpSize, ngwas, device)
    phet = phetORIGINAL

    while  (phet < pValIterativeThreshold) & (N > 3):
            beta = beta[d< d.max()]
            gamma = gamma[d< d.max()]
            D = D[:, d < d.max()]
            D = D[d< d.max()]
            rsname = rsname[d< d.max()]
            ngwas = ngwas[d < d.max()]
            alpha = Alpha(beta,gamma, D, pseudoInverse, device)
            se = SE_iterative(beta, gamma, ngwas, qtlExpSize, pseudoInverse, device)
            z = alpha/se
            pval = 2*(1 - scipy.stats.norm.cdf(torch.abs(z).numpy()))
            phet, d = CohrainQ(alpha, beta, gamma, se, qtlExpSize, ngwas, device)
            N = len(d)        
    
    result = revTWMRResult(
        Alpha=float(alphaOriginal.cpu().numpy()), 
        Se=float(seOriginal.cpu().numpy()), 
        Pval=float(pvalOriginal), 
        N=float(Nstart), 
        HetP=float(phetORIGINAL),
        AlphaIterative=float(alpha.cpu().numpy()), 
        SeIterative=float(se.cpu().numpy()), 
        PvalIterative=float(pval), 
        HetPIterative=float(phet), 
        rsname=list(rsname))

    return result