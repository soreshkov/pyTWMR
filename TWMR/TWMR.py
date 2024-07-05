from collections import namedtuple
from typing import Union, List
import warnings


import numpy as np
import scipy.stats
import torch


    
TWMRresult = namedtuple('TWMRresult', ['Alpha', 'Se', 'Pval', 'D', 'HetP'])

def TWMR(
        beta: Union[np.ndarray, torch.Tensor], 
        gamma: Union[np.ndarray, torch.Tensor], 
        nEQTLs:int, 
        NGwas:int, 
        ldMatrix = None, 
        pseudoInverse = False, 
        device='cpu'):
    '''
    Perform Two-Sample Mendelian Randomization (TWMR) analysis.

    Parameters
    ----------

    beta : Union[np.ndarray, torch.Tensor]
        Matrix of QTL effect sizes for the exposure.
    gamma : Union[np.ndarray, torch.Tensor]
        Vector of QTL effect sizes for the outcome.
    nEQTLs : int
        Number of EQTLs in QTL study.
    NGwas : int
        Number of samples in genome-wide association study.
    ldMatrix : Union[np.ndarray, torch.Tensor], optional
        Linkage disequilibrium matrix. Default is None.
    pseudoInverse : bool, optional
        If True, use Moore-Penrose pseudo-inverse for matrix inversion. Default is False.
    device : str, optional
        Device to use (e.g., 'cpu' or 'cuda'). Default is 'cpu'.

    Returns
    -------
    result: TWMRresult
        Result of TWMR analysis.
        Attributes:
        Alpha : np.ndarray
            Causal effect estimated by TWMR.
        Se : np.ndarray
            Standard error of Alpha.
        Pval : float
            P-value calculated from Alpha and Se.
        D : np.ndarray
            Cohrian Q statistics for QTLs.
        HetP : np.ndarray
            P-value for heterogeneity test.   
    '''

    if not isinstance(beta, np.ndarray) and not isinstance(beta, torch.Tensor):
        raise ValueError("Beta should be a numpy array or torch tensor")
    
    if not isinstance(gamma, np.ndarray) and not isinstance(gamma, torch.Tensor):
        raise ValueError("Gamma should be a numpy array or torch tensor")
    
    if len(beta.shape) != 2:
        raise ValueError("Beta should be a 2D matrix")
    
    if len(gamma.shape) != 1:
        raise ValueError("Gamma should be a vector")
    
    if len(beta) != len(gamma):
        raise ValueError("Beta and gamma should have the same length")

    if torch.cuda.is_available() and device == 'cuda':
        device = 'cuda'
    elif not torch.cuda.is_available():
        device = 'cpu'
        warnings.warn('CUDA is not available, reverting to CPU')

    if isinstance(beta, np.ndarray):
        beta = torch.from_numpy(beta)

    if isinstance(gamma, np.ndarray):
        gamma = torch.from_numpy(gamma)

    beta = beta.to(torch.float64).to(device)
    gamma = gamma.to(torch.float64).to(device)

    if ldMatrix is None:
        ldMatrix = torch.eye(len(beta), dtype=torch.float64).to(device)
    else:
        if isinstance(ldMatrix, np.ndarray):
            ldMatrix = torch.from_numpy(ldMatrix)
        ldMatrix = ldMatrix.to(torch.float64).to(device)
    
    inverse = torch.linalg.pinv if pseudoInverse else torch.linalg.inv
    tS = beta.T @ inverse(ldMatrix) @ beta
    H = (1 -1/np.sqrt(3781))* tS + 1/np.sqrt(3781) * torch.eye(tS.shape[0]).to(device) 
    alpha = inverse(H) @ (beta.T @ inverse(ldMatrix) @ gamma)
    D  = ldMatrix
    invD = inverse(D)
    GCG = beta.T @ (invD @ beta)
    GCG = (1-1/np.sqrt(3781)) * GCG + 1/np.sqrt(3781) * torch.eye(GCG.shape[0]).to(device)

    GCG_inv = inverse(GCG)
    df_dg = GCG_inv @ beta.T @ invD
    p1= gamma @ invD @ (beta @ GCG_inv @ beta.T @ invD + torch.eye(len(beta)).to(device))
    p1 = torch.kron(GCG_inv.contiguous() , p1)
    p22 = GCG_inv @ beta.T @ invD
    p21 = - gamma @ invD @ beta @ GCG_inv
    df_dG = p1 + torch.kron(p21,p22)
    SEs = torch.cat((torch.full((beta.shape[0] * beta.shape[1],), 1 / np.sqrt(nEQTLs)), torch.full((len(gamma),), 1/np.sqrt(NGwas)))).to(device)
    R = torch.eye(beta.shape[1]+1).to(device)

    sigma = (SEs[:, None]  @ SEs[:, None] .T) * (torch.kron(ldMatrix.contiguous(), R.contiguous()))
    J = torch.cat((df_dG, df_dg), axis=1)
    V = J @ sigma @ J.T
    se = torch.sqrt(torch.diag(V))

    z = alpha/se
    pval = 2*(1 - scipy.stats.norm.cdf(torch.abs(z).cpu().numpy()))

    #Cohrain Q
    d = gamma - alpha @ beta.T
    var_d = 1 / nEQTLs  + se * se * beta * beta + alpha * alpha * 1 / NGwas + (1 / NGwas) * se * se
    phet = np.array([ 1 - scipy.stats.chi2.cdf(float(d @ torch.diag(var_d[:, v]) @ d), len(d)-1) for v in range(var_d.shape[1])])

    result = TWMRresult(Alpha=alpha.cpu().numpy(), Se=se.cpu().numpy(), Pval=pval, D=d.cpu().numpy(), HetP=phet)
    return result


def QCorrectedTWMR(
        beta: np.ndarray, 
        gamma: np.ndarray, 
        nEQTLs:int, 
        NGwas:int, 
        rsnames:Union[np.ndarray, List[str]],
        ldMatrix: Union[np.ndarray, torch.Tensor] = None,
        threshold=0.05, 
        pseudoInverse = False, device='cpu'):
    """
    TWMR method with heterogenity outlier correction. Excludes SNPs by Cohrain Q heterogenity test until SNP heterogenity is no longer observed.

    Parameters
    ----------

    beta : Union[np.ndarray, torch.Tensor]
        Matrix of QTL effect sizes for the exposure.
    gamma : Union[np.ndarray, torch.Tensor]
        Vector of QTL effect sizes for the outcome.
    nEQTLs : Union[np.ndarray, List[str]]
        Number of EQTLs in QTL study.
    NGwas : int
        Number of GWAS (genome-wide association study) samples.
    ldMatrix : Union[np.ndarray, torch.Tensor], optional
        Linkage disequilibrium matrix. Default is None.
    pseudoInverse : bool, optional
        If True, use Moore-Penrose pseudo-inverse for matrix inversion. Default is False.
    device : str, optional
        Device to use (e.g., 'cpu' or 'cuda'). Default is 'cpu'.
        
    Returns
    -------
    result: TWMRresult
        Result of TWMR analysis with outliers removed by iterative Cohrain Q test.
        Attributes:
        Alpha : float
            Causal effect estimated by TWMR.
        Se : float
            Standard error of Alpha.
        Pval : float
            P-value calculated from Alpha and Se.
        D : float
            Cohrian Q statistics for Alphas.
        HetP : float
            P-value for heterogeneity test.   
    rsnames: List[str]
        List of QTLs after removing outliers.
    """

    result = TWMR(beta=beta, gamma=gamma, nEQTLs=nEQTLs, NGwas=NGwas, ldMatrix = ldMatrix, pseudoInverse = pseudoInverse, device=device)

    d = result.D
    phet = result.HetP
    N = len(gamma)
    corrected_result = result

    while (N> 2 and (phet < threshold).all()):
        mask = d < np.max(d)
        beta = beta[mask]
        gamma = gamma[mask]
        rsnames = rsnames[mask]
        corrected_result = TWMR(beta=beta, gamma=gamma, nEQTLs=nEQTLs, NGwas=NGwas, ldMatrix = ldMatrix, pseudoInverse = pseudoInverse, device=device)
        d = corrected_result.D
        phet = corrected_result.HetP
        N = len(gamma)
    
    return corrected_result, list(rsnames)
