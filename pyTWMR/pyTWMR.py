import numpy as np

try:
    import torch

    def TWMR(beta: np.ndarray, gamma: np.ndarray, nEQTLs:int, NGwas:int, ldMatrix = None, pseudoInverse = False, device='cpu'):
        '''
    Perform Two-Sample Mendelian Randomization (TWMR) analysis.

    Parameters:
    - beta (np.ndarray): Matrix of QTL effect sizes for the exposure.
    - gamma (np.ndarray): Vector of QTL effect sizes for the outcome.
    - nEQTLs (int): Number of EQTLs in QTL study.
    - NGwas (int): Number of GWAS (genome-wide association study) samples.
    - ldMatrix (np.ndarray, optional): Linkage disequilibrium matrix. Default is None.
    - pseudoInverse (bool, optional): If True, use Moore-Penrose pseudo-inverse for matrix inversion. Default is False.
    - device (str, optional): Device to use (e.g., 'cpu' or 'cuda'). Default is 'cpu'.

    Returns:
    - Result of TWMR analysis.

    Notes:
    - The TWMR function performs a Two-Sample Mendelian Randomization analysis using the specified parameters.
    - The ldMatrix parameter is optional. If not provided, it can be computed internally.
    - If pseudoInverse is set to True, pseudo-inverse is used for matrix inversion.
    - The device parameter specifies the computational device to use, such as 'cpu' or 'cuda'.
        '''
        
        if not isinstance(beta, np.ndarray):
            raise ValueError("Beta should be a numpy array")
        
        if not isinstance(gamma, np.ndarray):
            raise ValueError("Gamma should be a numpy array")
        
        if beta.ndim != 2:
            raise ValueError("Beta should be a matrix")
        
        if gamma.ndim != 1:
            raise ValueError("Gamma should be a vector")
        
        if len(beta) != len(gamma):
            raise ValueError("Beta and gamma should have the same length")

        
        beta = torch.from_numpy(beta).to(device)
        gamma = torch.from_numpy(gamma).to(device)
        if ldMatrix is None:
            ldMatrix = torch.eye(len(beta)).to(device)
        else:
            ldMatrix = torch.from_numpy(ldMatrix).to(device)
        inverse = torch.linalg.pinv if pseudoInverse else torch.linalg.inv

        tS = beta.T @ inverse(ldMatrix) @ beta
        alpha = inverse(tS) @ (beta.T @ inverse(ldMatrix) @ gamma)
        D  = ldMatrix
        invD = inverse(D)
        GCG = beta.T @ (invD @ beta)
        GCG_inv = inverse(GCG)
        df_dg = GCG_inv @ beta.T @ invD
        p1= gamma @ invD @ (beta @ GCG_inv @ beta.T @ invD + torch.eye(len(beta)))
        p1 = torch.kron(GCG_inv.contiguous() , p1)
        p22 = GCG_inv @ beta.T @ invD
        p21 = - gamma.T @ invD @ beta @ GCG_inv
        df_dG = p1 + torch.kron(p21,p22)
        SEs = torch.cat((torch.full((beta.shape[0] * beta.shape[1],), 1 / np.sqrt(nEQTLs)), torch.full((len(gamma),), 1/np.sqrt(NGwas))))
        R = torch.eye(beta.shape[1]+1)
        sigma = (SEs @ SEs.T) * (torch.kron(ldMatrix, R))
        J = torch.cat((df_dG, df_dg), axis=1)
        V = J @ sigma @ J.T
        se = torch.sqrt(V)
        return alpha.numpy(), se.numpy()
        
            


except ImportError:
    def TWMR(beta, gamma, nEQTLs, NGwas, ldMatrix = None, pseudoInverse = False):
        raise ImportError("Please install pytorch to use the GPU version of TWMR")