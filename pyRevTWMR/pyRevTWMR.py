import numpy as np

try:
    import torch

    def revTWMR(beta, gamma, NGwas, ldMatrix = None, pseudoInverse = False, device='cpu'):
        """
        Perform reverse Two-Step Mendelian Randomization (revTWMR) analysis.

        Args:
            beta (numpy.ndarray): Vector of QTL effect sizes for the outcome.
            gamma (numpy.ndarray): Array of QTL effect sizes for the exposure.
            NGwas (int): The sample size of the GWAS.
            ldMatrix (numpy.ndarray, optional): The linkage disequilibrium matrix. Defaults to None.
            pseudoInverse (bool, optional): Whether to use the pseudo-inverse or the inverse of the ldMatrix. Defaults to False.
            device (str, optional): The device to perform the calculations on. Defaults to 'cpu'.

        Returns:
            tuple: A tuple containing the alpha estimates and the standard errors.

        """
        
        inverse = torch.linalg.pinv if pseudoInverse else torch.linalg.inv

        beta = torch.from_numpy(beta).to(device, dtype=torch.float32)
        gamma = torch.from_numpy(gamma).to(device, dtype=torch.float32)
        ldMatrix = torch.eye(len(gamma)).to(device, dtype=torch.float32) if ldMatrix is None else torch.from_numpy(ldMatrix).to(device, dtype=torch.float32)

        tS = beta.T @ (inverse(ldMatrix) @  beta)
        alpha = torch.nan_to_num(beta.T  @ (inverse(ldMatrix) @ gamma)  / tS)
        D  = ldMatrix
        invD = inverse(D)
        GCG = beta.T @ (invD @ beta)

        df_dg =  beta.T @ invD / GCG
        p1 = gamma.T @ invD 
        p1 = p1 @  ( - (beta[:, None] / GCG @ beta[None,:]) @ invD  + torch.eye(len(beta)).to(device, dtype=torch.float32))
        p1 = torch.kron(1/ GCG , p1)
        
        p21 = - gamma.T @ invD @ beta / GCG        
        p22 =  beta.T @ invD /GCG

        df_dG = p1 + torch.kron(p21,p22).view(p1.shape[0], p1.shape[1])

        alp = (1/GCG) * (beta.T @ invD @ gamma)

        Vgam = torch.var(gamma - beta[:, None] @ alp[None, :], dim=0).to(device)
        SEs = torch.cat((torch.full_like(beta, 1/np.sqrt(NGwas)).to(device), torch.sqrt(Vgam))  )
        R = torch.eye(2).to(device, dtype=torch.float32)
        J = torch.cat((df_dG, df_dg[None, :].repeat(df_dG.shape[0],1)), axis=1)
        sigma = torch.nan_to_num((SEs @ SEs.T) * (torch.kron(ldMatrix, R)))
        V = J @ sigma @ J.T        
        se = torch.nan_to_num(torch.sqrt(V))
        return alpha.numpy(), torch.diagonal(se).numpy()


except ImportError:
    def revTWMR(beta, gamma, nEQTLs, NGwas, ldMatrix = None, pseudoInverse = False):
        raise ImportError("Please install pytorch to use the GPU version of TWMR")
    