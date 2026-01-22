import numpy as np
from scipy.stats import dirichlet

def LCA_py(Z, par):
    P_LCA = par[Z - 1, :, :]
    return P_LCA


def sim_data_LCA(N=5000, I=100, L=10, poly_value=5):

    if np.isscalar(poly_value):
        poly_value = np.repeat(poly_value, I)
    else:
        poly_value = np.array(poly_value)
        if len(poly_value) != I:
            poly_value = np.repeat(poly_value[0], I)

    poly_max = int(np.max(poly_value))

    Z = np.random.choice(np.arange(1, L + 1), size=N, replace=True)
    par = np.full((L, I, poly_max), np.nan)

    for i in range(I):
        for l in range(L):
            k = poly_value[i]
            if k > 0:
                probs = dirichlet.rvs(alpha=np.ones(k), size=1)[0]
                par[l, i, :k] = probs

    P_LCA = LCA_py(Z, par)

    P_LCA_cumulative = np.zeros((N, I, poly_max))
    P_random = np.zeros((N, I, poly_max))

    for p in range(N):
        for i in range(I):
            k = poly_value[i]
            for po in range(k):
                P_LCA_cumulative[p, i, po] = np.sum(P_LCA[p, i, :po + 1])

    P_random[:, :, 0] = np.random.uniform(0, 1, size=(N, I))
    for po in range(1, poly_max):
        P_random[:, :, po] = P_random[:, :, 0]

    response = np.zeros((N, I), dtype=int)
    for i in range(I):
        k = poly_value[i]
        temp = P_random[:, i, :k] >= P_LCA_cumulative[:, i, :k]
        response[:, i] = np.sum(temp, axis=1)

    return {
        'response': response,
        'par': par,
        'Z': Z,
        'poly_value': poly_value
    }
    
def adjust_data(response):
    device = response.device
    N, I = response.shape

    poly_value = np.zeros(I, dtype=int)
    poly_orig = []

    max_k = 0
    unique_list = []
    for i in range(I):
        unique_vals, inverse_indices = torch.unique(response[:, i], return_inverse=True)
        unique_list.append((unique_vals, inverse_indices))
        k = len(unique_vals)
        poly_value[i] = k
        if k > max_k:
            max_k = k
    poly_max = max_k
    
    poly_orig = np.full((I, poly_max), np.nan)
    response_adjusted = torch.zeros((N, I), dtype=torch.float32, device=device)
    for i in range(I):
        unique_vals, inverse_indices = unique_list[i]
        poly_orig[i, :poly_value[i]] = unique_vals.detach().cpu().numpy().astype(int)
        response_adjusted[:, i] = inverse_indices
    
    return {
        'poly_orig': poly_orig,
        'poly_value': poly_value,
        'poly_max': poly_max,
        'response': response_adjusted
    }

import torch
import numpy as np
from torch.distributions import MultivariateNormal, Gamma

def sim_data_LPA(N, I, L, constraint="VV", distribution="random", device="GPU", dtype=torch.float32):
    # Validate device
    device = torch.device(device if device is not None else 'cpu')
    
    # Validate constraint
    constraints_all = ["E", "V", "E0", "V0", "EE", "VE", "EV", "VV"]
    if constraint not in constraints_all:
        raise ValueError(f"Invalid constraint '{constraint}'. Must be one of {constraints_all}")
    
    # Constraint validation for univariate case
    is_univariate = (I == 1)
    if I > 1 and constraint in ["E", "V"]:
        raise ValueError("Single-variable constraints ('E','V') require I = 1")
    if I == 1 and constraint not in ["E", "V"]:
        if constraint in ["E0", "EE", "EV"]:
            constraint = "E"
            print(f"I=1: Mapping '{constraint}' to univariate constraint 'E'")
        else:
            constraint = "V"
            print(f"I=1: Mapping '{constraint}' to univariate constraint 'V'")
    
    # Generate means: L classes, each with I variables in [-3, 3]
    means = [
        torch.rand(I, device=device, dtype=dtype) * 6.0 - 3.0
        for _ in range(L)
    ]  # List of (I,) tensors
    
    # Helper: Generate positive-definite covariance matrix
    def generate_cov_matrix(I, var_range=(0.0, 3.0)):
        # Create correlation matrix
        R = torch.eye(I, device=device, dtype=dtype)
        if I > 1:
            # Fill off-diagonals with random correlations
            for j in range(1, I):
                for i in range(j):
                    R[i, j] = R[j, i] = torch.rand(1, device=device, dtype=dtype).item() * 2.0 - 1.0
            
            # Ensure positive semi-definiteness
            eigvals, eigvecs = torch.linalg.eigh(R)
            eigvals = torch.clamp(eigvals, min=1e-5)
            R = eigvecs @ torch.diag(eigvals) @ eigvecs.T
            # Convert to correlation matrix
            D_inv = torch.diag(1.0 / torch.sqrt(torch.diag(R)))
            R = D_inv @ R @ D_inv
        
        # Generate standard deviations
        sds = torch.sqrt(torch.rand(I, device=device, dtype=dtype) * (var_range[1] - var_range[0]) + var_range[0])
        D = torch.diag(sds)
        
        # Construct covariance matrix
        Sigma = D @ R @ D
        
        # Final positive-definiteness check
        eigvals, eigvecs = torch.linalg.eigh(Sigma)
        eigvals = torch.clamp(eigvals, min=1e-5)
        Sigma = eigvecs @ torch.diag(eigvals) @ eigvecs.T
        return Sigma
    
    var_range = (0.0, 3.0)
    covs = []
    
    # Generate covariance matrices based on constraint
    if is_univariate:
        if constraint == "E":
            shared_var = torch.rand(1, device=device, dtype=dtype).item() * (var_range[1] - var_range[0]) + var_range[0]
            covs = [torch.tensor([[shared_var]], device=device, dtype=dtype) for _ in range(L)]
        else:  # "V"
            class_vars = torch.rand(L, device=device, dtype=dtype) * (var_range[1] - var_range[0]) + var_range[0]
            covs = [torch.tensor([[v.item()]], device=device, dtype=dtype) for v in class_vars]
    else:
        if constraint in ["EE", "VE", "EV"]:
            # Generate base correlation matrix
            base_R = torch.eye(I, device=device, dtype=dtype)
            if I > 1:
                for j in range(1, I):
                    for i in range(j):
                        base_R[i, j] = base_R[j, i] = torch.rand(1, device=device, dtype=dtype).item() * 2.0 - 1.0
                
                eigvals, eigvecs = torch.linalg.eigh(base_R)
                eigvals = torch.clamp(eigvals, min=1e-5)
                base_R = eigvecs @ torch.diag(eigvals) @ eigvecs.T
                # Convert to correlation matrix
                D_inv = torch.diag(1.0 / torch.sqrt(torch.diag(base_R)))
                base_R = D_inv @ base_R @ D_inv
            
            # Generate SDs for each class/dimension
            sds_list = [
                torch.sqrt(torch.rand(I, device=device, dtype=dtype) * (var_range[1] - var_range[0]) + var_range[0])
                for _ in range(L)
            ]
            
            if constraint == "EE":
                shared_sds = torch.sqrt(torch.rand(I, device=device, dtype=dtype) * (var_range[1] - var_range[0]) + var_range[0])
                Sigma_shared = torch.diag(shared_sds) @ base_R @ torch.diag(shared_sds)
                covs = [Sigma_shared.clone() for _ in range(L)]
            elif constraint == "VE":
                lambdas = torch.rand(L, device=device, dtype=dtype) * 0.4 + 0.8  # [0.8, 1.2]
                covs = []
                for l in range(L):
                    D = torch.diag(sds_list[l])
                    Sigma = lambdas[l] * (D @ base_R @ D)
                    covs.append(Sigma)
            elif constraint == "EV":
                shared_sds = torch.sqrt(torch.rand(I, device=device, dtype=dtype) * (var_range[1] - var_range[0]) + var_range[0])
                Sigma_base = torch.diag(shared_sds) @ base_R @ torch.diag(shared_sds)
                lambdas = torch.rand(L, device=device, dtype=dtype) * 0.2 + 0.9  # [0.9, 1.1]
                covs = [lambdas[l] * Sigma_base for l in range(L)]
        else:  # E0, V0, VV
            if constraint == "E0":
                shared_vars = torch.rand(I, device=device, dtype=dtype) * (var_range[1] - var_range[0]) + var_range[0]
                covs = [torch.diag(shared_vars) for _ in range(L)]
            elif constraint == "V0":
                covs = [torch.diag(
                    torch.rand(I, device=device, dtype=dtype) * (var_range[1] - var_range[0]) + var_range[0]
                ) for _ in range(L)]
            elif constraint == "VV":
                covs = [generate_cov_matrix(I, var_range) for _ in range(L)]
    
    # Generate class sizes (Z assignments)
    if distribution == "random":
        # Dirichlet-like distribution using Gamma
        alpha = torch.full((L,), 3.0, device=device, dtype=dtype)
        gamma_dist = Gamma(alpha, torch.ones(L, device=device, dtype=dtype))
        p = gamma_dist.sample()
        p = p / p.sum()
        sizes_class = (p * N).round().to(torch.long)
        diff = N - sizes_class.sum().item()
        if diff != 0:
            max_idx = torch.argmax(sizes_class)
            sizes_class[max_idx] += diff
        Z = torch.cat([torch.full((size.item(),), l + 1, device=device, dtype=torch.long) 
                      for l, size in enumerate(sizes_class)])
    elif distribution == "uniform":
        Z = torch.randint(1, L + 1, (N,), device=device, dtype=torch.long)
    else:
        raise ValueError("Invalid 'distribution'. Choose 'random' or 'uniform'.")
    
    # Generate response data
    response = []
    for l in range(L):
        idx = (Z == (l + 1))
        n_l = idx.sum().item()
        if n_l == 0:
            continue
        # Create multivariate normal distribution for class l
        dist = MultivariateNormal(means[l], covs[l])
        samples = dist.sample((n_l,))
        response.append(samples)
    
    response = torch.cat(response, dim=0) if response else torch.empty(0, I, device=device, dtype=dtype)
    
    # Clamp values to [-10, 10]
    response = torch.clamp(response, min=-10.0, max=10.0)
    
    # Shuffle observations
    shuffle_idx = torch.randperm(N, device=device)
    response = response[shuffle_idx]
    Z = Z[shuffle_idx]
    
    # Prepare metadata
    feature_names = ["V"] if is_univariate else [f"V{i+1}" for i in range(I)]
    class_names = [f"Class{l+1}" for l in range(L)]
    
    # Convert means to tensor (L, I)
    means_tensor = torch.stack(means) if means else torch.empty(0, I, device=device, dtype=dtype)
    
    return {
        "response": response,
        "means": means_tensor,
        "covs": covs,
        "Z": Z,
        "true_constraint": constraint,
        "feature_names": feature_names,
        "class_names": class_names
    }
