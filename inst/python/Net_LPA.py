import torch
from torch import linalg as LA
import numpy as np
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.optim.lr_scheduler import ReduceLROnPlateau
import matplotlib.pyplot as plt
import math
from sklearn.cluster import KMeans
from itertools import combinations
from torch.nn import TransformerEncoder, TransformerEncoderLayer
from scipy.linalg import eigh
import warnings
import copy
from torch.distributions import Dirichlet
import random

class LPAnet(nn.Module):
    def __init__(self, response, L=5, par_ini=None, constraint="VV",
                 hidden_layers=[32, 32], activation_function='tanh', 
                 d_model=None, nhead=None, dim_feedforward=None, eps=1e-8):
        super(LPAnet, self).__init__()
        self.L = L
        self.response = response.float()
        self.device = self.response.device
        self.N, self.I = self.response.shape
        self.constraint = constraint
        self.npar = self._compute_npar()
        self.eps = eps
        self.shared_mask = None
        
        self._init_shared_mask()
        
        data_std = torch.std(self.response, dim=0)
        data_std = data_std + self.eps
        self.register_buffer('data_std', data_std)
        self.register_buffer('data_mean', torch.mean(self.response, dim=0))
        
        self.input_dim = self.response.shape[1]

        activation_dict = {
            'relu': nn.ReLU(),
            'sigmoid': nn.Sigmoid(),
            'tanh': nn.Tanh(),
            'elu': nn.ELU(),
            'softmax': nn.Softmax(dim=1)
        }
        if activation_function not in activation_dict:
            raise ValueError(f"Unsupported activation function: {activation_function}")
        activation = activation_dict[activation_function]

        if par_ini is not None and not isinstance(par_ini, str):
            par_ini["means"] = np.asarray(par_ini["means"], dtype=np.float32)
            means_raw = torch.tensor(par_ini["means"], dtype=torch.float32, device=self.device)
            means_raw = torch.atanh(torch.clamp(means_raw, min=-1.99, max=1.99)/2.0)
            self.means = nn.Parameter(means_raw)

            par_ini["covs"] = np.asarray(par_ini["covs"], dtype=np.float32)
            target_covs = torch.tensor(par_ini["covs"], dtype=torch.float32, device=self.device)
            raw_covs = self.inverse_activation_for_covs(target_covs, eps=self.eps)
            self.covs = nn.Parameter(raw_covs)
            
            self.P_Z_init = torch.from_numpy(np.asarray(par_ini["P.Z"], dtype=np.float32)).float().to(self.device).squeeze()

        elif par_ini == "kmeans":
            _, covs_np, means_np, P_Z_np = self.kmeans_classify(self.response, self.L, constraint=self.constraint)

            means_raw = torch.from_numpy(means_np).float().to(self.device)
            means_raw = torch.atanh(torch.clamp(means_raw, min=-1.99, max=1.99)/2.0)
            self.means = nn.Parameter(means_raw)
            
            target_covs = torch.from_numpy(covs_np).float().to(self.device)
            raw_covs = self.inverse_activation_for_covs(target_covs, eps=self.eps)
            self.covs = nn.Parameter(raw_covs)
            
            self.P_Z_init = torch.from_numpy(np.asarray(P_Z_np, dtype=np.float32)).float().to(self.device).squeeze()
        else:
            I = self.I
            L = self.L
            response_np = self.response.detach().cpu().numpy()
            q25 = np.percentile(response_np, 10, axis=0)
            q75 = np.percentile(response_np, 90, axis=0)
            means_np = np.zeros((L, I))
            for i in range(I):
                means_np[:, i] = np.random.uniform(q25[i], q75[i], size=L)
            covs_np = np.zeros((I, I, L))
            for l in range(L):
                covs_np[:, :, l] = np.eye(I)
            
            means_raw = torch.from_numpy(means_np).float().to(self.device)
            means_raw = torch.atanh(torch.clamp(means_raw, min=-1.99, max=1.99)/2.0)
            self.means = nn.Parameter(means_raw)
            
            target_covs = torch.from_numpy(covs_np).float().to(self.device)
            raw_covs = self.inverse_activation_for_covs(target_covs, eps=self.eps)
            self.covs = nn.Parameter(raw_covs)
            
            alpha_z = torch.ones(L) * 3.0
            dirichlet_z = Dirichlet(alpha_z)
            P_Z_sample = dirichlet_z.sample().cpu().numpy()
            self.P_Z_init = torch.from_numpy(P_Z_sample).float().to(self.device)

        layers = []
        input_dim = self.input_dim
        for hidden_units in hidden_layers:
            linear_layer = nn.Linear(input_dim, hidden_units, dtype=torch.float32)
            layers.append(linear_layer)
            layers.append(nn.BatchNorm1d(hidden_units, dtype=torch.float32))
            layers.append(copy.deepcopy(activation))
            input_dim = hidden_units
        final_layer = nn.Linear(input_dim, self.L, dtype=torch.float32)
        layers.append(final_layer)
        self.network = nn.Sequential(*layers)

        if d_model is None: 
            d_model = max(8, int(max(1, np.floor(np.log10(max(2, self.N))) *
                                        np.floor(np.log2(max(2, self.L))))))
        if nhead is None: 
            nhead = min(2, max(1, d_model // 4))
        if d_model % nhead != 0:
            d_model += (nhead - d_model % nhead)
        if dim_feedforward is None: 
            dim_feedforward = max(d_model, 16)
        
        self.embed_proj = nn.Linear(self.L, d_model, dtype=torch.float32)
        encoder_layer = TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_feedforward,
            batch_first=True,
            dtype=torch.float32
        )
        num_attn_layers = min(2, max(1, nhead))
        self.attn_layer = TransformerEncoder(encoder_layer, num_layers=num_attn_layers)
        self.output_proj = nn.Linear(d_model, self.L, dtype=torch.float32)

        for layer in self.network:
            if isinstance(layer, nn.Linear):
                nn.init.kaiming_normal_(layer.weight, nonlinearity='relu')
                nn.init.zeros_(layer.bias)
        
        nn.init.kaiming_normal_(self.embed_proj.weight, nonlinearity='relu')
        nn.init.zeros_(self.embed_proj.bias)
        nn.init.kaiming_normal_(self.output_proj.weight, nonlinearity='relu')
        nn.init.zeros_(self.output_proj.bias)

        with torch.no_grad():
            pz = torch.clamp(self.P_Z_init, min=self.eps)
            pz = pz / pz.sum()
            logits_init = torch.log(pz + self.eps)
            if logits_init.numel() == self.L:
                self.output_proj.bias.copy_(logits_init)
        
        self.to(self.device)
    
    def generate_random_correlation_matrix(self):
        I = self.I
        A = np.random.normal(size=(I, I))
        C = A @ A.T
        standard_deviations = np.sqrt(np.diag(C))
        D = np.outer(standard_deviations, standard_deviations)
        R = np.divide(C, D, out=np.zeros_like(C), where=D!=0)
        np.fill_diagonal(R, 1.0)
        R = np.clip(R, -1.0, 1.0)
        return R

    def _init_shared_mask(self):
        I = self.I
        constraint = self.constraint
        device = self.device
        
        shared_mask = torch.zeros((I, I), dtype=torch.bool, device=device)
        
        if I == 1:
            if constraint in ["E", "EE"]:
                shared_mask[0, 0] = True
            self.shared_mask = shared_mask
            return shared_mask
        
        constraint_used = constraint
        if constraint_used in ["E", "V"]:
            constraint_used = "EE" if constraint_used == "E" else "VV"
        
        if constraint_used == "EE":
            shared_mask[:] = True 
        elif constraint_used == "E0":
            diag_indices = torch.arange(I, device=device)
            shared_mask[diag_indices, diag_indices] = True
        elif constraint_used == "EV":
            diag_indices = torch.arange(I, device=device)
            shared_mask[diag_indices, diag_indices] = True
        elif constraint_used in ["VV", "V0"]:
            shared_mask[:] = False
        elif constraint_used == "VE":
            off_diag_mask = ~torch.eye(I, dtype=torch.bool, device=device)
            shared_mask[off_diag_mask] = True
        
        elif isinstance(constraint_used, list) or isinstance(constraint_used, tuple):
            for pair in constraint_used:
                if len(pair) == 2:
                    i, j = pair
                    j -= 1
                    i -= 1
                    if i < I and j < I:
                        shared_mask[i, j] = True
                        shared_mask[j, i] = True
        
        self.shared_mask = shared_mask
        return shared_mask

    def _compute_npar(self):
        I = self.I
        L = self.L
        constraint = self.constraint
        npar = L * I + (L - 1)
        if I == 1:
            if constraint in ["V", "VV", "V0"]:
                npar += L
            else:
                npar += 1
            return npar

        if constraint == "E":
            constraint_used = "EE"
        elif constraint == "V":
            constraint_used = "VV"
        else:
            constraint_used = constraint
        
        if constraint_used == "E0":
            npar += I
        elif constraint_used == "V0":
            npar += L * I
        elif constraint_used == "EE":
            npar += I * (I + 1) // 2
        elif constraint_used == "VE":
            npar += L * I + I * (I - 1) // 2
        elif constraint_used == "EV":
            npar += I + L * (I * (I - 1) // 2)
        elif constraint_used == "VV":
            npar += L * (I * (I + 1) // 2)
        else:
            npar += L * (I * (I + 1) // 2)
            npar -= (L-1) * len(constraint_used)
        
        return npar

    def _apply_means_activation(self, raw_means):
        activated_means = torch.tanh(raw_means) * 2.0
        return activated_means

    def _apply_covs_activation(self, raw_covs):
        I, _, L = raw_covs.shape
        device = raw_covs.device
        dtype = raw_covs.dtype

        tril_idx = torch.tril_indices(I, I, offset=0, device=device)
        num_tril = tril_idx.shape[1]
        diag_mask = (tril_idx[0] == tril_idx[1])
        off_diag_mask = ~diag_mask

        tril_vals = raw_covs[tril_idx[0], tril_idx[1]]
        
        if diag_mask.any():
            diag_vals = tril_vals[diag_mask]
            diag_pos = torch.clamp(F.softplus(diag_vals), min=0.01, max=3.99)
            tril_vals[diag_mask] = diag_pos
        
        if off_diag_mask.any():
            off_diag_vals = tril_vals[off_diag_mask]
            tril_vals[off_diag_mask] = torch.tanh(off_diag_vals) * 4.0
        
        Lmats = torch.zeros((L, I, I), device=device, dtype=dtype)
        Lmats[:, tril_idx[0], tril_idx[1]] = tril_vals.T
        
        covs = torch.bmm(Lmats, Lmats.transpose(1, 2))
        activated_covs = covs.permute(1, 2, 0)
        return activated_covs

    def inverse_activation_for_covs(self, target_covs, eps=1e-10):
        I, I2, L = target_covs.shape
        assert I == I2, "Covariance matrix must be square"
        device = target_covs.device
        dtype = target_covs.dtype
    
        tril_idx = torch.tril_indices(I, I, offset=0, device=device)
        num_tril = tril_idx.shape[1]
        diag_mask = (tril_idx[0] == tril_idx[1])
        off_diag_mask = ~diag_mask
    
        covs_Lfirst = target_covs.permute(2, 0, 1)
        Lmats = torch.zeros((L, I, I), device=device, dtype=dtype)
        
        if self.shared_mask is None:
            self._init_shared_mask()
        shared_mask = self.shared_mask
        non_shared_mask = ~shared_mask
    
        for c in range(L):
            Sigma = covs_Lfirst[c]
            try:
                Lmat = torch.linalg.cholesky(Sigma)
            except RuntimeError:
                Lmat, _ = self.adaptive_cholesky_constrained(Sigma, non_shared_mask)
            Lmats[c] = Lmat
    
        tril_vals = Lmats[:, tril_idx[0], tril_idx[1]]
    
        if diag_mask.any():
            diag_vals = tril_vals[:, diag_mask]
            inv_diag = torch.log(torch.expm1(torch.clamp(diag_vals, min=0.01, max=3.99)) + eps)
            tril_vals[:, diag_mask] = inv_diag
    
        if off_diag_mask.any():
            off_diag_vals = tril_vals[:, off_diag_mask] / 4.0
            off_diag_vals = torch.clamp(off_diag_vals, min=-0.99, max=0.99)
            tril_vals[:, off_diag_mask] = torch.atanh(off_diag_vals)
            
        raw_covs = torch.zeros_like(target_covs)
        raw_covs[tril_idx[0], tril_idx[1]] = tril_vals.T
        return raw_covs
    
    @staticmethod
    def kmeans_classify(features, Y, L, max_attempts=5, nstart=1, eps=1e-4, random_state=None, constraint="VV"):
        features_np = features.detach().cpu().numpy()
        Y_np = Y.detach().cpu().numpy()
        mean = np.mean(Y_np, axis=0)
        std = np.std(Y_np, axis=0)
        std = np.where(std == 0, 1.0, std)
        Y_np_normalized = (Y_np - mean) / std
        N, I = Y_np.shape
        
        kmeans_attempts = 0
        kmeans_result = None
        cluster_assignments = None
        while kmeans_attempts < max_attempts:
            try:
                kmeans = KMeans(
                    n_clusters=L, 
                    n_init=nstart, 
                    max_iter=500, 
                    algorithm='lloyd'
                )
                kmeans.fit(Y_np_normalized)
                kmeans_result = kmeans
                cluster_assignments = kmeans.labels_
                break
            except Exception as e:
                kmeans_attempts += 1
                if kmeans_attempts == max_attempts:
                    warnings.warn(f"KMeans failed after {max_attempts} attempts: {str(e)}")
        
        if kmeans_result is None:
            rng = np.random.default_rng(random_state)
            cluster_assignments = rng.integers(0, L, size=N)
            global_mean = np.mean(Y_np, axis=0)
            means = np.tile(global_mean, (L, 1))
            P_Z = np.full(L, 1.0 / L)
        else:
            means = kmeans_result.cluster_centers_
            sizes = np.bincount(cluster_assignments, minlength=L)
            P_Z = np.maximum(sizes / N, 1e-8)
            P_Z = P_Z / np.sum(P_Z)

        try:
            if I == 1:
                covs_global = np.array([[np.var(Y_np, ddof=1)]])
            else:
                covs_global = np.cov(Y_np, rowvar=False, ddof=1)
        except Exception:
            vars = np.var(Y_np, axis=0, ddof=1)
            if not np.all(np.isfinite(vars)):
                vars = np.ones(I)
            covs_global = np.diag(np.maximum(vars, eps))

        if not (np.isfinite(covs_global).all() and covs_global.ndim == 2 and covs_global.shape == (I, I)):
            vars = np.var(Y_np, axis=0, ddof=1)
            if not np.all(np.isfinite(vars)):
                vars = np.ones(I)
            covs_global = np.diag(np.maximum(vars, eps))

        try:
            eigvals, eigvecs = eigh(covs_global)
            if np.min(eigvals) < eps:
                covs_global = eigvecs @ np.diag(np.maximum(eigvals, eps)) @ eigvecs.T
                covs_global = (covs_global + covs_global.T) / 2
        except Exception:
            covs_global = np.diag(np.maximum(np.diag(covs_global), eps))
        covs = np.zeros((I, I, L))
        for l in range(L):
            idx = (cluster_assignments == l)
            means[l, :] = np.mean(Y_np[idx], axis=0)
            if np.sum(idx) > 1:
                covs[:, :, l] = np.cov(Y_np[idx].T, ddof=1) if I > 1 else np.array([[np.var(Y_np[idx], ddof=1)]])
            else:
                covs[:, :, l] = covs_global.copy()

        if I == 1 and constraint in ["E", "V"]:
            if constraint == "E":
                var_shared = max(np.var(Y_np, ddof=1), eps)
                for l in range(L):
                    covs[0, 0, l] = var_shared
            elif constraint == "V":
                for l in range(L):
                    idx = (cluster_assignments == l)
                    if np.sum(idx) > 1:
                        var_l = np.var(Y_np[idx], ddof=1)
                    else:
                        var_l = np.var(Y_np, ddof=1)
                    covs[0, 0, l] = max(var_l, eps)
        else:
            if constraint == "E0" or constraint == "EE" or constraint == "EV":
                diag_vals = np.maximum(np.diag(covs_global), eps)
                diag_cov = np.diag(diag_vals)
                for l in range(L):
                    covs[:, :, l] = diag_cov
            elif constraint == "V0" or constraint == "VE" or constraint == "VV":
                for l in range(L):
                    idx = (cluster_assignments == l)
                    if np.sum(idx) > 1:
                        covc = np.cov(Y_np[idx].T, ddof=1) if I > 1 else np.array([[np.var(Y_np[idx], ddof=1)]])
                    else:
                        covc = covs_global.copy()
                    diag_vals = np.diag(covc)
                    diag_vals = np.maximum(diag_vals, eps)
                    covs[:, :, l] = np.diag(diag_vals)
        return cluster_assignments, covs, means, P_Z
    
    def apply_constraint(self, covs):
        I, I2, L = covs.shape
        assert I == I2, "Covariance matrix must be square"
        device = covs.device
        dtype = covs.dtype
        eps = self.eps
        
        covs = (covs + covs.transpose(0, 1)) * 0.5
        
        constraint = self.constraint
        constraint_used = constraint
        if I == 1:
            diag_vals = covs[0, 0, :]
            if constraint in ["E", "EE"]:
                pooled_var = torch.clamp(diag_vals.mean(), min=eps)
                covs_constrained = torch.full((1, 1, L), pooled_var, device=device, dtype=dtype)
                return covs_constrained
            
            diag_vals = torch.clamp(diag_vals, min=eps)
            covs_constrained = diag_vals.view(1, 1, L)
            return covs_constrained
        
        covs_constrained = covs.clone()
        pooled_cov = covs.mean(dim=2)
        pooled_cov = (pooled_cov + pooled_cov.T) * 0.5
        
        if constraint == "E": 
            constraint = "EE"
        elif constraint == "V":
            constraint = "VV"
        
        # Vectorized constraint application
        diag_indices = torch.arange(I, device=device)
        if constraint == "E0":
            diag_vals = torch.diagonal(pooled_cov)
            diag_vals = torch.clamp(diag_vals, min=eps)
            D = torch.diag(diag_vals)
            covs_constrained = D.unsqueeze(-1).expand(I, I, L).clone()
            return covs_constrained
        
        if constraint == "V0":
            class_diags = torch.diagonal(covs, dim1=0, dim2=1)  # (I, L)
            class_diags = torch.clamp(class_diags, min=eps)
            covs_constrained = torch.zeros_like(covs)
            covs_constrained[diag_indices, diag_indices, :] = class_diags.T
            return covs_constrained
        
        off_diag_mask = ~torch.eye(I, dtype=torch.bool, device=device)
        
        if constraint == "EE":
            covs_constrained = pooled_cov.unsqueeze(-1).expand(I, I, L).clone()
            return covs_constrained
        
        if constraint == "VV":
            return covs_constrained
        
        if constraint == "VE":
            class_diags = torch.diagonal(covs, dim1=0, dim2=1)  # (I, L)
            covs_constrained = pooled_cov.unsqueeze(-1).expand(I, I, L).clone()
            covs_constrained[diag_indices, diag_indices, :] = class_diags.T
            return covs_constrained
        
        if constraint == "EV":
            pooled_diag = torch.diagonal(pooled_cov)
            pooled_diag = torch.clamp(pooled_diag, min=eps)
            covs_constrained = covs.clone()
            covs_constrained[diag_indices, diag_indices, :] = pooled_diag.unsqueeze(1).expand(I, L)
            return covs_constrained
        
        if isinstance(constraint_used, list) or isinstance(constraint_used, tuple):
            covs_constrained = pooled_cov.unsqueeze(-1).expand(I, I, L).clone()
            for pair in constraint_used:
                if len(pair) == 2:
                    i, j = pair
                    i -= 1
                    j -= 1
                    if i < I and j < I:
                        covs_constrained[i, j, :] = pooled_cov[i, j]
                        covs_constrained[j, i, :] = pooled_cov[i, j]

        return covs_constrained
    
    def _ensure_positive_definite(self, covs):
        I, I2, L = covs.shape
        assert I == I2, "Covariance matrix must be square"
        device = covs.device
        dtype = covs.dtype
        eps = self.eps
        if self.shared_mask is None:
            self._init_shared_mask()
        non_shared_mask = ~self.shared_mask
        covs_pd = torch.zeros_like(covs)
        
        def safe_fallback(original_mat):
            diag_vals = torch.abs(torch.diag(original_mat))
            safe_diag = torch.clamp(diag_vals, min=eps * 10) 
            return torch.diag(safe_diag)
        
        if not torch.any(non_shared_mask):
            identity = torch.eye(I, device=device, dtype=dtype).unsqueeze(-1) * (eps * 10)
            return covs + identity
        
        total_jitter = 0.0
        max_jitter_per_matrix = 0.0
        
        for l in range(L):
            Sigma = covs[:, :, l]
            Sigma = 0.5 * (Sigma + Sigma.T)
            
            try:
                lmat, applied_jitter = self.optimized_cholesky_constrained(Sigma.clone(), non_shared_mask)
                candidate = lmat @ lmat.T
                if self._is_positive_definite(candidate, eps * 0.1):
                    covs_pd[:, :, l] = candidate
                    total_jitter += applied_jitter
                    max_jitter_per_matrix = max(max_jitter_per_matrix, applied_jitter)
                    continue
            except Exception as e1:
                pass 
            
            try:
                Sigma_sym = 0.5 * (Sigma + Sigma.T)
                eigvals, eigvecs = torch.linalg.eigh(Sigma_sym)
                eigvals = torch.clamp(eigvals, min=eps * 5)
                candidate = eigvecs @ torch.diag(eigvals) @ eigvecs.T
                candidate = 0.5 * (candidate + candidate.T)
                
                if self._is_positive_definite(candidate, eps * 0.5):
                    covs_pd[:, :, l] = candidate
                    continue
            except Exception as e2:
                pass
            
            covs_pd[:, :, l] = safe_fallback(Sigma)
        
        for l in range(L):
            mat = covs_pd[:, :, l]
            if not self._is_positive_definite(mat, eps * 0.1):
                covs_pd[:, :, l] = safe_fallback(covs[:, :, l])
        
        return covs_pd
    
    def _is_positive_definite(self, mat, tol=1e-8):
        try:
            _ = torch.linalg.cholesky(mat)
            return True
        except RuntimeError:
            pass
        
        try:
            min_eig = torch.linalg.eigvalsh(mat).min()
            return min_eig > tol
        except:
            return False
    
    def optimized_cholesky_constrained(self, Sigma, non_shared_mask, max_jitter=1e-3):
        device = Sigma.device
        dtype = Sigma.dtype
        I = Sigma.shape[0]
        
        if I == 1:
            val = Sigma[0, 0].item()
            safe_val = max(val, self.eps * 10)
            return torch.tensor([[math.sqrt(safe_val)]], device=device, dtype=dtype), 0.0
        
        Sigma = 0.5 * (Sigma + Sigma.T)
        
        if not torch.any(non_shared_mask):
            try:
                L = torch.linalg.cholesky(Sigma)
                return L, 0.0
            except:
                jitter = max(self.eps * 10, max_jitter)
                jittered = Sigma + torch.eye(I, device=device, dtype=dtype) * jitter
                try:
                    L = torch.linalg.cholesky(jittered)
                    return L, jitter
                except:
                    eigvals, eigvecs = torch.linalg.eigh(Sigma)
                    eigvals = torch.clamp(eigvals, min=self.eps * 5)
                    recon = eigvecs @ torch.diag(eigvals) @ eigvecs.T
                    return torch.linalg.cholesky(0.5 * (recon + recon.T)), self.eps * 5
        
        try:
            L = torch.linalg.cholesky(Sigma)
            return L, 0.0
        except:
            pass
        
        try:
            min_eig = torch.linalg.eigvalsh(Sigma).min().item()
            adjustment = max((self.eps * 5) - min_eig, self.eps * 5)
        except:
            adjustment = max_jitter * 10
        
        diag_mask = non_shared_mask & torch.eye(I, dtype=torch.bool, device=device)
        num_non_shared_diag = diag_mask.sum().item()
        
        if num_non_shared_diag > 0:
            jitter_per_diag = adjustment / num_non_shared_diag
            perturbation = torch.zeros((I, I), device=device, dtype=dtype)
            perturbation[diag_mask] = jitter_per_diag
            adjusted = Sigma + perturbation
            adjusted = 0.5 * (adjusted + adjusted.T)
            
            try:
                L = torch.linalg.cholesky(adjusted)
                return L, adjustment
            except:
                pass
        
        try:
            eigvals, eigvecs = torch.linalg.eigh(Sigma)
            new_eigvals = torch.clamp(eigvals, min=self.eps * 5)  # 严格下限
            new_mat = eigvecs @ torch.diag(new_eigvals) @ eigvecs.T
            new_mat = 0.5 * (new_mat + new_mat.T)
            
            result = Sigma.clone()
            result[non_shared_mask] = new_mat[non_shared_mask]
            result = 0.5 * (result + result.T)
            
            diag_vals = torch.diag(result)
            torch.diagonal(result).copy_(torch.clamp(diag_vals, min=self.eps * 10))
            
            return torch.linalg.cholesky(result), self.eps * 5
        except:
            diag_vals = torch.clamp(torch.diag(Sigma), min=self.eps * 10)
            safe_mat = torch.diag(diag_vals)
            return torch.linalg.cholesky(safe_mat), self.eps * 10
    
    def compute_log_pdf(self, dev, covs_constrained):
        N, L_dev, I = dev.shape
        assert L_dev == self.L, "dev second dim must equal self.L"
        device = dev.device
        eps = getattr(self, 'eps', 1e-6)
        covs_batch = covs_constrained.permute(2, 0, 1)  
        dev_batch_T = dev.permute(1, 2, 0)  # (L, I, N)
        eye = torch.eye(I, device=device).expand(L_dev, I, I)
        covs_batch = covs_batch + eps * eye
        try:
            chol_batch = torch.linalg.cholesky(covs_batch)
        except RuntimeError:
            covs_batch = covs_batch + (eps * 100) * eye
            chol_batch = torch.linalg.cholesky(covs_batch)
        diag_chol = torch.diagonal(chol_batch, dim1=1, dim2=2)
        logdet = 2.0 * torch.sum(torch.log(diag_chol), dim=1)  # (L,)
        y = torch.linalg.solve_triangular(chol_batch, dev_batch_T, upper=False)
        quad = torch.sum(y**2, dim=1)  # (L, N)
        const_term = -0.5 * I * math.log(2 * math.pi)
        log_pdfs = const_term - 0.5 * (logdet.unsqueeze(1) + quad)  # (L, N)
        log_pdfs = torch.clamp(log_pdfs, min=-1e8, max=1e8)
        return log_pdfs.t()  # (N, L)
    
    def get_Log_lik(self, P_Z, means, covs_constrained):
        N, I = self.response.shape
        L = self.L
        eps = self.eps
        
        dev = self.response.unsqueeze(1) - means.unsqueeze(0)
        log_pdfs = self.compute_log_pdf(dev, covs_constrained)
        log_pz = torch.log(P_Z.repeat(self.N, 1) + eps)
        log_joint = log_pdfs + log_pz

        row_max = torch.max(log_joint, dim=1, keepdim=True).values
        row_max = torch.where(torch.isfinite(row_max), row_max, torch.zeros_like(row_max))

        exp_term = torch.exp(log_joint - row_max)
        sum_exp = torch.sum(exp_term, dim=1, keepdim=True)
        mask_all_zero = (sum_exp < eps) 
        sum_exp = torch.where(mask_all_zero, torch.ones_like(sum_exp) * eps, sum_exp)

        log_sum_exp = row_max + torch.log(sum_exp)
        log_sum_exp = torch.where(mask_all_zero, torch.full_like(log_sum_exp, -float('inf')), log_sum_exp)
        log_sum_exp = torch.where(torch.isfinite(log_sum_exp), log_sum_exp, torch.full_like(log_sum_exp, -1e8))

        ll = torch.sum(log_sum_exp)
        return ll

    def loss(self, P_Z, means, covs_constrained):
        n_log_likelihood = -self.get_Log_lik(P_Z, means, covs_constrained)
        loss_total = n_log_likelihood
        return loss_total, -n_log_likelihood

    def get_P_Z_Xn(self):
        P_Z, means, covs = self.forward()
        P_Z = P_Z.clamp(min=self.eps)

        dev = self.response.unsqueeze(1) - means.unsqueeze(0)
        log_pdfs = self.compute_log_pdf(dev, covs)
        log_pz = torch.log(P_Z.unsqueeze(0).repeat(self.N, 1) + self.eps)

        log_joint = log_pdfs + log_pz
        log_norm = torch.logsumexp(log_joint, dim=1, keepdim=True)
        log_post = log_joint - log_norm
        P_Z_Xn = torch.exp(log_post)

        row_sums = P_Z_Xn.sum(dim=1, keepdim=True)
        row_sums = torch.where(row_sums < self.eps, torch.ones_like(row_sums), row_sums)
        P_Z_Xn = P_Z_Xn / row_sums

        return P_Z_Xn

    def get_Z(self):
        P_Z_Xn = self.get_P_Z_Xn()
        Z = torch.argmax(P_Z_Xn, dim=1) + 1
        return Z

    def get_constrained_parameters(self):
        with torch.no_grad():
            P_Z, means, covs = self.forward()
        return P_Z, means, covs

    def forward(self, x=None):
        if x is None:
            x_feat = self.response
        else:
            x_feat = x.to(self.device)

        logits = self.network(x_feat)
        logits_embed = self.embed_proj(logits)
        logits_embed = logits_embed.unsqueeze(1)
        logits_attn = self.attn_layer(logits_embed)
        logits_attn = logits_attn.squeeze(1)
        logits_mapped = self.output_proj(logits_attn)

        P_Z_Xn = F.softmax(logits_mapped, dim=1)
        P_Z = torch.sum(P_Z_Xn, dim=0) + self.eps
        P_Z = P_Z / torch.sum(P_Z)

        activated_means = self._apply_means_activation(self.means)
        activated_covs = self._apply_covs_activation(self.covs)
        
        covs_constrained = self.apply_constraint(activated_covs)
        covs_constrained = self._ensure_positive_definite(covs_constrained)
        
        return P_Z, activated_means, covs_constrained

def check_device():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return device

def simulated_annealing_optimization_LPA(LPAnet_model, response, par_ini=None, cycle=1, 
                                         initial_temperature=2000, cooling_rate=0.95, 
                                         threshold_sa=1e-4, maxiter=2000, vis=True):

    N, I  = response.shape
    
    with torch.no_grad():
        network_params = [p for p in LPAnet_model.parameters() if p is not LPAnet_model.means and p is not LPAnet_model.covs]
        means = LPAnet_model.means
        covs = LPAnet_model.covs
    
        current_network_params = [p.clone() for p in network_params]
        current_means = means.clone()
        current_covs = covs.clone()
        best_network_params = [p.clone() for p in network_params]
        best_means = means.clone()
        best_covs = covs.clone()
    
        temperature = initial_temperature
    
        P_Z, means, covs = LPAnet_model()
        best_loss, best_ll = LPAnet_model.loss(P_Z, means, covs )
        best_loss = best_loss.item()
        best_ll = best_ll.item()
    
        for iteration in range(maxiter):
            noise_network = [torch.randn_like(orig) * (orig * (temperature / initial_temperature)) for orig in current_network_params]
            new_network_params = [current + n for current, n in zip(current_network_params, noise_network)]
    
            noise_means = torch.randn_like(current_means) * (current_means * (temperature / initial_temperature))
            new_means = current_means + noise_means
            
            noise_covs = torch.randn_like(current_covs) * (current_covs * (temperature / initial_temperature))
            new_covs = current_covs + noise_covs
    
            for p, new_p in zip(network_params, new_network_params):
                p.data.copy_(new_p)
            means.data.copy_(new_means)
            covs.data.copy_(new_covs)
    
            P_Z, means, covs = LPAnet_model()
            loss, ll = LPAnet_model.loss(P_Z, means, covs)
            loss_value = loss.item()
            ll_value = ll.item()
    
            if loss_value < best_loss or np.random.rand() < np.exp((best_loss - loss_value) / temperature):
                current_network_params = [p.clone() for p in network_params]
                current_means = means.clone()
                current_covs = covs.clone()
                if loss_value < best_loss: 
                    best_loss = loss_value
                    best_ll = ll_value
                    best_network_params = [p.clone() for p in network_params]
                    best_means = means.clone()
                    best_covs = covs.clone()
    
            temperature *= cooling_rate
    
            if vis:
                print(
                        f"Iter = {iteration:{int(math.log10(abs(maxiter))) + 1}}, ", 
                        f"Loss: {-loss_value:{int(math.log10(abs(N*I))) + 3}.2f}, ", 
                        f"BIC:  {2*loss_value+np.log(N)*LPAnet_model.npar:{int(math.log10(abs(N*I))) + 3}.5f}, ", 
                        f"bets BIC: {2*best_ll+np.log(N)*LPAnet_model.npar:{int(math.log10(abs(N*I))) + 3}.5f}, ", 
                        f"Temperature: {temperature:{int(math.log10(abs(temperature))) + 1}.5f}", 
                        end='\r'
                    )
                    
            if temperature < threshold_sa:
                break
    
        for p, best_p in zip(network_params, best_network_params):
            p.data.copy_(best_p)
        means.data.copy_(best_means)
        covs.data.copy_(best_covs)
    
        if vis:
            print("\n")
        return LPAnet_model, best_network_params, best_means, best_covs

def NN_LPA(response,
           L=5,
           par_ini=None,
           constraint="VV", 
           nrep=2, 
           starts=50, 
           maxiter_wa=20, 
           vis=True,
           hidden_layers=[32],
           activation_function='tanh', 
           d_model=None, nhead=None, dim_feedforward=None, eps=1e-8, Lambda=1e-5, 
           initial_temperature=2000,
           cooling_rate=0.95,
           maxiter_sa=2000,
           threshold_sa=1e-5,
           maxiter=2000,
           maxiter_early=10,
           maxcycle=10, 
           lr = 0.025, 
           scheduler_patience = 10, 
           scheduler_factor = 0.70, 
           plot_interval=10, 
           device="CPU"):
    
    seed = 9
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        
    if isinstance(constraint, str):
        constraint = constraint
    else:
        constraint = [[int(x) for x in inner] for inner in constraint]

    if isinstance(hidden_layers, (list, tuple)):
        hidden_layers = [int(x) for x in hidden_layers]
    else:
        hidden_layers = [int(hidden_layers)]
        
    if d_model is not None:
        d_model = int(d_model)
    if nhead is not None:
        nhead = int(nhead)
    if dim_feedforward is not None:
        dim_feedforward = int(dim_feedforward)
    
    N, I = response.shape
    nrep = int(nrep)
    L = int(L)
    maxiter_sa = int(maxiter_sa)
    maxiter = int(maxiter)
    maxiter_early = int(maxiter_early)
    scheduler_patience = int(scheduler_patience)
    maxcycle = int(maxcycle)
    starts = int(starts)
    maxiter_wa = int(maxiter_wa)
    if device == "GPU":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device("cpu")

    response_tensor = torch.tensor(response, dtype=torch.float32).to(device)
    
    best_models = []
    
    if starts > 0:
        for s in range(starts):
            current_seed = seed + s
            torch.manual_seed(current_seed)
            np.random.seed(current_seed)
            random.seed(current_seed)
            if torch.cuda.is_available():
                torch.cuda.manual_seed_all(current_seed)
            
            LPAnet_warmup = LPAnet(response=response_tensor,
                                  L=L,
                                  par_ini=par_ini,
                                  constraint=constraint,
                                  hidden_layers=hidden_layers,
                                  activation_function=activation_function, 
                                  d_model=d_model, nhead=nhead, 
                                  dim_feedforward=dim_feedforward, eps=eps)
            LPAnet_warmup.to(device)
            
            network_params = [p for p in LPAnet_warmup.parameters() 
                             if p is not LPAnet_warmup.means and p is not LPAnet_warmup.covs]
            optimizer = torch.optim.AdamW(
                [
                    {'params': network_params, 'lr': lr},
                    {'params': LPAnet_warmup.means, 'lr': lr},
                    {'params': LPAnet_warmup.covs, 'lr': lr}
                ],
                weight_decay=Lambda
            )
            
            best_ll = -float('inf')
            for epoch in range(maxiter_wa):
                optimizer.zero_grad()
                P_Z, means, covs = LPAnet_warmup()
                loss, loss_ll = LPAnet_warmup.loss(P_Z, means, covs)
                loss.backward()
                optimizer.step()
                
                current_ll = loss_ll.item()
                if current_ll > best_ll:
                    best_ll = current_ll
                    best_state = {k: v.cpu() for k, v in LPAnet_warmup.state_dict().items()}
            
            if len(best_models) < nrep:
                best_models.append((best_ll, best_state))
            else:
                min_idx = min(range(len(best_models)), key=lambda i: best_models[i][0])
                if best_ll > best_models[min_idx][0]:
                    best_models[min_idx] = (best_ll, best_state)
            
            if vis:
                current_min_ll = min(model[0] for model in best_models)
                print(f"Warm {s+1}/{starts} | Best log-likelihood: {best_ll:.4f} | "
                      f"Min in top {nrep}: {current_min_ll:.4f}", end='\r')
            
            del LPAnet_warmup, optimizer
            torch.cuda.empty_cache()
        
        best_models.sort(key=lambda x: x[0], reverse=True)
    else:
        best_models = []
    
    if starts > 0:
        print("\n")

    if nrep <= 5:
        colors = plt.cm.tab10(np.linspace(0, 1, nrep))
    elif nrep <= 10:
        colors = plt.cm.tab20(np.linspace(0, 1, nrep))
    else:
        colors = plt.cm.viridis(np.linspace(0.1, 0.9, nrep))
    
    Log_Lik_nrep = []
    all_results = []
    best_overall_loss = float('inf')
    best_rep_index = -1

    global_best_ll = float('inf')
    global_best_loss = float('inf')
    global_best_step = None
    
    if vis:
        plt.ion()
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title('Loss Curves')
        ax.set_xlabel('Iterations')
        ax.set_ylabel('Loss')
        plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)
        
        line_trainings = []
        line_annealings = []
        global_best_loss = float('inf')
        global_best_step = None
        global_best_marker = ax.plot([], [], '*', color='gold', markersize=15, 
                                     markeredgecolor='k', label='Global Best Point')[0]
        
        for rep in range(nrep):
            color = colors[rep]
            line_trainings.append(ax.plot([], [], '.', color=color, alpha=0.5, 
                                         markersize=3,  # 点的大小
                                         label=f'Training Rep {rep+1}')[0])
            line_annealings.append(ax.plot([], [], 'o', color='gray', alpha=0.6, 
                                          markersize=6)[0])
        
        ax.legend().set_visible(False)
        fig.canvas.draw()
        plt.pause(0.001)

    for rep in range(nrep):
        current_seed = seed + starts + rep
        torch.manual_seed(current_seed)
        np.random.seed(current_seed)
        random.seed(current_seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(current_seed)
        
        LPAnet_model = LPAnet(response=response_tensor,
                             L=L,
                             par_ini=par_ini,
                             constraint=constraint,
                             hidden_layers=hidden_layers,
                             activation_function=activation_function, 
                             d_model=d_model, nhead=nhead, 
                             dim_feedforward=dim_feedforward, eps=eps)
        LPAnet_model.to(device)
        
        if starts > 0 and rep < len(best_models):
            LPAnet_model.load_state_dict(best_models[rep][1])
        
        log_records = []
        local_step = 0

        P_Z, means, covs = LPAnet_model()
        init_loss, init_ll = LPAnet_model.loss(P_Z, means, covs)
        init_loss = init_loss.item()
        init_ll = init_ll.item()
        log_records.append((local_step, init_loss))

        best_model_state = {k: v.clone() for k, v in LPAnet_model.state_dict().items()}
        best_means = means.clone().detach()
        best_covs = covs.clone().detach()
        best_P_Z = P_Z.clone().detach()
        best_loss = init_loss
        best_ll = init_ll
        best_step = 0

        if vis:
            train_x = [0]
            train_y = [init_loss]
            anneal_x = []
            anneal_y = []

        improved = True
        cycle = 0
        
        network_params = [p for p in LPAnet_model.parameters() 
                          if p is not LPAnet_model.means and p is not LPAnet_model.covs]
        means = LPAnet_model.means
        covs = LPAnet_model.covs
        
        optimizer = torch.optim.AdamW(
            [
                {'params': network_params, 'lr': lr},
                {'params': means, 'lr': lr},
                {'params': covs, 'lr': lr}
            ],
            weight_decay=Lambda
        )
        scheduler = ReduceLROnPlateau(optimizer, 'min', patience=scheduler_patience, factor=scheduler_factor)
        
        while improved and cycle < maxcycle:
            cycle += 1

            patience = 0
            for epoch in range(maxiter):
                optimizer.zero_grad()
                P_Z, means, covs = LPAnet_model()
                loss, loss_ll = LPAnet_model.loss(P_Z, means, covs)
                loss.backward()
                optimizer.step()
                scheduler.step(loss)

                cur_loss = loss.item()
                cur_ll = loss_ll.item()
                local_step += 1
                log_records.append((local_step, cur_loss))

                if vis:
                    train_x.append(local_step)
                    train_y.append(cur_loss)
                    line_trainings[rep].set_data(train_x, train_y)
                    if local_step % plot_interval == 0:
                        ax.relim()
                        ax.autoscale_view()
                        fig.canvas.draw()
                        plt.pause(0.001)

                if cur_loss < best_loss:
                    best_loss = cur_loss
                    best_ll = cur_ll
                    best_step = local_step
                    best_means = means.clone().detach()
                    best_covs = covs.clone().detach()
                    best_P_Z = P_Z.clone().detach()
                    best_model_state = {k: v.clone() for k, v in LPAnet_model.state_dict().items()}
                    
                    if cur_loss < global_best_loss:
                        global_best_loss = cur_loss
                        global_best_step = local_step
                        global_best_ll = cur_ll
                        if vis and global_best_marker is not None:
                            global_best_marker.set_data([global_best_step], [global_best_loss])
                        
                    patience = 0
                else:
                    patience += 1
                    if patience % scheduler_patience == 0:
                        LPAnet_model.load_state_dict(best_model_state)

                if vis and local_step % 50 == 0:
                    print(
                        f"Rep  {rep+1}/{nrep} | Iter = {local_step:{int(math.log10(abs(maxcycle*maxiter))) + 1}}, ", 
                        f"Loss: {cur_loss:{int(math.log10(abs(N*I))) + 3}.2f}, ", 
                        f"BIC: {-2*cur_ll+np.log(N)*LPAnet_model.npar:{int(math.log10(abs(N*I))) + 3}.5f}, ", 
                        f"Best BIC: {-2*global_best_ll+np.log(N)*LPAnet_model.npar:{int(math.log10(abs(N*I))) + 3}.5f}, ", 
                        f"Patience: {patience:3d}, Cycle: {cycle:{int(math.log10(abs(maxcycle))) + 1}}", 
                        end='\r'
                    )

                if patience >= maxiter_early:
                    break

            if best_model_state is not None:
                LPAnet_model.load_state_dict(best_model_state)

            improved = False
            LPAnet_model, sa_params, sa_means, sa_covs = simulated_annealing_optimization_LPA(
                LPAnet_model, response_tensor,
                par_ini=par_ini,
                cycle=cycle,
                initial_temperature=initial_temperature,
                cooling_rate=cooling_rate,
                threshold_sa=threshold_sa,
                maxiter=maxiter_sa,
                vis=False
            )
            P_Z, means, covs = LPAnet_model()
            sa_loss, sa_ll = LPAnet_model.loss(P_Z, means, covs)
            sa_loss = sa_loss.item()
            sa_ll = sa_ll.item()
            local_step += 1
            log_records.append((local_step, sa_loss))

            if vis:
                anneal_x.append(local_step)
                anneal_y.append(sa_loss)
                line_annealings[rep].set_data(anneal_x, anneal_y)
                ax.relim()
                ax.autoscale_view()
                fig.canvas.draw()
                plt.pause(0.001)

            if sa_loss < best_loss:
                best_loss = sa_loss
                best_ll = sa_ll
                best_step = local_step
                best_means = means.clone().detach()
                best_covs = covs.clone().detach()
                best_P_Z = P_Z.clone().detach()
                best_model_state = {k: v.clone() for k, v in LPAnet_model.state_dict().items()}
                
                if sa_loss < global_best_loss:
                    global_best_loss = sa_loss
                    global_best_step = local_step
                    global_best_ll = sa_ll
                    if vis and global_best_marker is not None:
                        global_best_marker.set_data([global_best_step], [global_best_loss])
                
                improved = True
            else:
                if best_model_state is not None:
                    LPAnet_model.load_state_dict(best_model_state)
        
        if best_model_state is not None:
            LPAnet_model.load_state_dict(best_model_state)
        P_Z, means, covs = LPAnet_model()
        loss, ll = LPAnet_model.loss(P_Z, means, covs)
        ll = ll.item()
        rep_result = {
            'params': {
                'means': means.clone().detach().cpu().numpy(),
                'covs': covs.clone().detach().cpu().numpy(),
                'P.Z': P_Z.clone().detach().cpu().numpy()
            }, 
            'model': LPAnet_model, 
            'npar': LPAnet_model.npar, 
            'Log.Lik': ll, 
            'AIC': -2*ll + 2*LPAnet_model.npar, 
            'BIC': -2*ll + np.log(N)*LPAnet_model.npar, 
            'P.Z.Xn': LPAnet_model.get_P_Z_Xn().cpu().detach().numpy(),
            'P.Z': P_Z.clone().detach().cpu().numpy(), 
            'Z': LPAnet_model.get_Z().clone().detach().cpu().numpy(), 
            'Log.Lik.history': log_records
        }
        all_results.append(rep_result)
        
        if best_loss < best_overall_loss:
            best_overall_loss = best_loss
            best_rep_index = rep
        
        Log_Lik_nrep.append(-2*ll + np.log(N)*LPAnet_model.npar)
        
    final_result = all_results[best_rep_index]
    final_result['Log.Lik.nrep'] = Log_Lik_nrep
    
    if vis:
        best_color = colors[best_rep_index]
        line_trainings[best_rep_index].set_color("red")
        line_trainings[best_rep_index].set_markersize(8)
        line_trainings[best_rep_index].set_zorder(10)
        line_annealings[best_rep_index].set_markersize(6)
        line_annealings[best_rep_index].set_alpha(1.0)

        ax.set_title(f'Loss Curves')

        annealing_legend = ax.plot([], [], 'o', color='gray', alpha=1.0, 
                                  markersize=6, label='Annealing')[0]
        
        if global_best_step is not None and global_best_loss is not None:
            global_best_marker.set_data([global_best_step], [global_best_loss])
            global_best_marker.set_zorder(10)
        
        handles = [line_trainings[best_rep_index], global_best_marker, annealing_legend]
        labels = ['Best Nrep', 'Global Best', 'Annealing']
        ax.legend(handles=handles, labels=labels, 
                 loc='upper right', fontsize='small')
        
        ax.relim()
        ax.autoscale_view()
        fig.canvas.draw()
        plt.pause(0.001)
        plt.ioff()
        plt.show()
        
        print("\n")

    return final_result
