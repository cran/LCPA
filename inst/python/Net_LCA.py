import torch
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
import six
from torch.distributions import Dirichlet
import random
import os
import sys

class LCAnet(nn.Module):
    def __init__(self, response, L=5, par_ini=None, hidden_layers=[32, 32], 
                 activation_function='tanh', use_attention=True, 
                 d_model=None, nhead=None, dim_feedforward=None, eps=1e-8):
        super(LCAnet, self).__init__()
        self.L = L
        self.eps = eps
        self.use_attention = use_attention
        
        device = response.device if hasattr(response, 'device') else torch.device('cpu')
        adjust_response_obj = self.adjust_response(response)
        self.response = adjust_response_obj["response"].to(device)
        self.poly_orig = adjust_response_obj["poly_orig"]
        self.poly_value = torch.tensor(adjust_response_obj["poly_value"], dtype=torch.long, device=device)
        self.poly_max = adjust_response_obj["poly_max"]
        self.device = self.response.device
        self.N, self.I = self.response.shape
        
        self.response_arange = torch.arange(self.poly_max, dtype=torch.long, device=self.device).view(1, 1, self.poly_max).expand(self.N, self.I, -1)
        self.response_hot = F.one_hot(self.response.long(), num_classes=self.poly_max).float()
        self.register_buffer('response_hot_buf', self.response_hot)

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

        if par_ini is not None and not isinstance(par_ini, (str, six.string_types)):
            par_ini["par"] = np.asarray(par_ini["par"], dtype=np.float32)
            self.par_mask_np = np.isnan(par_ini["par"])  # keep for shape only
            par_ini["par"] = np.clip(par_ini["par"], self.eps, 1 - self.eps)
            par = torch.tensor(par_ini["par"], dtype=torch.float32, device=self.device)
            self.par_mask = torch.tensor(self.par_mask_np, device=self.device)  # GPU tensor
            self.par = nn.Parameter(self.softmax_inverse(par, self.poly_value, method='zero_sum'))
            self.P_Z = torch.from_numpy(par_ini["P.Z"]).float().to(self.device)

        elif par_ini == "kmeans": 
            par, P_Z = self.kmeans_classify(self.response, self.L, self.poly_max, self.poly_value, nstart=1)
            self.par_mask_np = np.isnan(par.cpu().numpy())  # only for shape/debug
            self.par_mask = torch.isnan(par)  # GPU mask
            par = torch.nan_to_num(par, nan=float('-inf'))
            par = self.softmax_inverse(par, self.poly_value, method='zero_sum')
            self.par = nn.Parameter(par)
            self.P_Z = torch.softmax(P_Z, dim=0)

        elif par_ini == "random": 
            par_tensor, P_Z = self._random_init_params_torch()
            self.par_mask = torch.isnan(par_tensor)
            par_tensor = torch.nan_to_num(par_tensor, nan=float('-inf'))
            self.par = nn.Parameter(self.softmax_inverse(par_tensor, self.poly_value, method='zero_sum'))
            self.P_Z = P_Z

        else: 
            rand_num = random.random()
            if rand_num > 0.5:
                par_tensor, P_Z = self._random_init_params_torch()
                self.par_mask = torch.isnan(par_tensor)
                par_tensor = torch.nan_to_num(par_tensor, nan=float('-inf'))
                self.par = nn.Parameter(self.softmax_inverse(par_tensor, self.poly_value, method='zero_sum'))
                self.P_Z = P_Z
            else:
                par, P_Z = self.kmeans_classify(
                    self.response, 
                    self.L,
                    self.poly_max, 
                    self.poly_value, 
                    nstart=1
                )
                self.par_mask = torch.isnan(par)
                par = torch.nan_to_num(par, nan=float('-inf'))
                par = self.softmax_inverse(par, self.poly_value, method='zero_sum')
                self.par = nn.Parameter(par)
                self.P_Z = torch.softmax(P_Z, dim=0)

        layers = []
        input_dim = self.input_dim
        for hidden_units in hidden_layers:
            layers.append(nn.Linear(input_dim, hidden_units, dtype=torch.float32))
            layers.append(nn.BatchNorm1d(hidden_units, dtype=torch.float32))
            layers.append(activation)
            input_dim = hidden_units

        final_layer = nn.Linear(input_dim, self.L, dtype=torch.float32)
        layers.append(final_layer)
        self.network = nn.Sequential(*layers)
        
        if self.use_attention:
            if d_model is None: 
                d_model = 8
            if nhead is None: 
                nhead = 2
            if d_model % nhead != 0:
                d_model += (nhead - d_model % nhead)
            if dim_feedforward is None: 
                dim_feedforward = max(d_model, 16)

            self.embed_proj = nn.Linear(self.L, d_model)
            encoder_layer = TransformerEncoderLayer(d_model=d_model, nhead=nhead,
                                                    dim_feedforward=dim_feedforward, batch_first=True)
            self.attn_layer = TransformerEncoder(encoder_layer, num_layers=1, enable_nested_tensor=False)
            self.output_proj = nn.Linear(d_model, self.L)

        self.to(self.device)
        self.npar = self._compute_npar()
         
    def _random_init_params_torch(self):
        par_array = torch.full((self.L, self.I, self.poly_max), float('nan'), device=self.device)
        for i in range(self.I):
            k = self.poly_value[i]
            if k <= 0:
                continue
            alpha = torch.ones(k, device=self.device) * 3.0
            dist = Dirichlet(alpha)
            samples = dist.sample((self.L,))
            par_array[:, i, :k] = samples
            
        alpha = torch.ones(self.L, device=self.device) * 3.0
        dist = Dirichlet(alpha)
        P_Z = dist.sample()
        return par_array, P_Z
    
    def _compute_npar(self):
        pv = self.poly_value.cpu().numpy() if isinstance(self.poly_value, torch.Tensor) else self.poly_value
        npar = np.sum(pv * self.L - 1) + self.L - 1
        return int(npar)
    
    @staticmethod
    def adjust_response(response):
        device = response.device if hasattr(response, 'device') else torch.device('cpu')
        response_cpu = response.cpu() if hasattr(response, 'cpu') else response
        N, I = response_cpu.shape

        poly_value = np.zeros(I, dtype=int)
        unique_list = []

        max_k = 0
        for i in range(I):
            unique_vals, inverse_indices = torch.unique(response_cpu[:, i], return_inverse=True)
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
            response_adjusted[:, i] = inverse_indices.to(device)
        
        return {
            'poly_orig': poly_orig,
            'poly_value': poly_value,
            'poly_max': poly_max,
            'response': response_adjusted
        }

    @staticmethod
    def kmeans_classify(Y, L, poly_max, poly_value, nstart=100):
        N, I = Y.shape
        Y_np = Y.detach().cpu().numpy().astype(int)
        mean = np.mean(Y_np, axis=0)
        std = np.std(Y_np, axis=0)
        std = np.where(std == 0, 1.0, std)
        Y_normalized = (Y_np - mean) / std
        
        km = KMeans(
            init='k-means++',
            n_clusters=int(L),
            max_iter=500,
            n_init=nstart,
            algorithm='lloyd',
            verbose=0
        )
        cluster_labels = km.fit_predict(Y_normalized)
        
        P_Z = np.bincount(cluster_labels, minlength=L) / N
        P_Z = torch.from_numpy(P_Z).float().to(Y.device)
        
        par = np.full((L, I, poly_max), np.nan)
        for l in range(L):
            l_posi = np.where(cluster_labels == l)[0]
            Y_l = Y_np[l_posi, :]
            for i in range(I):
                unique_vals, counts = np.unique(Y_l[:, i], return_counts=True)
                par[l, i, unique_vals] = counts / sum(counts)
                par[l, i, :] = par[l, i, :] + 1e-4
                
        for i in range(I):
            par[:, i, 0:poly_value[i]] = np.nan_to_num(par[:, i, 0:poly_value[i]], nan=1e-4)
        par = torch.from_numpy(par).float().to(Y.device)
        return par, P_Z

    def softmax_inverse(self, par, poly_value, method='zero_sum'):
        L, I, poly_max = par.shape
        device = par.device
        
        arange_k = torch.arange(poly_max, device=device).unsqueeze(0)  # (1, poly_max)
        mask = arange_k < poly_value.unsqueeze(1)  # (I, poly_max)
        mask = mask.unsqueeze(0).expand(L, -1, -1)  # (L, I, poly_max)
        
        log_par = torch.log(torch.clamp(par, min=1e-12))
        
        valid_log_par = log_par * mask.float()
        count = mask.sum(dim=2, keepdim=True).clamp(min=1.0)  # (L, I, 1)
        
        if method == 'zero_sum':
            mean_log = valid_log_par.sum(dim=2, keepdim=True) / count
            results = log_par - mean_log
        elif method == 'max_zero':
            safe_log = torch.where(mask, log_par, torch.tensor(-float('inf'), device=device))
            max_log = safe_log.max(dim=2, keepdim=True).values
            results = log_par - max_log
        elif method == 'first_zero':
            first_log = log_par[:, :, 0].unsqueeze(2)
            results = log_par - first_log
        else:
            results = log_par
            
        return results

    def LCA(self, Z):
        idx = Z - 1  # (N,)
        par_selected = self.par[idx]  # (N, I, poly_max)
        mask_selected = ~self.par_mask[idx]  # (N, I, poly_max)
        P_LCA = torch.softmax(par_selected, dim=2) * mask_selected.float()
        return P_LCA

    def get_P_Z_Xn(self):
        P_Z, par = self.forward()
        eps = self.eps
        p = par.unsqueeze(0)  # (1, L, I, poly_max)
        rt = self.response_hot_buf.unsqueeze(1)  # (N, 1, I, poly_max)
        probs = (p * rt).sum(dim=3)  # (N, L, I)
        log_probs = torch.log(probs + eps)  # (N, L, I)
        log_pxz = torch.sum(log_probs, dim=2)  # (N, L)
        log_pz = torch.log(P_Z + eps)  # (L,)
        log_pz = log_pz.repeat(self.N, 1)  # (N, L)
        log_joint = log_pz + log_pxz  # (N, L)
        log_joint_max = torch.max(log_joint, dim=1, keepdim=True).values
        log_joint = log_joint - log_joint_max
        P_Z_Xn = F.softmax(log_joint, dim=1)  # (N, L)
        return P_Z_Xn
    
    def get_Z(self):
        P_Z_Xn = self.get_P_Z_Xn()
        Z = torch.argmax(P_Z_Xn, dim=1) + 1
        return Z
    
    def log_lik(self, P_Z, par):
        eps = self.eps
        p = par.unsqueeze(0)  # (1, L, I, poly_max)
        rt = self.response_hot_buf.unsqueeze(1)  # (N, 1, I, poly_max)
        probs = (p * rt).sum(dim=3)  # (N, L, I)
        log_probs = torch.log(probs + eps)  # (N, L, I)
        log_pxz = torch.sum(log_probs, dim=2)  # (N, L)
        log_pz = torch.log(P_Z + eps)  # (L,)
        log_pz = log_pz.repeat(self.N, 1)  # (N, L)
        log_joint = log_pz + log_pxz  # (N, L)
        log_marginal_per_sample = torch.logsumexp(log_joint, dim=1)  # (N,)
        log_likelihood = torch.sum(log_marginal_per_sample)  # scalar
        return log_likelihood
    
    def loss(self, P_Z, par):
        n_log_likelihood = - self.log_lik(P_Z, par)
        total_loss = n_log_likelihood
        return total_loss, -n_log_likelihood

    def forward(self, x=None):
        if x is None:
            x_feat = self.response
        else:
            x_feat = x.to(self.device)

        logits = self.network(x_feat)
        
        if self.use_attention:
            logits_embed = self.embed_proj(logits)
            logits_attn = self.attn_layer(logits_embed.unsqueeze(1)).squeeze(1)
            logits_mapped = self.output_proj(logits_attn)
        else:
            logits_mapped = logits
        
        P_Z_Xn = F.softmax(logits_mapped, dim=1)
        P_Z = torch.sum(P_Z_Xn, dim=0, keepdim=True) + self.eps
        self.P_Z = P_Z / torch.sum(P_Z)

        par_temp = torch.where(self.par_mask, torch.tensor(float('-inf'), device=self.device), self.par)
        par = torch.softmax(par_temp, dim=2)
        
        return self.P_Z, par

def check_device():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return device

def simulated_annealing_optimization_LCA(LCAnet_model, response, par_ini=None, cycle=1,
                                         initial_temperature=2000, cooling_rate=0.95, 
                                         threshold_sa=1e-4, maxiter=2000, vis=True):

    N, I  = response.shape
    
    with torch.no_grad():
        network_params = [p for p in LCAnet_model.parameters() if p is not LCAnet_model.par]
        par = LCAnet_model.par
    
        current_network_params = [p.clone() for p in network_params]
        current_par = par.clone()
        best_network_params = [p.clone() for p in network_params]
        best_par = par.clone()
    
        temperature = initial_temperature
    
        P_Z, par = LCAnet_model()
        best_loss, best_ll = LCAnet_model.loss(P_Z, par)
        best_loss = best_loss.item()
        best_ll = best_ll.item()
    
        for iteration in range(maxiter):
            noise_network = [torch.randn_like(orig) * (orig * (temperature / initial_temperature)) for orig in current_network_params]
            new_network_params = [current + n for current, n in zip(current_network_params, noise_network)]
    
            noise_par = torch.randn_like(current_par) * (current_par * (temperature / initial_temperature))
            new_par = current_par + noise_par
    
            for p, new_p in zip(network_params, new_network_params):
                p.data.copy_(new_p)
            par.data.copy_(new_par)
    
            P_Z, par = LCAnet_model()
            loss, loss_ll = LCAnet_model.loss(P_Z, par)
            loss_value = loss.item()
    
            if loss_value < best_loss or np.random.rand() < np.exp((best_loss - loss_value) / temperature):
                current_network_params = [p.clone() for p in network_params]
                current_par = par.clone()
                if loss_value < best_loss: 
                    best_loss = loss_value
                    best_network_params = [p.clone() for p in network_params]
                    best_par = par.clone()
    
            temperature *= cooling_rate
    
            if vis:
                print(
                        f"Iter = {iteration:{int(math.log10(abs(maxiter))) + 1}}, ", 
                        f"Loss: {loss_value:{int(math.log10(abs(N*I))) + 1}.2f}, ", 
                        f"Temperature: {temperature:{int(math.log10(abs(temperature))) + 1}.5f}", 
                        end='\r'
                    )
                    
            if temperature < threshold_sa:
                break
    
        for p, best_p in zip(network_params, best_network_params):
            p.data.copy_(best_p)
        par.data.copy_(best_par)
    
        if vis:
            print("\n")
        return LCAnet_model, best_network_params, best_par

def NN_LCA(response,
           L=5,
           par_ini=None,
           nrep=2, 
           starts=50, 
           maxiter_wa=20, 
           vis=True,
           hidden_layers=[32],
           activation_function='tanh', use_attention=True, 
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
    
    os.environ["TORCH_COMPILE_DISABLE"] = "1"
    os.environ["TORCHDYNAMO_DISABLE"] = "1"
    os.environ["TORCHINDUCTOR_DISABLE_REPRODUCIBILITY"] = "1"
    os.environ["TORCHINDUCTOR_MAX_AUTOTUNE"] = "1"
    os.environ["TORCHINDUCTOR_MAX_AUTOTUNE_GEMM"] = "1"
    os.environ["TORCH_COMPILE_DEBUG"] = "0"
    
    is_windows = sys.platform.startswith('win')
    use_compile = not is_windows
    
    seed = 56756765
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)

    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
        torch.set_float32_matmul_precision('high')
        torch.backends.cudnn.benchmark = True
        torch.backends.cudnn.enabled = True
        torch.backends.cudnn.deterministic = False
        
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
            
            LCAnet_warmup = LCAnet(response=response_tensor,
                                  L=L,
                                  par_ini=par_ini,
                                  hidden_layers=hidden_layers,
                                  activation_function=activation_function, use_attention=use_attention, 
                                  d_model=d_model, nhead=nhead, 
                                  dim_feedforward=dim_feedforward, eps=eps)
            LCAnet_warmup.to(device)
            
            if device.type == 'cuda' and use_compile:
                LCAnet_warmup = torch.compile(
                    LCAnet_warmup,
                    mode='reduce-overhead',
                    dynamic=True
                )
            
            network_params = [p for p in LCAnet_warmup.parameters() if p is not LCAnet_warmup.par]
            optimizer = torch.optim.AdamW(
                [
                    {'params': network_params, 'lr': lr},
                    {'params': LCAnet_warmup.par, 'lr': lr}
                ],
                weight_decay=Lambda
            )
            
            scaler = torch.amp.GradScaler('cuda', enabled=(device.type == 'cuda'))
            
            best_ll = -float('inf')
            for epoch in range(maxiter_wa):
                optimizer.zero_grad()
                with torch.amp.autocast(device_type='cuda', dtype=torch.bfloat16, enabled=(device.type == 'cuda')):
                    P_Z, par = LCAnet_warmup()
                    loss, loss_ll = LCAnet_warmup.loss(P_Z, par)
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()
                
                current_ll = loss_ll.item()
                if current_ll > best_ll:
                    best_ll = current_ll
                    best_state = {k: v.cpu() for k, v in LCAnet_warmup.state_dict().items()}
            
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
            
            del LCAnet_warmup, optimizer
            torch.cuda.empty_cache()
        
        best_models.sort(key=lambda x: x[0], reverse=True)
    else:
        best_models = []
        
    if vis & starts > 0:
       print("\n\n")

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
        global_best_marker = ax.plot([], [], '*', color='gold', markersize=15, 
                                     markeredgecolor='k', label='Global Best Point')[0]
        
        for rep in range(nrep):
            color = colors[rep]
            line_trainings.append(ax.plot([], [], '.', color=color, alpha=0.5, 
                                         markersize=3,
                                         label=f'Training Rep {rep+1}', )[0])
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
        
        LCAnet_model = LCAnet(response=response_tensor,
                              L=L,
                              par_ini=par_ini,
                              hidden_layers=hidden_layers,
                              activation_function=activation_function, use_attention=use_attention, 
                              d_model=d_model, nhead=nhead, 
                              dim_feedforward=dim_feedforward, eps=eps)
        LCAnet_model.to(device)
        
        if device.type == 'cuda' and use_compile:
            LCAnet_model = torch.compile(
                LCAnet_model,
                mode='reduce-overhead',
                dynamic=True
            )
        
        if starts > 0 and rep < len(best_models):
            LCAnet_model.load_state_dict(best_models[rep][1])
        
        log_records = []
        local_step = 0

        P_Z, par = LCAnet_model()
        init_loss, init_ll = LCAnet_model.loss(P_Z, par)
        init_loss = init_loss.item()
        init_ll = init_ll.item()
        log_records.append((local_step, init_loss))

        best_model_state = {k: v.clone() for k, v in LCAnet_model.state_dict().items()}
        best_par = par.clone().detach()
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
        
        network_params = [p for p in LCAnet_model.parameters() if p is not LCAnet_model.par]
        par = LCAnet_model.par
        
        optimizer = torch.optim.AdamW(
            [
                {'params': network_params, 'lr': lr},
                {'params': par, 'lr': lr}
            ],
            weight_decay=Lambda
        )
        scheduler = ReduceLROnPlateau(optimizer, 'min', patience=scheduler_patience, factor=scheduler_factor)
        
        scaler = torch.amp.GradScaler('cuda', enabled=(device.type == 'cuda'))
        
        while improved and cycle < maxcycle:
            cycle += 1

            patience = 0
            for epoch in range(maxiter):
                optimizer.zero_grad()
                with torch.amp.autocast(device_type='cuda', dtype=torch.bfloat16, enabled=(device.type == 'cuda')):
                    P_Z, par = LCAnet_model()
                    loss, loss_ll = LCAnet_model.loss(P_Z, par)
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()
                scheduler.step(loss.item())

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
                    best_par = par.clone().detach()
                    best_P_Z = P_Z.clone().detach()
                    best_model_state = {k: v.clone() for k, v in LCAnet_model.state_dict().items()}
                    
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
                        LCAnet_model.load_state_dict(best_model_state)

                if vis and local_step % 50 == 0:
                    print(
                        f"Rep  {rep+1}/{nrep} | Iter = {local_step:{int(math.log10(abs(maxcycle*maxiter))) + 1}}, ", 
                        f"Loss: {cur_loss:{int(math.log10(abs(N*I))) + 1}.2f}, ", 
                        f"BIC: {-2*cur_ll+np.log(N)*LCAnet_model.npar:{int(math.log10(abs(N*I))) + 3}.5f}, ", 
                        f"Best BIC: {-2*global_best_ll+np.log(N)*LCAnet_model.npar:{int(math.log10(abs(N*I))) + 3}.5f}, ", 
                        f"Patience: {patience:3d}, Cycle: {cycle:{int(math.log10(abs(maxcycle))) + 1}}", 
                        end='\r'
                    )

                if patience >= maxiter_early:
                    break

            if best_model_state is not None:
                LCAnet_model.load_state_dict(best_model_state)

            improved = False
            LCAnet_model, sa_params, sa_par = simulated_annealing_optimization_LCA(
                LCAnet_model, response_tensor,
                par_ini=par_ini,
                cycle=cycle,
                initial_temperature=initial_temperature,
                cooling_rate=cooling_rate,
                threshold_sa=threshold_sa,
                maxiter=maxiter_sa,
                vis=False
            )
            P_Z, par = LCAnet_model()
            sa_loss, sa_ll = LCAnet_model.loss(P_Z, par)
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
                best_par = par.clone().detach()
                best_P_Z = P_Z.clone().detach()
                best_model_state = {k: v.clone() for k, v in LCAnet_model.state_dict().items()}
                
                if sa_loss < global_best_loss:
                    global_best_loss = sa_loss
                    global_best_step = local_step
                    global_best_ll = sa_ll
                    if vis and global_best_marker is not None:
                        global_best_marker.set_data([global_best_step], [global_best_loss])
                
                improved = True
            else:
                if best_model_state is not None:
                    LCAnet_model.load_state_dict(best_model_state)
            
        if best_model_state is not None:
            LCAnet_model.load_state_dict(best_model_state)
        P_Z, par = LCAnet_model()
        loss, ll = LCAnet_model.loss(P_Z, par)
        ll = ll.item()
        rep_result = {
            'params': {
                'par': par.clone().detach().cpu().numpy(),
                'P.Z': P_Z.clone().detach().cpu().numpy()
            }, 
            'model': LCAnet_model, 
            'npar': LCAnet_model.npar, 
            'Log.Lik': ll, 
            'AIC': -2*ll + 2*LCAnet_model.npar, 
            'BIC': -2*ll + np.log(N)*LCAnet_model.npar, 
            'best_BIC': -2*ll + np.log(N)*LCAnet_model.npar, 
            'P.Z.Xn': LCAnet_model.get_P_Z_Xn().cpu().detach().numpy(),
            'P.Z': P_Z.clone().detach().cpu().numpy(), 
            'Z': LCAnet_model.get_Z().clone().detach().cpu().numpy(), 
            'Log.Lik.history': log_records
        }
        all_results.append(rep_result)
        
        if best_loss < best_overall_loss:
            best_overall_loss = best_loss
            best_rep_index = rep
        
        Log_Lik_nrep.append(-2*ll + np.log(N)*LCAnet_model.npar)
            
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

        annealing_legend = ax.plot([], [], 'o', color='gray', alpha=0.8, 
                                  markersize=6, label='Annealing')[0]
        
        if global_best_step is not None and global_best_loss is not None:
            global_best_marker.set_data([global_best_step], [global_best_loss])
            global_best_marker.set_zorder(10)
        
        handles = [global_best_marker, line_trainings[best_rep_index], annealing_legend]
        labels = ['Global Best', 'Best Nrep', 'Annealing']
        ax.legend(handles=handles, labels=labels, 
                 loc='upper right', fontsize='large')
        
        ax.relim()
        ax.autoscale_view()
        fig.canvas.draw()
        plt.pause(0.001)
        plt.ioff()
        plt.show()
        
        print("\n")

    return final_result
