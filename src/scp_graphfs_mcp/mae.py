"""Masked-autoencoder transfer (EXPLORATION.md Amendment 15, workstream 4): the nonlinear version of reference projection.

A denoising/masked autoencoder is pretrained on the external reference (cells/samples x target proteins, z-scored per
protein within each source dataset, NaN = unmeasured). Inputs are (z * observed, observed-mask); during training a random
fraction of the *observed* entries is hidden from the encoder and the loss is computed on exactly those entries, so the
network learns to predict unmeasured proteins from measured ones (BERT-style masked modelling, MAE-style imputation).
On a small target, the encoder is applied to each cell's observed proteins (target-z-scored) and the decoder's output for
the unobserved proteins is mapped back with the target's per-protein mean/SD. Nothing is fitted on the target except
those per-protein moments, so the model cannot overfit 30 cells.
"""

from __future__ import annotations

import numpy as np


def _torch():
    import torch

    torch.set_num_threads(max(1, torch.get_num_threads()))
    return torch


def train_mae(Zref: np.ndarray, hidden: int = 512, latent: int = 64, epochs: int = 30, batch: int = 256,
              mask_frac: float = 0.3, lr: float = 1e-3, weight_decay: float = 1e-5, dropout: float = 0.1,
              min_obs: int = 20, seed: int = 0, verbose: bool = False):
    """Pretrain on the reference. Returns a dict with the model state and the warm-protein mask."""
    torch = _torch()
    import torch.nn as nn

    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    obs = ~np.isnan(Zref)
    warm = obs.sum(0) >= min_obs
    Z = np.where(obs, Zref, 0.0)[:, warm].astype(np.float32)
    O = obs[:, warm].astype(np.float32)
    keep = O.sum(1) >= 20
    Z, O = Z[keep], O[keep]
    p = Z.shape[1]
    model = nn.Sequential(
        nn.Linear(2 * p, hidden), nn.GELU(), nn.Dropout(dropout),
        nn.Linear(hidden, latent), nn.GELU(),
        nn.Linear(latent, hidden), nn.GELU(),
        nn.Linear(hidden, p))
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    Zt, Ot = torch.from_numpy(Z), torch.from_numpy(O)
    n = Z.shape[0]
    for ep in range(epochs):
        perm = rng.permutation(n)
        tot = 0.0
        model.train()
        for s in range(0, n, batch):
            idx = torch.from_numpy(perm[s:s + batch])
            z, o = Zt[idx], Ot[idx]
            hide = (torch.rand_like(o) < mask_frac) & (o > 0)
            vis = o * (~hide)
            pred = model(torch.cat([z * vis, vis], 1))
            loss = (((pred - z) ** 2) * hide).sum() / hide.sum().clamp(min=1)
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += float(loss) * len(idx)
        if verbose:
            print(f"epoch {ep + 1}/{epochs} masked-MSE {tot / n:.4f}", flush=True)
    model.eval()
    return {"model": model, "warm": warm}


def mae_impute(Xm: np.ndarray, fitted: dict, base: np.ndarray) -> np.ndarray:
    """Predict unobserved entries of Xm (cells x proteins, NaN = missing) for warm proteins; others keep `base`."""
    torch = _torch()
    model, warm = fitted["model"], fitted["warm"]
    obs = ~np.isnan(Xm)
    mu = np.nanmean(Xm, 0)
    sd = np.nanstd(Xm, 0)
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    ok = np.isfinite(mu)
    Z = np.where(obs, (Xm - np.nan_to_num(mu)) / sd, 0.0)[:, warm].astype(np.float32)
    O = (obs & ok[None, :])[:, warm].astype(np.float32)
    with torch.no_grad():
        P = model(torch.from_numpy(np.concatenate([Z * O, O], 1))).numpy()
    out = base.copy()
    wi = np.where(warm)[0]
    pred = mu[wi][None, :] + sd[wi][None, :] * P
    miss = ~obs[:, wi] & ok[wi][None, :]
    sub = out[:, wi]
    sub[miss] = pred[miss]
    out[:, wi] = sub
    return out
