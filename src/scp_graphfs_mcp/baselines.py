"""Comparator panel for the v9 validation: imputation methods researchers typically reach for, plus graph-prior methods
that are structurally distinct from the ones tried in Phases A-C.

Every method has the signature `fn(Xm, ctx) -> Y` (cells x proteins, log2, NaN = missing; Y is a full matrix and only
its unobserved entries are scored). ctx carries `names` (protein ids), `base` (the kNN fill used everywhere else) and
`seed`; graphs are built lazily from `names` and cached in ctx. Nothing here sees a reference dataset: these are the
methods a user would run on their own matrix (the graph methods add only a public knowledge graph).

Data-only:
  mean, mindet, minprob        per-protein mean; left-censored MinDet / MinProb (Lazar et al. 2016, imputeLCMD)
  knn                          the benchmark's kNN (phase2.loo_knn), for reference
  softimpute                   nuclear-norm soft-thresholded SVD (Mazumder, Hastie & Tibshirani 2010)
  itersvd                      rank-k iterative hard SVD (Troyanskaya et al. 2001)
  bpca (ppca kept, unused)     probabilistic PCA by EM with missing data; bpca adds ARD on the loadings (Bishop 1999),
                               the model behind pcaMethods' BPCA (Oba et al. 2003)
  mice                         chained equations, BayesianRidge, 50 nearest features (sklearn IterativeImputer)
  missforest                   chained equations with random forests (missForest, Stekhoven & Buehlmann 2012)
  dae                          masked/denoising autoencoder trained on the target alone (PIMMS-style, transductive)
Graph priors (STRING v12.5 combined >= 700 unless noted):
  grmf_string                  graph-Laplacian-regularised low-rank MF (the Phase-A E3 model, STRING graph)
  grmf_datagraph               the same with a data-derived co-variation graph of equal edge count (control)
  rwr_resid                    multi-hop network propagation (random walk with restart) of kNN residuals
  gcn                          2-layer graph convolutional network over proteins, trained by masked reconstruction
  imc_embed                    inductive MF: protein factors tied to spectral embeddings of the STRING graph
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from . import phase2 as P
from .methods import adjacency


# ------------------------------------------------------------------------------------------------- helpers
def _stats(Xm):
    obs = ~np.isnan(Xm)
    mu = np.nanmean(Xm, 0)
    sd = np.nanstd(Xm, 0)
    ok = np.isfinite(mu)
    mu = np.where(ok, mu, np.nanmean(Xm))
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    return obs, mu, sd, ok


def _keep_obs(Xm, Y, ctx):
    Y = np.where(np.isfinite(Y), Y, ctx["base"])
    return np.where(np.isnan(Xm), Y, Xm)


def _graph(ctx, kind="string"):
    key = f"_graph_{kind}"
    if key not in ctx:
        names = ctx["names"]
        if kind == "string":
            from .priors import build_prior

            ctx[key] = build_prior("string", names)
        elif kind == "data":
            n_edges = max(1, _graph(ctx, "string").number_of_edges())
            ctx[key] = P.data_graph(ctx["_Xm_for_graph"], names, n_edges)
        else:
            raise ValueError(kind)
    return ctx[key]


def _norm_adj(ctx, self_loops=True):
    key = f"_nadj_{self_loops}"
    if key not in ctx:
        A = adjacency(_graph(ctx), ctx["names"]).astype(float)
        if self_loops:
            A = A + sp.eye(A.shape[0])
        d = np.asarray(A.sum(1)).ravel()
        dinv = np.where(d > 0, 1.0 / np.sqrt(np.maximum(d, 1e-12)), 0.0)
        ctx[key] = (sp.diags(dinv) @ A @ sp.diags(dinv)).tocsr()
    return ctx[key]


# ------------------------------------------------------------------------------------------------- data-only
def mean(Xm, ctx):
    _, mu, _, _ = _stats(Xm)
    return _keep_obs(Xm, np.broadcast_to(mu, Xm.shape).copy(), ctx)


def mindet(Xm, ctx, q=0.01):
    """MinDet: every missing value in a cell = the q-quantile of that cell's observed values."""
    lo = np.nanquantile(Xm, q, axis=1)
    return _keep_obs(Xm, np.broadcast_to(lo[:, None], Xm.shape).copy(), ctx)


def minprob(Xm, ctx, q=0.01, tune_sigma=1.0):
    """MinProb: draws from N(q-quantile of the cell, (tune_sigma * median protein SD)^2)."""
    lo = np.nanquantile(Xm, q, axis=1)
    sdp = np.nanstd(Xm, 0)
    sigma = tune_sigma * np.nanmedian(sdp[np.isfinite(sdp) & (sdp > 0)])
    rng = np.random.default_rng(ctx["seed"] + 101)
    return _keep_obs(Xm, lo[:, None] + sigma * rng.normal(size=Xm.shape), ctx)


def knn(Xm, ctx):
    return ctx["base"]


def softimpute(Xm, ctx, shrink=None, iters=200, tol=1e-5):
    """Soft-Impute on per-protein standardised data; default shrinkage = max singular value / 50 (fancyimpute)."""
    obs, mu, sd, _ = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, 0.0)
    Zf = Z.copy()
    lam = shrink if shrink is not None else np.linalg.svd(Z, compute_uv=False)[0] / 50.0
    prev = np.zeros_like(Z)
    for _ in range(iters):
        U, s, Vt = np.linalg.svd(Zf, full_matrices=False)
        s = np.maximum(s - lam, 0.0)
        L = (U * s) @ Vt
        Zf = np.where(obs, Z, L)
        if np.linalg.norm(L - prev) <= tol * max(np.linalg.norm(prev), 1e-12):
            break
        prev = L
    return _keep_obs(Xm, mu + sd * L, ctx)


def itersvd(Xm, ctx, rank=None, iters=100, tol=1e-5):
    """rank defaults to 20% of cells (Troyanskaya et al.'s eigengene heuristic)."""
    obs, mu, sd, _ = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, 0.0)
    Zf = Z.copy()
    k = min(rank or max(2, Xm.shape[0] // 5), min(Z.shape) - 1)
    prev = np.zeros_like(Z)
    for _ in range(iters):
        U, s, Vt = np.linalg.svd(Zf, full_matrices=False)
        L = (U[:, :k] * s[:k]) @ Vt[:k]
        Zf = np.where(obs, Z, L)
        if np.linalg.norm(L - prev) <= tol * max(np.linalg.norm(prev), 1e-12):
            break
        prev = L
    return _keep_obs(Xm, mu + sd * L, ctx)


def _ppca(Xm, ctx, rank=10, iters=200, ard=False, tol=1e-6):
    """EM for probabilistic PCA with missing entries (Ilin & Raiko 2010, unit-scaled data); ard=True puts a
    per-component precision alpha_k on W (Bayesian PCA), re-estimated each iteration (Bishop 1999)."""
    obs, mu, sd, ok = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, 0.0)
    n, p = Z.shape
    k = min(rank, n - 1)
    rng = np.random.default_rng(ctx["seed"] + 7)
    W = 0.1 * rng.normal(size=(p, k))
    m = np.zeros(p)
    s2 = 1.0
    alpha = np.full(k, 1e-3)
    O = obs.astype(float)
    last = np.inf
    for _ in range(iters):
        # E-step: per-cell posterior over latent x_i using that cell's observed proteins
        Ex = np.zeros((n, k))
        Sx = np.zeros((n, k, k))
        for i in range(n):
            o = obs[i]
            Wo = W[o]
            M = Wo.T @ Wo + s2 * np.eye(k)
            Minv = np.linalg.inv(M)
            Ex[i] = Minv @ (Wo.T @ (Z[i, o] - m[o]))
            Sx[i] = s2 * Minv
        # M-step: per-protein W_j, m_j using cells that observe protein j
        R = np.eye(k) * 0.0
        err = 0.0
        for j in range(p):
            c = obs[:, j]
            if c.sum() <= k:  # too few cells to identify W_j: mean only (avoids singular solves)
                W[j] = 0.0
                m[j] = Z[c, j].mean() if c.any() else 0.0
                continue
            Xc = Ex[c]
            A = Xc.T @ Xc + Sx[c].sum(0)
            if ard:
                A = A + s2 * np.diag(alpha)
            b = Xc.T @ (Z[c, j] - m[j])
            W[j] = np.linalg.solve(A + 1e-3 * s2 * np.eye(k), b)
            m[j] = np.mean(Z[c, j] - Xc @ W[j])
        resid = np.where(obs, Z - m - Ex @ W.T, 0.0)
        tr = sum(np.einsum("ij,jk,ik->", W[obs[i]], Sx[i], W[obs[i]]) for i in range(n))
        s2 = max((np.sum(resid ** 2) + tr) / O.sum(), 1e-6)
        if ard:
            alpha = p / np.maximum((W ** 2).sum(0), 1e-12)
        err = float(np.sum(resid ** 2))
        if abs(last - err) <= tol * max(err, 1e-12):
            break
        last = err
    L = m + Ex @ W.T
    return _keep_obs(Xm, mu + sd * L, ctx)


def ppca(Xm, ctx):
    return _ppca(Xm, ctx, ard=False)


def bpca(Xm, ctx):
    return _ppca(Xm, ctx, ard=True)


def _iterative(Xm, ctx, estimator, max_iter, n_near=50):
    from sklearn.experimental import enable_iterative_imputer  # noqa: F401
    from sklearn.impute import IterativeImputer

    obs, mu, sd, ok = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, np.nan)[:, ok]
    imp = IterativeImputer(estimator=estimator, max_iter=max_iter, n_nearest_features=min(n_near, Z.shape[1] - 1),
                           initial_strategy="mean", random_state=ctx["seed"], skip_complete=True,
                           keep_empty_features=True)
    Zi = imp.fit_transform(Z)
    L = np.zeros(Xm.shape)
    L[:, ok] = Zi
    return _keep_obs(Xm, mu + sd * L, ctx)


def mice(Xm, ctx):
    from sklearn.linear_model import BayesianRidge

    return _iterative(Xm, ctx, BayesianRidge(), max_iter=5)


def missforest(Xm, ctx):
    from sklearn.ensemble import RandomForestRegressor

    # reduced from missForest's defaults (100 trees, <= 10 iterations) for runtime: ~6 min -> ~1 min per call
    rf = RandomForestRegressor(n_estimators=10, max_features=0.33, min_samples_leaf=1, random_state=ctx["seed"],
                               n_jobs=1)
    return _iterative(Xm, ctx, rf, max_iter=3)


def dae(Xm, ctx, epochs=300):
    from . import mae as MAE

    obs, mu, sd, ok = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, np.nan)
    fitted = MAE.train_mae(Z, hidden=256, latent=16, epochs=epochs, batch=32, min_obs=5, seed=ctx["seed"],
                           dropout=0.2, weight_decay=1e-4)
    return _keep_obs(Xm, MAE.mae_impute(Xm, fitted, ctx["base"]), ctx)


# ------------------------------------------------------------------------------------------------- graph priors
def grmf_string(Xm, ctx):
    return _keep_obs(Xm, P.grmf(Xm, _graph(ctx), ctx["names"], lam_g=1.0), ctx)


def grmf_datagraph(Xm, ctx):
    ctx = dict(ctx)
    ctx["_Xm_for_graph"] = Xm
    ctx.pop("_graph_data", None)
    return _keep_obs(Xm, P.grmf(Xm, _graph(ctx, "data"), ctx["names"], lam_g=1.0), ctx)


def rwr_resid(Xm, ctx, alpha=0.5, steps=20):
    """Random walk with restart over STRING applied to the kNN residuals of observed entries (multi-hop E1):
    r_hat = sum_{t>=1} alpha^t (R W^t) / sum_{t>=1} alpha^t (O W^t), self excluded, then Y = base + sd * r_hat."""
    obs, mu, sd, ok = _stats(Xm)
    base = ctx["base"]
    loo = P.loo_knn(Xm)  # leave-one-out kNN estimate of every entry (the base equals Xm on observed entries)
    Rz = np.where(obs & np.isfinite(loo), (Xm - loo) / sd, 0.0)
    O = (obs & np.isfinite(loo)).astype(float)
    W = _norm_adj(ctx, self_loops=False).T
    num = np.zeros_like(Rz)
    den = np.zeros_like(Rz)
    a, b = Rz, O
    for t in range(1, steps + 1):
        a = (W @ a.T).T
        b = (W @ b.T).T
        num += alpha ** t * a
        den += alpha ** t * b
    rhat = np.where(den > 1e-3, num / np.maximum(den, 1e-12), 0.0)
    shrink = den / (den + 0.5)  # few/weak network partners -> stay near kNN
    return _keep_obs(Xm, base + sd * shrink * rhat, ctx)


def gcn(Xm, ctx, hidden=64, epochs=400, lr=1e-2, hide=0.15):
    """Proteins are nodes; node features are the protein's standardised values across cells plus the observed mask.
    out = A relu(A [Z^T|M^T] W1) W2 (A = sym-normalised STRING with self loops). Trained to reconstruct randomly
    hidden observed entries (hidden from the input, scored in the loss)."""
    import torch

    torch.manual_seed(ctx["seed"])
    obs, mu, sd, ok = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, 0.0).T.astype(np.float32)  # p x n
    O = obs.T.astype(np.float32)
    A = _norm_adj(ctx, self_loops=True).tocoo()
    At = torch.sparse_coo_tensor(np.vstack([A.row, A.col]), A.data.astype(np.float32), A.shape).coalesce()
    n = Z.shape[1]
    W1 = torch.nn.Linear(2 * n, hidden)
    W2 = torch.nn.Linear(hidden, n)
    opt = torch.optim.Adam(list(W1.parameters()) + list(W2.parameters()), lr=lr, weight_decay=1e-4)
    Zt, Ot = torch.from_numpy(Z), torch.from_numpy(O)

    def fwd(z, o):
        h = torch.relu(torch.sparse.mm(At, W1(torch.cat([z * o, o], 1))))
        return torch.sparse.mm(At, W2(torch.nn.functional.dropout(h, 0.2, training=W1.training)))

    for _ in range(epochs):
        W1.train()
        hid = (torch.rand_like(Ot) < hide) & (Ot > 0)
        vis = Ot * (~hid)
        out = fwd(Zt, vis)
        loss = (((out - Zt) ** 2) * hid).sum() / hid.sum().clamp(min=1)
        opt.zero_grad()
        loss.backward()
        opt.step()
    W1.eval()
    with torch.no_grad():
        out = fwd(Zt, Ot).numpy().T
    return _keep_obs(Xm, mu + sd * out, ctx)


def imc_embed(Xm, ctx, rank=10, dim=32, lam=5.0, lam_d=10.0, iters=15):
    """Z ~ U V^T with V = E B + D: E = top-`dim` eigenvectors of the normalised STRING adjacency (proteins without
    edges get E = 0), B learned by ridge, D a ridge-penalised protein-specific deviation (lam_d)."""
    from scipy.sparse.linalg import eigsh

    obs, mu, sd, ok = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, 0.0)
    n, p = Z.shape
    if "_spec_E" not in ctx:
        A = _norm_adj(ctx, self_loops=False)
        k = min(dim, p - 2)
        try:
            vals, vecs = eigsh(A, k=k, which="LA")
        except Exception:
            vecs = np.zeros((p, k))
        deg = np.asarray(adjacency(_graph(ctx), ctx["names"]).sum(1)).ravel()
        vecs[deg == 0] = 0.0
        ctx["_spec_E"] = vecs * np.sqrt(p)
    E = ctx["_spec_E"]
    rng = np.random.default_rng(ctx["seed"])
    V = 0.1 * rng.normal(size=(p, rank))
    B = np.zeros((E.shape[1], rank))
    I = np.eye(rank)
    U = np.zeros((n, rank))
    for _ in range(iters):
        for i in range(n):
            o = obs[i] & ok
            U[i] = np.linalg.solve(V[o].T @ V[o] + lam * I, V[o].T @ Z[i, o])
        prior = E @ B
        for j in np.where(ok)[0]:
            o = obs[:, j]
            Uo = U[o]
            V[j] = np.linalg.solve(Uo.T @ Uo + lam_d * I, Uo.T @ Z[o, j] + lam_d * prior[j])
        B = np.linalg.solve(E.T @ E + 1.0 * np.eye(E.shape[1]), E.T @ V)
    return _keep_obs(Xm, mu + sd * (U @ V.T), ctx)


# ppca (ML, no ARD) is defined but not in the panel: at 30 cells it diverges (RMSE 8-13x the champion on dev split
# petrosius2023/MCAR/100 even with a variance floor); bpca is the variant used in practice at small n.
DATA_ONLY = {"mean": mean, "mindet": mindet, "minprob": minprob, "knn": knn, "softimpute": softimpute,
             "itersvd": itersvd, "bpca": bpca, "mice": mice, "missforest": missforest, "dae": dae}
GRAPH = {"grmf_string": grmf_string, "grmf_datagraph": grmf_datagraph, "rwr_resid": rwr_resid, "gcn": gcn,
         "imc_embed": imc_embed}
ALL = {**DATA_ONLY, **GRAPH}
