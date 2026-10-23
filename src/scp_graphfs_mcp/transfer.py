"""Transfer strategies for tiny targets (EXPLORATION.md Amendment 10).

Cross-domain analogies:
- **Reference projection** (ancestry PCA / recommender cold-start / a linear masked auto-encoder pretrained elsewhere):
  protein loadings V are learned on a large *external* reference (other single-cell datasets, optionally bulk
  proteomics panels); the few target cells only fit their own scores u_i. Proteins the reference never saw get V from
  their knowledge-graph neighbours (the curated prior enters exactly where data cannot).
- **Reference donors** (genotype imputation copies from reference haplotypes): external cells join the kNN donor pool.

All references exclude the target (leave-one-dataset-out) and never contain holdouts.
"""

from __future__ import annotations

import numpy as np

from .external import CORPUS, HOLDOUTS, _load
from .ids import group_accessions


def _target_accs(var_names):
    return [(group_accessions(f) or [""])[0] for f in var_names]


def _z(X):
    mu, sd = np.nanmean(X, 0), np.nanstd(X, 0)
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    return (X - mu) / sd


# Amendment 11: extra scpdata datasets built by scripts/corpus_ext/ (whatever built successfully; see PROVENANCE).
EXCLUDED_FROM_CORPUS = HOLDOUTS | {"leduc2022_plexDIA", "specht2019v2"}
NON_HUMAN = {"woo2022"}  # woo2022_macrophage is mouse: 0 accessions overlap human datasets (found Oct 23); never contributed
AMENDMENT13_HOLDOUT = "brunner2022"  # reference-only until Amendment 13; never in a reference from then on
PHASE_B_HOLDOUTS = {"brunner2022", "wu2026", "derks2022"}  # Amendment 15: never in any Phase-B reference


def extended_corpus() -> list[str]:
    from .external import CORPUS_DIR

    ext = sorted(f.stem for f in CORPUS_DIR.glob("*.npz") if f.stem not in CORPUS)
    ext = [e for e in ext if not e.endswith("_labels")]
    bad = set(ext) & EXCLUDED_FROM_CORPUS
    if bad:
        raise ValueError(f"excluded datasets present in the corpus dir: {bad}")
    return [d for d in list(CORPUS) + ext if d not in NON_HUMAN]


def corpus_reference(var_names, exclude: tuple[str, ...], datasets: list[str] | None = None) -> np.ndarray:
    """External single cells aligned to var_names: (N_cells x p), z-scored per protein *within* each dataset.
    datasets=None -> the Amendment-4 corpus; pass `extended_corpus()` for Amendment 11."""
    if set(exclude) & HOLDOUTS or set(datasets or []) & EXCLUDED_FROM_CORPUS:
        raise ValueError("holdouts are never part of the corpus")
    target = _target_accs(var_names)
    blocks = []
    for name in (datasets or CORPUS):
        if name in exclude:
            continue
        X, accs = _load(name)
        pos = {a: i for i, a in enumerate(accs)}
        B = np.full((X.shape[0], len(target)), np.nan)
        cols = [(j, pos[a]) for j, a in enumerate(target) if a in pos]
        if cols:
            tj, sj = zip(*cols)
            B[:, list(tj)] = X[:, list(sj)]
        blocks.append(_z(B))
    return np.vstack(blocks)


BULK_MODULES = ("ccle", "proteomehd")  # Phase-A default; Phase B passes more (procan, gtex_proteome, cptac, ...)


def bulk_reference(var_names, modules=BULK_MODULES) -> np.ndarray:
    """Bulk reference samples (one block per module in `reps`) aligned to var_names, (N_samples x p), z-scored per
    protein (NaN = absent). Modules that are not built yet are skipped."""
    import importlib

    out = []
    for name in modules:
        try:
            m = importlib.import_module(f"scp_graphfs_mcp.reps.{name}")
        except ModuleNotFoundError:
            continue
        R, has = m.reference(var_names)
        R = np.where(np.asarray(has, bool)[:, None], R, np.nan)
        out.append(_z(R.T))
    return np.vstack(out) if out else np.zeros((0, len(var_names)))


def fit_loadings(Zref: np.ndarray, rank: int = 20, lam: float = 5.0, iters: int = 15, min_obs: int = 20,
                 seed: int = 0, kg=None, kappa: float = 0.0):
    """Masked low-rank factorisation Zref ~ U V^T by alternating ridge. Returns (V, warm) with V (p x rank);
    proteins observed in < min_obs reference samples are 'cold' (warm=False, V row = 0).

    kg / kappa (Amendment 13): a biological prior on the loadings. kg is a per-protein list of (neighbour idx,
    weights) over the full protein axis (e.g. curated complex partners). Each warm protein's loading is shrunk
    toward the weight-normalised mean loading of its warm neighbours with strength kappa * n_obs_j, i.e. a
    graph-Laplacian penalty whose weight scales with the evidence for that protein (kappa = 0: plain ALS)."""
    obs = ~np.isnan(Zref)
    warm = obs.sum(0) >= min_obs
    Z = np.where(obs, Zref, 0.0)[:, warm]
    O = obs[:, warm]
    keep = O.sum(1) >= rank
    Z, O = Z[keep], O[keep]
    rng = np.random.default_rng(seed)
    V = 0.1 * rng.normal(size=(Z.shape[1], rank))
    U = np.zeros((Z.shape[0], rank))
    I = lam * np.eye(rank)
    nb = None
    if kg is not None and kappa > 0:
        col = -np.ones(Zref.shape[1], int)
        col[np.where(warm)[0]] = np.arange(warm.sum())
        nb = []
        for j in np.where(warm)[0]:
            idx, w = kg[j]
            ok = (col[idx] >= 0) if len(idx) else np.zeros(0, bool)
            ww = np.abs(np.asarray(w, float)[ok]) if len(idx) else np.zeros(0)
            nb.append((col[idx[ok]], ww / ww.sum()) if ww.sum() > 0 else (np.zeros(0, int), np.zeros(0)))
    nobs = O.sum(0)
    for _ in range(iters):
        for i in range(Z.shape[0]):
            o = O[i]
            Vo = V[o]
            U[i] = np.linalg.solve(Vo.T @ Vo + I, Vo.T @ Z[i, o])
        for j in range(Z.shape[1]):
            o = O[:, j]
            Uo = U[o]
            A, b = Uo.T @ Uo + I, Uo.T @ Z[o, j]
            if nb is not None and len(nb[j][0]):
                g = kappa * nobs[j]
                A = A + g * np.eye(rank)
                b = b + g * (nb[j][1] @ V[nb[j][0]])
            V[j] = np.linalg.solve(A, b)
    Vfull = np.zeros((Zref.shape[1], rank))
    Vfull[warm] = V
    return Vfull, warm


def coldstart(V: np.ndarray, warm: np.ndarray, nbrs) -> tuple[np.ndarray, np.ndarray]:
    """Give each cold protein the weighted mean loading of its *warm* knowledge neighbours (nbrs: list of
    (idx, weights) as in refpanel). Returns (V', usable mask)."""
    V2 = V.copy()
    usable = warm.copy()
    for j in np.where(~warm)[0]:
        idx, w = nbrs[j]
        if not len(idx):
            continue
        ok = warm[idx]
        if not ok.any():
            continue
        ww = np.abs(w[ok])
        V2[j] = (ww[:, None] * V[idx[ok]]).sum(0) / max(ww.sum(), 1e-9)
        usable[j] = True
    return V2, usable


def project(Xm: np.ndarray, V: np.ndarray, usable: np.ndarray, base: np.ndarray, lam: float = 5.0) -> np.ndarray:
    """Fit each target cell's scores on its observed, usable proteins; predict unobserved usable entries.
    Entries of unusable proteins keep `base` (the kNN fill)."""
    obs = ~np.isnan(Xm)
    mu = np.nanmean(Xm, 0)
    sd = np.nanstd(Xm, 0)
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    Z = np.where(obs, (Xm - np.nan_to_num(mu)) / sd, 0.0)
    out = base.copy()
    I = lam * np.eye(V.shape[1])
    ok_col = usable & np.isfinite(mu)
    for i in range(Xm.shape[0]):
        o = obs[i] & ok_col
        if o.sum() < 3:
            continue
        u = np.linalg.solve(V[o].T @ V[o] + I, V[o].T @ Z[i, o])
        miss = ~obs[i] & ok_col
        out[i, miss] = mu[miss] + sd[miss] * (V[miss] @ u)
    return out


def permute_rows(V: np.ndarray, warm: np.ndarray, seed: int):
    """Null: shuffle which protein gets which loading (among warm proteins). Keeps the reference's spectrum,
    breaks protein identity."""
    V2 = V.copy()
    idx = np.where(warm)[0]
    V2[idx] = V[np.random.default_rng(seed).permutation(idx)]
    return V2


def donor_knn(Xm: np.ndarray, Zref: np.ndarray, base: np.ndarray, k: int = 10, n_donors: int = 2000,
              seed: int = 0) -> np.ndarray:
    """Reference-donor kNN: donors = external cells (z-space, mapped to target scale with target mu/sd).
    Distance = mean squared z-difference over co-observed proteins (>= 20)."""
    obs = ~np.isnan(Xm)
    mu = np.nanmean(Xm, 0)
    sd = np.nanstd(Xm, 0)
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    Z = np.where(obs, (Xm - np.nan_to_num(mu)) / sd, np.nan)
    rng = np.random.default_rng(seed)
    R = Zref[rng.choice(len(Zref), min(n_donors, len(Zref)), replace=False)]
    Ro = ~np.isnan(R)
    R0 = np.where(Ro, R, 0.0)
    out = base.copy()
    for i in range(Xm.shape[0]):
        o = obs[i]
        zi = np.where(o, Z[i], 0.0)
        cnt = Ro[:, o].sum(1)
        d = ((R0[:, o] - zi[o]) ** 2 * Ro[:, o]).sum(1) / np.maximum(cnt, 1)
        d = np.where(cnt >= 20, d, np.inf)
        for j in np.where(~o & np.isfinite(mu))[0]:
            has = Ro[:, j] & np.isfinite(d)
            if has.sum() < k:
                continue
            nn = np.argsort(np.where(has, d, np.inf))[:k]
            out[i, j] = mu[j] + sd[j] * R[nn, j].mean()
    return out


def reference_similarity(Zref: np.ndarray, Xm: np.ndarray, min_overlap: int = 20, top: int = 3) -> np.ndarray:
    """Score each reference sample by how much it resembles the target cells (Amendment 15, reference selection).

    Both sides are z-scored per protein within their own dataset; the score of reference sample r is the mean of its
    `top` highest Pearson correlations (over co-observed proteins) with the target cells. Uses only observed target
    entries."""
    Zt = _z(Xm)
    ot, orf = ~np.isnan(Zt), ~np.isnan(Zref)
    A, B = np.where(orf, Zref, 0.0), np.where(ot, Zt, 0.0)
    Of, Ot = orf.astype(float), ot.astype(float)
    n = Of @ Ot.T
    sxy = A @ B.T
    sxx = (A ** 2) @ Ot.T
    syy = Of @ (B ** 2).T
    r = sxy / np.sqrt(np.maximum(sxx * syy, 1e-12))
    r[n < min_overlap] = -1.0
    k = min(top, r.shape[1])
    return np.sort(r, axis=1)[:, -k:].mean(1)


def select_reference(Zref: np.ndarray, Xm: np.ndarray, n: int, src: np.ndarray | None = None, mode: str = "cells",
                     n_datasets: int = 3, seed: int | None = None) -> np.ndarray:
    """Row indices of the selected reference. mode='cells': top-n samples by similarity; mode='datasets': all samples
    of the n_datasets sources with the highest mean similarity (needs src); seed given -> random n rows (control)."""
    if seed is not None:
        return np.sort(np.random.default_rng(seed).choice(len(Zref), min(n, len(Zref)), replace=False))
    s = reference_similarity(Zref, Xm)
    if mode == "datasets":
        srcs = np.unique(src)
        means = {d: s[src == d].mean() for d in srcs}
        best = sorted(srcs, key=lambda d: -means[d])[:n_datasets]
        return np.where(np.isin(src, best))[0]
    return np.sort(np.argsort(-s)[:min(n, len(Zref))])


def adapt_project(Xm: np.ndarray, Vref: np.ndarray, usable: np.ndarray, base: np.ndarray, lam_v: float = 10.0,
                  lam: float = 5.0, extra: int = 3, iters: int = 10, seed: int = 0) -> np.ndarray:
    """Domain-adapted reference projection (Amendment 15): MAP low-rank on the *target*, with the reference loadings as
    a Gaussian prior, min |P_obs(Z - U V^T)|^2 + lam|U|^2 + lam_v|V - V_ref|^2, plus `extra` free factors (ridge only).
    Unobserved entries are predicted from the adapted factors; proteins with no target observations keep `base`."""
    obs = ~np.isnan(Xm)
    mu = np.nanmean(Xm, 0)
    sd = np.nanstd(Xm, 0)
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    ok = np.isfinite(mu)
    Z = np.where(obs, (Xm - np.nan_to_num(mu)) / sd, 0.0)
    r0 = Vref.shape[1]
    rng = np.random.default_rng(seed)
    Vp = np.hstack([np.where(usable[:, None], Vref, 0.0), np.zeros((len(Vref), extra))])
    V = Vp + np.hstack([np.zeros_like(Vref), 0.01 * rng.normal(size=(len(Vref), extra))])
    r = r0 + extra
    U = np.zeros((Xm.shape[0], r))
    I = np.eye(r)
    for _ in range(iters):
        for i in range(Xm.shape[0]):
            o = obs[i] & ok
            U[i] = np.linalg.solve(V[o].T @ V[o] + lam * I, V[o].T @ Z[i, o])
        for j in np.where(ok)[0]:
            o = obs[:, j]
            Uo = U[o]
            V[j] = np.linalg.solve(Uo.T @ Uo + lam_v * I, Uo.T @ Z[o, j] + lam_v * Vp[j])
    out = base.copy()
    P = U @ V.T
    miss = ~obs & ok[None, :]
    out[miss] = (mu[None, :] + sd[None, :] * P)[miss]
    return out
