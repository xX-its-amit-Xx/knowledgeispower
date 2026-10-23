"""Phase-C idea `ebshrink` (queen: representation & math).

Empirical-Bayes per-protein adaptation strength lambda_v,j for the champion's adapted projection.
  A. Fit rank-20 loadings separately on each (subsampled) leave-target-out reference dataset, Procrustes-align
     (with scale) to the pooled V, and estimate each protein's between-dataset loading heterogeneity tau2_j by
     method of moments (raw spread minus ridge sampling noise), shrunk toward the median.
  B. lambda_j = c / tau2_j (median over warm proteins = 100, clipped to [10, 1000]); used as a diagonal per-protein
     penalty in a verbatim copy of transfer.adapt_project.
Null: same pipeline, lambda vector permuted among warm proteins (rng 1000 + seed).
No external data. Per-target cache: data/cache/phasec/ebshrink/<target>_tau.npz.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "ebshrink"
RANK, LAM, MIN_OBS, ITERS_D = 20, 5.0, 20, 8
MAX_CELLS, MIN_CELLS, MIN_S = 300, 40, 50
LAM_MED, LAM_LO, LAM_HI = 100.0, 10.0, 1000.0


def _fit_loadings_U(Zref, rank=RANK, lam=LAM, iters=ITERS_D, min_obs=MIN_OBS, seed=0):
    """Local copy of transfer.fit_loadings (no kg) that also returns U and the kept rows' Z / O."""
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
    for _ in range(iters):
        for i in range(Z.shape[0]):
            o = O[i]
            Vo = V[o]
            U[i] = np.linalg.solve(Vo.T @ Vo + I, Vo.T @ Z[i, o])
        for j in range(Z.shape[1]):
            o = O[:, j]
            Uo = U[o]
            V[j] = np.linalg.solve(Uo.T @ Uo + I, Uo.T @ Z[o, j])
    Vfull = np.zeros((Zref.shape[1], rank))
    Vfull[warm] = V
    return Vfull, warm, U, Z, O


def _noise(V, warm, U, Z, O, lam=LAM):
    """noise_j = s2_j * trace(inv(U_o'U_o + lam I)) / rank for warm proteins (full-length vector, nan elsewhere)."""
    p = len(warm)
    out = np.full(p, np.nan)
    Vw = V[warm]
    I = lam * np.eye(U.shape[1])
    for k, j in enumerate(np.where(warm)[0]):
        o = O[:, k]
        if o.sum() < 2:
            continue
        Uo = U[o]
        res = Z[o, k] - Uo @ Vw[k]
        s2 = float(np.mean(res ** 2))
        out[j] = s2 * np.trace(np.linalg.inv(Uo.T @ Uo + I)) / U.shape[1]
    return out


def compute_tau(Z, src, Vpool, warm, log=print):
    p = Z.shape[1]
    rng = np.random.default_rng(0)
    raw_sum, nz_sum, n = np.zeros(p), np.zeros(p), np.zeros(p, int)
    used = []
    for d in np.unique(src):
        rows = np.where(src == d)[0]
        if len(rows) < MIN_CELLS:
            continue
        if len(rows) > MAX_CELLS:
            rows = np.sort(rng.choice(rows, MAX_CELLS, replace=False))
        Vd, warm_d, U, Zd, Od = _fit_loadings_U(Z[rows])
        if U.shape[0] < MIN_CELLS:
            continue
        S = warm & warm_d
        if S.sum() < MIN_S:
            continue
        A, sig, Bt = np.linalg.svd(Vd[S].T @ Vpool[S])
        R = A @ Bt
        s = sig.sum() / np.sum(Vd[S] ** 2)
        Vt = s * Vd @ R
        nz = _noise(Vd, warm_d, U, Zd, Od) * s ** 2
        ok = S & np.isfinite(nz)
        raw_sum[ok] += np.sum((Vt[ok] - Vpool[ok]) ** 2, 1) / RANK
        nz_sum[ok] += nz[ok]
        n[ok] += 1
        used.append((str(d), int(len(rows)), int(S.sum()), float(s)))
        log(f"  ebshrink tau: {d} cells={len(rows)} |S|={S.sum()} scale={s:.3f}")
    has = n > 0
    raw = np.full(p, np.nan)
    nzm = np.full(p, np.nan)
    raw[has] = raw_sum[has] / n[has]
    nzm[has] = nz_sum[has] / n[has]
    med_raw = np.nanmedian(raw[has])
    tau = np.full(p, np.nan)
    tau[has] = np.maximum(raw[has] - nzm[has], 0.1 * med_raw)
    good = warm & (n >= 2)
    med = float(np.median(tau[good]))
    tau[good] = (n[good] * tau[good] + 3 * med) / (n[good] + 3)
    tau[~good] = med
    lam = np.full(p, LAM_MED)
    c = LAM_MED / np.median(1.0 / tau[warm])
    lam[warm] = np.clip(c / tau[warm], LAM_LO, LAM_HI)
    return {"lam": lam, "tau2": tau, "n": n, "raw": raw, "noise": nzm, "warm": warm,
            "used": np.array([u[0] for u in used]), "used_cells": np.array([u[1] for u in used]),
            "used_S": np.array([u[2] for u in used]), "used_scale": np.array([u[3] for u in used])}


def get_lambda(ctx):
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{ctx['target']}_tau.npz"
    if f.exists():
        return np.load(f)["lam"]
    Z, src = ctx["ref"](), ctx["ref_src"]()
    r = compute_tau(Z, src, ctx["V"], ctx["warm"])
    r["n_obs_ref"] = (~np.isnan(Z)).sum(0)
    tmp = f.with_suffix(f".{os.getpid()}.tmp.npz")
    np.savez_compressed(tmp, **r)
    os.replace(tmp, f)
    return r["lam"]


def adapt_project_diag(Xm, Vref, usable, base, lam_vec, lam=5.0, extra=10, iters=10, seed=0):
    """transfer.adapt_project verbatim except the V-step uses a per-protein penalty lam_vec[j]."""
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
            lj = lam_vec[j]
            V[j] = np.linalg.solve(Uo.T @ Uo + lj * I, Uo.T @ Z[o, j] + lj * Vp[j])
    out = base.copy()
    P = U @ V.T
    miss = ~obs & ok[None, :]
    out[miss] = (mu[None, :] + sd[None, :] * P)[miss]
    return out


def _run(Xm, ctx, lam_vec):
    return adapt_project_diag(Xm, ctx["V"], ctx["warm"], ctx["base"], lam_vec)


def impute(Xm, ctx):
    return _run(Xm, ctx, get_lambda(ctx))


def impute_null(Xm, ctx):
    lam = get_lambda(ctx).copy()
    w = np.where(ctx["warm"])[0]
    rng = np.random.default_rng(1000 + ctx["seed"])
    lam[w] = lam[rng.permutation(w)]
    return _run(Xm, ctx, lam)


def impute_nobs(Xm, ctx):
    """Secondary diagnostic (not the formal null): lambda assigned by reference observation-count rank only
    (more evidence -> stronger anchor), same sorted strength distribution."""
    lam = get_lambda(ctx).copy()
    f = CACHE / f"{ctx['target']}_tau.npz"
    nobs = np.load(f)["n_obs_ref"]
    w = np.where(ctx["warm"])[0]
    vals = np.sort(lam[w])
    lam[w[np.argsort(nobs[w], kind="stable")]] = vals
    return _run(Xm, ctx, lam)


# ---- Revision v2 (dose-response): same tau2 cache, contrast squared: lambda_j = c / tau2_j**2, median 100, clip ----
def get_lambda_v2(ctx):
    get_lambda(ctx)
    d = np.load(CACHE / f"{ctx['target']}_tau.npz")
    tau, warm = d["tau2"], d["warm"]
    lam = np.full(len(tau), LAM_MED)
    inv = 1.0 / tau[warm] ** 2
    lam[warm] = np.clip(LAM_MED * inv / np.median(inv), LAM_LO, LAM_HI)
    return lam


def impute_v2(Xm, ctx):
    return _run(Xm, ctx, get_lambda_v2(ctx))


def impute_null_v2(Xm, ctx):
    lam = get_lambda_v2(ctx).copy()
    w = np.where(ctx["warm"])[0]
    rng = np.random.default_rng(1000 + ctx["seed"])
    lam[w] = lam[rng.permutation(w)]
    return _run(Xm, ctx, lam)
