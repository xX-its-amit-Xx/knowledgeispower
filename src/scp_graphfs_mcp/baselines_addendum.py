"""v9 ADDENDUM comparators (post hoc, exploratory; written 2026-11-08 after PREREGISTRATION_v9 was frozen).

Why: the independent review of the frozen panel (docs/v9_review.md) found that `baselines.bpca` runs Bayesian PCA with the
~30 cells as observations, so ARD prunes every component and it equals per-protein mean imputation on 20/24 dev splits;
and the literature check found that several of the most-used proteomics imputers were missing. These methods are scored
by benchmarks/addendum_v9.py on the same splits as the frozen panel. They change no preregistered verdict.

Same signature as baselines.py: fn(Xm, ctx) -> full matrix; only unobserved entries are scored.
  bpca_oba      BPCA as in Oba et al. 2003 / pcaMethods / MSnbase impute("bpca"): PROTEINS are the observations, the
                cells are the dimensions, ARD on k = n-1 components, missing entries filled from the posterior.
  lls           local least squares (Kim, Golub & Park 2005): for each protein with missing values, regress it on its
                k most correlated proteins over the cells where it is observed (k = 20).
  knn_protein   kNN over proteins (impute.knn / Troyanskaya 2001 orientation): average of the k = 10 nearest proteins
                (Euclidean on co-observed cells), each standardised.
  qrilc         QRILC (Lazar et al. 2016, imputeLCMD): per cell, fit a normal to the observed quantiles assuming the
                missing fraction is the left tail, and draw missing values from the truncated left tail.
  perseus       Perseus down-shift imputation (Tyanova et al. 2016): per cell, draws from N(mean - 1.8 sd, (0.3 sd)^2).
  missforest_full  missForest closer to its defaults (100 trees, up to 10 iterations, 50 nearest features) - CPU heavy,
                dev check only.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from .baselines import _iterative, _keep_obs, _stats


def bpca_oba(Xm, ctx, iters=200, tol=1e-6):
    """Variational-free BPCA as in Oba's BPCAfill: rows = proteins (observations), columns = cells (dimensions)."""
    obs, mu, sd, ok = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, np.nan)[:, ok].T  # proteins x cells
    O = ~np.isnan(Z)
    N, D = Z.shape
    k = D - 1
    Y = np.where(O, Z, 0.0)
    m = np.nanmean(Z, 0)
    m = np.where(np.isfinite(m), m, 0.0)
    rng = np.random.default_rng(ctx["seed"] + 11)
    W = 0.1 * rng.normal(size=(D, k))
    tau = 1.0
    alpha = np.ones(k)
    last = np.inf
    for _ in range(iters):
        Ex = np.zeros((N, k))
        Sx_sum = np.zeros((k, k))
        Yfill = np.where(O, Z, 0.0)
        for i in range(N):
            o = O[i]
            Wo = W[o]
            M = np.eye(k) + tau * Wo.T @ Wo
            Minv = np.linalg.inv(M)
            Ex[i] = tau * Minv @ (Wo.T @ (Z[i, o] - m[o]))
            Sx_sum += Minv
            Yfill[i, ~o] = m[~o] + W[~o] @ Ex[i]
        m = Yfill.mean(0) - W @ Ex.mean(0)
        Yc = Yfill - m
        A = Ex.T @ Ex + Sx_sum + np.diag(alpha) / tau
        W = np.linalg.solve(A, Ex.T @ Yc).T
        resid = Yc - Ex @ W.T
        tau = (N * D) / max(np.sum(resid ** 2) + np.trace(W @ Sx_sum @ W.T), 1e-12)
        alpha = D / np.maximum((W ** 2).sum(0), 1e-12)
        err = float(np.sum(np.where(O, resid, 0.0) ** 2))
        if abs(last - err) <= tol * max(err, 1e-12):
            break
        last = err
    L = np.zeros(Xm.shape)
    L[:, ok] = (m + Ex @ W.T).T
    return _keep_obs(Xm, mu + sd * L, ctx)


def _protein_corr(Z, O, min_overlap=5):
    Zf = np.where(O, Z, 0.0)
    Of = O.astype(float)
    cnt = Of.T @ Of
    C = (Zf.T @ Zf) / np.maximum(cnt, 1)
    C[cnt < min_overlap] = 0.0
    np.fill_diagonal(C, 0.0)
    return C, cnt


def lls(Xm, ctx, k=10, ridge=1.0):
    """k neighbours by |correlation|, capped at a third of the protein's observed cells (Kim et al. choose k by the data;
    with ~30 cells an uncapped k interpolates noise); ridge on standardised values."""
    obs, mu, sd, ok = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, 0.0)
    C, _ = _protein_corr(Z, obs)
    out = np.zeros_like(Z)
    full = np.where(obs, Z, 0.0)  # neighbours' own missing values = 0 (their mean) when used as predictors
    for j in np.where((~obs).any(0) & ok)[0]:
        kk = max(2, min(k, int(obs[:, j].sum()) // 3))
        nb = np.argsort(-np.abs(C[j]))[:kk]
        nb = nb[np.abs(C[j, nb]) > 0]
        if len(nb) == 0:
            continue
        o = obs[:, j]
        A = full[o][:, nb]
        beta = np.linalg.solve(A.T @ A + ridge * np.eye(len(nb)), A.T @ Z[o, j])
        out[~o, j] = full[~o][:, nb] @ beta
    return _keep_obs(Xm, mu + sd * out, ctx)


def knn_protein(Xm, ctx, k=10, min_overlap=5):
    obs, mu, sd, ok = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, 0.0)
    Of = obs.astype(float)
    sq = (Z ** 2).T @ Of
    cross = Z.T @ Z
    cnt = Of.T @ Of
    D2 = (sq + sq.T - 2 * cross) / np.maximum(cnt, 1)  # mean squared distance over co-observed cells
    D2[cnt < min_overlap] = np.inf
    np.fill_diagonal(D2, np.inf)
    out = np.zeros_like(Z)
    for j in np.where((~obs).any(0) & ok)[0]:
        for i in np.where(~obs[:, j])[0]:
            cand = np.where(obs[i] & np.isfinite(D2[j]))[0]
            if len(cand) == 0:
                continue
            nb = cand[np.argsort(D2[j, cand])[:k]]
            out[i, j] = Z[i, nb].mean()
    return _keep_obs(Xm, mu + sd * out, ctx)


def qrilc(Xm, ctx, tune_sigma=1.0):
    """imputeLCMD::impute.QRILC: for each cell, the observed fraction q = 1 - pmiss; fit mean/sd of a normal by regressing
    the observed values' quantiles on the normal quantiles (pmiss + q * ppoints), then draw missing values from the
    normal truncated above at the pmiss quantile."""
    rng = np.random.default_rng(ctx["seed"] + 203)
    Y = np.array(Xm, float)
    for i in range(Xm.shape[0]):
        x = Xm[i]
        o = ~np.isnan(x)
        pm = 1 - o.mean()
        if o.sum() < 5 or pm == 0:
            continue
        xs = np.sort(x[o])
        pp = (np.arange(1, len(xs) + 1) - 0.5) / len(xs)
        zq = stats.norm.ppf(pm + (1 - pm) * pp)
        slope, icept = np.polyfit(zq, xs, 1)
        upper = stats.norm.ppf(max(pm, 1e-3))
        draws = stats.truncnorm.rvs(-np.inf, upper, loc=0, scale=1, size=(~o).sum(), random_state=rng)
        Y[i, ~o] = icept + slope * tune_sigma * draws
    return _keep_obs(Xm, Y, ctx)


def perseus(Xm, ctx, shift=1.8, width=0.3):
    rng = np.random.default_rng(ctx["seed"] + 307)
    mu_c = np.nanmean(Xm, 1)
    sd_c = np.nanstd(Xm, 1)
    Y = (mu_c - shift * sd_c)[:, None] + (width * sd_c)[:, None] * rng.normal(size=Xm.shape)
    return _keep_obs(Xm, Y, ctx)


def missforest_full(Xm, ctx):
    from sklearn.ensemble import RandomForestRegressor

    rf = RandomForestRegressor(n_estimators=100, max_features=0.33, random_state=ctx["seed"], n_jobs=1)
    return _iterative(Xm, ctx, rf, max_iter=10)


ADDENDUM = {"bpca_oba": bpca_oba, "lls": lls, "knn_protein": knn_protein, "qrilc": qrilc, "perseus": perseus}
HEAVY = {"missforest_full": missforest_full}
