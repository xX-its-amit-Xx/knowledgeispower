"""v9 ADDENDUM comparators (post hoc, exploratory; written 2026-11-08 after PREREGISTRATION_v9 was frozen).

Why: the independent review of the frozen panel (docs/v9_review.md) found that `baselines.bpca` runs Bayesian PCA with the
~30 cells as observations, so ARD prunes every component and it equals per-protein mean imputation on 20/24 dev splits;
and the literature check found that several of the most-used proteomics imputers were missing. These methods are scored
by benchmarks/addendum_v9.py on the same splits as the frozen panel. They change no preregistered verdict.

Same signature as baselines.py: fn(Xm, ctx) -> full matrix; only unobserved entries are scored.
  bpca_oba      EM-ARD Bayesian PCA in Oba et al.'s orientation (pcaMethods / MSnbase impute("bpca")): PROTEINS are the
                observations, the cells the dimensions, k = n-1 components with ARD; exact missing-data EM (see docstring).
  lls           local least squares (Kim, Golub & Park 2005): for each protein with missing values, regress it on its
                k most correlated proteins over the cells where it is observed (k = min(10, n_obs // 3), ridge 1),
                iterated up to 5 passes with neighbours' missing values re-filled from the previous pass.
  knn_protein   kNN over proteins (impute.knn / Troyanskaya 2001 orientation): average of the k = 10 nearest proteins
                (Euclidean on co-observed cells), each standardised.
  qrilc         QRILC (Lazar et al. 2016 model): per cell, fit a normal to the observed quantiles assuming the missing
                fraction is the left tail, and draw from the truncated left tail with SD = the fitted slope. (imputeLCMD's
                code passes the variance where an SD is expected, i.e. effectively SD = sqrt(slope); we follow the paper.)
  perseus       Perseus down-shift imputation (Tyanova et al. 2016): per cell, draws from N(mean - 1.8 sd, (0.3 sd)^2).
  missforest_full  missForest closer to its defaults (100 trees, up to 10 iterations, 50 nearest features) - CPU heavy,
                dev check only.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

from .baselines import _iterative, _keep_obs, _stats


def bpca_oba(Xm, ctx, iters=200, tol=1e-6):
    """EM-ARD probabilistic PCA (Bishop 1999 ARD) in Oba et al.'s orientation: PROTEINS are the observations and the
    cells are the dimensions (as in pcaMethods/MSnbase impute("bpca")), k = D-1 components, exact missing-data EM:
    per-protein posterior over the latent from its observed cells; per-cell M-step over the proteins that observe it
    (so no imputed value enters the sufficient statistics); noise precision from observed entries plus the posterior
    variance term. Not a line-by-line port of Oba's variational BPCAfill (no SigW term, Bishop's ARD update)."""
    obs, mu, sd, ok = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, np.nan)[:, ok].T  # proteins x cells
    O = ~np.isnan(Z)
    Of = O.astype(float)
    Z0 = np.where(O, Z, 0.0)
    N, D = Z.shape
    k = D - 1
    m = np.where(Of.sum(0) > 0, Z0.sum(0) / np.maximum(Of.sum(0), 1), 0.0)
    rng = np.random.default_rng(ctx["seed"] + 11)
    W = 0.1 * rng.normal(size=(D, k))
    tau = 1.0
    alpha = np.ones(k)
    last = np.inf
    Ex = np.zeros((N, k))
    Sx = np.zeros((N, k, k))
    for _ in range(iters):
        for i in range(N):  # E-step: posterior of x_i from protein i's observed cells
            o = O[i]
            Wo = W[o]
            Minv = np.linalg.inv(np.eye(k) + tau * Wo.T @ Wo)
            Ex[i] = tau * Minv @ (Wo.T @ (Z0[i, o] - m[o]))
            Sx[i] = Minv
        XX = Sx + np.einsum("ik,il->ikl", Ex, Ex)  # E[x x'] per protein
        Sj = np.einsum("ij,ikl->jkl", Of, XX)  # per cell: sum over proteins observing it
        for j in range(D):  # M-step per cell over observing proteins
            o = O[:, j]
            if not o.any():
                continue
            b = Ex[o].T @ (Z0[o, j] - m[j])
            W[j] = np.linalg.solve(Sj[j] + np.diag(alpha) / tau, b)
            m[j] = np.mean(Z0[o, j] - Ex[o] @ W[j])
        R = np.where(O, Z0 - m - Ex @ W.T, 0.0)
        var_term = np.einsum("jk,ikl,jl,ij->", W, Sx, W, Of)
        tau = Of.sum() / max(np.sum(R ** 2) + var_term, 1e-12)
        alpha = D / np.maximum((W ** 2).sum(0), 1e-12)
        err = float(np.sum(R ** 2))
        if abs(last - err) <= tol * max(err, 1e-12):
            break
        last = err
    for i in range(N):  # final E-step with the final parameters
        o = O[i]
        Wo = W[o]
        Minv = np.linalg.inv(np.eye(k) + tau * Wo.T @ Wo)
        Ex[i] = tau * Minv @ (Wo.T @ (Z0[i, o] - m[o]))
    L = np.zeros(Xm.shape)
    L[:, ok] = (m + Ex @ W.T).T
    return _keep_obs(Xm, mu + sd * L, ctx)


def _protein_corr(Z, O, min_overlap=None):
    """Co-observed cross-moment of per-protein standardised values (~ correlation), shrunk by overlap cnt/(cnt+5);
    pairs observed together in fewer than max(5, n/4) cells get 0."""
    Zf = np.where(O, Z, 0.0)
    Of = O.astype(float)
    cnt = Of.T @ Of
    C = (Zf.T @ Zf) / np.maximum(cnt, 1) * (cnt / (cnt + 5.0))
    C[cnt < (min_overlap or max(5, Z.shape[0] // 4))] = 0.0
    np.fill_diagonal(C, 0.0)
    return C, cnt


def lls(Xm, ctx, k=10, ridge=1.0, passes=5):
    """k neighbours by |correlation|, capped at a third of the protein's observed cells (Kim et al. choose k by the data;
    with ~30 cells an uncapped k interpolates noise); ridge on standardised values."""
    obs, mu, sd, ok = _stats(Xm)
    Z = np.where(obs, (Xm - mu) / sd, 0.0)
    C, _ = _protein_corr(Z, obs)
    out = np.zeros_like(Z)
    full = np.where(obs, Z, 0.0)  # pass 1: neighbours' missing values = 0 (their mean); then refined iteratively
    for _ in range(passes):
        out = _lls_pass(Z, obs, ok, C, full, k, ridge)
        new = np.where(obs, Z, out)
        if np.linalg.norm(new - full) <= 1e-4 * max(np.linalg.norm(full), 1e-12):
            full = new
            break
        full = new
    return _keep_obs(Xm, mu + sd * out, ctx)


def _lls_pass(Z, obs, ok, C, full, k, ridge):
    out = np.zeros_like(Z)
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
    return out


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
