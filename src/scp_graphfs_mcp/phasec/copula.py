"""Phase-C idea `copula` (queen: representation & mathematics).

Nonparanormal (Gaussian-copula) encoding around the champion's adapted projection (Liu, Lafferty & Wasserman 2009,
JMLR 10:2295-2328; Blom 1958 plotting positions).

  A. Reference: within every reference dataset, each protein with >= 5 observed cells is replaced by its Blom normal
     scores Phi^-1((r - 3/8)/(n + 1/4)) (average-tie ranks), re-standardised to mean 0 / sd 1 over its observed entries
     (location/scale identical to the champion's per-dataset z; only the SHAPE changes). Loadings are refitted with
     ttprior.fit_loadings_weighted (w = 1, cold start, 15 iters = the champion recipe).
  B. Target: champion z (nanmean / nanstd), then the same per-column Blom map (n_j >= 5), then the unchanged
     adapt_project ALS (lam 5, lam_v 100, extra 10, iters 10, seed 0) with the copula prior Vc, then a per-column
     monotone back-transform: np.interp of the predicted latent onto the target's own observed (latent, raw) pairs
     inside the observed latent range, linear extrapolation (OLS slope of raw on latent) outside it; n_j < 5 -> the
     champion's mu + sd * P. Only missing & ok entries are filled; the rest keep ctx['base'].
Both sides use ONE transform function g(z_j; knots of protein k) = Phi^-1(F_k(z_j)), F_k = Blom-smoothed
piecewise-linear ECDF of k's observed values clipped to [1/(2n_k), 1 - 1/(2n_k)]. Real run: k = j (this reproduces the
Blom scores exactly). Null: k = pi(j), a derangement of the eligible proteins of the same dataset (or of the target's
eligible columns): monotone in z_j (ranks / dependence kept), shape from the wrong protein.
Cache: data/cache/phasec/copula/<target>_<variant>.npz (atomic). Precompute:
    python -m scp_graphfs_mcp.phasec.copula --precompute [--targets a,b] [--variants real,null]
"""

from __future__ import annotations

import json
import os
import sys
import zlib
from pathlib import Path

import numpy as np
from scipy.stats import norm, rankdata

from scp_graphfs_mcp import transfer as T
from scp_graphfs_mcp.phasec import ttprior as TT

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "copula"
NMIN = 5
LAM, LAM_V, EXTRA, ITERS, SEED = 5.0, 100.0, 10, 10, 0
FAST = ["petrosius2023", "leduc2022", "weiss2026_liver", "gregoire2023_mixCTRL", "px063590_ADbrain", "sadiku2025"]


def _seed(*parts) -> int:
    return zlib.crc32("|".join(map(str, parts)).encode())


# --------------------------------------------------------------------------------------------- marginal transforms
def knots(v):
    """Blom-smoothed ECDF knots of observed values v: unique sorted values and their Blom probabilities
    (average-tie ranks, so a value's knot probability equals its own Blom position)."""
    n = len(v)
    r = rankdata(v, method="average")
    p = (r - 0.375) / (n + 0.25)
    u, inv = np.unique(v, return_inverse=True)
    pk = np.zeros(len(u))
    np.add.at(pk, inv, p)
    pk /= np.bincount(inv)
    return u, pk, n


def g(z, kn):
    """Phi^-1(F_k(z)), F_k piecewise linear through the knots, clipped to [1/(2n), 1 - 1/(2n)], then re-standardised
    to mean 0 / sd 1 over the given entries (falls back to centring if the sd collapses)."""
    u, pk, n = kn
    p = np.interp(z, u, pk) if len(u) > 1 else np.full(len(z), 0.5)
    p = np.clip(p, 1.0 / (2 * n), 1.0 - 1.0 / (2 * n))
    s = norm.ppf(p)
    s = s - s.mean()
    sd = s.std()
    return s / sd if sd > 1e-9 else s


def derangement(m, rng):
    if m < 2:
        return np.arange(m)
    while True:
        p = rng.permutation(m)
        if not np.any(p == np.arange(m)):
            return p


def encode_block(B, null=False, rng=None, return_map=False):
    """Column-wise copula encoding of a (cells x proteins) block with NaNs. Columns with < NMIN observed are returned
    unchanged. null=True: column j uses the knots of column pi(j) (derangement over eligible columns)."""
    out = B.copy()
    obs = ~np.isnan(B)
    nobs = obs.sum(0)
    elig = np.where(nobs >= NMIN)[0]
    src = elig.copy()
    if null and len(elig) >= 2:
        src = elig[derangement(len(elig), rng)]
    kn = {}
    for j, k in zip(elig, src):
        if k not in kn:
            kn[k] = knots(B[obs[:, k], k])
        o = obs[:, j]
        out[o, j] = g(B[o, j], kn[k])
    if return_map:
        return out, dict(zip(elig.tolist(), src.tolist()))
    return out


# --------------------------------------------------------------------------------------------- reference side
def encode_reference(Z, src, target, null=False):
    Zc = Z.copy()
    for d in sorted(set(src.tolist())):
        rows = np.where(src == d)[0]
        rng = np.random.default_rng(_seed(target, d)) if null else None
        Zc[rows] = encode_block(Z[rows], null=null, rng=rng)
    return Zc


def get_Vc(ctx, variant="real"):
    """variant: 'real' (Blom per protein) or 'null' (deranged per-protein shapes)."""
    target = ctx["target"]
    f = CACHE / f"{target}_{variant}.npz"
    if f.exists():
        return np.load(f)["V"]
    Z = ctx["ref"]()
    src = np.asarray(ctx["ref_src"]())
    Zc = encode_reference(Z, src, target, null=(variant == "null"))
    assert np.array_equal(np.isnan(Zc), np.isnan(Z))
    Vc, warm, _ = TT.fit_loadings_weighted(Zc, np.ones(len(Zc)), ctx["V"], ctx["warm"])
    CACHE.mkdir(parents=True, exist_ok=True)
    tmp = f.with_suffix(f".tmp{os.getpid()}.npz")
    np.savez_compressed(tmp, V=Vc)
    os.replace(tmp, f)
    return Vc


# --------------------------------------------------------------------------------------------- target side
def _standardise(Xm):
    obs = ~np.isnan(Xm)
    mu = np.nanmean(Xm, 0)
    sd = np.nanstd(Xm, 0)
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    ok = np.isfinite(mu)
    Zt = np.where(obs, (Xm - np.nan_to_num(mu)) / sd, np.nan)
    return obs, mu, sd, ok, Zt


def solve(S, obs, ok, Vref, usable, lam=LAM, lam_v=LAM_V, extra=EXTRA, iters=ITERS, seed=SEED):
    """The adapt_project ALS body, verbatim; S = encoded target (NaN where unobserved). Returns P = U V^T."""
    Z = np.where(obs, S, 0.0)
    r0 = Vref.shape[1]
    rng = np.random.default_rng(seed)
    Vp = np.hstack([np.where(usable[:, None], Vref, 0.0), np.zeros((len(Vref), extra))])
    V = Vp + np.hstack([np.zeros_like(Vref), 0.01 * rng.normal(size=(len(Vref), extra))])
    r = r0 + extra
    U = np.zeros((S.shape[0], r))
    I = np.eye(r)
    for _ in range(iters):
        for i in range(S.shape[0]):
            o = obs[i] & ok
            U[i] = np.linalg.solve(V[o].T @ V[o] + lam * I, V[o].T @ Z[i, o])
        for j in np.where(ok)[0]:
            o = obs[:, j]
            Uo = U[o]
            V[j] = np.linalg.solve(Uo.T @ Uo + lam_v * I, Uo.T @ Z[o, j] + lam_v * Vp[j])
    return U @ V.T


def decode_col(p, s, x):
    """Monotone back-transform of latent predictions p using observed (latent s, raw x) pairs of one column."""
    si, inv = np.unique(s, return_inverse=True)
    xi = np.zeros(len(si))
    np.add.at(xi, inv, x)
    xi /= np.bincount(inv)
    if len(si) < 2:
        return np.full(len(p), xi.mean())
    sc = s - s.mean()
    b = float(sc @ (x - x.mean()) / (sc @ sc))
    y = np.interp(p, si, xi)
    lo, hi = p < si[0], p > si[-1]
    y[lo] = xi[0] + b * (p[lo] - si[0])
    y[hi] = xi[-1] + b * (p[hi] - si[-1])
    return y


def project(Xm, ctx, Vref, encode="real"):
    """encode: 'identity' (champion z), 'real' (Blom per column), 'null' (deranged column shapes)."""
    obs, mu, sd, ok, Zt = _standardise(Xm)
    if encode == "identity":
        S = Zt
    else:
        rng = np.random.default_rng(_seed(ctx["target"], "target", ctx.get("seed"))) if encode == "null" else None
        S = encode_block(Zt, null=(encode == "null"), rng=rng)
    P = solve(S, obs, ok, Vref, ctx["warm"])
    out = ctx["base"].copy()
    nobs = obs.sum(0)
    for j in np.where(ok)[0]:
        m = ~obs[:, j]
        if not m.any():
            continue
        if nobs[j] >= NMIN:
            o = obs[:, j]
            out[m, j] = decode_col(P[m, j], S[o, j], Xm[o, j])
        else:
            out[m, j] = mu[j] + sd[j] * P[m, j]
    return out


# --------------------------------------------------------------------------------------------- entry points
def impute(Xm, ctx):
    return project(Xm, ctx, get_Vc(ctx, "real"), "real")


def impute_null(Xm, ctx):
    return project(Xm, ctx, get_Vc(ctx, "null"), "null")


def impute_refonly(Xm, ctx):
    return T.adapt_project(Xm, get_Vc(ctx, "real"), ctx["warm"], ctx["base"], lam_v=LAM_V, extra=EXTRA)


def impute_null_refonly(Xm, ctx):
    return T.adapt_project(Xm, get_Vc(ctx, "null"), ctx["warm"], ctx["base"], lam_v=LAM_V, extra=EXTRA)


def impute_tgtonly(Xm, ctx):
    return project(Xm, ctx, ctx["V"], "real")


def impute_null_tgtonly(Xm, ctx):
    return project(Xm, ctx, ctx["V"], "null")


def impute_identity(Xm, ctx):
    return project(Xm, ctx, ctx["V"], "identity")


# --------------------------------------------------------------------------------------------- checks / diagnostics
def _bench(target):
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from explore_phaseb import load_target, subsample
    from run import make_mask

    ctx = TT._bench_ctx(target)
    X, names, labels = load_target(target)
    return ctx, X, labels, subsample, make_mask


def unit_checks(target="petrosius2023", seed=100):
    from scp_graphfs_mcp import phase2 as P

    ctx, X, labels, subsample, make_mask = _bench(target)
    idx = subsample(X.shape[0], labels, 30, seed)
    Xs = X[idx]
    elig = ((~np.isnan(Xs)).mean(0) >= 0.5) & ~np.isnan(Xs).all(0)
    res = {"target": target, "seed": seed}
    for kind in ("MCAR", "MNAR"):
        mask = make_mask(Xs, elig, kind, seed)
        Xm = Xs.copy()
        Xm[mask] = np.nan
        b = P._fill(Xm, P.loo_knn(Xm))
        b = np.where(np.isnan(b), np.nanmean(Xm), b)
        champ = T.adapt_project(Xm, ctx["V"], ctx["warm"], b, lam_v=LAM_V, extra=EXTRA)
        c = dict(ctx, base=b, seed=seed)
        ident = impute_identity(Xm, c)
        res[f"i_identity_vs_champion_maxabs_{kind}"] = float(np.max(np.abs(ident - champ)))
    Z = ctx["ref"]()
    V1, _, _ = TT.fit_loadings_weighted(Z, np.ones(len(Z)), ctx["V"], ctx["warm"])
    res["ii_fitw1_vs_V_maxabs"] = float(np.max(np.abs(V1 - ctx["V"])))
    # encoding sanity: real g on own knots == Blom scores (restandardised)
    v = np.random.default_rng(0).gamma(2.0, size=50)
    v[:5] = v[5]
    s = norm.ppf((rankdata(v) - 0.375) / (len(v) + 0.25))
    s = (s - s.mean()) / s.std()
    res["blom_g_selfknots_maxabs"] = float(np.max(np.abs(g(v, knots(v)) - s)))
    # null map is monotone in z_j: Spearman of encoded vs input on the reference (one dataset)
    src = np.asarray(ctx["ref_src"]())
    d = sorted(set(src.tolist()))[0]
    B = Z[src == d]
    En = encode_block(B, null=True, rng=np.random.default_rng(_seed(target, d)))
    from scipy.stats import spearmanr

    rhos = []
    for j in np.where((~np.isnan(B)).sum(0) >= NMIN)[0][:300]:
        o = ~np.isnan(B[:, j])
        if np.unique(B[o, j]).size > 2:
            rr = spearmanr(B[o, j], En[o, j]).statistic
            if np.isfinite(rr):
                rhos.append(rr)
    res["null_monotone_spearman_min_mean"] = [float(np.min(rhos)), float(np.mean(rhos)), len(rhos), d]
    return res


def diagnostics(targets=FAST):
    """Mean |Blom score - z| over observed entries (columns with >= NMIN obs), per reference dataset (averaged over the
    targets' protein alignments) and per fast target (its full dev matrix, champion z)."""
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from explore_phaseb import load_target

    per_ds, per_tgt = {}, {}
    for t in targets:
        ctx = TT._bench_ctx(t)
        Z = ctx["ref"]()
        src = np.asarray(ctx["ref_src"]())
        for d in sorted(set(src.tolist())):
            B = Z[src == d]
            E = encode_block(B)
            m = ~np.isnan(B) & ((~np.isnan(B)).sum(0) >= NMIN)[None, :]
            if m.any():
                per_ds.setdefault(d, []).append(float(np.mean(np.abs(E[m] - B[m]))))
        X, _, _ = load_target(t)
        _, _, _, _, Zt = _standardise(X)
        E = encode_block(Zt)
        m = ~np.isnan(Zt) & ((~np.isnan(Zt)).sum(0) >= NMIN)[None, :]
        per_tgt[t] = {"mean_abs_blom_minus_z": float(np.mean(np.abs(E[m] - Zt[m]))),
                      "acq": TT.META.get(t, {}).get("acq", "unknown"), "n_cells": int(X.shape[0])}
        print(t, per_tgt[t], flush=True)
    ds = {d: {"mean_abs_blom_minus_z": float(np.mean(v)), "acq": TT.META.get(d, {}).get("acq", "unknown"),
              "n_target_alignments": len(v)} for d, v in per_ds.items()}
    by_acq = {}
    for d, r in ds.items():
        by_acq.setdefault(r["acq"], []).append(r["mean_abs_blom_minus_z"])
    by_acq = {a: {"mean": float(np.mean(v)), "n_datasets": len(v), "min": float(np.min(v)), "max": float(np.max(v))}
              for a, v in by_acq.items()}
    return {"reference_datasets": ds, "reference_by_acq": by_acq, "fast_targets": per_tgt}


def precompute(targets, variants):
    import time

    for t in targets:
        ctx = TT._bench_ctx(t)
        for v in variants:
            t0 = time.time()
            get_Vc(ctx, v)
            print(f"{t} {v} {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    out_dir = ROOT / "results" / "phaseC"
    if "--unit" in sys.argv:
        r = unit_checks()
        print(json.dumps(r, indent=1))
        CACHE.mkdir(parents=True, exist_ok=True)
        (CACHE / "unit_checks.json").write_text(json.dumps(r, indent=1))
    if "--precompute" in sys.argv:
        vs = sys.argv[sys.argv.index("--variants") + 1].split(",") if "--variants" in sys.argv else ["real", "null"]
        tg = sys.argv[sys.argv.index("--targets") + 1].split(",") if "--targets" in sys.argv else FAST
        precompute(tg, vs)
    if "--diag" in sys.argv:
        r = diagnostics()
        (out_dir / "copula_diag.json").write_text(json.dumps(r, indent=1) + "\n")
        print(json.dumps(r["reference_by_acq"], indent=1))
