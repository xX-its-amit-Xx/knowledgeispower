"""Phase-C idea `nfbase` (queen: innovation abduction).

Netflix-Prize baseline predictors for the target's per-protein standardisation (Bell & Koren 2007; Koren 2008).
The champion standardises the target with per-protein nanmean / nanstd of the 30 masked cells and puts predictions
back as mu_j + sd_j * (U V^T)_ij. Here mu_j and sd_j are replaced by empirical-Bayes posteriors that fuse the target's
own estimates with a hierarchical (protein + technology x protein) prior fitted to the UNSTANDARDISED dev reference.

  1. Reference baselines per dataset d (n_dj >= 10): relative location y = mean - median_j(mean), log-scale l = log sd,
     sampling variances s^2/n and 1/(2(n-1)). Datasets with < 200 such proteins are skipped. Cached per target.
  2. Hierarchical EB (separately for y and l): value_dj = a_j + b_{g(d),j} + e_dj, group = ttprior.META[d]['acq'].
  3. Fuse with the target: per-target OLS calibration ybar ~ alpha + beta * prior, tau2 by target-side method of
     moments, BLUP weights w = tau2 / (tau2 + psi * s^2/n). psi = 0 reproduces the champion exactly.
  4. adapt_project_mb = transfer.adapt_project with injected mu / sd.
Nulls: impute_null permutes each reference dataset's baseline vectors among its own covered proteins (rng 4242 + idx);
impute_technull deranges the reference datasets' acquisition labels (ttprior.null_meta_map); the target keeps its own.
No external data. Cache: data/cache/phasec/nfbase/<target>_ref.npz.
"""

from __future__ import annotations

import json
import os
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

from scp_graphfs_mcp import transfer as T
from scp_graphfs_mcp.external import _load
from scp_graphfs_mcp.phasec.ttprior import META, null_meta_map

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "nfbase"
FORBIDDEN = {"derks2022", "wu2026", "brunner2022", "khan2023"}
MIN_N_REF, MIN_PROT_REF, MIN_N_CAL = 10, 200, 5
LAM_V, EXTRA = 100.0, 10


def _reference_datasets(target):
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from explore_phaseb import reference_datasets

    ds = reference_datasets(exclude=(target,))
    assert not set(ds) & (FORBIDDEN | T.PHASE_B_HOLDOUTS | T.HOLDOUTS), "holdout in nfbase reference"
    assert target not in ds
    return ds


# ---------------------------------------------------------------------------------------------------------------
# 0. adapt_project with injectable baselines (verbatim copy of transfer.adapt_project otherwise)
def adapt_project_mb(Xm, Vref, usable, base, mu, sd, lam_v=10.0, lam=5.0, extra=3, iters=10, seed=0):
    obs = ~np.isnan(Xm)
    ok = np.isfinite(np.nanmean(Xm, 0))
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
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


# ---------------------------------------------------------------------------------------------------------------
# 1. reference baselines (independent of Xm; cached per target)
def ref_baselines(target, names):
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{target}_ref.npz"
    if f.exists():
        d = np.load(f, allow_pickle=False)
        return {k: d[k] for k in d.files}
    accs = T._target_accs(names)
    p = len(accs)
    Ys, Ls, VYs, VLs, kept, cs, nprot = [], [], [], [], [], [], []
    for d in _reference_datasets(target):
        X, a = _load(d)
        pos = {x: i for i, x in enumerate(a)}
        B = np.full((X.shape[0], p), np.nan)
        cols = [(j, pos[x]) for j, x in enumerate(accs) if x in pos]
        if cols:
            tj, sj = zip(*cols)
            B[:, list(tj)] = X[:, list(sj)]
        n = (~np.isnan(B)).sum(0)
        good = n >= MIN_N_REF
        if good.sum() < MIN_PROT_REF:
            continue
        with np.errstate(all="ignore"):
            m = np.nanmean(B, 0)
            s = np.nanstd(B, 0)
        good &= np.isfinite(s) & (s > 1e-6)
        c = np.median(m[good])
        y = np.where(good, m - c, np.nan)
        l = np.where(good, np.log(np.where(good, s, 1.0)), np.nan)
        vy = np.where(good, s ** 2 / np.maximum(n, 1), np.nan)
        vl = np.where(good, 1.0 / (2.0 * np.maximum(n - 1, 1)), np.nan)
        Ys.append(y), Ls.append(l), VYs.append(vy), VLs.append(vl), kept.append(d), cs.append(c)
        nprot.append(int(good.sum()))
    out = {"y": np.array(Ys), "l": np.array(Ls), "vy": np.array(VYs), "vl": np.array(VLs),
           "datasets": np.array(kept), "c": np.array(cs), "nprot": np.array(nprot)}
    tmp = CACHE / f".{target}_ref.tmp.npz"
    np.savez_compressed(tmp, **out)
    os.replace(tmp, f)
    return out


# ---------------------------------------------------------------------------------------------------------------
# 2. hierarchical EB fit
def eb_prior(vals, vars_, groups, gstar):
    """vals, vars_: (D x p) with NaN = not covered; groups: (D,) labels. Returns (prior p_j, info dict)."""
    G = sorted(set(groups))
    groups = np.asarray(groups)
    D, p = vals.shape
    cov = np.isfinite(vals)
    s2e = {g: None for g in G}
    theta = {g: np.full(p, np.nan) for g in G}
    Dg = {g: cov[groups == g].sum(0) for g in G}

    def group_means(s2e):
        for g in G:
            rows = groups == g
            v, vv, c = vals[rows], vars_[rows], cov[rows]
            w = np.where(c, 1.0 / ((s2e[g] if s2e[g] is not None else 0.0) + np.where(c, vv, 1.0)), 0.0)
            with np.errstate(all="ignore"):
                theta[g] = np.where(Dg[g] >= 1, (w * np.nan_to_num(v)).sum(0) / w.sum(0), np.nan)

    def within_var():
        out = {}
        for g in G:
            rows = groups == g
            v, vv, c = vals[rows], vars_[rows], cov[rows]
            sel = Dg[g] >= 2
            if not sel.any():
                out[g] = None
                continue
            r2 = np.where(c, (np.nan_to_num(v) - np.nan_to_num(theta[g]))[...] ** 2, 0.0)
            sv = r2.sum(0)[sel] / (Dg[g][sel] - 1)
            mv = np.nanmean(np.where(c, vv, np.nan)[:, sel])
            out[g] = max(1e-4, float(np.mean(sv) - mv))
        pooled = [x for x in out.values() if x is not None]
        fb = float(np.mean(pooled)) if pooled else 1e-2
        return {g: (out[g] if out[g] is not None else fb) for g in G}

    # init: unweighted group means; then iterate (b) -> (a) twice
    group_means({g: 1e12 for g in G})
    for _ in range(2):
        s2e = within_var()
        group_means(s2e)
    TH = np.array([theta[g] for g in G])  # (G x p)
    DG = np.array([Dg[g] for g in G])
    has = DG >= 1
    with np.errstate(all="ignore"):
        a = np.where(has.any(0), np.nansum(np.where(has, TH, 0.0), 0) / has.sum(0), np.nan)
    multi = has.sum(0) >= 2
    if multi.any() and len(G) >= 2:
        bv = np.array([np.var(TH[has[:, j], j], ddof=1) for j in np.where(multi)[0]])
        noise = np.array([np.mean([s2e[G[k]] / DG[k, j] for k in np.where(has[:, j])[0]]) for j in np.where(multi)[0]])
        s2b = max(1e-4, float(np.mean(bv) - np.mean(noise)))
    else:
        s2b = 1e-4
    if gstar in G and len(G) >= 2:
        k = G.index(gstar)
        Ds = DG[k]
        kj = Ds * s2b / (Ds * s2b + s2e[gstar])
        prior = np.where(Ds >= 1, a + kj * (TH[k] - a), a)
    else:
        prior = a.copy()
    info = {"sigma2_e": {g: float(s2e[g]) for g in G}, "sigma2_b": float(s2b),
            "n_prior": int(np.isfinite(prior).sum()), "n_gstar": int((Dg.get(gstar, np.zeros(p)) >= 1).sum())}
    return prior, info


@lru_cache(maxsize=64)
def _prior_cached(target, names_t, mode):
    """mode: 'true' | 'null' (protein-permuted) | 'technull' (deranged acq) | 'flat' | 'flat_null'."""
    R = ref_baselines(target, list(names_t))
    ds = [str(x) for x in R["datasets"]]
    y, l, vy, vl = R["y"].copy(), R["l"].copy(), R["vy"].copy(), R["vl"].copy()
    if mode in ("null", "flat_null"):
        for k, d in enumerate(ds):
            rng = np.random.default_rng(4242 + sorted(ds).index(d))
            cov = np.where(np.isfinite(y[k]))[0]
            perm = cov[rng.permutation(len(cov))]
            for A in (y, l, vy, vl):
                A[k, cov] = A[k, perm]
    if mode == "technull":
        mp = null_meta_map(target, ds)
        groups = [META.get(mp[d], {}).get("acq", "unknown") for d in ds]
    elif mode.startswith("flat"):
        groups = ["all"] * len(ds)
    else:
        groups = [META.get(d, {}).get("acq", "unknown") for d in ds]
    gstar = "all" if mode.startswith("flat") else META.get(target, {}).get("acq", "unknown")
    py, iy = eb_prior(y, vy, groups, gstar)
    pl, il = eb_prior(l, vl, groups, gstar)
    return py, pl, {"loc": iy, "scale": il, "datasets": ds, "groups": groups, "gstar": gstar}


# ---------------------------------------------------------------------------------------------------------------
# 3. fuse with the target
def _fuse(stat, var, prior, cal):
    """stat: target statistic (centred), var: its sampling variance, prior: p_j. Returns posterior, info."""
    sel = cal & np.isfinite(prior) & np.isfinite(stat)
    if sel.sum() < 10:
        return stat.copy(), {"alpha": None, "beta": None, "tau2": None, "mean_w": None, "n_cal": int(sel.sum())}
    A = np.c_[np.ones(sel.sum()), prior[sel]]
    alpha, beta = np.linalg.lstsq(A, stat[sel], rcond=None)[0]
    pm = alpha + beta * prior
    return pm, {"alpha": float(alpha), "beta": float(beta), "n_cal": int(sel.sum()), "sel": sel}


def baselines(Xm, py, pl, psi=1.0, use_loc=True, use_scale=True):
    """Posterior mu, sd for the target (champion's nanmean / nanstd where no prior or n < 2)."""
    obs = ~np.isnan(Xm)
    n = obs.sum(0)
    with np.errstate(all="ignore"):
        xbar = np.nanmean(Xm, 0)
        s = np.nanstd(Xm, 0)
    mu0, sd0 = xbar.copy(), s.copy()
    info = {}
    ok = np.isfinite(xbar)
    base_ok = ok & (n >= 2) & np.isfinite(s) & (s > 1e-6)
    cal = base_ok & (n >= MIN_N_CAL)
    mu, sd = mu0.copy(), sd0.copy()
    if use_loc:
        c_t = np.median(xbar[ok & (n >= 2)])
        ybar = xbar - c_t
        vy = np.where(base_ok, s ** 2 / np.maximum(n, 1), np.nan)
        pm, fi = _fuse(ybar, vy, py, cal)
        if fi.get("alpha") is not None:
            sel = fi.pop("sel")
            tau2 = max(1e-4, float(np.var(ybar[sel] - pm[sel]) - psi * np.mean(vy[sel])))
            w = tau2 / (tau2 + psi * vy)
            post = c_t + w * ybar + (1 - w) * pm
            use = base_ok & np.isfinite(py) & np.isfinite(post)
            mu = np.where(use, post, mu0)
            fi.update(tau2=tau2, mean_w=float(np.mean(w[use])), n_used=int(use.sum()), c_t=float(c_t))
        info["loc"] = fi
    if use_scale:
        lbar = np.where(base_ok, np.log(np.where(base_ok, s, 1.0)), np.nan)
        vl = np.where(base_ok, 1.0 / (2.0 * np.maximum(n - 1, 1)), np.nan)
        pm, fi = _fuse(lbar, vl, pl, cal)
        if fi.get("alpha") is not None:
            sel = fi.pop("sel")
            tau2 = max(1e-4, float(np.var(lbar[sel] - pm[sel]) - psi * np.mean(vl[sel])))
            w = tau2 / (tau2 + psi * vl)
            post = w * lbar + (1 - w) * pm
            use = base_ok & np.isfinite(pl) & np.isfinite(post)
            sd = np.where(use, np.exp(np.where(use, post, 0.0)), sd0)
            fi.update(tau2=tau2, mean_w=float(np.mean(w[use])), n_used=int(use.sum()))
        info["scale"] = fi
    return mu, sd, info


# ---------------------------------------------------------------------------------------------------------------
# 4. candidates
def _run(Xm, ctx, mode="true", psi=1.0, use_loc=True, use_scale=True):
    py, pl, _ = _prior_cached(ctx["target"], tuple(ctx["names"]), mode)
    mu, sd, _ = baselines(Xm, py, pl, psi, use_loc, use_scale)
    return adapt_project_mb(Xm, ctx["V"], ctx["warm"], ctx["base"], mu, sd, lam_v=LAM_V, extra=EXTRA)


def impute(Xm, ctx):
    return _run(Xm, ctx, "true", 1.0)


def impute_null(Xm, ctx):
    return _run(Xm, ctx, "null", 1.0)


def impute_technull(Xm, ctx):
    return _run(Xm, ctx, "technull", 1.0)


def impute_psi4(Xm, ctx):
    return _run(Xm, ctx, "true", 4.0)


def impute_psi4_null(Xm, ctx):
    return _run(Xm, ctx, "null", 4.0)


def impute_flat(Xm, ctx):
    return _run(Xm, ctx, "flat", 1.0)


def impute_flat_null(Xm, ctx):
    return _run(Xm, ctx, "flat_null", 1.0)


def impute_meanonly(Xm, ctx):
    return _run(Xm, ctx, "true", 1.0, use_scale=False)


def impute_meanonly_null(Xm, ctx):
    return _run(Xm, ctx, "null", 1.0, use_scale=False)


def impute_sdonly(Xm, ctx):
    return _run(Xm, ctx, "true", 1.0, use_loc=False)


def impute_sdonly_null(Xm, ctx):
    return _run(Xm, ctx, "null", 1.0, use_loc=False)


# ---------------------------------------------------------------------------------------------------------------
# unit check and go/no-go diagnostic
FAST = ["petrosius2023", "leduc2022", "weiss2026_liver", "gregoire2023_mixCTRL", "px063590_ADbrain", "sadiku2025"]


def _split(target, kind, seed):
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from explore_phaseb import load_target, subsample
    from run import make_mask

    X, names, labels = load_target(target)
    idx = subsample(X.shape[0], labels, 30, seed)
    Xs = X[idx]
    eligible = ((~np.isnan(Xs)).mean(0) >= 0.5) & ~np.isnan(Xs).all(0)
    mask = make_mask(Xs, eligible, kind, seed)
    Xm = Xs.copy()
    Xm[mask] = np.nan
    return Xs, Xm, mask, names


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items() if k != "sel"}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


def unit_check():
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H
    from scp_graphfs_mcp import phase2 as P

    out = {}
    for t in FAST:
        Xs, Xm, mask, names = _split(t, "MCAR", 100)
        V, warm = H.target_models(t, names)
        b = P._fill(Xm, P.loo_knn(Xm))
        b = np.where(np.isnan(b), np.nanmean(Xm), b)
        champ = T.adapt_project(Xm, V, warm, b, lam_v=LAM_V, extra=EXTRA)
        py, pl, _ = _prior_cached(t, tuple(names), "true")
        mu, sd, _ = baselines(Xm, py, pl, psi=0.0)
        mine = adapt_project_mb(Xm, V, warm, b, mu, sd, lam_v=LAM_V, extra=EXTRA)
        out[t] = {"max_abs_diff_psi0": float(np.nanmax(np.abs(mine - champ))),
                  "pass": bool(np.nanmax(np.abs(mine - champ)) < 1e-10)}
        print(t, out[t], flush=True)
    (CACHE / "unit_check.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def diagnostics():
    rows = []
    for t in FAST:
        for kind in ("MCAR", "MNAR"):
            for seed in (100, 101):
                Xs, Xm, mask, names = _split(t, kind, seed)
                r = {"target": t, "mask": kind, "seed": seed}
                ok = np.isfinite(np.nanmean(Xm, 0))
                M = mask & ok[None, :]
                truth = Xs[M]
                cols = np.nonzero(M)[1]
                mu0 = np.nanmean(Xm, 0)
                r["rmse_meanonly_champion"] = float(np.sqrt(np.mean((mu0[cols] - truth) ** 2)))
                pcols = np.where(M.any(0))[0]
                mt = np.array([np.mean(Xs[M[:, j], j]) for j in pcols])
                r["mae_mu_vs_masked_mean_champion"] = float(np.mean(np.abs(mu0[pcols] - mt)))
                for mode, psi in (("true", 1.0), ("true", 4.0), ("null", 1.0), ("technull", 1.0), ("flat", 1.0)):
                    py, pl, pinfo = _prior_cached(t, tuple(names), mode)
                    mu, sd, info = baselines(Xm, py, pl, psi)
                    tag = f"{mode}_psi{psi:g}"
                    r[f"rmse_meanonly_{tag}"] = float(np.sqrt(np.mean((mu[cols] - truth) ** 2)))
                    r[f"mae_mu_vs_masked_mean_{tag}"] = float(np.mean(np.abs(mu[pcols] - mt)))
                    r[f"info_{tag}"] = _jsonable(info)
                    if kind == "MCAR" and seed == 100:
                        r[f"eb_{tag}"] = _jsonable({k: pinfo[k] for k in ("loc", "scale", "groups", "gstar")})
                rows.append(r)
                print(t, kind, seed, {k: round(v, 4) for k, v in r.items() if k.startswith("rmse")}, flush=True)
    summ = {}
    keys = [k for k in rows[0] if k.startswith("rmse_") or k.startswith("mae_")]
    for t in FAST:
        rr = [r for r in rows if r["target"] == t]
        summ[t] = {k: float(np.mean([r[k] for r in rr])) for k in keys}
    summ["ALL"] = {k: float(np.mean([summ[t][k] for t in FAST])) for k in keys}
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "diagnostics.json").write_text(json.dumps({"summary": summ, "rows": rows}, indent=1) + "\n")
    return summ


def feasibility():
    """Pairwise Pearson r of per-protein relative location / log-sd between reference datasets (dev only)."""
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from explore_phaseb import load_target

    _, names, _ = load_target("petrosius2023")
    R = ref_baselines("petrosius2023", names)
    ds = [str(d) for d in R["datasets"]]
    out = {"datasets": ds, "acq": [META.get(d, {}).get("acq", "unknown") for d in ds], "pairs": []}
    for i in range(len(ds)):
        for j in range(i + 1, len(ds)):
            row = {"a": ds[i], "b": ds[j]}
            for key in ("y", "l"):
                s = np.isfinite(R[key][i]) & np.isfinite(R[key][j])
                row[f"r_{key}"] = float(np.corrcoef(R[key][i][s], R[key][j][s])[0, 1]) if s.sum() >= 50 else None
                row["n"] = int(s.sum())
            out["pairs"].append(row)
    out["sd_rel_location"] = {d: float(np.nanstd(R["y"][k])) for k, d in enumerate(ds)}
    (CACHE / "feasibility.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    a = sys.argv
    if "--unit" in a:
        unit_check()
    if "--diag" in a:
        print(json.dumps(diagnostics()["ALL"], indent=1))
    if "--feas" in a:
        feasibility()


# attribution at psi = 4 (added after the psi4 screen; reported alongside the a-priori psi = 1 attribution)
def impute_meanonly_psi4(Xm, ctx):
    return _run(Xm, ctx, "true", 4.0, use_scale=False)


def impute_meanonly_psi4_null(Xm, ctx):
    return _run(Xm, ctx, "null", 4.0, use_scale=False)


def impute_sdonly_psi4(Xm, ctx):
    return _run(Xm, ctx, "true", 4.0, use_loc=False)


def impute_sdonly_psi4_null(Xm, ctx):
    return _run(Xm, ctx, "null", 4.0, use_loc=False)
