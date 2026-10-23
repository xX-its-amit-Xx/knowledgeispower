"""Phase-C idea `tgtshrink` (queen: innovation abduction). Data-only revision of the champion's standardisation.

Fay-Herriot small-area estimation (Fay & Herriot 1979) + limma-trend (Smyth 2004): the target's per-protein 30-cell
nanmean and log nanstd are shrunk toward a WLS regression ("synthetic estimate") on the protein's own abundance and
observed fraction, with the between-protein variance tau2 by method of moments and the shrinkage strength psi chosen per
target on an inner selection mask built only from the input matrix. No external data, no reference baselines.

  0. stats from M: n_j, f_j = n_j / n_cells, xbar_j, s_j (ddof 0), base_ok, cal (n >= 5), ybar = xbar - c_t,
     vy = s^2/n, l = log s, vl = 1/(2(n-1)).
  1. trends (covariates standardised over cal, clipped to the cal range):
       scale    D_l = [1, ns(m), f, f^2, f*m]   (m = standardised xbar)
       location D_y = [1, ns(f)]                (fallback [1, f, f^2] if the f knots are not distinct)
     ns = natural cubic spline (ESL eq. 5.4-5.5 truncated-power form) with boundary knots at the cal min / max and
     interior knots at the 10/50/90% cal quantiles (4 non-constant columns: R's ns with 3 interior knots).
     WLS over cal with weights n - 1.
  2. FH: r = stat - t, tau2 = max(1e-4, var(r[cal]) - psi mean(v[cal])), w = tau2 / (tau2 + psi v), post = t + w r.
     psi = 0 -> the champion.
  3. adapt_project_mb(M, V, warm, base, mu, sd, lam_v=100, extra=10, iters=10).
  4. psi in {1,2,4,8,16} chosen on Msel (10% of eligible observed entries hidden: half uniform, half run.make_mask-MNAR
     weights computed on M; rng seed + 31337); score = mean(RMSE_MCAR-part, RMSE_MNAR-part); fits with iters = 5.
Null (impute_null): the (m, f) covariate rows are permuted jointly across cal proteins (rng seed + 4242) before the
trend fit; non-cal base_ok proteins get the covariates of a random cal protein; v_j stays with its own protein.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.special import expit

from scp_graphfs_mcp.phasec.nfbase import EXTRA, LAM_V, adapt_project_mb

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "tgtshrink"
PSI_GRID = (1.0, 2.0, 4.0, 8.0, 16.0)
MIN_N_CAL = 5


# ---------------------------------------------------------------------------------------------------------------
# 0. target statistics
def target_stats(M):
    obs = ~np.isnan(M)
    n = obs.sum(0)
    with np.errstate(all="ignore"):
        xbar = np.nanmean(M, 0)
        s = np.nanstd(M, 0)
    ok = np.isfinite(xbar)
    base_ok = ok & (n >= 2) & np.isfinite(s) & (s > 1e-6)
    cal = base_ok & (n >= MIN_N_CAL)
    c_t = float(np.median(xbar[ok & (n >= 2)]))
    nn = np.maximum(n, 1)
    return {"n": n, "f": n / M.shape[0], "xbar": xbar, "s": s, "ok": ok, "base_ok": base_ok, "cal": cal,
            "c_t": c_t, "ybar": xbar - c_t, "vy": np.where(base_ok, s ** 2 / nn, np.nan),
            "l": np.where(base_ok, np.log(np.where(base_ok, s, 1.0)), np.nan),
            "vl": np.where(base_ok, 1.0 / (2.0 * np.maximum(n - 1, 1)), np.nan)}


# ---------------------------------------------------------------------------------------------------------------
# 1. designs
def _ns_basis(x, knots):
    """Natural cubic spline basis (no intercept), knots sorted incl. boundary knots; K knots -> K-1 columns."""
    K = len(knots)
    xK, xK1 = knots[-1], knots[-2]

    def d(k):
        return (np.maximum(x - knots[k], 0) ** 3 - np.maximum(x - xK, 0) ** 3) / (xK - knots[k])

    dK1 = d(K - 2)
    cols = [x] + [d(k) - dK1 for k in range(K - 2)]
    B = np.column_stack(cols)
    return B


def _spline_knots(xc):
    q = np.quantile(xc, [0.10, 0.50, 0.90])
    lo, hi = float(np.min(xc)), float(np.max(xc))
    if not (lo < q[0] < q[1] < q[2] < hi):
        return None
    return np.array([lo, *q, hi])


def _designs(m, f, cal_ref, const=False):
    """m, f: raw covariates (per protein, any finite values for rows used); cal_ref: rows defining the standardisation,
    knots and the clipping range. Returns (D_l, D_y, info) for ALL rows."""
    if const:
        one = np.ones((len(m), 1))
        return one, one, {"ns_m": "const", "ns_f": "const"}
    mz_mu, mz_sd = np.mean(m[cal_ref]), np.std(m[cal_ref])
    fz_mu, fz_sd = np.mean(f[cal_ref]), np.std(f[cal_ref])
    mz = (np.nan_to_num(m) - mz_mu) / (mz_sd if mz_sd > 0 else 1.0)
    fz = (np.nan_to_num(f) - fz_mu) / (fz_sd if fz_sd > 0 else 1.0)
    mz = np.clip(mz, mz[cal_ref].min(), mz[cal_ref].max())
    fz = np.clip(fz, fz[cal_ref].min(), fz[cal_ref].max())
    info = {}
    km = _spline_knots(mz[cal_ref])
    if km is None:
        Bm = np.column_stack([mz, mz ** 2, mz ** 3])
        info["ns_m"] = "cubic_poly_fallback"
    else:
        Bm = _ns_basis(mz, km)
        info["ns_m"] = "ns"
    kf = _spline_knots(fz[cal_ref])
    if kf is None:
        Bf = np.column_stack([fz, fz ** 2])
        info["ns_f"] = "quad_fallback"
    else:
        Bf = _ns_basis(fz, kf)
        info["ns_f"] = "ns"
    # rescale spline columns for conditioning (does not change fitted values)
    def _sc(B):
        sd = B[cal_ref].std(0)
        return B / np.where(sd > 0, sd, 1.0)

    one = np.ones((len(m), 1))
    D_l = np.column_stack([one, _sc(Bm), fz, fz ** 2, fz * mz])
    D_y = np.column_stack([one, _sc(Bf)])
    return D_l, D_y, info


def _wls(D, y, w, sel):
    sw = np.sqrt(w[sel])
    beta = np.linalg.lstsq(D[sel] * sw[:, None], y[sel] * sw, rcond=None)[0]
    t = D @ beta
    r = y[sel] - t[sel]
    ybar_w = np.sum(w[sel] * y[sel]) / np.sum(w[sel])
    sst = np.sum(w[sel] * (y[sel] - ybar_w) ** 2)
    r2 = float(1 - np.sum(w[sel] * r ** 2) / sst) if sst > 0 else 0.0
    return t, r2


# ---------------------------------------------------------------------------------------------------------------
# 2. Fay-Herriot shrinkage
def shrink(M, psi, mode="trend", seed=0, cover=None):
    """mode: 'trend' | 'null' (permuted covariates) | 'const' (D = [1], unweighted, restricted to `cover`).
    Returns mu, sd, info."""
    S = target_stats(M)
    base_ok, cal, n = S["base_ok"], S["cal"], S["n"]
    mu, sd = S["xbar"].copy(), S["s"].copy()
    if psi == 0:
        return mu, sd, {"psi": 0.0}
    use = base_ok.copy()
    if cover is not None:
        use &= cover
        cal = cal & cover
    if cal.sum() < 10:
        return mu, sd, {"psi": psi, "n_cal": int(cal.sum())}
    m, f = S["xbar"].copy(), S["f"].astype(float).copy()
    if mode == "null":
        rng = np.random.default_rng(seed + 4242)
        ci = np.where(cal)[0]
        perm = ci[rng.permutation(len(ci))]
        m2, f2 = m.copy(), f.copy()
        m2[ci], f2[ci] = m[perm], f[perm]
        nc = np.where(base_ok & ~cal)[0]
        if len(nc):
            src = ci[rng.integers(0, len(ci), len(nc))]
            m2[nc], f2[nc] = m[src], f[src]
        m, f = m2, f2
    D_l, D_y, info = _designs(m, f, cal, const=(mode == "const"))
    wts = np.where(cal, (n - 1).astype(float), 0.0) if mode != "const" else np.where(cal, 1.0, 0.0)
    info["psi"] = psi
    for key, stat, var, D in (("y", S["ybar"], S["vy"], D_y), ("l", S["l"], S["vl"], D_l)):
        t, r2 = _wls(D, np.nan_to_num(stat), wts, cal)
        r = stat - t
        tau2 = max(1e-4, float(np.var(r[cal]) - psi * np.mean(var[cal])))
        w = tau2 / (tau2 + psi * var)
        post = t + w * r
        info.update({f"tau2_{key}": tau2, f"mean_w_{key}": float(np.mean(w[use])), f"r2_{key}": r2})
        if key == "y":
            mu = np.where(use & np.isfinite(post), S["c_t"] + post, mu)
            info["post_y"] = np.where(use, post, np.nan)
        else:
            sd = np.where(use & np.isfinite(post), np.exp(np.where(use, post, 0.0)), sd)
            info["post_l"] = np.where(use, post, np.nan)
    info["n_cal"] = int(cal.sum())
    info["n_used"] = int(use.sum())
    return mu, sd, info


# ---------------------------------------------------------------------------------------------------------------
# 4. inner selection mask (built only from M)
def selection_mask(M, seed):
    rng = np.random.default_rng(seed + 31337)
    obsM = ~np.isnan(M)
    elig = (obsM.mean(0) >= 0.5) & obsM.any(0)
    r, c = np.nonzero(obsM & elig[None, :])
    ntot = int(round(0.1 * len(r)))
    n1 = ntot // 2
    p1 = rng.choice(len(r), size=n1, replace=False)
    rest = np.setdiff1d(np.arange(len(r)), p1)
    with np.errstate(all="ignore"):
        q25 = np.nanpercentile(M, 25, axis=0)
        s = np.nanstd(M, axis=0)
    w = expit((q25[c[rest]] - M[r[rest], c[rest]]) / (0.5 * np.where(s > 0, s, 1.0)[c[rest]]))
    p2 = rest[rng.choice(len(rest), size=ntot - n1, replace=False, p=w / w.sum())]
    Msel = M.copy()
    Msel[r[p1], c[p1]] = np.nan
    Msel[r[p2], c[p2]] = np.nan
    return Msel, (r[p1], c[p1]), (r[p2], c[p2])


def _log(variant, rec):
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / f"log_{variant}.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def _clean(info):
    return {k: v for k, v in info.items() if not k.startswith("post_")}


def _run(M, ctx, mode, psi=None, variant="main", cover=None):
    V, warm, base, seed = ctx["V"], ctx["warm"], ctx["base"], ctx["seed"]
    scores = {}
    if psi is None:
        Msel, (r1, c1), (r2, c2) = selection_mask(M, seed)
        for p in PSI_GRID:
            mu, sd, _ = shrink(Msel, p, mode, seed, cover)
            Y = adapt_project_mb(Msel, V, warm, base, mu, sd, lam_v=LAM_V, extra=EXTRA, iters=5)
            e1 = np.sqrt(np.mean((Y[r1, c1] - M[r1, c1]) ** 2))
            e2 = np.sqrt(np.mean((Y[r2, c2] - M[r2, c2]) ** 2))
            scores[p] = 0.5 * (e1 + e2)
        psi = min(PSI_GRID, key=lambda p: (scores[p], p))
    mu, sd, info = shrink(M, psi, mode, seed, cover)
    rec = {"target": ctx.get("target"), "mask": "unknown", "seed": seed, "n_missing": int(np.isnan(M).sum()),
           "chosen_psi": psi, "scores": {str(k): float(v) for k, v in scores.items()}, **_clean(info)}
    _log(variant, rec)
    return adapt_project_mb(M, V, warm, base, mu, sd, lam_v=LAM_V, extra=EXTRA)


def impute(Xm, ctx):
    return _run(Xm, ctx, "trend", None, "main")


def impute_null(Xm, ctx):
    return _run(Xm, ctx, "null", None, "null")


def impute_psi4(Xm, ctx):
    return _run(Xm, ctx, "trend", 4.0, "psi4")


def impute_psi4_null(Xm, ctx):
    return _run(Xm, ctx, "null", 4.0, "psi4_null")


def _nfbase_cover(ctx):
    from scp_graphfs_mcp.phasec.nfbase import _prior_cached

    py, _, _ = _prior_cached(ctx["target"], tuple(ctx["names"]), "null")
    return np.isfinite(py)


def impute_const_psi4(Xm, ctx):
    return _run(Xm, ctx, "const", 4.0, "const_psi4", cover=_nfbase_cover(ctx))


# ---------------------------------------------------------------------------------------------------------------
# unit checks
def unit_check():
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H
    from scp_graphfs_mcp import phase2 as P
    from scp_graphfs_mcp import transfer as T
    from scp_graphfs_mcp.phasec import nfbase as NF

    out = {"a_psi0_vs_champion": {}, "b_const_vs_nfbase": {}}
    for t in NF.FAST:
        Xs, Xm, mask, names = NF._split(t, "MCAR", 100)
        V, warm = H.target_models(t, names)
        b = P._fill(Xm, P.loo_knn(Xm))
        b = np.where(np.isnan(b), np.nanmean(Xm), b)
        champ = T.adapt_project(Xm, V, warm, b, lam_v=LAM_V, extra=EXTRA)
        mu, sd, _ = shrink(Xm, 0.0, "trend", 100)
        mine = adapt_project_mb(Xm, V, warm, b, mu, sd, lam_v=LAM_V, extra=EXTRA)
        dmax = float(np.nanmax(np.abs(mine - champ)))
        out["a_psi0_vs_champion"][t] = {"max_abs_diff": dmax, "pass": bool(dmax < 1e-10)}
        # (b) constant trend vs nfbase.baselines with a constant prior, on nfbase-covered proteins
        py, _, _ = NF._prior_cached(t, tuple(names), "null")
        cover = np.isfinite(py)
        c = np.where(cover, 1.0, np.nan)
        mu_nf, sd_nf, _ = NF.baselines(Xm, c, c, psi=4.0)
        mu_c, sd_c, info = shrink(Xm, 4.0, "const", 100, cover)
        S = target_stats(Xm)
        u = cover & S["base_ok"]
        dmu = float(np.max(np.abs(mu_c[u] - mu_nf[u])))
        dsd = float(np.max(np.abs(np.log(sd_c[u]) - np.log(sd_nf[u]))))
        out["b_const_vs_nfbase"][t] = {"max_abs_diff_mu": dmu, "max_abs_diff_logsd": dsd, "n": int(u.sum()),
                                       "pass": bool(max(dmu, dsd) < 1e-8)}
        print(t, out["a_psi0_vs_champion"][t], out["b_const_vs_nfbase"][t], flush=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "unit_check.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    if "--unit" in sys.argv:
        unit_check()
