"""Phase-C idea `lwkrige` (queen: representation & mathematics, round 5).

Kriging with a locally re-estimated covariance, credited as an exact correction on top of the cshrink base.

  0. cshrink base, reproduced exactly: psi selected as cshrink.run does (same inner mask, same on-disk score cache),
     (mu, sd) = cshrink_stats(M, psi, psi), (Yc, U, Vad) = adapt_project_uv(...) (verbatim copy of
     nfbase.adapt_project_mb that also returns U and V). Yc == cshrink.impute(M, ctx) bit for bit (unit check).
  1. Block B = the <= 1,500 ok proteins with the most observed cells (stable argsort). Z = (M[:, B] - mu) / sd.
     Sigma_ref = Vad[B] Vad[B]^T + 5 I (the conditional-Gaussian mean under it is the U-step ridge posterior).
     S_target = pairwise-complete covariance of Z (pairwise-centred, nan-masked product formula); entries with < 8
     common cells are filled from c * Sigma_ref, c = sum of valid diag(S) / sum of the same diag(Sigma_ref).
  2. Sigma_d = (1 - d) S_target + d c Sigma_ref; for d < 1 eigenvalues clipped at 1e-3 * mean(eigenvalues).
  3. zhat_d[i, m] = Sigma_d[m, o] Sigma_d[o, o]^{-1} z[i, o] (dense Cholesky per cell; same path for d = 1).
  4. Y = Yc; Y[i, j] += sd_j (zhat_d - zhat_1)[i, j] for missing (i, j in B). d = 1 returns Yc identically.
  5. d chosen on tgtshrink.selection_mask(M, seed) (50/50 MCAR/MNAR inner mask) from D_GRID; on Msel, step 0 uses the
     psi already selected on M and iters = 5; ties go to the larger d (toward the base). Refit on M.

Entry points: impute (selected d), impute_base (d = 1, = cshrink), impute_null (S_target from Zp: every column of Z
permuted across cells independently, values and NaN together, rng(seed + 4242 + j), j = global protein index),
impute_d85 / impute_d85_null (d forced to 0.85), impute_tech (Sigma_ref from technology-matched raw reference loadings,
grid {base} U {(tech, d)}; built on the v2 estimator and grid), impute_allraw (same path with all-reference ctx['V']: content control for tech).

Pluggable: krige_correct(M, ctx, Yc, mu, sd, Vad, sigma_ref_fn=None, d_grid=D_GRID, S_fn=None) -> (Ys, info), where
Ys maps 'base' and every d in d_grid to a full matrix; sigma_ref_fn(B, ctx) -> |B| x |B| PSD matrix. lw_run(...) does
the whole select-and-refit pipeline with the same plug-ins (for mod's `linerna`).
No external data. Cache / logs: data/cache/phasec/lwkrige/.
"""

from __future__ import annotations

import fcntl
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from scp_graphfs_mcp.phasec import tgtshrink as TS
from scp_graphfs_mcp.phasec.cshrink import GRID, INF, _score, cshrink_stats  # noqa: F401
from scp_graphfs_mcp.phasec.nfbase import EXTRA, LAM_V
from scp_graphfs_mcp.phasec.ttprior import META

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "lwkrige"
LAM_U = 5.0
NB = 1500
MIN_PAIR = 8
D_GRID = (0.5, 0.7, 0.85, 0.95, 1.0)


# ---------------------------------------------------------------------------------------------------------------
# 0. adapt_project_mb that also returns U and V (verbatim copy of nfbase.adapt_project_mb otherwise)
def adapt_project_uv(Xm, Vref, usable, base, mu, sd, lam_v=10.0, lam=5.0, extra=3, iters=10, seed=0):
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
    return out, U, V


# ---------------------------------------------------------------------------------------------------------------
# step 0: the cshrink base
def select_psi(M, ctx):
    """Exactly cshrink.run's joint-psi selection (default grid, tag 'c', shared score cache)."""
    Msel, p1, p2 = TS.selection_mask(M, ctx["seed"])
    scores = {p: _score(M, Msel, (p1, p2), ctx, p, p, cshrink_stats, "c") for p in GRID}
    return min(GRID, key=lambda p: (scores[p], p))


def base_fit(M, ctx, psi, iters=10):
    mu, sd, _ = cshrink_stats(M, psi, psi)
    Yc, U, Vad = adapt_project_uv(M, ctx["V"], ctx["warm"], ctx["base"], mu, sd, lam_v=LAM_V, lam=LAM_U,
                                  extra=EXTRA, iters=iters)
    return Yc, mu, sd, Vad


# ---------------------------------------------------------------------------------------------------------------
# steps 1-4
def block(M, n_b=NB):
    with np.errstate(all="ignore"):
        ok = np.isfinite(np.nanmean(M, 0))
    nobs = (~np.isnan(M)).sum(0)
    idx = np.where(ok)[0]
    order = np.argsort(-nobs[idx], kind="stable")
    return idx[order[:n_b]]


def pairwise_cov(Z):
    """Pairwise-complete covariance (each pair centred on its common cells) and the common-cell counts."""
    O = (~np.isnan(Z)).astype(float)
    Z0 = np.where(O > 0, Z, 0.0)
    N = O.T @ O
    G = Z0.T @ Z0
    A = Z0.T @ O  # A[i, j] = sum of z_i over cells where i and j are both observed
    with np.errstate(all="ignore"):
        S = (G - A * A.T / N) / (N - 1)
    S[N < MIN_PAIR] = np.nan
    return S, N


def permute_columns(Z, cols, seed):
    Zp = np.empty_like(Z)
    for k, j in enumerate(cols):
        rng = np.random.default_rng(seed + 4242 + int(j))
        Zp[:, k] = Z[rng.permutation(Z.shape[0]), k]
    return Zp


def ss_dstar(Z, T):
    """Schafer-Strimmer analytic shrinkage intensity toward target T on the zero-filled standardised data
    (diagnostic only; never used for selection)."""
    X = np.where(np.isnan(Z), 0.0, Z)
    X = X - X.mean(0)
    n = X.shape[0]
    Wb = X.T @ X / n
    S = Wb * n / (n - 1)
    num = n / (n - 1) ** 3 * (float(np.sum(np.sum(X ** 2, 1) ** 2)) - n * float(np.sum(Wb ** 2)))
    den = float(np.sum((S - T) ** 2))
    return float(np.clip(num / den, 0.0, 1.0)) if den > 0 else 1.0


def _blend(S, Sr, d):
    Sd = (1.0 - d) * S + d * Sr
    nclip = 0
    if d < 1:
        w, Q = np.linalg.eigh(Sd)
        fl = 1e-3 * float(np.mean(w))
        nclip = int((w < fl).sum())
        Sd = (Q * np.maximum(w, fl)) @ Q.T
        Sd = 0.5 * (Sd + Sd.T)
    return Sd, nclip


def cond_mean(Z, Sd):
    """zhat[i, m] = Sd[m, o] Sd[o, o]^{-1} z[i, o]; observed entries 0."""
    O = ~np.isnan(Z)
    out = np.zeros_like(Z)
    for i in range(Z.shape[0]):
        o, m = O[i], ~O[i]
        if not m.any() or not o.any():
            continue
        L = cho_factor(Sd[np.ix_(o, o)], lower=True)
        out[i, m] = Sd[np.ix_(m, o)] @ cho_solve(L, Z[i, o])
    return out


def _sd_s(sd):
    return np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)


def main_sigma(Vad, B):
    Vb = Vad[B]
    return Vb @ Vb.T + LAM_U * np.eye(len(B))


def completed_cov(Zs, Sr_main):
    """v2 estimator: fill Zs's holes with the base's conditional mean under Sr_main (= the cshrink U-step ridge
    projection), then the ordinary sample covariance (PSD, rank <= n - 1). N = n everywhere."""
    Zf = np.where(np.isnan(Zs), cond_mean(Zs, Sr_main), Zs)
    X = Zf - Zf.mean(0)
    n = X.shape[0]
    return X.T @ X / (n - 1), np.full((Zs.shape[1], Zs.shape[1]), float(n))


def krige_correct(M, ctx, Yc, mu, sd, Vad, sigma_ref_fn=None, d_grid=D_GRID, S_fn=None, s_est="pairwise"):
    """Returns (Ys, info): Ys['base'] = Yc and Ys[d] = Yc + sd (zhat_d - zhat_1^main) on missing entries in B.
    sigma_ref_fn(B, ctx) -> |B| x |B| PSD Sigma_ref for the arm (default: the adapted model Vad Vad^T + 5 I).
    S_fn(Z, B, ctx) -> the matrix S_target is computed from (default: Z itself).
    s_est: 'pairwise' (v1, pairwise-complete covariance) | 'completed' (v2, base-completed sample covariance)."""
    B = block(M)
    sds = _sd_s(sd)
    Z = (M[:, B] - np.nan_to_num(mu)[B]) / sds[B]
    Sr_main = main_sigma(Vad, B)
    Sr = Sr_main if sigma_ref_fn is None else np.asarray(sigma_ref_fn(B, ctx), float)
    Zs = Z if S_fn is None else S_fn(Z, B, ctx)
    if s_est == "completed":
        S, N = completed_cov(Zs, Sr_main)
    else:
        S, N = pairwise_cov(Zs)
    dg = np.diag(S)
    vd = np.isfinite(dg)
    c = float(np.sum(dg[vd]) / np.sum(np.diag(Sr)[vd])) if vd.any() else 1.0
    fill = ~np.isfinite(S)
    S = np.where(fill, c * Sr, S)
    S = 0.5 * (S + S.T)
    cSr = c * Sr
    z1 = cond_mean(Z, c * Sr_main)
    miss = np.isnan(Z)
    Ys = {"base": Yc}
    nclip = {}
    for d in d_grid:
        if d == 1 and sigma_ref_fn is None:
            zd = z1
        else:
            Sd, nclip[d] = _blend(S, cSr, d)
            zd = cond_mean(Z, Sd)
        Y = Yc.copy()
        sub = Y[:, B]
        sub[miss] = (Yc[:, B] + sds[B][None, :] * (zd - z1))[miss]
        Y[:, B] = sub
        Ys[d] = Y
    iu = np.triu_indices(len(B), 1)
    info = {"c": c, "nB": int(len(B)), "frac_filled": float(fill[iu].mean()) if len(iu[0]) else 0.0,
            "median_pair_n": float(np.median(N[iu])) if len(iu[0]) else None,
            "n_clipped": {str(k): v for k, v in nclip.items()}, "ss_dstar": ss_dstar(Zs, cSr)}
    return Ys, info


# ---------------------------------------------------------------------------------------------------------------
# step 5: selection on the 50/50 inner mask, refit on M
def _log(variant, rec):
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / f"log_{variant}.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def _pref(k):
    return (0, 0.0) if k == "base" else (1, -float(k))  # ties: base first, then larger d


def lw_run(M, ctx, variant, sigma_ref_fn=None, d_grid=D_GRID, S_fn=None, forced=None, with_base=False, log=True,
           s_est="pairwise"):
    """Full pipeline. forced: a key of the candidate set ('base' or a d) used without selection."""
    seed = ctx["seed"]
    psi = select_psi(M, ctx)
    keys = (["base"] if with_base else []) + list(d_grid)
    scores = {}
    if forced is None and len(keys) > 1:
        Msel, (r1, c1), (r2, c2) = TS.selection_mask(M, seed)
        Yc_s, mu_s, sd_s, Vad_s = base_fit(Msel, ctx, psi, iters=5)
        Ys_s, _ = krige_correct(Msel, ctx, Yc_s, mu_s, sd_s, Vad_s, sigma_ref_fn, d_grid, S_fn, s_est)
        for k in keys:
            Y = Ys_s[k]
            e1 = np.sqrt(np.mean((Y[r1, c1] - M[r1, c1]) ** 2))
            e2 = np.sqrt(np.mean((Y[r2, c2] - M[r2, c2]) ** 2))
            scores[k] = float(0.5 * (e1 + e2))
        chosen = min(keys, key=lambda k: (scores[k], _pref(k)))
    else:
        chosen = forced if forced is not None else keys[0]
    Yc, mu, sd, Vad = base_fit(M, ctx, psi, iters=10)
    if chosen == "base":
        Ys, info = {"base": Yc}, {}
    else:
        Ys, info = krige_correct(M, ctx, Yc, mu, sd, Vad, sigma_ref_fn, (chosen,), S_fn, s_est)
    if log:
        rec = {"target": ctx.get("target"), "seed": seed, "n_missing": int(np.isnan(M).sum()),
               "psi": "inf" if np.isinf(psi) else psi, "chosen": str(chosen),
               "scores": {str(k): v for k, v in scores.items()}, **info, **ctx.get("_lwlog", {})}
        _log(variant, rec)
    return Ys[chosen]


def null_S_fn(Z, B, ctx):
    return permute_columns(Z, B, ctx["seed"])


def impute(Xm, ctx):
    return lw_run(Xm, ctx, "main")


def impute_base(Xm, ctx):
    return lw_run(Xm, ctx, "base", d_grid=(1.0,), forced=1.0)


def impute_null(Xm, ctx):
    return lw_run(Xm, ctx, "null", S_fn=null_S_fn)


def impute_d85(Xm, ctx):
    return lw_run(Xm, ctx, "d85", d_grid=(0.85,), forced=0.85)


def impute_d85_null(Xm, ctx):
    return lw_run(Xm, ctx, "d85_null", d_grid=(0.85,), forced=0.85, S_fn=null_S_fn)


# v2 (the one allowed revision): v1's pairwise-complete S_target is grossly indefinite at 30 cells (39-43% negative
# eigenvalues, min about -90 vs max about 840 on petrosius2023 / leduc2022), so every d < 1 blew up after clipping.
# v2 estimates S_target from the base-completed standardised block (PSD, no clipping needed) and moves the grid
# toward d = 1, where the inner-mask curve was found to bend (scratch diagnostic on petrosius2023 / leduc2022 seed 100).
D_GRID_V2 = (0.9, 0.95, 0.98, 0.99, 1.0)


def impute_v2(Xm, ctx):
    return lw_run(Xm, ctx, "v2", d_grid=D_GRID_V2, s_est="completed")


def impute_v2_null(Xm, ctx):
    return lw_run(Xm, ctx, "v2_null", d_grid=D_GRID_V2, S_fn=null_S_fn, s_est="completed")


def impute_v2_d98(Xm, ctx):
    return lw_run(Xm, ctx, "v2_d98", d_grid=(0.98,), forced=0.98, s_est="completed")


def impute_v2_d98_null(Xm, ctx):
    return lw_run(Xm, ctx, "v2_d98_null", d_grid=(0.98,), forced=0.98, S_fn=null_S_fn, s_est="completed")


# ---------------------------------------------------------------------------------------------------------------
# technology arm
def tech_loadings(ctx):
    """Raw reference loadings fitted on the reference rows of the target's acquisition class (cached, locked)."""
    from scp_graphfs_mcp import transfer as T

    t = ctx["target"]
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{t}_Vtech.npz"
    with open(CACHE / f".{t}_Vtech.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if f.exists():
            d = np.load(f, allow_pickle=False)
            return d["V"], d["warm"], {"fallback": bool(d["fallback"]), "datasets": [str(x) for x in d["datasets"]]}
        acq = META.get(t, {}).get("acq")
        src = np.asarray(ctx["ref_src"]())
        ds = sorted({s for s in set(src.tolist()) if acq is not None and META.get(s, {}).get("acq") == acq})
        if len(ds) < 2:
            V, warm, fb = ctx["V"], ctx["warm"], True
        else:
            Z = ctx["ref"]()[np.isin(src, ds)]
            V, warm = T.fit_loadings(Z)
            fb = False
        tmp = CACHE / f".{t}_Vtech.tmp{os.getpid()}.npz"
        np.savez_compressed(tmp, V=V, warm=warm, fallback=fb, datasets=np.array(ds, dtype=str), acq=str(acq))
        os.replace(tmp, f)
        return V, warm, {"fallback": fb, "datasets": ds}


def _raw_sigma(V, warm):
    def fn(B, ctx):
        Vb = np.where(warm[B, None], V[B], 0.0)
        rn = np.sum(Vb ** 2, 1)
        return Vb @ Vb.T + np.diag(np.maximum(1.0 - rn, 0.05))
    return fn


def impute_tech(Xm, ctx):
    V, warm, ti = tech_loadings(ctx)
    ctx = dict(ctx, _lwlog={"tech_fallback": ti["fallback"], "tech_datasets": ti["datasets"]})
    return lw_run(Xm, ctx, "tech", sigma_ref_fn=_raw_sigma(V, warm), with_base=True, d_grid=D_GRID_V2,
                  s_est="completed")


def impute_allraw(Xm, ctx):
    return lw_run(Xm, ctx, "allraw", sigma_ref_fn=_raw_sigma(ctx["V"], ctx["warm"]), with_base=True,
                  d_grid=D_GRID_V2, s_est="completed")


# ---------------------------------------------------------------------------------------------------------------
# unit checks
def unit_check(target="petrosius2023"):
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H
    from scp_graphfs_mcp import phase2 as P
    from scp_graphfs_mcp.phasec import cshrink as C
    from scp_graphfs_mcp.phasec import nfbase as NF

    out = {}
    Xs, Xm, mask, names = NF._split(target, "MCAR", 100)
    V, warm = H.target_models(target, names)
    b = P._fill(Xm, P.loo_knn(Xm))
    b = np.where(np.isnan(b), np.nanmean(Xm), b)
    ctx = {"V": V, "warm": warm, "base": b, "seed": 100, "target": target}
    psi = select_psi(Xm, ctx)
    mu, sd, _ = cshrink_stats(Xm, psi, psi)
    y_mb = NF.adapt_project_mb(Xm, V, warm, b, mu, sd, lam_v=LAM_V, extra=EXTRA, iters=10)
    Yc, _, Vad = adapt_project_uv(Xm, V, warm, b, mu, sd, lam_v=LAM_V, lam=LAM_U, extra=EXTRA, iters=10)
    out["U1_uv_vs_mb_maxdiff"] = float(np.max(np.abs(Yc - y_mb)))
    ref = C.run(Xm, ctx, "unit", log=False)
    out["U2_Yc_vs_cshrink_bitexact"] = bool(np.array_equal(Yc, ref))
    out["U2_maxdiff"] = float(np.max(np.abs(Yc - ref)))
    yb = impute_base(Xm, dict(ctx, _lwlog={"unit": True}))
    out["U3_base_vs_cshrink_bitexact"] = bool(np.array_equal(yb, ref))
    # U4: the conditional-Gaussian mean under Vad Vad^T + 5 I equals the ridge U-step projection on the block
    B = block(Xm)
    sds = _sd_s(sd)
    Z = (Xm[:, B] - mu[B]) / sds[B]
    z1 = cond_mean(Z, main_sigma(Vad, B))
    Vb = Vad[B]
    O = ~np.isnan(Z)
    err = 0.0
    for i in range(Z.shape[0]):
        o = O[i]
        u = np.linalg.solve(Vb[o].T @ Vb[o] + LAM_U * np.eye(Vb.shape[1]), Vb[o].T @ Z[i, o])
        err = max(err, float(np.max(np.abs((Vb[~o] @ u) - z1[i, ~o]))) if (~o).any() else 0.0)
    out["U4_condmean_eq_ridge_maxdiff"] = err
    out["pass"] = bool(out["U1_uv_vs_mb_maxdiff"] < 1e-10 and out["U2_Yc_vs_cshrink_bitexact"]
                       and out["U3_base_vs_cshrink_bitexact"] and err < 1e-8)
    out["target"], out["psi"] = target, ("inf" if np.isinf(psi) else psi)
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "unit_check.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    if "--unit" in sys.argv:
        unit_check()
