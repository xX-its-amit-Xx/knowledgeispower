"""Phase-C idea `otanalog` (queen: representation & mathematics, round 6).

Analog forecasting (Lorenz 1969) by semi-relaxed entropic optimal transport (Cuturi 2013; Chizat, Peyre, Schmitzer &
Vialard, arXiv:1607.05816) from the 30 target cells to real reference cells, as a residual correction on the cshrink base.

  0. Base (identical to cshrink): psi = lwkrige.select_psi(M, ctx); (mu, sd) = cshrink_stats(M, psi, psi);
     (Yc, U, Vad) = lwkrige.adapt_project_uv(..., lam_v=100, lam=5, extra=10, iters=10). Zhat = U Vad^T (z units).
  1. Reference subsample: ctx['ref']() rows (per-dataset z-scored), <= 400 cells per dataset, <= 5,000 in total,
     rng(12345 + crc32(target)); cached per target in a module dict and on disk.
  2. Embedding: W = Vad[:, :20] (reference-factor columns), J = ok & warm; e = (W_o^T W_o + 5 I)^{-1} W_o^T z_o on the
     cell's observed proteins in J, the same encoder for target (Zt) and reference rows; reference cells with |o| < 20
     dropped. C_ik = ||e_i - e_k||^2, c_med = median(C).
  3. Semi-relaxed entropic OT in the log domain (target marginal exact, reference marginal KL-relaxed with strength
     rho; rho = inf is balanced), a = 1/n, b = 1/m, <= 500 iterations or row-marginal L1 error < 1e-6.
     Grid eps in {0.05, 0.1, 0.5} c_med x rho in {0.1, 1, inf} c_med.
  4. Analog values yhat_ij = sum_k P_ik O_kj R_kj / sum_k P_ik O_kj where support_ij = sum_k P_ik O_kj / sum_k P_ik
     >= 0.2, else Zhat_ij. corr_ij = sd_j (yhat_ij - Zhat_ij). Y = Yc + gamma corr on missing & ok entries.
  5. (eps, rho, gamma) chosen on tgtshrink.selection_mask(M, seed) (50/50 MCAR/MNAR inner mask): steps 0-4 on Msel with
     the psi selected on M and iters = 5; score = 0.5 (RMSE_p1 + RMSE_p2); ties -> smaller gamma, larger eps, larger rho.
     Refit on M with the chosen plan; gamma = 0 returns Yc (= cshrink.impute bit for bit).
  6. Diagnostics -> data/cache/phasec/otanalog/log_<variant>.jsonl (one line per call).

Entry points: impute (main), impute_null (each row of C permuted over the reference columns, rng(seed + 4049), before
step 3), impute_base (gamma = 0, = cshrink), impute_bal (rho = inf only). numpy only (POT is not installed).
"""

from __future__ import annotations

import json
import os
import sys
import zlib
from pathlib import Path

import numpy as np

from scp_graphfs_mcp.phasec import tgtshrink as TS
from scp_graphfs_mcp.phasec.cshrink import cshrink_stats
from scp_graphfs_mcp.phasec.lwkrige import adapt_project_uv, select_psi
from scp_graphfs_mcp.phasec.nfbase import EXTRA, LAM_V
from scp_graphfs_mcp.phasec.ttprior import META

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "otanalog"
LAM_U = 5.0
R_REF = 20
PER_DS = 400
MAX_REF = 5000
MIN_OBS = 20
RIDGE_E = 5.0
EPS_GRID = (0.05, 0.1, 0.5)
RHO_GRID = (0.1, 1.0, float("inf"))
GAMMA_GRID = (0.0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0)
MIN_SUPPORT = 0.2
FIXED = (0.1, 1.0)
MAX_IT = 500
TOL = 1e-6
_REF = {}


# ---------------------------------------------------------------------------------------------------------------
# 1. reference subsample
def ref_subsample(ctx):
    t = ctx["target"]
    if t in _REF:
        return _REF[t]
    f = CACHE / f"{t}_refsub.npz"
    if f.exists():
        d = np.load(f, allow_pickle=False)
        out = (d["R"], np.array([str(x) for x in d["src"]]))
    else:
        R = ctx["ref"]()
        src = np.asarray(ctx["ref_src"]()).astype(str)
        rng = np.random.default_rng(12345 + zlib.crc32(t.encode()))
        keep = []
        for ds in sorted(set(src.tolist())):
            idx = np.where(src == ds)[0]
            if len(idx) > PER_DS:
                idx = np.sort(rng.choice(idx, PER_DS, replace=False))
            keep.append(idx)
        keep = np.concatenate(keep)
        if len(keep) > MAX_REF:
            keep = np.sort(rng.choice(keep, MAX_REF, replace=False))
        out = (R[keep].astype(np.float32), src[keep])
        CACHE.mkdir(parents=True, exist_ok=True)
        tmp = CACHE / f".{t}_refsub.tmp{os.getpid()}.npz"
        np.savez_compressed(tmp, R=out[0], src=out[1])
        os.replace(tmp, f)
    _REF[t] = out
    return out


# ---------------------------------------------------------------------------------------------------------------
# 2. embedding
def embed(Zrows, W, J, chunk=400):
    """Ridge projection of each row on W over its observed proteins in J. Returns (E, n_obs)."""
    WJ = W[J].astype(np.float64)
    r = WJ.shape[1]
    Wout = (WJ[:, :, None] * WJ[:, None, :]).reshape(len(WJ), r * r)
    E = np.zeros((Zrows.shape[0], r))
    nobs = np.zeros(Zrows.shape[0], int)
    I = RIDGE_E * np.eye(r)
    for s in range(0, Zrows.shape[0], chunk):
        Zc = np.asarray(Zrows[s:s + chunk][:, J], np.float64)
        O = np.isfinite(Zc)
        Z0 = np.where(O, Zc, 0.0)
        G = (O.astype(np.float64) @ Wout).reshape(-1, r, r) + I
        rhs = Z0 @ WJ
        E[s:s + chunk] = np.linalg.solve(G, rhs[:, :, None])[:, :, 0]
        nobs[s:s + chunk] = O.sum(1)
    return E, nobs


def sqdist(A, B):
    C = (A ** 2).sum(1)[:, None] + (B ** 2).sum(1)[None, :] - 2 * A @ B.T
    return np.maximum(C, 0.0)


# ---------------------------------------------------------------------------------------------------------------
# 3. semi-relaxed entropic OT (log domain)
def _lse(X, axis):
    m = X.max(axis=axis, keepdims=True)
    return (m + np.log(np.exp(X - m).sum(axis=axis, keepdims=True))).squeeze(axis)


def sinkhorn(C, eps, rho, max_it=MAX_IT, tol=TOL):
    n, m = C.shape
    la, lb = np.log(1.0 / n), np.log(1.0 / m)
    kappa = 1.0 if np.isinf(rho) else rho / (rho + eps)
    f, g = np.zeros(n), np.zeros(m)
    err, it = np.inf, 0
    for it in range(1, max_it + 1):
        f = eps * la - eps * _lse((g[None, :] - C) / eps, 1)
        g = kappa * (eps * lb - eps * _lse((f[:, None] - C) / eps, 0))
        rows = np.exp(_lse((f[:, None] + g[None, :] - C) / eps, 1))
        err = float(np.abs(rows - 1.0 / n).sum())
        if err < tol:
            break
    P = np.exp((f[:, None] + g[None, :] - C) / eps)
    return P, {"iters": it, "row_err": err}


# ---------------------------------------------------------------------------------------------------------------
# 0 + 2: base fit and embeddings for one input matrix
def _sd_s(sd):
    return np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)


def prepare(M, ctx, psi, iters):
    mu, sd, _ = cshrink_stats(M, psi, psi)
    Yc, U, Vad = adapt_project_uv(M, ctx["V"], ctx["warm"], ctx["base"], mu, sd, lam_v=LAM_V, lam=LAM_U,
                                  extra=EXTRA, iters=iters)
    with np.errstate(all="ignore"):
        ok = np.isfinite(np.nanmean(M, 0))
    sds = _sd_s(sd)
    Zt = (M - np.nan_to_num(mu)) / sds
    Zhat = U @ Vad.T
    R, src = ref_subsample(ctx)
    J = ok & np.asarray(ctx["warm"], bool)
    W = Vad[:, :R_REF]
    Et, _ = embed(Zt, W, J)
    Er, nr = embed(R, W, J)
    keep = nr >= MIN_OBS
    C = sqdist(Et, Er[keep])
    return {"Yc": Yc, "Zhat": Zhat, "sd": sds, "ok": ok, "C": C, "Rk": R[keep], "src": src[keep],
            "c_med": float(np.median(C)), "n_ref": int(keep.sum()), "n_ref_dropped": int((~keep).sum())}


def analog_corr(P, st, cols=None):
    """corr (n x p, 0 where not corrected) on missing & ok entries; plus support and covered mask."""
    R = st["Rk"]
    ok = st["ok"]
    if cols is None:
        cols = np.where(ok)[0]
    Rc = R[:, cols]
    O = np.isfinite(Rc).astype(np.float32)
    R0 = np.where(O > 0, Rc, 0.0).astype(np.float32)
    P32 = P.astype(np.float32)
    num = (P32 @ R0).astype(np.float64)
    den = (P32 @ O).astype(np.float64)
    tot = P.sum(1)[:, None]
    support = den / np.maximum(tot, 1e-300)
    with np.errstate(all="ignore"):
        yhat = num / den
    use = support >= MIN_SUPPORT
    Zh = st["Zhat"][:, cols]
    corr = np.zeros(st["Yc"].shape)
    corr[:, cols] = np.where(use, st["sd"][cols][None, :] * (np.nan_to_num(yhat) - Zh), 0.0)
    sup = np.zeros(st["Yc"].shape)
    sup[:, cols] = support
    return corr, sup


def apply(M, st, corr, gamma):
    Y = st["Yc"].copy()
    if gamma == 0:
        return Y
    miss = np.isnan(M) & st["ok"][None, :]
    Y[miss] = Y[miss] + gamma * corr[miss]
    return Y


def permute_rows(C, seed):
    rng = np.random.default_rng(seed + 4049)
    Cp = np.empty_like(C)
    for i in range(C.shape[0]):
        Cp[i] = C[i, rng.permutation(C.shape[1])]
    return Cp


# ---------------------------------------------------------------------------------------------------------------
# 6. diagnostics
def mass_diag(P, src, target):
    ds = sorted(set(src.tolist()))
    tot = P.sum()
    cm = P.sum(0)
    m = {d: float(cm[src == d].sum() / tot) for d in ds}
    s = {d: float(np.mean(src == d)) for d in ds}
    tv = 0.5 * sum(abs(m[d] - s[d]) for d in ds)
    top = sorted(ds, key=lambda d: -m[d] / s[d])[:3]
    me = META.get(target, {})
    out = {"tv": tv, "top3": [[d, round(m[d] / s[d], 3), round(m[d], 4)] for d in top],
           "mass": {d: round(v, 4) for d, v in m.items()}}
    for key, fn in (("acq", lambda x: me.get("acq") not in (None, "unknown") and x.get("acq") == me.get("acq")),
                    ("instrument", lambda x: me.get("instrument") not in (None, "unknown")
                     and x.get("instrument") == me.get("instrument")),
                    ("lineage", lambda x: bool(set(me.get("lineage", set())) & set(x.get("lineage", set()))))):
        sel = [d for d in ds if fn(META.get(d, {}))]
        out[f"{key}_mass"] = float(sum(m[d] for d in sel))
        out[f"{key}_share"] = float(sum(s[d] for d in sel))
    return out


def _log(variant, rec):
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / f"log_{variant}.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def _rs(x):
    return "inf" if np.isinf(x) else x


# ---------------------------------------------------------------------------------------------------------------
# 5. pipeline
def run(M, ctx, variant, null=False, rho_grid=RHO_GRID, base_only=False, log=True):
    seed = ctx["seed"]
    psi = select_psi(M, ctx)
    plans = [(e, r) for e in EPS_GRID for r in rho_grid]
    rec = {"target": ctx.get("target"), "seed": seed, "n_missing": int(np.isnan(M).sum()), "psi": _rs(psi)}
    if base_only:
        st = prepare(M, ctx, psi, 10) if log else None
        if log:
            _log(variant, {**rec, "chosen_gamma": 0.0})
        mu, sd, _ = cshrink_stats(M, psi, psi)
        Yc, _, _ = adapt_project_uv(M, ctx["V"], ctx["warm"], ctx["base"], mu, sd, lam_v=LAM_V, lam=LAM_U,
                                    extra=EXTRA, iters=10)
        del st
        return Yc

    # inner-mask selection
    Msel, (r1, c1), (r2, c2) = TS.selection_mask(M, seed)
    ss = prepare(Msel, ctx, psi, 5)
    Cs = permute_rows(ss["C"], seed) if null else ss["C"]
    scores, curves, pears = {}, {}, {}
    rr = np.concatenate([r1, r2])
    cc = np.concatenate([c1, c2])
    resid = M[rr, cc] - ss["Yc"][rr, cc]
    for e, r in plans:
        P, _ = sinkhorn(Cs, e * ss["c_med"], r * ss["c_med"])
        corr, _ = analog_corr(P, ss)
        cv = corr[rr, cc]
        pears[(e, r)] = float(np.corrcoef(cv, resid)[0, 1]) if np.std(cv) > 0 else None
        for gm in GAMMA_GRID:
            Y = ss["Yc"][rr, cc] + gm * cv
            e1 = np.sqrt(np.mean((Y[:len(r1)] - M[r1, c1]) ** 2))
            e2 = np.sqrt(np.mean((Y[len(r1):] - M[r2, c2]) ** 2))
            scores[(e, r, gm)] = float(0.5 * (e1 + e2))
        curves[(e, r)] = [scores[(e, r, gm)] for gm in GAMMA_GRID]
    ce, cr, cg = min(scores, key=lambda k: (scores[k], k[2], -k[0], -k[1]))

    # refit on M
    st = prepare(M, ctx, psi, 10)
    C = permute_rows(st["C"], seed) if null else st["C"]
    P, pinfo = sinkhorn(C, ce * st["c_med"], cr * st["c_med"])
    corr, sup = analog_corr(P, st)
    Y = apply(M, st, corr, cg)
    if log:
        miss = np.isnan(M) & st["ok"][None, :]
        fx = FIXED if FIXED[1] in rho_grid else (FIXED[0], rho_grid[0])
        Pf, finfo = sinkhorn(C, fx[0] * st["c_med"], fx[1] * st["c_med"])
        tv_all = {}
        for e, r in plans:
            Pp = Pf if (e, r) == fx else sinkhorn(C, e * st["c_med"], r * st["c_med"])[0]
            tv_all[f"{e}_{_rs(r)}"] = mass_diag(Pp, st["src"], ctx.get("target"))["tv"]
        rec.update({"chosen_eps": ce, "chosen_rho": _rs(cr), "chosen_gamma": cg,
                    "inner_best_score": scores[(ce, cr, cg)], "inner_base_score": scores[(ce, cr, 0.0)],
                    "gamma_curve": curves[(ce, cr)], "inner_pearson": pears[(ce, cr)],
                    "inner_pearson_all": {f"{e}_{_rs(r)}": v for (e, r), v in pears.items()},
                    "best_gamma_per_plan": {f"{e}_{_rs(r)}": GAMMA_GRID[int(np.argmin(curves[(e, r)]))]
                                            for (e, r) in plans},
                    "c_med": st["c_med"], "n_ref": st["n_ref"], "n_ref_dropped": st["n_ref_dropped"],
                    "sinkhorn": pinfo, "mean_support": float(sup[miss].mean()),
                    "frac_covered": float((sup[miss] >= MIN_SUPPORT).mean()),
                    "chosen_mass": mass_diag(P, st["src"], ctx.get("target")),
                    "fixed_plan": {"eps": fx[0], "rho": _rs(fx[1]), **finfo,
                                   **mass_diag(Pf, st["src"], ctx.get("target"))},
                    "tv_all_plans": tv_all})
        _log(variant, rec)
    return Y


def impute(Xm, ctx):
    return run(Xm, ctx, "main")


def impute_null(Xm, ctx):
    return run(Xm, ctx, "null", null=True)


def impute_base(Xm, ctx):
    return run(Xm, ctx, "base", base_only=True, log=False)


def impute_bal(Xm, ctx):
    return run(Xm, ctx, "bal", rho_grid=(float("inf"),))


# ---------------------------------------------------------------------------------------------------------------
# unit checks
def unit_check(target="petrosius2023"):
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H
    from scp_graphfs_mcp import phase2 as P_
    from scp_graphfs_mcp.phasec import cshrink as CS
    from scp_graphfs_mcp.phasec import nfbase as NF

    out = {"target": target}
    Xs, Xm, mask, names = NF._split(target, "MCAR", 100)
    V, warm = H.target_models(target, names)
    b = P_._fill(Xm, P_.loo_knn(Xm))
    b = np.where(np.isnan(b), np.nanmean(Xm), b)
    cache = {}

    def ref():
        if "Z" not in cache:
            cache["Z"], cache["src"] = H._ref_parts(target, names)
        return cache["Z"]

    def ref_src():
        ref()
        return cache["src"]

    ctx = {"V": V, "warm": warm, "base": b, "seed": 100, "target": target, "ref": ref, "ref_src": ref_src}
    # U2: gamma = 0 == cshrink
    yb = impute_base(Xm.copy(), ctx)
    yc = CS.run(Xm.copy(), ctx, "unit", log=False)
    out["U2_base_vs_cshrink_bitexact"] = bool(np.array_equal(yb, yc))
    out["U2_maxdiff"] = float(np.max(np.abs(yb - yc)))
    psi = select_psi(Xm, ctx)
    st = prepare(Xm, ctx, psi, 10)
    out["U2b_prepare_Yc_bitexact"] = bool(np.array_equal(st["Yc"], yc))
    # gamma = 0 through apply()
    corr, _ = analog_corr(sinkhorn(st["C"], 0.1 * st["c_med"], st["c_med"])[0], st)
    out["U2c_apply_gamma0_bitexact"] = bool(np.array_equal(apply(Xm, st, corr, 0.0), yc))
    # U1: balanced Sinkhorn marginals
    u1 = {}
    for e in EPS_GRID:
        P, info = sinkhorn(st["C"], e * st["c_med"], np.inf, max_it=5000, tol=1e-9)
        n, m = P.shape
        u1[str(e)] = {"iters": info["iters"], "row_err_L1": float(np.abs(P.sum(1) - 1 / n).sum()),
                      "col_err_L1": float(np.abs(P.sum(0) - 1 / m).sum()),
                      "row_err_max": float(np.abs(P.sum(1) - 1 / n).max()),
                      "col_err_max": float(np.abs(P.sum(0) - 1 / m).max())}
        P5, i5 = sinkhorn(st["C"], e * st["c_med"], np.inf)
        u1[str(e)]["default_500it"] = {"iters": i5["iters"], "row_err_L1": float(np.abs(P5.sum(1) - 1 / n).sum()),
                                       "col_err_L1": float(np.abs(P5.sum(0) - 1 / m).sum())}
    out["U1_balanced_marginals"] = u1
    out["U1_pass_1e-6"] = bool(all(v["row_err_L1"] < 1e-6 and v["col_err_L1"] < 1e-6 for v in u1.values()))
    # U3: TV of the null plan (balanced ~ 0, relaxed small) vs the real plan
    Cp = permute_rows(st["C"], 100)
    u3 = {}
    for e, r in [(0.1, np.inf), (0.1, 1.0), (0.1, 0.1), (0.05, 0.1)]:
        Pn = sinkhorn(Cp, e * st["c_med"], r * st["c_med"])[0]
        Pr = sinkhorn(st["C"], e * st["c_med"], r * st["c_med"])[0]
        u3[f"{e}_{_rs(r)}"] = {"null_tv": mass_diag(Pn, st["src"], target)["tv"],
                               "real_tv": mass_diag(Pr, st["src"], target)["tv"]}
    out["U3_tv"] = u3
    out["n_ref"], out["c_med"], out["psi"] = st["n_ref"], st["c_med"], _rs(psi)
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "unit_check.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    if "--unit" in sys.argv:
        unit_check()
