"""Phase-C idea `cshrink` (queen: innovation abduction). Champion revision, data-only.

Selected-psi James-Stein / Efron-Morris shrinkage of the target's per-protein 30-cell mean and log-sd toward CONSTANTS
(no trend, no permutation step), inside nfbase.adapt_project_mb. It is the explicit version of tgtshrink's null
pipeline (results/phaseC/tgtshrink_full.json null 0.9650), on a wider psi grid, with a one-global-sd arm.

  0. S = tgtshrink.target_stats(M); cal = base_ok & n >= 5 (cover = None, i.e. all base_ok proteins).
     cal.sum() < 10 -> champion-equivalent mu = xbar, sd = s.
  1. For (stat, var) in {(ybar, vy), (l, vl)}: t = mean(stat[cal]); r = stat - t;
     tau2 = max(1e-4, var(r[cal]) - psi mean(v[cal])); w = tau2 / (tau2 + psi v); post = t + w r.
     psi = inf -> w = 0 exactly; psi = 0 -> raw statistics (= the champion).
     mu = c_t + post_y, sd = exp(post_l) on base_ok; raw xbar / s elsewhere.
  2. Y = adapt_project_mb(M, V, warm, base, mu, sd, lam_v=100, extra=10, iters=10).
  3. psi selected on tgtshrink.selection_mask(M, seed) (50/50 MCAR/MNAR inner mask, rng seed + 31337), fits with
     iters = 5, score = 0.5 (RMSE_MCAR-part + RMSE_MNAR-part), argmin with ties to the smaller psi, refit on M.

Entry points: impute (joint psi on GRID), impute_null (psi uniform at random from GRID, rng seed + 9091 + n_missing),
impute_gsd (psi_l = inf, psi_y selected), impute_psi16, impute_psiinf, impute_oldgrid (joint psi on {1,2,4,8,16}).
Inner-mask scores are cached on disk keyed by sha1(Msel) + target + seed + (psi_y, psi_l) (exact; speed only).
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

from scp_graphfs_mcp.phasec import tgtshrink as TS
from scp_graphfs_mcp.phasec.nfbase import EXTRA, LAM_V, adapt_project_mb

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "cshrink"
INF = float("inf")
GRID = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, INF)
OLD_GRID = (1.0, 2.0, 4.0, 8.0, 16.0)


# ---------------------------------------------------------------------------------------------------------------
# 1. shrinkage toward constants
def cshrink_stats(M, psi_y, psi_l):
    S = TS.target_stats(M)
    base_ok, cal = S["base_ok"], S["cal"]
    mu, sd = S["xbar"].copy(), S["s"].copy()
    info = {"psi_y": psi_y, "psi_l": psi_l, "n_cal": int(cal.sum())}
    if cal.sum() < 10:
        info["fallback"] = "n_cal<10"
        return mu, sd, info
    for key, stat, var, psi in (("y", S["ybar"], S["vy"], psi_y), ("l", S["l"], S["vl"], psi_l)):
        if psi == 0:
            info.update({f"tau2_{key}": None, f"mean_w_{key}": 1.0})
            continue
        t = float(np.mean(stat[cal]))
        r = stat - t
        if np.isinf(psi):
            tau2 = 1e-4
            w = np.zeros_like(stat)
        else:
            tau2 = max(1e-4, float(np.var(r[cal]) - psi * np.mean(var[cal])))
            w = tau2 / (tau2 + psi * var)
        post = t + w * r
        info.update({f"tau2_{key}": tau2, f"mean_w_{key}": float(np.mean(w[base_ok]))})
        if key == "y":
            mu = np.where(base_ok & np.isfinite(post), S["c_t"] + post, mu)
        else:
            sd = np.where(base_ok & np.isfinite(post), np.exp(np.where(base_ok, post, 0.0)), sd)
    return mu, sd, info


# ---------------------------------------------------------------------------------------------------------------
# selection with an exact on-disk score cache
def _key(Msel, ctx, py, pl, tag):
    h = hashlib.sha1(np.ascontiguousarray(Msel).tobytes()).hexdigest()[:20]
    return f"{ctx.get('target')}_{ctx['seed']}_{h}_{tag}_{py}_{pl}"


def _score(M, Msel, parts, ctx, py, pl, shrink_fn, tag):
    d = CACHE / "scores"
    f = d / (_key(Msel, ctx, py, pl, tag) + ".json")
    if f.exists():
        try:
            return float(json.loads(f.read_text())["score"])
        except Exception:
            pass
    (r1, c1), (r2, c2) = parts
    mu, sd, _ = shrink_fn(Msel, py, pl)
    Y = adapt_project_mb(Msel, ctx["V"], ctx["warm"], ctx["base"], mu, sd, lam_v=LAM_V, extra=EXTRA, iters=5)
    e1 = np.sqrt(np.mean((Y[r1, c1] - M[r1, c1]) ** 2))
    e2 = np.sqrt(np.mean((Y[r2, c2] - M[r2, c2]) ** 2))
    s = float(0.5 * (e1 + e2))
    d.mkdir(parents=True, exist_ok=True)
    tmp = f.with_suffix(f".tmp{os.getpid()}")
    tmp.write_text(json.dumps({"score": s}))
    os.replace(tmp, f)
    return s


def _log(variant, rec):
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / f"log_{variant}.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def _ps(p):
    return "inf" if np.isinf(p) else str(p)


def run(M, ctx, variant, grid=GRID, fixed=None, gsd=False, random_psi=False, shrink_fn=None, tag="c", log=True):
    """fixed: joint psi (no selection). gsd: psi_l = inf, psi_y selected. random_psi: psi drawn from grid."""
    seed = ctx["seed"]
    if shrink_fn is None:
        shrink_fn = cshrink_stats
    scores = {}
    n_missing = int(np.isnan(M).sum())

    def pair(p):
        return (p, INF) if gsd else (p, p)

    if random_psi:
        rng = np.random.default_rng(seed + 9091 + n_missing)
        psi = grid[int(rng.integers(0, len(grid)))]
    elif fixed is not None:
        psi = fixed
    else:
        Msel, p1, p2 = TS.selection_mask(M, seed)
        for p in grid:
            scores[p] = _score(M, Msel, (p1, p2), ctx, *pair(p), shrink_fn, tag)
        psi = min(grid, key=lambda p: (scores[p], p))
    mu, sd, info = shrink_fn(M, *pair(psi))
    if log:
        rec = {"target": ctx.get("target"), "seed": seed, "n_missing": n_missing, "chosen_psi": _ps(psi),
               "scores": {_ps(k): v for k, v in scores.items()},
               **{k: (_ps(v) if isinstance(v, float) and np.isinf(v) else v) for k, v in info.items()
                  if not k.startswith("post_")}}
        _log(variant, rec)
    return adapt_project_mb(M, ctx["V"], ctx["warm"], ctx["base"], mu, sd, lam_v=LAM_V, extra=EXTRA, iters=10)


def impute(Xm, ctx):
    return run(Xm, ctx, "main")


def impute_null(Xm, ctx):
    return run(Xm, ctx, "null", random_psi=True)


def impute_gsd(Xm, ctx):
    return run(Xm, ctx, "gsd", gsd=True)


def impute_psi16(Xm, ctx):
    return run(Xm, ctx, "psi16", fixed=16.0)


def impute_psiinf(Xm, ctx):
    return run(Xm, ctx, "psiinf", fixed=INF)


def impute_oldgrid(Xm, ctx):
    return run(Xm, ctx, "oldgrid", grid=OLD_GRID)


# ---------------------------------------------------------------------------------------------------------------
# unit checks
def _tgt_null_shrink(seed):
    def f(M, py, pl):
        assert py == pl
        mu, sd, info = TS.shrink(M, py, "null", seed)
        return mu, sd, {k: v for k, v in info.items() if not k.startswith("post_")}
    return f


def unit_check():
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H  # noqa: F401
    from scp_graphfs_mcp import phase2 as P
    from scp_graphfs_mcp import transfer as T
    from scp_graphfs_mcp.phasec import nfbase as NF

    out = {"U1_codepath_vs_tgtshrink_null": {}, "U3_psi0_vs_adapt_project": {}}
    import tempfile

    global CACHE
    real = CACHE
    for t in ("petrosius2023", "leduc2022"):
        Xs, Xm, mask, names = NF._split(t, "MCAR", 100)
        V, warm = H.target_models(t, names)
        b = P._fill(Xm, P.loo_knn(Xm))
        b = np.where(np.isnan(b), np.nanmean(Xm), b)
        ctx = {"V": V, "warm": warm, "base": b, "seed": 100, "target": t}
        CACHE = Path(tempfile.mkdtemp())  # no score cache, no logs in the real dir
        TS.CACHE = CACHE / "ts"
        mine = run(Xm, ctx, "u1", grid=OLD_GRID, shrink_fn=_tgt_null_shrink(100), tag="u1")
        ref = TS.impute_null(Xm.copy(), ctx)
        miss = np.isnan(Xm)
        d = float(np.max(np.abs(mine[miss] - ref[miss])))
        out["U1_codepath_vs_tgtshrink_null"][t] = {"max_abs_diff": d, "pass": bool(d < 1e-9)}
        if t == "petrosius2023":
            mu, sd, _ = cshrink_stats(Xm, 0.0, 0.0)
            y = adapt_project_mb(Xm, V, warm, b, mu, sd, lam_v=LAM_V, extra=EXTRA, iters=10)
            c = T.adapt_project(Xm, V, warm, b, lam_v=LAM_V, extra=EXTRA)
            d3 = float(np.max(np.abs(y[miss] - c[miss])))
            out["U3_psi0_vs_adapt_project"][t] = {"max_abs_diff": d3, "pass": bool(d3 < 1e-9)}
        print(t, out, flush=True)
    CACHE = real
    TS.CACHE = ROOT / "data" / "cache" / "phasec" / "tgtshrink"
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "unit_check.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    if "--unit" in sys.argv:
        unit_check()
