"""Phase-D DEV screen: does stacking the H15 champion's components and the Phase-C revision with the strongest
data-only comparators found in v9 (bpca_oba, gcn_nograph) help?

    nice -n 15 env OMP_NUM_THREADS=1 python benchmarks/phased_dev.py            # compute missing splits, summarise
    python benchmarks/phased_dev.py --summarize                                 # summarise cached splits only

DEV ONLY. Targets = benchmarks/phasec.FAST, seeds = phasec.FAST_SEEDS (100, 101), masks MCAR and MNAR: 24 splits.
Holdouts (transfer.PHASE_B_HOLDOUTS | HOLDOUTS | V9_HOLDOUTS) are asserted out and never loaded.

Per split the recipe is exactly benchmarks/panel_v9.evaluate_split: 30 cells via explore_phaseb.subsample (with labels),
run.make_mask 10% MCAR/MNAR evaluation mask over eligible proteins (>= 50% observed), 10% inner mask of the remaining
observed eligible entries (rng seed + 7919) used only to choose stack weights. Components, each computed on both the
evaluation input and the inner input:
  knn   phase2 loo_knn fill (NaN -> global mean)             lr    phase2.grmf(lam_g=0)
  ad    transfer.adapt_project(lam_v=100, extra=10) with V / warm from phasec.target_models (leave-target-out)
  cs    phasec.cshrink.run(M, ctx, 'phaseD', log=False), psi-score cache redirected to data/cache/phaseD/cshrink
  bpca  baselines_addendum.bpca_oba                          gcn   baselines.gcn_nograph
Stacks: simplex weights on a 0.25 grid (explore_phaseb.GRID convention, itertools.product order, first argmin wins),
extended generically to 5 and 6 components (the 6-component stack also uses the 0.25 grid; an unconstrained
non-negative least-squares stack of the 6 components on the inner mask is reported as an extra row, all6_nnls).
Only the scored values (evaluation-mask entries, inner-mask entries) of each component are cached, one file per
(split, component) under data/cache/phaseD/parts/, so the run is resumable at component granularity.
Output: results/phaseD/dev_combos.json.
"""

from __future__ import annotations

import itertools
import json
import os
import sys
import time
import warnings

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from explore_phaseb import G2, GRID, load_target, subsample  # noqa: E402
from phasec import FAST, FAST_SEEDS, target_models  # noqa: E402
from run import make_mask  # noqa: E402

from scp_graphfs_mcp import baselines as B  # noqa: E402
from scp_graphfs_mcp import baselines_addendum as BA  # noqa: E402
from scp_graphfs_mcp import phase2 as P  # noqa: E402
from scp_graphfs_mcp import transfer as T  # noqa: E402
from scp_graphfs_mcp.complexes import in_any_complex  # noqa: E402
from scp_graphfs_mcp.ids import group_accessions  # noqa: E402
from scp_graphfs_mcp.phasec import cshrink as CS  # noqa: E402
from scp_graphfs_mcp.phasec import tgtshrink as TS  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache" / "phaseD"
PARTS = CACHE / "parts"
OUT = ROOT / "results" / "phaseD"
LAM_V, EXTRA = 100.0, 10
CS.CACHE = CACHE / "cshrink"  # never write into the frozen data/cache/phasec/cshrink
TS.CACHE = CACHE / "tgtshrink"  # (cshrink.run(log=False) does not write here; redirected for safety)
COMPONENTS = ("knn", "lr", "ad", "cs", "bpca", "gcn")
COMBOS = {  # name -> components, in stacking order
    "champion": ("knn", "lr", "ad"),
    "revised": ("knn", "lr", "ad", "cs"),
    "bpca_oba+stack": ("knn", "lr", "bpca"),
    "gcn_nograph+stack": ("knn", "lr", "gcn"),
    "knn_lr_ad_bpca": ("knn", "lr", "ad", "bpca"),
    "knn_lr_ad_gcn": ("knn", "lr", "ad", "gcn"),
    "knn_lr_cs_bpca": ("knn", "lr", "cs", "bpca"),
    "knn_lr_ad_cs_bpca": ("knn", "lr", "ad", "cs", "bpca"),
    "knn_lr_ad_cs_bpca_gcn": ("knn", "lr", "ad", "cs", "bpca", "gcn"),
}
NNLS6 = "knn_lr_ad_cs_bpca_gcn__nnls"
RAW = {"knn": "knn", "lowrank": "lr", "adapt": "ad", "cshrink": "cs", "bpca_oba": "bpca", "gcn_nograph": "gcn"}


def simplex_grid(n):
    if n in GRID:  # identical object (and order) to the frozen harnesses for n = 2..4
        return GRID[n]
    return [w for w in itertools.product(G2, repeat=n) if abs(sum(w) - 1) < 1e-9]


def split_setup(X, labels, kind, seed, names):
    """Exactly panel_v9.evaluate_split's split recipe."""
    idx = subsample(X.shape[0], labels, 30, seed)
    Xs = X[idx]
    eligible = ((~np.isnan(Xs)).mean(0) >= 0.5) & ~np.isnan(Xs).all(0)
    mask = make_mask(Xs, eligible, kind, seed)
    Xm = Xs.copy()
    Xm[mask] = np.nan
    truth = Xs[mask]
    sel = np.array([in_any_complex(group_accessions(f)) for f in names])[np.nonzero(mask)[1]]
    rng = np.random.default_rng(seed + 7919)
    c = np.argwhere(~np.isnan(Xm) & eligible[None, :])
    pick = c[rng.choice(len(c), max(1, int(0.1 * len(c))), replace=False)]
    Xi = Xm.copy()
    Xi[pick[:, 0], pick[:, 1]] = np.nan
    ti = Xm[pick[:, 0], pick[:, 1]]
    return Xm, Xi, mask, truth, sel, pick, ti


def compute_split(tname, kind, seed):
    tag = f"{tname}_{kind}_{seed}"
    files = {c: PARTS / f"{tag}_{c}.npz" for c in COMPONENTS}
    meta_f = PARTS / f"{tag}_meta.npz"
    if meta_f.exists() and all(f.exists() for f in files.values()):
        return
    X, names, labels = load_target(tname)
    V, warm = target_models(tname, names)
    Xm, Xi, mask, truth, sel, pick, ti = split_setup(X, labels, kind, seed, names)
    PARTS.mkdir(parents=True, exist_ok=True)
    if not meta_f.exists():
        tmp = PARTS / f".{tag}_meta.{os.getpid()}.npz"
        np.savez_compressed(tmp, truth=truth, sel=sel, ti=ti)
        os.replace(tmp, meta_f)
    memo = {}

    def base(key, M):
        if key not in memo:
            b = P._fill(M, P.loo_knn(M))
            memo[key] = np.where(np.isnan(b), np.nanmean(M), b)
        return memo[key]

    def ctx(key, M):
        return {"names": names, "target": tname, "seed": seed, "base": base(key, M), "V": V, "warm": warm}

    def fn(comp, key, M):
        b = base(key, M)
        if comp == "knn":
            return b
        if comp == "lr":
            return P.grmf(M, None, names, lam_g=0.0)
        if comp == "ad":
            return T.adapt_project(M, V, warm, b, lam_v=LAM_V, extra=EXTRA)
        if comp == "cs":
            return CS.run(M, ctx(key, M), "phaseD", log=False)
        if comp == "bpca":
            Y = np.asarray(BA.bpca_oba(M.copy(), ctx(key, M)), float)
        elif comp == "gcn":
            Y = np.asarray(B.gcn_nograph(M.copy(), ctx(key, M)), float)
        return np.where(np.isfinite(Y), Y, b)  # same guard as panel_v9 for comparator methods

    for comp in COMPONENTS:
        if files[comp].exists():
            continue
        t0 = time.time()
        ye = fn(comp, "eval", Xm)[mask]
        yi = fn(comp, "inner", Xi)[pick[:, 0], pick[:, 1]]
        tmp = PARTS / f".{tag}_{comp}.{os.getpid()}.npz"
        np.savez_compressed(tmp, eval=ye, inner=yi, seconds=time.time() - t0)
        os.replace(tmp, files[comp])
        print(f"[{time.strftime('%H:%M:%S')}] {tag} {comp} {time.time() - t0:.0f}s", flush=True)


def load_split(tname, kind, seed):
    tag = f"{tname}_{kind}_{seed}"
    m = np.load(PARTS / f"{tag}_meta.npz")
    comps = {}
    for c in COMPONENTS:
        d = np.load(PARTS / f"{tag}_{c}.npz")
        comps[c] = (d["eval"], d["inner"], float(d["seconds"]))
    return m["truth"], m["sel"].astype(bool), m["ti"], comps


def stack(comps, names, ti):
    g = simplex_grid(len(names))
    inn = [comps[c][1] for c in names]
    errs = [np.nanmean((sum(w * p for w, p in zip(ws, inn)) - ti) ** 2) for ws in g]
    ws = g[int(np.nanargmin(errs))]
    return sum(w * comps[c][0] for w, c in zip(ws, names)), [float(w) for w in ws]


def stack_nnls(comps, names, ti):
    from scipy.optimize import nnls

    A = np.column_stack([comps[c][1] for c in names])
    ok = np.isfinite(A).all(1) & np.isfinite(ti)
    w, _ = nnls(A[ok], ti[ok])
    return sum(wi * comps[c][0] for wi, c in zip(w, names)), [float(x) for x in w]


def split_rows(tname, kind, seed):
    truth, sel, ti, comps = load_split(tname, kind, seed)
    preds, wts = {}, {}
    for raw, c in RAW.items():
        preds[f"raw:{raw}"] = comps[c][0]
    for name, cs in COMBOS.items():
        preds[name], wts[name] = stack(comps, cs, ti)
    preds[NNLS6], wts[NNLS6] = stack_nnls(comps, COMBOS["knn_lr_ad_cs_bpca_gcn"], ti)
    rows = []
    for m, y in preds.items():
        err = (y - truth) ** 2
        rows.append({"target": tname, "mask": kind, "seed": seed, "method": m,
                     "rmse_all": float(np.sqrt(np.mean(err))),
                     "rmse_in_complex": float(np.sqrt(np.mean(err[sel]))) if sel.any() else None,
                     "weights": wts.get(m)})
    secs = {c: comps[c][2] for c in COMPONENTS}
    return rows, secs


def _t95(x):
    """Mean of per-split log RMSE ratios with a t 95% interval, back-transformed to ratios (descriptive: the 24
    splits share 6 targets, so this is optimistic; dev screen only)."""
    from scipy import stats

    x = x[np.isfinite(x)]
    if len(x) < 2 or np.allclose(x, 0):
        return {"geo_mean_ratio": round(float(np.exp(x.mean())), 4) if len(x) else None, "ci95": None}
    h = stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return {"geo_mean_ratio": round(float(np.exp(x.mean())), 4),
            "ci95": [round(float(np.exp(x.mean() - h)), 4), round(float(np.exp(x.mean() + h)), 4)]}


def summarize(rows):
    import pandas as pd

    df = pd.DataFrame(rows)
    out = {}
    for col in ("rmse_all", "rmse_in_complex"):
        piv = df.pivot_table(index=["target", "mask", "seed"], columns="method", values=col)
        rel_split = piv.div(piv["champion"], axis=0)
        tm = df.groupby(["target", "mask", "method"])[col].mean().unstack("method")
        rel_tm = tm.div(tm["champion"], axis=0)
        out[col] = {
            "mean_rel_to_champion_over_target_masks": rel_tm.mean().round(4).to_dict(),
            "mean_rel_to_champion_over_splits": rel_split.mean().round(4).to_dict(),
            "median_rel_to_champion_over_splits": rel_split.median().round(4).to_dict(),
            "n_splits_le_0.99": (rel_split <= 0.99).sum().astype(int).to_dict(),
            "n_splits_lt_1": (rel_split < 1 - 1e-12).sum().astype(int).to_dict(),
            "n_splits_gt_1": (rel_split > 1 + 1e-12).sum().astype(int).to_dict(),
            "n_splits_ge_1.01": (rel_split >= 1.01).sum().astype(int).to_dict(),
            "n_splits": int(len(rel_split)),
            "per_target_mask_rel_to_champion": rel_tm.round(4).reset_index().to_dict("records"),
        }
        rv_split = piv.div(piv["revised"], axis=0)  # does anything add to the Phase-C revision?
        rv_tm = tm.div(tm["revised"], axis=0)
        out[col]["vs_revised"] = {
            "mean_rel_to_revised_over_target_masks": rv_tm.mean().round(4).to_dict(),
            "mean_rel_to_revised_over_splits": rv_split.mean().round(4).to_dict(),
            "n_splits_le_0.99": (rv_split <= 0.99).sum().astype(int).to_dict(),
            "n_splits_lt_1": (rv_split < 1 - 1e-12).sum().astype(int).to_dict(),
            "n_splits_gt_1": (rv_split > 1 + 1e-12).sum().astype(int).to_dict(),
            "n_splits_ge_1.01": (rv_split >= 1.01).sum().astype(int).to_dict(),
            "split_t95_of_log_ratio": {m: _t95(np.log(rv_split[m].values)) for m in rv_split.columns},
        }
        out[col]["split_t95_of_log_ratio_vs_champion"] = {m: _t95(np.log(rel_split[m].values))
                                                          for m in rel_split.columns}
    wt = {}
    for name, cs in list(COMBOS.items()) + [(NNLS6, COMBOS["knn_lr_ad_cs_bpca_gcn"])]:
        W = np.array(df[df.method == name]["weights"].tolist())
        d = {c: round(float(W[:, i].mean()), 4) for i, c in enumerate(cs)}
        for c in ("bpca", "gcn"):
            if c in cs:
                i = cs.index(c)
                d[f"n_splits_{c}_weight_gt0"] = int((W[:, i] > 0).sum())
        wt[name] = d
    out["mean_stack_weights"] = wt
    return out


def splits():
    return [(t, k, s) for t in FAST for k in ("MCAR", "MNAR") for s in FAST_SEEDS]


def main():
    assert not set(FAST) & (T.PHASE_B_HOLDOUTS | T.HOLDOUTS | T.V9_HOLDOUTS), "holdouts are never Phase-D dev targets"
    t0 = time.time()
    if "--summarize" not in sys.argv:
        for t, k, s in splits():
            compute_split(t, k, s)
    done = [sp for sp in splits() if (PARTS / f"{sp[0]}_{sp[1]}_{sp[2]}_meta.npz").exists()
            and all((PARTS / f"{sp[0]}_{sp[1]}_{sp[2]}_{c}.npz").exists() for c in COMPONENTS)]
    rows, timing = [], []
    for sp in done:
        r, secs = split_rows(*sp)
        rows += r
        timing.append({"target": sp[0], "mask": sp[1], "seed": sp[2], "seconds_eval_plus_inner": secs})
    summ = summarize(rows)
    meta = {"phase": "D", "dev_only": True, "holdouts_touched": False, "targets": FAST, "seeds": FAST_SEEDS,
            "masks": ["MCAR", "MNAR"], "n_cells": 30, "n_splits_done": len(done), "n_splits_planned": len(splits()),
            "champion": {"lam_v": LAM_V, "extra": EXTRA},
            "components": {"knn": "phase2 loo_knn fill", "lr": "phase2.grmf lam_g=0",
                           "ad": "transfer.adapt_project lam_v=100 extra=10, V/warm phasec.target_models",
                           "cs": "phasec.cshrink.run(M, ctx, 'phaseD', log=False), cache data/cache/phaseD/cshrink",
                           "bpca": "baselines_addendum.bpca_oba", "gcn": "baselines.gcn_nograph"},
            "combos": {k: list(v) for k, v in COMBOS.items()},
            "stacking": "simplex grid step 0.25 chosen on the 10% inner mask for every combo (2-6 components; "
                        "first argmin in itertools.product order); extra row " + NNLS6 + " = unconstrained "
                        "non-negative least squares of the 6 components on the inner mask",
            "relative_metric": "RMSE / RMSE(champion); headline = mean over the 12 target x mask cells of the "
                               "seed-averaged RMSE ratio (phasec.summarize convention); split-level also given",
            "runtime_s_this_call": round(time.time() - t0, 1)}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "dev_combos.json").write_text(json.dumps({"meta": meta, "summary": summ, "timing": timing,
                                                     "rows": rows}, indent=1) + "\n")
    print(json.dumps({"n_splits": len(done), "rmse_all": summ["rmse_all"]["mean_rel_to_champion_over_target_masks"],
                      "weights": summ["mean_stack_weights"]}, indent=1))


if __name__ == "__main__":
    main()
