"""v9 comparator panel: the H15 champion and its Phase-C revision against a broad set of standard imputation methods
and graph-prior methods (src/scp_graphfs_mcp/baselines.py), on the same splits.

    python benchmarks/panel_v9.py --pool fast|full [--methods a,b,...] [--jobs 2] [--tag dev]
        -> results/v9/dev_panel_<tag>_<pool>.json

DEV ONLY: targets come from benchmarks/phasec.py (leave-target-out references, holdouts asserted out). The holdout
run reuses `evaluate_split` from benchmarks/confirm_v9.py after PREREGISTRATION_v9.md is committed.

Per split (30 cells, MCAR or MNAR 10% mask, 10% inner mask for stacking), every method m is reported
  raw      : m alone                                   (what a user of m would get)
  stacked  : simplex stack of (kNN, low-rank, m)       (m given the same stacking help as the champion)
and the reference methods:
  knn, stack_lowrank (best data-only stack), champion = stack(kNN, low-rank, adapt)  [H15, confirmed on wu2026],
  revised = stack(kNN, low-rank, adapt, cshrink)  [Phase-C revision, dev only so far].
"""

from __future__ import annotations

import json
import os
import sys
import time
import warnings

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):  # one thread per worker process
    os.environ.setdefault(_v, "1")
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from explore_phaseb import GRID, load_target, subsample  # noqa: E402
from run import make_mask  # noqa: E402

from scp_graphfs_mcp import baselines as B  # noqa: E402
from scp_graphfs_mcp import phase2 as P  # noqa: E402
from scp_graphfs_mcp import transfer as T  # noqa: E402
from scp_graphfs_mcp.complexes import in_any_complex  # noqa: E402
from scp_graphfs_mcp.ids import group_accessions  # noqa: E402
from scp_graphfs_mcp.phasec import cshrink as CS  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
LAM_V, EXTRA = 100.0, 10


def evaluate_split(X, names, labels, V, warm, kind, seed, target, methods, n_cells=30, Vnull=None):
    """Returns (dict method -> squared errors on masked entries, in-complex mask, dict method -> stack weights,
    dict method -> seconds). With Vnull, also the champion's controls: champion@null (reference loadings permuted
    across proteins) and champion@noref (zero prior, same free factors)."""
    idx = subsample(X.shape[0], labels, n_cells, seed)
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

    def pk(M):
        return M[pick[:, 0], pick[:, 1]]

    def stack(fulls, inners):
        g = GRID[len(fulls)]
        errs = [np.nanmean((sum(w * pk(p) for w, p in zip(ws, inners)) - ti) ** 2) for ws in g]
        ws = g[int(np.nanargmin(errs))]
        return sum(w * p for w, p in zip(ws, fulls)), list(ws)

    parts, secs = {}, {}
    graph_cache = {}  # STRING graph / spectra depend only on names: share across eval and inner fits
    for key, M in (("eval", Xm), ("inner", Xi)):
        t0 = time.time()
        b = P._fill(M, P.loo_knn(M))
        b = np.where(np.isnan(b), np.nanmean(M), b)
        lr = P.grmf(M, None, names, lam_g=0.0)
        ad = T.adapt_project(M, V, warm, b, lam_v=LAM_V, extra=EXTRA)
        ctx = {"names": names, "target": target, "seed": seed, "base": b, "V": V, "warm": warm, **graph_cache}
        cs = CS.run(M, ctx, "v9panel", log=False)
        secs["_core"] = secs.get("_core", 0.0) + time.time() - t0
        row = {"base": b, "lr": lr, "ad": ad, "cs": cs}
        if Vnull is not None:
            row["ad_null"] = T.adapt_project(M, Vnull, warm, b, lam_v=LAM_V, extra=EXTRA)
            row["ad_noref"] = T.adapt_project(M, np.zeros_like(V), np.zeros_like(warm), b, lam_v=LAM_V, extra=EXTRA)
        for m in methods:
            t1 = time.time()
            Y = np.asarray(B.ALL[m](M.copy(), ctx), float)
            row[m] = np.where(np.isfinite(Y), Y, b)
            secs[m] = secs.get(m, 0.0) + time.time() - t1
        graph_cache.update({k: v for k, v in ctx.items() if k.startswith("_graph") or k.startswith("_nadj")
                            or k == "_spec_E"})
        parts[key] = row
    e, i = parts["eval"], parts["inner"]
    preds, wts = {"knn": e["base"], "lowrank": e["lr"]}, {}
    preds["stack_lowrank"], wts["stack_lowrank"] = stack([e["base"], e["lr"]], [i["base"], i["lr"]])
    preds["champion"], wts["champion"] = stack([e["base"], e["lr"], e["ad"]], [i["base"], i["lr"], i["ad"]])
    preds["revised"], wts["revised"] = stack([e["base"], e["lr"], e["ad"], e["cs"]],
                                             [i["base"], i["lr"], i["ad"], i["cs"]])
    if Vnull is not None:
        for k in ("null", "noref"):
            preds[f"champion@{k}"], wts[f"champion@{k}"] = stack([e["base"], e["lr"], e[f"ad_{k}"]],
                                                                 [i["base"], i["lr"], i[f"ad_{k}"]])
    for m in methods:
        preds[f"{m}"] = e[m]
        preds[f"{m}+stack"], wts[f"{m}+stack"] = stack([e["base"], e["lr"], e[m]], [i["base"], i["lr"], i[m]])
    se = {k: (Y[mask] - truth) ** 2 for k, Y in preds.items()}
    return se, sel, wts, secs


def _dev_split(args):
    from phasec import target_models

    tname, kind, seed, methods, tag = args
    cache = ROOT / "results" / "v9" / "dev_splits" / f"{tag}_{tname}_{kind}_{seed}.json"
    if cache.exists():  # resumable: a finished split is never recomputed
        d = json.loads(cache.read_text())
        if set(methods) <= set(d["methods"]):
            return d["rows"], d["timing"]
    X, names, labels = load_target(tname)
    V, warm = target_models(tname, names)
    se, sel, wts, secs = evaluate_split(X, names, labels, V, warm, kind, seed, tname, methods)
    rows = [{"target": tname, "mask": kind, "seed": seed, "method": m, "rmse_all": float(np.sqrt(v.mean())),
             "rmse_in_complex": float(np.sqrt(v[sel].mean())) if sel.any() else None, "weights": wts.get(m)}
            for m, v in se.items()]
    print(f"[{time.strftime('%H:%M:%S')}] {tname} {kind} {seed} " +
          " ".join(f"{k}={v:.0f}s" for k, v in secs.items()), flush=True)
    timing = {"target": tname, "mask": kind, "seed": seed, "seconds": secs}
    cache.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_suffix(f".tmp{os.getpid()}")
    tmp.write_text(json.dumps({"methods": methods, "rows": rows, "timing": timing}))
    os.replace(tmp, cache)
    return rows, timing


def summarize(rows):
    import pandas as pd

    df = pd.DataFrame(rows)
    out = {}
    for col in ("rmse_all", "rmse_in_complex"):
        m = df.groupby(["target", "mask", "method"])[col].mean().unstack("method")
        out[col] = {"mean_rel_to_champion": m.div(m["champion"], axis=0).mean().round(4).to_dict(),
                    "mean_rel_to_stack_lowrank": m.div(m["stack_lowrank"], axis=0).mean().round(4).to_dict(),
                    "n_target_masks_beating_champion_by_1pct":
                        (m.div(m["champion"], axis=0) <= 0.99).sum().astype(int).to_dict()}
    return out


def main():
    from phasec import FAST, FAST_SEEDS, FULL_SEEDS, dev_targets, target_models

    a = sys.argv
    pool = a[a.index("--pool") + 1] if "--pool" in a else "fast"
    tag = a[a.index("--tag") + 1] if "--tag" in a else "dev"
    jobs = int(a[a.index("--jobs") + 1]) if "--jobs" in a else 2
    methods = a[a.index("--methods") + 1].split(",") if "--methods" in a else list(B.ALL)
    targets = a[a.index("--targets") + 1].split(",") if "--targets" in a else (FAST if pool == "fast" else dev_targets())
    seeds = FAST_SEEDS if pool == "fast" else FULL_SEEDS
    assert not set(targets) & (T.PHASE_B_HOLDOUTS | T.HOLDOUTS | T.V9_HOLDOUTS), "holdouts are never dev targets"
    for t in targets:
        X, names, _ = load_target(t)
        target_models(t, names)
    t0 = time.time()
    jl = [(t, k, s, methods, tag) for t in targets for k in ("MCAR", "MNAR") for s in seeds]
    with ProcessPoolExecutor(jobs) as ex:
        res = list(ex.map(_dev_split, jl))
    rows = [r for rs, _ in res for r in rs]
    timing = [t for _, t in res]
    summ = summarize(rows)
    out_dir = ROOT / "results" / "v9"
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = {"pool": pool, "targets": targets, "seeds": seeds, "methods": methods, "runtime_s": round(time.time() - t0)}
    (out_dir / f"dev_panel_{tag}_{pool}.json").write_text(
        json.dumps({"meta": meta, "summary": summ, "timing": timing, "rows": rows}, indent=1) + "\n")
    rel = summ["rmse_all"]["mean_rel_to_champion"]
    for k in sorted(rel, key=rel.get):
        print(f"{k:22s} {rel[k]:.4f}  (vs stack_lowrank {summ['rmse_all']['mean_rel_to_stack_lowrank'][k]:.4f})")


if __name__ == "__main__":
    main()
