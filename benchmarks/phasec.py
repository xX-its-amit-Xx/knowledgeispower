"""Phase-C DEV harness: plug-in screening of new ideas against the H15 champion (docs/phaseC/README.md).

    python benchmarks/phasec.py --cand phasec.ideas.my_idea:impute [--null phasec.ideas.my_idea:impute_null]
                                --id my_idea [--pool fast|full] [--jobs 2]

Champion = kNN + low-rank + domain-adapted reference projection (lam_v=100, 10 extra), stacked on an inner mask; i.e.
the method confirmed on wu2026 (H15). A candidate is a function `impute(Xm, ctx) -> Y` returning a full cells x proteins
matrix (observed entries are ignored). It is called twice per split: on the evaluation input and on the inner-mask input,
so it must only use what it is given. ctx keys: names (protein ids), target, labels (subsample labels or None), seed,
base (kNN fill), lowrank, champion (the adapted projection alone), V / warm (reference loadings, leave-target-out),
ref() -> the z-scored single-cell reference aligned to `names` (lazy, cached), ref_src() -> dataset of each ref row.

Reported per candidate (all relative to the champion stack, < 1 is better):
  cand_stack = stack(kNN, low-rank, adapted, candidate)   -> does the idea ADD to the champion?
  cand_solo  = stack(kNN, low-rank, candidate)            -> could it REPLACE the adapted projection?
  null_stack = cand_stack with the --null function        -> is the gain from the idea's biological content?
Holdouts (transfer.PHASE_B_HOLDOUTS, transfer.HOLDOUTS) are never loaded or put in a reference. Dev seeds only.
"""

from __future__ import annotations

import importlib
import json
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from explore_phaseb import GRID, load_target, reference_datasets, subsample  # noqa: E402
from explore_phaseb import dev_targets as _phaseb_targets  # noqa: E402
from run import make_mask  # noqa: E402

from scp_graphfs_mcp import phase2 as P  # noqa: E402
from scp_graphfs_mcp import transfer as T  # noqa: E402
from scp_graphfs_mcp.complexes import in_any_complex  # noqa: E402
from scp_graphfs_mcp.ids import group_accessions  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache" / "phasec"
OUT = ROOT / "results" / "phaseC"
LAM_V, EXTRA = 100.0, 10
FAST = ["petrosius2023", "leduc2022", "weiss2026_liver", "gregoire2023_mixCTRL", "px063590_ADbrain", "sadiku2025"]
FAST_SEEDS = [100, 101]
FULL_SEEDS = [100, 101, 102]


def dev_targets():
    return [t for t in _phaseb_targets() if t not in T.NON_HUMAN]


def _ref_parts(tname, names):
    parts, src = [], []
    for d in reference_datasets(exclude=(tname,)):
        zd = T.corpus_reference(names, (tname,), [d])
        parts.append(zd)
        src += [d] * len(zd)
    return np.vstack(parts), np.array(src)


def target_models(tname, names):
    """Leave-target-out single-cell reference loadings, cached (same recipe as the v8 confirmation)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{tname}_V.npz"
    if f.exists():
        d = np.load(f)
        return d["V"], d["warm"]
    Z, _ = _ref_parts(tname, names)
    V, warm = T.fit_loadings(Z)
    np.savez_compressed(f, V=V, warm=warm)
    return V, warm


def _resolve(spec):
    mod, fn = spec.split(":")
    return getattr(importlib.import_module(mod), fn)


def run_split(args):
    tname, kind, seed, cand_spec, null_spec = args
    X, names, labels = load_target(tname)
    V, warm = target_models(tname, names)
    cand = _resolve(cand_spec)
    null = _resolve(null_spec) if null_spec else None
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

    def pk(M):
        return M[pick[:, 0], pick[:, 1]]

    def stack(fulls, inners):
        g = GRID[len(fulls)]
        errs = [np.nanmean((sum(w * pk(p) for w, p in zip(ws, inners)) - ti) ** 2) for ws in g]
        ws = g[int(np.nanargmin(errs))]
        return sum(w * p for w, p in zip(ws, fulls)), list(ws)

    cache = {}

    def ref():
        if "Z" not in cache:
            cache["Z"], cache["src"] = _ref_parts(tname, names)
        return cache["Z"]

    def ref_src():
        ref()
        return cache["src"]

    parts = {}
    for key, M in (("eval", Xm), ("inner", Xi)):
        b = P._fill(M, P.loo_knn(M))
        b = np.where(np.isnan(b), np.nanmean(M), b)
        lr = P.grmf(M, None, names, lam_g=0.0)
        ad = T.adapt_project(M, V, warm, b, lam_v=LAM_V, extra=EXTRA)
        ctx = {"names": names, "target": tname, "labels": None if labels is None else labels[idx], "seed": seed,
               "base": b, "lowrank": lr, "champion": ad, "V": V, "warm": warm, "ref": ref, "ref_src": ref_src}
        cp = np.asarray(cand(M.copy(), ctx), float)
        cp = np.where(np.isfinite(cp), cp, b)
        row = {"base": b, "lr": lr, "ad": ad, "cand": cp}
        if null:
            npred = np.asarray(null(M.copy(), ctx), float)
            row["null"] = np.where(np.isfinite(npred), npred, b)
        parts[key] = row
    e, i = parts["eval"], parts["inner"]
    preds, wts = {"knn": e["base"]}, {}
    preds["stack_lowrank"], wts["stack_lowrank"] = stack([e["base"], e["lr"]], [i["base"], i["lr"]])
    preds["champion"], wts["champion"] = stack([e["base"], e["lr"], e["ad"]], [i["base"], i["lr"], i["ad"]])
    preds["cand_raw"] = e["cand"]
    preds["cand_solo"], wts["cand_solo"] = stack([e["base"], e["lr"], e["cand"]], [i["base"], i["lr"], i["cand"]])
    preds["cand_stack"], wts["cand_stack"] = stack([e["base"], e["lr"], e["ad"], e["cand"]],
                                                   [i["base"], i["lr"], i["ad"], i["cand"]])
    if null:
        preds["null_stack"], wts["null_stack"] = stack([e["base"], e["lr"], e["ad"], e["null"]],
                                                       [i["base"], i["lr"], i["ad"], i["null"]])
    rows = []
    for m, Y in preds.items():
        err = (Y[mask] - truth) ** 2
        rows.append({"target": tname, "mask": kind, "seed": seed, "method": m, "rmse_all": float(np.sqrt(err.mean())),
                     "rmse_in_complex": float(np.sqrt(err[sel].mean())) if sel.any() else None,
                     "weights": wts.get(m)})
    return rows


def summarize(rows):
    import pandas as pd

    df = pd.DataFrame(rows)
    out = {}
    for col in ("rmse_all", "rmse_in_complex"):
        m = df.groupby(["target", "mask", "method"])[col].mean().unstack("method")
        rel = m.div(m["champion"], axis=0)
        out[col] = {"mean_rel_to_champion": rel.mean().round(4).to_dict(),
                    "per_target_mask": rel.round(4).reset_index().to_dict("records")}
        if "cand_stack" in rel:
            out[col]["n_splits_cand_stack_le_0.99"] = int((rel["cand_stack"] <= 0.99).sum())
            out[col]["n_splits"] = int(len(rel))
    w = df[df.method == "cand_stack"]["weights"].dropna().tolist()
    out["mean_cand_weight_in_stack"] = float(np.mean([x[3] for x in w])) if w else None
    return out


def main():
    a = sys.argv
    cand = a[a.index("--cand") + 1]
    null = a[a.index("--null") + 1] if "--null" in a else None
    iid = a[a.index("--id") + 1]
    pool = a[a.index("--pool") + 1] if "--pool" in a else "fast"
    jobs = int(a[a.index("--jobs") + 1]) if "--jobs" in a else 2
    targets = a[a.index("--targets") + 1].split(",") if "--targets" in a else (FAST if pool == "fast" else dev_targets())
    seeds = FAST_SEEDS if pool == "fast" else FULL_SEEDS
    assert not set(targets) & (T.PHASE_B_HOLDOUTS | T.HOLDOUTS), "holdouts are never Phase-C dev targets"
    for t in targets:  # build reference caches serially (memory)
        X, names, _ = load_target(t)
        target_models(t, names)
    t0 = time.time()
    jobs_list = [(t, k, s, cand, null) for t in targets for k in ("MCAR", "MNAR") for s in seeds]
    with ProcessPoolExecutor(jobs) as ex:
        rows = [r for rs in ex.map(run_split, jobs_list) for r in rs]
    summ = summarize(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    meta = {"id": iid, "cand": cand, "null": null, "pool": pool, "targets": targets, "seeds": seeds,
            "champion": {"lam_v": LAM_V, "extra": EXTRA}, "runtime_s": round(time.time() - t0, 1)}
    (OUT / f"{iid}_{pool}.json").write_text(json.dumps({"meta": meta, "summary": summ, "rows": rows}, indent=1) + "\n")
    print(json.dumps({"id": iid, "pool": pool, "rmse_all": summ["rmse_all"]["mean_rel_to_champion"],
                      "rmse_in_complex": summ["rmse_in_complex"]["mean_rel_to_champion"],
                      "splits_improved": summ["rmse_all"].get("n_splits_cand_stack_le_0.99"),
                      "n_splits": summ["rmse_all"].get("n_splits"), "runtime_s": meta["runtime_s"]}, indent=1))


if __name__ == "__main__":
    main()
