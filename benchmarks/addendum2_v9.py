"""v9 ADDENDUM 2 (post hoc, exploratory): cheap literature-standard imputers added after the frozen run, scored on the
same splits as the preregistered evaluation, against the champion.

    python benchmarks/addendum2_v9.py <holdout|--dev target> --methods mod:fn,mod:fn,...   [--jobs 1]
        -> results/v9/addendum2/<holdout or dev_target>.json

Per split it recomputes only kNN, low-rank, the adapted projection and the champion stack (the preregistered recipe), then each
listed method raw and stacked with kNN + low-rank. The `revised` method and the frozen panel stay in results/confirm_v9/ and
results/v9/addendum/. Per-split errors are cached under data/cache/v9_addendum2/ so restarts resume.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import os
import sys
import time
import warnings

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from explore_phaseb import GRID, reference_datasets, subsample  # noqa: E402
from run import make_mask  # noqa: E402

from scp_graphfs_mcp import phase2 as P  # noqa: E402
from scp_graphfs_mcp import transfer as T  # noqa: E402
from scp_graphfs_mcp.complexes import in_any_complex  # noqa: E402
from scp_graphfs_mcp.ids import group_accessions  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache" / "v9_addendum2"
OUT = ROOT / "results" / "v9" / "addendum2"


def _resolve(spec):
    mod, fn = spec.split(":")
    return getattr(importlib.import_module(mod), fn)


def split_errors(X, names, labels, V, warm, kind, seed, specs, n_cells=30):
    idx = subsample(X.shape[0], labels, n_cells, seed)
    Xs = X[idx]
    eligible = ((~np.isnan(Xs)).mean(0) >= 0.5) & ~np.isnan(Xs).all(0)
    mask = make_mask(Xs, eligible, kind, seed)
    Xm = Xs.copy()
    Xm[mask] = np.nan
    truth = Xs[mask]
    r, c = np.nonzero(mask)
    sel = np.array([in_any_complex(group_accessions(f)) for f in names])[c]
    rng = np.random.default_rng(seed + 7919)
    cand = np.argwhere(~np.isnan(Xm) & eligible[None, :])
    pick = cand[rng.choice(len(cand), max(1, int(0.1 * len(cand))), replace=False)]
    Xi = Xm.copy()
    Xi[pick[:, 0], pick[:, 1]] = np.nan
    ti = Xm[pick[:, 0], pick[:, 1]]

    def pk(M):
        return M[pick[:, 0], pick[:, 1]]

    def stack(fulls, inners):
        g = GRID[len(fulls)]
        errs = [np.nanmean((sum(w * pk(p) for w, p in zip(ws, inners)) - ti) ** 2) for ws in g]
        ws = g[int(np.nanargmin(errs))]
        return sum(w * p for w, p in zip(ws, fulls))

    parts, secs = {}, {}
    for key, M in (("eval", Xm), ("inner", Xi)):
        b = P._fill(M, P.loo_knn(M))
        b = np.where(np.isnan(b), np.nanmean(M), b)
        row = {"base": b, "lr": P.grmf(M, None, names, lam_g=0.0),
               "ad": T.adapt_project(M, V, warm, b, lam_v=100.0, extra=10)}
        ctx = {"names": names, "seed": seed, "base": b}
        for spec in specs:
            t0 = time.time()
            try:
                Y = np.asarray(_resolve(spec)(M.copy(), ctx), float)
            except Exception as e:  # a failing method is recorded, not fatal
                print(f"  {spec} failed: {e}", flush=True)
                Y = b
            row[spec] = np.where(np.isfinite(Y), Y, b)
            secs[spec] = secs.get(spec, 0.0) + time.time() - t0
        parts[key] = row
    e, i = parts["eval"], parts["inner"]
    preds = {"knn": e["base"], "stack_lowrank": stack([e["base"], e["lr"]], [i["base"], i["lr"]]),
             "champion": stack([e["base"], e["lr"], e["ad"]], [i["base"], i["lr"], i["ad"]])}
    for spec in specs:
        name = spec.split(":")[1]
        preds[name] = e[spec]
        preds[f"{name}+stack"] = stack([e["base"], e["lr"], e[spec]], [i["base"], i["lr"], i[spec]])
    se = {k: (Y[mask] - truth) ** 2 for k, Y in preds.items()}
    return se, sel, idx[r], c, secs


def main():
    a = sys.argv
    specs = a[a.index("--methods") + 1].split(",")
    if "--dev" in a:
        from explore_phaseb import load_target

        target = a[a.index("--dev") + 1]
        assert target not in (T.PHASE_B_HOLDOUTS | T.HOLDOUTS | T.V9_HOLDOUTS)
        X, names, labels = load_target(target)
        seeds, refs, excl, tag = [100, 101], reference_datasets(exclude=(target,)), (target,), f"dev_{target}"
    else:
        import confirm_v9 as C

        target = a[1]
        assert target in C.PRIMARY + C.SECONDARY
        if not (ROOT / "results" / "confirm_v9" / f"{target}.json").exists():
            sys.exit(f"refusing: preregistered result for {target} not produced yet")
        X, names = C.load_holdout(target)
        labels, seeds, refs, excl, tag = None, C.SEEDS, reference_datasets(), (), target
    if "--seeds" in a:  # a subset of the preregistered seeds (runtime); recorded in the output
        seeds = [int(x) for x in a[a.index("--seeds") + 1].split(",")]
    assert not set(refs) & (T.PHASE_B_HOLDOUTS | T.HOLDOUTS | T.V9_HOLDOUTS)
    V, warm = T.fit_loadings(T.corpus_reference(names, excl, refs))
    out = {"dataset": target, "seeds": seeds, "post_hoc": True, "methods": specs, "masks": {}}
    for kind in ("MCAR", "MNAR"):
        SE, SEL, SPL = {}, [], []
        for seed in seeds:
            f = CACHE / f"{tag}_{kind}_{seed}_{hashlib.sha1(','.join(specs).encode()).hexdigest()[:10]}.npz"
            if f.exists():
                d = np.load(f, allow_pickle=True)
                se, sel = d["se"].item(), d["sel"]
            else:
                t0 = time.time()
                se, sel, cell, prot, secs = split_errors(X, names, labels, V, warm, kind, seed, specs)
                CACHE.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(f, se=np.array(se, dtype=object), sel=sel, cell=cell, prot=prot)
                print(f"[{time.strftime('%H:%M:%S')}] {target} {kind} {seed} ({time.time() - t0:.0f}s) "
                      + " ".join(f"{k.split(':')[1]}={v:.0f}s" for k, v in secs.items()), flush=True)
            for m, v in se.items():
                SE.setdefault(m, []).append(v)
            SEL.append(sel)
            SPL.append(np.full(len(sel), seed))
        SE = {m: np.concatenate(v) for m, v in SE.items()}
        SEL, SPL = np.concatenate(SEL), np.concatenate(SPL)
        rmse = {m: float(np.sqrt(v.mean())) for m, v in SE.items()}
        per_split = {m: [float(np.sqrt(SE[m][SPL == s].mean()) / np.sqrt(SE["champion"][SPL == s].mean()))
                         for s in seeds] for m in SE}
        out["masks"][kind] = {"rmse": rmse, "rel_to_champion": {m: v / rmse["champion"] for m, v in rmse.items()},
                              "rmse_in_complex": {m: float(np.sqrt(v[SEL].mean())) for m, v in SE.items()},
                              "per_split_rel_to_champion": per_split}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{tag}.json").write_text(json.dumps(out, indent=1) + "\n")
    for kind, b in out["masks"].items():
        print(kind, {m: round(v, 3) for m, v in sorted(b["rel_to_champion"].items(), key=lambda kv: kv[1])})


if __name__ == "__main__":
    main()
