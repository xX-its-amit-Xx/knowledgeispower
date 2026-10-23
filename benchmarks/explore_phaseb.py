"""Phase-B DEV harness (Amendment 15): leave-one-dataset-out over many single-cell targets.

    python benchmarks/explore_phaseb.py --tag <name> [--bulk ccle,proteomehd,...] [--methods ...]

For each dev target, the reference is every other non-holdout corpus dataset plus the chosen bulk atlases. Methods:
  knn, stack_lowrank (best data-only), stack3_plain (Phase-A best: kNN + low-rank + linear reference projection),
  mae / stack4_mae (kNN + low-rank + projection + masked autoencoder), stack4_mae@null (MAE trained on a
  protein-scrambled reference). 30-cell subsamples (stratified when labels exist), MCAR + MNAR, dev seeds 100-102.
Holdouts (transfer.PHASE_B_HOLDOUTS, khan2023, derks2022) are never loaded or put in a reference.
"""

from __future__ import annotations

import itertools
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run import make_mask  # noqa: E402

from scp_graphfs_mcp import mae as MAE  # noqa: E402
from scp_graphfs_mcp import phase2 as P  # noqa: E402
from scp_graphfs_mcp import transfer as T  # noqa: E402
from scp_graphfs_mcp.complexes import in_any_complex  # noqa: E402
from scp_graphfs_mcp.datasets import load_leduc2022, load_petrosius2023, load_specht2021  # noqa: E402
from scp_graphfs_mcp.external import CORPUS_DIR  # noqa: E402
from scp_graphfs_mcp.ids import group_accessions  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
SEEDS = [100, 101, 102]
G2 = (0.0, 0.25, 0.5, 0.75, 1.0)
GRID = {n: [w for w in itertools.product(G2, repeat=n) if abs(sum(w) - 1) < 1e-9] for n in (2, 3, 4)}
ANNDATA_TARGETS = {"petrosius2023": load_petrosius2023, "leduc2022": load_leduc2022, "specht2021": load_specht2021}
MIN_CELLS = 100


def reference_datasets(exclude=()):
    ds = [d for d in T.extended_corpus() if d not in T.PHASE_B_HOLDOUTS and d not in exclude]
    assert not set(ds) & (T.PHASE_B_HOLDOUTS | T.HOLDOUTS)
    return ds


def dev_targets():
    out = list(ANNDATA_TARGETS)
    for d in reference_datasets():
        if d in out:
            continue
        X = np.load(CORPUS_DIR / f"{d}.npz")["X"]
        if X.shape[0] >= MIN_CELLS:
            out.append(d)
    return out


def load_target(name):
    """(X cells x proteins log2, var_names, labels or None)."""
    if name in ANNDATA_TARGETS:
        a = ANNDATA_TARGETS[name]()
        return np.asarray(a.X, float), list(a.var_names), a.obs["cell_type"].astype(str).values
    d = np.load(CORPUS_DIR / f"{name}.npz")
    lab = CORPUS_DIR / f"{name}_labels.json"
    labels = np.array(json.loads(lab.read_text())) if lab.exists() else None
    return d["X"].astype(float), [str(a) for a in d["accs"]], labels


def subsample(n_total, labels, n_cells, seed):
    rng = np.random.default_rng(seed)
    if labels is not None:
        types = [t for t in sorted(set(labels)) if (labels == t).sum() >= n_cells // len(set(labels))]
        if len(types) >= 2:
            k = n_cells // len(types)
            return np.sort(np.concatenate([rng.choice(np.where(labels == t)[0], k, replace=False) for t in types]))
    return np.sort(rng.choice(n_total, n_cells, replace=False))


def evaluate(X, names, models, n_cells, kind, seed, labels, ceiling=False):
    idx = subsample(X.shape[0], labels, n_cells, seed)
    Xs = X[idx]
    keepcol = ~np.isnan(Xs).all(0)
    eligible = ((~np.isnan(Xs)).mean(0) >= 0.5) & keepcol
    mask = make_mask(Xs, eligible, kind, seed)
    Xm = Xs.copy()
    Xm[mask] = np.nan
    truth = Xs[mask]
    in_cx = np.array([in_any_complex(group_accessions(f)) for f in names])
    sel = in_cx[np.nonzero(mask)[1]]
    rng = np.random.default_rng(seed + 7919)
    cand = np.argwhere(~np.isnan(Xm) & eligible[None, :])
    pick = cand[rng.choice(len(cand), max(1, int(0.1 * len(cand))), replace=False)]
    Xi = Xm.copy()
    Xi[pick[:, 0], pick[:, 1]] = np.nan
    ti = Xm[pick[:, 0], pick[:, 1]]
    pk = lambda M: M[pick[:, 0], pick[:, 1]]  # noqa: E731

    def stack(fulls, inners):
        g = GRID[len(fulls)]
        errs = [np.nanmean((sum(w * pk(p) for w, p in zip(ws, inners)) - ti) ** 2) for ws in g]
        ws = g[int(np.nanargmin(errs))]
        return sum(w * p for w, p in zip(ws, fulls)), list(ws)

    base, bi = P._fill(Xm, P.loo_knn(Xm)), P._fill(Xi, P.loo_knn(Xi))
    base, bi = np.where(np.isnan(base), np.nanmean(Xm), base), np.where(np.isnan(bi), np.nanmean(Xi), bi)
    lr, lri = P.grmf(Xm, None, names, lam_g=0.0), P.grmf(Xi, None, names, lam_g=0.0)
    V, warm = models["V"], models["warm"]
    pj, pji = T.project(Xm, V, warm, base), T.project(Xi, V, warm, bi)
    preds, wts = {"knn": base}, {}
    preds["stack_lowrank"], wts["stack_lowrank"] = stack([base, lr], [bi, lri])
    preds["stack3_plain"], wts["stack3_plain"] = stack([base, lr, pj], [bi, lri, pji])
    if "Zref" in models:  # Amendment 15 workstream 3: reference selection (+ random-selection control of equal size)
        Zr, src = models["Zref"], models["src"]
        for lab, kw in (("sel_cells2000", {"n": 2000}), ("sel_rand2000", {"n": 2000, "seed": seed}),
                        ("sel_ds3", {"n": 0, "mode": "datasets", "src": src})):
            rows_sel = T.select_reference(Zr, Xm, **kw)
            Vs, ws = T.fit_loadings(Zr[rows_sel])
            ps, psi = T.project(Xm, Vs, ws, base), T.project(Xi, Vs, ws, bi)
            preds[f"stack3_{lab}"], wts[f"stack3_{lab}"] = stack([base, lr, ps], [bi, lri, psi])
    if ceiling:  # diagnostic, not a candidate: the ideal reference = the *other* cells of the same dataset
        rest = np.setdiff1d(np.arange(X.shape[0]), idx)
        Vin, win = T.fit_loadings(T._z(X[rest]))
        pin, pini = T.project(Xm, Vin, win, base), T.project(Xi, Vin, win, bi)
        preds["ceiling_indomain"], wts["ceiling_indomain"] = stack([base, lr, pin], [bi, lri, pini])
    for key in ("mae", "mae@null"):
        if key not in models:
            continue
        m, mi = MAE.mae_impute(Xm, models[key], base), MAE.mae_impute(Xi, models[key], bi)
        if key == "mae":
            preds["mae"] = m
        preds[f"stack4_{key}"], wts[f"stack4_{key}"] = stack([base, lr, pj, m], [bi, lri, pji, mi])
    rows = []
    for meth, Y in preds.items():
        e = (Y[mask] - truth) ** 2
        rows.append({"method": meth, "rmse_all": float(np.sqrt(e.mean())),
                     "rmse_in_complex": float(np.sqrt(e[sel].mean())) if sel.any() else None,
                     "weights": wts.get(meth), "n_cells": int(len(idx)), "mask": kind, "seed": seed})
    return rows


def main():
    tag = sys.argv[sys.argv.index("--tag") + 1] if "--tag" in sys.argv else "phaseb"
    bulk = tuple(b for b in sys.argv[sys.argv.index("--bulk") + 1].split(",") if b) if "--bulk" in sys.argv else T.BULK_MODULES
    no_sc = "--no-sc-ref" in sys.argv  # ablation: bulk-only reference
    only = sys.argv[sys.argv.index("--targets") + 1].split(",") if "--targets" in sys.argv else None
    epochs = int(sys.argv[sys.argv.index("--epochs") + 1]) if "--epochs" in sys.argv else 30
    out_dir = ROOT / "results" / "phaseB" / "dev"
    out_dir.mkdir(parents=True, exist_ok=True)
    targets = only or dev_targets()
    meta = {"tag": tag, "bulk": bulk, "targets": targets, "seeds": SEEDS, "epochs": epochs}
    rows, t0 = [], time.time()
    for tname in targets:
        X, names, labels = load_target(tname)
        refs = reference_datasets(exclude=(tname,))
        parts, src = [], []
        for d in ([] if no_sc else refs):
            if d == tname:
                continue
            zd = T.corpus_reference(names, (tname,), [d])
            parts.append(zd)
            src += [d] * len(zd)
        zb = T.bulk_reference(names, bulk)
        parts.append(zb)
        src += ["bulk"] * len(zb)
        Z, src = np.vstack(parts), np.array(src)
        V, warm = T.fit_loadings(Z)
        models = {"V": V, "warm": warm}
        if "--select" in sys.argv:
            models["Zref"], models["src"] = Z, src
        if "--no-mae" not in sys.argv:
            models["mae"] = MAE.train_mae(Z, epochs=epochs, seed=0)
            Zs = Z[:, np.random.default_rng(0).permutation(Z.shape[1])]
            models["mae@null"] = MAE.train_mae(Zs, epochs=epochs, seed=0)
        print(f"[{time.strftime('%H:%M:%S')}] {tname}: ref {Z.shape}, warm {warm.mean():.3f} ({time.time() - t0:.0f}s)",
              flush=True)
        for kind in ("MCAR", "MNAR"):
            for seed in SEEDS:
                for r in evaluate(X, names, models, 30, kind, seed, labels, ceiling="--ceiling" in sys.argv):
                    rows.append({"target": tname, **r})
        (out_dir / f"dev_{tag}.json").write_text(json.dumps({"meta": meta, "rows": rows}, indent=1) + "\n")
    summarize(rows)


def summarize(rows):
    import pandas as pd

    df = pd.DataFrame(rows)
    m = df.groupby(["target", "mask", "method"])["rmse_all"].mean().unstack("method")
    rel = m.div(m["stack_lowrank"], axis=0)
    print(rel.round(4).to_string())
    print("\nmean relative to stack_lowrank (all targets, both masks):")
    print(rel.mean().round(4).to_string())


if __name__ == "__main__":
    main()
