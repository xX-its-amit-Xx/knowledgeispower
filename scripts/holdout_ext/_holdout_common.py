"""Shared writer for the v9 holdout builders (scripts/holdout_ext/build_*.py).

Processing is the corpus save() logic of scripts/corpus_ext/build_woo2025_PC3.py (first canonical UniProt accession
per protein group; duplicate accessions -> keep the most-observed column; drop cells with < 50 observed proteins;
per-cell median-centring of log2 values; NaN = missing), plus aligned per-cell labels.

HOLDOUT INTEGRITY: outputs go ONLY to data/cache/holdout/<name>.npz and data/datasets/<name>_holdout.h5ad.
Nothing is ever written under data/cache/corpus/ (every .npz there is pulled into the external reference).
Only QC counts are printed: no imputation, correlation or structure statistic is computed here.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from scp_graphfs_mcp.datasets import from_log2_matrix

ROOT = Path(__file__).resolve().parents[2]
HOLDOUT = ROOT / "data" / "cache" / "holdout"
DATASETS = ROOT / "data" / "datasets"
CORPUS = ROOT / "data" / "cache" / "corpus"
RAW = ROOT / "data" / "cache" / "corpus_raw"


def _guard(path: Path) -> Path:
    p = path.resolve()
    if p == CORPUS.resolve() or CORPUS.resolve() in p.parents:
        raise SystemExit(f"refusing to write into the reference corpus: {p}")
    return p


def to_log2(X):
    X = np.asarray(X, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(X > 0, np.log2(X), np.nan)


def save_holdout(name, X, accs, cells, labels=None, extra_obs=None, uns=None, min_obs=50):
    """X: cells x groups (log2, NaN = missing); accs: first accession per group ('' = unmappable)."""
    X = np.asarray(X, float)
    cells = np.asarray(cells, dtype=str)
    accs = [str(a) for a in accs]
    keep_cols = [i for i, a in enumerate(accs) if a]
    X, accs = X[:, keep_cols], [accs[i] for i in keep_cols]
    order = np.argsort(-(~np.isnan(X)).sum(0), kind="stable")
    seen, cols = set(), []
    for i in order:
        if accs[i] not in seen:
            seen.add(accs[i]); cols.append(i)
    cols.sort()
    X, accs = X[:, cols], [accs[i] for i in cols]
    rows = (~np.isnan(X)).sum(1) >= min_obs
    n_dropped = int((~rows).sum())
    X, cells = X[rows], cells[rows]
    labels = None if labels is None else np.asarray(labels, dtype=str)[rows]
    X = X - np.nanmedian(X, axis=1, keepdims=True)
    pseen = ~np.isnan(X).all(0)  # proteins never observed in a kept cell
    X, accs = X[:, pseen], [a for a, k in zip(accs, pseen) if k]

    HOLDOUT.mkdir(parents=True, exist_ok=True)
    npz = _guard(HOLDOUT / f"{name}.npz")
    arrays = {"X": X.astype(np.float32), "accs": np.array(accs), "cells": cells}
    if labels is not None:
        arrays["labels"] = labels
    np.savez_compressed(npz, **arrays)

    obs = pd.DataFrame(index=pd.Index(cells, name="cell"))
    obs["cell_type"] = labels if labels is not None else name
    if extra_obs is not None:
        for k, v in extra_obs.items():
            obs[k] = np.asarray(v, dtype=str)[rows]
    a = from_log2_matrix(X, cells, accs, obs=obs)
    a.uns["dataset"] = {"name": f"{name}_holdout", "role": "HOLDOUT v9: unseen until preregistration is committed",
                        "value": "log2, per-cell median-centred; NaN = missing; var = first UniProt accession",
                        **(uns or {})}
    h5 = _guard(DATASETS / f"{name}_holdout.h5ad")
    a.write_h5ad(h5, compression="gzip")
    qc(name, X, accs, labels, n_dropped, npz, h5)


def qc(name, X, accs, labels, n_dropped, npz, h5):
    obsd = ~np.isnan(X)
    ref = set()
    for f in sorted(CORPUS.glob("*.npz")):
        ref |= set(np.load(f, allow_pickle=False)["accs"].astype(str).tolist())
    print(f"{name}: {X.shape[0]} cells x {X.shape[1]} proteins | missing {np.isnan(X).mean():.3f} | "
          f">=50% cells {(obsd.mean(0) >= 0.5).sum()} | cells dropped (<50 obs) {n_dropped} | "
          f"overlap w/ corpus union {len(set(accs) & ref)}/{len(accs)} (corpus union {len(ref)})")
    if labels is not None:
        print("  labels", pd.Series(labels).value_counts().sort_index().to_dict())
    print("  wrote", npz.relative_to(ROOT), "and", h5.relative_to(ROOT))
