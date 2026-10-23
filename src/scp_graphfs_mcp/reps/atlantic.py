"""ATLANTiC cell-line proteome reference panel (Frejno et al. 2020, Nat Commun, "Proteome activity landscapes of
tumor cell lines determine drug responses"): NCI60 (60 lines) + CRC65 (65 colorectal lines), label-free
MaxQuant, trypsin digests -> 125 cell lines.

Raw file (``scripts/fetch_features.sh``, Phase B section): ``41467_2020_17336_MOESM6_ESM.xlsx`` (Supplementary
Data 3). Sheets used (legend verified 2026-10-23): ``Annotation NCI60 Trypsin`` / ``Annotation CRC65 Trypsin``
(MaxQuant proteinGroups columns: ``ID, Protein IDs, Majority protein IDs, ..., Only identified by site,
Potential contaminant``) and ``NCI60 Trypsin`` / ``CRC65 Trypsin`` ("log10-transformed abundance [iBAQ]"; first
column ``ID`` = ``<MaxQuant id> <gene names>``, then one column per cell line; empty = missing).

Quirks handled: the gene-name part of the data-sheet ``ID`` is corrupted for ~200 rows (every "NA" substring
removed, e.g. ``98 GT3`` for GNAT3), so rows are joined to the annotation by the numeric MaxQuant id only; the
NCI60 data sheet also contains ~1,500 id-only rows without values (dropped). Contaminants and "only identified
by site" groups are dropped. log10 iBAQ -> log2 (x * log2(10)).

The two panels come from separate MaxQuant searches; a protein is matched across them by UniProt accession
(every canonical accession of the Majority protein IDs; within a panel an accession claimed by several groups
goes to the group with most observations). ``reference(..., center="panel")`` removes each panel's per-protein
mean first (the panels were measured separately).

Derived cache: ``atlantic_matrix.npz``. API: ``reference(var_names) -> (R, has)``, ``corr_graph``, ``block``.
"""

from __future__ import annotations

import warnings
from functools import lru_cache

import networkx as nx
import numpy as np

from ..ids import group_accessions
from . import _bulk
from . import _refpanel as rp

NAME = "atlantic"
URL = ("https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-020-17336-9/"
       "MediaObjects/41467_2020_17336_MOESM6_ESM.xlsx")
DIR = rp.REPS_DIR / NAME
RAW = DIR / "41467_2020_17336_MOESM6_ESM.xlsx"
CACHE = DIR / "atlantic_matrix.npz"
PANELS = ("NCI60", "CRC65")
MIN_SHARED = 30
MIN_OBS = 30


def _num(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def _read_panel(panel: str):
    ann = {}
    for r in _bulk.read_xlsx_rows(RAW, f"Annotation {panel} Trypsin"):
        if r.get(0, "ID") == "ID" or len(r) <= 1:
            continue
        if r.get(21) == "+" or r.get(22) == "+":            # only identified by site / contaminant
            continue
        ann[int(r[0].split()[0])] = group_accessions(r.get(2, ""))
    header, ids, vals = None, [], []
    for r in _bulk.read_xlsx_rows(RAW, f"{panel} Trypsin"):
        if r.get(0) == "ID":
            header = r
            cols = sorted(k for k in r if k >= 1)
            continue
        if header is None or len(r) <= 1:
            continue
        i = int(r[0].split()[0])
        if i in ann:
            ids.append(i)
            vals.append([_num(r.get(k)) for k in cols])
    X = np.array(vals, np.float32) * np.float32(np.log2(10.0))
    lines = [f"{panel}:{header[k]}" for k in cols]
    nobs = (~np.isnan(X)).sum(1)
    return X, _bulk.build_index([ann[i] for i in ids], nobs), lines


def _build() -> None:
    parts = [_read_panel(p) for p in PANELS]
    accs = sorted(set().union(*[set(ix) for _, ix, _ in parts]))
    combos: dict[tuple, int] = {}
    row_of_acc = []
    for a in accs:
        key = tuple(ix.get(a, -1) for _, ix, _ in parts)
        row_of_acc.append(combos.setdefault(key, len(combos)))
    X = np.full((len(combos), sum(len(l) for _, _, l in parts)), np.nan, np.float32)
    keys = sorted(combos, key=combos.get)
    j = 0
    for p, (Xp, _, lines) in enumerate(parts):
        rows = np.array([k[p] for k in keys])
        ok = rows >= 0
        X[np.flatnonzero(ok), j:j + len(lines)] = Xp[rows[ok]]
        j += len(lines)
    lines = sum((l for _, _, l in parts), [])
    np.savez_compressed(CACHE, accs=np.array(accs), row_of_acc=np.array(row_of_acc, np.int32), X=X,
                        lines=np.array(lines), raw_sha256=np.array(rp.sha256(RAW)), url=np.array(URL))


@lru_cache(maxsize=1)
def _panel():
    if not CACHE.exists() or CACHE.stat().st_mtime < RAW.stat().st_mtime:
        _build()
    d = np.load(CACHE)
    index = {str(a): int(r) for a, r in zip(d["accs"], d["row_of_acc"])}
    return d["X"].astype(float), index, [str(x) for x in d["lines"]], str(d["raw_sha256"])


def provenance() -> dict:
    return {"source": "ATLANTiC NCI60 + CRC65 proteomes (Frejno et al. 2020)",
            "version": "Supplementary Data 3, sheets 'NCI60 Trypsin' + 'CRC65 Trypsin'", "url": URL,
            "sha256": _panel()[3]}


def reference(var_names, center: str = "none"):
    """(p x 125) per-protein standardised log2 iBAQ aligned to var_names (NaN = absent), has.
    center="panel" subtracts each panel's (NCI60 / CRC65) per-protein mean first."""
    X, index, lines, _ = _panel()
    if center == "panel":
        X = X.copy()
        lab = np.array([s.split(":")[0] for s in lines])
        for p in PANELS:
            m = lab == p
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                X[:, m] -= np.nanmean(X[:, m], 1, keepdims=True)
    elif center != "none":
        raise ValueError(center)
    return rp.reference(var_names, X, index)


def corr_graph(var_names, top: int = 20, min_shared: int = MIN_SHARED) -> nx.Graph:
    """Each protein linked to its `top` most positively correlated partners (pairwise-complete Pearson across
    the 125 lines, >= `min_shared` shared lines); weight = r. Isolated nodes kept."""
    X, index, _, _ = _panel()
    return _bulk.corr_graph(var_names, X, index, top, min_shared, provenance(), "lines")


@lru_cache(maxsize=4)
def _loadings(d: int):
    return rp.svd_loadings(_panel()[0], d, MIN_OBS)


def block(var_names, d: int = 24) -> np.ndarray:
    """SVD loadings (d) of the standardised ATLANTiC panel (proteins with >= 30 observed lines, NaN -> 0) + has."""
    L, ok = _loadings(d)
    return _bulk.block(var_names, _panel()[1], L, ok)
