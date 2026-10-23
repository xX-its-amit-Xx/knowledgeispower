"""Sanger ProCan pan-cancer cell-line proteome reference panel (Gonçalves et al. 2022, Cancer Cell; ProCan-
DepMapSanger: 949 cancer cell lines, DIA-MS (DIA-NN, MaxLFQ), 8,498 proteins).

Raw files (``scripts/fetch_features.sh``, Phase B section), from figshare 10.6084/m9.figshare.19345397 (v1):

* ``ProCan-DepMapSanger_protein_matrix_8498_averaged.txt``: tab-separated, **cell lines as rows** (first column
  ``Project_Identifier`` = ``SIDMxxxxx;<cell line>``), **proteins as columns** named ``<UniProt acc>;<entry name>``
  (e.g. ``P37108;SRP14_HUMAN``; canonical accessions only, no isoforms, no protein groups). Values are log2
  MaxLFQ intensities averaged across replicates (paper STAR Methods: "Data were then log2 - transformed"); empty
  = missing. Kept as provided (already log2).
* ``ProCan-DepMapSanger_mapping_file_averaged.txt``: per-line metadata (``Cell_line, SIDM, Project_Identifier,
  Tissue_type, Cancer_type, Cancer_subtype``); only ``Tissue_type`` is cached (for reference).

The 8,498-protein matrix includes proteins quantified from a single peptide; the authors also publish a
6,692-protein (>1 peptide) subset. We use the full matrix for coverage.

Derived cache: ``procan_matrix.npz`` (accs, lines, tissue, X float32 proteins x lines, raw sha256).
API: ``reference(var_names) -> (R, has)``, ``corr_graph(var_names, top=20)``, ``block(var_names, d=24)``.
"""

from __future__ import annotations

from functools import lru_cache

import networkx as nx
import numpy as np
import pandas as pd

from ..ids import canonical_accession
from . import _bulk
from . import _refpanel as rp

NAME = "procan"
FIGSHARE = "https://doi.org/10.6084/m9.figshare.19345397.v1"
URL = "https://ndownloader.figshare.com/files/34411172"          # ..._protein_matrix_8498_averaged.txt
URL_META = "https://ndownloader.figshare.com/files/34411133"     # ..._mapping_file_averaged.txt
DIR = rp.REPS_DIR / NAME
RAW = DIR / "ProCan-DepMapSanger_protein_matrix_8498_averaged.txt"
RAW_META = DIR / "ProCan-DepMapSanger_mapping_file_averaged.txt"
CACHE = DIR / "procan_matrix.npz"
MIN_SHARED = 50   # corr_graph: minimum lines observed for both proteins
MIN_OBS = 50      # block: minimum observed lines for a protein to enter the SVD


def _build() -> None:
    df = pd.read_csv(RAW, sep="\t", index_col=0)
    X = df.to_numpy(np.float32).T                     # proteins x lines
    accs_raw = [canonical_accession(c.split(";")[0]) for c in df.columns]
    nobs = (~np.isnan(X)).sum(1)
    index = _bulk.build_index([[a] if a else [] for a in accs_raw], nobs)
    keys = sorted(index)
    rows = [index[k] for k in keys]
    lines = [str(s) for s in df.index]
    tissue = np.array([""] * len(lines), dtype=object)
    if RAW_META.exists():
        meta = pd.read_csv(RAW_META, sep="\t", dtype=str).set_index("Project_Identifier")
        tissue = np.array([meta["Tissue_type"].get(s, "") for s in lines], dtype=object)
    np.savez_compressed(CACHE, accs=np.array(keys), lines=np.array(lines), tissue=tissue.astype(str),
                        X=X[rows], raw_sha256=np.array(rp.sha256(RAW)), url=np.array(URL))


@lru_cache(maxsize=1)
def _panel():
    if not CACHE.exists() or CACHE.stat().st_mtime < RAW.stat().st_mtime:
        _build()
    d = np.load(CACHE)
    accs = [str(a) for a in d["accs"]]
    return d["X"].astype(float), {a: i for i, a in enumerate(accs)}, [str(x) for x in d["lines"]], str(d["raw_sha256"])


def provenance() -> dict:
    return {"source": "Sanger ProCan-DepMapSanger (Gonçalves et al. 2022)",
            "version": "figshare 19345397 v1, protein_matrix_8498_averaged", "url": URL, "sha256": _panel()[3]}


def reference(var_names):
    """(p x 949) per-protein standardised ProCan log2 abundances aligned to var_names (NaN where absent), has."""
    X, index, _, _ = _panel()
    return rp.reference(var_names, X, index)


def corr_graph(var_names, top: int = 20, min_shared: int = MIN_SHARED) -> nx.Graph:
    """Each protein linked to its `top` most positively correlated partners (pairwise-complete Pearson across
    ProCan lines, >= `min_shared` shared lines); weight = r. Isolated nodes kept."""
    X, index, _, _ = _panel()
    return _bulk.corr_graph(var_names, X, index, top, min_shared, provenance(), "lines")


@lru_cache(maxsize=4)
def _loadings(d: int):
    return rp.svd_loadings(_panel()[0], d, MIN_OBS)


def block(var_names, d: int = 24) -> np.ndarray:
    """SVD loadings (d) of the standardised ProCan panel (proteins with >= 50 observed lines, NaN -> 0) + has."""
    L, ok = _loadings(d)
    return _bulk.block(var_names, _panel()[1], L, ok)
