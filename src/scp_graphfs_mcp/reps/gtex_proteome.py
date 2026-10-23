"""GTEx tissue proteome reference panel (Jiang et al. 2020, Cell, "A Quantitative Proteome Map of the Human
Body"; 32 normal tissue types from GTEx donors, TMT 10-plex).

Raw file (``scripts/fetch_features.sh``, Phase B section): ``1-s2.0-S0092867420310783-mmc3.xlsx`` (Table S2 of
the paper, Elsevier supplementary CDN). Worksheet ``D protein relative abundance`` ("Cleaned relative protein
abundances (log2 transformed abundance to the average of the reference in the same run)"); layout verified
2026-10-23: row 1 title, row 2 tissue name per sample column, row 3 header ``gene.id.full, gene.id`` + 201 GTEx
sample ids (``GTEX-YEC4-1326-SM-9KI3A``), then one row per gene (``ENSG00000000003.10``, ``ENSG00000000003``) with
values or ``NA``. Values are already log2 ratios -> kept as provided.

Identifier mapping: Ensembl gene id -> HGNC symbol(s) (``ids.resolve_symbol``) -> UniProt accession(s)
(``ids.uniprot_for_symbol``); an accession claimed by several genes goes to the gene with most observations.

Derived cache: ``gtex_proteome_matrix.npz`` (accs, row_of_acc, genes, samples, tissue, X float32 genes x samples).
API: ``reference(var_names) -> (R, has)``, ``corr_graph(var_names, top=20)``, ``block(var_names, d=24)``.
"""

from __future__ import annotations

from functools import lru_cache

import networkx as nx
import numpy as np

from . import _bulk
from . import _refpanel as rp

NAME = "gtex_proteome"
URL = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867420310783-mmc3.xlsx"
DIR = rp.REPS_DIR / NAME
RAW = DIR / "1-s2.0-S0092867420310783-mmc3.xlsx"
SHEET = "D protein relative abundance"
CACHE = DIR / "gtex_proteome_matrix.npz"
MIN_SHARED = 30   # corr_graph: minimum samples observed for both proteins
MIN_OBS = 30      # block: minimum observed samples for a protein to enter the SVD


def _num(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def _build() -> None:
    tissue_row, header, genes, vals = None, None, [], []
    for row in _bulk.read_xlsx_rows(RAW, SHEET):
        if header is None:
            if row.get(0) == "gene.id.full":
                header = row
                cols = sorted(k for k in row if k >= 2)
            elif len(row) > 10:
                tissue_row = row
            continue
        if not row.get(1):
            continue
        genes.append(row[1])
        vals.append(np.array([_num(row.get(k)) for k in cols], np.float32))
    if header is None:
        raise ValueError(f"header row 'gene.id.full' not found in sheet {SHEET!r} of {RAW.name}")
    X = np.vstack(vals)
    samples = [header[k] for k in cols]
    tissue = [(tissue_row or {}).get(k, "") for k in cols]
    nobs = (~np.isnan(X)).sum(1)
    index = _bulk.build_index([_bulk.gene_accessions(g) for g in genes], nobs)
    keys = sorted(index)
    np.savez_compressed(CACHE, accs=np.array(keys), row_of_acc=np.array([index[k] for k in keys], np.int32),
                        genes=np.array(genes), samples=np.array(samples), tissue=np.array(tissue), X=X,
                        raw_sha256=np.array(rp.sha256(RAW)), url=np.array(URL))


@lru_cache(maxsize=1)
def _panel():
    if not CACHE.exists() or CACHE.stat().st_mtime < RAW.stat().st_mtime:
        _build()
    d = np.load(CACHE)
    index = {str(a): int(r) for a, r in zip(d["accs"], d["row_of_acc"])}
    return d["X"].astype(float), index, [str(x) for x in d["samples"]], str(d["raw_sha256"])


def provenance() -> dict:
    return {"source": "GTEx tissue proteome (Jiang et al. 2020)", "version": f"Table S2 sheet '{SHEET}'",
            "url": URL, "sha256": _panel()[3]}


def reference(var_names):
    """(p x 201) per-protein standardised GTEx log2 relative abundances aligned to var_names (NaN = absent), has."""
    X, index, _, _ = _panel()
    return rp.reference(var_names, X, index)


def corr_graph(var_names, top: int = 20, min_shared: int = MIN_SHARED) -> nx.Graph:
    """Each protein linked to its `top` most positively correlated partners (pairwise-complete Pearson across
    GTEx tissue samples, >= `min_shared` shared samples); weight = r. Isolated nodes kept."""
    X, index, _, _ = _panel()
    return _bulk.corr_graph(var_names, X, index, top, min_shared, provenance(), "samples")


@lru_cache(maxsize=4)
def _loadings(d: int):
    return rp.svd_loadings(_panel()[0], d, MIN_OBS)


def block(var_names, d: int = 24) -> np.ndarray:
    """SVD loadings (d) of the standardised GTEx panel (proteins with >= 30 observed samples, NaN -> 0) + has."""
    L, ok = _loadings(d)
    return _bulk.block(var_names, _panel()[1], L, ok)
