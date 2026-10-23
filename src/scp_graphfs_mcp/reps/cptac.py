"""CPTAC pan-cancer protein abundance reference panel: the Baylor College of Medicine (BCM) harmonised CPTAC
global-proteome tables (gene-level log2 "reference intensity normalised" TMT abundances), 10 cancer types,
tumour + normal adjacent tissue.

Raw files (``scripts/fetch_features.sh``, Phase B section), from Zenodo record 8394329 (doi:10.5281/zenodo.8394329,
CC-BY-4.0; the data store the ``cptac`` Python package (v1.5.x) downloads from), into ``data/cache/reps/cptac/``:

* ``bcm-<cohort>-proteomics-<COHORT>_proteomics_gene_abundance_log2_reference_intensity_normalized_<Tumor|Normal>.txt.gz``
  for brca, ccrcc, coad, gbm, hnscc, lscc, luad, ov, pdac, ucec (18 files; brca and gbm have no Normal file).
  Layout (verified 2026-10-23): tab-separated, first column ``idx`` = versioned Ensembl gene id, one column per
  sample (CPTAC case id, e.g. ``C3N-01799``; tumour and normal of one case share the id), ``NA`` = missing.
  Values are log2 abundances -> kept as provided.
* ``bcm-brca-mapping-gencode.v34.basic.annotation-mapping.txt.gz`` (``transcript, protein, gene, gene_name, coding,
  MANE, SwissProt``): used only as a fallback for gene ids that HGNC does not resolve.

Samples are labelled ``<cohort>:<T|N>:<case id>`` (1,549 samples in total, all kept: below the ~2,000 cap).
Identifier mapping: Ensembl gene id -> HGNC symbol(s) -> UniProt (``ids.resolve_symbol`` / ``uniprot_for_symbol``),
fallback: SwissProt accessions the BCM GENCODE v34 mapping lists for that gene.

Caveat: proteins are quantified per cohort (separate TMT studies, different protein coverage); pooling cohorts
means between-cohort (tissue + batch) differences contribute to correlations. ``reference(..., center="cohort")``
removes each cohort's per-protein mean first (cohort = cancer type x tumour/normal).

Derived cache: ``cptac_matrix.npz``. API: ``reference(var_names) -> (R, has)``, ``corr_graph``, ``block``.
"""

from __future__ import annotations

import gzip
import warnings
from functools import lru_cache

import networkx as nx
import numpy as np
import pandas as pd

from . import _bulk
from . import _refpanel as rp

NAME = "cptac"
ZENODO = "https://zenodo.org/api/records/8394329"
COHORTS = ("brca", "ccrcc", "coad", "gbm", "hnscc", "lscc", "luad", "ov", "pdac", "ucec")
DIR = rp.REPS_DIR / NAME
MAPPING = DIR / "bcm-brca-mapping-gencode.v34.basic.annotation-mapping.txt.gz"
CACHE = DIR / "cptac_matrix.npz"
MIN_SHARED = 50   # corr_graph: minimum samples observed for both proteins
MIN_OBS = 50      # block: minimum observed samples for a protein to enter the SVD


def raw_files() -> list:
    """The BCM proteomics files present locally, sorted (cohort, Tumor/Normal)."""
    return sorted(DIR.glob("bcm-*-proteomics-*_proteomics_gene_abundance_log2_reference_intensity_normalized_*.txt.gz"))


RAW = DIR / "bcm-brca-proteomics-BRCA_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt.gz"


def url(fname: str) -> str:
    return f"{ZENODO}/files/{fname}/content"


def _swissprot_fallback() -> dict[str, list[str]]:
    if not MAPPING.exists():
        return {}
    m = pd.read_csv(MAPPING, sep="\t", dtype=str, usecols=["gene", "SwissProt"]).dropna()
    out: dict[str, list[str]] = {}
    for g, a in zip(m["gene"].str.split(".").str[0], m["SwissProt"]):
        for x in str(a).replace(",", ";").split(";"):
            x = x.strip().split("-")[0]
            if x and x not in out.setdefault(g, []):
                out[g].append(x)
    return out


def _build() -> None:
    files = raw_files()
    frames, samples, cohort = [], [], []
    for f in files:
        parts = f.name.split("-")
        coh = parts[1]
        kind = "T" if f.name.endswith("_Tumor.txt.gz") else "N"
        with gzip.open(f, "rt") as fh:
            df = pd.read_csv(fh, sep="\t", index_col=0, na_values=["NA"])
        df.index = df.index.astype(str).str.split(".").str[0]
        df = df.groupby(level=0).mean()                      # (version-stripped duplicates, if any)
        df.columns = [f"{coh}:{kind}:{c}" for c in df.columns]
        frames.append(df.astype(np.float32))
        samples += list(df.columns)
        cohort += [f"{coh}:{kind}"] * df.shape[1]
    genes = sorted(set().union(*[set(d.index) for d in frames]))
    X = np.full((len(genes), len(samples)), np.nan, np.float32)
    pos = {g: i for i, g in enumerate(genes)}
    j = 0
    for d in frames:
        rows = np.array([pos[g] for g in d.index])
        X[rows, j:j + d.shape[1]] = d.to_numpy(np.float32)
        j += d.shape[1]
    fb = _swissprot_fallback()
    acc_lists = [(_bulk.gene_accessions(g) or fb.get(g, [])) for g in genes]
    nobs = (~np.isnan(X)).sum(1)
    index = _bulk.build_index(acc_lists, nobs)
    keys = sorted(index)
    shas = [rp.sha256(f) for f in files]
    np.savez_compressed(CACHE, accs=np.array(keys), row_of_acc=np.array([index[k] for k in keys], np.int32),
                        genes=np.array(genes), samples=np.array(samples), cohort=np.array(cohort), X=X,
                        files=np.array([f.name for f in files]), raw_sha256=np.array(shas))


@lru_cache(maxsize=1)
def _panel():
    files = raw_files()
    if not CACHE.exists() or any(CACHE.stat().st_mtime < f.stat().st_mtime for f in files):
        _build()
    d = np.load(CACHE)
    index = {str(a): int(r) for a, r in zip(d["accs"], d["row_of_acc"])}
    return d["X"].astype(float), index, [str(x) for x in d["samples"]], [str(x) for x in d["cohort"]]


def provenance() -> dict:
    d = np.load(CACHE)
    return {"source": "CPTAC pan-cancer proteome, BCM harmonised gene-level tables (via cptac package data)",
            "version": "Zenodo 8394329", "url": ZENODO,
            "sha256": {str(f): str(s) for f, s in zip(d["files"], d["raw_sha256"])}}


def _centered(X: np.ndarray, cohort: list[str]) -> np.ndarray:
    X = X.copy()
    lab = np.array(cohort)
    for c in np.unique(lab):
        m = lab == c
        with warnings.catch_warnings():          # proteins absent from a cohort: all-NaN slice
            warnings.simplefilter("ignore", RuntimeWarning)
            X[:, m] -= np.nanmean(X[:, m], 1, keepdims=True)
    return X


def reference(var_names, center: str = "none"):
    """(p x 1549) per-protein standardised CPTAC log2 abundances aligned to var_names (NaN = absent), has.
    center="cohort" subtracts each cohort's per-protein mean before standardising."""
    X, index, _, cohort = _panel()
    if center == "cohort":
        X = _centered(X, cohort)
    elif center != "none":
        raise ValueError(center)
    return rp.reference(var_names, X, index)


def corr_graph(var_names, top: int = 20, min_shared: int = MIN_SHARED) -> nx.Graph:
    """Each protein linked to its `top` most positively correlated partners (pairwise-complete Pearson across
    CPTAC samples, >= `min_shared` shared samples); weight = r. Isolated nodes kept."""
    X, index, _, _ = _panel()
    return _bulk.corr_graph(var_names, X, index, top, min_shared, provenance(), "samples")


@lru_cache(maxsize=4)
def _loadings(d: int):
    return rp.svd_loadings(_panel()[0], d, MIN_OBS)


def block(var_names, d: int = 24) -> np.ndarray:
    """SVD loadings (d) of the standardised CPTAC panel (proteins with >= 50 observed samples, NaN -> 0) + has."""
    L, ok = _loadings(d)
    return _bulk.block(var_names, _panel()[1], L, ok)
