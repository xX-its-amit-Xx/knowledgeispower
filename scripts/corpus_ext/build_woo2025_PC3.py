"""Corpus dataset woo2025_PC3 -> data/cache/corpus/woo2025_PC3.npz (cells x UniProt, log2, median-centred).

Source: PRIDE PXD059079 (Woo et al. 2025, Mol Cell Proteomics, doi:10.1016/j.mcpro.2025.100949; "Single-Cell
Proteomic Characterization of Drug-Resistant Prostate Cancer Cells ..."). Homo sapiens PC-3 prostate cancer cell line:
untreated control cells (CONT) and cisplatin-induced polyaneuploid cancer cells (PACC). Label-free DIA, Orbitrap
Ascend + FAIMS, 60 SPD; DIA-NN.
File: Result_pg_matrix.tsv (DIA-NN protein-group matrix, 1590 groups x 144 runs; linear, DIA-NN-normalised MaxLFQ
values; missing = empty/NaN, i.e. unimputed; 19.5% missing, so the authors' DIA-NN output appears already
restricted to well-detected groups).
Filters: every quantity column is a single-cell run (`..._PC3_{CONT,PACC}_SC_<well>.raw`, 72 + 72; SDRF
`cells per well` = 1). Scale: positive values -> log2. `Protein.Group` (`A;B`) -> first canonical UniProt accession
via ids.group_accessions. Then the build_corpus.save() logic.
Labels: `PC3_CONT` / `PC3_PACC` parsed from the run name.

    python scripts/corpus_ext/build_woo2025_PC3.py data/cache/corpus_raw/px059079/Result_pg_matrix.tsv
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from scp_graphfs_mcp.ids import group_accessions

OUT = Path(__file__).resolve().parents[2] / "data" / "cache" / "corpus"


def save(name, X, accs, labels=None):  # build_corpus.py save() + aligned per-cell labels
    OUT.mkdir(parents=True, exist_ok=True)
    X = np.asarray(X, float)
    keep_cols = [i for i, a in enumerate(accs) if a]
    X, accs = X[:, keep_cols], [accs[i] for i in keep_cols]
    order = np.argsort(-(~np.isnan(X)).sum(0))
    seen, cols = set(), []
    for i in order:
        if accs[i] not in seen:
            seen.add(accs[i]); cols.append(i)
    cols.sort()
    X, accs = X[:, cols], [accs[i] for i in cols]
    rows = (~np.isnan(X)).sum(1) >= 50
    X = X[rows]
    X = X - np.nanmedian(X, axis=1, keepdims=True)
    np.savez_compressed(OUT / f"{name}.npz", X=X.astype(np.float32), accs=np.array(accs))
    if labels is not None:
        (OUT / f"{name}_labels.json").write_text(json.dumps([str(l) for l, k in zip(labels, rows) if k]))
    print(f"{name}: {X.shape[0]} cells x {X.shape[1]} proteins, missing {np.isnan(X).mean():.3f}, "
          f"p1/p50/p99 {np.round(np.nanpercentile(X, [1, 50, 99]), 3).tolist()}")


d = pd.read_csv(sys.argv[1], sep="\t")
cols, labels = [], []
for c in d.columns:
    m = re.search(r"_PC3_(CONT|PACC)_SC_[A-P]\d+\.raw$", c)
    if m:
        cols.append(c); labels.append(f"PC3_{m.group(1)}")
assert len(cols) == 144, len(cols)
X = d[cols].to_numpy(float).T
with np.errstate(divide="ignore", invalid="ignore"):
    X = np.where(X > 0, np.log2(X), np.nan)
save("woo2025_PC3", X, [(group_accessions(str(g)) or [""])[0] for g in d["Protein.Group"]], labels)
