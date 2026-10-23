"""Corpus dataset inns2026_CCS -> data/cache/corpus/inns2026_CCS.npz (cells x UniProt, log2, median-centred).

Source: PRIDE PXD073250 (Inns et al. 2026, Life Sci Alliance, doi:10.26508/lsa.202603759; "Democratized single-cell
proteomics resolves cell state heterogeneity in skin tumors"). Homo sapiens primary cells from one CYLD cutaneous
syndrome (CCS) skin tumour, FACS-sorted into hexadecane-filled 96-well plates (6 plates), label-free diaPASEF on a
timsTOF HT (Evosep Whisper-Zoom 40 SPD), DIA-NN 2.0.1.
File: report.pg_matrix.tsv (DIA-NN protein-group matrix, linear normalised quantities, missing = empty/NaN, i.e.
unimputed; the authors' Seurat normalisation is NOT applied here).
Filters: only `..._plate<k>_singlecell_pop{1,2}_...` runs (529 single cells); true blanks, end/post blanks dropped;
protein groups with no value in any kept cell dropped. Scale: positive values -> log2. `Protein.Group` (`A;B`) ->
first canonical UniProt accession via ids.group_accessions. Then the build_corpus.save() logic (cells with < 50
proteins dropped; the paper's own QC kept 419 cells, a stricter rule not reproduced here).
Labels: FACS gate from the run name. The paper sorted CD45-/CD200+ tumour keratinocytes and CD45+/CD74+ immune
(myeloid) cells; pop1 = keratinocyte gate and pop2 = immune gate, inferred from the data (KRT5/KRT14/KRT17
detected in 80-85% of pop1 vs 38-62% of pop2 cells; PTPRC/CD74/HLA-DRA in 3-11% of pop1 vs 63-76% of pop2).
Labels: `CD45neg_CD200pos_keratinocyte_gate` / `CD45pos_CD74pos_immune_gate`.

    python scripts/corpus_ext/build_inns2026_CCS.py data/cache/corpus_raw/cyld2026/report.pg_matrix.tsv
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from scp_graphfs_mcp.ids import group_accessions

OUT = Path(__file__).resolve().parents[2] / "data" / "cache" / "corpus"
GATE = {"pop1": "CD45neg_CD200pos_keratinocyte_gate", "pop2": "CD45pos_CD74pos_immune_gate"}


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
    m = re.search(r"_plate\d+_singlecell_(pop[12])_", c)
    if m:
        cols.append(c); labels.append(GATE[m.group(1)])
X = d[cols].to_numpy(float).T
with np.errstate(divide="ignore", invalid="ignore"):
    X = np.where(X > 0, np.log2(X), np.nan)
accs = [(group_accessions(str(g)) or [""])[0] for g in d["Protein.Group"]]
seen_sc = ~np.isnan(X).all(0)
X, accs = X[:, seen_sc], [a for a, k in zip(accs, seen_sc) if k]
save("inns2026_CCS", X, accs, labels)
