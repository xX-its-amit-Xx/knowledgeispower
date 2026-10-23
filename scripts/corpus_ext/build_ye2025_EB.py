"""Corpus dataset ye2025_EB -> data/cache/corpus/ye2025_EB.npz (cells x UniProt, log2, median-centred).

Source: PRIDE PXD049181 (project title "Deep single-cell proteomics quantifies low-abundant embryonic transcription
factors and lineage markers during induced pluripotent stem cell differentiation"; paper: Ye et al. 2025, Nat Methods
22:499-509, doi:10.1038/s41592-024-02558-2, "Enhanced sensitivity and scalability with a Chip-Tip workflow enables
deep single-cell proteomics"). Homo sapiens single cells from iPSC-derived embryoid bodies (EBs), cellenONE isolation,
label-free DIA on an Orbitrap Astral (Evosep 40 SPD), Spectronaut.
File: iPSC_EBs_search/EBs/20231218_185013_20231207_AST5_Evo5_ML_SA_40SPD_DIA_SCP_EB_Report_pivot_Protein.tsv,
extracted with HTTP range reads from iPSC_EBs_search.zip (Spectronaut protein pivot; `<run>.PG.Quantity` linear).
The matching setup file (..._EB__Report_.setup.txt) states "Imputation Strategy: None" and
"Cross-Run Normalization: False", so values are unimputed and un-normalised (missing = NaN).
Filters: only `..._SCP_EB_1Cellc_<n>.raw.PG.Quantity` columns (53 single cells); the four 20-cell runs are dropped,
as are protein groups with no value in any kept cell. Scale: positive quantities -> log2. `PG.ProteinGroups` ->
first canonical UniProt accession via ids.group_accessions. Then the build_corpus.save() logic.
Labels: constant `iPSC_EB` (no per-cell annotation is deposited; the PRIDE community SDRF labels these runs
"HeLa", which contradicts the run names and the paper and is ignored).

    python scripts/corpus_ext/build_ye2025_EB.py data/cache/corpus_raw/ye2025/EB_Report_pivot_Protein.tsv
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


d = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)
cols = [c for c in d.columns if re.search(r"_SCP_EB_1Cellc_\d+\.raw\.PG\.Quantity$", c)]
assert len(cols) == 53, len(cols)
X = d[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float).T
with np.errstate(divide="ignore", invalid="ignore"):
    X = np.where(X > 0, np.log2(X), np.nan)
accs = [(group_accessions(str(g)) or [""])[0] for g in d["PG.ProteinGroups"]]
seen_sc = ~np.isnan(X).all(0)
X, accs = X[:, seen_sc], [a for a, k in zip(accs, seen_sc) if k]
save("ye2025_EB", X, accs, ["iPSC_EB"] * len(cols))
