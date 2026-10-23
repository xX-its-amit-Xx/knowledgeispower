"""Corpus dataset bubis2025_A549 -> data/cache/corpus/bubis2025_A549.npz (cells x UniProt, log2, median-centred).

Source: PRIDE PXD064564 (Bubis et al. 2025, Nat Methods, doi:10.1038/s41592-024-02559-1; "Challenging the Astral
mass analyzer - up to 5300 proteins per single-cell"). Homo sapiens A549 (lung adenocarcinoma cell line), label-free
DIA on an Orbitrap Astral (FAIMS CV -48, 50 SPD), cellenONE isolation. Spectronaut directDIA+ per-group MBR, strict FDR.
File: MM_Astral_paper_A549_Figure7_MBR_perGroup_directDIA_strictFDR_Report_Protein_JB_Pivot.tsv (Spectronaut protein
pivot report; `<run>.PG.Quantity` columns are linear PG quantities; missing = NaN / "Filtered", i.e. unimputed;
Spectronaut's own cross-run normalisation may have been applied - the .sne setup is not inspected).
Filters: only `..._50SPD_SC_A549_<well>.raw.PG.Quantity` columns (66 single cells); the 20xSC / 40xSC multi-cell
runs and the 3 blanks are dropped, as are protein groups with no value in any kept cell. Scale: positive quantities -> log2. Row IDs are Spectronaut protein groups
(`A;B`) -> first canonical UniProt accession via ids.group_accessions. Then the build_corpus.save() logic.
Labels: all cells are untreated A549 (labels file is the constant "A549").

    python scripts/corpus_ext/build_bubis2025_A549.py \
        data/cache/corpus_raw/bubis2025/MM_Astral_paper_A549_Figure7_MBR_perGroup_directDIA_strictFDR_Report_Protein_JB_Pivot.tsv
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
cols = [c for c in d.columns if re.search(r"_50SPD_SC_A549_[A-P]\d+\.raw\.PG\.Quantity$", c)]
assert len(cols) == 66, len(cols)
X = d[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float).T
with np.errstate(divide="ignore", invalid="ignore"):
    X = np.where(X > 0, np.log2(X), np.nan)
accs = [(group_accessions(str(g)) or [""])[0] for g in d["PG.ProteinGroups"]]
seen_sc = ~np.isnan(X).all(0)  # groups quantified only in the multi-cell runs / blanks carry no single-cell data
X, accs = X[:, seen_sc], [a for a, k in zip(accs, seen_sc) if k]
save("bubis2025_A549", X, accs, ["A549"] * len(cols))
