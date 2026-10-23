"""Corpus dataset amrute2026_MFS -> data/cache/corpus/amrute2026_MFS.npz (cells x UniProt, log2, median-centred).

Source: PRIDE PXD078860 (Amrute et al. 2026, Arterioscler Thromb Vasc Biol, doi:10.1161/atvbaha.125.324104;
"Single-Cell Proteomics Uncovers Cell-Specific Proteins in Vascular Health and Marfan Syndrome Disease").
Human part only: single cells from Marfan-syndrome (MFS) aneurysm aortic tissue, label-free DIA (Evosep Whisper
100 40 SPD, FAIMS, Orbitrap Eclipse/Ascend), DIA-NN with MBR. The mouse matrices are not used.
Files: report.human_1_pg_matrix.tsv (766 runs `Evosep_JHC1_*`) and report.human_2_pg_matrix.tsv (384 runs
`Evosep_JHC2_*`): DIA-NN protein-group matrices, linear normalised quantities, missing = empty/NaN (unimputed).
The two searches are separate DIA-NN runs; they are concatenated on the first-UniProt-accession key.
Filters: PRIDE has no per-run sample sheet (no SDRF), so empty wells cannot be identified from metadata; runs with
< 200 quantified protein groups (in their own matrix) are dropped as likely empty / failed wells (the paper applied
its own TIC / protein-count / keratin QC, not reproducible from the deposit). Protein groups with no value in any kept
run are dropped. Scale: positive values -> log2. `Protein.Group` -> first canonical UniProt accession via
ids.group_accessions. Then the build_corpus.save() logic.
Labels: batch only (`MFS_aorta_JHC1` / `MFS_aorta_JHC2`); no per-cell cell-type annotation is deposited. Cells are
dominated by smooth-muscle markers (MYH11/ACTA2 detected in ~90-95% of runs); keratin (KRT10) is detected in ~95%
of runs, i.e. contamination is present.

    python scripts/corpus_ext/build_amrute2026_MFS.py data/cache/corpus_raw/amrute2025/report.human_1_pg_matrix.tsv \
        data/cache/corpus_raw/amrute2025/report.human_2_pg_matrix.tsv
"""
import json
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


blocks, labels = [], []
for path, batch in zip(sys.argv[1:3], ["JHC1", "JHC2"]):
    d = pd.read_csv(path, sep="\t")
    runs = [c for c in d.columns if c.endswith(".raw")]
    assert all(f"Evosep_{batch}_" in c for c in runs)
    x = d[runs].to_numpy(float).T
    x = x[(~np.isnan(x) & (x > 0)).sum(1) >= 200]
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.where(x > 0, np.log2(x), np.nan)
    accs = [(group_accessions(str(g)) or [""])[0] for g in d["Protein.Group"]]
    # collapse duplicate accessions within the block (keep the most-detected group), as save() does
    det = (~np.isnan(x)).sum(0)
    best = {}
    for i, a in enumerate(accs):
        if a and (a not in best or det[i] > det[best[a]]):
            best[a] = i
    blocks.append(pd.DataFrame(x[:, list(best.values())], columns=list(best.keys())))
    labels += [f"MFS_aorta_{batch}"] * x.shape[0]
M = pd.concat(blocks, axis=0, join="outer", ignore_index=True)
M = M.loc[:, M.notna().any(axis=0)]
save("amrute2026_MFS", M.to_numpy(float), list(M.columns), labels)
