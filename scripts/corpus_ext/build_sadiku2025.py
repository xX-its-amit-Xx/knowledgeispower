"""Corpus dataset sadiku2025 -> data/cache/corpus/sadiku2025.npz (cells x UniProt, log2, median-centred).

Source: PRIDE PXD061065 (Sadiku et al. 2025, Nat Commun, doi:10.1038/s41467-025-67367-3; "Single cell proteomics of
human brain neutrophils from glioblastoma patients"). Homo sapiens primary neutrophils isolated from glioblastoma
tissue of 6 patients (57 wells per patient), label-free DIA on an Orbitrap Astral, Spectronaut 19.4.
Files: 20250204_112949_gbm_sc_full_fasta_19.4_Report.tsv (Spectronaut protein pivot report; `[i] <run>.raw.PG.Quantity`
columns are linear PG quantities, missing = NaN, i.e. unimputed; Spectronaut's default cross-run normalisation may
have been applied) and experimentalDesignTemplate.txt (raw file -> patient, cell type).
Filters: runs listed in experimentalDesignTemplate.txt whose name does not contain `0cells` (12 empty-well controls
dropped); the one `500x` run is not in the design and is dropped; protein groups with no value in any kept cell are
dropped. Scale: positive quantities -> log2. Protein groups (`A;B`) -> first canonical UniProt accession via
ids.group_accessions. Then the build_corpus.save() logic (cells with < 50 proteins dropped).
Labels: `GBM_patient<k>` (all cells are neutrophils; the patient is the only per-cell annotation).

    python scripts/corpus_ext/build_sadiku2025.py \
        data/cache/corpus_raw/sadiku2025/20250204_112949_gbm_sc_full_fasta_19.4_Report.tsv \
        data/cache/corpus_raw/sadiku2025/experimentalDesignTemplate.txt
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
design = pd.read_csv(sys.argv[2], sep="\t")
patient = dict(zip(design["raw file"], design["patient"]))
cols, labels = [], []
for c in d.columns:
    if not c.endswith(".raw.PG.Quantity"):
        continue
    run = re.sub(r"^\[\d+\] ", "", c)[: -len(".PG.Quantity")]
    if run in patient and "0cells" not in run:
        cols.append(c); labels.append(f"GBM_patient{patient[run]}")
X = d[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float).T
with np.errstate(divide="ignore", invalid="ignore"):
    X = np.where(X > 0, np.log2(X), np.nan)
accs = [(group_accessions(str(g)) or [""])[0] for g in d["PG.ProteinGroups"]]
seen_sc = ~np.isnan(X).all(0)
X, accs = X[:, seen_sc], [a for a, k in zip(accs, seen_sc) if k]
save("sadiku2025", X, accs, labels)
