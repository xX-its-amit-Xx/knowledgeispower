"""Corpus dataset weiss2026_liver -> data/cache/corpus/weiss2026_liver.npz (cells x UniProt, log2, median-centred).

Source: PRIDE PXD062231 (Weiss et al. 2026, Nat Metab, doi:10.1038/s42255-026-01459-2; "Single-cell spatial
proteomics maps human liver zonation patterns ..."). Homo sapiens hepatocytes from FFPE liver sections of 18
individuals, single-cell Deep Visual Proteomics (scDVP): one laser-microdissected hepatocyte contour per well
(44 per individual, 792 runs), label-free DIA on an Orbitrap Astral (Evosep Whisper 40), DIA-NN 1.8.1 with MBR.
File: search_engine_results/report.pg_matrix.tsv, extracted with HTTP range reads from
search_engine_results.zip (DIA-NN --matrices protein-group matrix; linear normalised MaxLFQ values; missing = empty/NaN,
i.e. unimputed).
Filters: all 792 runs are single-hepatocyte shapes (SDRF `cells per well` = 1); protein groups with no value in
any run dropped. Scale: positive values -> log2. `Protein.Group` -> first canonical UniProt accession via
ids.group_accessions. Then the build_corpus.save() logic (cells with < 50 proteins dropped).
Labels: `donor_<code>` from the run name (`..._Box<k>_<code>_shape_<n>`; 18 codes NML* / CHTL*, 44 shapes each).
The SDRF renames runs to `individual1..18` (14 normal, 4 fibrosis) but repeats (box, shape) keys across individuals,
so code -> SDRF individual / disease cannot be mapped unambiguously; only the donor code is used. Caveat: an FFPE
section contour is a partial-volume single cell, not a whole dissociated cell.

    python scripts/corpus_ext/build_weiss2026_liver.py data/cache/corpus_raw/weiss2026/report.pg_matrix.tsv
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
pat = re.compile(r"_box\d+_(.+?)_shape_\d+\.mzML$", re.I)
cols = [c for c in d.columns if pat.search(c)]
assert len(cols) == 792, len(cols)
labels = [f"donor_{pat.search(c).group(1)}" for c in cols]
X = d[cols].to_numpy(float).T
with np.errstate(divide="ignore", invalid="ignore"):
    X = np.where(X > 0, np.log2(X), np.nan)
accs = [(group_accessions(str(g)) or [""])[0] for g in d["Protein.Group"]]
seen_sc = ~np.isnan(X).all(0)
X, accs = X[:, seen_sc], [a for a, k in zip(accs, seen_sc) if k]
save("weiss2026_liver", X, accs, labels)
