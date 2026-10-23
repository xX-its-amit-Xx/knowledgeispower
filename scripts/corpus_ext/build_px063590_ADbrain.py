"""Corpus dataset px063590_ADbrain -> data/cache/corpus/px063590_ADbrain.npz (cells x UniProt, log2, median-centred).

Source: PRIDE PXD063590 ("Single cell proteomic analysis of human prefrontal cortex in Alzheimer's Disease"; no
publication linked in PRIDE at build time). Homo sapiens: cells dissociated from dorsolateral prefrontal cortex
(Brodmann area 9) of three female donors (control Braak 0, AD Braak III, AD Braak V per metadata), cellenONE
sorting, nPOP digestion, TMTpro 32/35-plex multiplexed DIA-free pSCoPE (prioritized DDA) on an Orbitrap Exploris 480;
Proteome Discoverer 3.2 + CHIMERYS. Carrier (126, 10 ng bulk tissue) and reference (127N, 1 ng) channels.
Files: results_pSCOPE_AD_plate1_matched.txt and results_pSCOPE_AD_plate2_matched.txt (long PSM x channel tables:
FileName, MasterProteinAccessions, AnnotatedSequence, Modifications, qValue, Channel, Abundance, matchingID = <well>_<channel>,
Type = donor ID / carrier / reference / negativectrl, Braak). Only quantified reporter ions have rows (missing =
absent row, i.e. unimputed).
Processing (SCoPE2/scp-style, no imputation; the PSM table carries no protein-level matrix):
  1. keep rows with Type == a donor ID (single-cell channels; carrier, reference and negative-control channels are
     dropped), qValue <= 0.01, Abundance > 0 and a single master protein accession (no ';');
  2. peptide = (FileName, AnnotatedSequence, Modifications); duplicate PSMs of a peptide in one channel -> median;
  3. within each TMT set (FileName) and peptide, divide each cell's value by the mean over that set's single-cell
     channels (removes peptide ionisation efficiency and set effects), then log2;
  4. protein value per cell = median of its peptides' log2 ratios; accession -> canonical UniProt via
     ids.group_accessions.
  The reference channel is not used: it is quantified for only ~4% of PSMs.
Cells are identified as plate<k>:<well>_<channel>. Then the build_corpus.save() logic (cells with < 50 proteins dropped).
Scale caveat: values are log2 ratios to the within-set single-cell mean, not log2 intensities (the corpus z-scores
each protein within a dataset, so only within-dataset covariation is used). TMT ratio compression applies.
Labels: `donor<ID>_Braak<stage>` from the PSM table's Type and Braak columns (no per-cell cell-type annotation).

    python scripts/corpus_ext/build_px063590_ADbrain.py \
        data/cache/corpus_raw/px063590/results_pSCOPE_AD_plate1_matched.txt \
        data/cache/corpus_raw/px063590/results_pSCOPE_AD_plate2_matched.txt
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from scp_graphfs_mcp.ids import group_accessions

OUT = Path(__file__).resolve().parents[2] / "data" / "cache" / "corpus"
USE = ["FileName", "MasterProteinAccessions", "AnnotatedSequence", "Modifications", "qValue", "Channel",
       "Abundance", "matchingID", "Type", "Braak"]


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


parts = []
for plate, path in enumerate(sys.argv[1:], start=1):
    for ch in pd.read_csv(path, sep="\t", usecols=USE, chunksize=2_000_000, dtype={"Type": str, "Braak": str},
                          low_memory=False):
        ch = ch[~ch["Type"].isin(["carrier", "reference", "negativectrl"]) & ch["Type"].notna()]
        ch = ch[(ch["qValue"] <= 0.01) & (ch["Abundance"] > 0)
                & ~ch["MasterProteinAccessions"].astype(str).str.contains(";")]
        ch = ch.assign(cell=f"plate{plate}:" + ch["matchingID"].astype(str),
                       label="donor" + ch["Type"].astype(str) + "_Braak" + ch["Braak"].astype(str),
                       pep=ch["AnnotatedSequence"].astype(str) + "|" + ch["Modifications"].astype(str).fillna(""),
                       FileName=ch["FileName"].astype(str))
        parts.append(ch[["FileName", "pep", "MasterProteinAccessions", "cell", "label", "Abundance"]])
        print(path, len(ch), flush=True)
d = pd.concat(parts, ignore_index=True)
del parts
for c in ["FileName", "pep", "MasterProteinAccessions", "cell", "label"]:
    d[c] = d[c].astype("category")
lab = d.groupby("cell", observed=True)["label"].agg(lambda s: s.iloc[0] if s.nunique() == 1 else "AMBIGUOUS")
assert (lab != "AMBIGUOUS").all()
# 2. median over duplicate PSMs of a peptide in one cell (each cell belongs to exactly one FileName)
pep = d.groupby(["FileName", "pep", "MasterProteinAccessions", "cell"], observed=True)["Abundance"].median().reset_index()
del d
# 3. ratio to the within-set mean of single-cell channels, log2
pep["lr"] = np.log2(pep["Abundance"] / pep.groupby(["FileName", "pep"], observed=True)["Abundance"].transform("mean"))
# 4. protein = median over peptides
prot = pep.groupby(["cell", "MasterProteinAccessions"], observed=True)["lr"].median()
M = prot.unstack("MasterProteinAccessions")
cells = list(M.index.astype(str))
accs = [(group_accessions(str(a)) or [""])[0] for a in M.columns]
save("px063590_ADbrain", M.to_numpy(float), accs, [lab[c] for c in cells])
