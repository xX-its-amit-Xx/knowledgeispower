"""Build data/datasets/derks2022_holdout.h5ad: a reserved HOLDOUT (Phase B) dataset, from UNIMPUTED precursors.

Source: scpdata `derks2022` (Derks et al. 2022, plexDIA with mTRAQ Δ0/Δ4/Δ8; PDAC, melanoma and U-937 monocyte
single cells prepared with nPOP, plus 100-cell bulk samples; Q-Exactive and timsTOF SCP runs). The raw file is
data/cache/holdout/derks2022.Rda; its SHA-256 is checked below before anything is read.

Choices, and why:
1. Never use the `proteins` assay. It is IMPUTED (0% NaN), so a missing-value benchmark cannot use it.
2. Input = `qe_prec_extracted` (8590 precursors x 125 run.channel columns). This is the authors' precursor x cell
   matrix for the Q-Exactive single-cell runs (wJD1146-1193), taken after their DIA-NN/channel-level filtering and
   before protein roll-up or imputation (78.7% NaN). rowData gives Protein.Group / Genes for every precursor. The
   per-run assays (`F..JD.plexDIA.*`) hold the same runs before that filtering. They also contain the carrier,
   negative and bulk channels, so we do not use them.
3. Instrument: the timsTOF SCP cells come in a separate matrix (`tims_prec_extracted`, 30 cells). They differ in
   structure from the Q-Exactive cells on every axis: a different instrument, LC column (AUR2-25075 vs AUR2-15075),
   sort day, LC-MS batch and digest batch, and their raw intensity scale is about 10x lower (median ~1.0e3 vs
   ~9.6e3). We judged them too different to pool, so only the larger Q-Exactive set (125 cells) is kept. This was
   decided from metadata and scale only; no correlation or other analysis was run.
4. Cells: keep Real_single_cell == "1" and Celltype in {PDAC, Melanoma, U-937}. This excludes the "_DB" 100-cell
   bulk samples (Real_single_cell 0, which sit in `bulk_prec_extracted` anyway), the "_t" timsTOF cells, and any
   carrier or negative-control channel. All 125 Q-Exactive columns pass.
5. Labels are joined by the key File.Name + "." + Label (mTRAQ channel). This key equals the assay column names.
   The colData `rownames` vector, read by our R-free readers, does NOT line up with the colData rows. For example,
   with rownames, Q-Exactive columns resolve to timsTOF rows. The File.Name/Label/Raw/Instrument fields agree
   with each other within each row, so the key is reliable and the rownames are not.
6. Zeros in a DIA-NN quantity mean "not quantified". So value <= 0 counts as missing, as does NaN. The rest
   become log2.
7. Precursor -> protein: rowData `Protein.Group` (UniProt accessions; isoform/ambiguous groups such as
   "P07954;P07954-2" stay as one protein-group ID, matching the brunner2022 holdout). For each cell, a protein
   = the median of its log2 precursor values, ignoring missing ones. With no observed precursor it is NaN.
8. Each cell's protein vector is median-centred (nanmedian over its observed proteins). This matches the
   corpus and the brunner2022 holdout. Proteins never observed in any kept cell are dropped. X = log2 with
   NaN = missing; layers['raw'] = 2**X (via from_log2_matrix).
Nothing is imputed, correlated or evaluated here. The script prints structure only.

    python scripts/build_derks2022_holdout.py [path/to/derks2022.Rda]
"""
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from scp_graphfs_mcp.datasets import from_log2_matrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "corpus_ext"))
from rda_stream import dataframe, load_rda, qf_assays, sce_matrix  # noqa: E402

RDA = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "data/cache/holdout/derks2022.Rda"
RDA_SHA = "92f157cbe77978507896ab0504c1c3e0afeda5cbf1072cc28ffd77ab9aef9a8c"
ASSAY = "qe_prec_extracted"
CELLTYPES = ("PDAC", "Melanoma", "U-937")

assert hashlib.sha256(RDA.read_bytes()).hexdigest() == RDA_SHA, "derks2022.Rda checksum mismatch"
qf = load_rda(str(RDA), max_str=2_000_000)[0][1]
sce = qf_assays(qf)[ASSAY]
m, dim = sce_matrix(sce)
P = np.asarray(m.value, dtype=float).reshape(dim[1], dim[0])  # cells x precursors
cols = np.array(sce.attr["colData"].attr["rownames"].value)
_, rd = dataframe(sce.attr["rowRanges"].attr["elementMetadata"])
prot = np.array(rd["Protein.Group"], dtype=str)
assert len(prot) == P.shape[1] and (prot != "").all()

_, cd = dataframe(qf.attr["colData"])
cd = pd.DataFrame({k: ["" if x is None else str(x) for x in v] for k, v in cd.items()})
cd.index = cd["File.Name"] + "." + cd["Label"]  # see docstring point 5
assert cd.index.is_unique and set(cols) <= set(cd.index)
meta = cd.loc[cols]
keep = ((meta["Real_single_cell"] == "1") & meta["Celltype"].isin(CELLTYPES)
        & (meta["Instrument"] == "Q-Exactive")).to_numpy()
P, meta = P[keep], meta[keep]

with np.errstate(divide="ignore", invalid="ignore"):
    L = np.where(P > 0, np.log2(P), np.nan)
codes, proteins = pd.factorize(prot, sort=True)
X = np.full((L.shape[0], len(proteins)), np.nan)
for j in range(len(proteins)):
    blk = L[:, codes == j]
    obs = ~np.isnan(blk).all(1)
    X[obs, j] = np.nanmedian(blk[obs], axis=1)
X = X - np.nanmedian(X, axis=1, keepdims=True)
seen = ~np.isnan(X).all(0)
X, proteins = X[:, seen], np.asarray(proteins)[seen]

cells = meta.index.to_numpy()
obs = pd.DataFrame({"cell_type": meta["Celltype"].to_numpy(), "instrument": meta["Instrument"].to_numpy(),
                    "mtraq_label": meta["Label"].to_numpy(), "run": meta["Raw"].to_numpy()}, index=cells)
a = from_log2_matrix(X, cells, proteins, obs=obs)
a.uns["dataset"] = {"name": "derks2022_holdout",
                    "source": f"scpdata derks2022.Rda `{ASSAY}` (unimputed precursors; sha256 {RDA_SHA})",
                    "value": "log2; protein = median of log2 precursors (Protein.Group); per-cell median-centred; "
                             "NaN = missing",
                    "role": "HOLDOUT (Phase B): confirmation only; Q-Exactive single cells only"}
out = ROOT / "data" / "datasets" / "derks2022_holdout.h5ad"
a.write_h5ad(out, compression="gzip")
Xa = np.asarray(a.X)
print("shape", a.shape, f"missing={np.isnan(Xa).mean():.3f}", f"range=[{np.nanmin(Xa):.2f}, {np.nanmax(Xa):.2f}]")
print("cell_type", a.obs.cell_type.value_counts().to_dict(), "instrument", a.obs.instrument.value_counts().to_dict())
print("sha256", hashlib.sha256(out.read_bytes()).hexdigest())
