"""HOLDOUT v9 karagach2025_S2 -> data/cache/holdout/karagach2025_S2.npz + data/datasets/karagach2025_S2_holdout.h5ad.

Source: PRIDE PXD058457 (Karagach et al. 2025, Mol Cell Proteomics, doi:10.1016/j.mcpro.2025.101018; CC0). HeLa single
cells, automated (Bravo) prep, label-free Slice-PASEF DIA on timsTOF SCP (ZOOM40), DIA-NN 1.9.1.
File: https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD058457/S2_Automation.pg_matrix.tsv
(SHA-256 4cdbf497294d88762ebf3b971bca088801c2bd553be858971a3b5377043ab5be; 40 runs, all SDRF assay "S2_Automation",
cell type HeLa, names `Hela_SC_`). FEWER THAN 60 CELLS.
Labels (acquisition variants, technical): number of Slice-PASEF frames (`1F`, `2F`, `4F`, `2F_Demichev`) from run name.

    python scripts/holdout_ext/build_karagach2025_S2.py [data/cache/corpus_raw/karagach2025]
"""
import re
import sys
from pathlib import Path

import pandas as pd

from _holdout_common import RAW, save_holdout, to_log2
from scp_graphfs_mcp.ids import group_accessions

src = Path(sys.argv[1]) if len(sys.argv) > 1 else RAW / "karagach2025"
d = pd.read_csv(src / "S2_Automation.pg_matrix.tsv", sep="\t")
sdrf = pd.read_csv(src / "Karagach_Smollich_Resubmision.sdrf.tsv", sep="\t")
sdrf.index = sdrf["comment[data file]"].str.strip()
qcols = list(d.columns[5:])
runs = [re.sub(r".*[\\/]", "", c) for c in qcols]
labels = []
for r in runs:
    assert sdrf.loc[r, "characteristics[cell type]"] == "HeLa" and r.startswith("Hela_SC_"), r
    m = re.search(r"_slice_pasef_(\dF(?:_Demichev)?)_", r)
    assert m, r
    labels.append(m.group(1))
assert len(runs) == 40
X = to_log2(d[qcols].to_numpy(float).T)
accs = [(group_accessions(str(g)) or [""])[0] for g in d["Protein.Group"]]
save_holdout("karagach2025_S2", X, accs, [r[:-2] for r in runs], labels,
             uns={"source": "PXD058457 S2_Automation.pg_matrix.tsv (DIA-NN 1.9.1)",
                  "sha256": "4cdbf497294d88762ebf3b971bca088801c2bd553be858971a3b5377043ab5be"})
