"""HOLDOUT v9 karagach2025_S1 -> data/cache/holdout/karagach2025_S1.npz + data/datasets/karagach2025_S1_holdout.h5ad.

Source: PRIDE PXD058457 (Karagach, Smollich, ..., Geiger 2025, "High throughput single-cell proteomics of in-vivo
cells", Mol Cell Proteomics, doi:10.1016/j.mcpro.2025.101018; licence CC0). HeLa single cells, CellenONE, label-free
Slice-PASEF DIA on timsTOF SCP (Evosep ZOOM40), DIA-NN 1.9.1.
File: https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD058457/S1_rehyd_cover_384.pg_matrix.tsv
(SHA-256 60715ab23c38a24bd0746fdd94dfb1d9c71082a42622b652048dcffb07641b1a; 283 runs). SDRF
Karagach_Smollich_Resubmision.sdrf.tsv (same directory) annotates every run as HeLa (assays "S1 rehydration 50
percent" and "S1_rehydration vs covered"); every run name carries `HeLa_SC_` / `Hela_SC_`. No blanks/bulk in the file.
Labels (sample-prep variants, technical not biological): master-mix volume + rehydration volume (`MM<v>_rehyd<v>`) or
plate covered/uncovered (`MM<v>_cover` / `MM<v>_open`), parsed from the run name.
Processing: _holdout_common.save_holdout (corpus save() logic).

    python scripts/holdout_ext/build_karagach2025_S1.py [data/cache/corpus_raw/karagach2025]
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _holdout_common import RAW, save_holdout, to_log2
from scp_graphfs_mcp.ids import group_accessions

src = Path(sys.argv[1]) if len(sys.argv) > 1 else RAW / "karagach2025"
d = pd.read_csv(src / "S1_rehyd_cover_384.pg_matrix.tsv", sep="\t")
sdrf = pd.read_csv(src / "Karagach_Smollich_Resubmision.sdrf.tsv", sep="\t")
sdrf.index = sdrf["comment[data file]"].str.strip()
qcols = list(d.columns[5:])
runs = [re.sub(r".*[\\/]", "", c) for c in qcols]
labels = []
for r in runs:
    assert sdrf.loc[r, "characteristics[cell type]"] == "HeLa", r
    assert re.match(r"^he[lL]a_SC_", r, re.I), r
    m = re.search(r"_(\d+)nlMM_rehyd_(\d+)nl_", r)
    if m:
        labels.append(f"MM{m.group(1)}_rehyd{m.group(2)}"); continue
    m = re.search(r"_(\d+)nlMM_(cover_)?Bravo2_", r)
    assert m, r
    labels.append(f"MM{m.group(1)}_{'cover' if m.group(2) else 'open'}")
assert len(runs) == 283
X = to_log2(d[qcols].to_numpy(float).T)
accs = [(group_accessions(str(g)) or [""])[0] for g in d["Protein.Group"]]
save_holdout("karagach2025_S1", X, accs, [r[:-2] for r in runs], labels,
             extra_obs={"assay": sdrf.loc[runs, "assay name"].str.strip().to_numpy()},
             uns={"source": "PXD058457 S1_rehyd_cover_384.pg_matrix.tsv (DIA-NN 1.9.1)",
                  "sha256": "60715ab23c38a24bd0746fdd94dfb1d9c71082a42622b652048dcffb07641b1a"})
