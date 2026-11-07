"""Shared loader for the PXD061710 (Spectronaut) holdout builders. Single cells are selected by joining the result
columns to PXD061710_community_annotated.sdrf.tsv (cells per well == 1); 2/5/10/50-cell and blank runs are dropped."""
import re

import numpy as np
import pandas as pd

from _holdout_common import save_holdout, to_log2
from scp_graphfs_mcp.ids import group_accessions


def build(name, tsv, sdrf_path, quant_suffix, label_map, n_expected, sha, uns_extra):
    d = pd.read_csv(tsv, sep="\t", low_memory=False)
    s = pd.read_csv(sdrf_path, sep="\t")
    s["run"] = s["comment[data file]"].str.replace(r"\.d$", "", regex=True)
    s["short"] = s["run"].str.replace(r"_[sS]\d-[a-pA-P]\d+_\d+_\d+$", "", regex=True)
    qcols, runs = [], []
    for c in d.columns:
        if not c.endswith(quant_suffix):
            continue
        r = re.sub(r"^\[\d+\] ", "", c)[: -len(quant_suffix)]
        r = r[:-2] if r.endswith(".d") else r
        qcols.append(c); runs.append(r)
    key = "run" if all(r in set(s["run"]) for r in runs if "_SC_" in r) else "short"
    s = s.drop_duplicates(key).set_index(key)
    keep = [i for i, r in enumerate(runs) if r in s.index and int(s.loc[r, "characteristics[cells per well]"]) == 1]
    assert all("_SC_" in runs[i] for i in keep) and len(keep) == n_expected, (len(keep), n_expected)
    dropped = sorted(r for i, r in enumerate(runs) if i not in set(keep))
    assert all(not re.search(r"_SC_", r) for r in dropped), dropped  # only blanks / multi-cell runs dropped
    runs = [runs[i] for i in keep]; qcols = [qcols[i] for i in keep]
    meta = s.loc[runs]
    labels = [label_map[t] for t in meta["characteristics[treatment]"]]
    X = to_log2(d[qcols].apply(pd.to_numeric, errors="coerce").to_numpy(float).T)
    accs = [(group_accessions(str(g)) or [""])[0] for g in d["PG.ProteinGroups"]]
    cells = meta["run"].to_numpy() if key == "short" else np.array(runs)
    print(f"  dropped non-single-cell runs: {dropped}")
    save_holdout(name, X, accs, cells, labels,
                 extra_obs={"cell_line": meta["characteristics[cell line]"].to_numpy(),
                            "treatment": meta["characteristics[treatment]"].to_numpy()},
                 uns={"source": f"PXD061710 {tsv.name} ({quant_suffix.strip('.')}, Spectronaut)", "sha256": sha,
                      **uns_extra})
