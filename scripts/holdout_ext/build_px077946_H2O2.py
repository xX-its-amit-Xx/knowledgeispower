"""HOLDOUT v9 px077946_H2O2 -> data/cache/holdout/px077946_H2O2.npz + data/datasets/px077946_H2O2_holdout.h5ad.

Source: PRIDE PXD077946 ("Refining liquid chromatography conditions across platforms to democratize single-cell
proteomics"; Wang ... Haiyuan Yu; CC0). PRIDE lists no reference; by title match the preprint is bioRxiv
doi:10.64898/2025.12.16.694739. K562, HEK293T and HCT116 single cells, untreated vs 300 uM H2O2 for 30 min, Evosep +
label-free dia-PASEF on timsTOF HT ("WhisperZoom"), DIA-NN.
File: https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD077946/diann_singlecell.zip
(zip SHA-256 0bd4c3b9bbd3ad8e94acf7e3ef9dba4330dd1c25667630f177539b5bc5714a05; members used:
report.unique_genes_matrix.tsv aeb2df163cedc56543f02ada9b99d6bfcda5edf3f36c7fbe418c72351c33d197,
report.pr_matrix.tsv 674317e5472390f1a9485ff8bd4bdcd46f5126eb3f3524f7d177126ab55c244d). FEWER THAN 60 CELLS (40).
The archive ships no pg_matrix. Protein values = DIA-NN's unique-genes matrix (MaxLFQ over proteotypic precursors,
linear, NaN = missing); each gene is mapped to the first canonical UniProt accession of the Protein.Group carried by its
proteotypic precursors in report.pr_matrix.tsv (most frequent group; 1 gene, TMPO, has two groups).
All 40 runs come from the archive named "singlecell" and the run directory `...h2o2sc`; there are no blank/bulk runs
(the 3 ng runs are in a separate archive, diann_3ng.zip, not used). Labels: <cell line>_<ctrl|H2O2>.

    python scripts/holdout_ext/build_px077946_H2O2.py [data/cache/corpus_raw/px077946/x]
"""
import re
import sys
from pathlib import Path

import pandas as pd

from _holdout_common import RAW, save_holdout, to_log2
from scp_graphfs_mcp.ids import group_accessions

src = Path(sys.argv[1]) if len(sys.argv) > 1 else RAW / "px077946" / "x"
g = pd.read_csv(src / "report.unique_genes_matrix.tsv", sep="\t")
p = pd.read_csv(src / "report.pr_matrix.tsv", sep="\t", usecols=["Protein.Group", "Genes", "Proteotypic"])
p = p[p["Proteotypic"] == 1]
cnt = p.groupby(["Genes", "Protein.Group"]).size().rename("n").reset_index()
cnt = cnt.sort_values(["Genes", "n", "Protein.Group"], ascending=[True, False, True], kind="stable")
pg_of = cnt.drop_duplicates("Genes").set_index("Genes")["Protein.Group"]
qcols = list(g.columns[1:])
runs = [re.sub(r".*[\\/]", "", c)[:-2] for c in qcols]
labels = []
for r in runs:
    m = re.search(r"_1p5_(293T|HCT116|K562)_(300UMH2O2_)?S5-[A-H]\d+_1_\d+$", r)
    assert m, r
    labels.append(f"{m.group(1)}_{'H2O2' if m.group(2) else 'ctrl'}")
assert len(runs) == 40
X = to_log2(g[qcols].to_numpy(float).T)
accs = [(group_accessions(str(pg_of.get(gene, ""))) or [""])[0] for gene in g["Genes"]]
save_holdout("px077946_H2O2", X, accs, runs, labels,
             uns={"source": "PXD077946 diann_singlecell.zip: report.unique_genes_matrix.tsv (DIA-NN)",
                  "sha256_zip": "0bd4c3b9bbd3ad8e94acf7e3ef9dba4330dd1c25667630f177539b5bc5714a05",
                  "instrument": "timsTOF HT, dia-PASEF"})
