"""HOLDOUT v9 px061710_GSTP1KO -> data/cache/holdout/px061710_GSTP1KO.npz + data/datasets/px061710_GSTP1KO_holdout.h5ad.

Source: PRIDE PXD061710 ("Overcoming preservation challenges to enable single-cell proteomics of fixed cell and tissue
samples with retained proteome integrity"; Makar ... von Kriegsheim; CC0). PRIDE lists no reference; by title match the
paper is J Proteome Res 2025, doi:10.1021/acs.jproteome.5c00268 (preprint doi:10.1101/2025.03.10.642380).
RKO cells, sgControl vs sgGSTP1 knockout (Figure 1), CellenONE, label-free diaPASEF on timsTOF SCP, Spectronaut.
File: https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD061710/20230613_GSTP1KO_Figure1_search_engine_results.tsv
(SHA-256 740e3c618cf5897c3099b71cffb08995806aa32a7c8874da0c34c49f9214abe7). Quantity: `PG.MS2Quantity` (linear; the
file also has PG.MS1Quantity, not used). Single cells via SDRF PXD061710_community_annotated.sdrf.tsv (cells per well 1,
78 = 39 + 39); blanks and 2/5/10/50-cell runs dropped. Labels: RKO_sgCTRL / RKO_sgGSTP1.

    python scripts/holdout_ext/build_px061710_GSTP1KO.py [data/cache/corpus_raw/px061710]
"""
import sys
from pathlib import Path

from _holdout_common import RAW
from _px061710 import build

src = Path(sys.argv[1]) if len(sys.argv) > 1 else RAW / "px061710"
build("px061710_GSTP1KO", src / "20230613_GSTP1KO_Figure1_search_engine_results.tsv",
      src / "PXD061710_community_annotated.sdrf.tsv", ".PG.MS2Quantity",
      {"sgControl (non-targeting sgRNA)": "RKO_sgCTRL", "sgGSTP1 knockout": "RKO_sgGSTP1"}, 78,
      "740e3c618cf5897c3099b71cffb08995806aa32a7c8874da0c34c49f9214abe7", {"instrument": "timsTOF SCP, diaPASEF"})
