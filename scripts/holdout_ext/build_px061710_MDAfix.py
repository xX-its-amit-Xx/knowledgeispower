"""HOLDOUT v9 px061710_MDAfix -> data/cache/holdout/px061710_MDAfix.npz + data/datasets/px061710_MDAfix_holdout.h5ad.

Source: PRIDE PXD061710 (Makar ... von Kriegsheim; CC0; paper by title match J Proteome Res 2025,
doi:10.1021/acs.jproteome.5c00268). MDA-MB-231 cells fresh vs 0.1 / 0.3 / 1% PFA fixed for 24 h (Figure 3), CellenONE,
label-free diaPASEF on timsTOF SCP, Spectronaut.
File: https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD061710/20240709MDA_fresh_fixed_24hrs_Figure3_search_engine_results.tsv
(SHA-256 72d750b7847a0af2978019f37e1002251e23d3f9fa2efd4e1c31cfc3d434e5fb). Quantity `PG.Quantity` (linear). Single
cells via SDRF (cells per well 1; 4 x 19 = 76); blanks and 2/5/10/50-cell runs dropped.
Labels: MDA_fresh, MDA_PFA0.1_24h, MDA_PFA0.3_24h, MDA_PFA1_24h.

    python scripts/holdout_ext/build_px061710_MDAfix.py [data/cache/corpus_raw/px061710]
"""
import sys
from pathlib import Path

from _holdout_common import RAW
from _px061710 import build

src = Path(sys.argv[1]) if len(sys.argv) > 1 else RAW / "px061710"
build("px061710_MDAfix", src / "20240709MDA_fresh_fixed_24hrs_Figure3_search_engine_results.tsv",
      src / "PXD061710_community_annotated.sdrf.tsv", ".PG.Quantity",
      {"fresh, no fixation": "MDA_fresh", "0.1% paraformaldehyde fixation, 24 h": "MDA_PFA0.1_24h",
       "0.3% paraformaldehyde fixation, 24 h": "MDA_PFA0.3_24h", "1% paraformaldehyde fixation, 24 h": "MDA_PFA1_24h"},
      76, "72d750b7847a0af2978019f37e1002251e23d3f9fa2efd4e1c31cfc3d434e5fb", {"instrument": "timsTOF SCP, diaPASEF"})
