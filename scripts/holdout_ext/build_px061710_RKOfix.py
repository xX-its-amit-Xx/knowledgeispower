"""HOLDOUT v9 px061710_RKOfix -> data/cache/holdout/px061710_RKOfix.npz + data/datasets/px061710_RKOfix_holdout.h5ad.

Source: PRIDE PXD061710 (Makar ... von Kriegsheim; CC0; paper by title match J Proteome Res 2025,
doi:10.1021/acs.jproteome.5c00268). RKO cells fresh vs 0.1 / 0.3 / 1% PFA fixed for 20 min (Figure 4), CellenONE,
label-free diaPASEF on timsTOF SCP, Spectronaut.
File: https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD061710/20240809RKOFreshFixed20min_Figure4_search_engine_results.tsv
(SHA-256 2bb0dd829c3137c2fa54c7d1d843169229aeff9ea9146d5764b989759b36a03d). Quantity `PG.Quantity` (linear). Single
cells via SDRF (cells per well 1; 19 + 19 + 19 + 18 = 75); blanks and 2/5/10/50-cell runs dropped.
Labels: RKO_fresh, RKO_PFA0.1_20min, RKO_PFA0.3_20min, RKO_PFA1_20min.

    python scripts/holdout_ext/build_px061710_RKOfix.py [data/cache/corpus_raw/px061710]
"""
import sys
from pathlib import Path

from _holdout_common import RAW
from _px061710 import build

src = Path(sys.argv[1]) if len(sys.argv) > 1 else RAW / "px061710"
build("px061710_RKOfix", src / "20240809RKOFreshFixed20min_Figure4_search_engine_results.tsv",
      src / "PXD061710_community_annotated.sdrf.tsv", ".PG.Quantity",
      {"fresh, no fixation": "RKO_fresh", "0.1% paraformaldehyde fixation, 20 min": "RKO_PFA0.1_20min",
       "0.3% paraformaldehyde fixation, 20 min": "RKO_PFA0.3_20min", "1% paraformaldehyde fixation, 20 min": "RKO_PFA1_20min"},
      75, "2bb0dd829c3137c2fa54c7d1d843169229aeff9ea9146d5764b989759b36a03d", {"instrument": "timsTOF SCP, diaPASEF"})
