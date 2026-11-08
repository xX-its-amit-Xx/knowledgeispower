"""POST HOC: the revision vs the R reference imputers (and the other added comparators) on exactly matched splits.

Addendum 2 (R references) ran 4 of the 10 preregistered seeds and did not score `revised`; the addendum (v3) scored `revised` on
all 10. Both cache per-entry squared errors per split, on identical masks (checked: the champion's errors are identical in both).
This joins them on seeds 9000-9003 and writes results/v9/addendum2/matched_vs_revised.json with RMSE ratios to `revised`, a
per-split t-interval (4 splits: indicative only) and the per-split ratios.

    .venv/bin/python scripts/matched_vs_revised_v9.py
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
C1 = ROOT / "data" / "cache" / "v9_addendum" / "v3"
C2 = ROOT / "data" / "cache" / "v9_addendum2"
OUT = ROOT / "results" / "v9" / "addendum2" / "matched_vs_revised.json"
HOLD = ["derks2022", "karagach2025_S1", "px061710_GSTP1KO", "px061710_MDAfix", "px061710_RKOfix",
        "khan2023", "brunner2022", "karagach2025_S2", "px077946_H2O2"]
SEEDS = [9000, 9001, 9002, 9003]
FROM1 = ["champion", "revised", "bpca_oba", "bpca_oba+stack", "gcn_nograph+stack", "gcn+stack", "pimms_vae+stack"]
FROM2 = ["r_bpca", "r_bpca+stack", "r_missforest", "r_missforest+stack", "rem_protein", "seqknn+stack"]


def main():
    out = {"note": "post hoc; 4 of the 10 preregistered seeds; matched splits", "seeds": SEEDS, "holdouts": {}}
    for h in HOLD:
        res = {}
        for kind in ("MCAR", "MNAR"):
            per = {m: [] for m in FROM1 + FROM2}
            same = True
            for s in SEEDS:
                f1 = C1 / f"holdout_{h}_{kind}_{s}.npz"
                f2 = glob.glob(str(C2 / f"{h}_{kind}_{s}_*.npz"))
                if not f1.exists() or not f2:
                    per = None
                    break
                a = np.load(f1, allow_pickle=True)["se"].item()
                b = np.load(f2[0], allow_pickle=True)["se"].item()
                same &= bool(np.allclose(a["champion"], b["champion"]))
                for m in FROM1:
                    per[m].append(a[m])
                for m in FROM2:
                    per[m].append(b[m])
            if per is None:
                continue
            rv = np.concatenate(per["revised"])
            row = {"champion_identical_in_both_caches": same, "methods": {}}
            for m, v in per.items():
                ratio = float(np.sqrt(np.concatenate(v).mean() / rv.mean()))
                ps = np.array([np.sqrt(x.mean() / y.mean()) for x, y in zip(v, per["revised"])])
                half = stats.t.ppf(0.975, len(ps) - 1) * ps.std(ddof=1) / np.sqrt(len(ps)) if len(ps) > 1 else np.nan
                row["methods"][m] = {"rmse_ratio_to_revised": ratio, "per_split": ps.round(4).tolist(),
                                     "split_t95": [float(ps.mean() - half), float(ps.mean() + half)]}
            res[kind] = row
        if res:
            out["holdouts"][h] = res
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    for h, r in out["holdouts"].items():
        for kind, row in r.items():
            print(h, kind, row["champion_identical_in_both_caches"],
                  {m: round(v["rmse_ratio_to_revised"], 3) for m, v in row["methods"].items()
                   if m in ("champion", "r_bpca", "r_missforest", "bpca_oba", "gcn_nograph+stack")})


if __name__ == "__main__":
    main()
