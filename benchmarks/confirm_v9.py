"""v9 holdout evaluation (PREREGISTRATION_v9.md): the H15 champion, its Phase-C revision and the full comparator panel
on holdouts that no method in this repository has been developed on.

    python benchmarks/confirm_v9.py <holdout> [--jobs 2]   -> results/confirm_v9/<holdout>.json   (single use each)
    python benchmarks/confirm_v9.py --dry-run-dev           -> scratch output on DEV target weiss2026_liver, 2 seeds

Fixed by the preregistration: holdout list, seeds 9000-9009, 30-cell random subsamples (no stratification), MCAR and
MNAR masks, the comparator set, and the dev-selected "best panel" and "best graph-prior" competitors, which are read
from results/v9/v9_selection.json (committed with the preregistration). Reference = every corpus dataset in
explore_phaseb.reference_datasets() (no holdout of any version, no non-human data), aligned to the holdout's proteins.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from explore_phaseb import reference_datasets  # noqa: E402
from panel_v9 import evaluate_split  # noqa: E402
from run import SIG_TOL, boot_rmse_diff  # noqa: E402

from scp_graphfs_mcp import baselines as B  # noqa: E402
from scp_graphfs_mcp import transfer as T  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SEEDS = list(range(9000, 9010))
PRIMARY = ["derks2022", "karagach2025_S1", "px061710_GSTP1KO", "px061710_MDAfix", "px061710_RKOfix"]
SECONDARY = ["khan2023", "brunner2022", "karagach2025_S2", "px077946_H2O2"]
H5AD = {"khan2023": "khan2023_proteins_unimputed.h5ad"}
SELECTION = ROOT / "results" / "v9" / "v9_selection.json"


def load_holdout(name):
    import anndata as ad

    a = ad.read_h5ad(ROOT / "data" / "datasets" / H5AD.get(name, f"{name}_holdout.h5ad"))
    return np.asarray(a.X, float), [str(v) for v in a.var_names]


def _split(args):
    X, names, V, warm, Vnull, kind, seed, target = args
    se, sel, wts, secs = evaluate_split(X, names, None, V, warm, kind, seed, target, list(B.ALL), Vnull=Vnull)
    print(f"[{time.strftime('%H:%M:%S')}] {target} {kind} {seed} ({sum(secs.values()):.0f}s)", flush=True)
    return kind, seed, se, sel, wts, secs


def tests(SE, SEL, a, b):
    out = []
    for endpoint, pick in (("all", None), ("in_complex", SEL)):
        sa = SE[a] if pick is None else SE[a][pick]
        sb = SE[b] if pick is None else SE[b][pick]
        c = boot_rmse_diff(sa, [sb], 1000, seed=17)
        ref = float(np.sqrt(sb.mean()))
        pct = 100 * c["delta_rmse"] / ref
        out.append({"method": a, "vs": b, "endpoint": endpoint, **c, "pct_of_ref_rmse": pct,
                    "pass": bool(c["ci95"][1] < -SIG_TOL and pct <= -1.0)})
    return out


def main():
    a = sys.argv
    dry = "--dry-run-dev" in a
    jobs = int(a[a.index("--jobs") + 1]) if "--jobs" in a else 2
    sel = json.loads(SELECTION.read_text())
    best_panel, best_graph = sel["best_panel"], sel["best_graph"]
    if dry:
        from explore_phaseb import load_target

        target, seeds = "weiss2026_liver", [100, 101]
        X, names, _ = load_target(target)
        out_path = Path("/tmp/claude-0/-home-user-knowledgeispower/19f0adc1-0e52-561e-b422-c31fec9b93e5/scratchpad/"
                        "confirm_v9_dryrun.json")
        refs = reference_datasets(exclude=(target,))
    else:
        target = a[1]
        if target not in PRIMARY + SECONDARY:
            sys.exit(f"refusing: {target} is not a v9 holdout")
        out_path = ROOT / "results" / "confirm_v9" / f"{target}.json"
        if out_path.exists():
            sys.exit(f"refusing: {out_path} exists; each v9 holdout is single-use")
        if subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "src", "benchmarks", "results/v9",
                           "PREREGISTRATION_v9.md", "benchmarks/PREREGISTRATION_v9.md"], cwd=ROOT).returncode != 0:
            sys.exit("refusing: uncommitted changes in src/, benchmarks/ or the v9 selection")
        if subprocess.run(["git", "ls-files", "--error-unmatch", str(SELECTION.relative_to(ROOT)),
                           "benchmarks/PREREGISTRATION_v9.md"], cwd=ROOT, capture_output=True).returncode != 0:
            sys.exit("refusing: the preregistration and selection file must be committed first")
        seeds = SEEDS
        X, names = load_holdout(target)
        refs = reference_datasets()
    assert not set(refs) & (T.PHASE_B_HOLDOUTS | T.HOLDOUTS | T.V9_HOLDOUTS)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    t0 = time.time()
    Z = T.corpus_reference(names, (target,) if dry else (), refs)
    V, warm = T.fit_loadings(Z)
    Vnull = T.permute_rows(V, warm, 0)
    jl = [(X, names, V, warm, Vnull, k, s, target) for k in ("MCAR", "MNAR") for s in seeds]
    with ProcessPoolExecutor(jobs) as ex:
        res = list(ex.map(_split, jl))
    out = {"code_commit": commit, "dataset": target + (" (DRY RUN)" if dry else ""), "seeds": seeds,
           "shape": list(X.shape), "missing_fraction": float(np.isnan(X).mean()),
           "reference": {"datasets": refs, "shape": list(Z.shape), "warm_fraction": float(warm.mean())},
           "best_panel": best_panel, "best_graph": best_graph, "masks": {}}
    for kind in ("MCAR", "MNAR"):
        rr = [r for r in res if r[0] == kind]
        SE = {m: np.concatenate([r[2][m] for r in rr]) for m in rr[0][2]}
        SEL = np.concatenate([r[3] for r in rr])
        rmse = {m: float(np.sqrt(v.mean())) for m, v in SE.items()}
        rmse_cx = {m: float(np.sqrt(v[SEL].mean())) for m, v in SE.items()}
        primary = []
        for comp in ("stack_lowrank", "champion@null", "champion@noref", best_panel, best_graph):
            primary += tests(SE, SEL, "champion", comp)
        revision = tests(SE, SEL, "revised", "champion")
        revised_vs = [t for c in ("stack_lowrank", best_panel, best_graph) for t in tests(SE, SEL, "revised", c)]
        gcn_controls = [t for c in ("gcn_null", "gcn_nograph") for t in tests(SE, SEL, "gcn", c)]
        descriptive = {m: tests(SE, SEL, "champion", m) for m in SE if m != "champion"}
        out["masks"][kind] = {"rmse": rmse, "rmse_in_complex": rmse_cx, "primary": primary, "revision": revision,
                              "revised_vs_competitors": revised_vs, "gcn_string_vs_controls": gcn_controls,
                              "champion_vs_every_method": descriptive,
                              "stack_weights": {m: [r[4].get(m) for r in rr] for m in ("champion", "revised")},
                              "seconds": [r[5] for r in rr]}
    pp = [t["pass"] for k in out["masks"].values() for t in k["primary"]]
    rp = [t["pass"] for k in out["masks"].values() for t in k["revision"]]
    hp = [t["pass"] for k in out["masks"].values() for t in k["revised_vs_competitors"]]
    out["H18_confirmed"], out["H18_passes"] = bool(all(hp)), f"{sum(hp)}/{len(hp)}"
    out["H16_confirmed"], out["H16_passes"] = bool(all(pp)), f"{sum(pp)}/{len(pp)}"
    out["H17_confirmed"], out["H17_passes"] = bool(all(rp)), f"{sum(rp)}/{len(rp)}"
    out["runtime_seconds"] = round(time.time() - t0, 1)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("dataset", "H16_passes", "H17_passes", "H18_passes", "runtime_seconds")}))
    for kind, v in out["masks"].items():
        for t in v["primary"] + v["revision"] + v["revised_vs_competitors"]:
            print(kind, t["method"], "vs", t["vs"], t["endpoint"], f"{t['pct_of_ref_rmse']:+.2f}%",
                  [round(x, 4) for x in t["ci95"]], "PASS" if t["pass"] else "fail")


if __name__ == "__main__":
    main()
