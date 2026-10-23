"""Single-use confirmation of H15 (adapted transfer) on the wu2026 holdout. PREREGISTRATION_v8.md.

    python benchmarks/confirm_v8.py                -> results/confirm_v8/confirm_v8.json   (refuses to run twice)
    python benchmarks/confirm_v8.py --dry-run-dev  -> scratch output on a DEV target (weiss2026_liver, 2 seeds)
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from explore_phaseb import GRID, reference_datasets  # noqa: E402
from run import SIG_TOL, boot_rmse_diff, make_mask  # noqa: E402

from scp_graphfs_mcp import phase2 as P  # noqa: E402
from scp_graphfs_mcp import transfer as T  # noqa: E402
from scp_graphfs_mcp.complexes import in_any_complex  # noqa: E402
from scp_graphfs_mcp.external import CORPUS_DIR  # noqa: E402
from scp_graphfs_mcp.ids import group_accessions  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
LV, EXTRA = 100.0, 10


def main():
    dry = "--dry-run-dev" in sys.argv
    target = "weiss2026_liver" if dry else "wu2026"
    if dry:
        out_path = Path("/tmp/claude-0/-home-user-knowledgeispower/19f0adc1-0e52-561e-b422-c31fec9b93e5/scratchpad/confirm_v8_dryrun.json")
        seeds = [100, 101]
    else:
        out_path = ROOT / "results" / "confirm_v8" / "confirm_v8.json"
        if out_path.exists():
            sys.exit(f"refusing: {out_path} exists; the v8 confirmation is single-use")
        if subprocess.run(["git", "diff", "--quiet", "--", "src", "benchmarks"], cwd=ROOT).returncode != 0:
            sys.exit("refusing: uncommitted changes in src/ or benchmarks/")
        seeds = list(range(8000, 8010))
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    t0 = time.time()
    d = np.load(CORPUS_DIR / f"{target}.npz")
    X, names = d["X"].astype(float), [str(a) for a in d["accs"]]
    refs = [r for r in reference_datasets(exclude=(target,)) if r != target]
    assert "wu2026" not in refs and "derks2022" not in refs and "brunner2022" not in refs
    Z = T.corpus_reference(names, (target,) if target not in T.PHASE_B_HOLDOUTS else (), refs)
    V, warm = T.fit_loadings(Z)
    Vn = T.permute_rows(V, warm, 0)
    in_cx = np.array([in_any_complex(group_accessions(f)) for f in names])
    out = {"code_commit": commit, "dataset": target + (" (DRY RUN)" if dry else ""), "seeds": seeds,
           "reference": {"datasets": refs, "shape": list(Z.shape), "warm_fraction": float(warm.mean())},
           "params": {"lam_v": LV, "extra": EXTRA}, "masks": {}}
    for kind in ["MCAR", "MNAR"]:
        se, sels = {}, []
        for seed in seeds:
            idx = np.sort(np.random.default_rng(seed).choice(X.shape[0], 30, replace=False))
            Xs = X[idx]
            eligible = (~np.isnan(Xs)).mean(0) >= 0.5
            mask = make_mask(Xs, eligible, kind, seed)
            Xm = Xs.copy()
            Xm[mask] = np.nan
            truth = Xs[mask]
            sels.append(in_cx[np.nonzero(mask)[1]])
            r2 = np.random.default_rng(seed + 7919)
            cand = np.argwhere(~np.isnan(Xm) & eligible[None, :])
            pick = cand[r2.choice(len(cand), max(1, int(0.1 * len(cand))), replace=False)]
            Xi = Xm.copy()
            Xi[pick[:, 0], pick[:, 1]] = np.nan
            ti = Xm[pick[:, 0], pick[:, 1]]
            pk = lambda M: M[pick[:, 0], pick[:, 1]]  # noqa: E731

            def stack(fulls, inners):
                g = GRID[len(fulls)]
                errs = [np.nanmean((sum(w * pk(p) for w, p in zip(ws, inners)) - ti) ** 2) for ws in g]
                ws = g[int(np.nanargmin(errs))]
                return sum(w * p for w, p in zip(ws, fulls))

            base, bi = P._fill(Xm, P.loo_knn(Xm)), P._fill(Xi, P.loo_knn(Xi))
            base, bi = np.where(np.isnan(base), np.nanmean(Xm), base), np.where(np.isnan(bi), np.nanmean(Xi), bi)
            lr, lri = P.grmf(Xm, None, names, lam_g=0.0), P.grmf(Xi, None, names, lam_g=0.0)

            def adapted(VV, ww):
                return (T.adapt_project(Xm, VV, ww, base, lam_v=LV, extra=EXTRA),
                        T.adapt_project(Xi, VV, ww, bi, lam_v=LV, extra=EXTRA))

            preds = {"knn": base, "stack_lowrank": stack([base, lr], [bi, lri])}
            for key, (VV, ww) in {"adapt": (V, warm), "adapt@null": (Vn, warm),
                                  "adapt@noref": (np.zeros_like(V), np.zeros_like(warm))}.items():
                a, ai = adapted(VV, ww)
                preds[key] = stack([base, lr, a], [bi, lri, ai])
            for m, Y in preds.items():
                se.setdefault(m, []).append((Y[mask] - truth) ** 2)
            print(f"[{time.strftime('%H:%M:%S')}] {kind} seed {seed} ({time.time() - t0:.0f}s)", flush=True)
        SE = {m: np.concatenate(v) for m, v in se.items()}
        sel = np.concatenate(sels)
        res = {"rmse_all": {m: float(np.sqrt(v.mean())) for m, v in SE.items()},
               "rmse_in_complex": {m: float(np.sqrt(v[sel].mean())) for m, v in SE.items()}, "comparisons": []}
        for ep in ("all", "in_complex"):
            f = (lambda v: v) if ep == "all" else (lambda v: v[sel])
            for comp in ("knn", "stack_lowrank", "adapt@null", "adapt@noref"):
                c = boot_rmse_diff(f(SE["adapt"]), [f(SE[comp])], 1000, seed=17)
                pct = 100 * c["delta_rmse"] / float(np.sqrt(f(SE[comp]).mean()))
                res["comparisons"].append({"vs": comp, "endpoint": ep, **c, "pct_of_ref_rmse": pct,
                                           "pass": bool(c["ci95"][1] < -SIG_TOL and pct <= -1.0)})
        out["masks"][kind] = res
    passes = [c["pass"] for k in out["masks"].values() for c in k["comparisons"]]
    out["H15_confirmed"] = bool(all(passes))
    out["passes"] = f"{sum(passes)}/{len(passes)}"
    out["runtime_seconds"] = round(time.time() - t0, 1)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ["H15_confirmed", "passes", "runtime_seconds"]}))
    for k, v in out["masks"].items():
        for c in v["comparisons"]:
            print(k, c["endpoint"], "adapt vs", c["vs"], f"{c['pct_of_ref_rmse']:+.2f}%", [round(x, 4) for x in c["ci95"]],
                  "PASS" if c["pass"] else "fail")


if __name__ == "__main__":
    main()
