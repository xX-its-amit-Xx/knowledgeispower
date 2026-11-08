"""v9 ADDENDUM (post hoc, exploratory; no preregistered verdict depends on it). Written 2026-11-08 after the
independent review (docs/v9_review.md). For a holdout whose preregistered result already exists, or for a dev target,
it recomputes the same splits (same seeds, masks, reference) and adds:

  1. exact-reproduction check of the frozen run's RMSEs (champion, revised, controls, comparators);
  2. the addendum comparators (src/scp_graphfs_mcp/baselines_addendum.py: correctly oriented BPCA, LLS, kNN over
     proteins, QRILC, Perseus down-shift), each raw and stacked;
  3. robustness of the H16/H17/H18 tests: per-entry errors are kept with (split, cell, protein) coordinates, so the CI is
     recomputed with (a) the preregistered entry bootstrap, (b) a protein-cluster bootstrap, (c) a split-level
     bootstrap, and each test is re-scored under the strict SESOI reading (CI upper bound <= -1%).

    python benchmarks/addendum_v9.py <holdout> [--jobs 1]    -> results/v9/addendum/<holdout>.json
    python benchmarks/addendum_v9.py --dev <target>          -> results/v9/addendum/dev_<target>.json (seeds 100-101)

Heavy panel methods (missForest, MICE, BPCA-old, ...) are not recomputed: their frozen numbers stay in
results/confirm_v9/. Per-split results are cached under data/cache/v9_addendum/ so a restart resumes.
"""

from __future__ import annotations

import json
import os
import sys
import time
import warnings

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
from concurrent.futures import ProcessPoolExecutor  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from explore_phaseb import reference_datasets, subsample  # noqa: E402
from run import SIG_TOL, make_mask  # noqa: E402

from scp_graphfs_mcp import baselines as B  # noqa: E402
from scp_graphfs_mcp import baselines_addendum as A  # noqa: E402
from scp_graphfs_mcp import transfer as T  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache" / "v9_addendum"
OUT = ROOT / "results" / "v9" / "addendum"
PANEL = ["gcn_nograph", "gcn"] + list(A.ADDENDUM)  # recomputed here; everything else stays frozen
COMPS16 = ["stack_lowrank", "champion@null", "champion@noref", "gcn_nograph+stack", "gcn+stack"]
COMPS18 = ["stack_lowrank", "gcn_nograph+stack", "gcn+stack"]


def _coords(X, labels, kind, seed, n_cells=30):
    """The exact masked-entry coordinates evaluate_split scores (same subsample, eligibility and mask)."""
    idx = subsample(X.shape[0], labels, n_cells, seed)
    Xs = X[idx]
    eligible = ((~np.isnan(Xs)).mean(0) >= 0.5) & ~np.isnan(Xs).all(0)
    mask = make_mask(Xs, eligible, kind, seed)
    r, c = np.nonzero(mask)
    return idx[r], c


def _split(args):
    import panel_v9 as PV

    X, names, V, warm, Vnull, kind, seed, target = args
    f = CACHE / f"{target}_{kind}_{seed}.npz"
    if f.exists():
        d = np.load(f, allow_pickle=True)
        return kind, seed, d["se"].item(), d["sel"], d["cell"], d["prot"]
    for k, fn in A.ADDENDUM.items():  # registered in-process only; the frozen baselines.py is unchanged
        B.ALL[k] = fn
    t0 = time.time()
    se, sel, wts, secs = PV.evaluate_split(X, names, None, V, warm, kind, seed, target, PANEL, Vnull=Vnull)
    cell, prot = _coords(X, None, kind, seed)
    CACHE.mkdir(parents=True, exist_ok=True)
    tmp = CACHE / f".{target}_{kind}_{seed}.{os.getpid()}.npz"
    np.savez_compressed(tmp, se=np.array(se, dtype=object), sel=sel, cell=cell, prot=prot)
    os.replace(tmp, f)
    print(f"[{time.strftime('%H:%M:%S')}] {target} {kind} {seed} ({time.time() - t0:.0f}s)", flush=True)
    return kind, seed, se, sel, cell, prot


def _rmse_diff_boot(sa, sb, groups, n_boot=1000, seed=17):
    """Cluster bootstrap of RMSE(a) - RMSE(b): resample whole groups (entries if groups is None)."""
    rng = np.random.default_rng(seed)
    point = float(np.sqrt(sa.mean()) - np.sqrt(sb.mean()))
    if groups is None:
        n = len(sa)
        d = [np.sqrt(sa[i].mean()) - np.sqrt(sb[i].mean()) for i in (rng.integers(0, n, n) for _ in range(n_boot))]
    else:
        u, inv = np.unique(groups, return_inverse=True)
        A_ = np.bincount(inv, sa, len(u))
        B_ = np.bincount(inv, sb, len(u))
        N_ = np.bincount(inv, minlength=len(u)).astype(float)
        d = []
        for _ in range(n_boot):
            k = np.bincount(rng.integers(0, len(u), len(u)), minlength=len(u))
            d.append(np.sqrt((k * A_).sum() / (k * N_).sum()) - np.sqrt((k * B_).sum() / (k * N_).sum()))
    lo, hi = np.percentile(d, [2.5, 97.5])
    return point, [float(lo), float(hi)]


def robust_tests(SE, SEL, split_id, prot_id, a, b):
    out = []
    for endpoint, pick in (("all", np.ones_like(SEL, bool)), ("in_complex", SEL)):
        sa, sb = SE[a][pick], SE[b][pick]
        ref = float(np.sqrt(sb.mean()))
        row = {"method": a, "vs": b, "endpoint": endpoint}
        for scheme, g in (("entry", None), ("protein", prot_id[pick]), ("split", split_id[pick])):
            point, ci = _rmse_diff_boot(sa, sb, g)
            row[scheme] = {"pct": 100 * point / ref, "ci95_pct": [100 * x / ref for x in ci],
                           "pass_prereg": bool(ci[1] < -SIG_TOL and 100 * point / ref <= -1.0),
                           "pass_strict": bool(100 * ci[1] / ref <= -1.0)}
        out.append(row)
    return out


def main():
    a = sys.argv
    jobs = int(a[a.index("--jobs") + 1]) if "--jobs" in a else 1
    if "--dev" in a:
        from explore_phaseb import load_target

        target = a[a.index("--dev") + 1]
        X, names, _ = load_target(target)
        seeds, refs, excl, frozen = [100, 101], reference_datasets(exclude=(target,)), (target,), None
        tag = f"dev_{target}"
    else:
        import confirm_v9 as C

        target = a[1]
        fz = ROOT / "results" / "confirm_v9" / f"{target}.json"
        if not fz.exists():
            sys.exit(f"refusing: the preregistered run for {target} has not produced its result yet")
        frozen = json.loads(fz.read_text())
        X, names = C.load_holdout(target)
        seeds, refs, excl, tag = C.SEEDS, reference_datasets(), (), target
    Z = T.corpus_reference(names, excl, refs)
    V, warm = T.fit_loadings(Z)
    Vnull = T.permute_rows(V, warm, 0)
    jl = [(X, names, V, warm, Vnull, k, s, target) for k in ("MCAR", "MNAR") for s in seeds]
    if jobs > 1:
        with ProcessPoolExecutor(jobs) as ex:
            res = list(ex.map(_split, jl))
    else:
        res = [_split(j) for j in jl]
    out = {"dataset": target, "seeds": seeds, "post_hoc": True, "masks": {}}
    for kind in ("MCAR", "MNAR"):
        rr = [r for r in res if r[0] == kind]
        SE = {m: np.concatenate([r[2][m] for r in rr]) for m in rr[0][2]}
        SEL = np.concatenate([r[3] for r in rr])
        split_id = np.concatenate([np.full(len(r[3]), r[1]) for r in rr])
        prot_id = np.concatenate([r[5] for r in rr])
        rmse = {m: float(np.sqrt(v.mean())) for m, v in SE.items()}
        block = {"rmse": rmse, "rel_to_champion": {m: v / rmse["champion"] for m, v in rmse.items()}}
        if frozen is not None:
            fr = frozen["masks"][kind]["rmse"]
            block["reproduction_max_abs_diff"] = max(abs(fr[m] - rmse[m]) for m in rmse if m in fr)
        block["H16"] = [t for c in COMPS16 for t in robust_tests(SE, SEL, split_id, prot_id, "champion", c)]
        block["H17"] = robust_tests(SE, SEL, split_id, prot_id, "revised", "champion")
        block["H18"] = [t for c in COMPS18 for t in robust_tests(SE, SEL, split_id, prot_id, "revised", c)]
        block["addendum_vs"] = {m: robust_tests(SE, SEL, split_id, prot_id, "champion", m)
                                for m in SE if m.split("+")[0] in A.ADDENDUM}
        block["revised_vs_addendum"] = {m: robust_tests(SE, SEL, split_id, prot_id, "revised", m)
                                        for m in SE if m.split("+")[0] in A.ADDENDUM}
        out["masks"][kind] = block
    for h in ("H16", "H17", "H18"):
        out[f"{h}_passes"] = {s: {r: f"{sum(t[s][r] for m in out['masks'].values() for t in m[h])}/"
                                     f"{sum(1 for m in out['masks'].values() for t in m[h])}"
                                  for r in ("pass_prereg", "pass_strict")} for s in ("entry", "protein", "split")}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{tag}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"dataset": target, **{f"{h}_passes": out[f"{h}_passes"] for h in ("H16", "H17", "H18")},
                      "repro": {k: v.get("reproduction_max_abs_diff") for k, v in out["masks"].items()}}, indent=1))


if __name__ == "__main__":
    main()
