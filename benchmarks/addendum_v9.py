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
VERSION = "v3"  # bump when addendum code changes: part of every cache key
CACHE = ROOT / "data" / "cache" / "v9_addendum" / VERSION
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

    from scp_graphfs_mcp.phasec import cshrink as CS

    X, names, labels, V, warm, Vnull, kind, seed, target, mode = args
    f = CACHE / f"{mode}_{target}_{kind}_{seed}.npz"
    if f.exists():
        d = np.load(f, allow_pickle=True)
        se = d["se"].item()
        if set(PANEL) <= set(se):
            return kind, seed, se, d["sel"], d["cell"], d["prot"]
    for k, fn in A.ADDENDUM.items():  # registered in-process only; the frozen baselines.py is unchanged
        B.ALL[k] = fn
    CS.CACHE = CACHE / "cshrink_independent"  # fresh psi-selection scores: the reproduction check is independent
    t0 = time.time()
    se, sel, wts, secs = PV.evaluate_split(X, names, labels, V, warm, kind, seed, target, PANEL, Vnull=Vnull)
    cell, prot = _coords(X, labels, kind, seed)
    assert len(cell) == len(sel) == len(se["champion"]), "coordinate / error length mismatch"
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


def _split_t(sa, sb, split_id):
    """Per-split RMSE differences (G splits); t-interval for their mean, in RMSE units (percentile bootstraps are
    degenerate with ~10 clusters)."""
    from scipy import stats

    g = np.unique(split_id)
    d = np.array([np.sqrt(sa[split_id == k].mean()) - np.sqrt(sb[split_id == k].mean()) for k in g])
    half = stats.t.ppf(0.975, len(g) - 1) * d.std(ddof=1) / np.sqrt(len(g)) if len(g) > 1 else np.inf
    return float(d.mean()), [float(d.mean() - half), float(d.mean() + half)], len(g), int((d < 0).sum())


def robust_tests(SE, SEL, split_id, prot_id, a, b):
    out = []
    for endpoint, pick in (("all", np.ones_like(SEL, bool)), ("in_complex", SEL)):
        sa, sb = SE[a][pick], SE[b][pick]
        ref = float(np.sqrt(sb.mean()))
        row = {"method": a, "vs": b, "endpoint": endpoint}
        for scheme, g in (("entry", None), ("protein", prot_id[pick])):
            point, ci = _rmse_diff_boot(sa, sb, g)
            row[scheme] = {"pct": 100 * point / ref, "ci95_pct": [100 * x / ref for x in ci],
                           "n_clusters": int(len(sa) if g is None else len(np.unique(g))),
                           "pass_prereg": bool(ci[1] < -SIG_TOL and 100 * point / ref <= -1.0),
                           "pass_strict": bool(100 * ci[1] / ref <= -1.0)}
        mean_d, ci, G, n_better = _split_t(sa, sb, split_id[pick])
        row["split_t"] = {"mean_pct": 100 * mean_d / ref, "ci95_pct": [100 * x / ref for x in ci], "n_splits": G,
                          "n_splits_better": n_better, "valid": G >= 5,
                          "pass_prereg": bool(G >= 5 and ci[1] < 0 and 100 * mean_d / ref <= -1.0),
                          "pass_strict": bool(G >= 5 and 100 * ci[1] / ref <= -1.0)}
        out.append(row)
    return out


def main():
    a = sys.argv
    jobs = int(a[a.index("--jobs") + 1]) if "--jobs" in a else 1
    if "--dev" in a:
        from explore_phaseb import load_target

        target = a[a.index("--dev") + 1]
        assert target not in (T.PHASE_B_HOLDOUTS | T.HOLDOUTS | T.V9_HOLDOUTS), "holdouts are never dev targets"
        X, names, labels = load_target(target)  # stratified like the frozen dev panel
        seeds, refs, excl, frozen = [100, 101], reference_datasets(exclude=(target,)), (target,), None
        tag, mode = f"dev_{target}", "dev"
    else:
        import confirm_v9 as C

        target = a[1]
        fz = ROOT / "results" / "confirm_v9" / f"{target}.json"
        if not fz.exists():
            sys.exit(f"refusing: the preregistered run for {target} has not produced its result yet")
        frozen = json.loads(fz.read_text())
        assert target in C.PRIMARY + C.SECONDARY
        X, names = C.load_holdout(target)
        labels, seeds, refs, excl, tag, mode = None, C.SEEDS, reference_datasets(), (), target, "holdout"
    assert not set(refs) & (T.PHASE_B_HOLDOUTS | T.HOLDOUTS | T.V9_HOLDOUTS)
    Z = T.corpus_reference(names, excl, refs)
    V, warm = T.fit_loadings(Z)
    Vnull = T.permute_rows(V, warm, 0)
    jl = [(X, names, labels, V, warm, Vnull, k, s, target, mode) for k in ("MCAR", "MNAR") for s in seeds]
    if jobs > 1:
        with ProcessPoolExecutor(jobs) as ex:
            res = list(ex.map(_split, jl))
    else:
        res = [_split(j) for j in jl]
    out = {"dataset": target, "seeds": seeds, "post_hoc": True, "addendum_version": VERSION, "masks": {}}
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
                                  for r in ("pass_prereg", "pass_strict")} for s in ("entry", "protein", "split_t")}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{tag}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"dataset": target, **{f"{h}_passes": out[f"{h}_passes"] for h in ("H16", "H17", "H18")},
                      "repro": {k: v.get("reproduction_max_abs_diff") for k, v in out["masks"].items()}}, indent=1))


if __name__ == "__main__":
    main()
