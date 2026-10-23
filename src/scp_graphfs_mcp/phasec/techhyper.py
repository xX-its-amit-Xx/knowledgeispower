"""Phase-C idea `techhyper` (queen: biological priors; round 4, built on the cshrink/tgtshrink-null base).

A technology/tissue-kernel hyperprior on the shrinkage strength psi of the target's per-protein mean and log-sd.
Hierarchical empirical Bayes on a hyperparameter: every human reference dataset with >= 40 cells is a pseudo-target
(30-cell subsample, harness mask, 50/50 MCAR/MNAR inner mask), which gives a psi-loss curve per dataset; a metadata-
kernel average of those curves is the prior on the target's own (noisy) inner-mask curve.

  GRID = (1, 2, 4, 8, 16, 32, 64, INF); INF = 1e6 (tau2 at its floor, w ~ 0); grid index g = 0..7.
  Scoring of a curve (identical to tgtshrink._run): Msel = tgtshrink.selection_mask(M, seed); for psi in GRID:
     mu, sd = tgtshrink.shrink(Msel, psi, "null", seed); Y = nfbase.adapt_project_mb(Msel, V, warm, base, mu, sd,
     lam_v=100, extra=10, iters=5); S = 0.5 (RMSE_MCAR-part + RMSE_MNAR-part). Relative curve l = log S - mean log S.
  Step 1 (precompute, target-independent): lbar_d = mean over 4 replicates (MCAR/MNAR x seeds 900/901); s_t = sqrt(mean
     over d, psi of var_r l_{d,r}(psi)).
  Step 2 (prior, LODO over pseudo-targets != target): w_d = exp(-dist/1); lp = sum w lbar / sum w;
     s_p = max(s_t/2, sqrt(mean_psi weighted var of lbar around lp)); pi ∝ exp(-lp/s_p) (FLAT: uniform).
     dist: FULL = tech' + tissue, TECH = tech', TISSUE = tissue; tech' = ttprior.tech_dist + carrier term
     (0 same, 0.5 differ, 0.25 if either unknown). NULL = FULL with ttprior.null_meta_map applied to reference rows.
  Step 3 (impute): target curve lt (cached on disk by matrix hash), P ∝ pi exp(-lt/s_t), gbar = sum P g,
     psi* = INF if gbar >= 6.5 else 2**gbar; return adapt_project_mb(Xm, V, warm, base, *shrink(Xm, psi*, "null", seed)).
Entry points: impute (FULL), impute_null (NULL), impute_tech, impute_tissue, impute_base (FLAT, posterior mean),
impute_argmin (argmin lt over GRID).
    python -m scp_graphfs_mcp.phasec.techhyper --precompute     # curves.json
    python -m scp_graphfs_mcp.phasec.techhyper --diag           # results/phaseC/techhyper_diag.json
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from functools import lru_cache
from pathlib import Path

import numpy as np

from scp_graphfs_mcp import transfer as T
from scp_graphfs_mcp.phasec.nfbase import EXTRA, LAM_V, adapt_project_mb
from scp_graphfs_mcp.phasec.tgtshrink import selection_mask, shrink
from scp_graphfs_mcp.phasec.ttprior import META, U_, null_meta_map, tech_dist, tissue_dist

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "techhyper"
CURVES = CACHE / "curves.json"
TCURVES = CACHE / "tcurves"
OUT = ROOT / "results" / "phaseC"
FORBIDDEN = {"derks2022", "wu2026", "brunner2022", "khan2023"}
INF = 1e6
GRID = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, INF)
GIDX = np.arange(len(GRID), dtype=float)
TAU = 1.0
PSEUDO_SEEDS = (900, 901)
MIN_CELLS_PSEUDO = 40

# Step 0: carrier channel (isobaric carrier/booster), only from the ttprior META source strings, the scpdata manual
# (data/cache/phasec/ttprior/meta/scpdata_data.R) and the cached PRIDE project JSONs those strings cite.
# 'no' only where label-free quantification is stated; 'yes' only where a carrier/booster channel is stated.
CARRIER = {
    "petrosius2023": U_,          # META: directDIA Spectronaut 17; no labelling statement
    "leduc2022": "yes",           # META 'TMT-18 + carrier'; scpdata leduc2022_pSCoPE 'carrier channels'
    "specht2021": "yes",          # META 'SCoPE2 TMT-11/16 carrier'; scpdata specht2019v3 'carrier channels'
    "schoof2021": "yes",          # META / scpdata schoof2021 '200 cells (booster channel)'
    "ai2025a": U_,                # META: DIA-NN report; no labelling statement
    "amrute2026_MFS": "no",       # META 'label-free'; PXD078860 'label-free single-cell mass spectrometry'
    "bubis2025_A549": "no",       # META 'label-free'; PXD064564 quantificationMethods 'label-free quantification'
    "cong2020AC": "no",           # META 'label-free'; scpdata cong2020AC 'label-free quantitative information'
    "gregoire2023_mixCTRL": "yes",  # scpdata gregoire2023_mixCTRL 'single-cell channels, carrier channels'
    "guise2024": U_,              # META: nanoPOTS, Exploris 480; no labelling statement
    "hu2023_K562": "no",          # META 'label-free'
    "inns2026_CCS": "no",         # PXD073250 description 'a label free SCP methodology'
    "krull2024": U_,              # META: diaPASEF timsTOF Pro; no labelling statement
    "liang2020_hela": U_,         # META / scpdata liang2020_hela: no labelling statement
    "px063590_ADbrain": "yes",    # META 'TMTpro + carrier'; PXD063590 '126 (carrier)'
    "sadiku2025": U_,             # PXD061065: 'Relative quantification'; no labelling statement
    "weiss2026_liver": "no",      # PXD062231 quantificationMethods 'label-free quantification'
    "williams2020_lfq": "no",     # META / scpdata 'label-free nanoPOTS'
    "williams2020_tmt": "yes",    # META / scpdata 'carrier (10 ng)'
    "woo2022_lung": "no",         # META 'label-free'
    "woo2025_PC3": "no",          # META 'LFQ DIA'; PXD059079 'label-free quantification'
    "ye2025_EB": "no",            # PXD049181 description '80 label-free SCP samples per day'
}
assert set(CARRIER) == set(META)


def meta(d):
    return {**META[d], "carrier": CARRIER[d]}


def carrier_term(a, b):
    if U_ in (a, b):
        return 0.25
    return 0.5 * (a != b)


def tech2(m, t):
    return tech_dist(m, t) + carrier_term(m["carrier"], t["carrier"])


def kdist(m, t, kernel):
    if kernel == "full":
        return tech2(m, t) + tissue_dist(m, t)
    if kernel == "tech":
        return tech2(m, t)
    if kernel == "tissue":
        return tissue_dist(m, t)
    raise ValueError(kernel)


# ---------------------------------------------------------------------------------------------------------------
def _bench():
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H  # noqa: E402

    return H


def pseudo_targets():
    H = _bench()
    out = []
    for d in H.reference_datasets():
        assert d not in FORBIDDEN | T.PHASE_B_HOLDOUTS | T.HOLDOUTS
        if d in T.NON_HUMAN:
            continue
        X, _, _ = H.load_target(d)
        if X.shape[0] >= MIN_CELLS_PSEUDO:
            out.append(d)
    return out


def score_curve(M, V, warm, base, seed):
    """S(psi) for psi in GRID, exactly tgtshrink._run's scoring with mode='null'."""
    Msel, (r1, c1), (r2, c2) = selection_mask(M, seed)
    S = []
    for p in GRID:
        mu, sd, _ = shrink(Msel, p, "null", seed)
        Y = adapt_project_mb(Msel, V, warm, base, mu, sd, lam_v=LAM_V, extra=EXTRA, iters=5)
        e1 = np.sqrt(np.mean((Y[r1, c1] - M[r1, c1]) ** 2))
        e2 = np.sqrt(np.mean((Y[r2, c2] - M[r2, c2]) ** 2))
        S.append(0.5 * (e1 + e2))
    return np.array(S)


def rel(S):
    L = np.log(np.asarray(S, float))
    return L - L.mean()


def pseudo_split(d, kind, seed):
    """Masked 30-cell input for pseudo-target d, exactly as benchmarks/phasec.run_split builds the eval input."""
    H = _bench()
    from scp_graphfs_mcp import phase2 as P

    X, names, labels = H.load_target(d)
    idx = H.subsample(X.shape[0], labels, 30, seed)
    Xs = X[idx]
    eligible = ((~np.isnan(Xs)).mean(0) >= 0.5) & ~np.isnan(Xs).all(0)
    mask = H.make_mask(Xs, eligible, kind, seed)
    Xm = Xs.copy()
    Xm[mask] = np.nan
    b = P._fill(Xm, P.loo_knn(Xm))
    b = np.where(np.isnan(b), np.nanmean(Xm), b)
    return Xm, b, names


def _pseudo_job(args):
    d, kind, seed, vfile = args
    t0 = time.time()
    Xm, b, names = pseudo_split(d, kind, seed)
    if vfile is None:
        V, warm = _bench().target_models(d, names)
    else:
        z = np.load(vfile)
        V, warm = z["V"], z["warm"]
    S = score_curve(Xm, V, warm, b, seed)
    return {"d": d, "kind": kind, "seed": seed, "S": S.tolist(), "l": rel(S).tolist(), "sec": round(time.time() - t0, 1)}


def _atomic_write(f, text):
    f.parent.mkdir(parents=True, exist_ok=True)
    tmp = f.with_suffix(f".tmp{os.getpid()}")
    tmp.write_text(text)
    os.replace(tmp, f)


def precompute(jobs=2):
    from concurrent.futures import ProcessPoolExecutor

    H = _bench()
    ds = pseudo_targets()
    print("pseudo-targets:", ds, flush=True)
    for d in ds:  # V_{-d} fits serially
        t0 = time.time()
        _, names, _ = H.load_target(d)
        H.target_models(d, names)
        print(f"V {d} {time.time() - t0:.1f}s", flush=True)
    todo = [(d, k, s, None) for d in ds for k in ("MCAR", "MNAR") for s in PSEUDO_SEEDS]
    with ProcessPoolExecutor(jobs) as ex:
        res = []
        for r in ex.map(_pseudo_job, todo):
            print(r["d"], r["kind"], r["seed"], r["sec"], np.round(r["l"], 3).tolist(), flush=True)
            res.append(r)
    out = summarise_curves(res, ds)
    _atomic_write(CURVES, json.dumps(out, indent=1) + "\n")
    return out


def summarise_curves(res, ds):
    reps = {d: [r for r in res if r["d"] == d] for d in ds}
    lbar = {d: np.mean([r["l"] for r in reps[d]], 0).tolist() for d in ds}
    var = [np.var([r["l"] for r in reps[d]], 0, ddof=1) for d in ds]
    s_t = float(np.sqrt(np.mean(var)))
    return {"grid": [float(g) for g in GRID], "datasets": ds, "seeds": list(PSEUDO_SEEDS), "replicates": res,
            "lbar": lbar, "s_t": s_t}


@lru_cache(maxsize=1)
def load_curves():
    c = json.loads(CURVES.read_text())
    return {"datasets": c["datasets"], "lbar": {d: np.array(v) for d, v in c["lbar"].items()}, "s_t": c["s_t"]}


# ---------------------------------------------------------------------------------------------------------------
def weights(target, D, kernel, null=False):
    if kernel == "flat":
        return {d: 1.0 for d in D}
    mp = null_meta_map(target, D) if null else {d: d for d in D}
    t = meta(target)
    return {d: float(np.exp(-kdist(meta(mp[d]), t, kernel) / TAU)) for d in D}


def prior_curve(target, kernel, null=False, lbar=None, D=None):
    c = load_curves()
    lbar = lbar if lbar is not None else c["lbar"]
    D = D if D is not None else [d for d in c["datasets"] if d != target]
    w = weights(target, D, kernel, null)
    W = np.array([w[d] for d in D])
    L = np.vstack([lbar[d] for d in D])
    lp = (W[:, None] * L).sum(0) / W.sum()
    s_p = float(np.sqrt(np.mean((W[:, None] * (L - lp) ** 2).sum(0) / W.sum())))
    s_p = max(s_p, c["s_t"] / 2)
    return lp, s_p, w


def prior_pi(target, kernel, null=False):
    if kernel == "flat":
        return np.full(len(GRID), 1.0 / len(GRID)), None
    lp, s_p, _ = prior_curve(target, kernel, null)
    a = -lp / s_p
    pi = np.exp(a - a.max())
    return pi / pi.sum(), s_p


def posterior(lt, pi, s_t):
    a = np.log(pi) - np.asarray(lt) / s_t
    P = np.exp(a - a.max())
    P /= P.sum()
    gbar = float(P @ GIDX)
    psi = INF if gbar >= 6.5 else float(2.0 ** gbar)
    return P, gbar, psi


def _mhash(M, seed, target):
    h = hashlib.sha1(np.ascontiguousarray(M).tobytes())
    h.update(f"{seed}|{target}|{M.shape}".encode())
    return h.hexdigest()[:20]


def target_curve(M, ctx):
    """Target's inner-mask curve S(psi), cached on disk (shared by every entry point and both harness calls)."""
    f = TCURVES / f"{ctx['target']}_{ctx['seed']}_{_mhash(M, ctx['seed'], ctx['target'])}.json"
    if f.exists():
        return np.array(json.loads(f.read_text())["S"])
    S = score_curve(M, ctx["V"], ctx["warm"], ctx["base"], ctx["seed"])
    _atomic_write(f, json.dumps({"S": S.tolist()}))
    return S


def _log(variant, rec):
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / f"log_{variant}.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def _run(Xm, ctx, kernel, null=False, variant="full", rule="mean"):
    target, seed = ctx["target"], ctx["seed"]
    assert target not in FORBIDDEN
    S = target_curve(Xm, ctx)
    lt = rel(S)
    s_t = load_curves()["s_t"]
    if rule == "argmin":
        g = int(np.argmin(lt))  # ties -> smaller psi
        pi, P, gbar, psi, s_p = None, None, float(g), GRID[g], None
    else:
        pi, s_p = prior_pi(target, kernel, null)
        P, gbar, psi = posterior(lt, pi, s_t)
    mu, sd, _ = shrink(Xm, psi, "null", seed)
    _log(variant, {"target": target, "seed": seed, "n_missing": int(np.isnan(Xm).sum()), "S": S.tolist(),
                   "lt": lt.tolist(), "pi": None if pi is None else pi.tolist(),
                   "P": None if P is None else P.tolist(), "s_p": s_p, "s_t": s_t, "gbar": gbar, "psi": psi})
    return adapt_project_mb(Xm, ctx["V"], ctx["warm"], ctx["base"], mu, sd, lam_v=LAM_V, extra=EXTRA)


def impute(Xm, ctx):
    return _run(Xm, ctx, "full", False, "full")


def impute_null(Xm, ctx):
    return _run(Xm, ctx, "full", True, "null")


def impute_tech(Xm, ctx):
    return _run(Xm, ctx, "tech", False, "tech")


def impute_tissue(Xm, ctx):
    return _run(Xm, ctx, "tissue", False, "tissue")


def impute_base(Xm, ctx):
    return _run(Xm, ctx, "flat", False, "base")


def impute_argmin(Xm, ctx):
    return _run(Xm, ctx, "flat", False, "argmin", rule="argmin")


# ---------------------------------------------------------------------------------------------------------------
# Step 4 diagnostics
def _refit_V_excl(t, d):
    """V_{-t,-d}: fit_loadings on every reference dataset except t and d, aligned to d's proteins (cached)."""
    H = _bench()
    f = CACHE / "leak" / f"V_{d}_minus_{t}.npz"
    if f.exists():
        return f
    _, names, _ = H.load_target(d)
    parts = []
    for dd in H.reference_datasets(exclude=(d, t)):
        parts.append(T.corpus_reference(names, (d, t), [dd]))
    V, warm = T.fit_loadings(np.vstack(parts))
    f.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(f, V=V, warm=warm)
    return f


def _target_eval_curves(t):
    """The target's harness eval-input curves for the fast pool (2 masks x seeds 100/101)."""
    H = _bench()
    out = []
    for kind in ("MCAR", "MNAR"):
        for seed in H.FAST_SEEDS:
            Xm, b, names = pseudo_split(t, kind, seed)
            V, warm = H.target_models(t, names)
            ctx = {"target": t, "seed": seed, "V": V, "warm": warm, "base": b}
            out.append((kind, seed, rel(target_curve(Xm, ctx))))
    return out


def diagnostics(jobs=2):
    from concurrent.futures import ProcessPoolExecutor

    c = load_curves()
    ds, lbar, s_t = c["datasets"], c["lbar"], c["s_t"]
    out = {"s_t": s_t, "grid": [float(g) for g in GRID], "datasets": ds}
    # (a) LODO curve prediction
    a = {}
    for d in ds:
        row = {}
        for name, kern, nl in (("FULL", "full", False), ("TECH", "tech", False), ("TISSUE", "tissue", False),
                               ("FLAT", "flat", False), ("NULL", "full", True)):
            lp, s_p, _ = prior_curve(d, kern, nl)
            row[name] = float(np.sqrt(np.mean((lbar[d] - lp) ** 2)))
        a[d] = row
    out["a_lodo_rmse"] = a
    out["a_mean"] = {k: float(np.mean([a[d][k] for d in ds])) for k in ("FULL", "TECH", "TISSUE", "FLAT", "NULL")}
    out["a_n_full_better_than_null"] = int(sum(a[d]["FULL"] < a[d]["NULL"] for d in ds))
    out["a_n_full_better_than_flat"] = int(sum(a[d]["FULL"] < a[d]["FLAT"] for d in ds))
    # (b) argmin psi of each dataset's mean curve next to its metadata
    out["b_argmin"] = {d: {"argmin_psi": float(GRID[int(np.argmin(lbar[d]))]), "lbar": np.round(lbar[d], 4).tolist(),
                           "acq": META[d]["acq"], "instrument": META[d]["instrument"], "carrier": CARRIER[d],
                           "sample_type": META[d]["sample_type"], "lineage": sorted(META[d]["lineage"])}
                       for d in ds}
    # prior summaries for the fast targets
    out["prior_fast"] = {}
    for t in ("petrosius2023", "leduc2022", "weiss2026_liver", "gregoire2023_mixCTRL", "px063590_ADbrain",
              "sadiku2025"):
        r = {}
        for name, kern, nl in (("FULL", "full", False), ("TECH", "tech", False), ("TISSUE", "tissue", False),
                               ("NULL", "full", True)):
            pi, s_p = prior_pi(t, kern, nl)
            r[name] = {"pi": np.round(pi, 4).tolist(), "s_p": s_p, "prior_gbar": float(pi @ GIDX)}
        _, _, w = prior_curve(t, "full")
        r["top3_full"] = sorted(w, key=w.get, reverse=True)[:3]
        r["w_full"] = w
        out["prior_fast"][t] = r
    # (c) leakage: V_{-d} contains the target's rows; refit V_{-t,-d} for the 3 top-weight d
    leak = {}
    for t in ("petrosius2023", "sadiku2025"):
        top = out["prior_fast"][t]["top3_full"]
        files = {d: str(_refit_V_excl(t, d)) for d in top}
        todo = [(d, k, s, files[d]) for d in top for k in ("MCAR", "MNAR") for s in PSEUDO_SEEDS]
        with ProcessPoolExecutor(jobs) as ex:
            res = list(ex.map(_pseudo_job, todo))
        lbar2 = dict(lbar)
        per_d = {}
        for d in top:
            l2 = np.mean([r["l"] for r in res if r["d"] == d], 0)
            lbar2[d] = l2
            per_d[d] = {"corr": float(np.corrcoef(lbar[d], l2)[0, 1]), "rmse": float(np.sqrt(np.mean((lbar[d] - l2) ** 2))),
                        "lbar_leakfree": np.round(l2, 4).tolist()}
        D = [d for d in ds if d != t]
        lp0, sp0, _ = prior_curve(t, "full", lbar=lbar, D=D)
        lp1, sp1, _ = prior_curve(t, "full", lbar=lbar2, D=D)
        pi0 = np.exp(-(lp0 - lp0.min()) / sp0); pi0 /= pi0.sum()
        pi1 = np.exp(-(lp1 - lp1.min()) / sp1); pi1 /= pi1.sum()
        calls = []
        for kind, seed, lt in _target_eval_curves(t):
            _, g0, p0 = posterior(lt, pi0, s_t)
            _, g1, p1 = posterior(lt, pi1, s_t)
            calls.append({"mask": kind, "seed": seed, "gbar": g0, "psi": p0, "gbar_leakfree": g1, "psi_leakfree": p1})
        leak[t] = {"top3": top, "curves": per_d, "prior_gbar": float(pi0 @ GIDX), "prior_gbar_leakfree": float(pi1 @ GIDX),
                   "eval_calls": calls}
    out["c_leakage"] = leak
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "techhyper_diag.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({"a_mean": out["a_mean"], "full<null": out["a_n_full_better_than_null"],
                      "full<flat": out["a_n_full_better_than_flat"], "n": len(ds)}, indent=1))
    return out


if __name__ == "__main__":
    j = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 2
    if "--precompute" in sys.argv:
        precompute(j)
    if "--diag" in sys.argv:
        diagnostics(j)
