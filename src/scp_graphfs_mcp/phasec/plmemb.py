"""Phase-C idea `plmemb` (queen: biological priors, round 5). ProtT5 embedding-kernel smoothing of weakly supported
reference loadings, inside the cshrink pipeline.

  0. E_t: UniProt 2026_03 ProtT5 per-protein embeddings (data/cache/reps/prott5/per-protein.h5, SHA-256 checked in
     data/cache/phasec/plmemb/sha256_check.json), PCA-128 (pca128.npz). Each target protein group -> mean of the
     L2-normalised Z rows of its mapped canonical accessions, L2-normalised again; no mapped accession -> has_emb=False.
  1. Support: n_obs_j, n_ds_j (# ref datasets with >= 5 obs in column j), nrm_j = ||Vref_j||;
     g_j = 1 - 0.5 (rankpct(n_ds) + rankpct(nrm)) over warm proteins (average ties); cold -> g = 1.
     Pool P_t = warm & has_emb & g <= 0.5.
  2. Kernel: k = 15 nearest pool neighbours (cosine, excluding j), w = exp((c - 1)/tau), tau = 0.1, normalised;
     Vs_j = sum w Vref_i; beta_j(b) = min(1, b g_j^2) (warm), 1 for cold & has_emb if b > 0;
     Vtil = (1 - beta) Vref + beta Vs; usable_b = warm | (cold & has_emb & b > 0).
  3. psi* from cshrink's selection (same inner mask, same on-disk score cache); b in {0, .25, .5, 1} selected on the same
     inner mask (b = 0 score = cshrink's own score at psi*); refit on M with iters = 10.

Entry points: impute (true E), impute_null (E rows permuted across has_emb proteins, rng seed + 7177, independent of
n_missing), impute_gauss (iid N(0,1) 128-d, normalised, rng seed + 4242), impute_base (b = 0 == cshrink),
impute_nullfix (permuted E, b pinned to the b* that true E selects for the same call).
`python -m scp_graphfs_mcp.phasec.plmemb diag|unit|report`.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from scp_graphfs_mcp.ids import group_accessions
from scp_graphfs_mcp.phasec import cshrink as CS
from scp_graphfs_mcp.phasec import tgtshrink as TS
from scp_graphfs_mcp.phasec.cshrink import GRID, cshrink_stats
from scp_graphfs_mcp.phasec.nfbase import EXTRA, LAM_V, adapt_project_mb

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "plmemb"
PCA = ROOT / "data" / "cache" / "reps" / "prott5" / "pca128.npz"
H5 = ROOT / "data" / "cache" / "reps" / "prott5" / "per-protein.h5"
H5_SHA = "9d713afe200557b2532aa8743886027772013c7140c6567e7084ee724af0877a"
K, TAU = 15, 0.1
B_GRID = (0.0, 0.25, 0.5, 1.0)
B_GRID_V2 = (0.0, 0.25, 0.5, 1.0, 2.0, 4.0)  # v2: b* hit the grid edge (1.0) in 25/48 fast calls
MIN_DS_OBS = 5


def _atomic_savez(f, **kw):
    f.parent.mkdir(parents=True, exist_ok=True)
    tmp = f.with_name(f.name + f".tmp{os.getpid()}.npz")
    np.savez_compressed(tmp, **kw)
    os.replace(tmp, f)


# ---------------------------------------------------------------------------------------------------------------
# 0. embeddings
@lru_cache(maxsize=1)
def _pca():
    d = np.load(PCA)
    Z = d["Z"].astype(np.float64)
    Z /= np.maximum(np.linalg.norm(Z, axis=1, keepdims=True), 1e-12)
    return {a: i for i, a in enumerate(d["acc"].tolist())}, Z


def target_embeddings(target, names):
    f = CACHE / f"{target}_E.npz"
    if f.exists():
        d = np.load(f)
        if len(d["has_emb"]) == len(names):
            return d["E"], d["has_emb"]
    pos, Z = _pca()
    E = np.zeros((len(names), Z.shape[1]))
    has = np.zeros(len(names), bool)
    for j, nm in enumerate(names):
        rows = [pos[a] for a in dict.fromkeys(a.split("-")[0] for a in group_accessions(str(nm))) if a in pos]
        if rows:
            e = Z[rows].mean(0)
            n = np.linalg.norm(e)
            if n > 1e-12:
                E[j], has[j] = e / n, True
    _atomic_savez(f, E=E, has_emb=has)
    return E, has


# ---------------------------------------------------------------------------------------------------------------
# 1. support diagnostic
def support(target, R, src, V, warm):
    f = CACHE / f"{target}_support.npz"
    if f.exists():
        d = np.load(f)
        if len(d["g"]) == len(warm):
            return {k: d[k] for k in d.files}
    obs = ~np.isnan(R)
    n_obs = obs.sum(0)
    n_ds = np.zeros(R.shape[1], int)
    for s in np.unique(src):
        n_ds += (obs[src == s].sum(0) >= MIN_DS_OBS).astype(int)
    nrm = np.linalg.norm(V, axis=1)
    g = np.ones(R.shape[1])
    w = np.where(warm)[0]
    if len(w):
        rp = lambda x: rankdata(x, method="average") / len(x)  # noqa: E731
        g[w] = 1.0 - 0.5 * (rp(n_ds[w]) + rp(nrm[w]))
    out = {"n_obs": n_obs, "n_ds": n_ds, "nrm": nrm, "g": g}
    _atomic_savez(f, **out)
    return out


def get_support(ctx):
    f = CACHE / f"{ctx['target']}_support.npz"
    if f.exists():
        d = np.load(f)
        if len(d["g"]) == len(ctx["warm"]):
            return {k: d[k] for k in d.files}
    return support(ctx["target"], ctx["ref"](), ctx["ref_src"](), ctx["V"], ctx["warm"])


# ---------------------------------------------------------------------------------------------------------------
# 2. kernel smoothing
def neighbours(E, pool, rows, k=K, exclude_self=True):
    """For each j in rows: indices (into the protein axis) and normalised weights of its k nearest pool proteins."""
    P = np.where(pool)[0]
    C = E[rows] @ E[P].T
    if exclude_self:
        same = rows[:, None] == P[None, :]
        C = np.where(same, -np.inf, C)
    kk = min(k, len(P) - (1 if exclude_self else 0))
    idx = np.argpartition(-C, kk - 1, axis=1)[:, :kk]
    c = np.take_along_axis(C, idx, 1)
    w = np.exp((c - 1.0) / TAU)
    w /= w.sum(1, keepdims=True)
    return P[idx], w


def smoothed(V, E, has, warm, g):
    pool = warm & has & (g <= 0.5)
    Vs = V.copy()
    rows = np.where(has)[0]
    if pool.sum() < K + 2 or len(rows) == 0:
        return Vs, pool
    nb, w = neighbours(E, pool, rows)
    Vs[rows] = np.einsum("rk,rkd->rd", w, V[nb])
    return Vs, pool


def vtilde(V, Vs, has, warm, g, b):
    if b <= 0:
        return V, warm, 0
    beta = np.where(warm & has, np.minimum(1.0, b * g ** 2), 0.0)
    beta = np.where(~warm & has, 1.0, beta)
    Vt = (1 - beta)[:, None] * V + beta[:, None] * Vs
    usable = warm | (~warm & has)
    return Vt, usable, int((beta > 0).sum())


def embedding_variant(E, has, mode, seed):
    if mode == "true":
        return E
    h = np.where(has)[0]
    E2 = E.copy()
    if mode == "perm":
        rng = np.random.default_rng(seed + 7177)
        E2[h] = E[h[rng.permutation(len(h))]]
    elif mode == "gauss":
        rng = np.random.default_rng(seed + 4242)
        G = rng.normal(size=(len(h), E.shape[1]))
        E2[h] = G / np.linalg.norm(G, axis=1, keepdims=True)
    else:
        raise ValueError(mode)
    return E2


# ---------------------------------------------------------------------------------------------------------------
# 3. pipeline
def _bscore(M, Msel, parts, ctx, psi, b, Vt, usable, mu, sd, tag):
    h = hashlib.sha1(np.ascontiguousarray(Msel).tobytes()).hexdigest()[:20]
    f = CACHE / "scores" / f"{ctx['target']}_{ctx['seed']}_{h}_{tag}_{CS._ps(psi)}_{b}.json"
    if f.exists():
        try:
            return float(json.loads(f.read_text())["score"])
        except Exception:
            pass
    (r1, c1), (r2, c2) = parts
    Y = adapt_project_mb(Msel, Vt, usable, ctx["base"], mu, sd, lam_v=LAM_V, extra=EXTRA, iters=5)
    s = float(0.5 * (np.sqrt(np.mean((Y[r1, c1] - M[r1, c1]) ** 2)) + np.sqrt(np.mean((Y[r2, c2] - M[r2, c2]) ** 2))))
    f.parent.mkdir(parents=True, exist_ok=True)
    tmp = f.with_suffix(f".tmp{os.getpid()}")
    tmp.write_text(json.dumps({"score": s}))
    os.replace(tmp, f)
    return s


def _select(M, ctx, mode, force_b=None, b_grid=B_GRID):
    seed = ctx["seed"]
    V, warm = ctx["V"], ctx["warm"]
    Msel, p1, p2 = TS.selection_mask(M, seed)
    pscores = {p: CS._score(M, Msel, (p1, p2), ctx, p, p, cshrink_stats, "c") for p in GRID}
    psi = min(GRID, key=lambda p: (pscores[p], p))
    E, has = target_embeddings(ctx["target"], ctx["names"])
    sup = get_support(ctx)
    g = np.where(warm, sup["g"], 1.0)
    Ev = embedding_variant(E, has, mode, seed)
    Vs, pool = smoothed(V, Ev, has, warm, g)
    bscores = {0.0: pscores[psi]}
    if force_b is None:
        mu_s, sd_s, _ = cshrink_stats(Msel, psi, psi)
        for b in b_grid[1:]:
            Vt, us, _ = vtilde(V, Vs, has, warm, g, b)
            bscores[b] = _bscore(M, Msel, (p1, p2), ctx, psi, b, Vt, us, mu_s, sd_s, mode)
        bstar = min(b_grid, key=lambda b: (bscores[b], b))
    else:
        bstar = force_b
    return {"psi": psi, "pscores": pscores, "bscores": bscores, "b": bstar, "Vs": Vs, "has": has, "g": g,
            "pool": pool}


def run(M, ctx, variant, mode="true", force_b=None, pin_from_true=False, b_grid=B_GRID):
    if pin_from_true:
        force_b = _select(M, ctx, "true", b_grid=b_grid)["b"]
    s = _select(M, ctx, mode, force_b, b_grid)
    V, warm = ctx["V"], ctx["warm"]
    Vt, usable, moved = vtilde(V, s["Vs"], s["has"], warm, s["g"], s["b"])
    mu, sd, info = cshrink_stats(M, s["psi"], s["psi"])
    rec = {"target": ctx.get("target"), "seed": ctx["seed"], "n_missing": int(np.isnan(M).sum()), "mode": mode,
           "chosen_psi": CS._ps(s["psi"]), "chosen_b": s["b"], "b_scores": {str(k): v for k, v in s["bscores"].items()},
           "rows_moved": moved, "n_cold_made_usable": int((usable & ~warm).sum()), "pool": int(s["pool"].sum()),
           "has_emb": int(s["has"].sum())}
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / f"log_{variant}.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    return adapt_project_mb(M, Vt, usable, ctx["base"], mu, sd, lam_v=LAM_V, extra=EXTRA, iters=10)


def impute(Xm, ctx):
    return run(Xm, ctx, "main", "true")


def impute_null(Xm, ctx):
    return run(Xm, ctx, "null", "perm")


def impute_gauss(Xm, ctx):
    return run(Xm, ctx, "gauss", "gauss")


def impute_base(Xm, ctx):
    return run(Xm, ctx, "base", "true", force_b=0.0)


def impute_nullfix(Xm, ctx):
    return run(Xm, ctx, "nullfix", "perm", pin_from_true=True)


# v2 (one revision): b grid widened to {0, .25, .5, 1, 2, 4}; beta = min(1, b g^2) unchanged
def impute_v2(Xm, ctx):
    return run(Xm, ctx, "v2_main", "true", b_grid=B_GRID_V2)


def impute_v2_null(Xm, ctx):
    return run(Xm, ctx, "v2_null", "perm", b_grid=B_GRID_V2)


def impute_v2_gauss(Xm, ctx):
    return run(Xm, ctx, "v2_gauss", "gauss", b_grid=B_GRID_V2)


# ---------------------------------------------------------------------------------------------------------------
# offline: sha check, diag, unit check
def sha_check():
    h = hashlib.sha256()
    with open(H5, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    d = h.hexdigest()
    out = {"file": str(H5.relative_to(ROOT)), "bytes": H5.stat().st_size, "sha256": d, "expected": H5_SHA,
           "match": d == H5_SHA}
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "sha256_check.json").write_text(json.dumps(out, indent=1) + "\n")
    link = CACHE / "per-protein.h5"
    if not link.exists():
        link.symlink_to(H5)
    return out


def _harness():
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H
    return H


def diag(targets=None):
    H = _harness()
    from explore_phaseb import load_target
    targets = targets or H.FAST
    out = {}
    for t in targets:
        X, names, _ = load_target(t)
        V, warm = H.target_models(t, names)
        R, src = H._ref_parts(t, names)
        sup = support(t, R, src, V, warm)
        del R
        E, has = target_embeddings(t, names)
        g = np.where(warm, sup["g"], 1.0)
        pool = warm & has & (g <= 0.5)
        rows = np.where(pool)[0]
        res = {}
        for lab, Ev in (("true", E), ("perm", embedding_variant(E, has, "perm", 0)),
                        ("gauss", embedding_variant(E, has, "gauss", 0))):
            nb, w = neighbours(Ev, pool, rows)
            Vh = np.einsum("rk,rkd->rd", w, V[nb])
            Va = V[rows]
            cos = np.sum(Vh * Va, 1) / np.maximum(np.linalg.norm(Vh, axis=1) * np.linalg.norm(Va, axis=1), 1e-12)
            r2 = 1 - np.sum((Va - Vh) ** 2) / np.sum((Va - Va.mean(0)) ** 2)
            res[lab] = {"mean_cos": round(float(np.mean(cos)), 4), "R2": round(float(r2), 4)}
        out[t] = {"p": int(len(names)), "n_cold": int((~warm).sum()), "n_g_gt_0.5": int((g > 0.5).sum()),
                  "n_has_emb": int(has.sum()), "n_cold_has_emb": int((~warm & has).sum()), "n_pool": int(pool.sum()),
                  "median_n_ds_warm": float(np.median(sup["n_ds"][warm])),
                  "loo_loading_prediction_pool": res}
        print(t, json.dumps(out[t]), flush=True)
    note = ("perm uses rng default_rng(0 + 7177), gauss default_rng(0 + 4242) (seed 0 = the diag); loo = leave-one-out "
            "kernel prediction (k=15, tau=0.1) of each pool protein's Vref row from its pool neighbours")
    f = ROOT / "results" / "phaseC" / "plmemb_diag.json"
    f.write_text(json.dumps({"note": note, "targets": out}, indent=1) + "\n")
    return out


def unit_check():
    H = _harness()
    from scp_graphfs_mcp import phase2 as P
    from scp_graphfs_mcp.phasec import nfbase as NF
    t, seed = "petrosius2023", 100
    Xs, Xm, mask, names = NF._split(t, "MCAR", seed)
    V, warm = H.target_models(t, names)
    b = P._fill(Xm, P.loo_knn(Xm))
    b = np.where(np.isnan(b), np.nanmean(Xm), b)
    cache = {}

    def ref():
        if "Z" not in cache:
            cache["Z"], cache["src"] = H._ref_parts(t, names)
        return cache["Z"]

    def ref_src():
        ref()
        return cache["src"]

    ctx = {"V": V, "warm": warm, "base": b, "seed": seed, "target": t, "names": names, "ref": ref, "ref_src": ref_src}
    mine = run(Xm.copy(), ctx, "unit", "true", force_b=0.0)
    refy = CS.run(Xm.copy(), ctx, "unit", log=False)
    miss = np.isnan(Xm)
    d = float(np.max(np.abs(mine[miss] - refy[miss])))
    # vtilde(b=0) identity and b>0 sanity
    E, has = target_embeddings(t, names)
    sup = get_support(ctx)
    g = np.where(warm, sup["g"], 1.0)
    Vs, pool = smoothed(V, E, has, warm, g)
    Vt0, us0, m0 = vtilde(V, Vs, has, warm, g, 0.0)
    out = {"impute_base_vs_cshrink_impute_petrosius2023_MCAR_100": {"max_abs_diff": d, "pass": bool(d < 1e-9)},
           "vtilde_b0_identity": bool(np.array_equal(Vt0, V) and np.array_equal(us0, warm) and m0 == 0)}
    (CACHE / "unit_check.json").write_text(json.dumps(out, indent=1) + "\n")
    print(out)
    return out


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "sha":
        print(sha_check())
    elif cmd == "diag":
        diag()
    elif cmd == "unit":
        unit_check()


def report(ids=("plmemb", "plmemb_base", "plmemb_gauss"), pool="fast"):
    """Per-mask means from the result JSONs and chosen-b histograms from the run logs."""
    from collections import Counter
    out = {"per_mask": {}, "per_target": {}, "b_hist": {}}
    for i in ids:
        f = ROOT / "results" / "phaseC" / f"{i}_{pool}.json"
        if not f.exists():
            continue
        rows = json.loads(f.read_text())["summary"]["rmse_all"]["per_target_mask"]
        pm = {}
        for m in ("MCAR", "MNAR"):
            rr = [r for r in rows if r["mask"] == m]
            pm[m] = {k: round(float(np.mean([r[k] for r in rr])), 4) for k in rr[0]
                     if k not in ("target", "mask") and rr[0][k] is not None}
        out["per_mask"][i] = pm
        out["per_target"][i] = {f"{r['target']}|{r['mask']}": {k: r[k] for k in ("cand_raw", "cand_stack", "null_stack")
                                                               if k in r} for r in rows}
    for v in ("main", "null", "gauss", "nullfix"):
        f = CACHE / f"{pool}_log_{v}.jsonl"
        if f.exists():
            recs = [json.loads(x) for x in f.read_text().splitlines() if x.strip()]
            out["b_hist"][v] = {"all": dict(Counter(str(r["chosen_b"]) for r in recs)),
                                "per_target": {t: dict(Counter(str(r["chosen_b"]) for r in recs if r["target"] == t))
                                               for t in sorted({r["target"] for r in recs})}}
    (ROOT / "results" / "phaseC" / f"plmemb_report_{pool}.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["per_mask"], indent=1))
    print(json.dumps(out["b_hist"], indent=1))
    return out


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "report":
    report(tuple(sys.argv[3].split(",")) if len(sys.argv) > 3 else ("plmemb", "plmemb_base", "plmemb_gauss"),
           sys.argv[2] if len(sys.argv) > 2 else "fast")
