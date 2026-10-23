"""Phase-C idea `dca` (queen: innovation abduction).

Mean-field direct-coupling analysis (Morcos et al., PNAS 2011) on the leave-target-out single-cell reference:
  STEP 1 (per target and variant, cached): S = the P=1500 target proteins most observed in the REAL reference;
    pairwise-complete correlation C on Z[:, S] (pairs with n<30 -> 0), eigen-clip PSD repair, unit diagonal;
    shrinkage lambda chosen inside the reference (largest dataset held out, 300 rows, 20% hidden, grid
    {0.3,0.5,0.7,0.85}); C_lam = (1-lam) C + lam I; Theta = inv(C_lam).
    Real variant also exports APC-corrected partial-correlation couplings (not used for imputation).
  STEP 2 (per split): standardise Xm like transfer.adapt_project and impute each cell's missing S-proteins by the
    Gaussian conditional mean E[z_m | z_o] = -Theta_mm^-1 Theta_mo z_o (= C_mo C_oo^-1 z_o), clipped to [-4, 4].
    Everything outside S keeps ctx['champion'].
Null: identical pipeline, but before the statistics each reference dataset's protein columns are permuted with an
independent permutation (rng 12345 + index of the dataset in sorted order). S comes from the real reference.
No external data. Cache: data/cache/phasec/dca/{target}_{real|null}.npz, {target}_couplings.npz.
"""

from __future__ import annotations

import fcntl
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "dca"
P_TOP, MIN_PAIR, LAM_GRID, HOLD_ROWS, HIDE_FRAC, TOPK = 1500, 30, (0.3, 0.5, 0.7, 0.85), 300, 0.2, 30
_MEM: dict = {}


def _corr(Zs):
    """Pairwise-complete correlation with min-pair cutoff, PSD repair by eigenvalue clipping, unit diagonal."""
    W = np.isfinite(Zs).astype(float)
    Z0 = np.nan_to_num(Zs)
    n = W.T @ W
    C = (Z0.T @ Z0) / np.maximum(n, 1)
    C[n < MIN_PAIR] = 0.0
    d = np.sqrt(np.clip(np.diag(C), 1e-8, None))
    C = C / d[:, None] / d[None, :]
    np.fill_diagonal(C, 1.0)
    w, U = np.linalg.eigh(C)
    n_neg = int((w < 0).sum())
    w = np.clip(w, 1e-3 * w.mean(), None)
    C = (U * w) @ U.T
    d = np.sqrt(np.diag(C))
    C = C / d[:, None] / d[None, :]
    C = (C + C.T) / 2
    return C, n_neg


def _condmean(Clam, Theta, z, o, m):
    """E[z_m | z_o] under N(0, Clam); uses whichever system is smaller."""
    if o.sum() <= m.sum():
        return Clam[np.ix_(m, o)] @ np.linalg.solve(Clam[np.ix_(o, o)], z[o])
    return -np.linalg.solve(Theta[np.ix_(m, m)], Theta[np.ix_(m, o)] @ z[o])


def _select_lambda(Zs, src):
    ds, cnt = np.unique(src, return_counts=True)
    dstar = ds[np.argmax(cnt)]
    C, _ = _corr(Zs[src != dstar])
    rows = np.where(src == dstar)[0]
    rng = np.random.default_rng(0)
    if len(rows) > HOLD_ROWS:
        rows = np.sort(rng.choice(rows, HOLD_ROWS, replace=False))
    H = Zs[rows]
    obs = np.argwhere(np.isfinite(H))
    hide = obs[rng.choice(len(obs), int(HIDE_FRAC * len(obs)), replace=False)]
    Hm = H.copy()
    Hm[hide[:, 0], hide[:, 1]] = np.nan
    hmask = np.zeros(H.shape, bool)
    hmask[hide[:, 0], hide[:, 1]] = True
    I = np.eye(C.shape[0])
    scores = {}
    for lam in LAM_GRID:
        Clam = (1 - lam) * C + lam * I
        Theta = np.linalg.inv(Clam)
        se, k = 0.0, 0
        for i in range(len(H)):
            o = np.isfinite(Hm[i])
            m = hmask[i]
            if o.sum() < 5 or not m.any():
                continue
            zm = np.clip(_condmean(Clam, Theta, Hm[i], o, m), -4, 4)
            se += float(((zm - H[i, m]) ** 2).sum())
            k += int(m.sum())
        scores[lam] = float(np.sqrt(se / max(k, 1)))
    # reference RMSE of the column-mean (z=0) predictor on the same hidden entries, for scale
    scores["mean0"] = float(np.sqrt(np.mean(H[hmask] ** 2)))
    best = min(LAM_GRID, key=lambda l: scores[l])
    return best, scores, str(dstar)


def _export_couplings(target, Theta, S, names, tag="couplings"):
    dg = np.sqrt(np.diag(Theta))
    pc = -Theta / dg[:, None] / dg[None, :]
    np.fill_diagonal(pc, 0.0)
    A = np.abs(pc)
    P = len(S)
    rmean = A.sum(1) / (P - 1)
    amean = A.sum() / (P * (P - 1))
    F = A - np.outer(rmean, rmean) / amean
    np.fill_diagonal(F, -np.inf)
    top = np.argsort(-F, axis=1)[:, :TOPK]
    Ftop = np.take_along_axis(F, top, 1)
    acc = np.array([str(names[j]) for j in S])
    np.savez_compressed(CACHE / f"{target}_{tag}.npz", S=S, accession=acc, partners=top, partner_F=Ftop,
                        partner_accession=acc[top], degree=Ftop.sum(1))


def _build(target, ctx, variant):
    t0 = time.time()
    Z = ctx["ref"]()
    src = ctx["ref_src"]()
    nobs = np.isfinite(Z).sum(0)
    S = np.sort(np.argsort(-nobs, kind="stable")[:P_TOP])
    if variant == "null":
        Z = Z.copy()
        for k, d in enumerate(sorted(np.unique(src))):
            r = src == d
            pi = np.random.default_rng(12345 + k).permutation(Z.shape[1])
            Z[r] = Z[r][:, pi]
    Zs = Z[:, S]
    lam, scores, dstar = _select_lambda(Zs, src)
    C, n_neg = _corr(Zs)
    Clam = (1 - lam) * C + lam * np.eye(len(S))
    Theta = np.linalg.inv(Clam)
    info = {"target": target, "variant": variant, "lambda": lam, "lambda_scores": {str(k): v for k, v in scores.items()},
            "heldout_dataset": dstar, "n_ref_rows": int(Z.shape[0]), "n_ref_datasets": int(len(np.unique(src))),
            "p": int(Z.shape[1]), "P": int(len(S)), "n_neg_eig_before_repair": n_neg,
            "build_s": round(time.time() - t0, 1)}
    np.savez(CACHE / f"{target}_{variant}.npz", S=S, lam=lam, Theta=Theta, Clam=Clam, info=json.dumps(info))
    _export_couplings(target, Theta, S, ctx["names"], "couplings" if variant == "real" else "null_couplings")
    print(f"[dca] built {target} {variant}: {json.dumps(info)}", flush=True)


def _stats(ctx, variant):
    target = ctx["target"]
    key = (target, variant)
    if key in _MEM:
        return _MEM[key]
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{target}_{variant}.npz"
    with open(CACHE / f"{target}.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            if not f.exists():
                _build(target, ctx, variant)
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)
    d = np.load(f)
    if any(k[0] != target for k in _MEM):
        _MEM.clear()
    _MEM[key] = (d["S"], d["Clam"], d["Theta"])
    return _MEM[key]


def _impute(Xm, ctx, null: bool):
    S, Clam, Theta = _stats(ctx, "null" if null else "real")
    mu = np.nanmean(Xm, 0)
    sd = np.nanstd(Xm, 0)
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    ok = np.isfinite(mu)
    Zt = (Xm[:, S] - mu[S]) / sd[S]
    okS = ok[S]
    out = np.array(ctx["champion"], float, copy=True)
    for i in range(Xm.shape[0]):
        o = np.isfinite(Zt[i])
        m = ~o & okS
        if o.sum() < 5 or not m.any():
            continue
        zm = np.clip(_condmean(Clam, Theta, Zt[i], o, m), -4, 4)
        cols = S[m]
        out[i, cols] = mu[cols] + sd[cols] * zm
    return out


def impute(Xm, ctx):
    return _impute(Xm, ctx, False)


def impute_null(Xm, ctx):
    return _impute(Xm, ctx, True)


def couplings_diag(targets=("petrosius2023", "leduc2022", "weiss2026_liver", "gregoire2023_mixCTRL",
                            "px063590_ADbrain", "sadiku2025")):
    """Optional diagnostic (not a stack member): fraction of each protein's top-30 partners that share a CORUM
    complex, for APC-corrected DCA couplings (real, null), raw |correlation| (real) and random partners.
    Writes results/phaseC/dca_couplings_diag.json.  Run: python -m scp_graphfs_mcp.phasec.dca"""
    from scp_graphfs_mcp.complexes import _membership_index
    from scp_graphfs_mcp.ids import group_accessions

    idx = _membership_index(("corum",))
    out = {}
    for t in targets:
        if not (CACHE / f"{t}_couplings.npz").exists():
            continue
        r = np.load(CACHE / f"{t}_couplings.npz", allow_pickle=True)
        acc = r["accession"]
        sets = [set(c for a in group_accessions(str(f)) for c in idx.get(a, [])) for f in acc]

        def frac(top):
            hits = tot = 0
            for i in range(len(acc)):
                if sets[i]:
                    tot += len(top[i])
                    hits += sum(bool(sets[i] & sets[j]) for j in top[i])
            return round(hits / max(tot, 1), 4)

        d = np.load(CACHE / f"{t}_real.npz")
        A = np.abs(d["Clam"])
        np.fill_diagonal(A, -np.inf)
        rng = np.random.default_rng(0)
        n = len(acc)
        res = {"P": int(n), "n_in_corum": int(sum(bool(s) for s in sets)), "lambda": float(d["lam"]),
               "real_top30_F": frac(r["partners"]), "real_top30_abs_corr": frac(np.argsort(-A, 1)[:, :TOPK]),
               "random_top30": frac([rng.choice(np.delete(np.arange(n), i), TOPK, replace=False) for i in range(n)])}
        if (CACHE / f"{t}_null_couplings.npz").exists():
            res["null_top30_F"] = frac(np.load(CACHE / f"{t}_null_couplings.npz")["partners"])
        out[t] = res
    f = ROOT / "results" / "phaseC" / "dca_couplings_diag.json"
    f.write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(couplings_diag(), indent=1))
