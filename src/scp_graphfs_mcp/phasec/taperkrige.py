"""Phase-C idea `taperkrige` (queen: biological priors, round 6).

Biology-tapered in-domain kriging: EnKF-style covariance localisation (Houtekamer & Mitchell 2001; Gaspari & Cohn 1999)
of the target's 30-cell base-completed covariance S (lwkrige v2) with a PSD heat kernel on a STRING + HPA graph, as an
exact correction on the cshrink base.

  0. psi = lwkrige.select_psi(M, ctx); (Yc, mu, sd, Vad) = lwkrige.base_fit(M, ctx, psi, iters=10) (= cshrink.impute).
  1. B = lwkrige.block(M) (<= 1,500 proteins); Z = (M[:, B] - mu[B]) / sd[B]; Sr = Vad_B Vad_B^T + 5 I;
     S = lwkrige.completed_cov(Z, Sr); c = tr(S) / tr(Sr); z1 = cond_mean(Z, c Sr).
  2. A on B = A_STRING (combined >= 700, weight combined / 1000, max over the group's accessions)
     + A_HPA (H diag(1 / n_l) H^T, diagonal zeroed; H = multi-hot HPA main locations).
     L = I - D^-1/2 A D^-1/2 (isolated rows 0); K_t = Q exp(-t Lambda) Q^T rescaled to unit diagonal; T_t = K_t.
  3. S_t = T_t o S; Sigma = lwkrige._blend(S_t, c Sr, delta); zhat = cond_mean(Z, Sigma);
     Y = Yc; Y[miss in B] += sd_B (zhat - z1). delta = 1 -> Yc.
  4. (t, delta) selected on tgtshrink.selection_mask (50/50 MCAR/MNAR inner mask, rng seed + 31337) from
     T_GRID x D_GRID plus 'base'; ties: base, then larger delta, then smaller t. Refit on M.

Entry points: impute (STRING+HPA graph), impute_null (same A, protein labels permuted within 10 weighted-degree
quantile bins, rng(seed + 7171)), impute_knn (k = 20 |corr| neighbours of c Sr, data control), impute_diag (T = I),
impute_base (delta = 1 forced, = cshrink).
Data: data/derived/string_v12.5_9606_uniprot.tsv.gz, HPA main locations via seqprior.protein_table(). No downloads.
Cache / logs: data/cache/phasec/taperkrige/.
"""

from __future__ import annotations

import fcntl
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

from scp_graphfs_mcp.phasec import tgtshrink as TS
from scp_graphfs_mcp.phasec.lwkrige import (_blend, _pref, _sd_s, base_fit, block, completed_cov,  # noqa: F401
                                            cond_mean, main_sigma, select_psi)
from scp_graphfs_mcp.phasec.ttprior import META

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "taperkrige"
STRING = ROOT / "data" / "derived" / "string_v12.5_9606_uniprot.tsv.gz"
STRING_SHA = "a933a71db0671f11d8ef86e0084ccfb3f375920a834e7cd3f3588914b47f1a78"
MIN_COMBINED = 700
T_GRID = (0.5, 1.0, 2.0)
D_GRID = (0.5, 0.7, 0.85, 0.95, 0.99)
K_NN = 20
N_BINS = 10
TMT = ("leduc2022", "gregoire2023_mixCTRL", "px063590_ADbrain")


# ---------------------------------------------------------------------------------------------------------------
# step 2: biological graph (built from protein names only, cached per target)
def _graph_all(target, names):
    """STRING edges (i, j, w) among all target proteins and the HPA multi-hot matrix H (p x L), cached."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{target}_graph.npz"
    with open(CACHE / f".{target}_graph.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if f.exists():
            d = np.load(f, allow_pickle=False)
            if int(d["p"]) == len(names):
                return d["ei"], d["ej"], d["ew"], d["H"]
        import pandas as pd

        from scp_graphfs_mcp.ids import group_accessions
        from scp_graphfs_mcp.phasec.seqprior import _strip, protein_table

        groups = [[_strip(a) for a in (group_accessions(n) or [])] for n in names]
        acc2idx = {}
        for j, g in enumerate(groups):
            for a in set(g):
                if a:
                    acc2idx.setdefault(a, []).append(j)
        df = pd.read_csv(STRING, sep="\t", usecols=["a", "b", "combined"])
        df = df[df["combined"] >= MIN_COMBINED]
        df = df[df["a"].isin(acc2idx) & df["b"].isin(acc2idx)]
        W = {}
        for a, b, w in zip(df["a"].tolist(), df["b"].tolist(), (df["combined"] / 1000.0).tolist()):
            for i in acc2idx[a]:
                for j in acc2idx[b]:
                    if i == j:
                        continue
                    k = (min(i, j), max(i, j))
                    if W.get(k, 0.0) < w:
                        W[k] = w
        ei = np.array([k[0] for k in W], dtype=np.int64)
        ej = np.array([k[1] for k in W], dtype=np.int64)
        ew = np.array(list(W.values()), dtype=float)
        tab = protein_table()
        S, Hd = tab["seq"], tab["hpa"]
        locs = sorted({x for v in Hd.values() for x in v})
        li = {x: k for k, x in enumerate(locs)}
        H = np.zeros((len(names), len(locs)), dtype=np.uint8)
        for j, g in enumerate(groups):
            for a in g:
                if a in S and S[a][0] in Hd:
                    for x in Hd[S[a][0]]:
                        H[j, li[x]] = 1
        tmp = CACHE / f".{target}_graph.tmp{os.getpid()}.npz"
        np.savez_compressed(tmp, ei=ei, ej=ej, ew=ew, H=H, p=len(names), locs=np.array(locs, dtype=str))
        os.replace(tmp, f)
        return ei, ej, ew, H


def bio_adj(B, ctx):
    """Dense |B| x |B| A = A_STRING + A_HPA and coverage diagnostics."""
    ei, ej, ew, H = _graph_all(ctx["target"], ctx["names"])
    pos = -np.ones(len(ctx["names"]), dtype=np.int64)
    pos[B] = np.arange(len(B))
    keep = (pos[ei] >= 0) & (pos[ej] >= 0)
    a, b, w = pos[ei[keep]], pos[ej[keep]], ew[keep]
    AS = np.zeros((len(B), len(B)))
    AS[a, b] = w
    AS[b, a] = w
    Hb = H[B].astype(float)
    n_l = Hb.sum(0)
    AH = (Hb / np.where(n_l > 0, n_l, 1.0)) @ Hb.T
    np.fill_diagonal(AH, 0.0)
    diag = {"frac_string": float((AS > 0).any(1).mean()), "frac_hpa": float((Hb.sum(1) > 0).mean()),
            "n_string_edges": int(keep.sum())}
    return AS + AH, diag


def permute_adj(A, seed):
    """Protein labels permuted within N_BINS quantile bins of weighted degree (isomorphic graph)."""
    deg = A.sum(1)
    rng = np.random.default_rng(seed + 7171)
    order = np.argsort(deg, kind="stable")
    perm = np.arange(len(A))
    for chunk in np.array_split(order, N_BINS):
        perm[chunk] = chunk[rng.permutation(len(chunk))]
    return A[np.ix_(perm, perm)]


def knn_adj(Sr, k=K_NN):
    d = np.sqrt(np.diag(Sr))
    R = np.abs(Sr / np.outer(d, d))
    np.fill_diagonal(R, -1.0)
    idx = np.argpartition(-R, k, axis=1)[:, :k]
    W = np.zeros_like(R)
    rows = np.repeat(np.arange(len(R)), k)
    W[rows, idx.ravel()] = R[rows, idx.ravel()]
    return np.maximum(W, W.T)


def heat_tapers(A, t_grid):
    """{t: K_t} with K_t = exp(-t L) (normalised Laplacian, isolated rows 0), rescaled to unit diagonal."""
    A = 0.5 * (A + A.T)
    np.fill_diagonal(A, 0.0)
    dg = A.sum(1)
    di = np.where(dg > 0, 1.0 / np.sqrt(np.where(dg > 0, dg, 1.0)), 0.0)
    L = np.diag((dg > 0).astype(float)) - di[:, None] * A * di[None, :]
    lam, Q = np.linalg.eigh(0.5 * (L + L.T))
    out = {}
    for t in t_grid:
        K = (Q * np.exp(-t * lam)) @ Q.T
        s = np.sqrt(np.clip(np.diag(K), 1e-300, None))
        K = K / np.outer(s, s)
        out[t] = 0.5 * (K + K.T)
    return out


# ---------------------------------------------------------------------------------------------------------------
# steps 1 and 3
def taper_correct(M, ctx, Yc, mu, sd, Vad, keys, adj_fn, ref_untapered=False):
    """Ys[k] for k in keys ('base' or (t, delta)); adj_fn(B, ctx, Sr) -> (A or None for T = I, diag)."""
    B = block(M)
    sds = _sd_s(sd)
    Z = (M[:, B] - np.nan_to_num(mu)[B]) / sds[B]
    Sr = main_sigma(Vad, B)
    S, _ = completed_cov(Z, Sr)
    c = float(np.trace(S) / np.trace(Sr))
    cSr = c * Sr
    z1 = cond_mean(Z, cSr)
    miss = np.isnan(Z)
    A, gdiag = adj_fn(B, ctx, cSr)
    ts = sorted({k[0] for k in keys if k != "base" and k[1] < 1})
    if A is None:
        tapers = {t: np.eye(len(B)) for t in ts}
    else:
        tapers = heat_tapers(A, ts) if ts else {}

    def corr(zd):
        Y = Yc.copy()
        sub = Y[:, B]
        sub[miss] = (Yc[:, B] + sds[B][None, :] * (zd - z1))[miss]
        Y[:, B] = sub
        return Y

    Ys = {}
    for k in keys:
        if k == "base" or k[1] >= 1:
            Ys[k] = Yc
            continue
        t, d = k
        Sd, _ = _blend(tapers[t] * S, cSr, d)
        Ys[k] = corr(cond_mean(Z, Sd))
    if ref_untapered:
        Sd, _ = _blend(S, cSr, 0.99)
        Ys["untapered_0.99"] = corr(cond_mean(Z, Sd))
    off = {}
    for t, T in tapers.items():
        off[str(t)] = float((np.abs(T).sum(1) - np.abs(np.diag(T))).mean())
    info = {"c": c, "nB": int(len(B)), "offdiag_mass": off, **gdiag}
    return Ys, info


# ---------------------------------------------------------------------------------------------------------------
# step 4: selection on the 50/50 inner mask, refit on M
def _log(variant, rec):
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / f"log_{variant}.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def _key_pref(k):
    return (0, 0.0, 0.0) if k == "base" else (1, -float(k[1]), float(k[0]))


def _kstr(k):
    return "base" if k == "base" else f"t{k[0]}_d{k[1]}"


def taper_run(M, ctx, variant, adj_fn, t_grid=T_GRID, d_grid=D_GRID, forced=None, log=True):
    seed = ctx["seed"]
    t0 = time.time()
    psi = select_psi(M, ctx)
    keys = ["base"] + [(t, d) for t in t_grid for d in d_grid]
    scores = {}
    sel_info = {}
    if forced is None:
        Msel, (r1, c1), (r2, c2) = TS.selection_mask(M, seed)
        Yc_s, mu_s, sd_s, Vad_s = base_fit(Msel, ctx, psi, iters=5)
        Ys_s, sel_info = taper_correct(Msel, ctx, Yc_s, mu_s, sd_s, Vad_s, keys, adj_fn, ref_untapered=log)
        for k, Y in Ys_s.items():
            e1 = np.sqrt(np.mean((Y[r1, c1] - M[r1, c1]) ** 2))
            e2 = np.sqrt(np.mean((Y[r2, c2] - M[r2, c2]) ** 2))
            scores[k] = float(0.5 * (e1 + e2))
        chosen = min(keys, key=lambda k: (scores[k], _key_pref(k)))
    else:
        chosen = forced
    Yc, mu, sd, Vad = base_fit(M, ctx, psi, iters=10)
    if chosen == "base" or chosen[1] >= 1:
        Y, info = Yc, {}
    else:
        Ys, info = taper_correct(M, ctx, Yc, mu, sd, Vad, [chosen], adj_fn)
        Y = Ys[chosen]
    if log:
        rec = {"target": ctx.get("target"), "acq": META.get(ctx.get("target"), {}).get("acq"), "seed": seed,
               "n_missing": int(np.isnan(M).sum()), "psi": "inf" if np.isinf(psi) else psi,
               "chosen": _kstr(chosen), "chosen_t": None if chosen == "base" else chosen[0],
               "chosen_delta": 1.0 if chosen == "base" else chosen[1],
               "scores": {(k if isinstance(k, str) else _kstr(k)): v for k, v in scores.items()},
               "score_untapered_0.99": scores.get("untapered_0.99"), "sel_info": sel_info, "fit_info": info,
               "seconds": round(time.time() - t0, 1)}
        _log(variant, rec)
    return Y


def _bio(B, ctx, cSr):
    return bio_adj(B, ctx)


def _bio_null(B, ctx, cSr):
    A, dg = bio_adj(B, ctx)
    return permute_adj(A, ctx["seed"]), dg


def _knn(B, ctx, cSr):
    A = knn_adj(cSr)
    return A, {"frac_string": None, "frac_hpa": None}


def _diag(B, ctx, cSr):
    return None, {}


def impute(Xm, ctx):
    return taper_run(Xm, ctx, "main", _bio)


def impute_null(Xm, ctx):
    return taper_run(Xm, ctx, "null", _bio_null)


def impute_knn(Xm, ctx):
    return taper_run(Xm, ctx, "knn", _knn)


def impute_diag(Xm, ctx):
    return taper_run(Xm, ctx, "diag", _diag, t_grid=(1.0,))


def impute_base(Xm, ctx):
    return taper_run(Xm, ctx, "base", _bio, forced="base")


# v2 (the one allowed revision, a scale fix): v1 chose t = 2, the top of T_GRID, in 37/37 tapered calls, and at t = 2
# the unit-diagonal heat kernel keeps only ~4.6 units of off-diagonal mass per protein (out of ~1,500; the untapered
# lwkrige v2 S keeps all of it). On petrosius2023 / leduc2022 the mass is ~27 at t = 4, ~230 at t = 8, ~930 at t = 16
# and ~1,390 at t = 32 (scratch diagnostic). v2 spans t in {2, 4, 8, 16}; delta 0.5 / 0.7 were never chosen in v1.
T_GRID_V2 = (2.0, 4.0, 8.0, 16.0)
D_GRID_V2 = (0.85, 0.9, 0.95, 0.98, 0.99)


def impute_v2(Xm, ctx):
    return taper_run(Xm, ctx, "v2", _bio, t_grid=T_GRID_V2, d_grid=D_GRID_V2)


def impute_v2_null(Xm, ctx):
    return taper_run(Xm, ctx, "v2_null", _bio_null, t_grid=T_GRID_V2, d_grid=D_GRID_V2)


def impute_v2_knn(Xm, ctx):
    return taper_run(Xm, ctx, "v2_knn", _knn, t_grid=T_GRID_V2, d_grid=D_GRID_V2)


# ---------------------------------------------------------------------------------------------------------------
# unit checks
def unit_check(target="petrosius2023"):
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H
    from scp_graphfs_mcp import phase2 as P
    from scp_graphfs_mcp.phasec import cshrink as C
    from scp_graphfs_mcp.phasec import nfbase as NF

    out = {}
    Xs, Xm, mask, names = NF._split(target, "MCAR", 100)
    V, warm = H.target_models(target, names)
    b = P._fill(Xm, P.loo_knn(Xm))
    b = np.where(np.isnan(b), np.nanmean(Xm), b)
    ctx = {"V": V, "warm": warm, "base": b, "seed": 100, "target": target, "names": names}
    ref = C.impute(Xm, ctx)
    yb = impute_base(Xm, dict(ctx))
    out["U1_base_vs_cshrink_maxdiff"] = float(np.max(np.abs(yb - ref)))
    psi = select_psi(Xm, ctx)
    Yc, mu, sd, Vad = base_fit(Xm, ctx, psi)
    B = block(Xm)
    Z = (Xm[:, B] - mu[B]) / _sd_s(sd)[B]
    Sr = main_sigma(Vad, B)
    t0 = time.time()
    S, _ = completed_cov(Z, Sr)
    out["T_cond_mean_s"] = round(time.time() - t0, 2)
    A, dg = bio_adj(B, ctx)
    out.update(dg)
    Ks = heat_tapers(A, T_GRID)
    for t, K in Ks.items():
        w = np.linalg.eigvalsh(K)
        ws = np.linalg.eigvalsh(K * S)
        out[f"U2_K{t}_min_eig"] = float(w.min())
        out[f"U2_KoS{t}_min_eig_rel"] = float(ws.min() / ws.max())
        out[f"U2_K{t}_diag_maxdev"] = float(np.max(np.abs(np.diag(K) - 1)))
        out[f"offdiag_mass_t{t}"] = float((np.abs(K).sum(1) - 1).mean())
    An = permute_adj(A, 100)
    out["U3_null_same_spectrum"] = bool(np.allclose(np.sort(An.sum(1)), np.sort(A.sum(1))))
    # delta = 1 key returns Yc exactly
    Ys, _ = taper_correct(Xm, ctx, Yc, mu, sd, Vad, ["base", (1.0, 1.0)], _bio)
    out["U4_delta1_is_base"] = bool(np.array_equal(Ys[(1.0, 1.0)], Yc))
    out["pass"] = bool(out["U1_base_vs_cshrink_maxdiff"] < 1e-10 and out["U3_null_same_spectrum"]
                       and out["U4_delta1_is_base"] and min(out[f"U2_K{t}_min_eig"] for t in T_GRID) > -1e-8)
    out["target"] = target
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "unit_check.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    if "--unit" in sys.argv:
        unit_check()
