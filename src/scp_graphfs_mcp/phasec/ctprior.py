"""Phase-C idea `ctprior` (queen: biological priors).

HPA-marker cell-type signature scores as soft per-CELL weights on the reference in a weighted fit_loadings.
  1. Markers (per target): HPA single-cell-type RNA (pivot.npz, cached by scrnaprior), L = log2(nCPM + 1); each target
     protein -> Ensembl genes (scrnaprior._acc2genes: proteinatlas.tsv first, HGNC fallback), RNA row = mean over genes;
     e_jk = L_jk - median_k' L_jk'; markers of type k = proteins with e_jk >= 2, top 50 by e_jk; keep types with >= 10.
  2. Raw aligned reference: corpus_reference WITHOUT the per-dataset z-score (same datasets, same order; asserted).
  3. Within-cell percentile rank r_ij = rank of x_ij among the observed entries of cell i / n_obs_i.
  4. Abundance-matched background (AddModuleScore style): a_j = mean r_ij over reference rows, 24 quantile bins, 5
     controls per marker from its bin (excluding that type's markers), rng default_rng(0).
  5. s_ik = mean r(observed markers) - mean r(observed controls), NaN if < 3 markers observed; standardised per type by
     the reference mean/SD.
  6. Per call: target scores from Xm (reference mean/SD), sbar_k (>= 5 finite cells), d_i = K_used * nanmean_k
     (s_ik - sbar_k)^2 (rows with < 10 finite types -> median d), w_i = exp(-d_i / tau).
     impute: tau = median(d); impute_half: tau = median(d)/2; impute_inf: w = 1.
  7. Vw = ttprior.fit_loadings_weighted(Z, w, V, warm, iters=3, warm_start=True); memoised by sha1(w).
  8. transfer.adapt_project(Xm, Vw, warm, base, lam_v=100, extra=10).
Null (impute_null*): every type's marker set is replaced by random non-marker proteins drawn one-for-one from the same
reference-observation-rate decile (and the same abundance bin where possible); controls redrawn by the same rule; the
rng is default_rng(12345 + index of target in FAST). Everything else identical.
Cache: data/cache/phasec/ctprior/<target>_markers.npz, <target>_S.npz. Precompute (serial, before screening):
    python -m scp_graphfs_mcp.phasec.ctprior --precompute [--targets a,b]
Diagnostics on the harness's eval inputs:  python -m scp_graphfs_mcp.phasec.ctprior --diag [--targets a,b]
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

from scp_graphfs_mcp import transfer as T

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "ctprior"
PIVOT = ROOT / "data" / "cache" / "phasec" / "scrnaprior" / "pivot.npz"
FAST = ["petrosius2023", "leduc2022", "weiss2026_liver", "gregoire2023_mixCTRL", "px063590_ADbrain", "sadiku2025"]
E_MIN, TOP, MIN_MARKERS = 2.0, 50, 10
N_BINS, N_CTRL, MIN_OBS_MARK = 24, 5, 3
MIN_CELLS_SBAR, MIN_TYPES_ROW = 5, 10
ITERS_WARM, LAM_V, EXTRA = 3, 100.0, 10

_MEM: dict = {}


# ---------------------------------------------------------------------------------------------------------------
def _atomic_savez(f, **kw):
    tmp = f.with_suffix(f".tmp{os.getpid()}.npz")
    np.savez_compressed(tmp, **kw)
    os.replace(tmp, f)


def within_cell_rank(X, chunk=2000):
    """r_ij = rank of x_ij among the observed entries of row i / n_obs_i (NaN where unobserved). float32."""
    out = np.full(X.shape, np.nan, np.float32)
    for a in range(0, X.shape[0], chunk):
        B = X[a:a + chunk]
        obs = ~np.isnan(B)
        Bi = np.where(obs, B, np.inf)
        rk = np.argsort(np.argsort(Bi, axis=1, kind="stable"), axis=1, kind="stable").astype(np.float32) + 1.0
        n = obs.sum(1, keepdims=True).astype(np.float32)
        out[a:a + chunk] = np.where(obs, rk / np.maximum(n, 1.0), np.nan)
    return out


def raw_reference(target, names):
    """corpus_reference(names, (target,), [d]) for each reference dataset, but WITHOUT _z. Returns (X, src)."""
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from explore_phaseb import reference_datasets

    from scp_graphfs_mcp.external import _load

    tacc = T._target_accs(names)
    parts, src = [], []
    for d in reference_datasets(exclude=(target,)):
        X, accs = _load(d)
        pos = {a: i for i, a in enumerate(accs)}
        B = np.full((X.shape[0], len(tacc)), np.nan)
        cols = [(j, pos[a]) for j, a in enumerate(tacc) if a in pos]
        if cols:
            tj, sj = zip(*cols)
            B[:, list(tj)] = X[:, list(sj)]
        parts.append(B)
        src += [d] * len(B)
    return np.vstack(parts), np.array(src)


# ---------------------------------------------------------------------------------------------------------------
def build_markers(names):
    """(M bool p x K, types list, e_markers) for types with >= MIN_MARKERS in-target markers."""
    from scp_graphfs_mcp.ids import group_accessions
    from scp_graphfs_mcp.phasec.scrnaprior import _acc2genes

    d = np.load(PIVOT)
    genes, types = d["genes"], d["celltypes"]
    L = np.log2(d["nCPM"].astype(float) + 1.0)
    gi = {g: i for i, g in enumerate(genes)}
    acc, fb = _acc2genes()
    p = len(names)
    R = np.full((p, L.shape[1]), np.nan)
    for j, nm in enumerate(names):
        accs = group_accessions(str(nm))
        hit = {g for a in accs for g in acc.get(a, ()) if g in gi}
        if not hit:
            hit = {g for a in accs for g in fb.get(a, ()) if g in gi}
        if hit:
            R[j] = L[[gi[g] for g in sorted(hit)]].mean(0)
    E = R - np.nanmedian(R, axis=1, keepdims=True)
    M = np.zeros((p, L.shape[1]), bool)
    for k in range(L.shape[1]):
        e = np.where(np.isfinite(E[:, k]), E[:, k], -np.inf)
        cand = np.where(e >= E_MIN)[0]
        cand = cand[np.argsort(-e[cand], kind="stable")][:TOP]
        M[cand, k] = True
    keep = M.sum(0) >= MIN_MARKERS
    return M[:, keep], [str(t) for t in types[keep]], int(np.isfinite(R[:, 0]).sum())


def draw_controls(M, bins, valid, rng):
    """C (p x K float counts): for each marker of type k, N_CTRL proteins from its abundance bin, excluding type-k
    markers (and proteins with no abundance, `valid`=False). With replacement across markers (AddModuleScore)."""
    p, K = M.shape
    C = np.zeros((p, K))
    for k in range(K):
        pool_ok = valid & ~M[:, k]
        for j in np.where(M[:, k])[0]:
            if not valid[j]:
                continue
            pool = np.where(pool_ok & (bins == bins[j]))[0]
            if len(pool) == 0:
                continue
            C[rng.choice(pool, N_CTRL, replace=len(pool) < N_CTRL), k] += 1
    return C


def null_markers(M, bins, valid, obsrate, rng):
    """Each real marker -> a random non-marker (of that type, not yet chosen) from the same observation-rate decile,
    and from the same abundance bin where possible. Markers with no abundance (never observed in the reference) are
    replaced by a random never-observed non-marker protein."""
    dec_edges = np.quantile(obsrate[valid], np.linspace(0, 1, 11)[1:-1]) if valid.any() else []
    dec = np.digitize(obsrate, dec_edges)
    p, K = M.shape
    Mn = np.zeros_like(M)
    n_bin_match = n_tot = 0
    for k in range(K):
        taken = M[:, k].copy()
        for j in np.where(M[:, k])[0]:
            n_tot += 1
            if valid[j]:
                pool = np.where(~taken & valid & (dec == dec[j]) & (bins == bins[j]))[0]
                if len(pool):
                    n_bin_match += 1
                else:
                    pool = np.where(~taken & valid & (dec == dec[j]))[0]
                if len(pool) == 0:
                    pool = np.where(~taken & valid)[0]
            else:
                pool = np.where(~taken & ~valid)[0]
                if len(pool) == 0:
                    pool = np.where(~taken)[0]
            q = rng.choice(pool)
            Mn[q, k] = True
            taken[q] = True
    return Mn, n_bin_match / max(n_tot, 1)


def scores(Rk, M, C):
    """s_ik from within-cell ranks Rk (n x p, NaN unobserved)."""
    obs = ~np.isnan(Rk)
    R0 = np.where(obs, Rk, 0.0).astype(float)
    O = obs.astype(float)
    Mf = M.astype(float)
    nm = O @ Mf
    sm = R0 @ Mf
    nc = O @ C
    sc = R0 @ C
    with np.errstate(invalid="ignore", divide="ignore"):
        s = sm / nm - sc / nc
    s[(nm < MIN_OBS_MARK) | (nc <= 0)] = np.nan
    return s


def precompute_target(target, names=None, src_check=None):
    """Markers, controls, null sets and reference S / S_null for one target (cached)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{target}_S.npz"
    if f.exists():
        return dict(np.load(f, allow_pickle=False))
    if names is None:
        sys.path.insert(0, str(ROOT / "benchmarks"))
        from explore_phaseb import load_target

        _, names, _ = load_target(target)
    M, types, n_mapped = build_markers(names)
    Xr, src = raw_reference(target, names)
    if src_check is not None:
        assert np.array_equal(np.asarray(src_check), src), "raw reference rows differ from ctx['ref_src']()"
    Rk = within_cell_rank(Xr)
    del Xr
    with np.errstate(invalid="ignore"):
        a = np.nanmean(Rk, 0)
    obsrate = (~np.isnan(Rk)).mean(0)
    valid = np.isfinite(a)
    edges = np.quantile(a[valid], np.linspace(0, 1, N_BINS + 1)[1:-1])
    bins = np.where(valid, np.digitize(np.nan_to_num(a), edges), -1)
    C = draw_controls(M, bins, valid, np.random.default_rng(0))
    ti = FAST.index(target) if target in FAST else len(FAST) + sorted(FAST + [target]).index(target)
    rng_n = np.random.default_rng(12345 + ti)
    Mn, frac_bin = null_markers(M, bins, valid, obsrate, rng_n)
    Cn = draw_controls(Mn, bins, valid, rng_n)
    out = {"types": np.array(types, dtype="U"), "M": M, "C": C, "Mn": Mn, "Cn": Cn, "src": src}
    for tag, (MM, CC) in {"": (M, C), "_null": (Mn, Cn)}.items():
        S = scores(Rk, MM, CC)
        with np.errstate(invalid="ignore"):
            mu, sd = np.nanmean(S, 0), np.nanstd(S, 0)
        sd = np.where(np.isfinite(sd) & (sd > 1e-9), sd, 1.0)
        out["S" + tag] = ((S - mu) / sd).astype(np.float32)
        out["mu" + tag], out["sd" + tag] = mu, sd
    _atomic_savez(CACHE / f"{target}_markers.npz", types=out["types"], M=M, Mn=Mn, n_mapped=n_mapped)
    _atomic_savez(f, **out)
    meta = {"target": target, "n_proteins": len(names), "n_mapped_to_rna": n_mapped, "n_types_kept": len(types),
            "markers_per_type_median": float(np.median(M.sum(0))),
            "observed_markers_per_ref_cell_median": float(np.median((~np.isnan(Rk)).astype(float) @ M.astype(float))),
            "null_frac_markers_bin_matched": frac_bin, "N_ref_rows": int(len(src))}
    (CACHE / f"{target}_precompute.json").write_text(json.dumps(meta, indent=1))
    return out


def _cached(target, ctx=None):
    if _MEM.get("S_target") != target:
        _MEM["S"] = precompute_target(target, ctx["names"] if ctx else None,
                                      ctx["ref_src"]() if ctx else None)
        _MEM["S_target"] = target
    return _MEM["S"]


# ---------------------------------------------------------------------------------------------------------------
def target_scores(Xm, d, null=False):
    tag = "_null" if null else ""
    M, C = (d["Mn"], d["Cn"]) if null else (d["M"], d["C"])
    s = scores(within_cell_rank(np.asarray(Xm, float)), M, C)
    return (s - d["mu" + tag]) / d["sd" + tag]


def distances(Sref, St):
    """d_i for each reference row; K_used = types with >= MIN_CELLS_SBAR finite target scores."""
    nfin = np.isfinite(St).sum(0)
    used = nfin >= MIN_CELLS_SBAR
    with np.errstate(invalid="ignore"):
        sbar = np.where(used, np.nanmean(np.where(used[None, :], St, np.nan), 0), np.nan)
    Ku = int(used.sum())
    D = (Sref[:, used].astype(float) - sbar[used][None, :]) ** 2
    nf = np.isfinite(D).sum(1)
    with np.errstate(invalid="ignore"):
        dd = Ku * np.nanmean(D, 1)
    bad = nf < MIN_TYPES_ROW
    med = np.nanmedian(dd[~bad]) if (~bad).any() else 1.0
    dd[bad] = med
    return dd, sbar, Ku, int(bad.sum())


def weights(Xm, ctx_or_target, mode, null=False, d=None):
    if mode == "inf":
        d = d if d is not None else _cached(ctx_or_target["target"], ctx_or_target)
        return np.ones(len(d["src"])), None
    d = d if d is not None else _cached(ctx_or_target["target"], ctx_or_target)
    St = target_scores(Xm, d, null)
    dd, sbar, Ku, nbad = distances(d["S_null" if null else "S"], St)
    tau = np.median(dd) * (0.5 if mode == "half" else 1.0)
    w = np.exp(-dd / max(tau, 1e-12))
    return w, {"sbar": sbar, "K_used": Ku, "n_rows_median_d": nbad, "tau": float(tau)}


def _ref_Z(ctx):
    if _MEM.get("Z_target") != ctx["target"]:
        _MEM.pop("Z", None)
        _MEM["Z"] = ctx["ref"]()
        _MEM["Z_target"] = ctx["target"]
        _MEM["V"] = {}
    return _MEM["Z"]


def get_Vw(ctx, w):
    from scp_graphfs_mcp.phasec.ttprior import fit_loadings_weighted

    Z = _ref_Z(ctx)
    key = hashlib.sha1(np.ascontiguousarray(w, dtype=np.float64).tobytes()).hexdigest()
    if key not in _MEM["V"]:
        Vw, _, _ = fit_loadings_weighted(Z, w, ctx["V"], ctx["warm"], iters=ITERS_WARM, warm_start=True)
        _MEM["V"][key] = Vw
    return _MEM["V"][key]


def _impute(Xm, ctx, mode, null=False):
    w, _ = weights(Xm, ctx, mode, null)
    Vw = get_Vw(ctx, w)
    return T.adapt_project(Xm, Vw, ctx["warm"], ctx["base"], lam_v=LAM_V, extra=EXTRA)


def impute(Xm, ctx):
    return _impute(Xm, ctx, "med", False)


def impute_null(Xm, ctx):
    return _impute(Xm, ctx, "med", True)


def impute_half(Xm, ctx):
    return _impute(Xm, ctx, "half", False)


def impute_null_half(Xm, ctx):
    return _impute(Xm, ctx, "half", True)


def impute_inf(Xm, ctx):
    return _impute(Xm, ctx, "inf", False)


# ---------------------------------------------------------------------------------------------------------------
def eval_inputs(target):
    """The harness's evaluation inputs Xm (same subsample + mask as benchmarks/phasec.run_split), FAST seeds."""
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from explore_phaseb import load_target, subsample
    from run import make_mask

    X, names, labels = load_target(target)
    out = []
    for kind in ("MCAR", "MNAR"):
        for seed in (100, 101):
            idx = subsample(X.shape[0], labels, 30, seed)
            Xs = X[idx]
            eligible = ((~np.isnan(Xs)).mean(0) >= 0.5) & ~np.isnan(Xs).all(0)
            mask = make_mask(Xs, eligible, kind, seed)
            Xm = Xs.copy()
            Xm[mask] = np.nan
            out.append((kind, seed, Xm))
    return out, names


def diag(target):
    from scipy.stats import spearmanr

    d = precompute_target(target)
    src = d["src"]
    types = [str(t) for t in d["types"]]
    # rows fit_loadings keeps (same rule as the weighted fit): >= rank observed warm proteins in the z-scored ref
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H
    from explore_phaseb import load_target

    _, names, _ = load_target(target)
    _, warm = H.target_models(target, names)
    Xr, _ = raw_reference(target, names)
    keep = (~np.isnan(Xr[:, warm])).sum(1) >= 20
    del Xr
    datasets = sorted(set(src.tolist()))
    n_unw = {ds: int(((src == ds) & keep).sum()) for ds in datasets}
    res = {"target": target, "N_kept": int(keep.sum()), "K_types": len(types), "per_call": []}
    acc = {m: {"neff": [], "share": [], "rho": [], "sbar": []} for m in ("med", "half")}
    for kind, seed, Xm in eval_inputs(target)[0]:
        for mode in ("med", "half"):
            w, info = weights(Xm, None, mode, False, d=d)
            wn, info_n = weights(Xm, None, mode, True, d=d)
            wk, wnk = w[keep], wn[keep]
            wk = wk * len(wk) / wk.sum()
            wnk = wnk * len(wnk) / wnk.sum()
            neff = float(wk.sum() ** 2 / (wk ** 2).sum())
            neff_n = float(wnk.sum() ** 2 / (wnk ** 2).sum())
            share = {ds: float(wk[src[keep] == ds].sum() / wk.sum()) for ds in datasets}
            rho = float(spearmanr(w[keep], wn[keep]).statistic)
            acc[mode]["neff"].append((neff, neff_n))
            acc[mode]["share"].append(share)
            acc[mode]["rho"].append(rho)
            acc[mode]["sbar"].append(info["sbar"])
            res["per_call"].append({"mask": kind, "seed": seed, "mode": mode, "n_eff": neff, "n_eff_null": neff_n,
                                    "spearman_w_true_vs_null": rho, "K_used": info["K_used"],
                                    "K_used_null": info_n["K_used"], "n_rows_median_d": info["n_rows_median_d"]})
    tot = sum(n_unw.values())
    unw = {ds: n_unw[ds] / tot for ds in datasets}
    for mode in ("med", "half"):
        a = acc[mode]
        sh = {ds: float(np.mean([s[ds] for s in a["share"]])) for ds in datasets}
        gain = sorted(datasets, key=lambda ds: sh[ds] - unw[ds], reverse=True)
        sbar = np.nanmean(np.array(a["sbar"]), 0)
        order = [types[i] for i in np.argsort(-np.nan_to_num(sbar, nan=-np.inf))]
        res[mode] = {"n_eff_mean": float(np.mean([x[0] for x in a["neff"]])),
                     "n_eff_null_mean": float(np.mean([x[1] for x in a["neff"]])),
                     "spearman_w_true_vs_null_mean": float(np.mean(a["rho"])),
                     "share_weighted": sh, "share_unweighted": unw,
                     "top5_gainers": [(ds, round(sh[ds] - unw[ds], 4)) for ds in gain[:5]],
                     "top5_losers": [(ds, round(sh[ds] - unw[ds], 4)) for ds in gain[::-1][:5]],
                     "sbar_top10": order[:10]}
        if "hepatocytes" in types:
            res[mode]["rank_hepatocytes_in_sbar"] = order.index("hepatocytes") + 1
    (CACHE / f"diag_{target}.json").write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    tg = sys.argv[sys.argv.index("--targets") + 1].split(",") if "--targets" in sys.argv else FAST
    if "--precompute" in sys.argv:
        import time

        sys.path.insert(0, str(ROOT / "benchmarks"))
        import phasec as H
        from explore_phaseb import load_target

        for t in tg:
            t0 = time.time()
            _, names, _ = load_target(t)
            _, src = H._ref_parts(t, names)
            precompute_target(t, names, src)
            print(t, f"{time.time() - t0:.1f}s", (CACHE / f"{t}_precompute.json").read_text(), flush=True)
    if "--diag" in sys.argv:
        for t in tg:
            r = diag(t)
            print(t, json.dumps({m: {k: r[m][k] for k in ("n_eff_mean", "n_eff_null_mean",
                                                            "spearman_w_true_vs_null_mean", "sbar_top10",
                                                            "top5_gainers")} for m in ("med",)}, indent=1),
                  "N_kept", r["N_kept"], flush=True)
