"""Phase-C idea `mnarjoint` (queen: representation & mathematics, round 4).

Missingness as a second channel of a joint low-rank encoding (Sportisse, Boyer & Josse 2020, Stat. Comput.
30:1629-1643, doi:10.1007/s11222-020-09963-5): the cell factors U of the adapted projection must explain both the
standardised observed values Z_obs (ridge-anchored loadings V, as in the champion) and the column-centred missing
indicators w * (M - mean) (free loadings B, Gaussian / least-squares surrogate). The M block never produces values;
it moves U, and through U the imputations, entry by entry. No global shift.

Base (round-4 rule, cshrink not landed): tgtshrink.impute_null's pipeline = shrink(M, psi, "null", seed) with psi
selected on tgtshrink.selection_mask (50/50 MCAR/MNAR inner mask, rng seed + 31337), then
nfbase.adapt_project_mb(lam_v=100, extra=10, lam=5, iters=10, seed=0). w = 0 reproduces it bit for bit.

  adapt_project_mj : adapt_project_mb + M block; per iteration (a) U-step with B, (b) champion ridge V-step,
                     (c) B = Wm^T U (U^T U + lam_b I)^-1.
  _run             : psi as tgtshrink._run (mode "null"); w from {0} U W_GRID on the same two inner hole sets
                     (score = 0.5 (RMSE_MCAR-part + RMSE_MNAR-part), iters = 5; ties -> smaller w); final iters = 10.
Null (impute_null): the missing indicators fed to the M block are permuted within each cell across ok proteins
(rng default_rng(seed + 4242), rows drawn in order), both for the inner-selection fits and the final fit.

    python -m scp_graphfs_mcp.phasec.mnarjoint unit     -> results/phaseC/mnarjoint_unit.json
    python -m scp_graphfs_mcp.phasec.mnarjoint report   -> per-mask means + chosen-w distribution (stdout / json)
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from scp_graphfs_mcp.phasec.nfbase import EXTRA, LAM_V, adapt_project_mb
from scp_graphfs_mcp.phasec.tgtshrink import PSI_GRID, selection_mask, shrink

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "mnarjoint"
OUT = ROOT / "results" / "phaseC"
W_GRID = (0.1, 0.3, 1.0, 3.0)
LAM, LAM_B = 5.0, 5.0


# ---------------------------------------------------------------------------------------------------------------
# 1. joint [Z_obs | w (M - mean)] factorisation
def adapt_project_mj(Xm, Vref, usable, base, mu, sd, w, Mind=None, lam_v=100.0, lam=5.0, lam_b=5.0, extra=10,
                     iters=10, seed=0, return_U=False, force_m=False):
    obs = ~np.isnan(Xm)
    ok = np.isfinite(np.nanmean(Xm, 0))
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    Z = np.where(obs, (Xm - np.nan_to_num(mu)) / sd, 0.0)
    r0 = Vref.shape[1]
    rng = np.random.default_rng(seed)
    Vp = np.hstack([np.where(usable[:, None], Vref, 0.0), np.zeros((len(Vref), extra))])
    V = Vp + np.hstack([np.zeros_like(Vref), 0.01 * rng.normal(size=(len(Vref), extra))])
    r = r0 + extra
    U = np.zeros((Xm.shape[0], r))
    I = np.eye(r)
    # missingness block (ok columns only)
    if Mind is None:
        Mind = np.isnan(Xm).astype(float)
    Mo = np.asarray(Mind, float)[:, ok]
    Mc = Mo - Mo.mean(0)
    Wm = w * Mc
    Bt = np.zeros((int(ok.sum()), r))
    use_m = (w != 0) or force_m
    for _ in range(iters):
        if use_m:
            BtB = Bt.T @ Bt
            BW = Wm @ Bt  # row i = Bt^T Wm[i]
        for i in range(Xm.shape[0]):
            o = obs[i] & ok
            if use_m:
                U[i] = np.linalg.solve(V[o].T @ V[o] + BtB + lam * I, V[o].T @ Z[i, o] + BW[i])
            else:
                U[i] = np.linalg.solve(V[o].T @ V[o] + lam * I, V[o].T @ Z[i, o])
        for j in np.where(ok)[0]:
            o = obs[:, j]
            Uo = U[o]
            V[j] = np.linalg.solve(Uo.T @ Uo + lam_v * I, Uo.T @ Z[o, j] + lam_v * Vp[j])
        if use_m:
            Bt = np.linalg.solve(U.T @ U + lam_b * I, U.T @ Wm).T
    out = base.copy()
    P = U @ V.T
    miss = ~obs & ok[None, :]
    out[miss] = (mu[None, :] + sd[None, :] * P)[miss]
    if return_U:
        return out, U
    return out


# ---------------------------------------------------------------------------------------------------------------
# 2. indicators (true or within-cell permuted)
def make_mind(M, null, seed):
    Mind = np.isnan(M).astype(float)
    if not null:
        return Mind
    ok = np.isfinite(np.nanmean(M, 0))
    okc = np.where(ok)[0]
    rng = np.random.default_rng(seed + 4242)
    P = Mind.copy()
    for i in range(M.shape[0]):
        P[i, okc] = Mind[i, okc][rng.permutation(len(okc))]
    return P


# ---------------------------------------------------------------------------------------------------------------
# 3. psi selection (exactly tgtshrink._run, mode "null"), memoised per input within a process
_PSI_MEMO: dict = {}


def _key(M, ctx):
    h = hashlib.sha1(np.ascontiguousarray(M).tobytes()).hexdigest()
    return (ctx.get("target"), ctx["seed"], h)


def _select_psi(M, ctx):
    k = _key(M, ctx)
    if k in _PSI_MEMO:
        return _PSI_MEMO[k]
    V, warm, base, seed = ctx["V"], ctx["warm"], ctx["base"], ctx["seed"]
    Msel, (r1, c1), (r2, c2) = selection_mask(M, seed)
    scores = {}
    for p in PSI_GRID:
        mu, sd, _ = shrink(Msel, p, "null", seed, None)
        Y = adapt_project_mb(Msel, V, warm, base, mu, sd, lam_v=LAM_V, extra=EXTRA, iters=5)
        e1 = np.sqrt(np.mean((Y[r1, c1] - M[r1, c1]) ** 2))
        e2 = np.sqrt(np.mean((Y[r2, c2] - M[r2, c2]) ** 2))
        scores[p] = 0.5 * (e1 + e2)
    psi = min(PSI_GRID, key=lambda p: (scores[p], p))
    res = {"psi": psi, "scores": scores, "Msel": Msel, "h1": (r1, c1), "h2": (r2, c2)}
    _PSI_MEMO.clear()  # keep at most one entry (cand and null are called back to back on the same input)
    _PSI_MEMO[k] = res
    return res


def _sig(M):
    return hashlib.sha1(np.packbits(np.isnan(M)).tobytes()).hexdigest()[:16]


def _log(variant, rec):
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / f"log_{variant}.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def _run(M, ctx, null=False, w_grid=W_GRID, fixed_w=None, variant="main"):
    V, warm, base, seed = ctx["V"], ctx["warm"], ctx["base"], ctx["seed"]
    S = _select_psi(M, ctx)
    psi, Msel = S["psi"], S["Msel"]
    (r1, c1), (r2, c2) = S["h1"], S["h2"]
    wscores = {0.0: float(S["scores"][psi])}
    if fixed_w is not None:
        w = float(fixed_w)
    else:
        if len(w_grid):
            mu_s, sd_s, _ = shrink(Msel, psi, "null", seed, None)
            Mind_s = make_mind(Msel, null, seed)
            for wv in w_grid:
                Y = adapt_project_mj(Msel, V, warm, base, mu_s, sd_s, wv, Mind_s, lam_v=LAM_V, lam=LAM, lam_b=LAM_B,
                                     extra=EXTRA, iters=5)
                e1 = np.sqrt(np.mean((Y[r1, c1] - M[r1, c1]) ** 2))
                e2 = np.sqrt(np.mean((Y[r2, c2] - M[r2, c2]) ** 2))
                wscores[float(wv)] = float(0.5 * (e1 + e2))
        w = min(wscores, key=lambda x: (wscores[x], x))
    mu, sd, _ = shrink(M, psi, "null", seed, None)
    Mind = make_mind(M, null, seed)
    out, U = adapt_project_mj(M, V, warm, base, mu, sd, w, Mind, lam_v=LAM_V, lam=LAM, lam_b=LAM_B, extra=EXTRA,
                              iters=10, return_U=True)
    if w != 0:
        _, U0 = adapt_project_mj(M, V, warm, base, mu, sd, 0.0, None, lam_v=LAM_V, lam=LAM, extra=EXTRA, iters=10,
                                 return_U=True)
        du = float(np.mean(np.abs(U - U0)) / max(np.mean(np.abs(U0)), 1e-12))
    else:
        du = 0.0
    _log(variant, {"target": ctx.get("target"), "seed": seed, "n_missing": int(np.isnan(M).sum()),
                   "nan_sig": _sig(M),
                   "null": bool(null), "chosen_psi": psi, "chosen_w": w,
                   "psi_scores": {str(k): float(v) for k, v in S["scores"].items()},
                   "w_scores": {str(k): v for k, v in wscores.items()}, "rel_dU": du})
    return out


def impute(Xm, ctx):
    return _run(Xm, ctx, null=False, variant="main")


def impute_null(Xm, ctx):
    return _run(Xm, ctx, null=True, variant="null")


def impute_base(Xm, ctx):
    return _run(Xm, ctx, null=False, w_grid=(), variant="base")


def impute_w1(Xm, ctx):
    return _run(Xm, ctx, null=False, fixed_w=1.0, variant="w1")


def impute_w1_null(Xm, ctx):
    return _run(Xm, ctx, null=True, fixed_w=1.0, variant="w1_null")


# ---------------------------------------------------------------------------------------------------------------
# unit check
def _harness_ctx(t, kind, seed):
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H
    from scp_graphfs_mcp import phase2 as P
    from scp_graphfs_mcp.phasec import nfbase as NF

    Xs, Xm, mask, names = NF._split(t, kind, seed)
    V, warm = H.target_models(t, names)
    b = P._fill(Xm, P.loo_knn(Xm))
    b = np.where(np.isnan(b), np.nanmean(Xm), b)
    return Xm, {"names": names, "target": t, "seed": seed, "base": b, "V": V, "warm": warm}


def unit_check():
    from scp_graphfs_mcp.phasec import nfbase as NF
    from scp_graphfs_mcp.phasec import tgtshrink as TS

    TS._log = lambda *a, **k: None  # never write into tgtshrink's cache from here
    out = {}
    for t in NF.FAST:
        Xm, ctx = _harness_ctx(t, "MCAR", 100)
        ref = TS.impute_null(Xm.copy(), ctx)
        mine = impute_base(Xm.copy(), ctx)
        d1 = float(np.max(np.abs(mine - ref)))
        mu, sd, _ = shrink(Xm, 4.0, "null", 100, None)
        a = adapt_project_mb(Xm, ctx["V"], ctx["warm"], ctx["base"], mu, sd, lam_v=LAM_V, extra=EXTRA)
        bmj = adapt_project_mj(Xm, ctx["V"], ctx["warm"], ctx["base"], mu, sd, 0.0, None, lam_v=LAM_V, extra=EXTRA)
        d2 = float(np.max(np.abs(a - bmj)))
        # the general (M-block) branch with w = 0 forced through it: Wm = 0 so Bt stays 0
        bfm = adapt_project_mj(Xm, ctx["V"], ctx["warm"], ctx["base"], mu, sd, 0.0, None, lam_v=LAM_V, extra=EXTRA,
                               force_m=True)
        d2 = max(d2, float(np.max(np.abs(a - bfm))))
        out[t] = {"max_abs_base_vs_tgtshrink_null": d1, "pass_base": bool(d1 < 1e-10),
                  "max_abs_mj_w0_vs_mb": d2, "pass_mj": bool(d2 < 1e-12)}
        print(t, out[t], flush=True)
    out["all_pass"] = bool(all(v["pass_base"] and v["pass_mj"] for k, v in out.items() if isinstance(v, dict)))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "mnarjoint_unit.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


# ---------------------------------------------------------------------------------------------------------------
# reporting: per-mask means of per_target_mask rows + chosen-w distribution per mask from the logs
def _n_missing_map():
    """NaN-pattern signature -> 'MCAR'/'MNAR' x 'eval'/'inner', rebuilt exactly as the harness builds its inputs
    (n_missing alone is ambiguous: MCAR and MNAR hide the same number of entries)."""
    from scp_graphfs_mcp.phasec import nfbase as NF

    f = CACHE / "nansig_map.json"
    if f.exists():
        return json.loads(f.read_text())
    mp = {}
    for t in NF.FAST:
        for kind in ("MCAR", "MNAR"):
            for seed in (100, 101):
                Xs, Xm, mask, names = NF._split(t, kind, seed)
                eligible = ((~np.isnan(Xs)).mean(0) >= 0.5) & ~np.isnan(Xs).all(0)
                rng = np.random.default_rng(seed + 7919)
                c = np.argwhere(~np.isnan(Xm) & eligible[None, :])
                pick = c[rng.choice(len(c), max(1, int(0.1 * len(c))), replace=False)]
                Xi = Xm.copy()
                Xi[pick[:, 0], pick[:, 1]] = np.nan
                mp[_sig(Xm)] = f"{kind}/eval"
                mp[_sig(Xi)] = f"{kind}/inner"
    CACHE.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(mp))
    return mp


def report(ids):
    mp = _n_missing_map()
    res = {}
    for iid in ids:
        p = OUT / f"{iid}_fast.json"
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        rows = d["summary"]["rmse_all"]["per_target_mask"]
        pm = {}
        for kind in ("MCAR", "MNAR"):
            rr = [r for r in rows if r["mask"] == kind]
            pm[kind] = {m: round(float(np.mean([r[m] for r in rr])), 4)
                        for m in ("cand_raw", "cand_solo", "cand_stack", "null_stack") if m in rr[0]}
            ws = [r["weights"][3] for r in d["rows"] if r["method"] == "cand_stack" and r["mask"] == kind]
            pm[kind]["cand_weight"] = round(float(np.mean(ws)), 4)
        res[iid] = {"summary": d["summary"]["rmse_all"]["mean_rel_to_champion"],
                    "in_complex": d["summary"]["rmse_in_complex"]["mean_rel_to_champion"],
                    "splits_le_0.99": d["summary"]["rmse_all"].get("n_splits_cand_stack_le_0.99"),
                    "mean_cand_weight": d["summary"]["mean_cand_weight_in_stack"],
                    "per_mask": pm, "runtime_s": d["meta"]["runtime_s"]}
    logs = {}
    for v in ("main", "null", "base", "w1", "w1_null"):
        f = CACHE / f"log_{v}.jsonl"
        if not f.exists():
            continue
        cnt, du, seen = {}, {}, set()
        for line in f.read_text().splitlines():
            rec = json.loads(line)
            sig = rec.get("nan_sig")
            if sig in seen:  # deterministic duplicates (e.g. impute_base as cand and null, impute_null in two runs)
                continue
            seen.add(sig)
            if "nan_sig" not in rec:
                continue
            key = mp.get(rec["nan_sig"], "unknown")
            cnt.setdefault(key, {}).setdefault(str(rec["chosen_w"]), 0)
            cnt[key][str(rec["chosen_w"])] += 1
            du.setdefault(key, []).append(rec["rel_dU"])
        logs[v] = {"chosen_w_counts": cnt, "mean_rel_dU": {k: round(float(np.mean(x)), 4) for k, x in du.items()}}
    res["logs"] = logs
    (OUT / "mnarjoint_report.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    if "unit" in sys.argv:
        unit_check()
    elif "report" in sys.argv:
        report(["mnarjoint", "mnarjoint_base", "mnarjoint_nullraw", "mnarjoint_w1"])
