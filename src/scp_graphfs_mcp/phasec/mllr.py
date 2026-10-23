"""Phase-C idea `mllr` (queen: innovation abduction, round 5). Data-only, on the cshrink base.

MLLR-style shared adaptation of the reference loadings (Leggetter & Woodland 1995, Computer Speech & Language
9(2):171-185, doi:10.1006/csla.1995.0010): instead of only pulling each protein's loading row toward Vref_j
(lam_v = 100, the champion's per-protein MAP), the prior mean of every row becomes Vref_j A, with ONE shared
r0 x r0 matrix A fitted on all observed entries of the target, ridge-pulled toward I with strength kappa * n_usable.

  adapt_project_mllr: nfbase.adapt_project_mb line for line, plus per iteration
    (a) U-step (unchanged);
    (b) A-step (kappa finite): D = V1 - Vr A; R = Z - U2 Ve^T - U1 D^T on obs & ok; R_ij ~ vr_j^T A u1_i;
        row-major vec(A), x_ij = kron(vr_j, u1_i); H = sum_i kron(G_i, u1_i u1_i^T), G_i = Vr[o_i]^T Vr[o_i];
        b = vec(Vr^T R^T U1); (H + kap I) vec(A) = b + kap vec(I). diag=True: diagonal of A only.
        Null: Vr[perm] (fixed permutation of the usable proteins) replaces Vr in G_i and b only.
    (c) Vp[:, :r0] = Vr A; V-step (unchanged) pulls toward the new prior mean.
  kappa = inf skips (b) and (c): the code path is exactly adapt_project_mb (A = I).

  run: psi = cshrink's own selection (C._score, tag 'c', same cache) -> mu, sd = cshrink_stats(M, psi, psi);
       kappa selected on TS.selection_mask(M, seed) (50/50 MCAR/MNAR, rng seed + 31337), iters = 5,
       score = 0.5 (RMSE_MCAR + RMSE_MNAR), kappa = inf score = C._score (asserted), ties -> larger kappa;
       refit on M with iters = 10.

Note (identifiability, see docs/phaseC/ideas/mllr.md): with free per-cell factors U the likelihood is invariant to
A -> A Q, u -> Q^T u, D -> D Q for orthogonal Q; the ridge on U is isotropic, so only A A^T (the target's factor
covariance in reference coordinates, relative to the U ridge and the D ridge) is identified by data; the rotation part
of A is held only by the kappa ridge toward I.

Entry points: impute, impute_null, impute_base (kappa = inf), impute_diag.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from scp_graphfs_mcp.phasec import cshrink as C
from scp_graphfs_mcp.phasec import tgtshrink as TS
from scp_graphfs_mcp.phasec.nfbase import EXTRA, LAM_V

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "mllr"
INF = float("inf")
KGRID = (0.1, 1.0, 10.0, 100.0, INF)


# ---------------------------------------------------------------------------------------------------------------
# A-step
def a_step(Z, obs, ok, usable, Vr_fit, Vr_true, A, U, V, r0, kappa, diag=False):
    """Closed-form ridge update of A. Vr_fit is used in G_i and b (= Vr_true, or Vr[perm] in the null)."""
    n_us = int(usable.sum())
    kap = kappa * max(n_us, 1)
    U1, U2 = U[:, :r0], U[:, r0:]
    V1, Ve = V[:, :r0], V[:, r0:]
    D = V1 - Vr_true @ A
    m = obs & ok[None, :]
    R = np.where(m, Z - U2 @ Ve.T - U1 @ D.T, 0.0)
    q = r0 * r0
    H = np.zeros((q, q))
    for i in range(Z.shape[0]):
        o = m[i] & usable
        Vo = Vr_fit[o]
        G = Vo.T @ Vo
        H += np.kron(G, np.outer(U1[i], U1[i]))
    Rm = np.where(m & usable[None, :], R, 0.0)
    b = (Vr_fit.T @ Rm.T @ U1).reshape(-1)  # row-major vec
    trH = float(np.trace(H) / q)
    I0 = np.eye(r0)
    if diag:
        idx = np.arange(r0) * r0 + np.arange(r0)
        Hd = H[np.ix_(idx, idx)]
        a = np.linalg.solve(Hd + kap * np.eye(r0), b[idx] + kap * np.ones(r0))
        return np.diag(a), trH
    a = np.linalg.solve(H + kap * np.eye(q), b + kap * I0.reshape(-1))
    return a.reshape(r0, r0), trH


def adapt_project_mllr(Xm, Vref, usable, base, mu, sd, kappa, lam_v=LAM_V, lam=5.0, extra=EXTRA, iters=10, seed=0,
                       perm=None, diag=False, return_info=False):
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
    # --- mllr additions
    A = np.eye(r0)
    Vr = np.where(usable[:, None], Vref, 0.0)
    if perm is not None:
        uidx = np.where(usable)[0]
        Vr_fit = Vr.copy()
        Vr_fit[uidx] = Vr[uidx[perm]]
    else:
        Vr_fit = Vr
    use_A = np.isfinite(kappa)
    trH = None
    for _ in range(iters):
        for i in range(Xm.shape[0]):
            o = obs[i] & ok
            U[i] = np.linalg.solve(V[o].T @ V[o] + lam * I, V[o].T @ Z[i, o])
        if use_A:
            A, trH = a_step(Z, obs, ok, usable, Vr_fit, Vr, A, U, V, r0, kappa, diag=diag)
            Vp[:, :r0] = Vr @ A
        for j in np.where(ok)[0]:
            o = obs[:, j]
            Uo = U[o]
            V[j] = np.linalg.solve(Uo.T @ Uo + lam_v * I, Uo.T @ Z[o, j] + lam_v * Vp[j])
    out = base.copy()
    P = U @ V.T
    miss = ~obs & ok[None, :]
    out[miss] = (mu[None, :] + sd[None, :] * P)[miss]
    if return_info:
        return out, {"A": A, "trH_per_entry": trH}
    return out


# ---------------------------------------------------------------------------------------------------------------
# diagnostics of A
def a_diag(A):
    r0 = A.shape[0]
    I0 = np.eye(r0)
    s = np.linalg.svd(A, compute_uv=False)
    ev = np.linalg.eigvals(A)
    W, sv, Vt = np.linalg.svd(A)
    Q = W @ Vt  # polar A = Q S
    S = Vt.T @ np.diag(sv) @ Vt
    return {"normAI": float(np.linalg.norm(A - I0)), "sv_top": float(s[0]), "sv_bot": float(s[-1]),
            "sv": [float(x) for x in s], "eig_real": sorted([float(x) for x in ev.real], reverse=True),
            "normQI": float(np.linalg.norm(Q - I0)), "normSI": float(np.linalg.norm(S - I0))}


def _ks(k):
    return "inf" if np.isinf(k) else str(k)


def _log(variant, rec):
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / f"log_{variant}.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def _cshrink_psi(M, ctx):
    """cshrink's own psi selection, exactly as C.run (joint psi on C.GRID, ties to the smaller psi)."""
    Msel, p1, p2 = TS.selection_mask(M, ctx["seed"])
    scores = {p: C._score(M, Msel, (p1, p2), ctx, p, p, C.cshrink_stats, "c") for p in C.GRID}
    psi = min(C.GRID, key=lambda p: (scores[p], p))
    return psi, scores, Msel, (p1, p2)


def run(M, ctx, variant, kgrid=KGRID, null=False, diag=False, fixed_kappa=None, log=True):
    seed = ctx["seed"]
    V, warm, base = ctx["V"], ctx["warm"], ctx["base"]
    psi, psi_scores, Msel, ((r1, c1), (r2, c2)) = _cshrink_psi(M, ctx)
    mu, sd, _ = C.cshrink_stats(M, psi, psi)
    perm = None
    if null:
        perm = np.random.default_rng(seed + 4242).permutation(int(np.asarray(warm, bool).sum()))
    kscores, kdiag = {}, {}
    if fixed_kappa is not None:
        kappa = fixed_kappa
    else:
        mus, sds, _ = C.cshrink_stats(Msel, psi, psi)
        for k in kgrid:
            if np.isinf(k):
                s = psi_scores[psi]
                Y = C.adapt_project_mb(Msel, V, warm, base, mus, sds, lam_v=LAM_V, extra=EXTRA, iters=5) \
                    if log == "verify" else None
                if Y is not None:
                    s2 = 0.5 * (np.sqrt(np.mean((Y[r1, c1] - M[r1, c1]) ** 2))
                                + np.sqrt(np.mean((Y[r2, c2] - M[r2, c2]) ** 2)))
                    assert abs(s2 - s) < 1e-9, (s2, s)
                kscores[k] = float(s)
                continue
            Y, info = adapt_project_mllr(Msel, V, warm, base, mus, sds, k, iters=5, perm=perm, diag=diag,
                                         return_info=True)
            e1 = np.sqrt(np.mean((Y[r1, c1] - M[r1, c1]) ** 2))
            e2 = np.sqrt(np.mean((Y[r2, c2] - M[r2, c2]) ** 2))
            kscores[k] = float(0.5 * (e1 + e2))
            kdiag[_ks(k)] = {**a_diag(info["A"]), "trH_per_entry": info["trH_per_entry"]}
        kappa = min(kgrid, key=lambda k: (kscores[k], -k))  # ties -> larger kappa (toward the base)
    Y, info = adapt_project_mllr(M, V, warm, base, mu, sd, kappa, iters=10, perm=perm, diag=diag, return_info=True)
    if log:
        rec = {"target": ctx.get("target"), "seed": seed, "n_missing": int(np.isnan(M).sum()), "psi": _ks(psi),
               "kappa_scores": {_ks(k): v for k, v in kscores.items()}, "chosen_kappa": _ks(kappa),
               "final": {**a_diag(info["A"]), "trH_per_entry": info["trH_per_entry"]}, "sel_fits": kdiag}
        _log(variant, rec)
    return Y


def impute(Xm, ctx):
    return run(Xm, ctx, "main")


def impute_null(Xm, ctx):
    return run(Xm, ctx, "null", null=True)


def impute_base(Xm, ctx):
    return run(Xm, ctx, "base", fixed_kappa=INF)


def impute_diag(Xm, ctx):
    return run(Xm, ctx, "diag", diag=True)


def impute_cshrink(Xm, ctx):
    """cshrink.impute exactly (C.run 'main'), without appending to cshrink's own log files."""
    return C.run(Xm, ctx, "main", log=False)


# ---------------------------------------------------------------------------------------------------------------
# unit checks
def unit_check():
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H  # noqa: F401
    from scp_graphfs_mcp import phase2 as P
    from scp_graphfs_mcp.phasec import nfbase as NF

    out = {"U1_base_vs_cshrink": {}, "U2_kappa1e12_vs_inf": {}, "U3_synthetic": {}}
    for t in ("petrosius2023", "leduc2022"):
        Xs, Xm, mask, names = NF._split(t, "MCAR", 100)
        V, warm = H.target_models(t, names)
        b = P._fill(Xm, P.loo_knn(Xm))
        b = np.where(np.isnan(b), np.nanmean(Xm), b)
        ctx = {"V": V, "warm": warm, "base": b, "seed": 100, "target": t}
        mine = run(Xm, ctx, "unit", fixed_kappa=INF, log=False)
        ref = C.run(Xm, ctx, "unit", log=False)
        miss = np.isnan(Xm)
        d = float(np.max(np.abs(mine[miss] - ref[miss])))
        out["U1_base_vs_cshrink"][t] = {"max_abs_diff": d, "pass": bool(d < 1e-9)}
        # the kappa = inf score inside selection equals C._score (assert path)
        run(Xm, ctx, "unit", log="verify")
        out["U1_base_vs_cshrink"][t]["inf_score_assert"] = "pass"
        if t == "petrosius2023":
            psi, *_ = _cshrink_psi(Xm, ctx)
            mu, sd, _ = C.cshrink_stats(Xm, psi, psi)
            y1 = adapt_project_mllr(Xm, V, warm, b, mu, sd, 1e12, iters=10)
            y2 = adapt_project_mllr(Xm, V, warm, b, mu, sd, INF, iters=10)
            d2 = float(np.max(np.abs(y1[miss] - y2[miss])))
            out["U2_kappa1e12_vs_inf"][t] = {"max_abs_diff": d2, "pass": bool(d2 < 1e-4)}
        print(t, out, flush=True)
    # synthetic: A-step normal equations with known U recover a known A0
    rng = np.random.default_rng(0)
    n, p, r0 = 30, 1500, 20
    Vr = rng.normal(size=(p, r0)) / np.sqrt(r0)
    Qr, _ = np.linalg.qr(rng.normal(size=(r0, r0)))
    A0 = Qr @ np.diag(np.linspace(0.5, 1.5, r0))
    U1 = rng.normal(size=(n, r0))
    Z = U1 @ (Vr @ A0).T + 0.1 * rng.normal(size=(n, p))
    obs = rng.random((n, p)) > 0.3
    Z = np.where(obs, Z, 0.0)
    ok = np.ones(p, bool)
    usable = np.ones(p, bool)
    U = np.hstack([U1, np.zeros((n, 1))])
    Vfull = np.hstack([Vr @ np.eye(r0), np.zeros((p, 1))])  # V1 = Vr A (A = I) -> D = 0
    Ahat, trH = a_step(Z, obs, ok, usable, Vr, Vr, np.eye(r0), U, Vfull, r0, 1e-6)
    rel = float(np.linalg.norm(Ahat - A0) / np.linalg.norm(A0))
    out["U3_synthetic"]["A_step_known_U"] = {"rel_frob_err": rel, "pass": bool(rel < 0.05), "trH_per_entry": trH}
    # full alternating fit on synthetic data: only A A^T is identified (U is free); report, not a pass criterion
    Xsyn = np.where(obs, Z, np.nan)
    _, info = adapt_project_mllr(Xsyn, Vr, usable, np.zeros((n, p)), np.zeros(p), np.ones(p), 1e-3, extra=1,
                                 iters=30, return_info=True)
    Af = info["A"]
    out["U3_synthetic"]["full_fit_info_only"] = {
        "rel_err_A": float(np.linalg.norm(Af - A0) / np.linalg.norm(A0)),
        "rel_err_AAt": float(np.linalg.norm(Af @ Af.T - A0 @ A0.T) / np.linalg.norm(A0 @ A0.T)),
        "note": "A identifiable only up to right-orthogonal Q (u -> Q^T u); AA^T also scaled by the U/D ridges"}
    print(out, flush=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "unit_check.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    if "--unit" in sys.argv:
        unit_check()
