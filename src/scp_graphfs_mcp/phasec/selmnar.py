"""selmnar (Phase C, queen: representation & mathematics): a Heckman-style probit selection model inside an EM U-step
of the domain-adapted projection. Encodes the missingness pattern (MNAR), not the values.

STEP 0  local copy of transfer.adapt_project (lam_v=100, extra=10, lam=5, iters=10, seed=0) -> U, V, mu, sd, ok, obs;
        asserted equal to ctx['champion']. Predictive N(m_ij, v_ij), m = mu + sd*(U V^T), v_j = sd_j^2 * s2_j.
STEP 1  probit selection P(obs | x) = Phi((x - c_ij)/sigma_d), c_ij = t_i + delta_j (delta ~ N(0, tau^2)),
        fitted by L-BFGS-B on the Heckman likelihood (observed: -log Phi((x-c)/sd_d); missing:
        -log Phi((c-m)/sqrt(sd_d^2+v))).
STEP 2  E-step: xhat = m - rho*(v/s)*phi(alpha)/Phi(alpha), s = sqrt(sd_d^2 + v), alpha = (c - m)/s.
STEP 3  ridge U-step with pseudo-observations on missing ok entries, V fixed at the STEP-0 adapted V.
4 EM iterations, final E-step. Heckman JJ (1979) Econometrica 47(1):153-161.

    python -m scp_graphfs_mcp.phasec.selmnar diag     # -> results/phaseC/selmnar_diag.json
See docs/phaseC/ideas/selmnar.md.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import log_ndtr, ndtr

ROOT = Path(__file__).resolve().parents[3]
LAM_V, EXTRA, LAM, ITERS = 100.0, 10, 5.0, 10
EM_ITERS = 4
RHOS = (0.0, 0.25, 0.5, 1.0)
LOG_SQRT_2PI = 0.5 * np.log(2 * np.pi)


def _mills(a):
    """phi(a)/Phi(a), stable."""
    return np.exp(-0.5 * a * a - LOG_SQRT_2PI - log_ndtr(a))


# ---------------------------------------------------------------- STEP 0
def adapt_state(Xm, Vref, usable, base, lam_v=LAM_V, lam=LAM, extra=EXTRA, iters=ITERS, seed=0):
    """Verbatim copy of transfer.adapt_project that also returns its state."""
    obs = ~np.isnan(Xm)
    mu = np.nanmean(Xm, 0)
    sd = np.nanstd(Xm, 0)
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    ok = np.isfinite(mu)
    Z = np.where(obs, (Xm - np.nan_to_num(mu)) / sd, 0.0)
    r0 = Vref.shape[1]
    rng = np.random.default_rng(seed)
    Vp = np.hstack([np.where(usable[:, None], Vref, 0.0), np.zeros((len(Vref), extra))])
    V = Vp + np.hstack([np.zeros_like(Vref), 0.01 * rng.normal(size=(len(Vref), extra))])
    r = r0 + extra
    U = np.zeros((Xm.shape[0], r))
    I = np.eye(r)
    for _ in range(iters):
        for i in range(Xm.shape[0]):
            o = obs[i] & ok
            U[i] = np.linalg.solve(V[o].T @ V[o] + lam * I, V[o].T @ Z[i, o])
        for j in np.where(ok)[0]:
            o = obs[:, j]
            Uo = U[o]
            V[j] = np.linalg.solve(Uo.T @ Uo + lam_v * I, Uo.T @ Z[o, j] + lam_v * Vp[j])
    out = base.copy()
    P = U @ V.T
    miss = ~obs & ok[None, :]
    out[miss] = (mu[None, :] + sd[None, :] * P)[miss]
    return {"out": out, "U": U, "V": V, "mu": mu, "sd": sd, "ok": ok, "obs": obs, "Z": Z, "lam": lam}


def _s2(st):
    """Per-protein predictive z-variance (in-sample residual, shrunk, x1.5, floor 0.05)."""
    obs, ok, Z = st["obs"], st["ok"], st["Z"]
    P = st["U"] @ st["V"].T
    R = np.where(obs, (Z - P) ** 2, 0.0)
    n = obs.sum(0).astype(float)
    r = np.where(n > 0, R.sum(0) / np.maximum(n, 1), np.nan)
    med = float(np.nanmedian(r[ok]))
    s2 = 1.5 * (n * np.nan_to_num(r) + 5 * med) / (n + 5)
    return np.maximum(s2, 0.05)


# ---------------------------------------------------------------- STEP 1
def fit_selection(X, M, v, obsI, theta0, tau, sd_med, fixed_sigma=None, marginal=False):
    """X, M: cells x p_ok raw values (X NaN where missing) and predictive means; v: p_ok variances;
    obsI: indicator used by the selection likelihood. theta = [t (n), delta (p), log sigma_d]."""
    n, p = M.shape
    N = float(n * p)
    lo, hi = np.log(0.05 * sd_med), np.log(100 * sd_med)
    Xf = np.where(obsI, np.nan_to_num(X), 0.0)

    def f(th):
        t, d, ls = th[:n], th[n:n + p], th[-1]
        sg = np.exp(ls)
        c = t[:, None] + d[None, :]
        s2 = sg * sg + v[None, :]
        s = np.sqrt(s2)
        # missing term (also the marginal observed term in the null, with the sign flipped)
        b = (c - M) / s
        bm = np.where(obsI, -b, b) if marginal else b
        if marginal:
            lm = log_ndtr(bm)
            Mb = _mills(bm)
            # d(-log Phi(bm))/dc = -Mb * dbm/dc ; dbm/dc = +-1/s
            sign = np.where(obsI, -1.0, 1.0)
            gc = -Mb * sign / s
            # dbm/dsigma = -bm*sigma/s2
            gls = Mb * bm * sg * sg / s2
            val = -lm.sum()
        else:
            a = (Xf - c) / sg
            la = log_ndtr(a)
            Ma = _mills(a)
            lb = log_ndtr(b)
            Mb = _mills(b)
            val = -(np.where(obsI, la, lb)).sum()
            gc = np.where(obsI, Ma / sg, -Mb / s)
            gls = np.where(obsI, Ma * a, Mb * b * sg * sg / s2)
        val += (d * d).sum() / (2 * tau * tau)
        g = np.concatenate([gc.sum(1), gc.sum(0) + d / (tau * tau), [gls.sum()]])
        return val / N, g / N

    th0 = theta0.copy()
    bounds = [(None, None)] * (n + p) + [(lo, hi)]
    if fixed_sigma is not None:
        th0[-1] = np.log(fixed_sigma)
        return th0
    res = minimize(f, th0, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": 300})
    return res.x


# ---------------------------------------------------------------- EM
def _em(st, s2, Xin, obs, theta, rho=1.0, fixed_sigma=None, sel_obs=None, U0=None):
    """Run STEPS 1-3 for EM_ITERS iterations on input Xin (cells x proteins) with V/mu/sd fixed.
    sel_obs: indicator used by the selection fit (None -> true obs, standard Heckman likelihood)."""
    ok, mu, sd, V, lam = st["ok"], st["mu"], st["sd"], st["V"], st["lam"]
    j = np.where(ok)[0]
    Vk = V[j]
    muk, sdk = mu[j], sd[j]
    vk = sd[j] ** 2 * s2[j]
    Xk = Xin[:, j]
    ob = obs[:, j]
    Zo = np.where(ob, (np.nan_to_num(Xk) - muk) / sdk, 0.0)
    U = st["U"].copy() if U0 is None else U0.copy()
    n, p = Xk.shape
    tau = float(np.median(sdk))
    sdm = tau
    if theta is None:
        t0 = np.array([np.nanpercentile(Xk[i, ob[i]], 5) if ob[i].any() else np.nanmin(muk) for i in range(n)])
        theta = np.concatenate([t0, np.zeros(p), [np.log(sdm)]])
    Ainv = np.linalg.inv(Vk.T @ Vk + lam * np.eye(Vk.shape[1]))  # all ok proteins carry an (pseudo)observation
    marginal = sel_obs is not None
    selI = ob if sel_obs is None else sel_obs[:, j]

    def estep(M, th):
        sg = np.exp(th[-1])
        c = th[:n, None] + th[None, n:n + p]
        s = np.sqrt(sg * sg + vk[None, :])
        alpha = (c - M) / s
        return M - rho * (vk[None, :] / s) * _mills(alpha)

    for _ in range(EM_ITERS):
        M = muk + sdk * (U @ Vk.T)
        theta = fit_selection(Xk, M, vk, selI, theta, tau, sdm, fixed_sigma=fixed_sigma, marginal=marginal)
        xh = estep(M, theta)
        Zf = np.where(ob, Zo, (xh - muk) / sdk)
        U = Zf @ Vk @ Ainv
    M = muk + sdk * (U @ Vk.T)
    xh = estep(M, theta)
    return {"xhat": xh, "M": M, "theta": theta, "U": U, "j": j, "ob": ob}


def _state(Xm, ctx):
    st = adapt_state(Xm, ctx["V"], ctx["warm"], ctx["base"])
    assert np.allclose(st["out"], ctx["champion"], atol=1e-8, equal_nan=True), "STEP-0 copy != champion"
    return st


def _assemble(ctx, st, em):
    out = np.array(ctx["champion"], float, copy=True)
    j = em["j"]
    miss = ~em["ob"]
    blk = out[:, j]
    blk[miss] = em["xhat"][miss]
    out[:, j] = blk
    return out


def _diag(st, em):
    th = em["theta"]
    sdm = float(np.median(st["sd"][em["j"]]))
    miss = ~em["ob"]
    return {"sigma_d_over_med_sd": float(np.exp(th[-1]) / sdm), "mean_tilt": float((em["M"] - em["xhat"])[miss].mean()),
            "obs_frac": float(em["ob"].mean())}


LAST = {}


def _run(Xm, ctx, rho=1.0, fixed_sigma=None, null=False):
    st = _state(Xm, ctx)
    s2 = _s2(st)
    sel = None
    if null:
        rng = np.random.default_rng(ctx["seed"] + 4242)
        sel = st["obs"].copy()
        j = np.where(st["ok"])[0]
        for i in range(sel.shape[0]):
            sel[i, j] = sel[i, j][rng.permutation(len(j))]
    em = _em(st, s2, Xm, st["obs"], None, rho=rho, fixed_sigma=fixed_sigma, sel_obs=sel)
    LAST.clear()
    LAST.update(_diag(st, em))
    return _assemble(ctx, st, em), st, s2, em


def impute(Xm, ctx):
    return _run(Xm, ctx)[0]


def impute_null(Xm, ctx):
    return _run(Xm, ctx, null=True)[0]


def impute_inf(Xm, ctx):
    return _run(Xm, ctx, fixed_sigma=1e6)[0]


def impute_cal(Xm, ctx):
    _, st, s2, em = _run(Xm, ctx, rho=1.0)
    j, ob = em["j"], em["ob"]
    th = em["theta"]
    n, p = ob.shape
    sg = np.exp(th[-1])
    c = th[:n, None] + th[None, n:n + p]
    Xk = Xm[:, j]
    elig = ob.mean(0) >= 0.5
    cand = np.argwhere(ob & elig[None, :])
    w = ndtr((c[cand[:, 0], cand[:, 1]] - Xk[cand[:, 0], cand[:, 1]]) / sg)
    w = np.maximum(w, 1e-300)
    rng = np.random.default_rng(ctx["seed"] + 1)
    k = max(1, int(0.1 * len(cand)))
    pick = cand[rng.choice(len(cand), k, replace=False, p=w / w.sum())]
    Xi = Xm.copy()
    cols = j[pick[:, 1]]
    truth = Xm[pick[:, 0], cols]
    Xi[pick[:, 0], cols] = np.nan
    obs_i = ~np.isnan(Xi)
    # one ridge U-step on the observed entries of the internal input, V fixed
    V, mu, sd, lam = st["V"], st["mu"], st["sd"], st["lam"]
    Zi = np.where(obs_i, (Xi - np.nan_to_num(mu)) / sd, 0.0)
    U0 = np.zeros_like(st["U"])
    I = np.eye(V.shape[1])
    for i in range(Xi.shape[0]):
        o = obs_i[i] & st["ok"]
        U0[i] = np.linalg.solve(V[o].T @ V[o] + lam * I, V[o].T @ Zi[i, o])
    errs = []
    for rho in RHOS:
        e = _em(st, s2, Xi, obs_i, th.copy(), rho=rho, U0=U0)
        errs.append(float(np.mean((e["xhat"][pick[:, 0], pick[:, 1]] - truth) ** 2)))
    rho = RHOS[int(np.argmin(errs))]
    LAST_CAL.append({"target": ctx.get("target"), "seed": ctx.get("seed"), "errs": errs, "rho": rho})
    if rho == 1.0:
        return _assemble(ctx, st, em)
    em2 = _em(st, s2, Xm, st["obs"], None, rho=rho)
    return _assemble(ctx, st, em2)


LAST_CAL: list = []


# ---------------------------------------------------------------- diag
def diag():
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H
    from explore_phaseb import load_target, subsample
    from run import make_mask

    from scp_graphfs_mcp import transfer as T

    out = {"targets": {}}
    for tname in H.FAST:
        X, names, labels = load_target(tname)
        V, warm = H.target_models(tname, names)
        seed = 100
        idx = subsample(X.shape[0], labels, 30, seed)
        Xs = X[idx]
        eligible = ((~np.isnan(Xs)).mean(0) >= 0.5) & ~np.isnan(Xs).all(0)
        rec = {}
        for kind in ("MCAR", "MNAR"):
            mask = make_mask(Xs, eligible, kind, seed)
            Xm = Xs.copy()
            Xm[mask] = np.nan
            from scp_graphfs_mcp import phase2 as P
            b = P._fill(Xm, P.loo_knn(Xm))
            b = np.where(np.isnan(b), np.nanmean(Xm), b)
            ad = T.adapt_project(Xm, V, warm, b, lam_v=LAM_V, extra=EXTRA)
            ctx = {"names": names, "target": tname, "seed": seed, "base": b, "champion": ad, "V": V, "warm": warm}
            yc = impute(Xm.copy(), ctx)
            dc = dict(LAST)
            yn = impute_null(Xm.copy(), ctx)
            dn = dict(LAST)
            yi = impute_inf(Xm.copy(), ctx)
            y2 = impute_v2(Xm.copy(), ctx)
            d2 = dict(LAST)
            yn2 = impute_null_v2(Xm.copy(), ctx)
            dn2 = dict(LAST)
            sd = np.nanstd(Xm, 0)
            sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
            dev = float(np.nanmax(np.abs(yi - ad) / sd[None, :]))
            tr = Xs[mask]

            def rm(Y):
                return float(np.sqrt(np.mean((Y[mask] - tr) ** 2)))

            rec[kind] = {"cand": dc, "null": dn, "v2": d2, "null_v2": dn2, "inf_max_abs_dev_over_sd": dev,
                         "mask_bias_true_minus_champion": float(np.mean(tr - ad[mask])),
                         "mask_bias_true_minus_cand": float(np.mean(tr - yc[mask])),
                         "rmse_rel_adapted_alone": {"cand": rm(yc) / rm(ad), "null": rm(yn) / rm(ad),
                                                    "inf": rm(yi) / rm(ad), "v2": rm(y2) / rm(ad),
                                                    "null_v2": rm(yn2) / rm(ad)}}
            print(tname, kind, json.dumps(rec[kind]), flush=True)
        out["targets"][tname] = rec
    f = ROOT / "results" / "phaseC" / "selmnar_inf_fast.json"
    if f.exists():
        out["selmnar_inf_cand_raw"] = json.loads(f.read_text())["summary"]["rmse_all"]["mean_rel_to_champion"]["cand_raw"]
    out["inf_unit_check_pass"] = all(r[k]["inf_max_abs_dev_over_sd"] < 1e-3 for r in out["targets"].values()
                                     for k in r)
    out["note"] = ("one seed (100) per fast target, eval input only (no inner mask); rmse_rel_adapted_alone is vs the "
                   "unstacked adapted projection, not the harness champion stack")
    (ROOT / "results" / "phaseC" / "selmnar_diag.json").write_text(json.dumps(out, indent=1) + "\n")


# ================================================================ v2 (the one revision): MAR + MNAR mixture selection
# P(miss | x) = pi + (1 - pi) * Phi((c - x)/sigma_d): a global random-dropout probability pi (MAR component) plus the
# probit detection limit. v1 forced the probit to explain every hole, including random dropouts (e.g. the MCAR eval
# mask), so every missing entry was tilted down (fast MCAR cand_raw 1.2721). Posterior mean of
# N(m, v) * [pi + (1-pi) Phi((c-x)/sigma_d)]:  xhat = m - (v/s) * (1-pi) phi(alpha) / (pi + (1-pi) Phi(alpha)).
_LOG_PHI = lambda a: -0.5 * a * a - LOG_SQRT_2PI  # noqa: E731


def _mix_ratio(b, lp, l1p):
    """(1-pi) phi(b) / (pi + (1-pi) Phi(b)) and log D, stable."""
    logD = np.logaddexp(lp, l1p + log_ndtr(b))
    return np.exp(l1p + _LOG_PHI(b) - logD), logD


def fit_selection_mix(X, M, v, obsI, theta0, tau, sd_med, marginal=False):
    """theta = [t (n), delta (p), log sigma_d, logit pi]."""
    n, p = M.shape
    N = float(n * p)
    Xf = np.where(obsI, np.nan_to_num(X), 0.0)

    def f(th):
        t, d, ls, eta = th[:n], th[n:n + p], th[n + p], th[n + p + 1]
        sg = np.exp(ls)
        pi = 1.0 / (1.0 + np.exp(-eta))
        lp, l1p = np.log(pi), np.log1p(-pi)
        c = t[:, None] + d[None, :]
        s2 = sg * sg + v[None, :]
        s = np.sqrt(s2)
        b = (c - M) / s
        # missing: -log(pi + (1-pi) Phi(b))
        Rb, logD = _mix_ratio(b, lp, l1p)
        Pb = ndtr(b)
        gc_m = -Rb / s
        gls_m = Rb * b * sg * sg / s2
        gpi_m = -(1 - Pb) / np.exp(logD)
        if marginal:  # observed (permuted): -log((1-pi) Phi(-b))
            a = -b
            la = log_ndtr(a)
            Ma = _mills(a)
            val_o = -(l1p + la)
            gc_o = Ma / s
            gls_o = Ma * a * sg * sg / s2
        else:
            a = (Xf - c) / sg
            la = log_ndtr(a)
            Ma = _mills(a)
            val_o = -(l1p + la)
            gc_o = Ma / sg
            gls_o = Ma * a
        gpi_o = 1.0 / (1 - pi)
        val = np.where(obsI, val_o, -logD).sum() + (d * d).sum() / (2 * tau * tau)
        gc = np.where(obsI, gc_o, gc_m)
        gls = np.where(obsI, gls_o, gls_m).sum()
        geta = np.where(obsI, gpi_o, gpi_m).sum() * pi * (1 - pi)
        g = np.concatenate([gc.sum(1), gc.sum(0) + d / (tau * tau), [gls, geta]])
        return val / N, g / N

    bounds = [(None, None)] * (n + p) + [(np.log(0.05 * sd_med), np.log(100 * sd_med)), (-9.0, 2.2)]
    res = minimize(f, theta0, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": 300})
    return res.x


def _em_mix(st, s2, Xin, obs, sel_obs=None):
    ok, mu, sd, V, lam = st["ok"], st["mu"], st["sd"], st["V"], st["lam"]
    j = np.where(ok)[0]
    Vk, muk, sdk = V[j], mu[j], sd[j]
    vk = sd[j] ** 2 * s2[j]
    Xk, ob = Xin[:, j], obs[:, j]
    Zo = np.where(ob, (np.nan_to_num(Xk) - muk) / sdk, 0.0)
    U = st["U"].copy()
    n, p = Xk.shape
    tau = float(np.median(sdk))
    t0 = np.array([np.nanpercentile(Xk[i, ob[i]], 5) if ob[i].any() else np.nanmin(muk) for i in range(n)])
    theta = np.concatenate([t0, np.zeros(p), [np.log(tau), np.log(0.1 / 0.9)]])
    Ainv = np.linalg.inv(Vk.T @ Vk + lam * np.eye(Vk.shape[1]))
    marginal = sel_obs is not None
    selI = ob if sel_obs is None else sel_obs[:, j]

    def estep(M, th):
        sg = np.exp(th[n + p])
        pi = 1.0 / (1.0 + np.exp(-th[n + p + 1]))
        c = th[:n, None] + th[None, n:n + p]
        s = np.sqrt(sg * sg + vk[None, :])
        R, _ = _mix_ratio((c - M) / s, np.log(pi), np.log1p(-pi))
        return M - (vk[None, :] / s) * R

    for _ in range(EM_ITERS):
        M = muk + sdk * (U @ Vk.T)
        theta = fit_selection_mix(Xk, M, vk, selI, theta, tau, tau, marginal=marginal)
        xh = estep(M, theta)
        U = np.where(ob, Zo, (xh - muk) / sdk) @ Vk @ Ainv
    M = muk + sdk * (U @ Vk.T)
    xh = estep(M, theta)
    return {"xhat": xh, "M": M, "theta": theta, "U": U, "j": j, "ob": ob}


def _run_v2(Xm, ctx, null=False):
    st = _state(Xm, ctx)
    s2 = _s2(st)
    sel = None
    if null:
        rng = np.random.default_rng(ctx["seed"] + 4242)
        sel = st["obs"].copy()
        j = np.where(st["ok"])[0]
        for i in range(sel.shape[0]):
            sel[i, j] = sel[i, j][rng.permutation(len(j))]
    em = _em_mix(st, s2, Xm, st["obs"], sel_obs=sel)
    d = _diag(st, em)
    n, p = em["ob"].shape
    d["sigma_d_over_med_sd"] = float(np.exp(em["theta"][n + p]) / np.median(st["sd"][em["j"]]))
    d["pi_mar"] = float(1.0 / (1.0 + np.exp(-em["theta"][n + p + 1])))
    LAST.clear()
    LAST.update(d)
    return _assemble(ctx, st, em)


def impute_v2(Xm, ctx):
    return _run_v2(Xm, ctx)


def impute_null_v2(Xm, ctx):
    return _run_v2(Xm, ctx, null=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "diag":
        diag()
