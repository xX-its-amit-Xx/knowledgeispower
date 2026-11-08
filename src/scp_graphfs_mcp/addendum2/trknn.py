"""trKNN / KNN-TN: truncation k-nearest-neighbour imputation (post hoc, exploratory addendum; written 2026-11-08 after
PREREGISTRATION_v9 was frozen).

Reference: Shah JS, Rai SN, DeFilippis AP, Hill BG, Bhatnagar A, Brock GN (2017) "Distribution based nearest neighbor
imputation for truncated high dimensional data with applications to pre-clinical and clinical metabolomics studies."
BMC Bioinformatics 18:114, doi:10.1186/s12859-017-1547-6 (title, authors, journal, volume, article number and date
checked against Crossref on 2026-11-08; the abstract, read via Europe PMC, says the method "assumes that the data follow
a truncated normal distribution with the truncation point at the detection limit (LOD)"). It is a METABOLOMICS paper;
the full methods text could not be fetched from this session (publisher redirects to a login, PMC/Europe PMC full text
unavailable), so the algorithm below follows the reference R code, not the paper prose.
Reference code: NAguideR (Wang et al. 2020, Nucleic Acids Res, doi:10.1093/nar/gkaa498), files
  inst/NAguideRapp/Trunc_KNN/Imput_funcs.r  (Shah et al.'s functions mklhood, NewtonRaphsonLike, EstimatesComputation,
                                              imputeKNN; sha256 ca9651f062d487a2bd57f465a114fab4494a6c8dd47cafbcf49387c60b090bd4)
  inst/NAguideRapp/app.R lines 1807-1814     (the NAguideR call)
read from raw.githubusercontent.com/wangshisheng/NAguideR/master on 2026-11-08. The NAguideR call is
    imputeKNN(proteins x samples, k = 10, distance = 'truncation', perc = 0)
i.e. rows = PROTEINS, columns = the ~30 cells (our Xm.T), on NAguideR's log2 table.

Algorithm (imputeKNN with distance = 'truncation'):
  1. Per protein, a location/scale (mu_i, s_i) from EstimatesComputation(perc):
       start from the sample mean and SD (R sd, n - 1);
       keep them if the protein's missing fraction >= perc            ("Case 1"),
       or if mu_i > 3 s_i + LOD, LOD = the minimum observed value of the WHOLE matrix ("Case 2");
       otherwise ("Case 3") replace them by the maximum-likelihood estimates of a normal truncated to
       [LOD, max(X) + 2 max_i s_i] fitted to the protein's observed values by Newton-Raphson (start = sample
       estimates, tol 1e-7, maxit 100), accepted only if it converged in fewer than 50 iterations.
  2. Standardise every protein: z = (x - mu_i) / s_i.
  3. For a missing entry (protein i, cell e): candidates = proteins observed in cell e; r = Pearson correlation of
     protein i with each candidate over their pairwise-complete cells (R cor(..., use = "pairwise.complete.obs"));
     distance d = 1 - |r|, with NA and |r| == 1 set to Inf.
  4. The k nearest candidates (R order(), ascending d) give z_i(e) = sum_j w_j sign(r_j) z_j(e),
     w_j = (1/d_j) / sum_l (1/d_l); back-transform x = mu_i + s_i z_i(e).
  Neighbour values are always observed ones (no reuse of imputed values), so the order of visits does not matter.

IMPORTANT, found while porting: with NAguideR's perc = 0, Case 1 (missing fraction >= 0) holds for EVERY protein, so the
truncated-normal step never runs and NAguideR's "trknn" is exactly KNN-CR (correlation kNN on sample-standardised
values). `trknn` below reproduces that call. `trknn_tn` uses perc = 1 (imputeKNN's own default), where the
truncated-normal estimates are actually used for proteins that are not fully missing and not > 3 SD above the LOD. Which
perc the paper's own experiments used could NOT be verified (full text not readable here).

Deviations from the R code (edge cases where R errors out or returns NA; none changes a well-posed entry):
  - R stops with "Fewer than K finite distances found" when an entry has < k candidates with finite distance; here the
    finite ones are used, and an entry with none gets ctx['base'] (via baselines._keep_obs).
  - Proteins with < 2 observed values or zero SD (R: NA/NaN standardisation, all distances Inf, then the stop above)
    are not imputed by this method (ctx['base']) and are not candidates.
  - A correlation over fewer than 3 co-observed cells is treated as NA (with 2 cells |r| is 1 mathematically, which R
    sets to Inf; computing it in floating point can give |r| = 1 - 1e-16 and a near-infinite weight).
  - Truncated-normal step: the integrals R evaluates with integrate() are computed in closed form (normal pdf/cdf
    moments); the upper limit uses max over proteins of the FINITE sample SDs (R's max() without na.rm returns NA
    whenever any protein has a single observation, which makes every Newton-Raphson fit fail and silently reverts
    to Case 1). The R Hessian has a typo in d2psi/dsigma2 (5*(y-mu) where 5*(y-mu)^2 is correct); it is reproduced by
    default (`hessian_typo=True`) because it changes Newton-Raphson's path and hence which fits pass the < 50 iteration
    rule; the gradient (which fixes the converged point) is correct in R and here. R's convergence test
    |ll_old - ll_new| / ll_old < tol divides by the signed negative log-likelihood, so it passes at once whenever that
    is negative; reproduced as written.
  - Exact distance ties at the k-th neighbour may be broken differently from R's stable order().
  - NaN payloads are canonicalised first (see addendum2/simple.py: leduc2022 carries R's NA_real_ bit pattern).

Same signature as baselines.py: fn(Xm, ctx) -> full matrix; only unobserved entries are scored.
"""

from __future__ import annotations

import numpy as np
from scipy.special import ndtr

from ..baselines_addendum import _keep_obs
from .simple import _canon

_SQ2PI = np.sqrt(2.0 * np.pi)


# --------------------------------------------------------------------------- truncated-normal MLE (EstimatesComputation)


def _moments(mu, s, a, b):
    """M_m = int_alpha^beta z^m phi(z) dz, m = 0..4, alpha = (a-mu)/s, beta = (b-mu)/s (vectorised)."""
    al, be = (a - mu) / s, (b - mu) / s
    pa, pb = np.exp(-0.5 * al * al) / _SQ2PI, np.exp(-0.5 * be * be) / _SQ2PI
    M0 = ndtr(be) - ndtr(al)
    M1 = pa - pb
    M2 = M0 + al * pa - be * pb
    M3 = 2.0 * M1 + al ** 2 * pa - be ** 2 * pb
    M4 = 3.0 * M2 + al ** 3 * pa - be ** 3 * pb
    return M0, M1, M2, M3, M4


def _negll(mu, s, n, sx, sxx, a, b):
    """R ll.tnorm2 (returns MINUS the truncated-normal log-likelihood)."""
    ssq = sxx - 2.0 * mu * sx + n * mu * mu
    Z = ndtr((b - mu) / s) - ndtr((a - mu) / s)
    with np.errstate(divide="ignore", invalid="ignore"):
        ll = -n * np.log(Z) - n * np.log(_SQ2PI * s) - ssq / (2.0 * s * s)
    return -ll


def _grad_hess(mu, s, n, sx, sxx, a, b, typo):
    """R grad.tnorm and hessian.tnorm with the integrals int_a^b dpsi in closed form."""
    M0, M1, M2, M3, M4 = _moments(mu, s, a, b)
    ssq = sxx - 2.0 * mu * sx + n * mu * mu
    dm = n * mu - sx
    Z = M0
    Zm = M1 / s
    Zs = (M2 - M0) / s
    Zmm = (M2 - M0) / s ** 2
    Zms = (M3 - 3.0 * M1) / s ** 2
    Zss = (2.0 * M0 + M4) / s ** 2 - 5.0 * M1 / s ** 3 if typo else (2.0 * M0 - 5.0 * M2 + M4) / s ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        g1 = -n * Zm / Z - dm / s ** 2
        g2 = -n * Zs / Z - n / s + ssq / s ** 3
        h1 = -n * (Z * Zmm - Zm ** 2) / Z ** 2 - n / s ** 2
        h3 = -n * (Z * Zms - Zm * Zs) / Z ** 2 + 2.0 * dm / s ** 3
        h2 = -n * (Z * Zss - Zs ** 2) / Z ** 2 + n / s ** 2 - 3.0 * ssq / s ** 4
    return g1, g2, h1, h2, h3


def estimates(P, perc, iter_cap=50, typo=True):
    """EstimatesComputation(missingdata = P (proteins x cells), perc). Returns (mu, sd, used_tn: bool per protein)."""
    O = ~np.isnan(P)
    n = O.sum(1).astype(float)
    P0 = np.where(O, P, 0.0)
    sx = P0.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = sx / n
        Dv = np.where(O, P - mu[:, None], 0.0)  # centred observed values (two-pass, as R's sd)
        sxxc_all = (Dv * Dv).sum(1)
        sd = np.sqrt(sxxc_all / (n - 1))  # R sd (n - 1); NaN when n < 2
    sd = np.where(n >= 2, sd, np.nan)
    mu = np.where(n >= 1, mu, np.nan)
    nsamp = P.shape[1]
    case1 = (P.shape[1] - n) / nsamp >= perc
    lod = np.nanmin(P)
    with np.errstate(invalid="ignore"):
        case2 = mu > 3.0 * sd + lod  # NaN -> False, as R which()
    nr = ~case1 & ~case2 & (n >= 2) & np.isfinite(sd) & (sd > 0)
    used = np.zeros(len(mu), bool)
    if nr.any():
        upp = np.nanmax(P) + 2.0 * np.nanmax(sd)  # R: max(ParamEstim[,2]) (no na.rm); see module docstring
        # fit on values centred at the sample mean (the likelihood is translation invariant when the limits move
        # with the data), for accuracy of sum((x - mu)^2); undone afterwards
        idx = np.where(nr)[0]
        c = mu[idx]
        sxc = Dv[idx].sum(1)
        sxxc = sxxc_all[idx]
        m_hat, s_hat, it, ok = _tn_fit(np.zeros(len(idx)), sd[idx], n[idx], sxc, sxxc, lod - c, upp - c, typo)
        acc = ok & (it < iter_cap)
        mu[idx[acc]] = m_hat[acc] + c[acc]
        sd[idx[acc]] = s_hat[acc]
        used[idx[acc]] = True
    return mu, sd, used


def _tn_fit(mu0, s0, n, sx, sxx, a, b, typo=True, tol=1e-7, maxit=100):
    """Vectorised R NewtonRaphsonLike over proteins (per-protein limits a, b). Returns (mu, s, iter, ok);
    iter = maxit when not converged; ok = False on any R error path (singular Hessian, sigma <= 0, NaN)."""
    m = len(mu0)
    mu, s = mu0.astype(float).copy(), s0.astype(float).copy()
    it = np.full(m, maxit)
    ok = np.zeros(m, bool)
    active = np.isfinite(mu) & np.isfinite(s) & (s > 0) & (n >= 2)
    g1, g2, *_ = _grad_hess(mu, s, n, sx, sxx, a, b, typo)
    with np.errstate(invalid="ignore"):
        done0 = active & (np.abs(g1) + np.abs(g2) < tol)
    it[done0], ok[done0] = 0, True
    active &= ~done0
    for i in range(1, maxit + 1):
        if not active.any():
            break
        idx = np.where(active)[0]
        c_mu, c_s, nn, sxi, sxxi, ai, bi = mu[idx], s[idx], n[idx], sx[idx], sxx[idx], a[idx], b[idx]
        g1, g2, h1, h2, h3 = _grad_hess(c_mu, c_s, nn, sxi, sxxi, ai, bi, typo)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            det = h1 * h2 - h3 * h3
            sing = ~np.isfinite(det) | (np.abs(det) <= 2.2e-16 * (np.abs(h1 * h2) + h3 * h3))
            n_mu = c_mu - (h2 * g1 - h3 * g2) / det  # new = cur - g %*% solve(H), H symmetric
            n_s = c_s - (-h3 * g1 + h1 * g2) / det
            fail = sing | ~np.isfinite(n_mu) | ~np.isfinite(n_s) | (n_s <= 0)
            llc = _negll(c_mu, c_s, nn, sxi, sxxi, ai, bi)
            lln = _negll(n_mu, np.where(n_s > 0, n_s, np.nan), nn, sxi, sxxi, ai, bi)
            rel = np.abs(llc - lln) / llc  # R: signed denominator, reproduced
        fail |= ~np.isfinite(rel)
        conv = ~fail & (rel < tol)
        keep = ~fail & ~conv
        active[idx[fail]] = False
        c_idx = idx[conv]
        mu[c_idx], s[c_idx], it[c_idx], ok[c_idx] = n_mu[conv], n_s[conv], i, True
        active[c_idx] = False
        mu[idx[keep]], s[idx[keep]] = n_mu[keep], n_s[keep]
    return mu, s, it, ok


# --------------------------------------------------------------------------- correlation kNN (imputeKNN)


def _impute_knn(P, mu, sd, k, block=128):
    """imputeKNN core on P (proteins x cells) standardised by (mu, sd). Returns proteins x cells with imputed values at
    missing entries (NaN where nothing could be imputed) and NaN elsewhere."""
    n_prot, n_cell = P.shape
    O = ~np.isnan(P)
    good = np.isfinite(mu) & np.isfinite(sd) & (sd > 0) & (O.sum(1) >= 2)
    with np.errstate(invalid="ignore"):
        Z = (P - mu[:, None]) / sd[:, None]
    Z[~good] = np.nan
    Oc = ~np.isnan(Z)  # candidate observations
    Zc = np.where(Oc, Z, 0.0)
    Ofc = Oc.astype(float)
    Zc2 = Zc * Zc
    ZcT, OfcT, Zc2T = Zc.T.copy(), Ofc.T.copy(), Zc2.T.copy()
    out = np.full(P.shape, np.nan)
    targets = np.where(good & (~O).any(1))[0]
    kk = min(k, n_prot)
    for s0 in range(0, len(targets), block):
        blk = targets[s0 : s0 + block]
        T, To, T2 = Zc[blk], Ofc[blk], Zc2[blk]
        cnt = To @ OfcT  # co-observed cells, b x n_prot
        sx, sy = T @ OfcT, To @ ZcT
        sxx, syy, sxy = T2 @ OfcT, To @ Zc2T, T @ ZcT
        with np.errstate(divide="ignore", invalid="ignore"):
            vx = sxx - sx * sx / cnt
            vy = syy - sy * sy / cnt
            r = (sxy - sx * sy / cnt) / np.sqrt(vx * vy)
        bad = (cnt < 3) | ~np.isfinite(r) | (vx <= 1e-12 * cnt) | (vy <= 1e-12 * cnt)
        r = np.clip(r, -1.0, 1.0)  # R's cor() clamps to [-1, 1]
        r[bad] = np.nan
        r[np.arange(len(blk)), blk] = np.nan  # self (never a candidate: it is missing wherever it is imputed)
        dist = 1.0 - np.abs(r)
        dist[np.isnan(r) | (np.abs(r) == 1.0)] = np.inf
        sgn = np.sign(np.nan_to_num(r))
        for e in range(n_cell):
            rows = np.where(~O[blk, e])[0]
            if len(rows) == 0:
                continue
            D = dist[rows]
            D[:, ~Oc[:, e]] = np.inf  # candidates: proteins observed in cell e
            nb = np.argpartition(D, kk - 1, axis=1)[:, :kk] if kk < n_prot else np.tile(np.arange(n_prot), (len(rows), 1))
            d = np.take_along_axis(D, nb, 1)
            with np.errstate(divide="ignore"):
                w = np.where(np.isfinite(d), 1.0 / d, 0.0)
            tot = w.sum(1)
            w = w / np.where(tot > 0, tot, 1.0)[:, None] * np.take_along_axis(sgn[rows], nb, 1)
            zval = (w * np.nan_to_num(Zc[nb, e])).sum(1)
            ti = blk[rows]
            out[ti, e] = np.where(tot > 0, mu[ti] + sd[ti] * zval, np.nan)
    return out


def _run(Xm, ctx, k, perc, typo):
    X = _canon(Xm)
    P = X.T  # proteins x cells (NAguideR / imputeKNN orientation)
    mu, sd, used = estimates(P, perc, typo=typo)
    out = _impute_knn(P, mu, sd, k)
    ctx.setdefault("_trknn_diag", {})[f"perc={perc}"] = {
        "n_tn_fits_used": int(used.sum()),
        "n_proteins": int(P.shape[0]),
    }
    return _keep_obs(Xm, out.T, ctx)


def trknn(Xm, ctx, k=10):
    """NAguideR's trknn call: imputeKNN(k = 10, distance = 'truncation', perc = 0). With perc = 0 every protein takes the
    sample mean/SD (EstimatesComputation Case 1), so this equals KNN-CR; see the module docstring."""
    return _run(Xm, ctx, k, perc=0.0, typo=True)


def trknn_tn(Xm, ctx, k=10):
    """imputeKNN(k = 10, distance = 'truncation', perc = 1) (imputeKNN's default perc): truncated-normal mean/SD from
    Newton-Raphson for proteins not fully missing and not > 3 SD above the global LOD (= minimum observed value)."""
    return _run(Xm, ctx, k, perc=1.0, typo=True)


TRKNN = {"trknn": trknn, "trknn_tn": trknn_tn}
