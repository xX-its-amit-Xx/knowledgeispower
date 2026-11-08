"""impSeq / impSeqRob: sequential imputation (post hoc, exploratory addendum; written 2026-11-08 after PREREGISTRATION_v9
was frozen). They change no preregistered verdict.

References
  impSeq:    Verboven S, Vanden Branden K, Goos P (2007) "Sequential imputation for missing values." Computational Biology
             and Chemistry 31(5-6):320-327, doi:10.1016/j.compbiolchem.2007.07.001 (title, authors, journal, volume, issue
             and pages checked against Crossref on 2026-11-08).
  impSeqRob: rrcovNA's impSeqRob() also cites only Verboven et al. 2007 (its "robust seqimpute" reference line is left
             as "......" in the source). The robust variant is described in Vanden Branden K, Verboven S (2009) "Robust
             data imputation." Comput Biol Chem 33(1):7-13, doi:10.1016/j.compbiolchem.2008.07.019 (Crossref record found
             2026-11-08; the paper's text was NOT read, so this port follows the package code, not the paper).
  Code ported: CRAN rrcovNA 0.5-3 (published 2025-04-22; source tarball rrcovNA_0.5-3.tar.gz, SHA-256
             ad5445c2f8ef5dfe092d4e2ca0faabf3a587347e7cd3a72caa0399b4dfddca80, read 2026-11-08): R/impSeq.R,
             R/impSeqRob.R (incl. .covSD, .extradir, .update_invcov). Helpers read from robustbase 0.99-7
             (h.alpha.n in R/covMcd.R, rankMM in R/classPC.R) and rrcov 1.7-7 (vecnorm in R/utils.R); the EM fallback
             mirrors norm::prelim.norm/em.norm/imp.norm (norm 1.0-11.1 as installed here; em.norm defaults
             maxits = 1000, criterion = 1e-4 on the max absolute parameter change, on norm's centred-and-scaled data).
             rrcovNA could NOT be installed in this session (installation was not permitted), so R was not run on the
             benchmark matrices. Check done instead: the example in the headers of R/impSeq.R and R/impSeqRob.R
             (robustbase's `phosphor` data, columns 1:2, x[10,2] and x[15,1] set to NA; values documented there:
             impSeq 41.02312 / 16.60647, impSeqRob 37.51664 / 26.38606) is reproduced to all printed digits by this
             port (41.023122 / 16.606471 and 37.516642 / 26.386060). That example only exercises the n >> p path
             (16 complete rows, no ridge, no fallback, all-pairs directions); the ridge path, the 500-direction LCG
             path and the EM fallback are NOT checked against R and their equality with R is unverified.

Orientation (NAguideR, verified): NAguideR (Wang et al. 2020, Nucleic Acids Res, doi:10.1093/nar/gkaa498), file
inst/NAguideRapp/app.R of github.com/wangshisheng/NAguideR (commit ae26a5344c5a042e066ea3e127132a6be4bfb7a4, read
2026-11-08), calls `impSeq(df1)` and `impSeqRob(df1, alpha = 0.9)$x` on the data frame whose ROWS are proteins and
COLUMNS are samples (the same df1 is passed unchanged to impute.knn, i.e. genes-as-rows, and transposed for llsImpute).
So rrcovNA's "observations" are the PROTEINS and its "variables" are the cells. Here that is P = Xm.T: n = 1000-5000
observations in p ~ 30 dimensions, the regime the method was designed for (n >> p). The covariance is the p x p
cell-by-cell covariance, so NO extra regularisation is needed in this orientation; the only ridge is rrcovNA's own
(impSeq adds 0.01 * I while the running complete set has <= p rows). The cells-as-observations orientation (p =
1000-5000 protein variables, 30 observations) is not implemented: it is not how NAguideR uses the method and would need a
heavily regularised 5000 x 5000 covariance, i.e. a different method.

impseq (rrcovNA::impSeq, Verboven et al. 2007)
  1. Rows (proteins) are sorted by their number of missing cells (stable sort, as R's sort.int(method = "radix")).
  2. Mean and sample covariance (ddof 1) of the complete rows start the set.
  3. Each incomplete row, in that order, gets the conditional Gaussian mean of its missing cells given its observed
     cells, x_m = mu_m + S_mo S_oo^-1 (x_o - mu_o) (R writes it with the precision matrix,
     mu_m - Lam_mm^-1 Lam_mo (x_o - mu_o), an algebraic identity), where S = cov + 0.01 I while the set has <= p rows
     (R: `if (p >= length(complobs)) solve(covx + 0.01 * diag(p))`).
  4. The completed row joins the set; the mean and covariance are updated with R's recursion, which equals the exact
     sample covariance of the enlarged set (checked algebraically: (n-2)/(n-1) S + d d'/n with d = x - old mean).
  Rows with no observed cell get the current mean, as in R (in the benchmark they are never scored).

impseq_rob (rrcovNA::impSeqRob, alpha = 0.9 as in NAguideR)
  1. Same sort. h = robustbase::h.alpha.n(alpha, n_complete, p); at least ceiling(p / alpha) complete rows are needed
     and their robust covariance must have full rank p (rankMM), unless n_complete >= 5p.
  2. Stahel-Donoho outlyingness of the complete rows over directions through pairs of complete rows (all pairs if
     there are <= 500, else 500 pairs drawn by the package's deterministic LCG `uniran`, seed 0, ported exactly):
     max over directions of |y - median(y)| / mad(y) (mad with R's 1.4826 constant). The h least outlying complete
     rows give the initial classical mean and covariance.
  3. Each incomplete row gets the conditional mean as in impseq (no ridge; R inverts covx directly), then its
     outlyingness on the SAME directions/medians/MADs. If it is <= the h-th smallest outlyingness so far, it joins the
     estimation set (h += 1) and the mean/covariance are updated with the package's recursion
     cov <- (n-2)/(n-1) cov + (x - old mean)(x - old mean)' / (n-1) (n = new set size; note: NOT the exact sample
     covariance, which would divide the rank-one term by n; ported as written). R keeps the inverse with a
     Sherman-Morrison update of exactly this matrix; here the covariance is kept and solved directly (same quantity).
  4. Not enough complete rows (or rank < p): R imputes the whole matrix once under the multivariate normal MLE
     (norm::em.norm, then one random draw by imp.norm with rngseed(1234567)) and copies the least-missing incomplete
     rows of that draw into the starting set until ceiling(p / alpha) rows are complete, adding one more row at a time
     until the robust covariance has rank p. Ported with a Python EM (same model and stopping rule) and one draw from
     numpy's default_rng(ctx['seed'] + 1234567), so the fallback rows cannot match R's draws.

Deviations (all at edge cases rrcovNA does not handle)
  - impseq with < 2 complete rows: R's cov() of 0 or 1 rows is NA and impSeq() fails. Here the (p + 1 - n_complete)
    least-missing incomplete rows are first completed with their conditional MEAN under the EM normal MLE (the same
    model impSeqRob falls back on, but deterministic) and added to the starting set, so it has p + 1 rows (full rank).
    With 30-60% MCAR and ~30 cells there is essentially never a complete protein, so on such data impseq is "EM-seeded
    sequential imputation"; the number of seeded rows is in the `info` dict (see _impseq_core).
  - Cells with no observed protein are excluded from the analysis (R would produce NA covariances); their entries get
    ctx['base'] via _keep_obs.
  - Directions whose MAD is 0 are dropped from the outlyingness (R would produce Inf/NaN and then mis-sort).
  - norm's prelim.norm encodes missingness patterns as sum(2^j), which overflows R integers for p > 31 variables
    (norm warns it "does not work reliably when the number of variables exceeds 30"), so R's impSeqRob fallback cannot
    run on matrices with more than ~31 cells. The Python EM has no such limit.
  - If even the fallback cannot reach rank p (fewer than ~p/alpha proteins with an observation), impseq_rob returns
    the EM conditional-mean imputation for all rows (flagged in info['fallback'] = 'em_only').
  - NaN payloads are canonicalised first (addendum2/simple.py: leduc2022 carries R's NA_real_ bit pattern).

Same signature as baselines.py: fn(Xm, ctx) -> full matrix; only unobserved entries are scored. Deterministic for a
given ctx['seed'] (impseq does not use the seed at all).
"""

from __future__ import annotations

import bisect
import math

import numpy as np

from ..baselines_addendum import _keep_obs
from .simple import _canon

_MAD_CONST = 1.4826  # R stats::mad default constant


# ----------------------------------------------------------------------------------------------- shared helpers


def _prep(Xm):
    """Canonical NaNs, proteins x usable-cells matrix P (rrcovNA orientation), the usable-cell mask, row sort order."""
    X = _canon(Xm)
    keep = ~np.isnan(X).all(1)  # cells (rows of X) with >= 1 observed protein
    P = X[keep].T.copy()  # proteins (observations) x cells (variables)
    nmiss = np.isnan(P).sum(1)
    order = np.argsort(nmiss, kind="stable")
    return X, keep, P, nmiss, order


def _finish(Xm, ctx, X, keep, Ps, order):
    out = np.empty_like(Ps)
    out[order] = Ps  # undo the sort
    Y = np.full(X.shape, np.nan)
    Y[keep] = out.T
    return _keep_obs(Xm, Y, ctx)


def _cond_mean(mx, S, x, miss):
    """mu_m + S_mo S_oo^-1 (x_o - mu_o); all-missing row -> mu."""
    o = ~miss
    if not o.any():
        return mx[miss].copy()
    rhs = np.linalg.solve(S[np.ix_(o, o)], x[o] - mx[o])
    return mx[miss] + S[np.ix_(miss, o)] @ rhs


def h_alpha_n(alpha, n, p):
    """robustbase::h.alpha.n (covMcd.R, robustbase 0.99-7)."""
    n2 = (n + p + 1) // 2
    return math.floor(2 * n2 - n + 2 * (n - n2) * alpha)


def rank_mm(A):
    """robustbase::rankMM: number of singular values >= max(dim) * eps * sv[0]."""
    sv = np.linalg.svd(A, compute_uv=False)
    if sv.size == 0 or sv[0] == 0:
        return 0
    return int(np.sum(sv >= max(A.shape) * np.finfo(float).eps * abs(sv[0])))


def em_norm(P, maxits=1000, criterion=1e-4, chunk=256):
    """Multivariate-normal MLE by EM for rows of P with NaNs (the model of norm::em.norm, MLE without prior).

    As norm does, the data are centred and scaled per column first and convergence is max |theta_new - theta_old|
    <= criterion on that scale (theta = mean and covariance). Rows with no observation are left out. Start: mean 0,
    identity covariance on the standardised scale. Returns (mu, Sigma) on the original scale and the iteration count.
    """
    use = ~np.isnan(P).all(1)
    Q = P[use]
    O = ~np.isnan(Q)
    m = np.nanmean(Q, 0)
    s = np.nanstd(Q, 0, ddof=1)
    s = np.where(np.isfinite(s) & (s > 0), s, 1.0)
    Z = np.where(O, (Q - m) / s, 0.0)
    n, p = Z.shape
    mu, S = np.zeros(p), np.eye(p)
    eye = np.eye(p)
    it = 0
    for it in range(1, maxits + 1):
        T1 = np.zeros(p)
        T2 = np.zeros((p, p))
        for a in range(0, n, chunk):
            Oc, Zc = O[a:a + chunk], Z[a:a + chunk]
            Of = Oc.astype(float)
            # A_i = D S D + (I - D), D = diag(o_i): its inverse holds S_oo^-1 in the observed block
            A = S[None] * (Of[:, :, None] * Of[:, None, :]) + eye[None] * (1.0 - Of)[:, None, :]
            Ainv = np.linalg.inv(A) * (Of[:, :, None] * Of[:, None, :])  # D A^-1 D (zero outside observed block)
            SK = np.einsum("jk,ikl->ijl", S, Ainv)  # S D A^-1 D  (p x p per row)
            xh = mu + np.einsum("ijl,il->ij", SK, Zc - mu)  # observed coordinates reproduce z_o exactly
            C = S[None] - np.einsum("ijl,lk->ijk", SK, S)  # conditional covariance (0 on observed block)
            T1 += xh.sum(0)
            T2 += np.einsum("ij,ik->jk", xh, xh) + C.sum(0)
        mu_new = T1 / n
        S_new = T2 / n - np.outer(mu_new, mu_new)
        S_new = 0.5 * (S_new + S_new.T)
        delta = max(np.max(np.abs(mu_new - mu)), np.max(np.abs(S_new - S)))
        mu, S = mu_new, S_new
        if delta <= criterion:
            break
    return m + s * mu, np.outer(s, s) * S, it


def _cond_draw(mx, S, x, miss, rng):
    """One draw of the missing coordinates from the conditional normal (norm::imp.norm semantics)."""
    o = ~miss
    mean = _cond_mean(mx, S, x, miss)
    if o.any():
        C = S[np.ix_(miss, miss)] - S[np.ix_(miss, o)] @ np.linalg.solve(S[np.ix_(o, o)], S[np.ix_(o, miss)])
    else:
        C = S[np.ix_(miss, miss)]
    C = 0.5 * (C + C.T)
    w, V = np.linalg.eigh(C)
    return mean + V @ (np.sqrt(np.clip(w, 0, None)) * rng.standard_normal(len(w)))


# ----------------------------------------------------------------------------------------------- impSeq


def _impseq_core(Ps, nmiss_s, ridge=0.01):
    """Ps: rows sorted by number of missing values. Returns completed Ps and an info dict."""
    Ps = Ps.copy()
    n, p = Ps.shape
    miss = np.isnan(Ps)
    nc = int(np.sum(nmiss_s == 0))
    info = {"n_complete": nc, "n_seeded": 0, "p": p}
    if nc == n:
        return Ps, info
    if nc < 2:  # R's impSeq fails here (cov of < 2 rows is NA); seed with EM conditional means (documented deviation)
        mu, S, its = em_norm(Ps)
        target = min(n, p + 1)
        for r in range(nc, target):
            Ps[r, miss[r]] = _cond_mean(mu, S, Ps[r], miss[r])
        info.update(n_seeded=target - nc, em_iterations=its)
        nc = target
    start = Ps[:nc]
    mx = start.mean(0)
    covx = np.cov(start, rowvar=False, ddof=1).reshape(p, p)
    ncompl = nc
    eye = np.eye(p)
    for r in range(nc, n):
        S = covx + ridge * eye if p >= ncompl else covx
        Ps[r, miss[r]] = _cond_mean(mx, S, Ps[r], miss[r])
        # R's update (applied at the start of its next iteration), == exact sample mean/covariance of the enlarged set
        ncompl += 1
        d = Ps[r] - mx
        mx = mx + d / ncompl
        covx = (ncompl - 2) / (ncompl - 1) * covx + np.outer(d, d) / ncompl
    return Ps, info


def impseq(Xm, ctx, ridge=0.01, return_info=False):
    X, keep, P, nmiss, order = _prep(Xm)
    Ps, info = _impseq_core(P[order], nmiss[order], ridge=ridge)
    Y = _finish(Xm, ctx, X, keep, Ps, order)
    return (Y, info) if return_info else Y


# ----------------------------------------------------------------------------------------------- impSeqRob


def _uniran(seed):
    seed = math.floor(seed * 5761) + 999
    quot = math.floor(seed / 65536)
    seed = math.floor(seed) - math.floor(quot * 65536)
    return seed, seed / 65536


def _randomset(n, k, seed):
    ranset = []
    for j in range(k):
        seed, u = _uniran(seed)
        num = math.floor(u * n) + 1
        if j > 0:
            while num in ranset:
                seed, u = _uniran(seed)
                num = math.floor(u * n) + 1
        ranset.append(num)
    return ranset, seed


def _extradir(data, ndirect, all_pairs):
    """rrcovNA:::.extradir: directions through pairs of rows (all pairs, or ndirect pairs from the package LCG)."""
    n = data.shape[0]
    if all_pairs:
        i, j = np.triu_indices(n, 1)  # same pairs as combn(n, 2), same order
        return data[i] - data[j]
    B = np.zeros((ndirect, data.shape[1]))
    seed = 0
    for r in range(ndirect):
        (a, b), seed = _randomset(n, 2, seed)
        B[r] = data[a - 1] - data[b - 1]
    return B


def _cov_sd(xc, h):
    """rrcovNA:::.covSD: Stahel-Donoho outlyingness; classical mean/cov of the h least outlying rows."""
    n = xc.shape[0]
    nrich1 = n * (n - 1) // 2
    ndirect = min(500, nrich1)
    B = _extradir(xc, ndirect, ndirect == nrich1)
    bn = np.sqrt((B ** 2).sum(1))
    B = B[bn > 1e-12]
    A = B / bn[bn > 1e-12][:, None]
    Y = xc @ A.T
    med = np.median(Y, 0)
    mad = _MAD_CONST * np.median(np.abs(Y - med), 0)
    good = mad > 0  # deviation: zero-MAD directions dropped
    A, med, mad = A[good], med[good], mad[good]
    d = (np.abs(Y[:, good] - med) / mad).max(1)
    ds = np.argsort(d, kind="stable")
    sub = xc[ds[:h]]
    return {"cov": np.cov(sub, rowvar=False, ddof=1).reshape(xc.shape[1], xc.shape[1]), "mean": sub.mean(0),
            "d": d, "A": A, "med": med, "mad": mad}


def _impseqrob_core(Ps, nmiss_s, alpha, rng):
    Ps = Ps.copy()
    n, p = Ps.shape
    miss = np.isnan(Ps)
    nc = int(np.sum(nmiss_s == 0))
    info = {"n_complete": nc, "n_seeded": 0, "p": p, "fallback": None}
    if nc == n:
        return Ps, info
    h = h_alpha_n(alpha, nc, p)
    mincomplete = math.ceil(p / alpha)
    r = None
    if nc >= 5 * p:
        r = p
    elif nc >= mincomplete:
        r = rank_mm(_cov_sd(Ps[:nc], h)["cov"])
    if nc < mincomplete or r < p:
        mu, S, its = em_norm(Ps)
        info.update(em_iterations=its, n_complete_before_seed=nc)
        while True:
            new = max(nc + 1, mincomplete)  # R: x[(ncomplobs+1):mincomplete, ] <- ximp[...] (one row if descending)
            if new > n:  # cannot reach rank p: EM conditional means for everything (deviation)
                for i in range(nc, n):
                    Ps[i, miss[i]] = _cond_mean(mu, S, Ps[i], miss[i])
                info["fallback"] = "em_only"
                return Ps, info
            for i in range(nc, new):
                Ps[i, miss[i]] = _cond_draw(mu, S, Ps[i], miss[i], rng)
            nc = new
            h = h_alpha_n(alpha, nc, p)
            r = rank_mm(_cov_sd(Ps[:nc], h)["cov"])
            if r >= p:
                break
            mincomplete += 1
        info["n_seeded"] = nc - info["n_complete"]
        info["fallback"] = "norm_draw"
    cs = _cov_sd(Ps[:nc], h)
    covx, mx, A, med, mad = cs["cov"], cs["mean"], cs["A"], cs["med"], cs["mad"]
    dsorted = sorted(cs["d"].tolist())
    ncu = h
    entered = 0
    for i in range(nc, n):
        Ps[i, miss[i]] = _cond_mean(mx, covx, Ps[i], miss[i])
        dx = float((np.abs(Ps[i] @ A.T - med) / mad).max())
        if dx <= dsorted[h - 1]:  # not an outlier: joins the estimation set
            ncu += 1
            h += 1
            entered += 1
            dlt = Ps[i] - mx
            mx = mx + dlt / ncu
            covx = (ncu - 2) / (ncu - 1) * covx + np.outer(dlt, dlt) / (ncu - 1)  # package recursion, as written
        bisect.insort(dsorted, dx)
    info["n_entered"] = entered
    return Ps, info


def impseq_rob(Xm, ctx, alpha=0.9, return_info=False):
    X, keep, P, nmiss, order = _prep(Xm)
    rng = np.random.default_rng(int(ctx.get("seed", 0)) + 1234567)
    Ps, info = _impseqrob_core(P[order], nmiss[order], alpha, rng)
    Y = _finish(Xm, ctx, X, keep, Ps, order)
    return (Y, info) if return_info else Y


IMPSEQ = {"impseq": impseq, "impseq_rob": impseq_rob}
