"""Regularised EM (RegEM) imputation (post hoc, exploratory addendum; written 2026-11-08 after PREREGISTRATION_v9 was frozen).

References (checked 2026-11-08):
  - Schneider T (2001) "Analysis of incomplete climate data: Estimation of mean values and covariance matrices and
    imputation of missing values." J Climate 14(5):853-871, doi:10.1175/1520-0442(2001)014<0853:AOICDE>2.0.CO;2
    (title, author, journal, volume, issue, pages checked against Crossref).
  - Webb-Robertson B-JM et al. (2015) "Review, evaluation, and discussion of the challenges of missing value imputation
    for mass spectrometry-based label-free global proteomics." J Proteome Res 14(5):1993-2001, doi:10.1021/pr501138h
    (checked against Crossref; full text read at PMC4776766). Their "REM" is Schneider's RegEM: "an iterative process of
    linear regression of variables (peptide intensities) with missing values on peptides without missing values.
    Regression coefficients are estimated by ridge regression ... determined by generalized cross-validation", run with
    Schneider's MATLAB code (their Table 2 cites http://www.clidyn.ethz.ch/imputation/index.html, ref. 19 = Schneider
    2001). Which orientation of their peptide x sample table they passed to regem.m is NOT stated in the text we read
    (their wording suggests peptides were the variables); unverified.
  - Reference code: Schneider's RegEM MATLAB package, github.com/tapios/RegEM (master branch, read 2026-11-08 via
    raw.githubusercontent.com): regem.m (driver), mridge.m (multiple ridge regression), gcvridge.m + gcvfctn.m (GCV),
    peigs.m, center.m, rescale.m, missingness_patterns.m. The functions below port regem.m with its DEFAULT options:
    regress='mridge' (one ridge parameter per missingness pattern, chosen by GCV), stagtol=1e-2, maxit=30,
    inflation=1, relvar_res=5e-2, minvarfrac=0, regularisation on the correlation matrix (variables rescaled to unit
    variance every iteration), initial fill = mean substitution, degrees of freedom n-1.

The algorithm (regem.m): centre each variable; fill missing with 0 (= its mean); C = X'X/(n-1). Each iteration: rescale
to unit variance; for every missingness pattern regress its missing variables on its available ones with
B = (C_aa + h^2 I)^-1 C_am (eigendecomposition of C_aa, h by GCV), impute X_m = X_a B and accumulate the ridge residual
covariance S = S0 + F' diag(h^4/(d+h^2)^2) F into CovRes; then re-centre and C = (X'X + CovRes)/(n-1). Stop when the RMS
change of the imputed values relative to their RMS (uncentred) value falls below stagtol, or after maxit iterations.

Which orientation (the task asked to justify it). Xm is cells x proteins with ~30 cells and 1000-5000 proteins.
  rem (primary)       records = PROTEINS, variables = the ~30 CELLS, i.e. regem.m on Xm.T. C is 30 x 30 and is estimated
                      from 1000-5000 records, so the EM is well posed, the port is exact (residual-covariance term
                      included) and it is cheap. A protein's missing cells are predicted by a GCV-ridge regression on its
                      observed cells, with coefficients learnt from all proteins: a regularised, learnt version of
                      "borrow from similar cells" (the kNN-over-cells reference does the same with a fixed 10-NN mean).
                      This is the same orientation as baselines_addendum.bpca_oba (proteins as observations).
  rem_protein_exact   records = CELLS, variables = PROTEINS (regem.m on Xm as given; the n << p case RegEM was designed
                      for). Needs the dense p x p covariance (5000^2 doubles = 200 MB) and one top-(n-1) eigensolve of a
                      p_a x p_a block per cell per iteration: ~30 x 30 large eigensolves per call. Exact but only
                      practical for small p; it refuses p > max_p (default 1500) and is used here for the synthetic
                      check only.
  rem_protein         the protein-space variant made tractable for p = 5000 with the dual (Gram-matrix) form: without
                      the residual-covariance term, C = X'X/(n-1) has rank <= n-1 and every quantity mridge/gcvridge
                      needs (eigenvalues, Fourier coefficients F = U'Y/sqrt(n-1), trace S0, the prediction
                      X_a B = U diag(d/(d+h^2)) U'Y) comes from the n x n Gram matrix of the available block.
                      DEVIATION: dropping CovRes makes it "regularised iterated regression" rather than EM (Schneider
                      2001 notes that omitting the conditional-covariance term biases the covariance estimate low).
                      The synthetic test compares it with rem_protein_exact to measure what that costs.
Choice: rem (cell orientation) is the primary method because it is the only one that is both an exact RegEM and cheap at
p = 5000; the two protein-space variants are kept so the orientation question is answered by data, not assumed.

Deviations common to all three:
  - Each protein is first centred on its observed-cell mean (baselines_addendum._stats) and the mean is added back at
    the end. In protein space this is a per-variable shift, to which RegEM is invariant EXCEPT for its stopping rule:
    regem.m measures the change of the imputed values relative to their uncentred RMS, so on log2 intensities (~20)
    the default stagtol = 1e-2 stopped after 3-4 iterations on synthetic data with N(20, 2^2) protein offsets
    (rem_protein_exact RMSE 0.88 / 1.06 at 30% / 60% MCAR vs 0.58 / 0.71 without offsets, one seed). In the cell
    orientation the centring also stops between-protein abundance differences dominating the 30 x 30 cell covariance.
    Per-protein SCALING was tried for rem and rejected on synthetic data (RMSE 0.61 / 0.81 vs 0.58 / 0.70 centring
    only, one seed): it inflates low-signal proteins' noise.
  - Records with no available value make regem.m stop with an error; here they are left out of the EM and get
    ctx['base'] via baselines._keep_obs (for rem: proteins observed in no cell; for the protein variants: cells with no
    observed protein). Variables observed in fewer than 2 records are left out the same way.
  - MATLAB's fminbnd is replaced by scipy.optimize.minimize_scalar(method='bounded') (the same Brent bounded minimiser,
    xatol = regem's TolX h_tol); peigs' Arnoldi branch (eigs) by scipy.linalg.eigh(subset_by_index) for the same top
    eigenpairs.
  - NaN payloads are canonicalised first (addendum2/simple.py: leduc2022 carries R's NA_real_ bit pattern).

Checks run 2026-11-08 (nice -n 19, one thread; scratch scripts, not part of the repo):
  - rem_protein (dual) equals rem_protein_exact run with inflation=0 (CovRes off) to 1e-13 on a 30 x 200 example.
  - regem on a well-posed case (2000 records x 4 variables, 20% MCAR) agrees with plain unregularised MVN EM to
    0.028 (max abs imputed difference; covariance to 0.003): the difference is the small GCV ridge.
  - Synthetic, rank 3 + N(0, 0.5^2) noise, 30 x 600, 3 seeds, RMSE on masked entries (signal-only oracle ~0.50):
        MCAR 30%: mean 1.78, kNN-over-cells 1.32, rem 0.59, rem_protein 0.58, rem_protein_exact 0.58
        MCAR 60%: mean 1.82, kNN-over-cells 1.67, rem 0.71, rem_protein 0.72, rem_protein_exact 0.74
    (same ranking with N(20, 2^2) protein offsets added). CPU per call: rem 2-4 s, rem_protein 0.2-0.5 s,
    rem_protein_exact ~5 s at p = 600.
  - One dev split (leduc2022, 30 cells, seed 100, MCAR 10% of entries of proteins observed in >= 50% of cells;
    30 x 2844, 2140 masked): RMSE relative to phase2 kNN: rem 1.003, rem_protein 1.045, per-protein mean 1.047.
    rem 9.9 s CPU. rem_protein_exact not run (2453 usable proteins > max_p).

Same signature as baselines.py: fn(Xm, ctx) -> full matrix; only unobserved entries are scored.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar

from ..baselines_addendum import _keep_obs, _stats
from .simple import _canon

_EPS = np.finfo(float).eps


def _peigs(A, rmax):
    """peigs.m: positive eigenpairs of symmetric A, descending; top-rmax only when rmax < min(size)/10 (eigs branch)."""
    m = A.shape[0]
    rmax = min(rmax, m)
    if rmax < m / 10:
        d, V = eigh(A, subset_by_index=[m - rmax, m - 1])
    else:
        d, V = np.linalg.eigh(A)
    o = np.argsort(-d, kind="stable")
    d, V = d[o], V[:, o]
    r = int(np.sum(d > d.max() * m * _EPS)) if d.size else 0
    return V[:, :r], d[:r], r


def _gcvridge(F, d, trS0, n, r, trSmin, minvarfrac=0.0):
    """gcvridge.m + gcvfctn.m: GCV-optimal ridge parameter h for the regularised normal equations."""
    fc2 = np.sum(F ** 2, 1)
    h_tol = 0.2 / np.sqrt(n)
    varfrac = np.cumsum(d) / np.sum(d)
    if minvarfrac > varfrac.min():
        h_max = np.sqrt(np.interp(minvarfrac, varfrac, d))
    else:
        h_max = np.sqrt(d.max()) / h_tol
    if trS0 > trSmin:
        h_min = np.sqrt(_EPS)
    else:
        rtsvd = np.zeros(r)
        rtsvd[r - 1] = trS0
        for j in range(r - 2, -1, -1):
            rtsvd[j] = rtsvd[j + 1] + fc2[j + 1]
        rmin = int(np.argmin(np.abs(rtsvd - trSmin)))
        h_min = np.sqrt(max(d[rmin], d.min() / n))
    if not h_min < h_max:
        return h_min
    dd, ff, dof0 = d[:r], fc2[:r], n - r

    def g(h):
        fil = h * h / (dd + h * h)
        return (np.sum(fil ** 2 * ff) + trS0) / (dof0 + np.sum(fil)) ** 2

    return float(minimize_scalar(g, bounds=(h_min, h_max), method="bounded", options={"xatol": h_tol}).x)


def _mridge(Cxx, Cyy, Cxy, dof, relvar_res=5e-2, minvarfrac=0.0):
    """mridge.m: B = (Cxx + h^2 I)^-1 Cxy with h by GCV; residual covariance S; peff = trace(Mxx_h Cxx)."""
    px, py = Cxx.shape[0], Cyy.shape[0]
    V, d, r = _peigs(Cxx, min(dof, px))
    if r == 0:
        return np.zeros((px, py)), Cyy.copy(), np.nan, 0.0
    F = (V / np.sqrt(d)).T @ Cxy
    S0 = Cyy - F.T @ F if dof > r else np.zeros((py, py))
    h = _gcvridge(F, d, np.trace(S0), dof, r, relvar_res * np.trace(Cyy), minvarfrac)
    B = V @ ((np.sqrt(d) / (d + h * h))[:, None] * F)
    S = S0 + F.T @ ((h ** 4 / (d + h * h) ** 2)[:, None] * F)
    return B, S, h, float(np.sum(d / (d + h * h)))


def _patterns(avail):
    pats, inv = np.unique(avail, axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    order = np.argsort(inv, kind="stable")
    bounds = np.searchsorted(inv[order], np.arange(len(pats) + 1))
    return [(pats[j], order[bounds[j]:bounds[j + 1]]) for j in range(len(pats))]


def regem(X, maxit=30, stagtol=1e-2, inflation=1.0, relvar_res=5e-2, minvarfrac=0.0):
    """Port of regem.m (regress='mridge', GCV) on a records x variables matrix with NaN = missing.
    Every record needs >= 1 available value and every variable >= 1. Returns (completed X, mean M, C, n_iter)."""
    X = np.array(X, float)
    n, p = X.shape
    dofC = n - 1
    miss = np.isnan(X)
    nmis = int(miss.sum())
    if nmis == 0:
        return X, X.mean(0), np.cov(X, rowvar=False), 0
    kmis = np.nonzero(miss)[1]
    pats = [(~pa, rows) for pa, rows in _patterns(~miss) if (~pa).any()]
    M = np.nanmean(X, 0)
    X = X - M
    X[miss] = 0.0
    C = X.T @ X / dofC
    it, rd = 0, np.inf
    while it < maxit and rd > stagtol:
        it += 1
        CovRes = np.zeros((p, p))
        D = np.sqrt(np.diag(C)).copy()
        D[np.abs(D) < _EPS] = 1.0
        X = X / D
        C = C / D[None, :] / D[:, None]
        Xmis = X.copy()
        for km, rows in pats:
            ka = ~km
            B, S, _, _ = _mridge(C[np.ix_(ka, ka)], C[np.ix_(km, km)], C[np.ix_(ka, km)], n - 1, relvar_res, minvarfrac)
            Xmis[np.ix_(rows, km)] = X[np.ix_(rows, ka)] @ B
            S = inflation * S * D[km][None, :] * D[km][:, None]
            CovRes[np.ix_(km, km)] += len(rows) * S
        X = X * D
        Xmis = Xmis * D
        dX = np.linalg.norm(Xmis[miss] - X[miss]) / np.sqrt(nmis)
        nX = np.linalg.norm(X[miss] + M[kmis]) / np.sqrt(nmis)
        rd = np.inf if nX < _EPS else dX / nX
        X[miss] = Xmis[miss]
        Mup = X.mean(0)
        X = X - Mup
        M = M + Mup
        C = (X.T @ X + CovRes) / dofC
    return X + M, M, C, it


def _usable(X, rec_axis):
    """Masks (records, variables) that regem can use: variables with >= 2 observed records, records with >= 1
    available value among those variables. rec_axis = 0 if records are rows of X, 1 if columns."""
    obs = ~np.isnan(X)
    var_ok = obs.sum(rec_axis) >= 2
    rec_ok = (obs[:, var_ok] if rec_axis == 0 else obs[var_ok, :]).sum(1 - rec_axis) >= 1
    return rec_ok, var_ok


def _centre(Xm):
    """Per-protein centring on the observed-cell mean (baselines_addendum._stats); see module docstring."""
    X = _canon(Xm)
    _, mu, _, _ = _stats(X)
    return X - mu, mu


def rem(Xm, ctx, **kw):
    """RegEM with records = proteins, variables = cells (regem.m on Xm.T), default regem.m options."""
    X, mu = _centre(Xm)
    prot_ok, cell_ok = _usable(X, rec_axis=1)
    T = X[np.ix_(cell_ok, prot_ok)].T  # proteins x cells
    out = np.full(X.shape, np.nan)
    if T.shape[0] >= 3 and T.shape[1] >= 2:
        Y, _, _, _ = regem(T, **kw)
        out[np.ix_(cell_ok, prot_ok)] = Y.T
    return _keep_obs(Xm, mu + out, ctx)


def rem_protein_exact(Xm, ctx, max_p=1500, **kw):
    """RegEM with records = cells, variables = proteins (regem.m on Xm), dense p x p covariance."""
    X, mu = _centre(Xm)
    cell_ok, prot_ok = _usable(X, rec_axis=0)
    if prot_ok.sum() > max_p:
        raise ValueError(f"rem_protein_exact: {int(prot_ok.sum())} proteins > max_p={max_p}; use rem_protein")
    out = np.full(X.shape, np.nan)
    Y, _, _, _ = regem(X[np.ix_(cell_ok, prot_ok)], **kw)
    out[np.ix_(cell_ok, prot_ok)] = Y
    return _keep_obs(Xm, mu + out, ctx)


def rem_protein(Xm, ctx, maxit=30, stagtol=1e-2, relvar_res=5e-2, minvarfrac=0.0):
    """Protein-space RegEM in dual (n x n Gram) form WITHOUT the residual-covariance term (see module docstring).
    Same per-pattern GCV ridge, rescaling, convergence rule and defaults as regem.m."""
    X0, mu = _centre(Xm)
    cell_ok, prot_ok = _usable(X0, rec_axis=0)
    X = X0[np.ix_(cell_ok, prot_ok)]
    n, p = X.shape
    dof = n - 1
    miss = np.isnan(X)
    nmis = int(miss.sum())
    out = np.full(X0.shape, np.nan)
    if nmis and n >= 3:
        kmis = np.nonzero(miss)[1]
        pats = [(~pa, rows) for pa, rows in _patterns(~miss) if (~pa).any()]
        M = np.nanmean(X, 0)
        X = X - M
        X[miss] = 0.0
        it, rd = 0, np.inf
        while it < maxit and rd > stagtol:
            it += 1
            D = np.sqrt(np.sum(X ** 2, 0) / dof)
            D[np.abs(D) < _EPS] = 1.0
            Xs = X / D
            Xmis = Xs.copy()
            for km, rows in pats:
                ka = ~km
                A, Yb = Xs[:, ka], Xs[:, km]
                K = A @ A.T / dof
                dK, U = np.linalg.eigh(K)
                o = np.argsort(-dK, kind="stable")
                dK, U = dK[o], U[:, o]
                r = int(np.sum(dK > dK.max() * ka.sum() * _EPS))
                r = min(r, dof, int(ka.sum()))  # peigs keeps at most rmax = min(dof, p_a) pairs
                d, U = dK[:r], U[:, :r]
                F = U.T @ Yb / np.sqrt(dof)
                trCyy = np.sum(Yb ** 2) / dof
                trS0 = trCyy - np.sum(F ** 2) if dof > r else 0.0
                h = _gcvridge(F, d, trS0, dof, r, relvar_res * trCyy, minvarfrac)
                Xmis[np.ix_(rows, km)] = (U[rows] * (d / (d + h * h))) @ (U.T @ Yb)
            Xmis = Xmis * D
            dX = np.linalg.norm(Xmis[miss] - X[miss]) / np.sqrt(nmis)
            nX = np.linalg.norm(X[miss] + M[kmis]) / np.sqrt(nmis)
            rd = np.inf if nX < _EPS else dX / nX
            X[miss] = Xmis[miss]
            Mup = X.mean(0)
            X = X - Mup
            M = M + Mup
        out[np.ix_(cell_ok, prot_ok)] = X + M
    return _keep_obs(Xm, mu + out, ctx)
