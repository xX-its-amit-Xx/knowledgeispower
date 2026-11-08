"""Exact R reference imputers through scripts/r_bridge.R (post hoc, exploratory; written 2026-11-08 after PREREGISTRATION_v9).

These call the actual R/Bioconductor packages rather than Python re-implementations, so the Python ports in
baselines_addendum.py (bpca_oba, lls, knn_protein, qrilc) can be checked against the reference code. Orientation follows
MsCoreUtils::impute_matrix (the wrapper behind MSnbase::impute): the matrix handed to R is proteins x cells.

  r_bpca            pcaMethods::pca(method = "bpca", nPcs = ncells - 1) on proteins x cells (MsCoreUtils::impute_bpca)
  r_lls             pcaMethods::llsImpute(k = 10) on cells x proteins (proteins are the "genes"/variables)
  r_lls_onepass     NOT a pcaMethods default: llsImpute(k = 10, center = TRUE, maxSteps = 1). The defaults (center =
                    FALSE, maxSteps = 100, unregularised ginv) diverge on ~30-cell matrices (see the addendum report)
  r_knn             impute::impute.knn(k = 10) on proteins x cells, rng.seed = ctx["seed"] (MsCoreUtils::impute_knn)
  r_qrilc           imputeLCMD::impute.QRILC(tune.sigma = 1) on proteins x cells (MsCoreUtils::impute_QRILC)
  r_missforest      missForest::missForest defaults on proteins x cells (MsCoreUtils::impute_RF orientation:
                    proteins are the observations), ranger backend, num.threads = 1
  r_missforest_cells  missForest::missForest defaults with cells as observations and proteins as variables
  r_msimpute        msImpute::msImpute(method = "v2-mnar") with ONE group (the benchmark has no labels at imputation
                    time; `group` is required by v2-mnar, so this is an assumption, not a package default)
  r_msimpute_v2     msImpute::msImpute(method = "v2") (MAR low-rank part only)

Packages tested with (Ubuntu 24.04, R 4.3.3): pcaMethods 1.94.0, impute 1.76.0 (Debian binaries); imputeLCMD 2.1,
norm 1.0-11.1, missForest 1.6.1, ranger 0.16.0 (CRAN); msImpute 1.12.0 (Bioconductor 3.18).

Preprocessing done here, not in R: proteins observed in fewer than `min_obs` cells (method-specific; pcaMethods refuses
all-NA rows/columns, msImpute refuses proteins with < 4 observations) and cells with no observation are left out of the
R call and get ctx["base"] (the kNN fill) via _keep_obs. R runs single-threaded (OMP/OpenBLAS threads = 1). The R
executable is $RSCRIPT or `Rscript` on PATH.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from ..baselines import _keep_obs

R_SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "r_bridge.R"
_MIN_OBS = {"msimpute": 4, "msimpute_v2": 4, "lls": 2, "lls_onepass": 2, "missforest_cells": 1}


def r_available() -> bool:
    return shutil.which(os.environ.get("RSCRIPT", "Rscript")) is not None


def _run_r(Xm, ctx, method, timeout=3600, min_obs=None):
    Xm = np.asarray(Xm, float)
    obs = ~np.isnan(Xm)
    cols = obs.sum(0) >= (min_obs or _MIN_OBS.get(method, 1))
    rows = obs[:, cols].any(1)
    Y = np.full(Xm.shape, np.nan)
    if cols.sum() < 2 or rows.sum() < 3:
        return _keep_obs(Xm, Y, ctx)
    sub = Xm[np.ix_(rows, cols)]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with tempfile.TemporaryDirectory() as d:
        src, dst = Path(d) / "in.csv", Path(d) / "out.csv"
        np.savetxt(src, sub, delimiter=",", fmt="%.17g")  # NaN is written as "nan"; r_bridge.R reads it as NA
        cmd = [os.environ.get("RSCRIPT", "Rscript"), str(R_SCRIPT), str(src), method, str(int(ctx["seed"])), str(dst)]
        p = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout)
        if p.returncode != 0:
            raise RuntimeError(f"r_bridge.R {method} failed:\n{p.stderr[-2000:]}")
        out = np.genfromtxt(dst, delimiter=",", missing_values="NA", filling_values=np.nan)
    out = np.atleast_2d(out)
    if out.shape != sub.shape:
        raise RuntimeError(f"r_bridge.R {method}: shape {out.shape} != {sub.shape}")
    Y[np.ix_(rows, cols)] = out
    return _keep_obs(Xm, Y, ctx)


def r_bpca(Xm, ctx):
    return _run_r(Xm, ctx, "bpca")


def r_lls(Xm, ctx):
    return _run_r(Xm, ctx, "lls")


def r_lls_onepass(Xm, ctx):
    return _run_r(Xm, ctx, "lls_onepass")


def r_knn(Xm, ctx, rowmax=0.5):
    """impute.knn defaults. impute's knnimp() mean-imputes (per-cell mean over the other proteins) every protein missing in
    more than trunc(rowmax * n) cells and runs kNN on the rest; we do that split here, in exactly that way, before calling
    R, because impute.knn's colmax = 0.8 guard is evaluated on ALL rows and stops on SCP data (cells > 80% missing). So
    the result equals impute.knn(k = 10) except that the colmax check sees only the proteins kNN actually uses."""
    Xm = np.asarray(Xm, float)
    obs = ~np.isnan(Xm)
    n = Xm.shape[0]
    keep = (n - obs.sum(0)) <= int(rowmax * n)
    Y = _run_r(Xm, ctx, "knn", min_obs=max(1, n - int(rowmax * n)))
    xbar = np.array([np.nanmean(r) if np.isfinite(r).any() else np.nan for r in np.where(keep[None, :], Xm, np.nan)])
    drop = ~keep[None, :] & ~obs
    Y = np.where(drop, xbar[:, None], Y)
    return _keep_obs(Xm, Y, ctx)


def r_qrilc(Xm, ctx):
    return _run_r(Xm, ctx, "qrilc")


def r_missforest(Xm, ctx):
    return _run_r(Xm, ctx, "missforest")


def r_missforest_cells(Xm, ctx):
    return _run_r(Xm, ctx, "missforest_cells")


def r_msimpute(Xm, ctx):
    return _run_r(Xm, ctx, "msimpute")


def r_msimpute_v2(Xm, ctx):
    return _run_r(Xm, ctx, "msimpute_v2")


RBRIDGE = {"r_bpca": r_bpca, "r_lls": r_lls, "r_knn": r_knn, "r_qrilc": r_qrilc, "r_missforest": r_missforest,
           "r_msimpute": r_msimpute, "r_msimpute_v2": r_msimpute_v2}
VARIANTS = {"r_lls_onepass": r_lls_onepass}
HEAVY = {"r_missforest_cells": r_missforest_cells}
