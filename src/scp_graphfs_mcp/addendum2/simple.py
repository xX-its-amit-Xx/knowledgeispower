"""Simple single-value imputers (post hoc, exploratory addendum; written 2026-11-08 after PREREGISTRATION_v9 was frozen).

These are the trivial fills that most proteomics imputation benchmarks list as reference points, e.g. NAguideR
(Wang et al. 2020, Nucleic Acids Res, doi:10.1093/nar/gkaa498: "zero", "minimum", "colmedian", "rowmedian", ...) and
Harris et al. 2023 (J Proteome Res, doi:10.1021/acs.jproteome.3c00205). They are included so the panel shows how much any
model beats a constant. Same signature as baselines.py: fn(Xm, ctx) -> full matrix; only unobserved entries are scored.

Orientation: Xm is cells x proteins, log2 intensity, NaN = missing. "Per protein" = per column (over the ~30 cells),
"per cell" = per row ("sample" in the bulk-proteomics literature, where samples are columns of a protein x sample table).

  zero          fill 0 on the log2 scale (i.e. linear intensity 1). This is the literal "zero" baseline; on log2 data it
                is far below every observed value and is expected to be the worst comparator. It exists only as a
                reference row; it does NOT violate the repo invariant (data are never stored with 0 for NaN).
  median        per-protein median of the observed cells.
  protein_min   per-protein minimum of the observed cells.
  half_min      per-protein half-minimum: half the minimum linear intensity, i.e. min - log2(2) = min - 1 on log2.
  sample_min    per-cell minimum of that cell's observed proteins.

Proteins (or cells, for sample_min) with no observed value get NaN here and are replaced by ctx['base'] (the kNN fill)
via baselines._keep_obs, which also restores every observed entry.
"""

from __future__ import annotations

import warnings

import numpy as np

from ..baselines import _keep_obs


def _col(Xm, fn):
    with warnings.catch_warnings():  # all-NaN columns -> NaN -> ctx['base'] in _keep_obs
        warnings.simplefilter("ignore", RuntimeWarning)
        return fn(Xm, axis=0)


def zero(Xm, ctx):
    return _keep_obs(Xm, np.zeros_like(Xm, dtype=float), ctx)


def median(Xm, ctx):
    v = _col(Xm, np.nanmedian)
    return _keep_obs(Xm, np.broadcast_to(v, Xm.shape), ctx)


def protein_min(Xm, ctx):
    v = _col(Xm, np.nanmin)
    return _keep_obs(Xm, np.broadcast_to(v, Xm.shape), ctx)


def half_min(Xm, ctx):
    v = _col(Xm, np.nanmin) - np.log2(2.0)
    return _keep_obs(Xm, np.broadcast_to(v, Xm.shape), ctx)


def sample_min(Xm, ctx):
    v = _col(Xm.T, np.nanmin)
    return _keep_obs(Xm, np.broadcast_to(v[:, None], Xm.shape), ctx)


SIMPLE = {
    "zero": zero,
    "median": median,
    "protein_min": protein_min,
    "half_min": half_min,
    "sample_min": sample_min,
}
