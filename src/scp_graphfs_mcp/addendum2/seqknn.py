"""SeqKNN: sequential kNN imputation (post hoc, exploratory addendum; written 2026-11-08 after PREREGISTRATION_v9 was frozen).

Reference: Kim K-Y, Kim B-J, Yi G-S (2004) "Reuse of imputed data in microarray analysis increases imputation
efficiency." BMC Bioinformatics 5:160, doi:10.1186/1471-2105-5-160 (title, authors, journal, volume and article number
checked against Crossref on 2026-11-08). Reference code: CRAN package SeqKnn 1.0.1 (archived; Kim & Yi, 2008),
functions SeqKNN() and nnmiss(), read from cran.r-project.org/src/contrib/Archive/SeqKnn/SeqKnn_1.0.1.tar.gz on
2026-11-08. NAguideR (Wang et al. 2020, Nucleic Acids Res, doi:10.1093/nar/gkaa498) lists SeqKNN among its methods with
proteins as rows; the exact NAguideR call (k = 10 on the log2 protein x sample table) is taken from the task brief and
was NOT verified against NAguideR's source (its GitHub repository could not be fetched from this session).

Algorithm (a port of SeqKnn 1.0.1, rows = PROTEINS, columns = the ~30 cells, i.e. Xm.T):
  1. Rows with no missing value form the candidate pool ("xcomplete").
  2. Rows with missing values are visited in increasing order of their number of missing cells (R's order(), stable,
     so ties keep the original protein order).
  3. For the current row, the distance to every pool row is the SUM of squared differences over the cells the current
     row observes (nnmiss: scale(x, xmiss, FALSE)[, !ismiss]^2 summed; unscaled, not averaged, no square root).
  4. The k nearest pool rows (all of them if the pool has fewer than k) are averaged with weights
     1 / (d + 1e-15), normalised to sum to 1, at the cells the current row misses.
  5. The completed row is appended to the pool and so becomes a candidate neighbour for later rows ("reuse").
  No standardisation is applied, as in the R code: neighbours are matched on raw log2 values.

Deviations from SeqKnn 1.0.1 (edge cases the R code handles badly):
  - Proteins with no observed cell are not imputed here (R would average the first k pool rows, since all distances are
    0); they get ctx['base'] via baselines._keep_obs and are not added to the pool.
  - If the pool is empty when a row is reached (no complete protein), R's weighted mean of zero rows returns 0; here
    that row's missing cells get the row's own observed mean (per-protein mean over observed cells) and the row then
    enters the pool as usual.
  - NaN payloads are canonicalised first (see addendum2/simple.py: leduc2022 carries R's NA_real_ bit pattern).

Same signature as baselines.py: fn(Xm, ctx) -> full matrix; only unobserved entries are scored.
"""

from __future__ import annotations

import numpy as np

from ..baselines_addendum import _keep_obs
from .simple import _canon


def seqknn(Xm, ctx, k=10, eps=1e-15):
    X = _canon(Xm)
    P = X.T.copy()  # proteins x cells (SeqKnn / NAguideR orientation)
    miss = np.isnan(P)
    nmiss = miss.sum(1)
    n_prot, n_cell = P.shape
    complete = np.where(nmiss == 0)[0]
    todo = np.where((nmiss > 0) & (nmiss < n_cell))[0]
    todo = todo[np.argsort(nmiss[todo], kind="stable")]  # R: order(missing) over the incomplete rows

    pool = np.empty((len(complete) + len(todo), n_cell))
    pool[: len(complete)] = P[complete]
    n_pool = len(complete)
    out = np.full_like(P, np.nan)
    for r in todo:
        o = ~miss[r]
        x = P[r]
        if n_pool == 0:
            fill = np.full((~o).sum(), x[o].mean())
        else:
            diff = pool[:n_pool][:, o] - x[o]
            d = np.einsum("ij,ij->i", diff, diff)
            nb = np.argsort(d, kind="stable")[:k]
            w = 1.0 / (d[nb] + eps)
            w /= w.sum()
            fill = w @ pool[nb][:, ~o]
        row = x.copy()
        row[~o] = fill
        out[r] = row
        pool[n_pool] = row
        n_pool += 1
    return _keep_obs(Xm, out.T, ctx)


SEQKNN = {"seqknn": seqknn}
