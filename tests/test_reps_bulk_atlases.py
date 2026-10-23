"""Phase B bulk-atlas reference panels (reps.procan, reps.gtex_proteome, reps.cptac, reps.atlantic).

Skipped unless the raw downloads (scripts/fetch_features.sh, "Phase B bulk atlases" section) are in
data/cache/reps/. Expected values are derived from the snapshots in-test; the proteasome check uses PSMA1 vs
PSMA2 where the panel has PSMA2 (GTEx Table S2 and the CPTAC BCM tables do not), else PSMA1 vs PSMA4.
"""

from __future__ import annotations

import networkx as nx
import numpy as np
import pytest

from scp_graphfs_mcp.reps import atlantic, cptac, gtex_proteome, procan

PSMA1, PSMA2, PSMA4 = "P25786", "P25787", "P25789"
VARS = [f"A0A024R1R8;{PSMA1}", PSMA2, PSMA4, "P25788", "P28066", "P04637", "Q00987", "NOT_A_PROTEIN"]

PANELS = [
    pytest.param(procan, marks=pytest.mark.skipif(not procan.RAW.exists(), reason="ProCan download absent")),
    pytest.param(gtex_proteome, marks=pytest.mark.skipif(not gtex_proteome.RAW.exists(), reason="GTEx download absent")),
    pytest.param(cptac, marks=pytest.mark.skipif(len(cptac.raw_files()) == 0, reason="CPTAC download absent")),
    pytest.param(atlantic, marks=pytest.mark.skipif(not atlantic.RAW.exists(), reason="ATLANTiC download absent")),
]


def _r(a, b):
    m = ~np.isnan(a) & ~np.isnan(b)
    return m.sum(), np.corrcoef(a[m], b[m])[0, 1]


@pytest.mark.parametrize("mod", PANELS)
def test_reference(mod):
    R, has = mod.reference(VARS)
    X, index, samples, _ = mod._panel()
    assert X.shape[1] == len(samples) and len(samples) >= 100
    assert R.shape == (len(VARS), X.shape[1]) and has.dtype == bool
    assert has[0] and not has[-1] and np.isnan(R[-1]).all()
    assert np.isnan(R[~has]).all()
    assert np.allclose(np.nanmean(R[has], 1), 0, atol=1e-6)
    sd = np.nanstd(R[has], 1)
    assert np.allclose(sd[sd > 0], 1, atol=1e-6)
    partner = 1 if has[1] else 2
    assert has[partner]
    n, r = _r(R[0], R[partner])
    assert n >= 50 and r > 0.5                                  # proteasome subunits co-vary across samples


@pytest.mark.parametrize("mod", PANELS)
def test_graph_and_block(mod):
    g = mod.corr_graph(VARS, top=3)
    assert isinstance(g, nx.Graph) and set(g.nodes) == set(VARS)
    assert nx.number_of_selfloops(g) == 0
    assert all(w > 0 for _, _, w in g.edges(data="weight"))
    assert g.degree(VARS[0]) > 0 and g.degree("NOT_A_PROTEIN") == 0
    B = mod.block(VARS, d=8)
    assert B.shape == (len(VARS), 9) and np.isfinite(B).all()
    assert set(np.unique(B[:, -1])) <= {0.0, 1.0}
    assert B[-1, -1] == 0 and np.all(B[-1, :-1] == 0)


@pytest.mark.skipif(len(cptac.raw_files()) == 0, reason="CPTAC download absent")
def test_cptac_cohort_centering():
    R, has = cptac.reference(VARS, center="cohort")
    R0, has0 = cptac.reference(VARS)
    assert np.array_equal(has, has0) and R.shape == R0.shape
    _, _, _, cohort = cptac._panel()
    assert len(set(cohort)) == len(cptac.raw_files())
