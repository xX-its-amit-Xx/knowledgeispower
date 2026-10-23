"""Synthetic check for masked-autoencoder transfer (no real data)."""
import numpy as np
import pytest

torch = pytest.importorskip("torch")
from scp_graphfs_mcp import mae  # noqa: E402


def test_mae_transfers_structure_and_beats_scrambled_reference():
    rng = np.random.default_rng(0)
    p, r = 80, 4
    V = rng.normal(size=(p, r))
    ref = np.tanh(rng.normal(size=(1500, r)) @ V.T) + 0.1 * rng.normal(size=(1500, p))
    ref = (ref - ref.mean(0)) / ref.std(0)
    ref[rng.random(ref.shape) < 0.4] = np.nan
    X = np.tanh(rng.normal(size=(30, r)) @ V.T) + 0.1 * rng.normal(size=(30, p))
    M = rng.random(X.shape) < 0.3
    Xm = X.copy()
    Xm[M] = np.nan
    base = np.where(np.isnan(Xm), np.nanmean(Xm, 0), Xm)
    rm = lambda Y: np.sqrt(np.mean((Y[M] - X[M]) ** 2))  # noqa: E731
    good = rm(mae.mae_impute(Xm, mae.train_mae(ref, hidden=64, latent=8, epochs=40), base))
    scr = rm(mae.mae_impute(Xm, mae.train_mae(ref[:, rng.permutation(p)], hidden=64, latent=8, epochs=40), base))
    assert good < 0.8 * rm(base) and good < 0.8 * scr
