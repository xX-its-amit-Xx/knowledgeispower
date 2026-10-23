"""Harness sanity checks: the champion itself (cand_stack must be ~1.0) and a pure-noise candidate (weight ~0)."""

import numpy as np


def impute(Xm, ctx):
    return ctx["champion"]


def impute_null(Xm, ctx):
    rng = np.random.default_rng(ctx["seed"])
    return ctx["base"] + rng.normal(0, 1, Xm.shape)
