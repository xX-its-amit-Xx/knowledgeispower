"""pepev (Phase C, queen: new data modalities): per-entry precursor/peptide-count precision weights in a
heteroscedastic U-step of the domain-adapted projection.

Build (once per target, cached under data/cache/phasec/pepev/<target>_N.npz):

    python -m scp_graphfs_mcp.phasec.pepev build <target> [<target> ...]

Each build replays the corpus builder on the protein-level source to get X_rebuilt and carries the evidence count
matrix N through the same column / row selection, then asserts X_rebuilt == load_target(target)[0] (equal_nan) after
the builder's per-cell median centring. N = number of precursors / stripped peptides used for the (cell, protein)
value; where X is observed but N is 0/NaN (MBR / match transfer) N is set to 1.

Screen: impute / impute_null (within-cell permutation of the weights), plus variants impute_a1 (alpha = 1),
impute_uv (U- and V-step weighted), impute_prot (per-protein column-mean N) and impute_null_prot (protein labels of
N permuted within the dataset). The one revision: impute_v2 (= impute_prot) / impute_null_v2 (column means permuted
across proteins). Diagnostics + W = ones unit check:

    python -m scp_graphfs_mcp.phasec.pepev diag      # -> results/phaseC/pepev_diag.json

px063590_ADbrain has no N (PSM tables 2.16 GB, over the cap): every impute_* returns ctx['champion'] there.
See docs/phaseC/ideas/pepev.md.
"""
from __future__ import annotations

import hashlib
import io
import os
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "pepev"
RAW = ROOT / "data" / "cache" / "corpus_raw"
sys.path.insert(0, str(ROOT / "benchmarks"))
sys.path.insert(0, str(ROOT / "scripts" / "corpus_ext"))

ALIGNED = ("sadiku2025", "weiss2026_liver", "petrosius2023", "gregoire2023_mixCTRL", "leduc2022")
W_CLIP = (0.25, 4.0)


def _sha(X):
    a = np.ascontiguousarray(np.asarray(X, np.float32))
    return hashlib.sha256(a.tobytes()).hexdigest()


def _replay_save(X, accs, N):
    """build_corpus.save() selection logic, carrying N through the same columns / rows."""
    X = np.asarray(X, float)
    N = np.asarray(N, float)
    keep_cols = [i for i, a in enumerate(accs) if a]
    X, N, accs = X[:, keep_cols], N[:, keep_cols], [accs[i] for i in keep_cols]
    order = np.argsort(-(~np.isnan(X)).sum(0))
    seen, cols = set(), []
    for i in order:
        if accs[i] not in seen:
            seen.add(accs[i])
            cols.append(i)
    cols.sort()
    X, N, accs = X[:, cols], N[:, cols], [accs[i] for i in cols]
    rows = (~np.isnan(X)).sum(1) >= 50
    X, N = X[rows], N[rows]
    X = X - np.nanmedian(X, axis=1, keepdims=True)
    return X, accs, N


def _finish(target, Xr, accs, N, note):
    from explore_phaseb import load_target

    X, names, _ = load_target(target)
    assert list(accs) == list(names), f"{target}: protein order differs"
    assert Xr.shape == X.shape, (Xr.shape, X.shape)
    # the corpus stores float32; compare at float32 precision
    eq = np.allclose(np.asarray(Xr, np.float32), X, equal_nan=True, atol=1e-5, rtol=0)
    assert eq, f"{target}: X_rebuilt != X_corpus"
    obs = ~np.isnan(X)
    N = np.where(obs, N, np.nan)
    raw_pos = np.isfinite(N) & (N > 0)
    frac_obs_with_ev = float(raw_pos[obs].mean())
    N = np.where(obs & ~raw_pos, 1.0, N)
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(CACHE / f"{target}_N.npz", N=N.astype(np.float32), sha=_sha(X), names=np.array(names),
                        frac_obs_with_evidence=frac_obs_with_ev, note=note)
    print(f"{target}: aligned {X.shape}, X equality OK, observed entries with N>0 in source = {frac_obs_with_ev:.4f}")


# --- (a) sadiku2025: Spectronaut pivot, per-run PG.NrOfStrippedSequencesUsedForQuantification -----------------
def build_sadiku2025():
    import pandas as pd
    from scp_graphfs_mcp.ids import group_accessions

    rep = RAW / "sadiku2025" / "20250204_112949_gbm_sc_full_fasta_19.4_Report.tsv"
    d = pd.read_csv(rep, sep="\t", low_memory=False)
    design = pd.read_csv(RAW / "sadiku2025" / "experimentalDesignTemplate.txt", sep="\t")
    patient = dict(zip(design["raw file"], design["patient"]))
    cols, ncols = [], []
    for c in d.columns:
        if not c.endswith(".raw.PG.Quantity"):
            continue
        run = re.sub(r"^\[\d+\] ", "", c)[: -len(".PG.Quantity")]
        if run in patient and "0cells" not in run:
            cols.append(c)
            nc = c[: -len("PG.Quantity")] + "PG.NrOfStrippedSequencesUsedForQuantification"
            assert nc in d.columns, nc
            ncols.append(nc)
    X = d[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float).T
    N = d[ncols].apply(pd.to_numeric, errors="coerce").to_numpy(float).T
    with np.errstate(divide="ignore", invalid="ignore"):
        X = np.where(X > 0, np.log2(X), np.nan)
    accs = [(group_accessions(str(g)) or [""])[0] for g in d["PG.ProteinGroups"]]
    seen_sc = ~np.isnan(X).all(0)
    X, N, accs = X[:, seen_sc], N[:, seen_sc], [a for a, k in zip(accs, seen_sc) if k]
    Xr, accs, N = _replay_save(X, accs, N)
    _finish("sadiku2025", Xr, accs, N, "Spectronaut PG.NrOfStrippedSequencesUsedForQuantification per run")


# --- (b) weiss2026_liver: DIA-NN report.pr_matrix.tsv, N = # non-NaN precursor rows per (run, Protein.Group) ---
def build_weiss2026_liver():
    import pandas as pd
    from scp_graphfs_mcp.ids import group_accessions

    d = pd.read_csv(RAW / "weiss2026" / "report.pg_matrix.tsv", sep="\t")
    pat = re.compile(r"_box\d+_(.+?)_shape_\d+\.mzML$", re.I)
    cols = [c for c in d.columns if pat.search(c)]
    assert len(cols) == 792, len(cols)
    pr = pd.read_csv(CACHE / "dl_weiss" / "report.pr_matrix.tsv", sep="\t", usecols=["Protein.Group"] + cols)
    cnt = pr[cols].notna().groupby(pr["Protein.Group"]).sum()
    miss = set(d["Protein.Group"]) - set(cnt.index)
    print(f"weiss2026_liver: {len(miss)} pg_matrix protein groups without precursor rows")
    N = cnt.reindex(d["Protein.Group"]).to_numpy(float).T  # runs x PG, NaN where no precursor row
    X = d[cols].to_numpy(float).T
    with np.errstate(divide="ignore", invalid="ignore"):
        X = np.where(X > 0, np.log2(X), np.nan)
    accs = [(group_accessions(str(g)) or [""])[0] for g in d["Protein.Group"]]
    seen_sc = ~np.isnan(X).all(0)
    X, N, accs = X[:, seen_sc], N[:, seen_sc], [a for a, k in zip(accs, seen_sc) if k]
    Xr, accs, N = _replay_save(X, accs, N)
    _finish("weiss2026_liver", Xr, accs, N, "DIA-NN report.pr_matrix.tsv: non-NaN precursors per run and Protein.Group")


# --- (c) petrosius2023: Spectronaut PEPQuant (Normal) export, N = nunique(EG.PrecursorId) at EG.Qvalue <= 0.01 -------
def build_petrosius2023():
    import pandas as pd
    from scp_graphfs_mcp.datasets import load_petrosius2023

    a = load_petrosius2023()
    runs = [r[: -len(".raw")] if r.endswith(".raw") else r for r in a.obs["run"].astype(str)]
    names = list(a.var_names)
    ri = {r: i for i, r in enumerate(runs)}
    ci = {c: j for j, c in enumerate(names)}
    parts = []
    for ch in pd.read_csv(CACHE / "dl_petro" / "orig.tsv", sep="\t", chunksize=2_000_000,
                          usecols=["R.FileName", "PG.ProteinAccessions", "EG.PrecursorId", "EG.Qvalue", "PEP.Quantity"]):
        ch = ch[(pd.to_numeric(ch["EG.Qvalue"], errors="coerce") <= 0.01)
                & (pd.to_numeric(ch["PEP.Quantity"], errors="coerce") > 0)]
        parts.append(ch[["R.FileName", "PG.ProteinAccessions", "EG.PrecursorId"]].drop_duplicates())
    p = pd.concat(parts).drop_duplicates()
    n = p.groupby(["R.FileName", "PG.ProteinAccessions"])["EG.PrecursorId"].nunique().reset_index()
    print(f"petrosius2023: runs in table {n['R.FileName'].nunique()}, in h5ad {len(runs)}, matched "
          f"{len(set(runs) & set(n['R.FileName']))}; PGs matched {len(set(names) & set(n['PG.ProteinAccessions']))}"
          f"/{len(names)}")
    N = np.zeros((len(runs), len(names)))
    r = n["R.FileName"].map(ri)
    c = n["PG.ProteinAccessions"].map(ci)
    k = r.notna() & c.notna()
    N[r[k].astype(int).values, c[k].astype(int).values] = n["EG.PrecursorId"][k].values
    X = np.asarray(a.X, float)  # load_target uses exactly this matrix; no corpus rebuild applies
    _finish("petrosius2023", X, names, N, "Spectronaut PEPQuant: nunique EG.PrecursorId, EG.Qvalue<=0.01, PEP.Quantity>0")


# --- (d) gregoire2023_mixCTRL: scpdata peptides_<batch>_norm_log (parent of proteins_<batch>_norm_log) --------------
def _sce_rows(sce, m):
    names = sce.attr.get("NAMES")
    if names is not None and names.type == 16 and names.value is not None:
        return names.value
    if m.attr.get("dimnames") is not None and m.attr["dimnames"].value[0] is not None:
        return m.attr["dimnames"].value[0].value
    return sce.attr["rowRanges"].attr["partitioning"].attr["NAMES"].value


def build_gregoire2023_mixCTRL():
    import pandas as pd
    from rda_stream import dataframe, load_rda, qf_assays, sce_matrix
    from scp_graphfs_mcp.ids import group_accessions

    batches = ["CBIO680", "CBIO681", "CBIO703", "CBIO715", "CBIO725", "CBIO754", "GIGA"]
    qf = load_rda(RAW / "gregoire2023_mixCTRL" / "gregoire2023_mixCTRL.Rda", max_str=10_000_000)[0][1]
    assays = qf_assays(qf)
    _, sm = dataframe(qf.attr["sampleMap"])
    cd_rn, cd = dataframe(qf.attr["colData"])
    ok = {r for r, t, n in zip(cd_rn, cd["sample_type"], np.array(cd["cell_number"]).astype(str)) if t == "SC" and n == "1"}
    fx, fn = [], []
    for b in batches:
        out = {}
        for kind in ("proteins", "peptides"):
            a = f"{kind}_{b}_norm_log"
            sce = assays[a]
            m, dim = sce_matrix(sce)
            M = np.asarray(m.value, dtype=float).reshape(dim[1], dim[0]).T
            rows = _sce_rows(sce, m)
            cols = sce.attr["colData"].attr["rownames"].value
            k = np.array(sm["assay"]) == a
            c2p = dict(zip(np.array(sm["colname"])[k], np.array(sm["primary"])[k]))
            prim = np.array([c2p.get(c, "") for c in cols])
            prot = None
            if kind == "peptides":
                _, rd = dataframe(sce.attr["rowRanges"].attr["elementMetadata"])
                prot = np.array(rd["proteins"])
            out[kind] = (M, np.array(rows), prim, prot)
        Xp, prow, pprim, _ = out["proteins"]
        Mpep, _, eprim, eprot = out["peptides"]
        cnt = pd.DataFrame((~np.isnan(Mpep)).astype(float), columns=eprim).groupby(eprot).sum()
        miss = len(set(prow) - set(cnt.index))
        assert list(pprim) == list(eprim), f"{b}: peptide / protein assay columns differ"
        Nb = cnt.reindex(prow).to_numpy(float)
        keep = np.array([p in ok for p in pprim])
        fx.append(pd.DataFrame(Xp[:, keep].T, columns=prow, index=pprim[keep]))
        fn.append(pd.DataFrame(Nb[:, keep].T, columns=prow, index=pprim[keep]))
        print(f"  {b}: {keep.sum()} cells, {len(prow)} protein groups ({miss} without peptides in the peptide assay)")
    P = pd.concat(fx, axis=0, join="outer")
    Q = pd.concat(fn, axis=0, join="outer").reindex(index=P.index, columns=P.columns)
    accs = [(group_accessions(str(r)) or [""])[0] for r in P.columns]
    Xr, accs, N = _replay_save(P.values, accs, Q.values)
    _finish("gregoire2023_mixCTRL", Xr, accs, N, "scpdata peptides_<batch>_norm_log: non-NaN peptides per protein group")


# --- (e) leduc2022: scpdata peptides_log (parent of proteins_norm2), N = non-NaN peptides per Leading.razor.protein --
def build_leduc2022():
    import pandas as pd
    from rda_stream import dataframe, load_rda, qf_assays, sce_matrix
    from explore_phaseb import load_target

    ex = CACHE / "leduc2022_extract.npz"
    if not ex.exists():
        qf = load_rda(CACHE / "dl_leduc" / "leduc2022_pSCoPE.Rda", max_str=10_000_000)[0][1]
        assays = qf_assays(qf)
        res = {}
        for a in ("proteins_norm2", "peptides_log"):
            sce = assays[a]
            m, dim = sce_matrix(sce)
            res[f"M_{a}"] = np.asarray(m.value, dtype=float).reshape(dim[1], dim[0]).T
            res[f"rows_{a}"] = np.array(list(_sce_rows(sce, m)))
            res[f"cols_{a}"] = np.array(list(sce.attr["colData"].attr["rownames"].value))
            _, rd = dataframe(sce.attr["rowRanges"].attr["elementMetadata"])
            for k in ("Leading.razor.protein", "Leading.razor.protein.id", "Leading.razor.protein.symbol"):
                res[f"{k}_{a}"] = np.array(rd[k]).astype(str)
        np.savez_compressed(ex, **res)
    d = np.load(ex)
    out = {a: (d[f"M_{a}"], list(d[f"rows_{a}"]), list(d[f"cols_{a}"]), d[f"Leading.razor.protein.id_{a}"])
           for a in ("proteins_norm2", "peptides_log")}
    Xp, prow, pcol, prz = out["proteins_norm2"]
    Mpep, _, ecol, erz = out["peptides_log"]
    assert pcol == ecol, "peptides_log / proteins_norm2 columns differ"
    cnt = pd.DataFrame((~np.isnan(Mpep)).astype(float), columns=ecol).groupby(erz).sum()
    # proteins_norm2 rows are named by Leading.razor.protein.id
    assert list(prz) == list(prow), "proteins_norm2 row names != Leading.razor.protein.id"
    X, names, _ = load_target("leduc2022")
    from scp_graphfs_mcp.datasets import load_leduc2022

    cells = list(load_leduc2022().obs_names)
    Xr = pd.DataFrame(Xp.T, index=pcol, columns=prow).loc[cells, names].to_numpy(float)
    N = cnt.T.reindex(index=cells, columns=names).to_numpy(float)
    _finish("leduc2022", Xr, names, N, "scpdata peptides_log: non-NaN peptides per Leading.razor.protein.id")


BUILDERS = {"sadiku2025": build_sadiku2025, "weiss2026_liver": build_weiss2026_liver,
            "petrosius2023": build_petrosius2023, "gregoire2023_mixCTRL": build_gregoire2023_mixCTRL,
            "leduc2022": build_leduc2022}


# ======================================================================================================================
# Model: heteroscedastic U-step of the domain-adapted projection
# ======================================================================================================================
_TGT = {}


def _target(target):
    """(X, names, labels, N or None), cached per process; N is checked against the sha of the current X."""
    if target not in _TGT:
        from explore_phaseb import load_target

        X, names, labels = load_target(target)
        f = CACHE / f"{target}_N.npz"
        N = None
        if f.exists():
            d = np.load(f)
            assert str(d["sha"]) == _sha(X), f"{target}: N cache is stale (X sha differs)"
            assert list(d["names"]) == list(names), f"{target}: N cache protein order differs"
            N = d["N"].astype(float)
        _TGT[target] = (X, names, labels, N)
    return _TGT[target]


def _aligned_counts(Xm, ctx):
    """N for the 30 subsampled cells of this split, after asserting Xm is a masked copy of X[idx]; None if the
    target has no N (non-aligned)."""
    from explore_phaseb import subsample

    X, names, labels, N = _target(ctx["target"])
    if N is None:
        print(f"[pepev] {ctx['target']}: not aligned, returning the champion unchanged", file=sys.stderr)
        return None, None
    assert list(names) == list(ctx["names"])
    idx = subsample(X.shape[0], labels, 30, ctx["seed"])
    Xs = X[idx]
    om, os_ = ~np.isnan(Xm), ~np.isnan(Xs)
    assert Xm.shape == Xs.shape and not (om & ~os_).any(), "obs(Xm) is not a subset of obs(X[idx])"
    assert np.allclose(Xm[om], Xs[om]), "Xm values differ from X[idx]"
    assert om.sum() >= 0.7 * os_.sum(), "Xm keeps < 70% of X[idx]'s observed entries"
    return N[idx], idx


def _weights(Xm, Nsub, alpha=0.5, mode="entry", seed=0, Nfull=None):
    """Per-cell normalised precision weights; used only at o_i = obs[i] & ok (ok = column observed in Xm).
    mode: entry (real), perm_cell (null: permuted among each cell's o_i), protein (column mean of N over the
    dataset's cells outside this split -- Nfull is passed with the split's rows removed -- constant across cells)."""
    obs = ~np.isnan(Xm)
    ok = np.isfinite(np.nanmean(Xm, 0))
    if mode == "protein":
        cm = np.nanmean(Nfull, 0)
        cm = np.where(np.isfinite(cm) & (cm > 0), cm, 1.0)
        P = np.broadcast_to(cm ** alpha, Xm.shape).copy()
    else:
        P = np.where(np.isfinite(Nsub) & (Nsub > 0), Nsub, 1.0) ** alpha
    W = np.ones_like(P)
    for i in range(Xm.shape[0]):
        o = obs[i] & ok
        if not o.any():
            continue
        p = P[i, o]
        w = p / p.mean()
        w = np.clip(w, *W_CLIP)
        w = w / w.mean()
        if mode == "perm_cell":
            w = np.random.default_rng(seed + 1000 * i).permutation(w)
        W[i, o] = w
    return W


def adapt_project_w(Xm, Vref, usable, base, W, lam_v=100.0, lam=5.0, extra=10, iters=10, seed=0, weight_v=False):
    """transfer.adapt_project with a weighted U-step (and, if weight_v, a weighted V-step). W = ones reproduces
    transfer.adapt_project bit-exactly."""
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
            Vo = V[o]
            Vw = Vo * W[i, o][:, None]
            U[i] = np.linalg.solve(Vw.T @ Vo + lam * I, Vw.T @ Z[i, o])
        for j in np.where(ok)[0]:
            o = obs[:, j]
            Uo = U[o]
            if weight_v:
                Uw = Uo * W[o, j][:, None]
                V[j] = np.linalg.solve(Uw.T @ Uo + lam_v * I, Uw.T @ Z[o, j] + lam_v * Vp[j])
            else:
                V[j] = np.linalg.solve(Uo.T @ Uo + lam_v * I, Uo.T @ Z[o, j] + lam_v * Vp[j])
    out = base.copy()
    P = U @ V.T
    miss = ~obs & ok[None, :]
    out[miss] = (mu[None, :] + sd[None, :] * P)[miss]
    return out


def _run(Xm, ctx, alpha=0.5, mode="entry", weight_v=False, permute_proteins=False, permute_colmeans=False):
    Nsub, idx = _aligned_counts(Xm, ctx)
    if Nsub is None:
        return ctx["champion"]
    Nfull = _target(ctx["target"])[3]
    Nrest = np.delete(Nfull, idx, axis=0)
    if permute_proteins:  # weaker null: each column takes another protein's count profile (dataset-wide permutation)
        perm = np.random.default_rng(12345).permutation(Nfull.shape[1])
        Np = Nfull[idx][:, perm]
        fill = np.round(np.nanmean(Nrest[:, perm], 0))  # source entry unobserved -> its column's mean outside the split
        Nsub = np.where(np.isfinite(Np), Np, fill[None, :])
    if permute_colmeans:  # v2 null: protein-level evidence profiles reassigned to random other proteins
        Nrest = Nrest[:, np.random.default_rng(ctx["seed"] + 7).permutation(Nrest.shape[1])]
    W = _weights(Xm, Nsub, alpha=alpha, mode=mode, seed=ctx["seed"], Nfull=Nrest)
    return adapt_project_w(Xm, ctx["V"], ctx["warm"], ctx["base"], W, lam_v=100.0, extra=10, weight_v=weight_v)


def impute(Xm, ctx):
    return _run(Xm, ctx)


def impute_null(Xm, ctx):
    return _run(Xm, ctx, mode="perm_cell")


def impute_a1(Xm, ctx):
    return _run(Xm, ctx, alpha=1.0)


def impute_null_a1(Xm, ctx):
    return _run(Xm, ctx, alpha=1.0, mode="perm_cell")


def impute_uv(Xm, ctx):
    return _run(Xm, ctx, weight_v=True)


def impute_null_uv(Xm, ctx):
    return _run(Xm, ctx, mode="perm_cell", weight_v=True)


def impute_prot(Xm, ctx):
    return _run(Xm, ctx, mode="protein")


def impute_null_prot(Xm, ctx):
    return _run(Xm, ctx, permute_proteins=True)


# v2 (the one revision): protein-level evidence weights (column-mean N outside the split), since the per-entry part of
# N tracks the entry's own value (informative weights -> biased U under MNAR). Null: column means permuted.
impute_v2 = impute_prot


def impute_null_v2(Xm, ctx):
    return _run(Xm, ctx, mode="protein", permute_colmeans=True)


def diagnostics(out=ROOT / "results" / "phaseC" / "pepev_diag.json"):
    """N distribution, Spearman(N, value) pooled / within cell / within protein, and the W = ones unit check
    (one MNAR split, seed 100, per aligned target; same split construction as benchmarks/phasec.py)."""
    import json

    from scipy.stats import spearmanr

    from explore_phaseb import subsample
    from run import make_mask
    from scp_graphfs_mcp import phase2 as P2
    from scp_graphfs_mcp import transfer as T
    import phasec as H

    res = {}
    for t in ALIGNED:
        X, names, labels, N = _target(t)
        obs = ~np.isnan(X)
        n, x = N[obs], X[obs]
        d = dict(np.load(CACHE / f"{t}_N.npz"))
        r = {"cells": int(X.shape[0]), "proteins": int(X.shape[1]), "observed_entries": int(obs.sum()),
             "frac_observed_with_source_evidence": float(d["frac_obs_with_evidence"]), "source": str(d["note"]),
             "N_quantiles_p10_p25_p50_p75_p90": np.percentile(n, [10, 25, 50, 75, 90]).tolist(),
             "N_mean": float(n.mean()), "frac_N_eq_1": float((n == 1).mean()),
             "spearman_N_value_pooled": float(spearmanr(n, x)[0])}
        rc = [spearmanr(N[i, obs[i]], X[i, obs[i]])[0] for i in range(X.shape[0]) if obs[i].sum() > 10]
        cc = [spearmanr(N[obs[:, j], j], X[obs[:, j], j])[0] for j in range(X.shape[1]) if obs[:, j].sum() > 10]
        cc = [c for c in cc if np.isfinite(c)]
        r["spearman_within_cell_median"] = float(np.nanmedian(rc))
        r["spearman_within_protein_median"] = float(np.median(cc))
        cm = np.nanmean(N, 0)
        r["frac_var_logN_between_proteins"] = float(np.nanvar(np.log(np.broadcast_to(cm, N.shape)[obs]))
                                                    / np.var(np.log(n)))
        V, warm = H.target_models(t, names)
        seed = 100
        idx = subsample(X.shape[0], labels, 30, seed)
        Xs = X[idx]
        elig = ((~np.isnan(Xs)).mean(0) >= 0.5) & ~np.isnan(Xs).all(0)
        mask = make_mask(Xs, elig, "MNAR", seed)
        Xm = Xs.copy()
        Xm[mask] = np.nan
        b = P2._fill(Xm, P2.loo_knn(Xm))
        b = np.where(np.isnan(b), np.nanmean(Xm), b)
        ad = T.adapt_project(Xm, V, warm, b, lam_v=100, extra=10)
        ad1 = adapt_project_w(Xm, V, warm, b, np.ones_like(Xm), lam_v=100, extra=10)
        r["unit_check_W_ones_MNAR_seed100"] = {"array_equal": bool(np.array_equal(ad, ad1)),
                                               "max_abs_diff": float(np.max(np.abs(ad - ad1)))}
        Ns, _ = _aligned_counts(Xm, {"target": t, "seed": seed, "names": names})
        W = _weights(Xm, Ns)
        o = ~np.isnan(Xm) & np.isfinite(np.nanmean(Xm, 0))[None, :]
        r["W_alpha05_split_MNAR_seed100"] = {"min": float(W[o].min()), "max": float(W[o].max()),
                                             "p10_p50_p90": np.percentile(W[o], [10, 50, 90]).tolist()}
        # MNAR mask removes low values: are the masked entries' counts lower than the kept entries' counts?
        Nss = N[idx]
        r["N_median_masked_vs_kept_MNAR_seed100"] = [float(np.median(Nss[mask])), float(np.median(Nss[o]))]
        res[t] = r
        print(t, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k != "source"})
    res["px063590_ADbrain"] = "not aligned (size): PSM tables 1.22 GB + 0.94 GB = 2.16 GB exceed the 2 GB cap"
    Path(out).write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "build":
        for t in sys.argv[2:]:
            BUILDERS[t]()
    elif sys.argv[1] == "diag":
        diagnostics()
