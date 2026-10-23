"""Phase-C idea `scrnaprior` (queen: new data modalities).

HPA single-cell-type RNA consensus (gene x 154 cell types, nCPM) -> genome-wide log1p, drop max nCPM < 1, z-score each
gene across cell types, truncated SVD, keep k=8 gene loadings L = A_k Sigma_k / sqrt(n_celltypes). Map each target
protein to genes (proteinatlas.tsv Ensembl->Uniprot, fallback HGNC ensembl_gene_id->uniprot_ids; mean over several
genes), giving S (p x 8), zero for unmapped proteins, globally rescaled to the RMS of the warm proteomic loadings.
S enters a local copy of transfer.adapt_project as an extra prior-anchored block between the 20 proteomic and the
10 free factors, with its own ridge: D = diag([100]*20 + [lam_s]*8 + [100]*10), lam_s = 30 (fixed a priori).
Null: identical, but the rows of S are permuted among mapped proteins (rng seed + 7).
Data: https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip (CC BY-SA 4.0), cached in
data/cache/phasec/scrnaprior/; genome-wide SVD cached as svd_k8.npz; per-names S cached as S_<hash>.npz.
"""

from __future__ import annotations

import hashlib
import os
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "scrnaprior"
RAW = CACHE / "rna_single_cell_type.tsv.zip"
HPA_FEAT = ROOT / "data" / "cache" / "features" / "proteinatlas.tsv.zip"
K, LAM_S, LAM_V, LAM, EXTRA, ITERS = 8, 30.0, 100.0, 5.0, 10, 10
LAM_S = float(os.environ.get("SCRNAPRIOR_LAM_S", LAM_S))  # diagnostic only; default fixed a priori

_MEM: dict = {}


def _svd_cache():
    f = CACHE / f"svd_k{K}.npz"
    if f.exists():
        d = np.load(f, allow_pickle=False)
        return d["genes"], d["L"]
    import pandas as pd

    with zipfile.ZipFile(RAW) as z:
        df = pd.read_csv(z.open(z.namelist()[0]), sep="\t", usecols=["Gene", "Cell type", "nCPM"])
    M = df.pivot_table(index="Gene", columns="Cell type", values="nCPM", aggfunc="mean")
    M = M.fillna(0.0)
    np.savez_compressed(CACHE / "pivot.npz", genes=np.array(M.index.tolist(), dtype="U"), celltypes=np.array(M.columns.tolist(), dtype="U"),
                        nCPM=M.values.astype(np.float32))
    keep = M.values.max(1) >= 1.0
    R = np.log1p(M.values[keep])
    Rz = (R - R.mean(1, keepdims=True)) / np.maximum(R.std(1, keepdims=True), 1e-9)
    A, s, _ = np.linalg.svd(Rz, full_matrices=False)
    L = A[:, :K] * s[:K] / np.sqrt(Rz.shape[1])
    genes = np.array(M.index.tolist(), dtype="U")[keep]
    np.savez_compressed(f, genes=genes, L=L, sv=s, n_celltypes=Rz.shape[1])
    return genes, L


def _acc2genes():
    if "acc" in _MEM:
        return _MEM["acc"]
    import pandas as pd

    acc: dict[str, set[str]] = {}
    with zipfile.ZipFile(HPA_FEAT) as z:
        d = pd.read_csv(z.open(z.namelist()[0]), sep="\t", usecols=["Ensembl", "Uniprot"], dtype=str).dropna()
    for g, u in zip(d["Ensembl"], d["Uniprot"]):
        for a in str(u).replace(";", ",").split(","):
            a = a.strip()
            if a:
                acc.setdefault(a, set()).add(g)
    from scp_graphfs_mcp import ids

    h = ids.hgnc_table()
    fb: dict[str, set[str]] = {}
    for g, u in zip(h["ensembl_gene_id"], h["uniprot_ids"]):
        if not g:
            continue
        for a in str(u).replace(";", "|").replace(",", "|").split("|"):
            a = a.strip()
            if a:
                fb.setdefault(a, set()).add(g)
    _MEM["acc"] = (acc, fb)
    return acc, fb


def build_S(names, Vw):
    """S (p x K) scaled to the warm proteomic loadings, and smask (mapped proteins). Cached per names list."""
    key = hashlib.sha1("\n".join(map(str, names)).encode()).hexdigest()[:16]
    f = CACHE / f"S_{key}.npz"
    if f.exists():
        d = np.load(f)
        S, smask = d["S"], d["smask"]
    else:
        from scp_graphfs_mcp.ids import group_accessions

        genes, L = _svd_cache()
        gi = {g: i for i, g in enumerate(genes)}
        acc, fb = _acc2genes()
        S = np.zeros((len(names), K))
        smask = np.zeros(len(names), bool)
        for j, nm in enumerate(names):
            accs = group_accessions(str(nm))
            hit = {g for a in accs for g in acc.get(a, ()) if g in gi}
            if not hit:
                hit = {g for a in accs for g in fb.get(a, ()) if g in gi}
            if hit:
                S[j] = L[[gi[g] for g in sorted(hit)]].mean(0)
                smask[j] = True
        np.savez_compressed(f, S=S, smask=smask)
    nz = np.abs(Vw).sum(1) > 0
    rms_v = np.sqrt(np.mean(Vw[nz] ** 2))
    nzs = np.abs(S).sum(1) > 0
    rms_s = np.sqrt(np.mean(S[nzs] ** 2)) if nzs.any() else 1.0
    return S * (rms_v / rms_s), smask


def _adapt(Xm, Vref, usable, S, base, lam_s=LAM_S, lam_v=LAM_V, lam=LAM, extra=EXTRA, iters=ITERS, seed=0):
    """Copy of transfer.adapt_project with an extra S block and a diagonal per-block ridge in the V-step."""
    obs = ~np.isnan(Xm)
    mu = np.nanmean(Xm, 0)
    sd = np.nanstd(Xm, 0)
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    ok = np.isfinite(mu)
    Z = np.where(obs, (Xm - np.nan_to_num(mu)) / sd, 0.0)
    r0, k = Vref.shape[1], S.shape[1]
    rng = np.random.default_rng(seed)
    Vp = np.hstack([np.where(usable[:, None], Vref, 0.0), S, np.zeros((len(Vref), extra))])
    V = Vp + np.hstack([np.zeros((len(Vref), r0 + k)), 0.01 * rng.normal(size=(len(Vref), extra))])
    r = r0 + k + extra
    U = np.zeros((Xm.shape[0], r))
    I = np.eye(r)
    dvec = np.array([lam_v] * r0 + [lam_s] * k + [lam_v] * extra)
    D = np.diag(dvec)
    for _ in range(iters):
        for i in range(Xm.shape[0]):
            o = obs[i] & ok
            U[i] = np.linalg.solve(V[o].T @ V[o] + lam * I, V[o].T @ Z[i, o])
        for j in np.where(ok)[0]:
            o = obs[:, j]
            Uo = U[o]
            V[j] = np.linalg.solve(Uo.T @ Uo + D, Uo.T @ Z[o, j] + dvec * Vp[j])
    out = base.copy()
    P = U @ V.T
    miss = ~obs & ok[None, :]
    out[miss] = (mu[None, :] + sd[None, :] * P)[miss]
    return out


def _S_for(ctx):
    return build_S(ctx["names"], ctx["V"][ctx["warm"]])


def impute(Xm, ctx):
    S, _ = _S_for(ctx)
    return _adapt(Xm, ctx["V"], ctx["warm"], S, ctx["base"], seed=ctx["seed"])


def impute_null(Xm, ctx):
    S, smask = _S_for(ctx)
    idx = np.where(smask)[0]
    perm = np.random.default_rng(ctx["seed"] + 7).permutation(idx)
    Sn = S.copy()
    Sn[idx] = S[perm]
    return _adapt(Xm, ctx["V"], ctx["warm"], Sn, ctx["base"], seed=ctx["seed"])


def impute_zero(Xm, ctx):
    """Secondary sanity check: S = 0 (8 extra ridge-to-zero factors at lam_s)."""
    S, _ = _S_for(ctx)
    return _adapt(Xm, ctx["V"], ctx["warm"], np.zeros_like(S), ctx["base"], seed=ctx["seed"])
