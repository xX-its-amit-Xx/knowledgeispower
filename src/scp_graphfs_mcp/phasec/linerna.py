"""Phase-C idea `linerna` (queen: new data modalities). Round 5, credited against the cshrink base.

Within-cell-line single-cell RNA programmes as a ridge-weighted second-stage factor block on the cshrink residual.

A. Per line (U2OS: GEO GSE146773 FUCCI U-2 OS Smart-seq2; U937 / PC9: GEO GSM3972674 Drop-seq PC9 + U-937 mix, split
   by marker sums with the 18 marker genes removed before PCA): x = log1p(1e4 count / total); genes detected in >= 5%
   of the line's cells; z-score per gene; truncated SVD k = 20; Lg = v_k s_k / sqrt(n - 1), hg = clip(1 - sum Lg^2,
   0.05, 1). Genes -> target proteins (HPA Ensembl/Gene -> UniProt, HGNC fallback); a protein row is the mean of its
   genes' rows (Lp, h); unmapped proteins: Lp = 0, h = 1, mapped = False.
B. _fit = a copy of nfbase.adapt_project_mb that also returns Z (standardised input) and P = U V^T.
C. impute: psi selected exactly as cshrink (shared score cache, tag 'c'); on the same 50/50 inner mask, residual
   E = Z - P on observed entries, scaled per column by its observed sd c_j; per cell the posterior of 20 latent
   programme scores f_i = (Lp_o' H_o^-1 Lp_o + I)^-1 Lp_o' H_o^-1 E_io (|o| >= 20), ehat_m = Lp_m f_i;
   Y = out + gamma sd c ehat at missing mapped entries, gamma in {0, .25, .5, 1, 2} on the inner mask (ties -> smaller);
   refit on the input with iters = 10.
Entry points: impute (matched line), impute_null (PC9 for every target), impute_perm (matched line, Lp/h rows
permuted among mapped proteins, rng seed + 7), impute_base (gamma = 0; equals cshrink.impute).
v2 (impute_v2 / impute_null_v2 / impute_perm_v2): the residual column scale c_j uses the observed sd only for
columns with >= 5 observed residuals (else 1), floored at 0.1; v1's 1e-3 floor let 1-2-observation columns explode.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import sys
import zipfile
from pathlib import Path

import numpy as np

from scp_graphfs_mcp.phasec import cshrink as CS
from scp_graphfs_mcp.phasec import tgtshrink as TS
from scp_graphfs_mcp.phasec.nfbase import EXTRA, LAM_V

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "linerna"
HPA_FEAT = ROOT / "data" / "cache" / "features" / "proteinatlas.tsv.zip"
U2OS_RAW = CACHE / "GSE146773_Counts.csv.gz"
DS_RAW = CACHE / "GSM3972674_PC9U937.dge.txt.gz"
K = 20
GAMMAS = (0.0, 0.25, 0.5, 1.0, 2.0)
MYELOID = ("LCP1", "CORO1A", "TYROBP", "SPI1", "CD33", "LYZ", "SRGN", "CD52", "ARHGDIB", "LAPTM5")
EPITHELIAL = ("EPCAM", "KRT8", "KRT18", "KRT19", "CLDN4", "S100A6", "TACSTD2", "CDH1")
LINE = {"krull2024": "U2OS", "specht2021": "U937", "gregoire2023_mixCTRL": "U937", "leduc2022": "U937"}
NULL_LINE = "PC9"
_MEM: dict = {}


# ---------------------------------------------------------------------------------------------------------------
# A. RNA factor loadings per line
def _pca(C, genes):
    """C: cells x genes raw counts. Returns dict with genes, L (g x K), h, var_explained."""
    C = np.asarray(C, float)
    tot = C.sum(1, keepdims=True)
    X = np.log1p(1e4 * C / np.maximum(tot, 1e-9))
    det = (C > 0).mean(0) >= 0.05
    X, genes = X[:, det], np.asarray(genes)[det]
    sd = X.std(0, ddof=1)
    good = sd > 1e-9
    X, genes, sd = X[:, good], genes[good], sd[good]
    Zg = (X - X.mean(0)) / sd
    n = Zg.shape[0]
    _, s, vt = np.linalg.svd(Zg, full_matrices=False)
    L = vt[:K].T * s[:K] / np.sqrt(n - 1)
    h = np.clip(1.0 - (L ** 2).sum(1), 0.05, 1.0)
    ve = (s ** 2) / np.sum(s ** 2)
    return {"genes": genes.astype("U"), "L": L, "h": h, "var_explained": ve[:K], "n_cells": n, "n_genes": len(genes)}


def _save(line, d):
    CACHE.mkdir(parents=True, exist_ok=True)
    tmp = CACHE / f".{line}_pcs.tmp{os.getpid()}.npz"
    np.savez_compressed(tmp, **d)
    os.replace(tmp, CACHE / f"{line}_pcs.npz")


def _build_u2os():
    import pandas as pd

    df = pd.read_csv(U2OS_RAW, index_col=0)
    cols = [c for c in df.columns if str(c).startswith("ENSG")]  # drop ERCC spike-ins / non-gene columns
    d = _pca(df[cols].values, [str(c).split(".")[0] for c in cols])
    d["n_dropped_nongene_cols"] = df.shape[1] - len(cols)
    _save("U2OS", d)


def _build_dropseq():
    import pandas as pd

    df = pd.read_csv(DS_RAW, sep="\t", index_col=0)
    df = df.astype(np.int32)
    g = df.index.astype(str)
    Mk = df.loc[[x for x in MYELOID if x in g]].sum(0).values
    Ek = df.loc[[x for x in EPITHELIAL if x in g]].sum(0).values
    u937 = (Mk > Ek) & (Mk >= 2)
    pc9 = (Ek > Mk) & (Ek >= 2)
    keep = ~np.isin(g, MYELOID + EPITHELIAL)
    C = df.values[keep].T  # cells x genes
    genes = np.asarray(g[keep])
    for line, sel in (("U937", u937), ("PC9", pc9)):
        d = _pca(C[sel], genes)
        d.update({"n_u937": int(u937.sum()), "n_pc9": int(pc9.sum()), "n_dropped": int((~u937 & ~pc9).sum()),
                  "n_markers_removed": int((~keep).sum())})
        _save(line, d)


def line_pcs(line):
    if line in _MEM:
        return _MEM[line]
    f = CACHE / f"{line}_pcs.npz"
    if not f.exists():
        (_build_u2os if line == "U2OS" else _build_dropseq)()
    d = np.load(f, allow_pickle=False)
    _MEM[line] = {k: d[k] for k in d.files}
    return _MEM[line]


def _maps():
    """(ensembl map, symbol map): UniProt accession -> set of gene ids, HPA primary with HGNC fallback dicts."""
    if "maps" in _MEM:
        return _MEM["maps"]
    import pandas as pd

    from scp_graphfs_mcp import ids

    with zipfile.ZipFile(HPA_FEAT) as z:
        d = pd.read_csv(z.open(z.namelist()[0]), sep="\t", usecols=["Gene", "Ensembl", "Uniprot"], dtype=str)
    h = ids.hgnc_table()

    def build(keys, unis):
        m: dict[str, set[str]] = {}
        for g, u in zip(keys, unis):
            if not isinstance(g, str) or not g or not isinstance(u, str):
                continue
            for a in u.replace(";", "|").replace(",", "|").split("|"):
                a = a.strip()
                if a:
                    m.setdefault(a, set()).add(g)
        return m

    ens = (build(d["Ensembl"], d["Uniprot"]), build(h["ensembl_gene_id"], h["uniprot_ids"]))
    sym = (build(d["Gene"], d["Uniprot"]), build(h["symbol"], h["uniprot_ids"]))
    _MEM["maps"] = {"ens": ens, "sym": sym}
    return _MEM["maps"]


def protein_loadings(names, line):
    """Lp (p x K), h (p,), mapped (p,) for the target's proteins. Cached per sha1(names) + line."""
    key = hashlib.sha1("\n".join(map(str, names)).encode()).hexdigest()[:16]
    f = CACHE / f"Lp_{line}_{key}.npz"
    if f.exists():
        d = np.load(f)
        return d["Lp"], d["h"], d["mapped"]
    from scp_graphfs_mcp.ids import group_accessions

    P = line_pcs(line)
    gi = {g: i for i, g in enumerate(P["genes"])}
    prim, fb = _maps()["ens" if line == "U2OS" else "sym"]
    p = len(names)
    Lp, h, mapped = np.zeros((p, K)), np.ones(p), np.zeros(p, bool)
    for j, nm in enumerate(names):
        accs = group_accessions(str(nm))
        hit = {g for a in accs for g in prim.get(a, ()) if g in gi}
        if not hit:
            hit = {g for a in accs for g in fb.get(a, ()) if g in gi}
        if hit:
            rows = [gi[g] for g in sorted(hit)]
            Lp[j] = P["L"][rows].mean(0)
            h[j] = P["h"][rows].mean()
            mapped[j] = True
    CACHE.mkdir(parents=True, exist_ok=True)
    tmp = CACHE / f".Lp_{line}_{key}.tmp{os.getpid()}.npz"
    np.savez_compressed(tmp, Lp=Lp, h=h, mapped=mapped)
    os.replace(tmp, f)
    return Lp, h, mapped


# ---------------------------------------------------------------------------------------------------------------
# B. the cshrink projection with its internals exposed (copy of nfbase.adapt_project_mb)
def _fit(Xm, Vref, usable, base, mu, sd, lam_v=10.0, lam=5.0, extra=3, iters=10, seed=0):
    obs = ~np.isnan(Xm)
    ok = np.isfinite(np.nanmean(Xm, 0))
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
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
            U[i] = np.linalg.solve(V[o].T @ V[o] + lam * I, V[o].T @ Z[i, o])
        for j in np.where(ok)[0]:
            o = obs[:, j]
            Uo = U[o]
            V[j] = np.linalg.solve(Uo.T @ Uo + lam_v * I, Uo.T @ Z[o, j] + lam_v * Vp[j])
    out = base.copy()
    P = U @ V.T
    miss = ~obs & ok[None, :]
    out[miss] = (mu[None, :] + sd[None, :] * P)[miss]
    return out, Z, P, obs, ok, mu, sd


# ---------------------------------------------------------------------------------------------------------------
# C. second stage
def _ehat(Z, P, obs, ok, Lp, h, mapped, cfix=False):
    """Returns (ehat in scaled-residual units, c (column scale), mean |o|, n cells used).

    v1 (cfix=False): c_j = observed residual sd, floor 1e-3 (as specified; columns with 1-2 observed residuals get
    c_j ~ 1e-3 and blow up E / c). v2 (cfix=True): c_j = observed sd only if the column has >= 5 observed residuals,
    else 1; floor 0.1."""
    E = np.where(obs & ok[None, :], Z - P, np.nan)
    with np.errstate(all="ignore"):
        c = np.nanstd(E, 0)
    if cfix:
        nobs = np.isfinite(E).sum(0)
        c = np.where(np.isfinite(c) & (nobs >= 5), c, 1.0)
        c = np.maximum(c, 0.1)
    else:
        c = np.where(np.isfinite(c), np.maximum(c, 1e-3), 1.0)
    Es = np.nan_to_num(E / c)
    eh = np.zeros_like(Z)
    Ik = np.eye(Lp.shape[1])
    no, used = [], 0
    for i in range(Z.shape[0]):
        o = obs[i] & ok & mapped
        m = ~obs[i] & ok & mapped
        no.append(int(o.sum()))
        if o.sum() < 20 or not m.any():
            continue
        Lo = Lp[o] / h[o][:, None]
        A = Lp[o].T @ Lo + Ik
        fi = np.linalg.solve(A, Lo.T @ Es[i, o])
        eh[i, m] = Lp[m] @ fi
        used += 1
    return eh, c, float(np.mean(no)), used


def _log(variant, rec):
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / f"log_{variant}.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def _select_psi(M, ctx):
    Msel, p1, p2 = TS.selection_mask(M, ctx["seed"])
    scores = {p: CS._score(M, Msel, (p1, p2), ctx, p, p, CS.cshrink_stats, "c") for p in CS.GRID}
    psi = min(CS.GRID, key=lambda p: (scores[p], p))
    return psi, scores, Msel, p1, p2


def run(M, ctx, variant, line=None, perm=False, force_gamma=None, cfix=False):
    seed = ctx["seed"]
    fit_kw = dict(lam_v=LAM_V, extra=EXTRA)
    psi, pscores, Msel, p1, p2 = _select_psi(M, ctx)
    if line is None and ctx.get("target") not in LINE:
        mu, sd, _ = CS.cshrink_stats(M, psi, psi)
        return _fit(M, ctx["V"], ctx["warm"], ctx["base"], mu, sd, iters=10, **fit_kw)[0]
    line = line or LINE[ctx["target"]]
    Lp, h, mapped = protein_loadings(ctx["names"], line)
    if perm:
        idx = np.where(mapped)[0]
        pr = np.random.default_rng(seed + 7).permutation(idx)
        Lp, h = Lp.copy(), h.copy()
        Lp[idx], h[idx] = Lp[pr], h[pr]
    rec = {"target": ctx.get("target"), "seed": seed, "n_missing": int(np.isnan(M).sum()), "line": line,
           "psi": CS._ps(psi), "frac_mapped": float(mapped.mean())}
    if force_gamma is not None:
        gamma = force_gamma
    else:
        mu_s, sd_s, _ = CS.cshrink_stats(Msel, psi, psi)
        out_s, Z_s, P_s, obs_s, ok_s, mu2, sd2 = _fit(Msel, ctx["V"], ctx["warm"], ctx["base"], mu_s, sd_s, iters=5,
                                                      **fit_kw)
        eh, c, mo, used = _ehat(Z_s, P_s, obs_s, ok_s, Lp, h, mapped, cfix)
        corr = sd2[None, :] * c[None, :] * eh
        (r1, c1), (r2, c2) = p1, p2
        curve = {}
        for g in GAMMAS:
            Y = out_s + g * corr
            e1 = np.sqrt(np.mean((Y[r1, c1] - M[r1, c1]) ** 2))
            e2 = np.sqrt(np.mean((Y[r2, c2] - M[r2, c2]) ** 2))
            curve[g] = float(0.5 * (e1 + e2))
        gamma = min(GAMMAS, key=lambda g: (curve[g], g))
        rr = np.concatenate([r1, r2])
        cc = np.concatenate([c1, c2])
        sel = mapped[cc] & (eh[rr, cc] != 0)
        res = (M[rr, cc] - out_s[rr, cc])[sel]
        pr_ = corr[rr, cc][sel]
        pear = float(np.corrcoef(pr_, res)[0, 1]) if sel.sum() > 2 and np.std(pr_) > 0 else None
        rec.update({"gamma_curve": {str(k): v for k, v in curve.items()}, "mean_n_obs_mapped": mo,
                    "n_cells_used_inner": used, "inner_pearson_ehat_resid": pear, "n_inner_pairs": int(sel.sum()),
                    "psi_scores": {CS._ps(k): v for k, v in pscores.items()}})
    rec["gamma"] = gamma
    mu, sd, _ = CS.cshrink_stats(M, psi, psi)
    out, Z, P, obs, ok, _, sd_ = _fit(M, ctx["V"], ctx["warm"], ctx["base"], mu, sd, iters=10, **fit_kw)
    if gamma != 0:
        eh, c, mo, used = _ehat(Z, P, obs, ok, Lp, h, mapped, cfix)
        miss = ~obs & ok[None, :] & mapped[None, :]
        out = out.copy()
        out[miss] = (out + gamma * sd_[None, :] * c[None, :] * eh)[miss]
        rec.update({"mean_n_obs_mapped_final": mo, "n_cells_used_final": used})
    _log(variant, rec)
    return out


def impute(Xm, ctx):
    return run(Xm, ctx, "main")


def impute_null(Xm, ctx):
    if ctx.get("target") not in LINE:
        return run(Xm, ctx, "null")
    return run(Xm, ctx, "null", line=NULL_LINE)


def impute_perm(Xm, ctx):
    return run(Xm, ctx, "perm", perm=True)


def impute_base(Xm, ctx):
    if ctx.get("target") not in LINE:
        return run(Xm, ctx, "base")
    return run(Xm, ctx, "base", force_gamma=0.0)


# v2 (one revision, column-scale fix only; base unchanged since gamma = 0 never touches c)
def impute_v2(Xm, ctx):
    return run(Xm, ctx, "v2_main", cfix=True)


def impute_null_v2(Xm, ctx):
    if ctx.get("target") not in LINE:
        return run(Xm, ctx, "v2_null", cfix=True)
    return run(Xm, ctx, "v2_null", line=NULL_LINE, cfix=True)


def impute_perm_v2(Xm, ctx):
    return run(Xm, ctx, "v2_perm", perm=True, cfix=True)


# ---------------------------------------------------------------------------------------------------------------
# unit checks
def unit_check():
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H  # noqa: F401
    from scp_graphfs_mcp import phase2 as P2
    from scp_graphfs_mcp.phasec import nfbase as NF

    out = {}
    t = "krull2024"
    Xs, Xm, mask, names = NF._split(t, "MCAR", 100)
    V, warm = H.target_models(t, names)
    b = P2._fill(Xm, P2.loo_knn(Xm))
    b = np.where(np.isnan(b), np.nanmean(Xm), b)
    ctx = {"V": V, "warm": warm, "base": b, "seed": 100, "target": t, "names": names}
    mu, sd, _ = CS.cshrink_stats(Xm, 4.0, 4.0)
    a = NF.adapt_project_mb(Xm, V, warm, b, mu, sd, lam_v=LAM_V, extra=EXTRA, iters=10)
    m = _fit(Xm, V, warm, b, mu, sd, lam_v=LAM_V, extra=EXTRA, iters=10)[0]
    out["U1_fit_vs_adapt_project_mb_bitwise"] = bool(np.array_equal(a, m))
    miss = np.isnan(Xm)
    yb = impute_base(Xm.copy(), ctx)
    yc = CS.run(Xm.copy(), ctx, "unit", log=False)
    d = float(np.max(np.abs(yb[miss] - yc[miss])))
    out["U2_base_vs_cshrink_impute"] = {"max_abs_diff": d, "pass": bool(d < 1e-3)}
    for line in ("U2OS", "U937", "PC9"):
        P = line_pcs(line)
        out[f"pcs_{line}"] = {k: (P[k].tolist() if k == "var_explained" else int(P[k]))
                              for k in P if k not in ("genes", "L", "h")}
        out[f"pcs_{line}"]["mean_h"] = float(P["h"].mean())
    print(json.dumps(out, indent=1))
    (CACHE / "unit_check.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def mapping_report():
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from explore_phaseb import load_target

    rep = {}
    for t, line in LINE.items():
        _, names, _ = load_target(t)
        rep[t] = {}
        for ln in (line, NULL_LINE):
            _, _, mp = protein_loadings(names, ln)
            rep[t][ln] = {"n_proteins": len(names), "n_mapped": int(mp.sum()), "frac_mapped": round(float(mp.mean()), 4)}
    print(json.dumps(rep, indent=1))
    (CACHE / "mapping_report.json").write_text(json.dumps(rep, indent=1) + "\n")
    return rep


if __name__ == "__main__":
    if "--unit" in sys.argv:
        unit_check()
    if "--map" in sys.argv:
        mapping_report()
