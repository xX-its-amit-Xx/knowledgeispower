"""Phase-C idea `varprior` (queen: biological priors).

Protein turnover (measured half-lives, Mathieson et al. 2018 Nat Commun 9:689) and complex buffering (CORUM
membership / size) as the log-sd shrinkage target in the target's per-protein standardisation. Machinery from
`nfbase` (imported, not copied): reference baseline tables, `baselines` (EB fusion), `adapt_project_mb`.

  1. Covariates per target protein j (first accession): h = log2 half-life (median-imputed + indicator mh),
     cx = in any CORUM complex, sz = log median member count of j's complexes; h, sz z-scored over the target.
  2. Reference model (LODO: nfbase.ref_baselines excludes the target and all holdouts): response r_dj = l_dj -
     median_j l_dj (l = log sd), design x = [1, h, mh, cx, sz, f, f^2] with f_dj = n_dj / N_d; pooled ridge
     (lambda = 1, intercept unpenalised).
  3. impute: pl_j = x_j . beta with the target's own observed fraction f_j; location prior = 0 (shrink to a constant);
     nfbase.baselines(psi = 4) -> adapt_project_mb(lam_v = 100, extra = 10).
Variants: impute_const (pl = 0), impute_fonly (x = [1, f, f^2]).
Null: impute_null permutes the covariate rows (h, mh, cx, sz) jointly across the target's p proteins (rng 4242 +
index of the target in sorted FAST), in both the reference fit and the target prediction; f is kept.
Cache: data/cache/phasec/varprior/.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import zipfile
import xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path

import numpy as np

from scp_graphfs_mcp import transfer as T
from scp_graphfs_mcp.phasec.nfbase import (EXTRA, FORBIDDEN, LAM_V, _load, adapt_project_mb, baselines,
                                           ref_baselines)

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "varprior"
XLSX = CACHE / "mathieson2018_MOESM5.xlsx"
XLSX_URL = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-03106-1/MediaObjects/"
            "41467_2018_3106_MOESM5_ESM.xlsx")
XLSX_SHA = "e0503f49230680cbf86fbc86db145ba439c03382e584fde8f0d83241c9826515"
PSI = 4.0
RIDGE = 1.0
FAST = ["petrosius2023", "leduc2022", "weiss2026_liver", "gregoire2023_mixCTRL", "px063590_ADbrain", "sadiku2025"]
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


# ---------------------------------------------------------------------------------------------------------------
# 0. data
def _col(ref):
    s = "".join(c for c in ref if c.isalpha())
    n = 0
    for c in s:
        n = n * 26 + ord(c) - 64
    return n - 1


@lru_cache(maxsize=1)
def half_lives():
    """{gene symbol: nanmedian human half-life (h) over the 8 human half_life columns} from sheet1."""
    js = CACHE / "halflife_by_symbol.json"
    if js.exists():
        return json.loads(js.read_text())
    sha = hashlib.sha256(XLSX.read_bytes()).hexdigest()
    assert sha == XLSX_SHA, f"Mathieson xlsx hash mismatch: {sha}"
    z = zipfile.ZipFile(XLSX)
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    names = [s.get("name") for s in wb.iter(NS + "sheet")]
    assert names[0] == "protein half lives high qual", names
    ss = [("".join(t.text or "" for t in si.iter(NS + "t"))) for si in ET.fromstring(z.read("xl/sharedStrings.xml"))]
    rows = []
    for row in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).iter(NS + "row"):
        r = {}
        for c in row.iter(NS + "c"):
            v = c.find(NS + "v")
            t = c.get("t")
            if t == "e" or v is None:
                val = None
            elif t == "s":
                val = ss[int(v.text)]
            elif t == "str":
                val = v.text
            else:
                val = float(v.text)
            r[_col(c.get("r"))] = val
        rows.append(r)
    hdr = rows[0]
    cols = [k for k, h in hdr.items() if isinstance(h, str) and "half_life" in h and "Mouse" not in h]
    assert len(cols) == 8, [hdr[k] for k in cols]
    out = {}
    for r in rows[1:]:
        g = r.get(0)
        if not isinstance(g, str) or not g:
            continue
        v = np.array([r.get(k) if isinstance(r.get(k), float) else np.nan for k in cols], float)
        if np.isfinite(v).any():
            out.setdefault(g, []).append(float(np.nanmedian(v)))
    out = {g: float(np.median(v)) for g, v in out.items()}
    meta = {"headers_used": [hdr[k] for k in cols], "n_genes": len(out), "median_h": float(np.median(list(out.values())))}
    CACHE.mkdir(parents=True, exist_ok=True)
    js.write_text(json.dumps(out) + "\n")
    (CACHE / "halflife_meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    return out


@lru_cache(maxsize=1)
def hl_by_acc():
    from scp_graphfs_mcp.ids import uniprot_for_symbol

    acc = {}
    n_res = 0
    for g, h in half_lives().items():
        us = uniprot_for_symbol(g)
        n_res += bool(us)
        for u in us:
            acc.setdefault(u, []).append(h)
    info = {"n_symbols": len(half_lives()), "n_resolved": n_res, "n_accessions": len(acc)}
    (CACHE / "mapping_meta.json").write_text(json.dumps(info, indent=1) + "\n")
    return {u: float(np.median(v)) for u, v in acc.items()}


@lru_cache(maxsize=1)
def _complex_sizes():
    from scp_graphfs_mcp.complexes import load_complexes

    df = load_complexes(("corum",))
    idx = {}
    for mem in df["members"]:
        for a in set(mem):
            idx.setdefault(a, []).append(len(mem))
    return idx


# ---------------------------------------------------------------------------------------------------------------
# 1. protein covariates
def covariates(target, names):
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{target}_cov.npz"
    if f.exists():
        d = np.load(f)
        return d["C"]
    accs = T._target_accs(names)
    hl = hl_by_acc()
    cs = _complex_sizes()
    h = np.array([np.log2(hl[a]) if a in hl and hl[a] > 0 else np.nan for a in accs])
    mh = (~np.isfinite(h)).astype(float)
    h = np.where(np.isfinite(h), h, np.nanmedian(h))
    cx = np.array([1.0 if a in cs else 0.0 for a in accs])
    sz = np.array([np.log(np.median(cs[a])) if a in cs else 0.0 for a in accs])

    def zs(v):
        s = v.std()
        return (v - v.mean()) / (s if s > 0 else 1.0)

    C = np.c_[zs(h), mh, cx, zs(sz)]
    np.savez_compressed(f, C=C)
    return C


# ---------------------------------------------------------------------------------------------------------------
# 2. reference log-sd model
def ncells(datasets):
    f = CACHE / "ncells.json"
    have = json.loads(f.read_text()) if f.exists() else {}
    new = False
    for d in datasets:
        assert d not in FORBIDDEN and d not in (T.PHASE_B_HOLDOUTS | T.HOLDOUTS), d
        if d not in have:
            X, _ = _load(d)
            have[d] = int(X.shape[0])
            new = True
    if new:
        CACHE.mkdir(parents=True, exist_ok=True)
        f.write_text(json.dumps(have, indent=1, sort_keys=True) + "\n")
    return have


def _null_perm(target, p):
    if target in FAST:
        k = sorted(FAST).index(target)
    else:  # full pool: deterministic per-target offset
        k = 100 + int(hashlib.sha256(target.encode()).hexdigest()[:6], 16) % 10000
    return np.random.default_rng(4242 + k).permutation(p)


def _design(C, f, mode):
    if mode == "fonly":
        return np.c_[np.ones(len(f)), f, f ** 2]
    if mode in ("nof", "nof_null"):  # v2: biological covariates only, no observed-fraction trend
        return np.c_[np.ones(len(f)), C]
    return np.c_[np.ones(len(f)), C, f, f ** 2]


def ref_table(target, names, mode):
    """Stacked (X, r, dataset index) over finite (d, j) for the reference."""
    R = ref_baselines(target, list(names))
    ds = [str(d) for d in R["datasets"]]
    N = ncells(ds)
    C = covariates(target, names)
    if mode in ("null", "nof_null"):
        C = C[_null_perm(target, len(C))]
    Xs, rs, gs = [], [], []
    for k, d in enumerate(ds):
        assert d not in FORBIDDEN and d != target
        l, vl = R["l"][k], R["vl"][k]
        fin = np.isfinite(l) & np.isfinite(vl)
        r = l - np.nanmedian(l[fin])
        n = 1.0 + 1.0 / (2.0 * vl)
        fd = n / N[d]
        Xs.append(_design(C[fin], fd[fin], mode))
        rs.append(r[fin])
        gs.append(np.full(fin.sum(), k))
    return np.vstack(Xs), np.concatenate(rs), np.concatenate(gs), ds


def ridge(X, y, lam=RIDGE):
    P = lam * np.eye(X.shape[1])
    P[0, 0] = 0.0
    A = X.T @ X + P
    b = np.linalg.solve(A, X.T @ y)
    res = y - X @ b
    s2 = float(res @ res) / max(1, len(y) - X.shape[1])
    Ai = np.linalg.inv(A)
    se = np.sqrt(np.diag(s2 * Ai @ (X.T @ X) @ Ai))
    return b, se


@lru_cache(maxsize=32)
def _beta(target, names_t, mode):
    X, r, _, _ = ref_table(target, names_t, mode)
    return ridge(X, r)[0]


# ---------------------------------------------------------------------------------------------------------------
# 3. candidates
def _pl(Xm, ctx, mode):
    names = ctx["names"]
    p = Xm.shape[1]
    if mode == "const":
        return np.zeros(p)
    beta = _beta(ctx["target"], tuple(names), mode)
    C = covariates(ctx["target"], names)
    if mode in ("null", "nof_null"):
        C = C[_null_perm(ctx["target"], p)]
    f = (~np.isnan(Xm)).sum(0) / Xm.shape[0]
    return _design(C, f, mode) @ beta


def _run(Xm, ctx, mode, psi=PSI):
    p = Xm.shape[1]
    pl = _pl(Xm, ctx, mode)
    mu, sd, _ = baselines(Xm, np.zeros(p), pl, psi=psi)
    return adapt_project_mb(Xm, ctx["V"], ctx["warm"], ctx["base"], mu, sd, lam_v=LAM_V, extra=EXTRA)


def impute(Xm, ctx):
    return _run(Xm, ctx, "true")


def impute_null(Xm, ctx):
    return _run(Xm, ctx, "null")


def impute_const(Xm, ctx):
    return _run(Xm, ctx, "const")


def impute_fonly(Xm, ctx):
    return _run(Xm, ctx, "fonly")


# v2 (one revision, after the main screen): the reference f-trend hurt weiss2026_liver (fonly 0.9873 vs const
# 0.9659 MCAR), so drop f / f^2 from the design; the prior is then the biological covariates alone.
def impute_v2(Xm, ctx):
    return _run(Xm, ctx, "nof")


def impute_v2_null(Xm, ctx):
    return _run(Xm, ctx, "nof_null")


# ---------------------------------------------------------------------------------------------------------------
# 5. diagnostics and unit check
def _spearman_partial(a, b, c):
    from scipy.stats import rankdata

    ra, rb, rc = (rankdata(v) for v in (a, b, c))
    D = np.c_[np.ones(len(rc)), rc]
    ea = ra - D @ np.linalg.lstsq(D, ra, rcond=None)[0]
    eb = rb - D @ np.linalg.lstsq(D, rb, rcond=None)[0]
    return float(np.corrcoef(ea, eb)[0, 1])


def diagnostics():
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from explore_phaseb import load_target

    out = {"cols_full": ["1", "h", "mh", "cx", "sz", "f", "f2"], "targets": {}}
    for t in FAST:
        _, names, _ = load_target(t)
        C = covariates(t, names)
        row = {"coverage_h": float(1 - C[:, 1].mean()), "coverage_cx": float(C[:, 2].mean()), "p": int(len(C))}
        for mode in ("true", "null", "fonly"):
            X, r, g, ds = ref_table(t, tuple(names), mode)
            b, se = ridge(X, r)
            # LODO R^2 within the reference (pooled out-of-dataset prediction)
            pred = np.full(len(r), np.nan)
            for k in np.unique(g):
                tr = g != k
                pred[g == k] = X[g == k] @ ridge(X[tr], r[tr])[0]
            r2 = float(1 - np.sum((r - pred) ** 2) / np.sum((r - r.mean()) ** 2))
            row[mode] = {"beta": b.tolist(), "se": se.tolist(), "lodo_r2": r2, "n_rows": int(len(r)),
                         "n_datasets": len(ds)}
            if mode in ("true", "null"):
                row[mode]["partial_spearman_h_r_given_f"] = _spearman_partial(X[:, 1], r, X[:, 5])
                row[mode]["partial_spearman_h_r_given_f_observed_h_only"] = _spearman_partial(
                    X[X[:, 2] == 0, 1], r[X[:, 2] == 0], X[X[:, 2] == 0, 5])
                row[mode]["partial_spearman_cx_r_given_f"] = _spearman_partial(X[:, 3], r, X[:, 5])
        row["incremental_lodo_r2_bio"] = row["true"]["lodo_r2"] - row["fonly"]["lodo_r2"]
        row["incremental_lodo_r2_null"] = row["null"]["lodo_r2"] - row["fonly"]["lodo_r2"]
        out["targets"][t] = row
        print(t, {k: (round(row[m]["lodo_r2"], 4)) for k, m in (("true", "true"), ("null", "null"), ("fonly", "fonly"))},
              "beta_true", np.round(row["true"]["beta"], 4).tolist(), flush=True)
    (CACHE / "diagnostics.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def unit_check():
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H
    from scp_graphfs_mcp import phase2 as P
    from scp_graphfs_mcp.phasec.nfbase import _split

    out = {}
    for t in FAST:
        Xs, Xm, mask, names = _split(t, "MCAR", 100)
        V, warm = H.target_models(t, names)
        b = P._fill(Xm, P.loo_knn(Xm))
        b = np.where(np.isnan(b), np.nanmean(Xm), b)
        champ = T.adapt_project(Xm, V, warm, b, lam_v=LAM_V, extra=EXTRA)
        ctx = {"names": names, "target": t, "V": V, "warm": warm, "base": b}
        mine = _run(Xm, ctx, "true", psi=0.0)
        d = float(np.nanmax(np.abs(mine - champ)))
        out[t] = {"max_abs_diff_psi0": d, "pass": bool(d < 1e-10)}
        print(t, out[t], flush=True)
    (CACHE / "unit_check.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    if "--diag" in sys.argv:
        diagnostics()
    if "--unit" in sys.argv:
        unit_check()
