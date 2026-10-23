"""Phase-C idea `seqprior` (queen: new data modalities).

Protein-intrinsic sequence and localisation features (UniProt human reviewed FASTA; HPA subcellular 'Main location')
predict each target protein's relative location, log-sd and observed fraction from the LODO reference. The
out-of-fold predictions are protein-specific shrinkage targets for the target's 30-cell mean / log-sd at psi = 4,
inside nfbase.adapt_project_mb.

  A. Data (cache data/cache/phasec/seqprior/): uniprot_human_sp.fasta.gz, subcellular_location.tsv.zip; parsed as
     data only (no code from the cache is executed).
  B. Features per target protein (first accession, isoform suffix stripped): log length, log(1 + unique fully tryptic
     7-30 aa peptides, Keil, 0 missed cleavages), log average MW, GRAVY (Kyte-Doolittle), pI (EMBOSS pKa, bisection),
     frac K+R, frac C, HPA multi-hot over the 20 most frequent main locations, has_hpa, has_seq; median fill for
     missing; columns standardised.
  C. Responses (LODO reference = nfbase.ref_baselines): protein-wise nanmean over the target's acquisition group (fall
     back to all datasets if < 2 datasets or < 500 covered proteins) of y (relative mean), l (log sd) and
     logit(clip(observed fraction, 0.02, 0.98)); weights = number of datasets covering the protein.
  D. Sub-models: HistGradientBoostingRegressor per response (ridge variant: RidgeCV). Covered proteins get a 5-fold
     out-of-fold prediction, uncovered ones the full-fit prediction.
  E. Fusion: OLS of ybar on [1, py, po] and lbar on [1, pl, po] over the cal proteins; nfbase tau2 / BLUP at psi = 4;
     applied to every base_ok protein; then nfbase.adapt_project_mb(lam_v = 100, extra = 10).
Null (impute_null): feature rows permuted across the target's proteins (rng 777) before training and prediction.
Baseline (impute_const): design [1] only (shrink to a constant at psi = 4).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import sys
import zipfile
from collections import Counter
from functools import lru_cache
from pathlib import Path

import numpy as np

from scp_graphfs_mcp import transfer as T
from scp_graphfs_mcp.phasec.nfbase import (EXTRA, FORBIDDEN, LAM_V, MIN_N_CAL, _load, _prior_cached,
                                           adapt_project_mb, ref_baselines)
from scp_graphfs_mcp.phasec.ttprior import META

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "seqprior"
FASTA = CACHE / "uniprot_human_sp.fasta.gz"
FASTA_SHA = "dc8ef8cbf0fc1b857b3be9a559613533ec49824079fb057e31fdfcf578cb229a"
HPA = CACHE / "subcellular_location.tsv.zip"
HPA_SHA = "a9466882bd3328416bcde77c32a6bcb2a17d70183af90e98d032f6f5f781d504"
PSI = 4.0
N_LOC = 20
FAST = ["petrosius2023", "leduc2022", "weiss2026_liver", "gregoire2023_mixCTRL", "px063590_ADbrain", "sadiku2025"]

# average residue masses (Da), Kyte-Doolittle hydropathy, EMBOSS pKa
AVG = {"A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388, "E": 129.1155, "Q": 128.1307,
       "G": 57.0519, "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766,
       "P": 97.1167, "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326, "U": 150.0388,
       "O": 237.3018}
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
PKA_POS = {"Nterm": 8.6, "K": 10.8, "R": 12.5, "H": 6.5}
PKA_NEG = {"Cterm": 3.6, "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
SEQ_COLS = ["log_len", "log1p_tryptic", "log_mw", "gravy", "pI", "frac_KR", "frac_C"]


# ---------------------------------------------------------------------------------------------------------------
# A/B. sequence and localisation features
def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def tryptic_peptides(seq, lo=7, hi=30):
    peps, start = set(), 0
    for i, a in enumerate(seq):
        if a in "KR" and (i + 1 == len(seq) or seq[i + 1] != "P"):
            pep = seq[start:i + 1]
            if lo <= len(pep) <= hi:
                peps.add(pep)
            start = i + 1
    pep = seq[start:]
    if lo <= len(pep) <= hi:
        peps.add(pep)
    return peps


def isoelectric_point(seq):
    cnt = Counter(seq)

    def charge(ph):
        pos = 1.0 / (1 + 10 ** (ph - PKA_POS["Nterm"])) + sum(
            cnt[a] / (1 + 10 ** (ph - k)) for a, k in PKA_POS.items() if a != "Nterm")
        neg = 1.0 / (1 + 10 ** (PKA_NEG["Cterm"] - ph)) + sum(
            cnt[a] / (1 + 10 ** (k - ph)) for a, k in PKA_NEG.items() if a != "Cterm")
        return pos - neg

    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def seq_features(seq):
    n = len(seq)
    kd = [KD[a] for a in seq if a in KD]
    mw = sum(AVG.get(a, 0.0) for a in seq) + 18.015
    return [np.log(n), np.log1p(len(tryptic_peptides(seq))), np.log(mw), float(np.mean(kd)) if kd else np.nan,
            isoelectric_point(seq), (seq.count("K") + seq.count("R")) / n, seq.count("C") / n]


@lru_cache(maxsize=1)
def protein_table():
    """{accession: (gene name, [7 sequence features])} plus HPA {gene name: set(main locations)}, cached as JSON."""
    js = CACHE / "features_by_acc.json"
    if js.exists():
        return json.loads(js.read_text())
    assert _sha(FASTA) == FASTA_SHA, "UniProt FASTA hash mismatch"
    assert _sha(HPA) == HPA_SHA, "HPA zip hash mismatch"
    seqs, cur, acc, gn = {}, [], None, None
    with gzip.open(FASTA, "rt") as f:
        for line in list(f) + [">end"]:
            if line.startswith(">"):
                if acc is not None:
                    seqs[acc] = (gn, "".join(cur))
                hdr = line[1:].split()
                parts = hdr[0].split("|")
                acc = parts[1] if len(parts) >= 2 else None
                gn = next((t[3:] for t in hdr if t.startswith("GN=")), None)
                cur = []
            else:
                cur.append(line.strip())
    feats = {a: [g, seq_features(s)] for a, (g, s) in seqs.items() if s}
    z = zipfile.ZipFile(HPA)
    lines = z.read("subcellular_location.tsv").decode("utf8").splitlines()
    hdr = lines[0].split("\t")
    ig, il = hdr.index("Gene name"), hdr.index("Main location")
    loc = {}
    for ln in lines[1:]:
        r = ln.split("\t")
        if len(r) > il and r[ig]:
            loc[r[ig]] = sorted({x.strip() for x in r[il].split(";") if x.strip()})
    out = {"seq": feats, "hpa": loc,
           "meta": {"n_fasta": len(seqs), "n_hpa_genes": len(loc), "fasta_sha256": FASTA_SHA,
                    "hpa_sha256": HPA_SHA}}
    tmp = CACHE / ".features_by_acc.tmp.json"
    tmp.write_text(json.dumps(out))
    os.replace(tmp, js)
    return out


def _strip(a):
    return a.split("-")[0] if a else a


def features(target, names):
    """(p x k) standardised feature matrix and column names; cached per target."""
    f = CACHE / f"{target}_feat.npz"
    if f.exists():
        d = np.load(f, allow_pickle=False)
        return d["F"], [str(c) for c in d["cols"]]
    tab = protein_table()
    S, H = tab["seq"], tab["hpa"]
    top = [k for k, _ in Counter(x for v in H.values() for x in v).most_common(N_LOC)]
    accs = [_strip(a) for a in T._target_accs(names)]
    p = len(accs)
    Fs = np.full((p, len(SEQ_COLS)), np.nan)
    L = np.zeros((p, N_LOC))
    has_seq, has_hpa = np.zeros(p), np.zeros(p)
    for j, a in enumerate(accs):
        if a in S:
            g, v = S[a]
            Fs[j] = v
            has_seq[j] = 1.0
            if g in H:
                has_hpa[j] = 1.0
                for x in H[g]:
                    if x in top:
                        L[j, top.index(x)] = 1.0
    med = np.nanmedian(Fs, 0)
    Fs = np.where(np.isfinite(Fs), Fs, med[None, :])
    F = np.c_[Fs, L, has_hpa, has_seq]
    sd = F.std(0)
    F = (F - F.mean(0)) / np.where(sd > 0, sd, 1.0)
    cols = SEQ_COLS + [f"hpa:{x}" for x in top] + ["has_hpa", "has_seq"]
    CACHE.mkdir(parents=True, exist_ok=True)
    tmp = CACHE / f".{target}_feat.tmp.npz"
    np.savez_compressed(tmp, F=F, cols=np.array(cols), has_seq=has_seq, has_hpa=has_hpa)
    os.replace(tmp, f)
    return F, cols


# ---------------------------------------------------------------------------------------------------------------
# C. responses
def observed_fraction(target, names, datasets):
    f = CACHE / f"{target}_of.npz"
    if f.exists():
        d = np.load(f, allow_pickle=False)
        if [str(x) for x in d["datasets"]] == list(datasets):
            return d["of"]
    accs = T._target_accs(names)
    OF = np.full((len(datasets), len(accs)), np.nan)
    for k, d in enumerate(datasets):
        assert d not in FORBIDDEN and d not in (T.PHASE_B_HOLDOUTS | T.HOLDOUTS) and d != target, d
        X, a = _load(d)
        pos = {x: i for i, x in enumerate(a)}
        cols = [(j, pos[x]) for j, x in enumerate(accs) if x in pos]
        if cols:
            tj, sj = zip(*cols)
            OF[k, list(tj)] = (~np.isnan(X[:, list(sj)])).sum(0) / X.shape[0]
    tmp = CACHE / f".{target}_of.tmp.npz"
    np.savez_compressed(tmp, of=OF, datasets=np.array(datasets))
    os.replace(tmp, f)
    return OF


def responses(target, names):
    R = ref_baselines(target, list(names))
    ds = [str(d) for d in R["datasets"]]
    OF = observed_fraction(target, names, ds)
    gstar = META.get(target, {}).get("acq", "unknown")
    rows = np.array([META.get(d, {}).get("acq", "unknown") == gstar for d in ds])
    covered = np.isfinite(R["y"][rows]).any(0).sum() if rows.any() else 0
    group = "acq"
    if rows.sum() < 2 or covered < 500:
        rows, group = np.ones(len(ds), bool), "all"
    lo = np.log(np.clip(OF, 0.02, 0.98) / (1 - np.clip(OF, 0.02, 0.98)))
    out = {}
    for key, A in (("y", R["y"]), ("l", R["l"]), ("o", lo)):
        A = A[rows]
        cnt = np.isfinite(A).sum(0)
        with np.errstate(all="ignore"):
            out[key] = (np.where(cnt > 0, np.nanmean(np.where(cnt[None, :] > 0, A, 0.0), 0), np.nan), cnt)
    info = {"group": group, "gstar": gstar, "datasets": [d for d, r in zip(ds, rows) if r]}
    return out, info


# ---------------------------------------------------------------------------------------------------------------
# D. sub-models
def _model(kind):
    if kind == "ridge":
        from sklearn.linear_model import RidgeCV

        return RidgeCV(alphas=np.logspace(-2, 3, 12))
    from sklearn.ensemble import HistGradientBoostingRegressor

    return HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=40,
                                         l2_regularization=1.0, random_state=0)


def _oof_predict(F, r, w, kind):
    from sklearn.model_selection import KFold

    cov = np.isfinite(r) & (w > 0)
    idx = np.where(cov)[0]
    pred = np.full(len(r), np.nan)
    oof = np.full(len(idx), np.nan)
    for tr, te in KFold(5, shuffle=True, random_state=0).split(idx):
        m = _model(kind).fit(F[idx[tr]], r[idx[tr]], sample_weight=w[idx[tr]])
        oof[te] = m.predict(F[idx[te]])
    pred[idx] = oof
    m = _model(kind).fit(F[idx], r[idx], sample_weight=w[idx])
    pred[~cov] = m.predict(F[~cov])
    ww, rr = w[idx], r[idx]
    mu = np.average(rr, weights=ww)
    r2 = 1 - np.sum(ww * (rr - oof) ** 2) / np.sum(ww * (rr - mu) ** 2)
    return pred, float(r2), int(cov.sum())


def null_perm(p):
    return np.random.default_rng(777).permutation(p)


@lru_cache(maxsize=64)
def predictions(target, names_t, mode):
    """mode: 'true' | 'null' | 'ridge' | 'ridge_null'. Returns (py, pl, po, info); disk-cached."""
    f = CACHE / f"{target}_pred_{mode}.npz"
    if f.exists():
        d = np.load(f, allow_pickle=False)
        return d["py"], d["pl"], d["po"], json.loads(str(d["info"]))
    names = list(names_t)
    F, cols = features(target, names)
    if mode.endswith("null"):
        F = F[null_perm(len(F))]
    kind = "ridge" if mode.startswith("ridge") else "gbm"
    resp, info = responses(target, names)
    preds = {}
    for key in ("y", "l", "o"):
        r, cnt = resp[key]
        preds[key], r2, ncov = _oof_predict(F, r, cnt.astype(float), kind)
        info[f"oof_r2_{key}"] = r2
        info[f"n_cov_{key}"] = ncov
    info.update(mode=mode, model=kind, p=int(len(F)), n_features=len(cols))
    tmp = CACHE / f".{target}_pred_{mode}.tmp.npz"
    np.savez_compressed(tmp, py=preds["y"], pl=preds["l"], po=preds["o"], info=json.dumps(info))
    os.replace(tmp, f)
    return preds["y"], preds["l"], preds["o"], info


# ---------------------------------------------------------------------------------------------------------------
# E. fusion (local copy of nfbase.baselines with a multi-regressor calibration)
def _blup(stat, var, D, cal, psi):
    sel = cal & np.isfinite(stat) & np.isfinite(D).all(1)
    if sel.sum() < 10 + D.shape[1]:
        return None, None, {"n_cal": int(sel.sum())}
    coef = np.linalg.lstsq(D[sel], stat[sel], rcond=None)[0]
    pm = D @ coef
    tau2 = max(1e-4, float(np.var(stat[sel] - pm[sel]) - psi * np.mean(var[sel])))
    w = tau2 / (tau2 + psi * var)
    return pm, w, {"coef": coef.tolist(), "tau2": tau2, "n_cal": int(sel.sum())}


def baselines(Xm, dloc, dscale, psi=PSI):
    obs = ~np.isnan(Xm)
    n = obs.sum(0)
    with np.errstate(all="ignore"):
        xbar = np.nanmean(Xm, 0)
        s = np.nanstd(Xm, 0)
    mu, sd = xbar.copy(), s.copy()
    ok = np.isfinite(xbar)
    base_ok = ok & (n >= 2) & np.isfinite(s) & (s > 1e-6)
    cal = base_ok & (n >= MIN_N_CAL)
    info = {}
    c_t = np.median(xbar[ok & (n >= 2)])
    ybar = xbar - c_t
    vy = np.where(base_ok, s ** 2 / np.maximum(n, 1), np.nan)
    pm, w, info["loc"] = _blup(ybar, vy, dloc, cal, psi)
    if pm is not None:
        post = c_t + w * ybar + (1 - w) * pm
        use = base_ok & np.isfinite(post)
        mu = np.where(use, post, mu)
    lbar = np.where(base_ok, np.log(np.where(base_ok, s, 1.0)), np.nan)
    vl = np.where(base_ok, 1.0 / (2.0 * np.maximum(n - 1, 1)), np.nan)
    pm, w, info["scale"] = _blup(lbar, vl, dscale, cal, psi)
    if pm is not None:
        post = w * lbar + (1 - w) * pm
        use = base_ok & np.isfinite(post)
        sd = np.where(use, np.exp(np.where(use, post, 0.0)), sd)
    return mu, sd, info


def designs(ctx, p, mode):
    one = np.ones((p, 1))
    if mode == "const":
        return one, one
    pmode = {"true": "true", "null": "null", "ridge": "ridge", "ridge_null": "ridge_null",
             "nf": "true", "nf_null": "null"}[mode]
    py, pl, po, _ = predictions(ctx["target"], tuple(ctx["names"]), pmode)
    dl, ds = np.c_[one, py, po], np.c_[one, pl, po]
    if mode.startswith("nf"):  # optional: add nfbase's reference prior where finite (0 + indicator elsewhere)
        ny, nl, _ = _prior_cached(ctx["target"], tuple(ctx["names"]), "true")
        hy, hl = np.isfinite(ny).astype(float), np.isfinite(nl).astype(float)
        dl = np.c_[dl, np.nan_to_num(ny), hy]
        ds = np.c_[ds, np.nan_to_num(nl), hl]
    return dl, ds


def _run(Xm, ctx, mode, psi=PSI):
    dl, ds = designs(ctx, Xm.shape[1], mode)
    mu, sd, _ = baselines(Xm, dl, ds, psi)
    return adapt_project_mb(Xm, ctx["V"], ctx["warm"], ctx["base"], mu, sd, lam_v=LAM_V, extra=EXTRA)


# F. entry points
def impute(Xm, ctx):
    return _run(Xm, ctx, "true")


def impute_null(Xm, ctx):
    return _run(Xm, ctx, "null")


def impute_const(Xm, ctx):
    return _run(Xm, ctx, "const")


def impute_ridge(Xm, ctx):
    return _run(Xm, ctx, "ridge")


def impute_ridge_null(Xm, ctx):
    return _run(Xm, ctx, "ridge_null")


def impute_nf(Xm, ctx):
    return _run(Xm, ctx, "nf")


def impute_nf_null(Xm, ctx):
    return _run(Xm, ctx, "nf_null")


# ---------------------------------------------------------------------------------------------------------------
# G. diagnostics and unit checks
def diagnostics(targets=FAST):
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from explore_phaseb import load_target

    out = {}
    for t in targets:
        _, names, _ = load_target(t)
        F, cols = features(t, names)
        d = np.load(CACHE / f"{t}_feat.npz")
        row = {"p": int(len(F)), "n_features": len(cols), "coverage_seq": float(d["has_seq"].mean()),
               "coverage_hpa": float(d["has_hpa"].mean())}
        for mode in ("true", "null", "ridge", "ridge_null"):
            row[mode] = predictions(t, tuple(names), mode)[3]
        out[t] = row
        print(t, {m: {k: round(row[m][k], 4) for k in ("oof_r2_y", "oof_r2_l", "oof_r2_o")}
                  for m in ("true", "null", "ridge", "ridge_null")}, row[mode]["group"], flush=True)
    res = ROOT / "results" / "phaseC" / "seqprior_diag.json"
    res.write_text(json.dumps({"columns": cols, "targets": out}, indent=1) + "\n")
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
        dd = float(np.nanmax(np.abs(mine - champ)))
        out[t] = {"max_abs_diff_psi0": dd, "pass": bool(dd < 1e-10)}
        print(t, out[t], flush=True)
    (CACHE / "unit_check.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    if "--build" in sys.argv:
        print(protein_table()["meta"])
    if "--diag" in sys.argv:
        diagnostics()
    if "--unit" in sys.argv:
        unit_check()
