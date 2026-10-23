"""subthr (Phase C round 4, queen: new data modalities): "forced photometry" -- sub-threshold detection evidence
added as low-precision observations at naturally missing entries, on the tgtshrink.impute_null base (the cshrink
stand-in: selected-psi Fay-Herriot shrinkage of the target's mean / log-sd toward a constant-like permuted trend).

Evidence sources (no downloads; everything is already cached):
  leduc2022  unfiltered per-batch PSM assays of leduc2022_pSCoPE.Rda: PSMs with 0.01 < dart_qval <= 0.1 (not
             carried into proteins_norm2), not Reverse / contaminant, Leading.razor.protein.id in the target's
             names. Per (cell, protein) at X-NaN entries:
               r = log2(reporter intensity), zeros -> NaN
               dev = r - nanmedian over the batch's target cells (per PSM)
               dev -= cell offset (median over the cell's PASSING PSMs, dart_qval <= 0.01, of the same dev;
                      X is per-cell median-centred, so the evidence is put on the same per-cell footing)
               D = nanmedian over the entry's sub-threshold PSMs of dev;  P = 1 - min dart_PEP, clipped to [0, 1]
  sadiku2025 Spectronaut pivot '[k] <run>.raw.PG.PEP (Run-Wise)', aligned by replaying pepev.build_sadiku2025;
             at X-NaN entries P = max(0, 1 - PEP / 0.5), D = NaN (detection score, no quantity).
  gregoire2023_mixCTRL: checked (see `build`), included only if >= 1 evidence entry per 30-cell split on average.

Build:   python -m scp_graphfs_mcp.phasec.subthr build      -> data/cache/phasec/subthr/<target>_E.npz
Diag:    python -m scp_graphfs_mcp.phasec.subthr diag       -> results/phaseC/subthr_diag.json (kappa = 0 unit check)

Model (per call): pseudo-values from Xm only (leduc v = colmean_j(Xm) + D; sadiku v = q10_j(Xm), cell q10 if
n_j < 5), weights pi = kappa * P; weighted target stats (n_eff = n + sum pi) feed tgtshrink's 'null' FH shrinkage;
weighted U-step over observed + evidence entries in nfbase.adapt_project_mb (V-step on observed only).
Selection on tgtshrink.selection_mask (50/50 MCAR/MNAR, seed + 31337): psi as the base (kappa = 0), then kappa in
{0, 0.1, 0.3, 1} at that psi. Null: within each cell, the (P, D) records are permuted among that cell's evidence
positions (seed + 777 + i) and the pseudo-value is recomputed for the protein each record lands on.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

from scp_graphfs_mcp.phasec import tgtshrink as TS
from scp_graphfs_mcp.phasec.nfbase import EXTRA, LAM_V, adapt_project_mb

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "subthr"
RAW = ROOT / "data" / "cache" / "corpus_raw"
PEPEV = ROOT / "data" / "cache" / "phasec" / "pepev"
OUT = ROOT / "results" / "phaseC"
sys.path.insert(0, str(ROOT / "benchmarks"))
sys.path.insert(0, str(ROOT / "scripts" / "corpus_ext"))

KAPPA_GRID = (0.0, 0.1, 0.3, 1.0)
EVIDENCE_TARGETS = ("leduc2022", "sadiku2025", "gregoire2023_mixCTRL")


def _sha(X):
    a = np.ascontiguousarray(np.asarray(X, np.float32))
    return hashlib.sha256(a.tobytes()).hexdigest()


def _save(target, X, names, P, D, meta):
    obs = ~np.isnan(X)
    Pf = np.isfinite(P)
    assert np.isnan(X[Pf]).all(), "evidence placed at an observed entry"
    assert not (np.isfinite(D) & ~Pf).any()
    nan = ~obs
    meta = {**meta, "n_evidence": int(Pf.sum()), "n_evidence_P_gt0": int((Pf & (P > 0)).sum()),
            "frac_X_nan_covered": float(Pf[nan].mean()), "frac_X_nan_covered_P_gt0": float((Pf & (P > 0))[nan].mean()),
            "sum_P": float(np.nansum(P)), "n_cells": int(X.shape[0]), "n_proteins": int(X.shape[1])}
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(CACHE / f"{target}_E.npz", P=P.astype(np.float32), D=D.astype(np.float32), sha=_sha(X),
                        names=np.array(names), meta=json.dumps(meta))
    print(target, json.dumps(meta), flush=True)
    return meta


# ======================================================================================================================
# (A) build
# ======================================================================================================================
def _leduc_extract():
    """Per-batch PSM assays (dart_qval <= 0.1 only), cached as one npz."""
    f = CACHE / "leduc2022_psm_extract.npz"
    if f.exists():
        return dict(np.load(f, allow_pickle=False))
    from rda_stream import dataframe, load_rda, qf_assays, sce_matrix

    qf = load_rda(PEPEV / "dl_leduc" / "leduc2022_pSCoPE.Rda", max_str=10_000_000)[0][1]
    assays = qf_assays(qf)
    Ms, cols, batch, q, pep, prot, bad = [], [], [], [], [], [], []
    n_all = 0
    for a in assays:
        if not re.match(r"^[ew]AL\d+$", a):
            continue
        sce = assays[a]
        m, dim = sce_matrix(sce)
        M = np.asarray(m.value, float).reshape(dim[1], dim[0]).T
        cn = list(sce.attr["colData"].attr["rownames"].value)
        assert cn == [f"{a}RI{k}" for k in range(1, 19)], (a, cn)
        _, rd = dataframe(sce.attr["rowRanges"].attr["elementMetadata"])
        dq = np.asarray(rd["dart_qval"], float)
        n_all += len(dq)
        keep = np.isfinite(dq) & (dq <= 0.1)
        rev = np.array([str(x) == "+" for x in rd["Reverse"]])
        con = np.array([str(x) == "+" for x in rd["Potential.contaminant"]])
        Ms.append(M[keep].astype(np.float32))
        batch += [a] * int(keep.sum())
        q.append(dq[keep])
        pep.append(np.asarray(rd["dart_PEP"], float)[keep])
        prot += list(np.asarray(rd["Leading.razor.protein.id"]).astype(str)[keep])
        bad.append((rev | con)[keep])
    out = {"M": np.vstack(Ms), "batch": np.array(batch), "dart_qval": np.concatenate(q),
           "dart_PEP": np.concatenate(pep), "prot": np.array(prot), "rev_or_con": np.concatenate(bad),
           "n_psm_all": np.array(n_all)}
    np.savez_compressed(f, **out)
    return out


def build_leduc2022():
    from scipy.stats import spearmanr

    from explore_phaseb import load_target
    from scp_graphfs_mcp.datasets import load_leduc2022

    X, names, _ = load_target("leduc2022")
    cells = list(load_leduc2022().obs_names)
    assert len(cells) == X.shape[0]
    ci = {c: i for i, c in enumerate(cells)}
    pj = {p: j for j, p in enumerate(names)}
    ex = _leduc_extract()
    M, batch, dq, pep, prot, bad = ex["M"], ex["batch"], ex["dart_qval"], ex["dart_PEP"], ex["prot"], ex["rev_or_con"]
    P = np.full(X.shape, np.nan)
    D = np.full(X.shape, np.nan)
    diag = {"n_psm_all_batches": int(ex["n_psm_all"]), "n_psm_qval_le_0.1": int(len(dq)),
            "n_psm_subthr_0.01_0.1": int(((dq > 0.01) & (dq <= 0.1)).sum()),
            "n_psm_subthr_not_rev_con": int(((dq > 0.01) & (dq <= 0.1) & ~bad).sum())}
    pass_pairs_dev, pass_pairs_x, pass_pairs_dev_raw = [], [], []
    n_sub_used, n_batches = 0, 0
    cell_off_all = []
    # accumulate per (cell, protein) lists
    acc_dev, acc_pep = {}, {}
    for b in np.unique(batch):
        rows = np.where(batch == b)[0]
        chans = [(k, ci[f"{b}RI{k + 1}"]) for k in range(18) if f"{b}RI{k + 1}" in ci]
        if not chans:
            continue
        n_batches += 1
        ks = np.array([k for k, _ in chans])
        cis = np.array([i for _, i in chans])
        R = M[rows][:, ks].astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            R = np.where(R > 0, np.log2(R), np.nan)
        with np.errstate(all="ignore"):
            med = np.nanmedian(R, 1, keepdims=True)
        dev_raw = R - med
        q_b, pep_b, prot_b, bad_b = dq[rows], pep[rows], prot[rows], bad[rows]
        passing = (q_b <= 0.01) & ~bad_b
        with np.errstate(all="ignore"):
            off = np.nanmedian(dev_raw[passing], 0)  # per cell offset from passing PSMs (raw data only)
        off = np.where(np.isfinite(off), off, 0.0)
        cell_off_all += off.tolist()
        dev = dev_raw - off[None, :]
        # diagnostic on passing PSMs whose protein is observed in X
        inp = np.array([p in pj for p in prot_b])
        Xb = X[cis]
        with np.errstate(all="ignore"):
            xmed = np.nanmedian(Xb, 0)
        for r in np.where(passing & inp)[0]:
            j = pj[prot_b[r]]
            xd = Xb[:, j] - xmed[j]
            ok = np.isfinite(xd) & np.isfinite(dev[r])
            if ok.any():
                pass_pairs_dev.append(dev[r, ok])
                pass_pairs_dev_raw.append(dev_raw[r, ok])
                pass_pairs_x.append(xd[ok])
        sub = (q_b > 0.01) & (q_b <= 0.1) & ~bad_b & inp
        for r in np.where(sub)[0]:
            j = pj[prot_b[r]]
            for c, i in enumerate(cis):
                if not np.isnan(X[i, j]) or not np.isfinite(dev[r, c]):
                    continue
                acc_dev.setdefault((i, j), []).append(dev[r, c])
                acc_pep.setdefault((i, j), []).append(pep_b[r])
            n_sub_used += 1
    for (i, j), v in acc_dev.items():
        D[i, j] = float(np.median(v))
        P[i, j] = float(np.clip(1.0 - np.min(acc_pep[(i, j)]), 0.0, 1.0))
    dv, dvr, xv = (np.concatenate(pass_pairs_dev), np.concatenate(pass_pairs_dev_raw), np.concatenate(pass_pairs_x))
    diag.update({"n_batches_with_target_cells": n_batches, "n_subthr_psms_in_names": n_sub_used,
                 "passing_psm_scale_check": {
                     "n_pairs": int(len(xv)),
                     "spearman_dev_offset_corrected_vs_X_within_batch_dev": float(spearmanr(dv, xv)[0]),
                     "spearman_dev_raw_vs_X_within_batch_dev": float(spearmanr(dvr, xv)[0]),
                     "sd_dev_offset_corrected": float(np.std(dv)), "sd_dev_raw": float(np.std(dvr)),
                     "sd_X_within_batch_dev": float(np.std(xv)),
                     "mad_ratio_dev_over_X": float(np.median(np.abs(dv)) / np.median(np.abs(xv)))},
                 "cell_offset_quantiles_p10_p50_p90": np.percentile(cell_off_all, [10, 50, 90]).tolist(),
                 "D_quantiles_p10_p50_p90": np.nanpercentile(D, [10, 50, 90]).tolist(),
                 "P_quantiles_p10_p50_p90": np.nanpercentile(P, [10, 50, 90]).tolist()})
    return _save("leduc2022", X, names, P, D, diag)


def build_sadiku2025():
    import pandas as pd

    from scp_graphfs_mcp.phasec.pepev import _replay_save as replay
    from explore_phaseb import load_target
    from scp_graphfs_mcp.ids import group_accessions

    rep = RAW / "sadiku2025" / "20250204_112949_gbm_sc_full_fasta_19.4_Report.tsv"
    d = pd.read_csv(rep, sep="\t", low_memory=False)
    design = pd.read_csv(RAW / "sadiku2025" / "experimentalDesignTemplate.txt", sep="\t")
    patient = dict(zip(design["raw file"], design["patient"]))
    cols, pcols = [], []
    for c in d.columns:
        if not c.endswith(".raw.PG.Quantity"):
            continue
        run = re.sub(r"^\[\d+\] ", "", c)[: -len(".PG.Quantity")]
        if run in patient and "0cells" not in run:
            cols.append(c)
            pc = c[: -len("PG.Quantity")] + "PG.PEP (Run-Wise)"
            assert pc in d.columns, pc
            pcols.append(pc)
    X = d[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float).T
    PEP = d[pcols].apply(pd.to_numeric, errors="coerce").to_numpy(float).T
    with np.errstate(divide="ignore", invalid="ignore"):
        X = np.where(X > 0, np.log2(X), np.nan)
    accs = [(group_accessions(str(g)) or [""])[0] for g in d["PG.ProteinGroups"]]
    seen_sc = ~np.isnan(X).all(0)
    X, PEP, accs = X[:, seen_sc], PEP[:, seen_sc], [a for a, k in zip(accs, seen_sc) if k]
    Xr, accs, PEP = replay(X, accs, PEP)
    X0, names, _ = load_target("sadiku2025")
    assert list(accs) == list(names) and Xr.shape == X0.shape
    assert np.allclose(np.asarray(Xr, np.float32), X0, equal_nan=True, atol=1e-5, rtol=0), "X_rebuilt != X_corpus"
    nan = np.isnan(X0)
    obs = ~nan
    diag = {"frac_entries_with_PEP": float(np.isfinite(PEP).mean()),
            "max_PEP_observed": float(np.nanmax(PEP[obs])),
            "PEP_quantiles_missing_p1_p10_p50": np.nanpercentile(PEP[nan], [1, 10, 50]).tolist(),
            "frac_missing_PEP_lt_0.05_0.15_0.5": [float((PEP[nan] < t).mean()) for t in (0.05, 0.15, 0.5)]}
    from scipy.stats import spearmanr

    cc = []
    for j in range(X0.shape[1]):
        o = obs[:, j] & np.isfinite(PEP[:, j])
        if o.sum() > 10 and np.std(PEP[o, j]) > 0:
            cc.append(spearmanr(PEP[o, j], X0[o, j])[0])
    diag["spearman_PEP_value_within_protein_median"] = float(np.nanmedian(cc))
    P = np.where(nan & np.isfinite(PEP), np.maximum(0.0, 1.0 - PEP / 0.5), np.nan)
    D = np.full(X0.shape, np.nan)
    return _save("sadiku2025", X0, names, P, D, diag)


def check_gregoire2023_mixCTRL():
    """Sage PSM assays: PSMs with protein_fdr > 0.01 at proteins absent from proteins_<batch>_norm_log for that cell.
    Returns the diagnostic dict; builds an E cache only if >= 1 evidence entry per 30-cell split on average."""
    from rda_stream import dataframe, load_rda, qf_assays, sce_matrix

    from explore_phaseb import load_target
    from scp_graphfs_mcp.ids import group_accessions

    qf = load_rda(RAW / "gregoire2023_mixCTRL" / "gregoire2023_mixCTRL.Rda", max_str=10_000_000)[0][1]
    assays = qf_assays(qf)
    names_all = list(assays)
    _, sm = dataframe(qf.attr["sampleMap"])
    X, names, _ = load_target("gregoire2023_mixCTRL")
    out = {"assays": names_all[:12] + ["..."] + names_all[-12:], "n_assays": len(names_all)}
    psm_assays = [a for a in names_all if not re.match(r"^(proteins|peptides)_", a)]
    out["psm_assays"] = psm_assays[:5] + ["..."] if len(psm_assays) > 5 else psm_assays
    # map assay columns -> primary sample ids, and primary -> rows of the corpus matrix via the pepev replay
    # (gregoire's corpus rows are the protein-assay single-cell columns in batch order; we only count evidence entries
    # against each batch's protein assay here, which is what decides inclusion)
    n_ev, n_cells_sc, n_rows_fdr, n_absent = 0, 0, 0, 0
    per_assay = {}
    prot_cols = {}
    for a in names_all:
        mm = re.match(r"^proteins_(.+)_norm_log$", a)
        if mm:
            sce = assays[a]
            m, dim = sce_matrix(sce)
            Mp = np.asarray(m.value, float).reshape(dim[1], dim[0]).T
            from scp_graphfs_mcp.phasec.pepev import _sce_rows

            rows = list(_sce_rows(sce, m))
            cols = list(sce.attr["colData"].attr["rownames"].value)
            k = np.array(sm["assay"]) == a
            c2p = dict(zip(np.array(sm["colname"])[k], np.array(sm["primary"])[k]))
            prot_cols[mm.group(1)] = (Mp, {r: i for i, r in enumerate(rows)}, {c2p.get(c, ""): ci for ci, c in enumerate(cols)})
    for a in psm_assays:
        sce = assays[a]
        try:
            m, dim = sce_matrix(sce)
            _, rd = dataframe(sce.attr["rowRanges"].attr["elementMetadata"])
        except Exception as e:  # noqa: BLE001
            per_assay[a] = f"unreadable: {e!r}"
            continue
        keys = list(rd.keys())
        if "protein_fdr" not in rd:
            per_assay[a] = {"keys": keys[:40]}
            continue
        M = np.asarray(m.value, float).reshape(dim[1], dim[0]).T
        cols = list(sce.attr["colData"].attr["rownames"].value)
        k = np.array(sm["assay"]) == a
        c2p = dict(zip(np.array(sm["colname"])[k], np.array(sm["primary"])[k]))
        pf = np.asarray(rd["protein_fdr"], float)
        protcol = "proteins" if "proteins" in rd else keys[0]
        prots = np.asarray(rd[protcol]).astype(str)
        sub = np.where(pf > 0.01)[0]
        n_rows_fdr += len(sub)
        bmatch = [b for b in prot_cols if b in a]
        rec = {"n_psm": int(len(pf)), "n_protein_fdr_gt_0.01": int(len(sub)), "batch": bmatch}
        if bmatch and len(sub):
            Mp, rix, cix = prot_cols[bmatch[0]]
            ne = 0
            for r in sub:
                for c in range(M.shape[1]):
                    if not np.isfinite(M[r, c]) or M[r, c] <= 0:
                        continue
                    pc = cix.get(c2p.get(cols[c], ""))
                    if pc is None:
                        continue
                    acc = prots[r]
                    pr = rix.get(acc)
                    if pr is None or np.isnan(Mp[pr, pc]):
                        ne += 1
            rec["evidence_entries_upper_bound"] = ne
            n_ev += ne
            # PSM-channel entries whose leading accession is absent from the batch's protein assay altogether
            batch_accs = {x for rr in rix for x in (group_accessions(str(rr)) or [])}
            na = 0
            for r in sub:
                acc = (group_accessions(prots[r]) or [""])[0]
                if acc not in batch_accs and acc in set(names):
                    na += int((np.isfinite(M[r]) & (M[r] > 0)).sum())
            rec["entries_at_proteins_absent_from_batch_assay_in_names"] = na
            n_absent += na
        per_assay[a] = rec
    out["per_psm_assay"] = per_assay
    out["n_psm_protein_fdr_gt_0.01"] = n_rows_fdr
    out["n_evidence_entries_upper_bound"] = n_ev
    out["n_entries_at_proteins_absent_from_batch_assay_in_names"] = n_absent
    n_ev = n_ev + n_absent
    out["n_cells_corpus"] = int(X.shape[0])
    out["evidence_per_30cell_split_upper_bound"] = float(n_ev * 30 / X.shape[0])
    out["include"] = bool(out["evidence_per_30cell_split_upper_bound"] >= 1.0)
    return out


def build():
    res = {}
    res["leduc2022"] = build_leduc2022()
    res["sadiku2025"] = build_sadiku2025()
    try:
        res["gregoire2023_mixCTRL"] = check_gregoire2023_mixCTRL()
    except Exception as e:  # noqa: BLE001
        res["gregoire2023_mixCTRL"] = {"error": repr(e)}
    print(json.dumps(res.get("gregoire2023_mixCTRL"), default=str)[:3000])
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "subthr_build.json").write_text(json.dumps(res, indent=1, default=str) + "\n")


# ======================================================================================================================
# (B) call-time alignment
# ======================================================================================================================
_TGT = {}


def _target(target):
    if target not in _TGT:
        from explore_phaseb import load_target

        X, names, labels = load_target(target)
        f = CACHE / f"{target}_E.npz"
        E = None
        if f.exists():
            d = np.load(f)
            assert str(d["sha"]) == _sha(X), f"{target}: E cache is stale"
            assert list(d["names"]) == list(names)
            E = (d["P"].astype(float), d["D"].astype(float))
        _TGT[target] = (X, names, labels, E)
    return _TGT[target]


def _aligned(Xm, ctx):
    from explore_phaseb import subsample

    X, names, labels, E = _target(ctx["target"])
    if E is None:
        return None
    assert list(names) == list(ctx["names"])
    idx = subsample(X.shape[0], labels, 30, ctx["seed"])
    Xs = X[idx]
    om, os_ = ~np.isnan(Xm), ~np.isnan(Xs)
    assert Xm.shape == Xs.shape and not (om & ~os_).any(), "obs(Xm) is not a subset of obs(X[idx])"
    assert np.allclose(Xm[om], Xs[om]), "Xm values differ from X[idx]"
    P, D = E[0][idx], E[1][idx]
    assert np.isnan(Xm[np.isfinite(P)]).all()
    return P, D


def _permute_records(P, D, ok, seed):
    """Null: within each cell, permute (P, D) records among that cell's evidence positions (on usable columns)."""
    P2, D2 = P.copy(), D.copy()
    for i in range(P.shape[0]):
        pos = np.where(np.isfinite(P[i]) & ok)[0]
        if len(pos) < 2:
            continue
        perm = np.random.default_rng(seed + 777 + i).permutation(len(pos))
        P2[i, pos] = P[i, pos[perm]]
        D2[i, pos] = D[i, pos[perm]]
    return P2, D2


def _pseudo(M, P, D, kind):
    """Pseudo-values on the matrix scale, from M (the matrix being fitted) only."""
    obs = ~np.isnan(M)
    n = obs.sum(0)
    with np.errstate(all="ignore"):
        if kind == "quant":
            m = np.nanmean(M, 0)
            V = m[None, :] + D
        else:
            q = np.nanpercentile(M, 10, axis=0)
            qc = np.nanpercentile(M, 10, axis=1)
            V = np.where((n >= 5)[None, :], q[None, :], qc[:, None]) * np.ones_like(M)
    return np.where(np.isfinite(P), V, np.nan)


# ======================================================================================================================
# (D) base with evidence
# ======================================================================================================================
def target_stats_w(M, Ev, Wev):
    """tgtshrink.target_stats with evidence entries (values Ev, weights Wev, both NaN/0 elsewhere)."""
    S = TS.target_stats(M)
    obs = ~np.isnan(M)
    ev = np.isfinite(Ev) & (Wev > 0) & ~obs
    if not ev.any():
        return S
    W = np.where(obs, 1.0, np.where(ev, Wev, 0.0))
    Xv = np.where(obs, M, np.where(ev, Ev, 0.0))
    n = obs.sum(0)
    sw = W.sum(0)
    with np.errstate(all="ignore"):
        xbar = (W * Xv).sum(0) / sw
        s = np.sqrt((W * (Xv - xbar[None, :]) ** 2).sum(0) / sw)
    ok = S["ok"]  # columns usable = observed in M (unchanged)
    xbar = np.where(ok, xbar, np.nan)
    s = np.where(ok, s, np.nan)
    n_eff = n + np.where(ev, Wev, 0.0).sum(0)
    base_ok = ok & (n >= 2) & np.isfinite(s) & (s > 1e-6)
    cal = base_ok & (n >= TS.MIN_N_CAL)
    c_t = float(np.median(xbar[ok & (n >= 2)]))
    return {"n": n, "f": n / M.shape[0], "xbar": xbar, "s": s, "ok": ok, "base_ok": base_ok, "cal": cal,
            "c_t": c_t, "ybar": xbar - c_t, "vy": np.where(base_ok, s ** 2 / np.maximum(n_eff, 1e-9), np.nan),
            "l": np.where(base_ok, np.log(np.where(base_ok, s, 1.0)), np.nan),
            "vl": np.where(base_ok, 1.0 / (2.0 * np.maximum(n_eff - 1, 1)), np.nan)}


def shrink_w(M, psi, seed, S):
    """tgtshrink.shrink(mode='null', cover=None) on the precomputed stats S."""
    base_ok, cal, n = S["base_ok"], S["cal"], S["n"]
    mu, sd = S["xbar"].copy(), S["s"].copy()
    if psi == 0 or cal.sum() < 10:
        return mu, sd, {"psi": psi}
    use = base_ok.copy()
    m, f = S["xbar"].copy(), S["f"].astype(float).copy()
    rng = np.random.default_rng(seed + 4242)
    ci = np.where(cal)[0]
    perm = ci[rng.permutation(len(ci))]
    m2, f2 = m.copy(), f.copy()
    m2[ci], f2[ci] = m[perm], f[perm]
    nc = np.where(base_ok & ~cal)[0]
    if len(nc):
        src = ci[rng.integers(0, len(ci), len(nc))]
        m2[nc], f2[nc] = m[src], f[src]
    D_l, D_y, info = TS._designs(m2, f2, cal)
    wts = np.where(cal, (n - 1).astype(float), 0.0)
    info["psi"] = psi
    for key, stat, var, Dd in (("y", S["ybar"], S["vy"], D_y), ("l", S["l"], S["vl"], D_l)):
        t, _ = TS._wls(Dd, np.nan_to_num(stat), wts, cal)
        r = stat - t
        tau2 = max(1e-4, float(np.var(r[cal]) - psi * np.mean(var[cal])))
        w = tau2 / (tau2 + psi * var)
        post = t + w * r
        if key == "y":
            mu = np.where(use & np.isfinite(post), S["c_t"] + post, mu)
        else:
            sd = np.where(use & np.isfinite(post), np.exp(np.where(use, post, 0.0)), sd)
    return mu, sd, info


def adapt_project_mb_w(Xm, Vref, usable, base, mu, sd, Ev, Wev, lam_v=10.0, lam=5.0, extra=3, iters=10, seed=0):
    """nfbase.adapt_project_mb with a weighted U-step over observed + evidence entries (V-step on observed only)."""
    obs = ~np.isnan(Xm)
    ok = np.isfinite(np.nanmean(Xm, 0))
    ev = np.isfinite(Ev) & (Wev > 0) & ~obs & ok[None, :]
    if not ev.any():
        return adapt_project_mb(Xm, Vref, usable, base, mu, sd, lam_v=lam_v, lam=lam, extra=extra, iters=iters,
                                seed=seed)
    sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
    Z = np.where(obs, (Xm - np.nan_to_num(mu)) / sd, 0.0)
    Ze = np.where(ev, (np.nan_to_num(Ev) - np.nan_to_num(mu)) / sd, 0.0)
    Zu = np.where(ev, Ze, Z)
    W = np.where(obs, 1.0, np.where(ev, Wev, 0.0))
    r0 = Vref.shape[1]
    rng = np.random.default_rng(seed)
    Vp = np.hstack([np.where(usable[:, None], Vref, 0.0), np.zeros((len(Vref), extra))])
    V = Vp + np.hstack([np.zeros_like(Vref), 0.01 * rng.normal(size=(len(Vref), extra))])
    r = r0 + extra
    U = np.zeros((Xm.shape[0], r))
    I = np.eye(r)
    for _ in range(iters):
        for i in range(Xm.shape[0]):
            o = (obs[i] & ok) | ev[i]
            Vo = V[o]
            Vw = Vo * W[i, o][:, None]
            U[i] = np.linalg.solve(Vw.T @ Vo + lam * I, Vw.T @ Zu[i, o])
        for j in np.where(ok)[0]:
            o = obs[:, j]
            Uo = U[o]
            V[j] = np.linalg.solve(Uo.T @ Uo + lam_v * I, Uo.T @ Z[o, j] + lam_v * Vp[j])
    out = base.copy()
    Pm = U @ V.T
    miss = ~obs & ok[None, :]
    out[miss] = (mu[None, :] + sd[None, :] * Pm)[miss]
    return out


def _fit(M, ctx, psi, kappa, P, D, kind, iters=10):
    if kappa == 0 or P is None:
        Ev = Wev = None
        S = TS.target_stats(M)
        mu, sd, _ = shrink_w(M, psi, ctx["seed"], S)
        return adapt_project_mb(M, ctx["V"], ctx["warm"], ctx["base"], mu, sd, lam_v=LAM_V, extra=EXTRA, iters=iters)
    Ev = _pseudo(M, P, D, kind)
    Wev = np.where(np.isfinite(P), kappa * P, 0.0)
    S = target_stats_w(M, Ev, Wev)
    mu, sd, _ = shrink_w(M, psi, ctx["seed"], S)
    return adapt_project_mb_w(M, ctx["V"], ctx["warm"], ctx["base"], mu, sd, Ev, Wev, lam_v=LAM_V, extra=EXTRA,
                              iters=iters)


def _log(variant, rec):
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / f"log_{variant}.jsonl", "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def _run(Xm, ctx, variant="main", null=False, kappa_fixed=None):
    al = _aligned(Xm, ctx)
    ok = np.isfinite(np.nanmean(Xm, 0))
    if al is None:  # no evidence cache: the base pipeline (= tgtshrink.impute_null) unchanged
        P = D = None
        kind = "none"
        kappa_fixed = 0.0
    else:
        P, D = al
        kind = "quant" if np.isfinite(D).any() else "detect"
    if null and P is not None:
        P, D = _permute_records(P, D, ok, ctx["seed"])
    seed = ctx["seed"]
    Msel, (r1, c1), (r2, c2) = TS.selection_mask(Xm, seed)

    def score(Y):
        e1 = np.sqrt(np.mean((Y[r1, c1] - Xm[r1, c1]) ** 2))
        e2 = np.sqrt(np.mean((Y[r2, c2] - Xm[r2, c2]) ** 2))
        return 0.5 * (e1 + e2)

    s_psi = {p: score(_fit(Msel, ctx, p, 0.0, None, None, kind, iters=5)) for p in TS.PSI_GRID}
    psi = min(TS.PSI_GRID, key=lambda p: (s_psi[p], p))
    s_k = {}
    if kappa_fixed is not None:
        kappa = kappa_fixed
    else:
        s_k = {0.0: s_psi[psi]}
        for k in KAPPA_GRID[1:]:
            s_k[k] = score(_fit(Msel, ctx, psi, k, P, D, kind, iters=5))
        kappa = min(KAPPA_GRID, key=lambda k: (s_k[k], k))
    Y = _fit(Xm, ctx, psi, kappa, P, D, kind)
    evm = (np.isfinite(P) & ok[None, :] & (P > 0)) if P is not None else np.zeros(Xm.shape, bool)
    P = P if P is not None else np.zeros(Xm.shape)
    _log(variant, {"target": ctx["target"], "seed": seed, "n_missing": int(np.isnan(Xm).sum()), "psi": psi,
                   "kappa": kappa, "scores_psi": {str(k): float(v) for k, v in s_psi.items()},
                   "scores_kappa": {str(k): float(v) for k, v in s_k.items()}, "n_evidence": int(evm.sum()),
                   "sum_P": float(np.nansum(np.where(evm, P, 0.0))), "sum_pi": float(kappa * np.nansum(np.where(evm, P, 0.0)))})
    return Y


def impute(Xm, ctx):
    return _run(Xm, ctx, "main")


def impute_null(Xm, ctx):
    return _run(Xm, ctx, "null", null=True)


def impute_base(Xm, ctx):
    return _run(Xm, ctx, "base", kappa_fixed=0.0)


def impute_k1(Xm, ctx):
    return _run(Xm, ctx, "k1", kappa_fixed=1.0)


def impute_k1_null(Xm, ctx):
    return _run(Xm, ctx, "k1_null", null=True, kappa_fixed=1.0)


# ======================================================================================================================
# diagnostics / unit check
# ======================================================================================================================
def diag():
    import phasec as H

    from explore_phaseb import load_target, subsample
    from run import make_mask
    from scp_graphfs_mcp import phase2 as P2

    res = {"unit_check_kappa0_vs_tgtshrink_impute_null": {}, "weighted_code_pi0_vs_unweighted": {},
           "evidence_per_split": {}}
    for t in ("leduc2022", "sadiku2025"):
        X, names, labels = load_target(t)
        V, warm = H.target_models(t, names)
        for kind in ("MCAR", "MNAR"):
            for seed in (100, 101):
                idx = subsample(X.shape[0], labels, 30, seed)
                Xs = X[idx]
                elig = ((~np.isnan(Xs)).mean(0) >= 0.5) & ~np.isnan(Xs).all(0)
                mask = make_mask(Xs, elig, kind, seed)
                Xm = Xs.copy()
                Xm[mask] = np.nan
                ctx = {"names": names, "target": t, "seed": seed, "V": V, "warm": warm}
                P, D = _aligned(Xm, ctx)
                ok = np.isfinite(np.nanmean(Xm, 0))
                evm = np.isfinite(P) & (P > 0) & ok[None, :]
                res["evidence_per_split"][f"{t}_{kind}_{seed}"] = {
                    "n_evidence_P_gt0": int(evm.sum()), "n_missing_usable": int((np.isnan(Xm) & ok[None, :]).sum()),
                    "n_natural_missing_usable": int((np.isnan(Xs) & ok[None, :]).sum()),
                    "sum_P": float(P[evm].sum()), "cells_with_evidence": int(evm.any(1).sum())}
                if kind == "MNAR" and seed == 100:
                    b = P2._fill(Xm, P2.loo_knn(Xm))
                    b = np.where(np.isnan(b), np.nanmean(Xm), b)
                    ctx["base"] = b
                    _orig_log = TS._log
                    TS._log = lambda *a, **k: None  # do not write into tgtshrink's own log
                    try:
                        ref = TS.impute_null(Xm.copy(), ctx)
                    finally:
                        TS._log = _orig_log
                    mine = impute_base(Xm.copy(), ctx)
                    dmax = float(np.max(np.abs(ref - mine)))
                    res["unit_check_kappa0_vs_tgtshrink_impute_null"][t] = {"max_abs_diff": dmax,
                                                                            "pass": bool(dmax == 0.0)}
                    # weighted code paths with all evidence weights = 0 vs unweighted (stats + projection)
                    kindq = "quant" if np.isfinite(D).any() else "detect"
                    Ev = _pseudo(Xm, P, D, kindq)
                    Wz = np.zeros_like(Xm)
                    S0, S1 = TS.target_stats(Xm), target_stats_w(Xm, Ev, Wz)
                    d_stats = float(np.nanmax(np.abs(S0["xbar"] - S1["xbar"])))
                    # weighted projection with tiny weights -> close to unweighted
                    mu, sd, _ = shrink_w(Xm, 4.0, seed, S0)
                    a0 = adapt_project_mb(Xm, V, warm, b, mu, sd, lam_v=LAM_V, extra=EXTRA)
                    a1 = adapt_project_mb_w(Xm, V, warm, b, mu, sd, Ev, np.where(np.isfinite(P), 1e-12, 0.0),
                                            lam_v=LAM_V, extra=EXTRA)
                    res["weighted_code_pi0_vs_unweighted"][t] = {"stats_xbar_max_abs_diff": d_stats,
                                                                 "projection_pi_1e-12_max_abs_diff":
                                                                     float(np.max(np.abs(a0 - a1)))}
                print(t, kind, seed, res["evidence_per_split"][f"{t}_{kind}_{seed}"], flush=True)
    print(json.dumps(res["unit_check_kappa0_vs_tgtshrink_impute_null"]), json.dumps(res["weighted_code_pi0_vs_unweighted"]))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "subthr_diag.json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "build":
        build()
    elif sys.argv[1] == "build_leduc":
        build_leduc2022()
    elif sys.argv[1] == "build_sadiku":
        build_sadiku2025()
    elif sys.argv[1] == "check_gregoire":
        print(json.dumps(check_gregoire2023_mixCTRL(), indent=1, default=str))
    elif sys.argv[1] == "diag":
        diag()
