"""Phase-C idea `ttprior` (queen: biological priors).

Tissue/technology-kernel-weighted reference loadings as the adapt_project prior.
  1. META: a priori metadata (species, sample_type, lineage set, acquisition, instrument) for every reference dataset,
     each field sourced (PRIDE API JSON, scpdata manual, Cellosaurus, data/PROVENANCE.md); 'unknown' when unstated.
  2. A kernel distance between the target's metadata and each reference dataset's (FULL = tissue + tech, or one axis).
  3. Row weights w_d = exp(-dist/tau), renormalised to mean 1 over the kept reference rows.
  4. A weighted V-step in a (vectorised) copy of transfer.fit_loadings, warm-started from the champion's V, 6 ALS iters.
  5. impute = transfer.adapt_project(Xm, Vw, warm, base, lam_v=100, extra=10).
Null: same machinery, but the reference datasets' metadata rows are deranged (no dataset keeps its own row); the
target's row stays true.
Cache: data/cache/phasec/ttprior/<target>_<variant>[_null].npz (written atomically). Precompute:
    python -m scp_graphfs_mcp.phasec.ttprior --precompute
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

from scp_graphfs_mcp import transfer as T

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "phasec" / "ttprior"
RANK, LAM, MIN_OBS, ITERS = 20, 5.0, 20, 15
LAM_V, EXTRA = 100.0, 10

U_ = "unknown"
# fields: species, sample_type, lineage (set), acq, instrument, source (short citation; full table in the idea doc)
META = {
    "petrosius2023": dict(species="human", sample_type="cell_line", lineage={"myeloid"}, acq="DIA",
                          instrument="Orbitrap_Astral",
                          source="scpdata manual petrosius2023_AstralAML (OCI-AML8227 FACS; Astral; directDIA "
                                 "Spectronaut 17); 'cultured AML 8227 cells' per scpdata schoof2021"),
    "leduc2022": dict(species="human", sample_type=U_, lineage={"melanocytic", "myeloid"}, acq="TMT_DDA",
                      instrument="Orbitrap_QExactive",
                      source="scpdata manual leduc2022_pSCoPE (melanoma cells vs monocytes; TMT-18 + carrier; "
                             "prioritized DDA; Q-Exactive)"),
    "specht2021": dict(species="human", sample_type=U_, lineage={"myeloid"}, acq="TMT_DDA",
                       instrument="Orbitrap_QExactive",
                       source="scpdata manual specht2019v3 (macrophages vs monocytes; SCoPE2 TMT-11/16 carrier; "
                              "Q-Exactive)"),
    "schoof2021": dict(species="human", sample_type="cell_line", lineage={"myeloid"}, acq="TMT_DDA", instrument=U_,
                       source="scpdata manual schoof2021 (cultured OCI-AML8227; TMTpro + 200-cell booster; "
                              "instrument model not stated)"),
    "ai2025a": dict(species="human", sample_type=U_, lineage={"cardiac_muscle"}, acq="DIA", instrument=U_,
                    source="scpdata manual ai2025a (adult-heart cardiomyocytes; DIA-NN report); PROVENANCE.md "
                           "(heart-region runs)"),
    "amrute2026_MFS": dict(species="human", sample_type=U_, lineage={"vascular_mesenchymal"}, acq="DIA",
                           instrument="Orbitrap_Fusion_Lumos/Eclipse",
                           source="PRIDE PXD078860 (human MFS aneurysm aortic tissue; label-free; sample protocol: "
                                  "Orbitrap Eclipse, DIA; project instrument field lists Orbitrap Ascend)"),
    "bubis2025_A549": dict(species="human", sample_type="cell_line", lineage={"epithelial_carcinoma", "lung"},
                           acq="DIA", instrument="Orbitrap_Astral",
                           source="PRIDE PXD064564 (A549; Orbitrap Astral; label-free; directDIA+ Spectronaut); "
                                  "Cellosaurus CVCL_0023 (lung adenocarcinoma)"),
    "cong2020AC": dict(species="human", sample_type="cell_line", lineage={"epithelial_carcinoma"}, acq="LFQ_DDA",
                       instrument="Orbitrap_Fusion_Lumos/Eclipse",
                       source="scpdata manual cong2020AC (HeLa; label-free; Orbitrap Fusion Eclipse; MaxQuant "
                              "1.6.3.3); Cellosaurus CVCL_0030"),
    "gregoire2023_mixCTRL": dict(species="human", sample_type="cell_line", lineage={"myeloid"}, acq="TMT_DDA",
                                 instrument="Orbitrap_Fusion_Lumos/Eclipse",
                                 source="scpdata manual gregoire2023_mixCTRL (THP1/U937 + macrophage-like; SCoPE2 "
                                        "TMT-16; Lumos mainly, 6/7 corpus batches CBIO=Lumos); Cellosaurus "
                                        "CVCL_0006/CVCL_0007"),
    "guise2024": dict(species="human", sample_type="tissue_LCM", lineage={"neural"}, acq="LFQ_DDA",
                      instrument="Orbitrap_Exploris",
                      source="scpdata manual guise2024 (LCM spinal motor neurons; nanoPOTS; Exploris 480; MS2 "
                             "dynamic exclusion 90 s)"),
    "hu2023_K562": dict(species="human", sample_type="cell_line", lineage={"myeloid"}, acq=U_,
                        instrument="Orbitrap_Fusion_Lumos/Eclipse",
                        source="scpdata manual hu2023_K562 (K562; label-free; Orbitrap Eclipse; acquisition mode not "
                               "stated); Cellosaurus CVCL_0004 (CML)"),
    "inns2026_CCS": dict(species="human", sample_type="tissue_dissociated", lineage={"skin_adnexal", "myeloid"},
                         acq="DIA", instrument="timsTOF",
                         source="PRIDE PXD073250 (dissociated CYLD skin tumours, FACS CD45/CD74/CD200; parts incl. "
                                "keratinocyte/macrophage/DC/monocyte; timsTOF HT diaPASEF; DIA-NN)"),
    "krull2024": dict(species="human", sample_type="cell_line", lineage={"other"}, acq="DIA", instrument="timsTOF",
                      source="scpdata manual krull2024 (U-2 OS; diaPASEF timsTOF Pro); Cellosaurus CVCL_0042 "
                             "(osteosarcoma)"),
    "liang2020_hela": dict(species="human", sample_type="cell_line", lineage={"epithelial_carcinoma"},
                           acq="LFQ_DDA", instrument="Orbitrap_Exploris",
                           source="scpdata manual liang2020_hela (HeLa; autoPOTS; Exploris 480; exclusion duration "
                                  "90 s; MaxQuant 1.6.7.0); Cellosaurus CVCL_0030"),
    "px063590_ADbrain": dict(species="human", sample_type="tissue_dissociated", lineage={"neural"}, acq="TMT_DDA",
                             instrument="Orbitrap_Exploris",
                             source="PRIDE PXD063590 (dissociated DLPFC; nPOP TMTpro + carrier; pSCoPE; Exploris "
                                    "480)"),
    "sadiku2025": dict(species="human", sample_type="tissue_dissociated", lineage={"myeloid"}, acq="DIA",
                       instrument="Orbitrap_Astral",
                       source="PRIDE PXD061065 (neutrophils sorted after brain-tumour dissociation; Astral DIA; "
                              "Spectronaut 19.4)"),
    "weiss2026_liver": dict(species="human", sample_type="tissue_LCM", lineage={"hepatic"}, acq="DIA",
                            instrument="Orbitrap_Astral",
                            source="PRIDE PXD062231 (laser-microdissected hepatocytes; Astral; DIA only; DIA-NN "
                                   "1.8.1)"),
    "williams2020_lfq": dict(species="human", sample_type="cell_line", lineage={"other"}, acq="LFQ_DDA",
                             instrument="Orbitrap_Fusion_Lumos/Eclipse",
                             source="scpdata manual williams2020_lfq (MCF10A; label-free nanoPOTS; Fusion Lumos; "
                                    "MaxQuant 1.6.2.10); Cellosaurus CVCL_0598 (non-cancer, immortalized)"),
    "williams2020_tmt": dict(species="human", sample_type="cell_line", lineage={"myeloid"}, acq="TMT_DDA",
                             instrument="Orbitrap_Fusion_Lumos/Eclipse",
                             source="scpdata manual williams2020_tmt (MOLM-14/K562/CMK 'AML cell lines'; TMT + 10 ng "
                                    "carrier; Fusion Lumos)"),
    "woo2022_lung": dict(species="human", sample_type="tissue_dissociated", lineage={"lung"}, acq="LFQ_DDA",
                         instrument="Orbitrap_Fusion_Lumos/Eclipse",
                         source="scpdata manual woo2022_lung (dissociated primary lung cells; label-free; Fusion "
                                "Lumos FAIMS; MaxQuant 1.6.2.10)"),
    "woo2025_PC3": dict(species="human", sample_type="cell_line", lineage={"epithelial_carcinoma"}, acq="DIA",
                        instrument=U_,
                        source="PRIDE PXD059079 (PC-3 prostate adenocarcinoma; LFQ DIA; two instruments listed, "
                               "Orbitrap Ascend and timsTOF HT); Cellosaurus CVCL_0035"),
    "ye2025_EB": dict(species="human", sample_type="cell_line", lineage={"pluripotent_embryonic"}, acq="DIA",
                      instrument="Orbitrap_Astral",
                      source="PRIDE PXD049181 (iPSC-derived embryoid bodies, 'Cell culture'; Astral; directDIA+ "
                             "Spectronaut 18)"),
}

# "label-multiplexed vs label-free" and "DDA vs DIA" axes for the acquisition term
_ACQ_AX = {"TMT_DDA": ("label", "DDA"), "LFQ_DDA": ("free", "DDA"), "DIA": ("free", "DIA"),
           "plexDIA": ("label", "DIA")}


def _acq_term(a, b):
    if U_ in (a, b):
        return 0.5
    if a == b:
        return 0.0
    shared = sum(x == y for x, y in zip(_ACQ_AX[a], _ACQ_AX[b]))
    return 0.5 if shared == 1 else 1.0


def tissue_dist(m, t):
    if U_ in m["lineage"] or U_ in t["lineage"] or not m["lineage"] or not t["lineage"]:
        lin = 0.5
    else:
        lin = 1.0 - len(m["lineage"] & t["lineage"]) / len(m["lineage"] | t["lineage"])
    st = 0.25 if U_ in (m["sample_type"], t["sample_type"]) else 0.5 * (m["sample_type"] != t["sample_type"])
    return lin + st


def tech_dist(m, t):
    ins = 0.25 if U_ in (m["instrument"], t["instrument"]) else 0.5 * (m["instrument"] != t["instrument"])
    return _acq_term(m["acq"], t["acq"]) + ins


def dist(m, t, kernel):
    if kernel == "full":
        return tissue_dist(m, t) + tech_dist(m, t)
    if kernel == "tissue":
        return tissue_dist(m, t)
    if kernel == "tech":
        return tech_dist(m, t)
    raise ValueError(kernel)


def null_meta_map(target, datasets):
    """Fixed derangement of the reference datasets' metadata rows (rng 12345 + index of target in sorted(META))."""
    ds = sorted(datasets)
    rng = np.random.default_rng(12345 + sorted(META).index(target))
    while True:
        p = rng.permutation(len(ds))
        if not np.any(p == np.arange(len(ds))):
            break
    return {d: ds[p[k]] for k, d in enumerate(ds)}  # dataset d is given the metadata row of ds[p[k]]


def dataset_weights(target, datasets, kernel, tau, null=False):
    t = META[target]
    mp = null_meta_map(target, datasets) if null else {d: d for d in datasets}
    return {d: float(np.exp(-dist(META[mp[d]], t, kernel) / tau)) for d in datasets}


# ---------------------------------------------------------------------------------------------------------------
def fit_loadings_weighted(Zref, w_rows, V0, warm0, rank=RANK, lam=LAM, iters=ITERS, min_obs=MIN_OBS, seed=0,
                          warm_start=False):
    """Vectorised copy of transfer.fit_loadings (kg=None) with per-row weights in the V-step only. `warm` is computed
    unweighted exactly as in the original. Default: COLD start with the original's init (rng seed 0) and 15 iterations,
    so w = 1 reproduces the champion V exactly. warm_start=True starts from V0[warm] (the spec's speed-up; it fails the
    w = 1 unit check because the champion's 15-iteration ALS is not converged: see the idea doc).
    w_rows: weight per row of Zref (renormalised here to mean 1 over the kept rows)."""
    obs = ~np.isnan(Zref)
    warm = obs.sum(0) >= min_obs
    assert np.array_equal(warm, warm0), "warm mask differs from the champion's"
    Z = np.where(obs, Zref, 0.0)[:, warm]
    O = obs[:, warm]
    keep = O.sum(1) >= rank
    Z, O, w = Z[keep], O[keep].astype(float), np.asarray(w_rows, float)[keep]
    w = w * (len(w) / w.sum())
    n, p = Z.shape
    V = V0[warm].copy() if warm_start else 0.1 * np.random.default_rng(seed).normal(size=(p, rank))
    I = lam * np.eye(rank)
    for _ in range(iters):
        VV = (V[:, :, None] * V[:, None, :]).reshape(p, rank * rank)
        A = (O @ VV).reshape(n, rank, rank) + I
        U = np.linalg.solve(A, (Z @ V)[:, :, None])[:, :, 0]
        UU = (U[:, :, None] * U[:, None, :]).reshape(n, rank * rank)
        Ow = O * w[:, None]
        A = (Ow.T @ UU).reshape(p, rank, rank) + I
        V = np.linalg.solve(A, ((Z * w[:, None]).T @ U)[:, :, None])[:, :, 0]
    Vfull = np.zeros((Zref.shape[1], rank))
    Vfull[warm] = V
    return Vfull, warm, keep


VARIANTS = {"full_t2": ("full", 2.0), "full_t1": ("full", 1.0), "full_t05": ("full", 0.5),
            "tissue_t1": ("tissue", 1.0), "tech_t1": ("tech", 1.0),
            "tissue_t05": ("tissue", 0.5), "tech_t05": ("tech", 0.5)}


def _diag(target, datasets, src_kept, kernel, tau, null):
    wd = dataset_weights(target, datasets, kernel, tau, null)
    n = {d: int((src_kept == d).sum()) for d in datasets}
    tot = sum(n[d] * wd[d] for d in datasets)
    neff = tot ** 2 / sum(n[d] ** 2 * wd[d] ** 2 for d in datasets)
    share = {d: n[d] * wd[d] / tot for d in datasets}
    top = sorted(share, key=share.get, reverse=True)[:5]
    return {"w_dataset": wd, "n_rows": n, "neff_datasets": neff, "neff_unweighted":
            sum(n.values()) ** 2 / sum(v ** 2 for v in n.values()), "top5_share": [(d, round(share[d], 3)) for d in top]}


def get_Vw(ctx, variant, null=False):
    target = ctx["target"]
    kernel, tau = VARIANTS[variant]
    f = CACHE / f"{target}_{variant}{'_null' if null else ''}.npz"
    if f.exists():
        return np.load(f)["V"]
    Z = ctx["ref"]()
    src = np.asarray(ctx["ref_src"]())
    datasets = sorted(set(src.tolist()))
    wd = dataset_weights(target, datasets, kernel, tau, null)
    w_rows = np.array([wd[s] for s in src])
    Vw, warm, keep = fit_loadings_weighted(Z, w_rows, ctx["V"], ctx["warm"])
    CACHE.mkdir(parents=True, exist_ok=True)
    tmp = f.with_suffix(f".tmp{os.getpid()}.npz")
    np.savez_compressed(tmp, V=Vw)
    os.replace(tmp, f)
    dg = _diag(target, datasets, src[keep], kernel, tau, null)
    dfile = CACHE / f"{target}_{variant}{'_null' if null else ''}_diag.json"
    tmpj = dfile.with_suffix(f".tmp{os.getpid()}")
    tmpj.write_text(json.dumps(dg, indent=1))
    os.replace(tmpj, dfile)
    return Vw


def _impute(Xm, ctx, variant, null=False):
    Vw = get_Vw(ctx, variant, null)
    return T.adapt_project(Xm, Vw, ctx["warm"], ctx["base"], lam_v=LAM_V, extra=EXTRA)


def _twoblock(Xm, ctx, null=False):
    Vref = np.hstack([get_Vw(ctx, "tissue_t05", null), get_Vw(ctx, "tech_t05", null)])
    return T.adapt_project(Xm, Vref, ctx["warm"], ctx["base"], lam_v=LAM_V, extra=EXTRA)


def _mk(variant):
    return (lambda Xm, ctx: _impute(Xm, ctx, variant, False)), (lambda Xm, ctx: _impute(Xm, ctx, variant, True))


impute_full_t2, impute_null_full_t2 = _mk("full_t2")
impute_full_t1, impute_null_full_t1 = _mk("full_t1")
impute_full_t05, impute_null_full_t05 = _mk("full_t05")
impute_tissue_t1, impute_null_tissue_t1 = _mk("tissue_t1")
impute_tech_t1, impute_null_tech_t1 = _mk("tech_t1")
impute = impute_full_t1          # primary variant (screened as id `ttprior`)
impute_null = impute_null_full_t1


def impute_twoblock(Xm, ctx):
    return _twoblock(Xm, ctx, False)


def impute_null_twoblock(Xm, ctx):
    return _twoblock(Xm, ctx, True)


# ---------------------------------------------------------------------------------------------------------------
def _bench_ctx(target):
    sys.path.insert(0, str(ROOT / "benchmarks"))
    import phasec as H
    from explore_phaseb import load_target

    _, names, _ = load_target(target)
    V, warm = H.target_models(target, names)
    Z, src = H._ref_parts(target, names)
    return {"target": target, "names": names, "V": V, "warm": warm, "ref": lambda: Z, "ref_src": lambda: src}


def unit_check(target="petrosius2023"):
    ctx = _bench_ctx(target)
    Z = ctx["ref"]()
    out = {}
    for ws, it in ((False, ITERS), (True, 6)):
        Vw, warm, _ = fit_loadings_weighted(Z, np.ones(len(Z)), ctx["V"], ctx["warm"], iters=it, warm_start=ws)
        out[f"warm_start={ws},iters={it}"] = float(np.linalg.norm(Vw - ctx["V"]) / np.linalg.norm(ctx["V"]))
    return out


def precompute(targets, variants):
    import time
    from scipy.stats import spearmanr

    out = {}
    for t in targets:
        ctx = _bench_ctx(t)
        src = np.asarray(ctx["ref_src"]())
        datasets = sorted(set(src.tolist()))
        for v in variants:
            for nl in (False, True):
                t0 = time.time()
                get_Vw(ctx, v, nl)
                print(f"{t} {v} null={nl} {time.time() - t0:.1f}s", flush=True)
            k, tau = VARIANTS[v]
            a = dataset_weights(t, datasets, k, tau, False)
            b = dataset_weights(t, datasets, k, tau, True)
            rho = spearmanr([a[d] for d in datasets], [b[d] for d in datasets]).statistic
            out.setdefault(t, {})[v] = {"spearman_true_vs_null": None if np.isnan(rho) else float(rho)}
        out[t]["null_map"] = null_meta_map(t, datasets)
        (CACHE / f"{t}_spearman.json").write_text(json.dumps(out[t], indent=1))


if __name__ == "__main__":
    FAST = ["petrosius2023", "leduc2022", "weiss2026_liver", "gregoire2023_mixCTRL", "px063590_ADbrain", "sadiku2025"]
    if "--unit" in sys.argv:
        res = unit_check()
        print("rel Frobenius |Vw(w=1) - V| / |V| =", res)
        CACHE.mkdir(parents=True, exist_ok=True)
        (CACHE / "unit_check.json").write_text(json.dumps(res, indent=1))
    if "--precompute" in sys.argv:
        vs = sys.argv[sys.argv.index("--variants") + 1].split(",") if "--variants" in sys.argv else list(VARIANTS)
        tg = sys.argv[sys.argv.index("--targets") + 1].split(",") if "--targets" in sys.argv else FAST
        precompute(tg, vs)
