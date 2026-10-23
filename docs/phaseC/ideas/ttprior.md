# ttprior: tissue/technology-kernel-weighted reference loadings as the adapt_project prior

- **Queen:** biological priors (round-1 directive in `KING_LOG.md`)
- **Status (worker's recommendation):** **dead** at the fast screen.
  - Best variant: FULL kernel, tau = 1 (`ttprior`). cand_stack is 0.9951 (all) / 0.9953 (in-complex). That misses the
    0.99 bar.
  - Margin over its null: 0.0040 (all) / 0.0052 (in-complex). The all-protein margin is below the 0.005 bar, though the
    sign agrees on both subsets.
  - None of the six screened variants reached 0.99, so the full pool was not run. No v2 revision was made: there was no
    bug or scale problem to fix (see "Honest read").
- **Code:** `src/scp_graphfs_mcp/phasec/ttprior.py`.
- **Results:** `results/phaseC/ttprior_fast.json` (FULL tau=1), `ttprior_full_t2_fast.json`,
  `ttprior_full_t05_fast.json`, `ttprior_tissue_t1_fast.json`, `ttprior_tech_t1_fast.json` and
  `ttprior_twoblock_fast.json`.
- **Diagnostics** (weights, effective number of datasets, null maps, Spearman): `data/cache/phasec/ttprior/diagnostics.json`
  and `unit_check.json`. These are cache files and are not committed.

## Idea
The champion's prior V is `fit_loadings` on every leave-target-out reference row, so each dataset counts in proportion to
its cell count. Four large sets dominate it: px063590_ADbrain, schoof2021, leduc2022 and specht2021, three of them TMT
carrier DDA. This idea re-fits V with **a priori metadata weights** instead,
w_d = exp(−dist(meta_d, meta_target)/tau), where the distance is built from lineage, sample type, acquisition and
instrument. The re-fitted V is then used as the `adapt_project` prior. It amounts to covariate-shift importance weighting
at the dataset level (Shimodaira 2000). As tau → ∞ it falls back to the champion.

## Mechanism (as built)
1. **Metadata table.** `META` in the module has one row per reference dataset. The full table, with sources, is below.
2. **Distance.**
   - tissue = (1 − Jaccard(lineage)) + 0.5·[sample_type differs].
   - tech = acq_term + 0.5·[instrument differs].
     - acq_term is 0 when the acquisition is equal.
     - It is 0.5 when exactly one of the axes (label-multiplexed vs label-free) and (DDA vs DIA) is shared.
     - It is 1 otherwise.
   - 'unknown' on either side contributes half the field maximum: lineage 0.5, sample_type 0.25, acq 0.5,
     instrument 0.25.
   - Kernels: FULL = tissue + tech, TISSUE = tissue only, TECH = tech only.
3. **Weights.** Each row gets its dataset's w_d. Weights are renormalised to mean 1 over the rows that `fit_loadings`
   keeps, so the ridge lam = 5 stays on the champion's scale.
4. **Weighted `fit_loadings`.**
   - This is a vectorised copy of `transfer.fit_loadings` (kg = None): rank 20, lam 5, min_obs 20.
   - `warm` and `keep` are computed **unweighted**, exactly as in the original. The code asserts that `warm` equals the
     champion's.
   - The U-step is unchanged. The V-step solves A = Uoᵀdiag(w)Uo + 5I, b = Uoᵀdiag(w)Z.
5. **Candidate.** `transfer.adapt_project(Xm, Vw, warm, base, lam_v=100, extra=10)`.
   - Vw does not depend on Xm, so it is cached per (target, variant, null) in
     `data/cache/phasec/ttprior/{target}_{variant}[_null].npz`, written with tmp + `os.replace`.
   - `python -m scp_graphfs_mcp.phasec.ttprior --precompute [--targets ...] [--variants ...]` builds the caches.
6. **Two-block variant.** Vref = [V_TISSUE,tau=0.5 | V_TECH,tau=0.5] (p × 40), passed to the same `adapt_project`
   call. Its null concatenates the two deranged-metadata fits.

### Deviation from the spec: a cold start instead of a warm start (unit check)
The spec asked for a warm start from the champion V plus 6 ALS iterations, *and* required that w = 1 reproduce
`ctx['V']` to a relative Frobenius error below 1e-2. The two requirements are incompatible. The champion's 15-iteration
ALS is not converged, so continuing it moves V.

Measured on petrosius2023 (`data/cache/phasec/ttprior/unit_check.json`):

| Start | w = 1 vs champion V (relative Frobenius error) |
|---|---|
| Warm start, 6 iterations | 0.1243 (fails) |
| Cold start: original init, rng seed 0, 15 iterations | 3.8e-15 (passes) |

Before the switch I checked the warm-start runs further:
- The vectorised step matches a loop copy of the original to 4.5e-15 after 1 iteration.
- Drift relative to the champion V grows with the number of extra iterations: 0.025 after 1, 0.048 after 2 and 0.124
  after 6. The relative change in VVᵀ is similar, so the drift is not just a rotation.

I therefore used the cold start (15 iterations, same init as the original). With it, tau = ∞ *is* the champion
exactly, and any difference comes from the weights alone. Each fit took about 20–46 s.

## Null control
- The same kernel, tau, renormalisation, weighted ALS and `adapt_project` call are used.
- The only change: the reference datasets' metadata rows are permuted by a fixed **derangement**.
  - The rng is `default_rng(12345 + index of target in sorted(META))`, and the derangement is checked to have no fixed
    point.
  - The target's own metadata row stays true.
- Spearman correlations between the true and null per-dataset weight vectors are listed below. They show the null
  really decorrelates the weights. They are the same for every tau of a given kernel, because the weights are monotone
  in the distance.

| Target | FULL | TISSUE | TECH |
|---|---|---|---|
| petrosius2023 | 0.21 | 0.09 | −0.26 |
| leduc2022 | 0.18 | 0.54 | 0.24 |
| weiss2026_liver | 0.08 | −0.07 | 0.14 |
| gregoire2023_mixCTRL | 0.34 | 0.01 | 0.29 |
| px063590_ADbrain | 0.06 | −0.17 | 0.18 |
| sadiku2025 | −0.25 | −0.17 | 0.09 |

(leduc2022 TISSUE at 0.54 is the most preserving permutation. Its TISSUE weights are flat anyway, because almost every
dataset is "unknown" or disjoint in lineage with it.)

## Weight diagnostics (FULL kernel, tau = 1)
- **Effective number of datasets** = (Σ n_d w_d)² / Σ n_d² w_d².
- **Top-5 share** is each dataset's share of total row weight.

| Target | n_eff (unweighted → weighted → null) | Top-5 share (weighted) |
|---|---|---|
| petrosius2023 | 7.10 → 6.52 → 6.90 | schoof2021 0.324, specht2021 0.120, sadiku2025 0.088, amrute2026_MFS 0.079, weiss2026_liver 0.079 |
| leduc2022 | 6.33 → 3.84 → 3.25 | schoof2021 0.371, specht2021 0.290, px063590_ADbrain 0.189, amrute2026_MFS 0.026, weiss2026_liver 0.019 |
| weiss2026_liver | 6.63 → 9.42 → 5.20 | amrute2026_MFS 0.171, schoof2021 0.157, px063590_ADbrain 0.132, leduc2022 0.099, specht2021 0.095 |
| gregoire2023_mixCTRL | 7.20 → 3.47 → 3.53 | schoof2021 0.483, specht2021 0.178, leduc2022 0.112, px063590_ADbrain 0.091, williams2020_tmt 0.039 |
| px063590_ADbrain | 6.93 → 4.77 → 4.04 | schoof2021 0.339, leduc2022 0.214, specht2021 0.206, amrute2026_MFS 0.050, weiss2026_liver 0.030 |
| sadiku2025 | 7.00 → 8.37 → 5.83 | schoof2021 0.225, specht2021 0.137, px063590_ADbrain 0.114, amrute2026_MFS 0.090, weiss2026_liver 0.090 |

The per-variant and per-null top-5 lists are in `diagnostics.json`.
- For DIA targets the weighting *spreads* weight away from the TMT giants, so n_eff rises (weiss2026_liver, sadiku2025).
- For TMT/myeloid targets it *concentrates* weight on schoof2021, specht2021 and leduc2022.

## Results (fast pool: 6 targets × MCAR/MNAR × seeds 100/101; ratios vs champion, < 1 better)
All values are `mean_rel_to_champion` from `results/phaseC/<id>_fast.json`, given as all / in-complex.

| id | kernel, tau | cand_stack | cand_solo | null_stack | cand_raw (all) | splits ≤ 0.99 | mean cand weight |
|---|---|---|---|---|---|---|---|
| ttprior | FULL, 1 | **0.9951 / 0.9953** | 0.9944 / 0.9947 | 0.9991 / 1.0005 | 0.9939 | 1/12 | 0.510 |
| ttprior_full_t2 | FULL, 2 | 0.9975 / 0.9971 | 0.9965 / 0.9962 | 0.9998 / 1.0002 | 0.9954 | 0/12 | 0.490 |
| ttprior_full_t05 | FULL, 0.5 | 0.9953 / 0.9970 | 0.9961 / 0.9989 | 0.9987 / 0.9994 | 0.9957 | 1/12 | 0.406 |
| ttprior_tissue_t1 | TISSUE, 1 | 0.9992 / 0.9992 | 0.9990 / 0.9987 | 0.9997 / 1.0000 | 0.9990 | 0/12 | 0.354 |
| ttprior_tech_t1 | TECH, 1 | 0.9970 / 0.9969 | 0.9970 / 0.9971 | 0.9980 / 0.9983 | 0.9955 | 0/12 | 0.500 |
| ttprior_twoblock | [TISSUE 0.5, TECH 0.5] | 0.9965 / 0.9962 | 0.9992 / 0.9982 | 0.9984 / 0.9980 | 1.0031 | 1/12 | 0.271 |
| sanity (tau = ∞) | the champion | 1.000 (by construction; `sanity_fast.json`, not re-run) | | | | | |

### Per-target cand_stack | null_stack (rmse_all)

| Target / mask | ttprior (FULL 1) | full_t05 | tech_t1 | tissue_t1 |
|---|---|---|---|---|
| petrosius2023 MCAR | 0.9905 \| 0.9984 | 0.9920 \| 0.9952 | 0.9915 \| 0.9954 | 0.9985 \| 1.0002 |
| petrosius2023 MNAR | 0.9943 \| 0.9955 | 0.9938 \| 0.9953 | 0.9942 \| 0.9974 | 0.9980 \| 0.9982 |
| sadiku2025 MCAR | 0.9797 \| 0.9986 | 0.9894 \| 0.9985 | 0.9910 \| 0.9997 | 0.9998 \| 1.0005 |
| sadiku2025 MNAR | 0.9927 \| 0.9984 | 0.9920 \| 0.9968 | 0.9953 \| 0.9867 | 1.0023 \| 1.0000 |
| weiss2026_liver MCAR | 0.9932 \| 0.9992 | 0.9907 \| 0.9969 | 0.9941 \| 0.9989 | 0.9986 \| 1.0000 |
| weiss2026_liver MNAR | 0.9956 \| 1.0058 | 0.9930 \| 1.0074 | 0.9967 \| 1.0026 | 0.9986 \| 1.0000 |
| gregoire2023_mixCTRL MCAR | 0.9977 \| 0.9979 | 0.9971 \| 0.9984 | 1.0031 \| 0.9999 | 1.0000 \| 0.9997 |
| gregoire2023_mixCTRL MNAR | 0.9967 \| 0.9967 | 0.9969 \| 0.9964 | 0.9978 \| 0.9965 | 0.9968 \| 0.9985 |
| leduc2022 MCAR | 1.0012 \| 1.0002 | 0.9994 \| 1.0009 | 1.0013 \| 1.0002 | 0.9991 \| 1.0001 |
| leduc2022 MNAR | 0.9999 \| 1.0000 | 0.9996 \| 0.9999 | 1.0000 \| 1.0000 | 1.0001 \| 1.0002 |
| px063590_ADbrain MCAR | 0.9992 \| 0.9988 | 0.9997 \| 0.9992 | 0.9989 \| 0.9993 | 0.9986 \| 0.9992 |
| px063590_ADbrain MNAR | 1.0000 \| 0.9996 | 0.9998 \| 0.9999 | 1.0007 \| 0.9998 | 0.9995 \| 1.0000 |

The full_t2 and twoblock per-target rows are in their JSON files.

## Honest read
- **Below the bar.** The best variant (FULL, tau = 1) gives about −0.5% vs the champion, against about −0.1% for its
  null. The direction is right and agrees on all proteins and in-complex. But it is an order of magnitude short of what
  the colony needs, and it fails both promotion bars (≤ 0.99, and a margin ≥ 0.005 on all proteins).
- **The signal is technology, not tissue.**
  - Every TISSUE-only variant is indistinguishable from its null: 0.9992 vs 0.9997.
  - TECH-only captures most of the FULL gain on the DIA targets.
  - The whole gain sits on the **label-free DIA, Orbitrap Astral targets**: sadiku2025 (best split 0.9797),
    petrosius2023 and weiss2026_liver.
  - The TMT-carrier DDA targets (leduc2022, px063590_ADbrain, gregoire2023_mixCTRL) are at 0.997–1.001. For them, the
    weighting mostly re-concentrates the prior on the same TMT giants that already dominate the pooled V, so the prior
    barely changes.
  - The answer to the scientific question is therefore: the transferable co-regulation signal that reweighting can
    exploit tracks the **acquisition platform**. Lineage overlap between public datasets and these targets is too thin
    to matter: most DIA targets have no lineage-matched reference with many cells.
- **The TECH null is close to TECH.** It reaches 0.9980, and one null split hits 0.9867 (sadiku2025 MNAR). So part of
  any reweighting effect is "any departure from cell-count weighting", not the right metadata.
- **cand_raw ≈ cand_solo ≈ cand_stack.** Mean candidate weight in the stack is about 0.5. The weighted prior is a
  slightly better *replacement* for the champion's adapted projection, rather than an additive component.
- **Caveat.** Several metadata fields are 'unknown' (see table). Some of these could be resolved from the papers, but
  given the TISSUE result this is unlikely to change the conclusion.

## Next step (for the queen)
- Do not spend the full pool on this.
- If anything, take the technology finding into a cheaper form: a **DIA-only / TMT-only split prior** (two V blocks, one
  fitted on each acquisition class, each with its own ridge) for DIA targets.
- Better still, treat acquisition as a *noise model*: per-dataset/per-technology protein offsets and scales, the
  Netflix-baseline idea already queued by the abduction queen. The gain here is concentrated where the reference
  platform differs most from the target, which suggests that technology-specific variance, not co-regulation content,
  is what the weighting partly corrects.

## Metadata table (META) and sources
Vocabularies are exactly as in the directive, with no extra values. The only conflict was Orbitrap Ascend: for
amrute2026_MFS it is resolved from the sample protocol, and for woo2025_PC3 the instrument is 'unknown'. Instrument and
sample-type differences enter the distance only through equality checks.

| Dataset | sample_type | lineage | acq | instrument | Source |
|---|---|---|---|---|---|
| petrosius2023 | cell_line | myeloid | DIA | Orbitrap_Astral | scpdata manual `petrosius2023_AstralAML` (FACS-sorted OCI-AML8227; Astral; directDIA Spectronaut 17); "cultured AML 8227 cells" in the scpdata `schoof2021` entry |
| leduc2022 | unknown | melanocytic, myeloid | TMT_DDA | Orbitrap_QExactive | scpdata `leduc2022_pSCoPE` (melanoma cells vs monocytes; TMT-18 + carrier; prioritized DDA; Q-Exactive) |
| specht2021 | unknown | myeloid | TMT_DDA | Orbitrap_QExactive | scpdata `specht2019v3` (macrophages vs monocytes; SCoPE2 TMT-11/16 with carrier; Q-Exactive) |
| schoof2021 | cell_line | myeloid | TMT_DDA | unknown | scpdata `schoof2021` (cultured OCI-AML8227; TMTpro + 200-cell booster; Orbitrap model not stated) |
| ai2025a | unknown | cardiac_muscle | DIA | unknown | scpdata `ai2025a` (adult-heart cardiomyocytes; generated from a DIA-NN report); PROVENANCE.md (heart-region runs) |
| amrute2026_MFS | unknown | vascular_mesenchymal | DIA | Orbitrap_Fusion_Lumos/Eclipse | PRIDE PXD078860 (human MFS aneurysm aortic tissue, label-free). The sample protocol says Orbitrap Eclipse with DIA; the project instrument field says Orbitrap Ascend. I used the protocol |
| bubis2025_A549 | cell_line | epithelial_carcinoma, lung | DIA | Orbitrap_Astral | PRIDE PXD064564 (A549; Orbitrap Astral; LFQ; directDIA+ Spectronaut); Cellosaurus CVCL_0023 (lung adenocarcinoma) |
| cong2020AC | cell_line | epithelial_carcinoma | LFQ_DDA | Orbitrap_Fusion_Lumos/Eclipse | scpdata `cong2020AC` (HeLa; label-free; Orbitrap Fusion Eclipse; MaxQuant 1.6.3.3 †); Cellosaurus CVCL_0030 |
| gregoire2023_mixCTRL | cell_line | myeloid | TMT_DDA | Orbitrap_Fusion_Lumos/Eclipse | scpdata `gregoire2023_mixCTRL` (THP1/U937 + macrophage-like; SCoPE2 TMT-16; "Lumos (mainly)", and 6/7 corpus batches are CBIO = Lumos); Cellosaurus CVCL_0006/0007 |
| guise2024 | tissue_LCM | neural | LFQ_DDA | Orbitrap_Exploris | scpdata `guise2024` (LCM spinal motor neurons; nanoPOTS; Exploris 480; MS2 dynamic-exclusion 90 s ⇒ DDA) |
| hu2023_K562 | cell_line | myeloid | unknown | Orbitrap_Fusion_Lumos/Eclipse | scpdata `hu2023_K562` (K562; label-free; Orbitrap Eclipse; mode not stated); Cellosaurus CVCL_0004 (CML) |
| inns2026_CCS | tissue_dissociated | skin_adnexal, myeloid | DIA | timsTOF | PRIDE PXD073250 (dissociated CYLD skin tumours, FACS CD45/CD74/CD200; keratinocyte/macrophage/DC/monocyte; timsTOF HT diaPASEF) |
| krull2024 | cell_line | other | DIA | timsTOF | scpdata `krull2024` (U-2 OS; diaPASEF on timsTOF Pro); Cellosaurus CVCL_0042 (osteosarcoma, no vocabulary tag) |
| liang2020_hela | cell_line | epithelial_carcinoma | LFQ_DDA | Orbitrap_Exploris | scpdata `liang2020_hela` (HeLa; autoPOTS; Exploris 480; exclusion duration 90 s; MaxQuant 1.6.7.0 †); Cellosaurus CVCL_0030 |
| px063590_ADbrain | tissue_dissociated | neural | TMT_DDA | Orbitrap_Exploris | PRIDE PXD063590 (dissociated DLPFC; nPOP TMTpro with carrier; pSCoPE; Exploris 480) |
| sadiku2025 | tissue_dissociated | myeloid | DIA | Orbitrap_Astral | PRIDE PXD061065 (neutrophils sorted after brain-tumour dissociation; Astral DIA; Spectronaut 19.4) |
| weiss2026_liver | tissue_LCM | hepatic | DIA | Orbitrap_Astral | PRIDE PXD062231 (laser-microdissected hepatocytes; Astral; "exclusively DIA"; DIA-NN 1.8.1) |
| williams2020_lfq | cell_line | other | LFQ_DDA | Orbitrap_Fusion_Lumos/Eclipse | scpdata `williams2020_lfq` (MCF10A; label-free nanoPOTS; Fusion Lumos; MaxQuant 1.6.2.10 †); Cellosaurus CVCL_0598 (non-cancer, so not epithelial_carcinoma) |
| williams2020_tmt | cell_line | myeloid | TMT_DDA | Orbitrap_Fusion_Lumos/Eclipse | scpdata `williams2020_tmt` ("three acute myeloid leukemia cell lines"; TMT + 10 ng carrier; Fusion Lumos) |
| woo2022_lung | tissue_dissociated | lung | LFQ_DDA | Orbitrap_Fusion_Lumos/Eclipse | scpdata `woo2022_lung` (dissociated primary lung cells; label-free; Fusion Lumos + FAIMS; MaxQuant 1.6.2.10 †) |
| woo2025_PC3 | cell_line | epithelial_carcinoma | DIA | unknown | PRIDE PXD059079 (PC-3 prostate adenocarcinoma; LFQ + DIA; lists both Orbitrap Ascend and timsTOF HT); Cellosaurus CVCL_0035 |
| ye2025_EB | cell_line | pluripotent_embryonic | DIA | Orbitrap_Astral | PRIDE PXD049181 (iPSC-derived embryoid bodies, "Cell culture"; Astral; directDIA+ Spectronaut 18) |

All rows are species = human.

† These LFQ rows are marked DDA because they were processed with label-free MaxQuant 1.6.x, which has no DIA
workflow (MaxDIA arrived with MaxQuant 2.x). Where a dynamic-exclusion setting is stated, it is cited.

**Fields left 'unknown' on purpose:**
- sample_type for leduc2022, specht2021, ai2025a and amrute2026_MFS: the sources consulted do not say whether these were
  cell lines or how the cells were isolated.
- Instrument for schoof2021, ai2025a and woo2025_PC3.
- acq for hu2023_K562.

## Data provenance (for the king to merge into `data/PROVENANCE.md`)
- All files are metadata only, fetched 2026-10-23 into `data/cache/phasec/ttprior/meta/` (git-ignored).
- No reference data were downloaded. The reference matrices come from the existing corpus via `ctx['ref']()`.

| File | URL | Licence | SHA-256 |
|---|---|---|---|
| PXD049181.json | https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD049181 | PRIDE dataset CC0 (as stated in the JSON); EMBL-EBI terms of use | `a71a78a1a0ac6e32897a8fe95a0a6293deebcaf2a44881f809d5666f65445247` |
| PXD059079.json | …/projects/PXD059079 | CC0 (JSON) | `4519f58fe147c1e43a200aa8475e1997c5aec6fbb58ffbcf93cd8b391ee40228` |
| PXD061065.json | …/projects/PXD061065 | CC0 (JSON) | `e0b44bc5449d60a9569c4e675caaa5581bdaf12442e20aaa320096f910509fbe` |
| PXD062231.json | …/projects/PXD062231 | CC0 (JSON) | `e7bf6880cb63a5bd865278d8d0b7df7f2bf070af55748a152a0c45f29aab21cf` |
| PXD063590.json | …/projects/PXD063590 | CC0 (JSON) | `f9cc5b5d32fbb027b26622a3f6eff7a518592ea389f8988c1f3d6039767107a2` |
| PXD064564.json | …/projects/PXD064564 | CC0 (JSON) | `1b54d796da87377208a7fb37bd0e7be6a4fbc12cc203cad3f18d30efd2d32acd` |
| PXD073250.json | …/projects/PXD073250 | CC0 (JSON) | `5a89c5ccd4a7c37f040caecae29c70ceddb441eb58668afa17f57bf97f655640` |
| PXD078860.json | …/projects/PXD078860 | CC0 (JSON) | `8bc1378265ae9e243895875c2e99ad3d91a4bd2424281ad20b466664758c5b5f` |
| scpdata.pdf | https://bioconductor.org/packages/release/data/experiment/manuals/scpdata/man/scpdata.pdf | scpdata package licence (Artistic-2.0 per Bioconductor; not verified by fetching) | `7301fe8448f374acbbeb45077032c232940b279bd273f3f8f244811025d3a76c` |
| scpdata_data.R | https://raw.githubusercontent.com/UCLouvain-CBIO/scpdata/devel/R/data.R | same | `a8ccbb5869dade773b9b18ac65144204846fec863ee851477d424b564031a4bd` |
| cello_CVCL_0004.txt (K-562) | https://api.cellosaurus.org/cell-line/CVCL_0004?format=txt | Cellosaurus CC BY 4.0 (not verified by fetching) | `165109f455cf14cc18ff2d64c85bba3c769b7170a7076f7bb16c2f5a9d63cb44` |
| cello_CVCL_0006.txt (THP-1) | …/CVCL_0006?format=txt | same | `2fa5c42d87000478c6fd84b59aa00927a1a780d3809e6c58e2ef6f57acb59245` |
| cello_CVCL_0007.txt (U-937) | …/CVCL_0007?format=txt | same | `50452e5fbe4fd0e01613af32d4c711f966aa77e926dc0bf1c87be057bf6081e8` |
| cello_CVCL_0023.txt (A-549) | …/CVCL_0023?format=txt | same | `0c80ca0343d95c306a55db433ad3dae695f3396fa5566dfa547e46dc2227815c` |
| cello_CVCL_0030.txt (HeLa) | …/CVCL_0030?format=txt | same | `cf62d43758f4917a6731d79d9170aeece775b690522dd6da8809df10f1e96391` |
| cello_CVCL_0035.txt (PC-3) | …/CVCL_0035?format=txt | same | `00e4299f6c8f5b68b9d3459d3ba483551e2204536db6509688f1d88c962435c5` |
| cello_CVCL_0042.txt (U2OS) | …/CVCL_0042?format=txt | same | `5bc783397ba47ec8addcb782c07696c03c5188ebed2257716df428bb27cf5b5c` |
| cello_CVCL_0598.txt (MCF-10A) | …/CVCL_0598?format=txt | same | `d4d534f1b16ed58784bf4dc5f2883e4133463091bc19ece40b0666b7006b5bb5` |

**Reproduce:**
```
python -m scp_graphfs_mcp.phasec.ttprior --unit
python -m scp_graphfs_mcp.phasec.ttprior --precompute
python benchmarks/phasec.py --cand scp_graphfs_mcp.phasec.ttprior:impute_<v> --null scp_graphfs_mcp.phasec.ttprior:impute_null_<v> --id ttprior_<v> --pool fast --jobs 2
```
`impute` / `impute_null` are FULL tau = 1, screened as id `ttprior`.
