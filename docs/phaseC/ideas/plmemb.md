# plmemb: ProtT5 embedding-kernel smoothing of weakly supported reference loadings, inside the cshrink pipeline

- **Queen:** biological priors (round-5 directive; abd's word2vec / collaborative-filtering cold-start analogy)
- **Status (worker's recommendation):** **dead.**
  - v1 fast cand_stack is 0.9529: 0.0010 below the permuted-embedding null (0.9539), 0.0006 below the Gaussian arm
    (0.9535), and 0.0025 below the base (0.9554).
  - v2 (wider b grid) is 0.9526: 0.0012 below the null (0.9538) and 0.0004 below its Gaussian arm (0.9530).
  - Every margin is < 0.005, so the full pool was not run and `impute_nullfix` was not run (the spec runs it only on a
    pass).
- **Code:** `src/scp_graphfs_mcp/phasec/plmemb.py`.
  - Entry points: `impute`, `impute_null`, `impute_gauss`, `impute_base`, `impute_nullfix` (not run), `impute_v2`,
    `impute_v2_null`, `impute_v2_gauss`.
  - Commands: `python -m scp_graphfs_mcp.phasec.plmemb sha|diag|unit` and `... report fast <ids>`.
- **Results:**
  - Fast runs: `results/phaseC/plmemb_fast.json`, `plmemb_base_fast.json`, `plmemb_gauss_fast.json`,
    `plmemb_v2_fast.json`, `plmemb_v2_gauss_fast.json`.
  - `plmemb_diag.json`: support and leave-one-out loading-prediction check.
  - `plmemb_report_fast.json`: per-mask means and b* histograms.
  - `plmemb_smoke_fast.json`: smoke test on petrosius2023 only; it does not count.
  - Run logs: `data/cache/phasec/plmemb/fast_log_*.jsonl`.
  - Unit check: `data/cache/phasec/plmemb/unit_check.json`.

## Data provenance (row copied from `data/PROVENANCE.md`; no new download)
| file | URL | version / licence | date | SHA-256 |
|---|---|---|---|---|
| `data/cache/reps/prott5/per-protein.h5` (symlinked as `data/cache/phasec/plmemb/per-protein.h5`) | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/UP000005640_9606/per-protein.h5 | UniProt release 2026_03 ProtT5 per-protein embeddings, human UP000005640 (147,489 × 1024, float16); CC BY 4.0 (https://www.uniprot.org/help/license); MD5 checked against `RELEASE.metalink` | downloaded 2026-10-22; re-hashed 2026-10-23 | `9d713afe200557b2532aa8743886027772013c7140c6567e7084ee724af0877a` |
| `data/cache/reps/prott5/pca128.npz` | derived from the above (PCA to 128 dimensions, 147,489 accessions) | — | 2026-10-22 | (derived) |

The re-hash matches (`data/cache/phasec/plmemb/sha256_check.json`: 356,339,800 bytes, `match: true`).

Per-target derived files are `data/cache/phasec/plmemb/<target>_E.npz` (p × 128) and `<target>_support.npz`, a few
MB in total.

## Idea
Proteins with weak reference support (few reference datasets, low ‖Vref_j‖) or none (cold) have noisy or zero prior
loadings. Pull only those rows toward a ProtT5-kernel average of the loadings of strongly supported embedding
neighbours. That changes the prior centre Vp_j in the target's MAP adaptation (`adapt_project_mb`, lam_v = 100), which
the target data can still overrule.

## Mechanism (as built, following the spec)
1. **Embeddings.**
   - Each target protein group is mapped with `ids.group_accessions`, isoforms stripped.
   - Its vector is the mean of the L2-normalised PCA-128 rows of its accessions, L2-normalised again.
   - A group with no mapped accession gets has_emb = False.
2. **Support and weakness score.**
   - Per protein: n_ds_j = number of leave-target-out reference datasets with ≥ 5 observations of the protein;
     nrm_j = ‖Vref_j‖.
   - g_j = 1 − ½(rankpct(n_ds) + rankpct(nrm)) over warm proteins (average ties); cold proteins get g = 1.
   - Neighbour pool P_t = warm & has_emb & g ≤ 0.5.
3. **Kernel smoothing.**
   - Each has_emb protein takes its k = 15 nearest pool proteins by cosine, excluding itself.
   - Weights are w = exp((c − 1)/0.1), normalised, and Vs_j = Σ w Vref_i.
   - β_j = min(1, b·g_j²) for warm proteins and 1 for cold has_emb proteins (when b > 0).
   - Vtil = (1 − β) Vref + β Vs; usable = warm | (cold & has_emb & b > 0).
4. **Pipeline.**
   - ψ* is cshrink's own selection: the same `tgtshrink.selection_mask` (50/50 MCAR/MNAR, rng seed + 31337) and
     cshrink's exact on-disk score cache.
   - b is selected from {0, .25, .5, 1} on the same inner mask with iters = 5. The b = 0 score is cshrink's score at
     ψ*. Ties go to the smaller b.
   - Refit on M with iters = 10.
5. **Arms.**
   - `impute_null`: E rows permuted among the has_emb proteins with rng seed + 7177, independent of n_missing, so both
     harness calls of a split see the same permutation.
   - `impute_gauss`: iid N(0,1) 128-d vectors, normalised, rng seed + 4242.
   - `impute_base`: b = 0.
   - Permuting only among has_emb proteins keeps the has_emb set, the pool, g and the usable set identical to the
     candidate, so only the neighbour identities change.
6. **v2 (the one revision).**
   - Reason: b* sat at the grid edge (b = 1.0) in 25/48 v1 calls, a scale limit.
   - Change: the b grid is widened to {0, .25, .5, 1, 2, 4}. Everything else is unchanged; the null and Gaussian arms
     get the same grid.

**Unit check** (`unit_check.json`):
- `impute_base` vs `cshrink.run` on petrosius2023 MCAR seed 100: max |diff| = 0.0 (pass).
- vtilde(b = 0) is an exact identity.
- Fast `plmemb_base` reproduces `cshrink_fast.json` exactly: 0.9554 / cand_raw 0.9539 / cand_solo 0.9548.

## Go/no-go diagnostic (`results/phaseC/plmemb_diag.json`)
Leave-one-out kernel prediction of each pool protein's Vref row from its pool neighbours. The diagnostic uses seed 0,
so the perm and gauss arms use rng 7177 and 4242.

| target | p | cold | g > 0.5 | has_emb | pool | median n_ds | mean cos true / perm / gauss | R² true / perm / gauss |
|---|---|---|---|---|---|---|---|---|
| petrosius2023 | 2813 | 106 | 1296 | 2812 | 1516 | 11 | 0.4673 / 0.2441 / 0.2654 | 0.1086 / −0.1420 / −0.0754 |
| leduc2022 | 2844 | 297 | 1359 | 2842 | 1484 | 11 | 0.4556 / 0.2691 / 0.2957 | 0.0609 / −0.1374 / −0.0729 |
| weiss2026_liver | 4520 | 360 | 2143 | 4518 | 2375 | 8 | 0.4202 / 0.2245 / 0.2629 | 0.0569 / −0.1504 / −0.0714 |
| gregoire2023_mixCTRL | 1878 | 75 | 869 | 1875 | 1008 | 13 | 0.4750 / 0.2364 / 0.2696 | 0.1152 / −0.1432 / −0.0666 |
| px063590_ADbrain | 5140 | 1157 | 2928 | 4991 | 2209 | 7 | 0.4995 / 0.3095 / 0.3307 | 0.1034 / −0.1390 / −0.0802 |
| sadiku2025 | 3039 | 208 | 1357 | 3036 | 1680 | 10 | 0.4607 / 0.2498 / 0.2790 | 0.0825 / −0.1342 / −0.0711 |

**Go signal: yes.**
- ProtT5 neighbours share reference loadings beyond chance on every target (true R² 0.06–0.12 vs negative for both
  nulls).
- But only about 6–12% of loading variance is predictable.

## Results (fast pool: 6 targets × MCAR/MNAR × seeds 100, 101; rmse_all `mean_rel_to_champion`)
| run | cand_raw | cand_solo | **cand_stack** | **null_stack** | in-cx cand_stack / null | splits ≤ 0.99 | mean cand weight | runtime |
|---|---|---|---|---|---|---|---|---|
| plmemb (true E, b ∈ {0,.25,.5,1}) | 0.9514 | 0.9524 | **0.9529** | **0.9539** (perm E) | 0.952 / 0.9517 | 12 / 12 | 0.7917 | 238.9 s |
| plmemb_gauss | 0.9522 | 0.9531 | 0.9535 | — | 0.9522 / — | 12 / 12 | 0.7917 | 393.9 s |
| plmemb_base (b = 0 = cshrink) | 0.9539 | 0.9548 | 0.9554 | — | 0.9533 / — | 12 / 12 | 0.7917 | 309.0 s |
| plmemb_v2 (b ∈ {0,…,4}) | 0.951 | 0.9521 | **0.9526** | **0.9538** (perm E) | 0.952 / 0.9513 | 12 / 12 | 0.78125 | 443.5 s |
| plmemb_v2_gauss | 0.952 | 0.9527 | 0.953 | — | 0.9517 / — | — | — | 367.0 s |

**Per-mask means** (raw / stack / null_stack; `plmemb_report_fast.json`):

| run | MCAR | MNAR |
|---|---|---|
| plmemb | 0.9682 / 0.9641 / 0.9651 | 0.9347 / 0.9417 / 0.9428 |
| plmemb_gauss | 0.9686 / 0.9645 / — | 0.9358 / 0.9425 / — |
| plmemb_base | 0.9713 / 0.9670 / — | 0.9365 / 0.9437 / — |
| plmemb_v2 | 0.9667 / 0.9630 / 0.9645 | 0.9354 / 0.9423 / 0.9432 |
| plmemb_v2_gauss | 0.9676 / 0.9632 / — | 0.9364 / 0.9428 / — |

MCAR cand_raw is better than the base (0.9682 vs 0.9713), so the MCAR guard passes. Part of that MCAR gain is shared
by the Gaussian arm (0.9686).

**Per target-mask, cand_stack:** v1 candidate / perm null / gauss / base, then v2 / v2 null / v2 gauss.

| target-mask | v1 cand | v1 null | v1 gauss | base | v2 | v2 null | v2 gauss |
|---|---|---|---|---|---|---|---|
| gregoire2023_mixCTRL MCAR | 0.9668 | 0.9668 | 0.9668 | 0.9668 | 0.9668 | 0.9668 | 0.9668 |
| gregoire2023_mixCTRL MNAR | 0.9336 | 0.9332 | 0.9325 | 0.9325 | 0.9351 | 0.9335 | 0.9345 |
| leduc2022 MCAR | 0.9778 | 0.9768 | 0.9766 | 0.9778 | 0.9762 | 0.977 | 0.9751 |
| leduc2022 MNAR | 0.9688 | 0.9683 | 0.9682 | 0.9688 | 0.9694 | 0.969 | 0.9683 |
| petrosius2023 MCAR | 0.9641 | 0.9655 | 0.9636 | 0.9692 | 0.9613 | 0.9655 | 0.9622 |
| petrosius2023 MNAR | 0.9652 | 0.9684 | 0.9674 | 0.9668 | 0.9656 | 0.9684 | 0.9674 |
| px063590_ADbrain MCAR | 0.9808 | 0.9795 | 0.9811 | 0.9812 | 0.9788 | 0.9791 | 0.9766 |
| px063590_ADbrain MNAR | 0.9788 | 0.9814 | 0.9795 | 0.98 | 0.9787 | 0.9814 | 0.9793 |
| sadiku2025 MCAR | 0.9676 | 0.9679 | 0.9676 | 0.9676 | 0.9693 | 0.9679 | 0.9676 |
| sadiku2025 MNAR | 0.933 | 0.9311 | 0.9311 | 0.9311 | 0.9355 | 0.9331 | 0.9311 |
| weiss2026_liver MCAR | 0.9277 | 0.9339 | 0.9314 | 0.9397 | 0.9257 | 0.9305 | 0.9307 |
| weiss2026_liver MNAR | 0.8709 | 0.8743 | 0.8765 | 0.8828 | 0.8693 | 0.8737 | 0.8765 |

**Chosen b** (48 calls per arm, eval and inner calls pooled):

| arm | b = 0 | b = 0.25 | b = 0.5 | b = 1 | b = 2 | b = 4 |
|---|---|---|---|---|---|---|
| v1 true E | 18 | 4 | 1 | 25 | | |
| v1 perm E | 21 | 10 | 0 | 17 | | |
| v1 gauss E | 21 | 5 | 1 | 21 | | |
| v2 true E | 14 | 1 | 0 | 7 | 13 | 13 |
| v2 perm E | 20 | 10 | 0 | 8 | 7 | 3 |
| v2 gauss E | 20 | 2 | 1 | 8 | 8 | 9 |

- **By target (v1, true E):** weiss2026_liver chose b > 0 in 8/8 calls, leduc2022 in 6/8, gregoire2023_mixCTRL in
  2/8, sadiku2025 in 4/8.
- **Margin against cold/weak counts:** the cold-rich targets do not lead.
  - px063590_ADbrain (1157 cold): margin over min(null, gauss) is +0.0013 MCAR and −0.0007 MNAR in v1. Positive means
    worse than a null.
  - leduc2022 (297 cold): worse than both nulls on both masks.
  - The only consistent real-over-null gain is weiss2026_liver (360 cold, median n_ds 8): v1 −0.0037 MCAR and −0.0034
    MNAR vs min(null, gauss); v2 −0.0048 and −0.0044.
  - A smaller real-over-null gain is on petrosius2023 (v1 MNAR −0.0022; v2 MCAR −0.0009, MNAR −0.0018).

### Pass bar (from the directive)
| bar | v1 | v2 |
|---|---|---|
| cand_stack ≤ 0.99 | pass (0.9529) | pass (0.9526) |
| ≥ 0.005 below the null | **fail** (0.9539; margin 0.0010) | **fail** (0.9538; margin 0.0012) |
| ≥ 0.005 below the base | **fail** (0.9554; margin 0.0025) | **fail** (margin 0.0028) |
| margin also holds vs gauss | **fail** (0.9535; margin 0.0006) | **fail** (0.9530; margin 0.0004) |
| MCAR cand_raw no worse than base + 0.005 | pass (0.9682 vs 0.9713) | pass (0.9667 vs 0.9713) |

## Honest read
- **The premise is true but too weak to matter.** ProtT5 neighbours share reference loadings beyond chance: LOO R² is
  0.06–0.12 vs ≤ −0.07 for the permuted and Gaussian embeddings, on all six targets.
- **About 75% of the gain over the base is regularisation, not biology.**
  - Inside the pipeline the idea gains 0.0025 over the cshrink base.
  - Of that, 0.0019 is reproduced by Gaussian embeddings (0.9535) and 0.0015 by permuted embeddings (0.9539).
  - Shrinking weak loading rows toward a kernel average of strong rows, any strong rows, already helps. The CORUM
    lesson from Phase A ("regularisation, not curation") repeats.
- **The real-over-null part is concentrated on weiss2026_liver (DIA, liver tissue)**, about 0.004–0.005 on both masks.
  This is the same target where earlier technology-specific signals appeared.
- **Cold proteins did not help.**
  - px063590_ADbrain has 1157 cold proteins (1011 with an embedding) and no gain over the nulls.
  - A likely reason, which was not measured here: the harness masks only proteins observed in ≥ 50% of the 30 cells,
    and those tend to be well covered by the reference.
  - If so, the support-graded prior acts where there is little to score.
- **Widening the b grid (v2) moved b* to 2–4** in 26/48 calls but changed cand_stack by only 0.0003. The scale is not
  the limit.

## Next step (for the queen)
- Close the ProtT5-as-loading-prior route on this harness.
- The embedding does carry loading information (LOO R² ≈ 0.1). It could matter only where weak or cold proteins are
  actually scored, for example an evaluation stratified to low-n_ds proteins. Without it, the effect is diluted to
  ~0.001.
- The weiss2026_liver-only gain is not enough to justify a target-specific arm.
- If anything is kept, the data-only part is reusable: smoothing weak loading rows toward the strong-row kernel mean,
  which the Gaussian arm reaches at 0.9535 vs base 0.9554. Credit it as math or regularisation, not bio.
