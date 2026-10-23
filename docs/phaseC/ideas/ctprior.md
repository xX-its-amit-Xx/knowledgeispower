# ctprior: HPA-marker cell-type scores as soft per-cell weights in a weighted `fit_loadings`

- **Queen:** biological priors (round 2; directive in `KING_LOG.md`, "Directive · bio: `ctprior`").
- **Status:** **dead** on the fast pool. cand_stack 0.9990 = null_stack 0.9990, so the margin is 0.0000 (bar: ≤ 0.99
  *and* ≥ 0.005 below the null). The full pool was not run.
- **Code:** `src/scp_graphfs_mcp/phasec/ctprior.py` (`impute`, `impute_null`, `impute_half`, `impute_null_half`,
  `impute_inf`). `transfer.py` was not edited.
- **Results:** `results/phaseC/ctprior_fast.json`, `ctprior_half_fast.json`, `ctprior_inf_fast.json`.
- **Diagnostics (cache, not committed):** `data/cache/phasec/ctprior/diag_<target>.json` and
  `<target>_precompute.json`.

## Idea
The champion's prior V weights every reference cell equally. This idea re-weights reference cells by how close
their curated cell-type identity is to the target's cells. The identity is scored from HPA RNA marker knowledge, not
from proteomic similarity. This is covariate-shift importance weighting (Shimodaira 2000) at cell level. The scores
follow single-cell practice: within-cell rank signatures as in AUCell/UCell, minus an expression-bin-matched control
set as in Seurat's AddModuleScore (Tirosh et al. 2016).

**Prediction (from the queen):**
- cand_stack < null_stack on at least 4 of 6 fast targets.
- The largest gains on weiss2026_liver and px063590_ADbrain.
- The largest n_eff drop on those same targets.

## Mechanism (as built)
1. **Markers, per target.**
   - Source: HPA `pivot.npz` (20,162 genes × 154 types, nCPM), with L = log2(nCPM + 1).
   - Mapping: each target protein → Ensembl genes, using `scrnaprior._acc2genes` (proteinatlas.tsv first, then the
     HGNC fallback). A protein's RNA row is the mean over its genes.
   - Enrichment: e_jk = L_jk − median_k′ L_jk′.
   - Markers of type k: proteins with e_jk ≥ 2, top 50 by e_jk. Types with ≥ 10 markers are kept.
   - Outcome: all 154 types are kept on all 6 targets, with a median of 50 markers per type.
   - Proteins mapped to RNA: 2804/2813, 2823/2844, 4507/4520, 1855/1878, 4616/5140 and 3024/3039 (target order as
     in the results table).
2. **Raw aligned reference.**
   - Built as in `corpus_reference` (same datasets, same order) but without `_z`.
   - The row sources were asserted equal to the harness's `_ref_parts` src.
   - Note: the corpus `.npz` files are log2 and row-median-centred. Within-cell ranks are invariant to that.
3. **Within-cell percentile rank:** r_ij = rank among cell i's observed entries / n_obs_i.
4. **Abundance-matched background.**
   - a_j = mean r_ij over the reference rows, cut into 24 quantile bins.
   - 5 controls per marker from the marker's bin, excluding that type's markers, with `default_rng(0)`.
5. **Scores.**
   - s_ik = mean r over the observed markers − mean r over the observed controls.
   - s_ik = NaN if fewer than 3 markers are observed.
   - Each type is standardised by the reference mean and SD.
   - Finite reference scores: 0.888, 0.831, 0.635, 0.967, 0.600 and 0.767 of entries (same target order).
6. **Per call (on Xm).**
   - Target scores use the same marker and control sets and the reference mean and SD.
   - sbar_k is computed over target cells, needing ≥ 5 finite cells (K_used = 154 on every eval input).
   - d_i = K_used · nanmean_k (s_ik − sbar_k)². Rows with < 10 finite types get median(d).
   - w_i = exp(−d_i/τ), with τ = median(d) for `impute` and median(d)/2 for `impute_half`. `impute_inf` uses w = 1.
7. **Weighted V.** `ttprior.fit_loadings_weighted(Z, w, V, warm, iters=3, warm_start=True)` on the harness's
   z-scored reference. It is memoised by sha1 of w, and the weights are renormalised to mean 1 over the kept rows.
8. **Projection.** `transfer.adapt_project(Xm, Vw, warm, base, lam_v=100, extra=10)`.

**Runtime:**
- Precompute (markers, controls, null sets, reference S and S_null): about 8–16 s per target. Run it with
  `python -m scp_graphfs_mcp.phasec.ctprior --precompute`.
- Screens: 958 s (ctprior), 606 s (ctprior_half) and 439 s (ctprior_inf).

## Null control
- `impute_null` keeps every step and swaps each type's marker set for random non-marker proteins.
  - The draw is one-for-one per real marker, from the same decile of reference observation rate, and from the same
    abundance bin where possible.
  - That bin match held for 0.95, 0.85, 0.92, 0.91, 0.85 and 0.95 of markers (same target order as above).
  - Its own controls are redrawn by the same rule.
  - The rng is `default_rng(12345 + index of target in FAST)`.
  - No null marker equals a real marker of the same type (overlap 0 on all targets).
- `impute_inf` is the procedural control: w = 1 with the identical 3-iteration warm start, so it measures warm-start
  drift alone.

## Validity checks (on the harness's eval inputs: 2 masks × 2 seeds; from `diag_<target>.json`)
- **(a) weiss2026_liver:** hepatocytes rank **4** of 154 in sbar. This passes (top 10). The top 5 are proximal tubule
  cells, enteric stem cells, pituicytes/fscs, hepatocytes and enteric transient amplifying cells.
- **(b) px063590_ADbrain:** the top 5 are brain inhibitory neurons, early primary spermatocytes, pituitary stem cells,
  hepatocytes and retinal amacrine cells.
  - This is a partial pass: neurons rank first, and oligodendrocyte progenitor cells and oligodendrocytes are 8th and
    10th.
  - There is also off-target noise. Hepatocytes rank 4th on this brain target, and on leduc2022 as well.
- **Unrequested sanity check, sadiku2025 (sorted neutrophils):** the top 4 are erythrocyte progenitors, erythrocytes,
  neutrophil progenitors and monocyte progenitors, and neutrophils rank 7th. So the scores carry real
  cell-type signal.

## Weight diagnostics (`impute`, τ = median d; means over the 4 eval inputs)
| target | N kept | n_eff (true) | n_eff (null) | Spearman w true vs null | top gainers (share Δ) | top losers (share Δ) |
|---|---|---|---|---|---|---|
| petrosius2023 | 11793 | 10391 | 10333 | 0.613 | schoof2021 +0.065, weiss2026_liver +0.017 | specht2021 −0.041, leduc2022 −0.037 |
| leduc2022 | 10435 | 9356 | 9543 | 0.393 | schoof2021 +0.055, williams2020_tmt +0.006 | amrute2026_MFS −0.027, px063590_ADbrain −0.015 |
| weiss2026_liver | 11245 | 10159 | 10370 | 0.417 | px063590_ADbrain +0.034, schoof2021 +0.030 | leduc2022 −0.046, specht2021 −0.029 |
| gregoire2023_mixCTRL | 11889 | 10550 | 10594 | 0.525 | schoof2021 +0.044, amrute2026_MFS +0.007 | px063590_ADbrain −0.023, specht2021 −0.017 |
| px063590_ADbrain | 9366 | 8044 | 8465 | 0.543 | schoof2021 +0.048, weiss2026_liver +0.026 | leduc2022 −0.047, specht2021 −0.041 |
| sadiku2025 | 11694 | 10177 | 10550 | 0.507 | amrute2026_MFS +0.030, px063590_ADbrain +0.027 | specht2021 −0.039, leduc2022 −0.036 |

With τ = median/2, n_eff falls further, to 6313–8649 (null 6636–8889). The Spearman correlation is unchanged,
because the ranks do not depend on τ.

**Read:**
- Weights are mild: n_eff is 86–90% of N.
- The true and null weights correlate 0.39–0.61, so the null only partly decorrelates them.
- The dataset shifts are largely the *same* on every target: schoof2021 gains, and the TMT-ratio sets leduc2022 and
  specht2021 lose. That pattern is a generic score-scale or technology effect (those sets are nearly column-centred,
  so their within-cell ranks are relative, not absolute). It is not target-specific cell-type identity.
- The n_eff drop is not concentrated on weiss2026_liver or px063590_ADbrain. Its ratio is 0.86–0.90 on all targets.

## Results (fast pool, `rmse_all` mean_rel_to_champion; in-complex in brackets)
| id | cand_stack | null_stack | cand_solo | cand weight in stack | splits ≤ 0.99 |
|---|---|---|---|---|---|
| `ctprior` (τ = median d) | **0.999** [0.9986] | 0.999 [0.9987] | 0.9988 | 0.5521 | 0 / 12 |
| `ctprior_half` (τ = median/2) | 0.998 [0.9982] | 0.9982 [0.9978] | 0.998 | 0.5312 | 0 / 12 |
| `ctprior_inf` (w = 1, warm-start drift only) | 0.9998 [0.9996] | 0.9998 [0.9996] | 0.9997 | 0.4062 | 0 / 12 |

**Per target, `ctprior`** (cand_stack / null_stack, per-target means of the MCAR and MNAR values in the JSON):
- gregoire2023_mixCTRL: 0.9994 / 0.9992
- leduc2022: 0.9999 / 0.9999
- petrosius2023: 0.9989 / 0.9984
- px063590_ADbrain: 0.9998 / 0.9999
- sadiku2025: 0.999 / 0.9989
- weiss2026_liver: 0.9971 / 0.9978. This is the best split: MNAR 0.9961 vs null 0.9974.

**Per target, `ctprior_half`:**
- gregoire2023_mixCTRL: 0.9993 / 0.9989
- leduc2022: 0.9997 / 0.9995
- petrosius2023: 0.9977 / 0.9965
- px063590_ADbrain: 0.9993 / 0.9998
- sadiku2025: 0.996 / 0.9979
- weiss2026_liver: 0.9959 / 0.9967

**Kill rule:** dead unless cand_stack ≤ 0.99 *and* ≥ 0.005 below null_stack *and* below ctprior_inf's cand_stack.
- `ctprior` meets only the third condition (0.999 < 0.9998).
- Its margin over the null is 0.0000.
- `ctprior_half`'s margin is 0.0002. → **dead.**

**Prediction check:**
- cand < null on ≥ 4/6 targets: **failed.** It held on 2/6 targets with `ctprior` (weiss2026_liver, px063590_ADbrain)
  and on 3/6 with `ctprior_half` (add sadiku2025).
- Largest gain on weiss2026_liver: yes, but tiny (0.0007).
- Largest n_eff drop on the tissue targets: no.

## Honest read
- The cell-type scores work as *annotation*: hepatocytes rank 4th on liver, neurons 1st on brain, and neutrophil
  lineages lead on sadiku2025.
- As *weights on the reference*, they do nothing that random marker sets do not also do.
- Weighting itself gives a ~0.001 improvement over the w = 1 warm-start control: 0.9990 vs 0.9998, and 0.9980 with
  τ halved. But the null gets the same improvement, so it is a generic effect of down-weighting some cells and
  datasets, not biology.
- The only target-specific hint matches the prediction direction: weiss2026_liver MNAR 0.9961 vs 0.9974, and
  px063590_ADbrain 0.9993 vs 0.9998 with τ/2. It is an order of magnitude below the bar.
- **Likely reasons:**
  1. The reference has few cells of the target's types. Most reference cells are cell lines or TMT-ratio data, and
     soft weights cannot create hepatocyte co-regulation that is absent.
  2. The weighted shares are dominated by a dataset-level score-scale effect common to every target. The TMT
     datasets' ranks are relative, not absolute.
  3. Only 3–12 markers are observed per reference cell and type (median), so per-cell scores are noisy.
  4. Round 1 already showed that the prior on V is saturated: adapt_project with lam_v = 100 and 10 free factors
     re-learns most of what a re-weighted V would carry.
- No revision was made (`ctprior_v2` not run). There was no bug or scale error to fix: two τ values and the
  procedural control were all screened, and all sit within 0.002 of the champion.

## Data provenance (row for the king; copied from `data/PROVENANCE.md` / `ideas/scrnaprior.md`, no new downloads)
| file | URL | licence | downloaded | SHA-256 | notes |
|---|---|---|---|---|---|
| `data/cache/phasec/scrnaprior/rna_single_cell_type.tsv.zip` | https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip | CC BY-SA 4.0 (Human Protein Atlas) | 2026-10-23 | `2ddd3d90fb050a8efea508792ba7741c510fb5aa3d89f352aab99794ec2c9d64` | already in PROVENANCE (scrnaprior); used via its derived `pivot.npz` |
| `data/cache/features/proteinatlas.tsv.zip` | https://www.proteinatlas.org/download/proteinatlas.tsv.zip | Human Protein Atlas full table ("current") | 2026-10-22 | `ad401f8519ecdee5b67b0ad3bf28175d402d1ea0ad6c542b90fcab2bee39b470` | already in PROVENANCE; Ensembl → UniProt map |

The HGNC fallback uses `scp_graphfs_mcp.ids.hgnc_table()` (existing).

**Derived caches** (git-ignored and regenerable, in `data/cache/phasec/ctprior/`):
- `<target>_markers.npz`
- `<target>_S.npz`
- `<target>_precompute.json`
- `diag_<target>.json`
- `logs/`

## Next step
- Do not take this to the full pool. Cell-level marker weighting joins dataset-level metadata weighting as a dead
  re-weighting of the reference prior.
- If biological priors come back to cell identity, the useful part is the *annotation*. One option is to label each
  reference cell with a type and fit type-specific loadings, then choose or blend them for the target. That is a
  structural change rather than a weight on a shared V. Given the round-1 saturation, its ceiling still looks low.
- Information that the reference lacks (new modalities, or the target's own standardisation, as in the `nfbase` line)
  remains the better bet.
