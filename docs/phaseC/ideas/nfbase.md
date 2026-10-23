# nfbase: Netflix-Prize baseline predictors for the target's per-protein location and scale

- **Queen:** innovation abduction (round-2 directive in `KING_LOG.md`)
- **Status (worker's recommendation):** **dead as a biological-prior idea** at the fast screen. The protein-specific
  reference content clears 0.99 but not the 0.005 margin over its null.
  - Main run (psi = 1): cand_stack is 0.9868 (all) / 0.9859 (in-complex). The protein-permuted null is 0.9867 / 0.9863,
    so the margin is -0.0001 / 0.0004.
  - Best a-priori variant (psi = 4): cand_stack is 0.9663 / 0.9647, null 0.9691 / 0.9688. The margin is 0.0028 /
    0.0041, below the 0.005 bar.
  - The full pool was not run, because no variant met both bars.
  - **Flag for the king: the null itself is a large effect.** Shrinking the target's 30-cell per-protein mean and SD
    toward a calibrated constant, with psi = 4, gives null_stack 0.9691. That is the largest fast-pool effect in
    Phase C so far, and it is present on every target. It is a *data-only* change to the champion's standardisation,
    not new biological information (see "Honest read").
- **Code:** `src/scp_graphfs_mcp/phasec/nfbase.py`
- **Results** (all in `results/phaseC/`):
  - `nfbase_fast.json`, `nfbase_technull_fast.json`, `nfbase_psi4_fast.json`, `nfbase_flat_fast.json`
  - `nfbase_meanonly_fast.json`, `nfbase_sdonly_fast.json`, `nfbase_meanonly_psi4_fast.json`,
    `nfbase_sdonly_psi4_fast.json`
- **Cache and diagnostics** (in `data/cache/phasec/nfbase/`, not committed):
  - `unit_check.json`, `diagnostics.json` (go/no-go), `feasibility.json` (cross-dataset correlations)
  - `<target>_ref.npz`, `log_*.txt`

## Idea
The champion standardises the target with per-protein `np.nanmean` / `np.nanstd` over the 30 masked cells. It puts each
prediction back as mu_j + sd_j * (U V^T)_ij, so any error in mu_j or sd_j goes straight into every imputed entry of
protein j.

The Netflix-Prize analogy (Bell & Koren, *SIGKDD Explorations* 2007; Koren, *KDD* 2008) maps as follows:
- The protein is the "item".
- Its relative abundance and log-dispersion are the "item bias".
- The baselines get shrunk before the factors.

The shrinkage prior is a hierarchical protein + technology x protein model, fitted to the *unstandardised* dev
reference.

## Mechanism (as built, following the abd spec)
0. **`adapt_project_mb`.** This is a verbatim copy of `transfer.adapt_project` with mu and sd injected as arguments.
   - `ok` and the sd guard are unchanged.
   - Unit check: psi = 0 reproduces `transfer.adapt_project` on MCAR seed 100 of all six fast targets. The max-abs
     difference is 1.3e-15 to 6.2e-15 (`unit_check.json`, all pass < 1e-10).
1. **Reference baselines.**
   - Datasets: `explore_phaseb.reference_datasets(exclude=(target,))`. The code asserts that none of derks2022, wu2026,
     brunner2022 or khan2023 is present. `extended_corpus()` is never called directly.
   - Raw `external._load` matrices are aligned to `names` by the first accession.
   - Per protein with n >= 10:
     - relative location y = mean - median_j(mean)
     - log-scale l = log(sd)
     - sampling variances s^2/n and 1/(2(n-1))
   - Datasets with < 200 such proteins are dropped: cong2020AC, liang2020_hela. With petrosius2023 as the target, 19
     datasets remain.
   - Group = `ttprior.META[d]['acq']`, imported, not copied. hu2023_K562 falls in 'unknown'.
   - Cached atomically per target.
2. **Hierarchical EB.** This is fitted separately for y and l.
   - Precision-weighted group means.
   - sigma2_e[g] by method of moments over proteins with D_gj >= 2, iterated twice. A group with no such protein gets
     the pooled mean.
   - a_j = mean of the group means.
   - sigma2_b by method of moments, floored at 1e-4.
   - Target prior: p_j = a_j + k_j (theta_{g*j} - a_j), with k_j = D sigma2_b / (D sigma2_b + sigma2_e[g*]).
3. **Fusion** (inside `impute`, using Xm only).
   - Centring: c_t = median of xbar_j.
   - Calibration: OLS ybar ~ alpha + beta p over proteins with a prior and n_j >= 5.
   - tau2 = max(1e-4, var(resid) - psi * mean(s^2/n)).
   - Weights: w = tau2 / (tau2 + psi s^2/n).
   - Posterior: mu = c_t + w ybar + (1 - w)(alpha + beta p).
   - The log-sd gets the same recipe with variance 1/(2(n-1)).
   - The champion's mu and sd are kept where n_j < 2, sd <= 1e-6, or there is no prior.
4. **Prediction.** `adapt_project_mb(Xm, V, warm, base, mu, sd, lam_v=100, extra=10)`.

Variants:
- **A priori:** psi4 (psi = 4) and flat (a single group).
- **Triggered** (main run <= 0.995): meanonly and sdonly at psi = 1.
- **Added after seeing psi4:** meanonly_psi4 and sdonly_psi4. These are attribution only and are labelled post hoc.

## Null controls
- **`impute_null` (protein-permuted).** Each reference dataset's (y, l, v_y, v_l) vectors are permuted jointly among
  that dataset's *own covered proteins*, with rng `default_rng(4242 + sorted index of d)`.
  - I deviated from the spec here by permuting within covered proteins only, not over all p. This keeps the coverage
    pattern, and so the set of proteins that have a prior, identical to the real run.
  - The EB refit, calibration, weights and `adapt_project_mb` are unchanged.
  - As intended, the calibration then learns beta ≈ 0: for example -0.077 on petrosius2023 location and -0.105 on
    weiss2026_liver, against 1.153 and 0.845 with the true prior (`diagnostics.json`).
- **`impute_technull`.** `ttprior.null_meta_map` deranges the acq labels across the kept reference datasets (rng
  12345 + index of the target in sorted(META), no fixed point). The target keeps its true group. It was run as its own
  harness id with `--null impute_null`.
- Every variant has a matching protein-permuted null: `impute_<v>_null`.

## Feasibility and go/no-go diagnostic (dev reference only, recomputed)
**Cross-dataset feasibility.** `feasibility.json` gives pairwise Pearson r between kept reference datasets, computed
from the petrosius2023-aligned baselines.

| Pair type | Location r: pairs, min / median / max | log-sd r: min / median / max |
|---|---|---|
| DIA-DIA | 36 pairs, -0.00 / 0.378 / 0.687 | -0.041 / 0.283 / 0.573 |
| TMT-TMT | 15 pairs, -0.378 / -0.027 / 0.09 | -0.131 / 0.131 / 0.484 |
| LFQ-LFQ | 3 pairs, 0.345 / 0.457 / 0.589 | not reported |
| Cross-group | median location r 0.068 | not reported |

The SD of relative location is 0.051 (leduc2022), 0.033 (specht2021) and 0.105 (px063590_ADbrain): these datasets
are ratio-normalised. This confirms the abd feasibility probe.

**Go/no-go diagnostic.** `diagnostics.json` covers 6 targets x MCAR/MNAR x seeds 100/101. It measures the masked-entry
RMSE of the mean-only predictor (mu_j alone), averaged over the 6 targets:

| Predictor | Mean-only RMSE |
|---|---|
| champion nanmean | 0.8828 |
| nfbase psi = 1 | 0.8743 |
| psi = 4 | 0.8763 |
| protein-permuted null | 0.8734 |
| technull | 0.8743 |
| flat | 0.8745 |

The null was already as good as the real prior on location, so I expected little protein-specific signal. The
harness was run anyway, because the diagnostic does not measure the scale.

Per-target EB fit (seed 100):
- **DIA targets:** the calibration slope on the location prior is strong. beta is 1.153 (petrosius2023), 0.845
  (weiss2026_liver) and 0.476 (sadiku2025). But tau2 is about 1.8, so the data weight stays at about 0.94 to 0.95. The
  target's own location is already precise relative to the prior's residual spread.
- **TMT targets:** beta on location ≈ 0, with tau2 0.001 (leduc2022) to 0.007 (px063590_ADbrain).
- **Scale weights:** mean w is about 0.48 to 0.73 on all targets, true or null.

## Results (fast pool, `rmse_all` / `rmse_in_complex` mean_rel_to_champion, exactly as in the JSONs)
| id | cand_raw | cand_solo | cand_stack (all / in-cx) | null_stack (all / in-cx) | margin (all) | splits <= 0.99 | stack weight |
|---|---|---|---|---|---|---|---|
| nfbase (psi 1) | 0.9846 | 0.9864 | 0.9868 / 0.9859 | 0.9867 / 0.9863 | -0.0001 | 6/12 | 0.7188 |
| nfbase_technull | 0.9848 | 0.9866 | 0.9868 / 0.9860 | 0.9867 / 0.9863 | — | 6/12 | 0.7188 |
| nfbase_flat | 0.9848 | 0.9866 | 0.9869 / 0.9859 | 0.9867 / 0.9863 | -0.0002 | 6/12 | 0.7188 |
| **nfbase_psi4** | 0.9644 | 0.9663 | **0.9663 / 0.9647** | 0.9691 / 0.9688 | 0.0028 | 12/12 | 0.7604 |
| nfbase_meanonly (psi 1) | 0.9918 | 0.9923 | 0.9916 / 0.9910 | 0.9915 / 0.9909 | -0.0001 | 4/12 | 0.4375 |
| nfbase_sdonly (psi 1) | 0.9934 | 0.9947 | 0.9947 / 0.9957 | 0.9950 / 0.9960 | 0.0003 | 2/12 | 0.7396 |
| nfbase_meanonly_psi4 (post hoc) | 1.0021 | 1.0004 | 0.9920 / 0.9916 | 0.9923 / 0.9922 | 0.0003 | 4/12 | 0.3125 |
| nfbase_sdonly_psi4 (post hoc) | 0.9832 | 0.9836 | 0.9832 / 0.9845 | 0.9829 / 0.9844 | -0.0003 | 7/12 | 0.6979 |

nfbase_psi4 per target x mask (cand_stack vs null_stack, from `nfbase_psi4_fast.json`):

| Target | MCAR | MNAR |
|---|---|---|
| weiss2026_liver | 0.949 vs 0.9645 | 0.9281 vs 0.9312 |
| petrosius2023 | 0.9751 vs 0.9827 | 0.979 vs 0.979 |
| sadiku2025 | 0.9767 vs 0.9815 | 0.9678 vs 0.9672 |
| gregoire2023_mixCTRL | 0.968 vs 0.9691 | 0.9425 vs 0.944 |
| leduc2022 | 0.9784 vs 0.9789 | 0.9683 vs 0.968 |
| px063590_ADbrain | 0.9791 vs 0.9792 | 0.9834 vs 0.9834 |

Runtime per screen: 126 to 358 s (`meta.runtime_s`).

## Honest read
1. **The protein-specific reference content is real, but small.** It appears only where abd predicted: DIA targets,
   with heavy shrinkage (psi = 4).
   - In nfbase_psi4, the true prior beats the protein-permuted null by 0.0155 on weiss2026_liver MCAR, 0.0076 on
     petrosius2023 MCAR and 0.0048 on sadiku2025 MCAR.
   - On the TMT targets the gap is about 0, which matches the beta ≈ 0 the calibration learns there.
   - The six-target mean margin is 0.0028, below the 0.005 bar.
   - At psi = 1 the margin is nil. The data weights on the location stay at about 0.95, so the prior barely moves the
     baseline.
2. **The technology x protein term carries nothing beyond a single pooled group.** technull (0.9868) and flat (0.9869)
   equal the main run (0.9868).
   - Per-target OLS calibration (alpha, beta) absorbs most of what the hierarchy distinguishes.
   - The bio finding that the technology axis carries signal does not transfer to baselines in this form.
3. **The big effect is the machinery, not the biology.** The protein-permuted null at psi = 4 is 0.9691, and is
   <= 0.99 on all 12 splits. What it does is James-Stein shrinkage of the target's own per-protein mean and log-SD
   toward a target-level constant (beta ≈ 0).
   - The effect is super-additive. The location-only arm gives 0.9920 and the scale-only arm 0.9832, but joint
     shrinkage gives 0.9663 (null 0.9691).
   - The champion's raw 30-cell standardisation is thus a measurable weak point that a *data-only* fix improves by
     about 3% on dev.
   - By the colony's rules this does not count as a biological gain, because the null destroys the biology and keeps
     the gain. It is still the largest effect seen in Phase C.
   - I did not tune psi beyond the two a-priori values. A grid would be fishing on the dev split.
4. **No v2 was made.** There was no bug or scale error to fix. The unit check is exact, and the null behaves as designed
   (beta ≈ 0).

## Next step (for the king)
- **Do not promote nfbase as a biological prior.**
- **Consider a separate, preregistered rep/abd idea, `tgtshrink`.** This would be target-only empirical-Bayes
  shrinkage of the per-protein location and log-scale inside `adapt_project`, with no reference baselines at all.
  - It needs its own null. One option is a shrink-to-permuted-target-constants control. A cleaner framing is to treat
    it as a champion revision and judge it against the champion on the full pool.
  - Its shrinkage strength psi should be chosen on the inner mask, not on the outer metric.
  - If it holds on the full pool, re-test nfbase's DIA-only protein prior *on top of* it. The margin there (0.0048 to
    0.0155 on DIA MCAR splits) is the only biological signal this idea found.
- **For rep's `copula`:** it composes with this split. nfbase covers location/scale; copula covers shape.

## Data provenance
No external data was downloaded. The idea uses only:
- the existing dev corpus npz files, through `scp_graphfs_mcp.external._load`, restricted to
  `explore_phaseb.reference_datasets(exclude=(target,))`;
- the metadata table `scp_graphfs_mcp.phasec.ttprior.META`, whose sources are in `ideas/ttprior.md`.

There are no new `PROVENANCE.md` rows.

abd's disclosed earlier probe called `transfer.extended_corpus()`, which returns brunner2022 and wu2026. None of its
numbers are used here. This worker never called that function directly. The king may want `extended_corpus()` to
exclude `PHASE_B_HOLDOUTS` itself, or assert that it does.
