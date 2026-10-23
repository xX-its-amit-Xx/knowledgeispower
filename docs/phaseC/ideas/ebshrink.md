# ebshrink: empirical-Bayes per-protein adaptation strength lambda_v,j

- **Queen:** representation & mathematics (round-1 directive)
- **Status (worker's recommendation):** **dead** at the fast screen.
  - cand_stack is 0.9986, which misses the 0.99 bar.
  - The margin over the null is 0.0014, below the 0.005 bar.
  - The one revision (v2) did not help. The full pool was not run.
- **Code:** `src/scp_graphfs_mcp/phasec/ebshrink.py`, with `impute` and `impute_null`, plus `impute_v2`, `impute_null_v2`
  and the diagnostic `impute_nobs`, which was not run.
- **Results:** `results/phaseC/ebshrink_fast.json` and `results/phaseC/ebshrink_v2_fast.json`.

## Idea
The champion's `adapt_project` uses MAP estimation with the prior V ~ N(V_ref, I/lam_v) and one global lam_v = 100.
In BLUP / James–Stein terms, the right shrinkage of each protein toward the reference is sigma2/tau2_j, where tau2_j
is that protein's between-context loading variance. We have about 17–21 dev reference datasets. Each one gives a
replicate estimate of each protein's loading, so tau2_j can be estimated from the data, and lambda_j ∝ 1/tau2_j can
replace the scalar.

## Mechanism (as built, following the spec)
1. **Per-dataset fits.** This works per target, on the leave-target-out reference from `ctx['ref']` and
   `ctx['ref_src']`. For each reference dataset with at least 40 cells:
   - Subsample at most 300 cells (rng seed 0).
   - Run a local copy of `fit_loadings` (rank 20, lam 5, 8 iterations, min_obs 20) that also returns U.
   - Compute a ridge sampling-noise term for each protein:
     noise_dj = s2_dj · tr((U_o'U_o + 5I)^-1) / 20.
2. **Alignment.** Procrustes with scale onto the pooled V, over S = warm ∩ warm_d, skipping a dataset when |S| < 50.
   The noise term is scaled by s_d².
   - All 6 fast targets used 17 datasets.
   - The Procrustes scales s_d ranged from 0.40 to 1.07; the petrosius2023 log is an example.
3. **Heterogeneity estimate.**
   - tau2_j = max(raw_j − noise_j, 0.1·median(raw)).
   - It is then shrunk toward the median with prior weight 3. Proteins with n_j < 2 get the median.
4. **Strengths.** lambda_j = c/tau2_j, with the median over warm proteins set to 100, clipped to [10, 1000]. Proteins
   that are not warm get 100.
5. **Adapted projection.** `transfer.adapt_project` is copied verbatim (lam 5, extra 10, iters 10, seed 0). The only
   change is the per-protein lam_j in the V-step, which applies to all 30 columns.
6. **Sanity check (C).**
   - Setup: petrosius2023, seed-100 subsample, 20% random mask.
   - Result: with lambda_j = 100 everywhere, the copy matches `T.adapt_project` exactly (max |diff| = 0.0).
   - This is a run-time check only and was not saved to a JSON file.
7. **Cache.** `data/cache/phasec/ebshrink/<target>_tau.npz` holds lam, tau2, n, raw, noise, warm, n_obs_ref and the
   datasets used. Building it takes 11–18 s per target.

**v2 revision (dose-response).** It uses the same tau2 cache with the contrast squared: lambda_j ∝ 1/tau2_j², median
100, clipped to [10, 1000].
- **Reason:** in v1 the 5–95% range of lambda was only about 56–151, so the per-protein contrast might have been too
  weak to register.
- **Expected if the signal were real:** a stronger contrast should improve the result.

## Null control
`impute_null` / `impute_null_v2` run the same pipeline, with the lambda vector permuted among warm proteins
(`default_rng(1000 + seed)`).
- **Kept:** the distribution of strengths, the clipping, the median and all of the machinery.
- **Destroyed:** which protein gets which strength.

## Diagnostics (from the cache files, D)
| target | warm | Spearman(tau2, ref n_obs) | lambda q0 / 5 / 25 / 50 / 75 / 95 / 100% | at clip 10 / 1000 | median noise/raw |
|---|---|---|---|---|---|
| gregoire2023_mixCTRL | 1803 | −0.032 | 18.4 / 56.2 / 84.8 / 100 / 118.6 / 151.0 / 226.1 | 0 / 0 | 0.050 |
| leduc2022 | 2547 | 0.119 | 31.8 / 62.9 / 86.7 / 100 / 113.8 / 144.6 / 205.5 | 0 / 0 | 0.058 |
| petrosius2023 | 2707 | 0.193 | 19.4 / 59.1 / 84.3 / 100 / 117.2 / 148.1 / 202.2 | 0 / 0 | 0.060 |
| px063590_ADbrain | 3983 | 0.318 | 32.3 / 64.7 / 87.8 / 100 / 112.7 / 141.8 / 184.6 | 0 / 0 | 0.083 |
| sadiku2025 | 2831 | 0.065 | 26.6 / 59.9 / 84.9 / 100 / 115.4 / 147.5 / 193.2 | 0 / 0 | 0.056 |
| weiss2026_liver | 4160 | 0.281 | 31.5 / 63.7 / 86.9 / 100 / 114.9 / 144.8 / 202.4 | 0 / 0 | 0.062 |

- **Not a missingness proxy:** tau2 has only a weak correlation with the reference observation count (Spearman −0.03
  to 0.32).
- **No clipping:** no protein hit either clip bound in v1. In v2 at most 0.4% sat at 10 and none at 1000.
- **Small noise correction:** the subtracted noise is about 5–8% of the raw spread. So tau2 is mostly real
  between-dataset spread, or alignment error, rather than ridge noise.

These diagnostics were computed by a scratch script from the cache files. They are not in a `results/` JSON; they can
be regenerated from `data/cache/phasec/ebshrink/*_tau.npz`.

## Results (fast pool: 6 targets × MCAR/MNAR × seeds 100, 101; rmse_all mean_rel_to_champion)
| run | cand_raw | cand_solo | **cand_stack** | **null_stack** | in-complex cand_stack / null | splits ≤ 0.99 | mean cand weight | runtime |
|---|---|---|---|---|---|---|---|---|
| ebshrink (v1) | 0.9996 | 0.9991 | **0.9986** | **1.0000** | 0.9983 / 0.9999 | 1 / 12 | 0.292 | 351.8 s |
| ebshrink_v2 | 1.0028 | 1.0014 | **0.9988** | **1.0002** | 0.9985 / 1.0001 | 1 / 12 | 0.167 | 194.8 s |

- **Per split, v1:** the only split at ≤ 0.99 is gregoire2023_mixCTRL MNAR, with cand_stack 0.9876 against a null of
  0.9988. Every other split lies between 0.9981 and 1.0001.
- **Per split, v2:** gregoire2023_mixCTRL MNAR is again the only split at ≤ 0.99 (0.9883 against 0.9993).
- **Pre-registered predictions:** the spec predicted cand_stack ≤ 0.99, null ≈ 1.00 and cand_solo < 1.00. We got
  cand_stack 0.9986 (fail), null 1.0000 (as predicted), and cand_solo 0.9991 (directionally right, but negligible).

## Honest read
- **Statistical status.** The null control behaves, and the candidate sits about 0.0014 below it. That is consistent in
  sign but tiny, and it comes mostly from a single split (gregoire MNAR).
- **No dose-response.** Squaring the contrast (v2) made the solo projection *worse* (cand_raw 0.9996 → 1.0028,
  cand_solo 0.9991 → 1.0014). The stack gain stayed at about the same 0.001, because the stack simply down-weights the
  candidate (weight 0.29 → 0.17). If cross-dataset stability really told the projection how hard to anchor each protein,
  a stronger contrast should have helped.
- **Conclusion.** Per-protein adaptation strength looks like the scalar lam_v: the adapted projection is insensitive to
  how strongly it is anchored around lam_v ≈ 100. The between-dataset loading spread that we measured (mostly real
  spread, since noise is only about 6%) is not information the 30-cell target can exploit through the V-step prior.
- **Recommendation.** Mark the idea dead. The "per-protein lambda" branch of the seed backlog is closed in the
  1/tau2 form.

## Next step
Switch encoding, following the queen's fallback: RBF/GP kernel ridge over the reference loadings (a kernel across
proteins) as a fourth stack member.
- A cheap side-use of this idea's cache: tau2_j could serve as a per-protein *feature*, for example to gate the kNN or
  low-rank weight per protein, rather than as a prior strength.
- That is a different mechanism and needs its own null.

## Data provenance
No external data. Everything comes from the dev single-cell reference already exposed through `ctx['ref']` and
`ctx['ref_src']` (leave-target-out, holdouts excluded by the harness). No row is needed in `data/PROVENANCE.md`.
