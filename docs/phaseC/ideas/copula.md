# copula: nonparanormal (rank-Gaussian, Blom normal-score) encoding around `adapt_project`

- **Queen:** representation & mathematics (round-2 directive)
- **Status (worker's recommendation):** **dead** at the fast screen.
  - The full method fails: cand_stack is 1.0005, against a null of 1.0004.
  - The reference-only variant has the only real signal: cand_stack 0.9980, null 0.9991. That is a margin of 0.0011,
    well short of both bars (≤ 0.99, and ≥ 0.005 below the null).
  - The full pool was not run. No v2 was run either: the variant screens already locate the failure, and the refonly
    configuration a v2 would use has already been screened (see the honest read).
- **Code:** `src/scp_graphfs_mcp/phasec/copula.py`.
  - Entry points: `impute` / `impute_null` (full method), `impute_refonly` / `impute_null_refonly`,
    `impute_tgtonly` / `impute_null_tgtonly`, and `impute_identity` (unit check).
  - Command-line flags: `--unit`, `--precompute` and `--diag`.
- **Results:**
  - `results/phaseC/copula_fast.json`
  - `results/phaseC/copula_refonly_fast.json`
  - `results/phaseC/copula_tgtonly_fast.json`
  - `results/phaseC/copula_diag.json`
  - Unit checks: `data/cache/phasec/copula/unit_checks.json`

## Idea
ttprior's residual signal was on the technology axis. The hypothesis was that platforms differ in the **shape** of each
protein's marginal distribution: TMT ratio compression, intensity-dependent DIA/LFQ noise. Per-protein z-scoring
removes only location and scale.

A Gaussian-copula ("nonparanormal") encoding keeps only rank dependence (Liu, Lafferty & Wasserman 2009, *JMLR*
10:2295–2328). It replaces every protein's observed values in every reference dataset, and in the target, by Blom
normal scores (Blom 1958). The loadings should then carry co-regulation only, and the target's empirical quantiles
restore its own shape on the way back.

## Mechanism (as built, following the spec)
**One transform function** is used on both sides.
- Definition: g(z; knots_k) = Φ⁻¹(F_k(z)), where F_k is the Blom-smoothed, piecewise-linear ECDF of protein k's observed
  values.
  - Knots: the unique values, each with its average-tie Blom position (r − 3/8)/(n + 1/4).
  - F_k is clipped to [1/(2n_k), 1 − 1/(2n_k)].
- The output is re-standardised to mean 0 and sd 1 over the observed entries.
- In the real run, k = j. On its own knots g returns exactly the Blom scores: max |diff| = 0.0 in the unit check.
- Columns with fewer than 5 observed values are left unchanged.

**A. Reference**
1. `ctx['ref']()` is split by `ctx['ref_src']()`, and every dataset is encoded column-wise.
2. The NaN pattern is unchanged; this is asserted.
3. Loadings are refitted with `ttprior.fit_loadings_weighted(Zc, 1, V, warm)`: a cold start with 15 iterations, the
   champion recipe.
4. Cache: `data/cache/phasec/copula/<target>_{real,null}.npz`, written atomically. Building one takes 30–61 s per
   target and variant. All 12 caches were precomputed.

**B. Target**
1. Champion standardisation (nanmean/nanstd with the same guard), followed by the same column-wise encoding.
2. The adapt_project ALS body, copied verbatim (lam 5, lam_v 100, extra 10, iters 10, seed 0, usable = warm), with the
   prior Vc. It gives P = UVᵀ.
3. Back-transform per column, using the observed (latent s, raw x) pairs, sorted with duplicate s averaged:
   - Inside the observed latent range: `np.interp`.
   - Outside it: linear extrapolation with the OLS slope of x on s.
   - Columns with n_j < 5: mu + sd·P.
4. Only missing entries with finite mu are filled. Everything else keeps `ctx['base']`.

**Unit checks** (`data/cache/phasec/copula/unit_checks.json`; petrosius2023, seed 100):
| check | result |
|---|---|
| (i) identity h, prior `ctx['V']`, steps 4–8 vs `ctx['champion']`, max abs diff | 1.78e-15 (MCAR), 1.78e-15 (MNAR) |
| (ii) `fit_loadings_weighted(ref, w = 1)` vs `ctx['V']`, max abs diff | 1.94e-14 |
| g on its own knots vs Blom scores (synthetic, with ties) | 0.0 |
| null map monotone: Spearman(input, null-encoded), ai2025a, 300 columns | min 0.9938, mean 0.9998 |

## Null control
`impute_null` uses the identical code path, except that column j uses the knots of column π(j).
- **Reference side:** π_d is a derangement of the eligible proteins of dataset d, with rng `crc32(target|d)`.
- **Target side:** π is a derangement of the target's eligible columns, with rng `crc32(target|target|seed)`.
- **Kept:** the map is monotone in z_j, so ranks and dependence are kept; the clipping makes it only weakly monotone in
  the tails (Spearman ≥ 0.994, above).
- **Destroyed:** the shape now comes from the wrong protein.
- **Decode:** it interpolates the observed (latent, raw) pairs exactly as in the real run.
- The variants use their own nulls: the refonly null swaps only the reference shapes, and the tgtonly null swaps only
  the target shapes.

## Diagnostics (E): non-Gaussianity per acquisition (`results/phaseC/copula_diag.json`)
The statistic is the mean |Blom score − z| over observed entries (columns with ≥ 5 observed values).
- **Reference datasets:** each value is averaged over the 5–6 fast-target protein alignments in which the dataset is
  part of the reference.
- **Targets:** each value is computed on the target's full dev matrix.

| acquisition | n datasets | mean | min | max |
|---|---|---|---|---|
| TMT_DDA | 6 | 0.1201 | 0.0755 (leduc2022) | 0.1664 (px063590_ADbrain) |
| DIA | 10 | 0.1660 | 0.1151 (amrute2026_MFS) | 0.2074 (sadiku2025) |
| LFQ_DDA | 3 | 0.2110 | 0.1690 (guise2024) | 0.2352 (williams2020_lfq) |
| unknown (hu2023_K562) | 1 | 0.1619 | | |

| fast target | acquisition | mean \|Blom − z\| (full matrix) |
|---|---|---|
| leduc2022 | TMT_DDA | 0.0769 |
| weiss2026_liver | DIA | 0.1389 |
| gregoire2023_mixCTRL | TMT_DDA | 0.1410 |
| petrosius2023 | DIA | 0.1572 |
| px063590_ADbrain | TMT_DDA | 0.1613 |
| sadiku2025 | DIA | 0.2009 |

- **Shapes do differ by platform, but in the opposite direction to the hypothesis.** TMT datasets are the *most*
  Gaussian. Label-free DDA and DIA are more skewed, which fits the intensity-dependent low-tail noise of label-free
  data. Reporter-ion compression does not show up as extra non-Gaussianity.
- **The spread within each platform is about as large as the spread between platforms.**

## Results (fast pool: 6 targets × MCAR/MNAR × seeds 100, 101; rmse_all mean_rel_to_champion)
| run | cand_raw | cand_solo | **cand_stack** | **null_stack** | in-complex cand_stack / null | splits ≤ 0.99 | mean cand weight | runtime |
|---|---|---|---|---|---|---|---|---|
| copula (ref + target encoding) | 1.0274 | 1.0193 | **1.0005** | **1.0004** | 1.0011 / 1.0004 | 0 / 12 | 0.156 | 416.7 s |
| copula_refonly (copula Vc, champion target path) | 0.9983 | 0.9983 | **0.9980** | **0.9991** | 0.9976 / 0.9988 | 1 / 12 | 0.323 | 355.8 s |
| copula_tgtonly (champion V, copula target path) | 1.0284 | 1.0202 | **1.0010** | **1.0001** | 1.0014 / 1.0002 | 0 / 12 | 0.146 | 382.3 s |

### Per split (rmse_all; cand_raw / cand_stack / null_stack)

**copula (full method)**
- **Worst raw:** sadiku2025 (1.0578 / 1.0594 raw for MCAR / MNAR) and weiss2026_liver (1.0364 / 1.0682 raw).
- **weiss2026_liver MNAR:** the stack *worsens* to 1.0101, against a null of 1.0052.
- **Best stack:** px063590_ADbrain MCAR at 0.9977, against a null of 0.9999.
- **petrosius2023**, the Astral target where a gain was predicted: raw 1.0368 / 1.0074, stack 1.0000 on both masks.

**copula_refonly**
- **Only split ≤ 0.99:** gregoire2023_mixCTRL MNAR at 0.9881, against a null of 0.9998.
- **weiss2026_liver:** 0.9961 / 0.9967, against nulls of 0.9964 / 0.9987.
- **petrosius2023:** raw 0.9988 / 0.9958, stack 1.0000 / 0.9990.
- **sadiku2025 MNAR:** raw 0.9703, but the stack is 1.0002, against a null of 0.9971.
- **leduc2022**, the TMT target where a gain was predicted: 0.9988 / 1.0000.

### Pre-registered predictions vs outcome
| prediction | outcome |
|---|---|
| cand_stack < 1 and below the null | 1.0005 vs 1.0004 (fail) |
| largest gains on platform-rare targets (leduc2022, gregoire2023, petrosius2023/Astral) | not observed: the full method is at or above 1.0 on all three; the refonly gain on gregoire MNAR is the only one of note |
| cand_stack beats copula with h = identity (= champion, 1.000) | no |

## Honest read
1. **The target-side encoding is what hurts.**
   - tgtonly alone gives cand_raw 1.0284, which matches the full method's 1.0274.
   - The cause is the decode. With 30 cells, the target's empirical quantile function has only about 15–30 knots per
     column. Interpolating the latent prediction through it is a noisy, piecewise-linear map, while the champion's
     affine map averages that noise away.
   - The Blom scores of 30 points also discard the within-target magnitude information that z keeps: two cells far
     apart in z can be adjacent in rank.
   - In short, a nonparanormal decode needs many more observations per protein than 30 single cells provide. This is a
     structural limit of the setting, not a bug: the unit check shows that the code reproduces the champion exactly
     with h = identity.
2. **The reference-side encoding is a small, real but sub-threshold effect.**
   - refonly beats its null by 0.0011 (0.9980 vs 0.9991). Most of that margin comes from one split (gregoire MNAR,
     0.9881 vs 0.9998).
   - This is the same order of magnitude as every round-1 idea (ttprior 0.0040, ebshrink 0.0014, scrnaprior 0.0014).
     It does not come close to the 0.005 bar.
   - The null (wrong-protein shapes) is itself slightly below 1 (0.9991). So part of the refonly effect is simply
     "refit V on a re-encoded, slightly regularised reference", not protein-specific shape.
3. **The motivating mechanism is not supported by the diagnostics.** TMT data are the most Gaussian, not the most
   distorted. Within-platform variation in shape is as large as between-platform variation.
4. **No v2 was run.**
   - The one obvious revision, dropping the harmful target side, is exactly `copula_refonly`, and it has already been
     screened.
   - Its margin of 0.0011 is a factor of about 5 short of the bar. No scale or bug fix is in sight that would close
     that gap.

## Next step
- Close the marginal-shape branch. Affine per-protein standardisation is not the bottleneck at 30 cells.
- If the queen wants to keep the reference-side piece, stack it with ttprior's tech-weighted V (both are cheap cached
  V's) as a combined-prior test. Expect ≤ 0.004 at best.
- The more promising direction, consistent with the Phase-B ceiling diagnostic, is new information (mod's `pepev`
  per-entry precision) rather than a re-encoding of the same reference values.

## Data provenance
No external data was used. Everything comes from the dev single-cell reference exposed through `ctx['ref']` and
`ctx['ref_src']` (leave-target-out, holdouts excluded by the harness), plus each fast target's own dev matrix for the
diagnostics. No row is needed in `data/PROVENANCE.md`.

Method citations:
- Liu, H., Lafferty, J. & Wasserman, L. (2009). The Nonparanormal: Semiparametric Estimation of High Dimensional
  Undirected Graphs. *JMLR* 10:2295–2328.
- Blom, G. (1958). *Statistical Estimates and Transformed Beta-Variables.* Wiley.
