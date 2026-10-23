# cshrink: selected-psi shrinkage of the target's protein mean and log-sd toward constants (+ one-global-sd arm)

- **Queen:** innovation abduction (round 4, directive in `KING_LOG.md` round-3 entry)
- **Status (worker's recommendation):** **screened-full; new best data-only pipeline, but it misses the
  preregistration-draft bar.** It passes the fast bar against its null and beats tgtshrink's null (the previous best) on
  both pools. The main result is an encoding finding, not a selection finding:
  - **The champion's per-protein standardisation is the harm.** Removing all per-protein centring *and* scaling
    (fixed psi = inf: one constant mean and one constant sd for every protein) gives fast 0.9559, against 0.9554 for
    the selected main arm.
  - Selection itself adds only 0.0015 (full) to 0.0022 (fast) over a fixed psi = 16.
  - Decision rule for the v9 draft: full MCAR mean is **0.9714 > 0.97**, and woo2025_PC3 MCAR is **1.0006 > 1.00**.
    So **no `PREREGISTRATION_v9_draft.md` was written.**
- **Code:** `src/scp_graphfs_mcp/phasec/cshrink.py`
- **Results** (in `results/phaseC/`): `cshrink_fast.json`, `cshrink_gsd_fast.json`, `cshrink_inf_fast.json`,
  `cshrink_full.json`, `cshrink_gsd_full.json`
- **Cache** (in `data/cache/phasec/cshrink/`, not committed):
  - `unit_check.json`;
  - per-call logs `fast_log_<variant>.jsonl` and `full_log_<variant>.jsonl` (variants main, null, gsd, psi16, psiinf,
    oldgrid);
  - `scores/`, an exact on-disk cache of inner-mask scores keyed by sha1(Msel), target, seed and (psi_y, psi_l). It is
    used for speed only and changes no number.

## Idea
The champion standardises every target protein with its raw 30-cell nanmean and nanstd. It then writes predictions back
as mu_j + sd_j (U V^T)_ij. Round 3 found that tgtshrink's null pipeline is in effect selected-psi James-Stein
shrinkage of these two statistics toward constants. That null scored full 0.9650. cshrink builds it explicitly, with
three changes:
- a wider grid {1, 2, 4, 8, 16, 32, 64, inf};
- a random-psi null that isolates the value of selection;
- a one-global-sd arm that answers the genomic-prediction question "should markers be standardised?".

## Mechanism (as built, as specified)
1. **Statistics.** `S = tgtshrink.target_stats(M)`, using all base_ok proteins, with cal = base_ok & n >= 5. If
   cal.sum() < 10, the champion's xbar and s are used; this never happened (median n_cal: 1750 fast, 1529 full).
2. **Shrinkage.** `cshrink_stats(M, psi_y, psi_l)`. For (ybar, vy) and (l = log s, vl):
   - t = mean over cal, r = stat - t;
   - tau2 = max(1e-4, var(r[cal]) - psi mean(v[cal]));
   - w = tau2 / (tau2 + psi v), post = t + w r.
   - psi = inf gives w = 0 exactly; psi = 0 returns the raw statistics.
   - mu = c_t + post_y and sd = exp(post_l) on base_ok; raw values elsewhere.
3. **Projection.** `nfbase.adapt_project_mb(M, V, warm, base, mu, sd, lam_v=100, extra=10, iters=10)`.
4. **Selection.** On `tgtshrink.selection_mask(M, seed)`: the 50/50 MCAR/MNAR inner mask with rng seed + 31337, built
   only from M.
   - Each psi in the grid is fitted with iters = 5.
   - Score = 0.5 (RMSE on the MCAR part + RMSE on the MNAR part).
   - Take the argmin, breaking ties toward the smaller psi, then refit on M with iters = 10.

**Entry points:**

| function | what it does |
|---|---|
| `impute` | joint psi_y = psi_l, selected on the full grid |
| `impute_null` | psi drawn uniformly from the grid; rng = default_rng(seed + 9091 + n_missing) |
| `impute_gsd` | psi_l = inf (one global sd), psi_y selected on the grid |
| `impute_psi16` | fixed joint psi = 16 |
| `impute_psiinf` | fixed joint psi = inf |
| `impute_oldgrid` | joint psi selected on {1, 2, 4, 8, 16} |

The `<id>_base` for this idea is the champion itself. psi = 0 is the champion's own standardisation (U3 below), so the
harness's 1.0 is the base.

### Unit checks
| check | result |
|---|---|
| U1 (`unit_check.json`): cshrink's selection/refit code path with `tgtshrink.shrink(mode='null')` on the old grid vs `tgtshrink.impute_null`, MCAR seed 100 | petrosius2023 max abs diff **0.0**, leduc2022 **0.0** (pass, < 1e-9) |
| U2: `cshrink_oldgrid` (null slot of `cshrink_inf_fast.json`) vs tgtshrink null (`tgtshrink_fast.json`) | all 0.957556 vs 0.957536 (gap **0.00002**); in-complex 0.95689 vs 0.956348 (gap **0.00054**); both < 0.001 (pass) |
| U3 (`unit_check.json`): psi = 0 through `cshrink_stats` + `adapt_project_mb` vs `transfer.adapt_project`, petrosius2023 MCAR 100 | max abs diff **0.0** (pass) |

The U2 values are unrounded means recomputed from the JSON rows. The rounded summaries are 0.9576 / 0.9569 vs
0.9575 / 0.9563.

## Null control
`impute_null` is the same shrink-to-constant code, the same grid and the same projection, with psi drawn at random per
call. This destroys only the selection information.

**Caveat found while reading the logs:** the harness calls the candidate twice, on the eval input and on the stacking
inner-mask input. n_missing differs between the two, so the null draws a different psi for each call. The stacking
weights are then fitted on predictions from a different psi than the ones they are applied to. Part of the null's loss
may therefore be a stacking-calibration effect, not the value of selection.

Fixed psi = 16 (`impute_psi16`) is the cleaner "selection vs no selection" control. It is reported next to the null.

## Results
All values are rmse_all `mean_rel_to_champion` from the JSONs (in-complex in brackets). Per-mask values are means of
`per_target_mask` over targets.

### Fast pool (6 targets x 2 masks x 2 seeds)
| run (cand / null slot) | cand_stack | cand_solo | cand_raw | null slot stack | cand_raw MCAR / MNAR | cand_stack MCAR / MNAR | null MCAR / MNAR | splits <= 0.99 | weight (cand / null) | worst cand target-mask |
|---|---|---|---|---|---|---|---|---|---|---|
| `cshrink` (main / random psi) | **0.9554** (0.9533) | 0.9548 | 0.9539 | 0.9637 (0.9641) | 0.9713 / 0.9365 | 0.9670 / 0.9437 | 0.9717 / 0.9558 | 12/12 (null 12/12) | 0.7917 / 0.7292 | px063590_ADbrain MCAR 0.9812 |
| `cshrink_gsd` (gsd / psi = 16) | **0.9554** (0.9533) | 0.9544 | 0.9539 | 0.9576 (0.9569) | 0.9715 / 0.9364 | 0.9667 / 0.9441 | 0.9679 / 0.9473 | 12/12 (12/12) | 0.75 / 0.75 | px063590_ADbrain MNAR 0.9830 |
| `cshrink_inf` (psi = inf / old grid) | **0.9559** (0.9535) | 0.9558 | 0.9553 | 0.9576 (0.9569) | 0.9730 / 0.9375 | 0.9679 / 0.9438 | 0.9678 / 0.9473 | 12/12 (12/12) | 0.7188 / 0.8021 | px063590_ADbrain MNAR 0.9830 |

Unrounded rmse_all cand_stack values: main 0.955352, gsd 0.955384, psi = inf 0.955854. The gsd null slot (psi = 16,
0.957601) and the inf null slot (old grid, 0.957556) differ only beyond the 4th decimal; their per-target values differ.

Worst null target-masks: random psi gregoire2023_mixCTRL MCAR 0.9835; psi = 16 px063590_ADbrain MNAR 0.9830; old grid
px063590_ADbrain MCAR 0.9812.

### Full pool (15 dev targets x 2 masks x 3 seeds)
| run | cand_stack | cand_solo | cand_raw | null slot stack | cand_raw MCAR / MNAR | cand_stack MCAR / MNAR | null MCAR / MNAR | target-masks <= 0.99 | weight (cand / null) | worst target-mask (cand / null) |
|---|---|---|---|---|---|---|---|---|---|---|
| `cshrink` (main / random psi) | **0.9609** (0.9610) | 0.9606 | 0.9605 | 0.9716 (0.9728) | 0.9790 / 0.9420 | **0.9714 / 0.9504** | 0.9793 / 0.9640 | 26/30 (null 25/30) | 0.6833 / 0.6806 | woo2025_PC3 MCAR **1.0006** / woo2025_PC3 MCAR 1.0015 |
| `cshrink_gsd` (gsd / psi = 16) | **0.9607** (0.9608) | 0.9611 | 0.9631 | 0.9624 (0.9629) | 0.9811 / 0.9452 | 0.9711 / 0.9503 | 0.9717 / 0.9531 | 26/30 (25/30) | 0.6778 / 0.675 | woo2025_PC3 MCAR 0.9994 / woo2025_PC3 MCAR 0.9986 |

The full main run took 680.5 s and full gsd 640.2 s. For reference, tgtshrink's null on the full pool scored 0.9650
(MCAR 0.9736, MNAR 0.9564).

Full main target-masks above 0.99: woo2025_PC3 MCAR 1.0006, schoof2021 MCAR 0.9972, williams2020_tmt MCAR 0.9989 and
MNAR 0.9947. The best were weiss2026_liver MNAR 0.8642 and amrute2026_MFS MNAR 0.8909.

### Chosen psi (from the logs; each split makes 2 calls, eval and stacking-inner)
| arm | pool | 1 | 2 | 4 | 8 | 16 | 32 | 64 | inf |
|---|---|---|---|---|---|---|---|---|---|
| main | fast (48) | 2 | 8 | 1 | 3 | 11 | 7 | 2 | 14 |
| gsd psi_y | fast (48) | 2 | 3 | 0 | 4 | 12 | 6 | 2 | 19 |
| old grid | fast (48) | 2 | 8 | 3 | 3 | 32 | – | – | – |
| main | full (180) | 17 | 30 | 15 | 16 | 29 | 15 | 10 | 48 |
| gsd psi_y | full (180) | 20 | 14 | 9 | 20 | 33 | 17 | 9 | 58 |
| random null | full (180, drawn) | 16 | 26 | 22 | 26 | 26 | 16 | 24 | 24 |

The choice is target-specific (full main):
- The two TMT targets select little shrinkage: williams2020_tmt picked psi = 1 in 11/12 calls, and woo2025_PC3 picked
  psi = 1 or 2 in 12/12. These are also the worst targets.
- guise2024 (64/inf in 12/12) and gregoire2023_mixCTRL (inf 8/12) select full collapse.

### tau2 and shrinkage weights (medians over calls)
| arm | pool | tau2_l | mean_w_l | tau2_y | mean_w_y | tau2_l at floor |
|---|---|---|---|---|---|---|
| main | fast | 1e-4 | 0.000113 | 1e-4 | 0.0057 | 0.79 |
| main | full | 1e-4 | 0.000162 | 1e-4 | 0.0117 | 0.74 |
| gsd | fast / full | 1e-4 | 0 | 1e-4 | 0.00097 / 0.0077 | 1.00 |
| psi = 16 | fast / full | 1e-4 | 0.000142 / 0.000139 | 0.305 / 0.442 | 0.210 / 0.417 | 1.00 |
| random null | full | 1e-4 | 0.000213 | 0.034 | 0.166 | 0.71 |

The log-sd shrinkage weight is about 1e-4 at every psi in the grid. tau2_l is at its floor because the between-protein
variance of log-sd is smaller than psi times its sampling variance. So every arm with psi >= 1 already uses (almost) one
global sd. That is why gsd = main.

## Predictions
- **P1 (explicit constant version on the wide grid matches or beats 0.9575 fast / 0.9650 full): holds.** It scored
  0.9554 fast and 0.9609 full. Full MCAR 0.9714 vs 0.9736, MNAR 0.9504 vs 0.9564; 26/30 vs 25/30 target-masks
  <= 0.99.
- **P2 (gain comes from selection; random-psi null worse by >= 0.005 fast): holds as specified, but read with care.**
  - Margin over the random null: **0.0083 fast** (0.9554 vs 0.9637) and **0.0107 full** (0.9609 vs 0.9716).
  - Against fixed psi = 16 the selection margin is only **0.0022 fast** (0.9554 vs 0.9576) and **0.0015 full**
    (0.9607 gsd vs 0.9624; same psi_y selection with one global sd).
  - Fixed psi = inf (no selection at all) is within **0.0005** of main on fast.
  - The random null is hurt by low-psi draws and by the eval/inner psi mismatch (see Null control).
  - Conclusion: *more shrinkage* is most of the effect; per-target selection adds about 0.15–0.2%.
- **P3 (scale gain comes from dropping per-protein unit-variance scaling; gsd within 0.002 of main and <= 0.9691
  fast): holds.**
  - gsd 0.9554 fast vs main 0.9554 (difference 0.00003); full 0.9607 vs 0.9609.
  - Mechanistically, main's own w_l is about 1e-4, so main already uses one global sd.
  - Stronger than predicted: psi = inf also removes per-protein *centring* and loses only 0.0005 fast. The low-rank
    projection, with 10 free extra factors, re-learns protein offsets better than a 30-cell nanmean supplies them.
  - **The champion's z-scoring is the harm.** In LDAK terms, alpha = 0 beats alpha = -1 on these data.
- **psi = inf "expected to do badly": refuted.** It scored 0.9559 fast, 12/12 target-masks <= 0.99, MNAR cand_raw
  0.9375.

## Decision rule (v9 draft)
The rule requires full cand_stack <= 0.97 on both the MCAR and MNAR means, and no full target-mask > 1.00.

| arm | MCAR mean | MNAR mean | worst target-mask | rule |
|---|---|---|---|---|
| main | **0.9714** (fail) | 0.9504 | woo2025_PC3 MCAR **1.0006** (fail) | **not met** |
| gsd | 0.9711 (fail) | 0.9503 | woo2025_PC3 MCAR 0.9994 | not met (and gsd is a secondary arm) |

No preregistration draft was written. No holdout was touched, and no `confirm_v*.py` was run.

## Honest read
- cshrink is the best data-only champion revision so far: full 0.9609 vs the champion, 0.0041 better than tgtshrink's
  null and better on both masks.
- The story is simpler than "selected EB shrinkage". The champion's per-protein standardisation of a 30-cell target
  injects noise, and almost any strong shrinkage toward one constant (psi >= 16 or inf) recovers about 4% on the fast
  pool.
- Per-target selection buys about 0.2%. Its main job is to protect the TMT targets (williams2020_tmt, woo2025_PC3),
  which prefer psi = 1–2 and are the only places where the revision does not help.
- The MCAR side (full 0.9714) is where the bar is missed. MCAR gains on TMT/carrier targets (schoof2021,
  williams2020_tmt, woo2025_PC3) are nil.
- **Recommendation to the king:** adopt `cshrink:impute` as the round-4+ base, since it is better than
  `tgtshrink.impute_null` on both pools and is the explicit code.
  - Freeze on the full grid; the psi = inf and gsd arms tie it.
  - If a single fixed rule is preferred for a preregistration, psi = inf is the simplest (no hyperparameter). It has
    only been run on the fast pool, so it needs a full-pool run before anyone relies on it.
  - techhyper (bio) is directly relevant: the technology axis (TMT wants low psi, DIA wants high psi) is visible in the
    selection logs.

## Data provenance
None. The method is data-only on the dev targets via `benchmarks/phasec.py`; no external file was downloaded. There are
no PROVENANCE rows.

## Analogies (mapping; citations as given by the queen, not re-verified online in this session)
1. **James-Stein / Efron-Morris empirical Bayes** (Efron & Morris 1975, *JASA* 70:311-319, baseball batting averages;
   Netflix-prize "global effects", Bell & Koren 2007).
   - Each protein's 30-cell mean and log-sd is a small-sample estimate, like one player's 45 at-bats.
   - Shrinking toward the grand mean, with strength from the ensemble (method of moments, psi inflation selected out of
     sample), dominates the raw estimates.
   - Result here: the shrinkage the data choose is near-total.
2. **Marker standardisation in genomic prediction** (Speed et al. 2012, *AJHG* 91:1011-1021, LDAK; verify the DOI
   before citing).
   - Shrinking log-sd toward a constant with weight w is the power scaling sd_j^w c^(1-w): w = 1 is the champion's
     z-score (alpha = -1), and w = 0 is no per-protein scaling (alpha = 0).
   - Fitted per target, w is about 1e-4, i.e. alpha = 0, as LDAK found for rare, noisy SNP columns.

## Round-5 analogies for the remaining bottleneck (covariance in the target's own domain)
The ideal in-domain reference reached -14.5% vs data-only in Phase B. Both proposals below act on covariance, not on
the standardisation, and both would be credited against the cshrink base.

**(a) MLLR speaker adaptation** (Leggetter & Woodland 1995, *Computer Speech and Language* 9:171-185; verify before
citing).
- *Historical move:* a speaker-independent acoustic model is adapted to a new speaker from seconds of speech. It does
  not re-estimate each Gaussian; it fits one shared affine transform of all the means by maximum likelihood.
- *Mapping:* the champion pulls each V_j toward Vref_j separately (lam_v = 100), spending 30 cells on p x r
  parameters. Instead fit V = Vref A (+ the free extra factors), with A an r x r matrix shared by all proteins and
  ridge-pulled toward I. The target's covariance is re-rotated and re-weighted inside the reference subspace with
  about r^2 = 400 parameters, which 30 cells x ~1,500 proteins can support.
- *Selection:* the ridge weight toward I is chosen on the 50/50 inner mask.
- *Null:* A fitted on a protein-permuted target.
- *Differs from ebshrink:* that was per-protein lam_v, not a shared transform.

**(b) Kriging with a locally re-estimated variogram = Ledoit-Wolf covariance shrinkage** (Ledoit & Wolf 2004,
*J. Multivariate Analysis* 88:365-411; verify before citing).
- *Historical move:* geostatisticians blend the local empirical variogram with a regional prior; they do not krige a
  new field with another region's variogram.
- *Mapping:* Sigma_t = (1 - delta) S_target + delta Sigma_ref.
  - S_target is the pairwise-complete covariance on the cshrink-standardised target.
  - Sigma_ref = V V^T + diag, from the reference.
- *Imputation:* the conditional Gaussian mean given each cell's observed proteins (Woodbury on the low-rank part), with
  delta selected on the 50/50 inner mask.
- *Null:* S_target computed after permuting cells independently within each protein. This keeps the marginals and
  destroys the in-domain covariance.
- *Differs from dca (dead):* dca inverted a target-only precision with no reference blend and no selection.

**Cross-pollination.**
- rep could treat (b) as the covariance counterpart of mnarjoint.
- bio's techhyper could supply a technology-aware prior on delta or on the A-ridge. The cshrink logs show that TMT
  targets want little shrinkage and DIA targets want a lot.
- New observation for the king: because psi = inf ties main, the base could drop per-protein standardisation entirely,
  which simplifies every round-5 encoding.

## Next step
- Promote `cshrink:impute` to the base for rounds 4+ (king's call).
- Run `cshrink_inf` on the full pool before choosing between selected psi and the parameter-free psi = inf for any
  preregistration.
- Look at the MCAR side on TMT/carrier targets, where neither arm helps.
