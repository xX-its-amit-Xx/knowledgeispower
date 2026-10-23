# techhyper: a technology/tissue-kernel hyperprior on the shrinkage strength psi

- **Queen:** biological priors (round 4, from the king's round-3 directive)
- **Status (worker's recommendation):** **dead.** The metadata hyperprior moves psi\* a little, but the harness cannot
  tell the result apart from its derangement null, its flat-prior base or the plain argmin selection. Every margin is
  ≤ 0.0004, against a bar of 0.005.
  - FULL kernel: fast cand_stack **0.9568** (all) / 0.9544 (in-complex).
  - Its metadata-deranged null: **0.9569** / 0.9543. Margin 0.0001.
  - Flat-prior base (`techhyper_base`, idea switched off): **0.9567** / 0.9542. FULL is worse than the base by 0.0001.
  - Argmin selection on the extended grid (cshrink-style): **0.9571** / 0.9548.
  - FULL clears the 0.99 bar vs the champion (12/12 splits ≤ 0.99) only through the shared shrinkage machinery, so the
    full pool was not run.
- **Code:** `src/scp_graphfs_mcp/phasec/techhyper.py`
- **Results** (all in `results/phaseC/`): `techhyper_fast.json`, `techhyper_techtissue_fast.json`,
  `techhyper_base_fast.json`, `techhyper_diag.json`
- **Cache** (all in `data/cache/phasec/techhyper/`, not committed):
  - `curves.json`: the pseudo-target curves.
  - `tcurves/`: the target curves, keyed by a hash of the input matrix and shared by every entry point.
  - `leak/`: the V_{-t,-d} refits.
  - `log_<variant>.jsonl`: one row per impute call, for full / null / tech / tissue / base / argmin.
  - Precompute: `python -m scp_graphfs_mcp.phasec.techhyper --precompute`, then `--diag`.
- **Revisions:** none. The first version ran as specified and failed for a substantive reason (below), not a fixable
  bug, so no `techhyper_v2` exists.

## Idea
With 30 target cells, choosing psi (the shrinkage strength of the per-protein mean and log-sd toward a constant) on the
target's own 50/50 MCAR/MNAR inner mask is noisy. The plan was to learn from the reference datasets how much each kind
of dataset should be shrunk:
- Every human reference dataset with ≥ 40 cells serves as a pseudo-target under the exact evaluation protocol.
- A ttprior-style metadata kernel (technology + carrier + tissue) averages their psi-loss curves into a prior.
- The prior is combined with the target's own curve, and psi is set to the posterior mean.

This is hierarchical empirical Bayes on a hyperparameter, i.e. transfer of tuning-parameter priors across related tasks.

## Mechanism (as built, following the spec)
- **Grid:** GRID = (1, 2, 4, 8, 16, 32, 64, INF), with INF = psi 1e6. The grid index is g = 0..7.
- **Scoring:** identical to `tgtshrink._run`.
  - Msel = `tgtshrink.selection_mask(M, seed)`: 50/50 MCAR/MNAR, rng seed + 31337.
  - For each psi: `tgtshrink.shrink(Msel, psi, "null", seed)`, then
    `nfbase.adapt_project_mb(..., lam_v=100, extra=10, iters=5)`.
  - S = ½ (RMSE on the MCAR part + RMSE on the MNAR part).
  - Relative curve: l = log S − mean_grid log S.
- **Step 1 (precompute):** pseudo-targets are every `explore_phaseb.reference_datasets()` entry that is human and has
  ≥ 40 cells. That gives 18 datasets:
  - petrosius2023, leduc2022, specht2021, schoof2021, ai2025a, amrute2026_MFS, bubis2025_A549, gregoire2023_mixCTRL,
    guise2024, hu2023_K562, inns2026_CCS, krull2024, px063590_ADbrain, sadiku2025, weiss2026_liver,
    williams2020_tmt, woo2025_PC3, ye2025_EB.
  - Holdouts are asserted absent.
  - Each pseudo-target gets 4 replicates: {MCAR, MNAR} × seeds {900, 901}, built with the harness's own subsample /
    `make_mask` / kNN base and V_{-d} from `benchmarks/phasec.target_models`.
  - Three V_{-d} were newly fitted: bubis2025_A549 (71.5 s), hu2023_K562 (26.0 s) and ye2025_EB (82.1 s). The other 15
    were already cached.
  - Global single-realisation noise: s_t = **0.00682**.
- **Step 2 (prior, LODO):** D_t = the pseudo-targets other than t.
  - Weights: w_d = exp(−dist/1).
  - Prior curve: lp = Σ w lbar / Σ w.
  - Prior scale: s_p = max(s_t/2, sqrt(mean_psi weighted var)).
  - Prior: pi ∝ exp(−lp/s_p).
  - Kernels: FULL = tech′ + tissue, TECH = tech′, TISSUE = tissue, FLAT = uniform. NULL = FULL with
    `ttprior.null_meta_map(t, D_t)` applied to the reference rows; the target's own row stays true.
  - tech′ = `ttprior.tech_dist` + a carrier term (0 if the same, 0.5 if different, 0.25 if either side is unknown).
- **Step 3 (impute):**
  - Target curve lt, then P ∝ pi · exp(−lt/s_t) and gbar = Σ P g.
  - psi\* = INF if gbar ≥ 6.5, else 2^gbar.
  - Return `adapt_project_mb(Xm, V, warm, base, *shrink(Xm, psi*, "null", seed)[:2], lam_v=100, extra=10)`.
- **Entry points:**
  - `impute`: FULL
  - `impute_null`: NULL
  - `impute_tech`, `impute_tissue`
  - `impute_base`: FLAT prior, same posterior-mean rule (the `<id>_base` run)
  - `impute_argmin`: argmin of lt over GRID, ties to the smaller psi
- **Deviation:** the carrier field lives in `techhyper.CARRIER`, not in `ttprior.META`, so that another idea's file is
  not edited. `techhyper.meta(d)` merges the two.

### Step 0: carrier field (evidence; nothing guessed)
The rule:
- `yes` only where a carrier or booster channel is stated.
- `no` only where label-free quantification is stated.
- Otherwise `unknown`.

The sources are:
- the ttprior META `source` strings;
- the scpdata manual cached by ttprior (`data/cache/phasec/ttprior/meta/scpdata_data.R` / `scpdata.pdf`);
- the PRIDE project JSONs those strings cite (`data/cache/phasec/ttprior/meta/PXD*.json`).

| dataset | carrier | evidence |
|---|---|---|
| petrosius2023 | unknown | META: directDIA Spectronaut 17; no labelling statement (not in the cached scpdata manual version) |
| leduc2022 | yes | META "TMT-18 + carrier"; scpdata `leduc2022_pSCoPE` "carrier channels" |
| specht2021 | yes | META "SCoPE2 TMT-11/16 carrier"; scpdata `specht2019v3` "carrier channels" |
| schoof2021 | yes | META / scpdata `schoof2021` "200 cells (booster channel)" |
| ai2025a | unknown | META: DIA-NN report; no labelling statement |
| amrute2026_MFS | no | META "label-free"; PXD078860 "label-free single-cell mass spectrometry" |
| bubis2025_A549 | no | META "label-free"; PXD064564 quantificationMethods "label-free quantification" |
| cong2020AC | no | scpdata `cong2020AC` "label-free quantitative information" |
| gregoire2023_mixCTRL | yes | scpdata `gregoire2023_mixCTRL` "single-cell channels, carrier channels" |
| guise2024 | unknown | META: nanoPOTS, Exploris 480; no labelling statement |
| hu2023_K562 | no | META "label-free" |
| inns2026_CCS | no | PXD073250 "a label free SCP methodology" |
| krull2024 | unknown | META: diaPASEF timsTOF Pro; no labelling statement |
| liang2020_hela | unknown | META / scpdata `liang2020_hela`: no labelling statement |
| px063590_ADbrain | yes | META "TMTpro + carrier"; PXD063590 "126 (carrier)" |
| sadiku2025 | unknown | PXD061065 quantificationMethods "Relative quantification"; no labelling statement |
| weiss2026_liver | no | PXD062231 quantificationMethods "label-free quantification" |
| williams2020_lfq | no | META / scpdata "label-free" nanoPOTS |
| williams2020_tmt | yes | META / scpdata "carrier (10 ng)" |
| woo2022_lung | no | META "label-free" |
| woo2025_PC3 | no | META "LFQ DIA"; PXD059079 "label-free quantification" |
| ye2025_EB | no | PXD049181 "80 label-free SCP samples per day" |

## Null controls
- **`impute_null`:** identical machinery (same curves, kernel form, tau = 1, s_p, s_t and posterior-mean rule). The
  reference datasets' metadata rows are deranged by `ttprior.null_meta_map(target, D_t)` (rng 12345 + the target's
  index in sorted(META)); the target's own row stays true.
- **`impute_base`:** the same code with the prior switched off (FLAT).
- **`impute_argmin`:** cshrink's selection rule on the same extended grid.

## Data provenance
No external data was downloaded. Internal sources:
- the dev corpus (`explore_phaseb.reference_datasets` / `load_target`);
- the leave-d-out loadings in `data/cache/phasec/{d}_V.npz` (three new ones, fitted with the harness's own function);
- ttprior's cached metadata files (above).

There is nothing to add to `data/PROVENANCE.md`.

## Diagnostics before the harness (`results/phaseC/techhyper_diag.json`)
- **Within-dataset reproducibility.** The curves are reproducible inside each dataset (s_t = 0.0068) and differ between
  datasets:
  - argmin psi of lbar_d is low (1–4) for leduc2022 (2), krull2024 (2), px063590_ADbrain (2), woo2025_PC3 (2),
    hu2023_K562 (4) and williams2020_tmt (1);
  - it is high (32 to INF) for sadiku2025 (32), weiss2026_liver (32), ye2025_EB (32) and amrute2026_MFS (INF);
  - every other dataset sits at 8–64.
  - The split does not follow technology. DIA appears in both groups (krull2024 and woo2025_PC3 low; amrute2026_MFS
    and weiss2026_liver high), and TMT-carrier datasets range from psi 1 (williams2020_tmt) to psi 64 (specht2021).
  - Most curves are almost flat beyond psi ≈ 16, where tau² reaches its floor.
- **(a) LODO curve prediction.** Mean RMSE between lbar_d and the prior built from the other 17 datasets:

  | FULL | TECH | TISSUE | FLAT | NULL |
  |---|---|---|---|---|
  | 0.01439 | 0.01449 | 0.01457 | 0.01464 | 0.01451 |

  - FULL beats NULL on 11/18 datasets and FLAT on 10/18. The advantage is about 1–2% of a curve error that is itself
    about twice s_t, so the metadata barely predicts the curve. Per the spec, the fast screen was run anyway.
- **Prior strength.** The prior scale s_p is 0.011–0.018 for the fast targets, against s_t = 0.0068, so the target's
  own curve dominates the posterior. Prior-mean gbar for the six fast targets:

  | kernel | prior-mean gbar range |
  |---|---|
  | FULL | 3.71–4.39 |
  | NULL | 3.75–4.42 |

  Every prior is centred on psi ≈ 16, whatever the kernel.
- **(c) Leakage.** V_{-d} contains the target's rows, so V_{-t,-d} was refitted for the 3 top-FULL-weight d of each
  test target. The curves are essentially unchanged and the effect on psi\* is negligible.

  | target | top-3 d | curve corr. (orig vs leak-free) | prior gbar orig → leak-free | eval-call psi\* shifts |
  |---|---|---|---|---|
  | petrosius2023 | sadiku2025, bubis2025_A549, hu2023_K562 | 0.9969–0.9999 | 4.179 → 4.207 | 12.71 → 12.80 to 8.67 → 8.72 |
  | sadiku2025 | petrosius2023, inns2026_CCS, ai2025a | 0.9974–0.9998 | 3.969 → 3.974 | ≤ 0.05 |

## Results (fast pool: 6 targets × 2 masks × 2 seeds; `rmse_all` mean_rel_to_champion)

| run (cand / null) | cand_raw | cand_solo | cand_stack | null_stack | in-cx cand / null | splits ≤ 0.99 | cand weight | runtime |
|---|---|---|---|---|---|---|---|---|
| `techhyper` (FULL / NULL) | 0.9554 | 0.9559 | **0.9568** | **0.9569** | 0.9544 / 0.9543 | 12/12 | 0.7396 | 460.5 s |
| `techhyper_techtissue` (TECH / TISSUE) | 0.9554 | 0.9558 | 0.9567 | 0.9567 | 0.9543 / 0.9543 | 12/12 | 0.7396 | 342.9 s |
| `techhyper_base` (FLAT / ARGMIN) | 0.9552 | 0.9554 | **0.9567** | **0.9571** | 0.9542 / 0.9548 | 12/12 | 0.7292 | 321.0 s |
| reference: `tgtshrink` null (5-point grid, argmin) | — | — | — | 0.9575 / 0.9563 | | | | |

Per mask (mean over targets):

| run | MCAR cand_raw / cand_solo / cand_stack / null_stack | MNAR cand_raw / cand_solo / cand_stack / null_stack |
|---|---|---|
| techhyper (FULL / NULL) | 0.9719 / 0.9681 / 0.9674 / 0.9674 | 0.9390 / 0.9438 / 0.9462 / 0.9464 |
| techtissue (TECH / TISSUE) | 0.9719 / 0.9681 / 0.9674 / 0.9674 | 0.9390 / 0.9436 / 0.9460 / 0.9459 |
| base (FLAT / ARGMIN) | 0.9718 / 0.9680 / 0.9674 / 0.9672 | 0.9385 / 0.9428 / 0.9461 / 0.9470 |

- **Per target-mask.** FULL and NULL cand_stack are identical to 4 decimals on 10/12 target-masks. The two exceptions:
  - sadiku2025 MNAR: 0.9381 vs 0.9390
  - weiss2026_liver MNAR: 0.8869 vs 0.8871
- **psi\* distribution (from `log_*.jsonl`).** 48 calls per variant (eval + inner input per split). Median psi\* per
  target, with the [min, max] range for FULL:

  | target | FULL | NULL | FLAT base | argmin |
  |---|---|---|---|---|
  | gregoire2023_mixCTRL | 27.6 [25.5, 31.1] | 29.7 | 30.6 | INF |
  | leduc2022 | 17.5 [16.5, 19.0] | 13.2 | 11.1 | 2 |
  | petrosius2023 | 13.3 [8.7, 23.2] | 13.6 | 12.5 | 16 |
  | px063590_ADbrain | 20.0 [16.4, 22.3] | 19.9 | 10.8 | 2 |
  | sadiku2025 | 26.7 [17.1, 48.1] | 30.8 | 32.7 | 24 |
  | weiss2026_liver | 54.2 [17.0, 61.1] | 54.3 | 58.6 | ≈ 5·10⁵, median of 16 and INF |

  - FULL gbar ranges from 3.12 to 5.93 and never reaches INF.
  - Mean |Δgbar| of FULL against the other variants:

    | vs | mean \|Δgbar\| | max \|Δgbar\| |
    |---|---|---|
    | NULL | 0.122 | 0.422 |
    | TECH | 0.065 | 0.239 |
    | TISSUE | 0.051 | 0.101 |
    | FLAT | 0.362 | 1.347 |

  - The argmin rule spreads over the whole grid: 1 ×2, 2 ×9, 4 ×1, 8 ×4, 16 ×12, 32 ×4, 64 ×1, INF ×15.

## Honest read
- **No credit.** FULL is 0.0001 below its null and 0.0001 above the flat-prior base. TECH vs TISSUE is a tie, as is
  FLAT vs argmin (0.0004).
- **The lever is flat around psi ≥ 16.** The posterior-mean rule pulls psi\* from argmin's 2-or-INF extremes toward
  13–55, yet even that does not change the harness RMSE by more than 0.0004. On the pseudo-targets the mean curves
  span 0.003–0.107 in log-RMSE across the whole grid, almost all of it between psi = 1 and psi = 8. Across psi 16–INF
  they span ≤ 0.024, and ≤ 0.007 for 12/18 datasets. The base pipeline is therefore already near the best point of a
  shallow, saturating curve. A better psi choice cannot be worth 0.005.
- **Technology does not predict the curve shape.** DIA datasets appear at both ends: krull2024 and woo2025_PC3 prefer
  psi 2, amrute2026_MFS and weiss2026_liver prefer psi ≥ 32. The LODO curve error of FULL is only 1% below the
  deranged null. ttprior's "technology carries signal" finding does not transfer to the shrinkage hyperparameter.
- **Leakage is not the explanation.** Refitting V without the target changes the curve correlations by < 0.004 and
  psi\* by < 1%.
- **A side observation, not a credited result.** The extended grid with any reasonable rule gives 0.9567–0.9571, a
  hair below tgtshrink's null on the 5-point grid (0.9575). That difference is within noise.
- **Conclusion:** technology-aware shrinkage hyperpriors are dead, as the hypothesis stated they would be if this
  failed.

## Next step (suggestion to the queen)
- Stop tuning psi; the shrinkage lever is saturated.
- The remaining gap to the ideal in-domain reference (−14.5%) concerns covariance in the target's domain, not its
  per-protein mean/sd. A biological-prior idea would have to change the *loadings* seen by the target, not its
  standardisation.
- The ttprior lesson was that a dataset-level V reweighting is too weak. One option is protein-pair-level evidence that
  is specific to the target's tissue/technology and enters `adapt_project` as an extra covariance block, with a
  protein-permuted null.
