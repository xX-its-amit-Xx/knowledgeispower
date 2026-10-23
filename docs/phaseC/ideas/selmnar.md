# selmnar: Heckman probit selection model inside an EM U-step of the adapted projection

- **Queen:** representation & mathematics (round-3 directive)
- **Status (worker's recommendation):** **dead under the README criterion** (fast cand_stack 1.0000, stack weight
  0.0000; v2 1.0000). **The preregistered secondary readout passes** (MNAR cand_raw 0.8774 vs the null's 0.9615; v2
  0.8587 vs 0.9674), so the missingness pattern carries real information. The harness stack cannot use it because
  the inner mask is random. Full pool not run.
- **Code:** `src/scp_graphfs_mcp/phasec/selmnar.py`: `impute`, `impute_null`, `impute_inf`, `impute_cal`, `impute_v2`,
  `impute_null_v2`, and `python -m scp_graphfs_mcp.phasec.selmnar diag`.
- **Results:** `results/phaseC/selmnar_fast.json`, `selmnar_cal_fast.json`, `selmnar_inf_fast.json`,
  `selmnar_nullraw_fast.json`, `selmnar_v2_fast.json`, `selmnar_v2_nullraw_fast.json`, `selmnar_diag.json`.
- **Data provenance:** no external data. The only input is the target's own NaN pattern plus ctx. No PROVENANCE row.
  The analogy is from Heckman JJ (1979), "Sample selection bias as a specification error", *Econometrica*
  47(1):153-161.

## Preregistration (written before the first screen, 2026-10-23)
- `impute`: the spec's STEPS 0-3, 4 EM iterations, tilt factor rho = 1 (full Heckman posterior mean).
- `impute_null`: same machinery; the selection fit sees the missing indicators permuted within each cell among ok
  proteins (`default_rng(seed + 4242)`), with the marginal probit term for every entry; E-step on the TRUE missing
  entries.
- `impute_inf`: sigma_d fixed at 1e6 (no selection) -> machinery-only control; required max |inf - champion| / sd_j
  < 1e-3 on every fast target; harness cand_raw expected 1.000 +/- 0.002.
- `impute_cal` (preregistered variant): rho in {0, 0.25, 0.5, 1} multiplies the tilt in every E-step,
  xhat = m - rho (v/s) lam. rho is chosen inside `impute` by self-simulation: using the fitted selection curve of the
  rho = 1 fit on the input, hide 10% of the observed eligible entries (eligible = ok proteins with observed fraction
  >= 0.5 in the input) with probability proportional to Phi((c_ij - x_ij)/sigma_d) (rng `default_rng(seed + 1)`),
  rerun STEPS 1-3 on that input with mu, sd, V and s2 fixed (U restarted by one ridge U-step on the observed entries),
  and pick the rho with the lowest RMSE on the hidden entries. Ties -> smaller rho. Screened as `selmnar_cal`
  (impute_cal vs impute_null).
- Promotion test (fixed): README criterion (fast cand_stack <= 0.99 and >= 0.005 below null_stack). Secondary readout
  (not a pass): MNAR-mean cand_raw <= 0.98 and >= 0.01 below the null's MNAR cand_raw (`selmnar_nullraw`, cand =
  impute_null). Credit only the margin below both null_stack and selmnar_inf.
- **Revision v2 (added after the v1/cal/inf/nullraw screens, BEFORE running v2).** v1 tilted every missing entry,
  including random dropouts, because the selection model had no MAR component (fast MCAR cand_raw 1.2721). v2 adds a
  global random-dropout probability pi: P(miss | x) = pi + (1 - pi) Phi((c - x)/sigma_d), fitted jointly (logit pi in
  [-9, 2.2], init 0.1). Closed-form E-step: xhat = m - (v/s) (1-pi) phi(alpha) / (pi + (1-pi) Phi(alpha)); rho = 1;
  everything else as v1. Null `impute_null_v2`: the same within-cell permutation null with the mixture marginal
  likelihood. Screened as `selmnar_v2` and `selmnar_v2_nullraw` (cand = impute_null_v2, for the secondary readout).

## Idea
The champion treats every NaN as missing at random. In single-cell MS, an entry is missing partly *because* its value
is low. A probit detection curve P(obs | x) = Phi((x - c_ij)/sigma_d), with c_ij = t_i + delta_j, is fitted from
the target's own NaN pattern against the champion's low-rank predictive Gaussian N(m_ij, v_ij). Each missing entry is
then replaced by the selection-tilted posterior mean E[x | missing] = m - (v/s) phi(alpha)/Phi(alpha). This is the
inverse-Mills correction of the Heckman model.

## Mechanism (as built, following the spec)
- **STEP 0.** A verbatim copy of `transfer.adapt_project` (lam_v 100, extra 10, lam 5, iters 10, seed 0) that also
  returns U, V, mu, sd, ok and obs. Every call asserts `np.allclose(out, ctx['champion'], atol=1e-8)`; this passed on
  every call of every run. m = mu + sd·UV^T, v_j = sd_j² · s2_j, s2_j = 1.5 (n_j r_j + 5 median r)/(n_j + 5), floored
  at 0.05.
- **STEP 1.** L-BFGS-B (analytic gradient, checked against finite differences to 5e-7) on the Heckman likelihood over
  all (i, j) with ok_j. delta has a N(0, tau²) penalty with tau = median sd_j, and sigma_d is bounded to
  [0.05, 100]·median sd_j. Initialisation: t_i = 5th percentile of the observed values in the cell, delta = 0,
  sigma_d = median sd_j. The objective is divided by n·p, which does not change the minimiser.
- **STEPS 2-3.** E-step on the missing ok entries, then the ridge U-step with V fixed. Every ok protein carries an
  observation or a pseudo-observation, so the U-step is the closed form U = Z_full V (V^T V + lam I)^-1. There are
  4 EM iterations, each warm-starting STEP 1, followed by a final E-step. Proteins without observations keep the
  champion/base value.
- **impute_cal** is as preregistered above.
- **v2** is as preregistered above: a MAR + MNAR mixture selection model.

## Null controls
- `impute_null` / `impute_null_v2`: the same machinery. The selection fit sees the missing indicators permuted within
  each cell among ok proteins (`default_rng(seed + 4242)`; per-cell counts are kept), and the marginal probit term is
  used for every entry. The E-step runs on the TRUE missing entries.
- `impute_inf`: sigma_d = 1e6, which is the machinery-only control.
- The null's own cand_raw was measured in separate runs: `selmnar_nullraw` and `selmnar_v2_nullraw`.

## Results (fast pool: 6 targets × MCAR/MNAR × seeds 100, 101; rmse_all `mean_rel_to_champion`)
| run | cand_raw | cand_solo | **cand_stack** | **null_stack** | in-cx cand_stack / null | splits ≤ 0.99 | mean cand weight (null) | runtime |
|---|---|---|---|---|---|---|---|---|
| selmnar (v1, rho = 1) | 1.0748 | 1.0192 | **1.0000** | **0.9985** | 1.0000 / 0.9986 | 0 / 12 | 0.0000 (0.2604) | 724.2 s |
| selmnar_cal (rho by self-simulation) | 0.9992 | 0.9932 | **0.9999** | **0.9985** | 1.0001 / 0.9986 | 0 / 12 | 0.0104 | 1063.0 s |
| selmnar_inf (sigma_d = 1e6, no null) | 0.9996 | 1.0000 | 1.0000 | — | 1.0000 / — | 0 / 12 | 0.3854 | 306.7 s |
| selmnar_nullraw (cand = impute_null) | 0.9874 | 0.9937 | 0.9985 | — | 0.9986 / — | 1 / 12 | 0.2604 | 334.6 s |
| selmnar_v2 (MAR + MNAR mixture) | 1.0133 | 0.9962 | **1.0000** | **0.9989** | 1.0001 / 0.9990 | 0 / 12 | 0.0208 (0.2708) | 892.9 s |
| selmnar_v2_nullraw (cand = impute_null_v2) | 0.9898 | 0.9923 | 0.9989 | — | 0.9990 / — | 0 / 12 | 0.2708 | 379.0 s |

**Per-mask means.** These are the mean over the 6 targets of the JSON `per_target_mask` rows, computed from the JSON;
they are not a stored field.

| run | MCAR raw / solo / stack / null_stack | MNAR raw / solo / stack / null_stack | stack weight MCAR / MNAR |
|---|---|---|---|
| selmnar | 1.2721 / 1.0604 / 1.0000 / 0.9998 | **0.8774** / 0.9779 / 1.0000 / 0.9971 | 0.0000 / 0.0000 |
| selmnar_cal | 1.1106 / 1.0450 / 0.9998 / 0.9998 | 0.8879 / 0.9415 / 1.0000 / 0.9971 | 0.0208 / 0.0000 |
| selmnar_inf | 1.0025 / 1.0000 / 1.0000 / — | 0.9967 / 1.0001 / 1.0000 / — | 0.3958 / 0.3750 |
| selmnar_nullraw | 1.0132 / 1.0037 / 0.9998 / — | 0.9615 / 0.9837 / 0.9971 / — | 0.2292 / 0.2917 |
| selmnar_v2 | 1.1680 / 1.0374 / 1.0001 / 0.9998 | **0.8587** / 0.9550 / 1.0000 / 0.9980 | 0.0417 / 0.0000 |
| selmnar_v2_nullraw | 1.0123 / 1.0034 / 0.9998 / — | 0.9674 / 0.9811 / 0.9980 / — | 0.2292 / 0.3125 |

**Per target (rmse_all, `per_target_mask`), cand_raw MCAR / MNAR.**

| target | selmnar | null (nullraw) | inf | v2 | v2 null (nullraw) |
|---|---|---|---|---|---|
| gregoire2023_mixCTRL | 1.3577 / 0.7964 | 1.0041 / 0.9897 | 1.0044 / 1.0049 | 1.2252 / 0.7962 | 1.0041 / 0.9960 |
| leduc2022 | 1.3014 / 0.7760 | 1.0561 / 0.8686 | 1.0064 / 1.0035 | 1.3012 / 0.7758 | 1.0531 / 0.8945 |
| petrosius2023 | 1.1444 / 0.9097 | 0.9994 / 0.9935 | 1.0002 / 1.0003 | 1.0825 / 0.8972 | 0.9990 / 0.9935 |
| px063590_ADbrain | 1.2223 / 0.8303 | 1.0171 / 0.9530 | 1.0031 / 1.0007 | 1.1534 / 0.8303 | 1.0157 / 0.9540 |
| sadiku2025 | 1.4206 / 0.9954 | 1.0015 / 0.9662 | 1.0010 / 0.9707 | 1.1998 / 0.9265 | 1.0011 / 0.9679 |
| weiss2026_liver | 1.1860 / 0.9568 | 1.0013 / 0.9979 | 0.9999 / 1.0001 | 1.0459 / 0.9261 | 1.0009 / 0.9985 |

Reading the inf column: cand_raw is relative to the champion *stack*, so inf ≈ the unstacked adapted projection.
For example, sadiku2025 MNAR 0.9707 means the adapted projection alone beats the stack on that split. The fair
raw baseline is therefore inf (MNAR mean 0.9967), not 1.

### Promotion tests (fixed in advance)
- **Main (README):** fails for v1, cal and v2. cand_stack is 1.0000 / 0.9999 / 1.0000, against the bar of ≤ 0.99 and
  at least 0.005 below null_stack (0.9985 / 0.9985 / 0.9989). The candidate is never below its null in the stack.
- **Secondary (MNAR cand_raw ≤ 0.98 and ≥ 0.01 below the null's MNAR cand_raw):**
  - **v1 passes:** 0.8774 vs 0.9615 (margin 0.0841) and vs inf 0.9967.
  - **v2 passes:** 0.8587 vs 0.9674 (margin 0.1087).
  - cal: 0.8879. Its null is the v1 null (0.9615).
  - Per target, v1 is below its null on MNAR on 6 of 6 targets: gregoire 0.7964 vs 0.9897, leduc 0.7760 vs 0.8686,
    petrosius 0.9097 vs 0.9935, px063590 0.8303 vs 0.9530, sadiku 0.9954 vs 0.9662 (**worse**), weiss 0.9568 vs
    0.9979. That is 5 of 6 better, with sadiku the exception. v2 fixes sadiku (0.9265 vs 0.9679), so v2 is 6 of 6.

### Diagnostics (`results/phaseC/selmnar_diag.json`)
These come from one seed (100) per fast target, on the eval input only. The rmse ratios are vs the unstacked adapted
projection.

| target | mask | sigma_d / med sd | mean tilt (log2) | obs frac | v2 pi | null sigma_d / med sd | inf max dev / sd | mask bias (true − champion / true − cand) | rmse vs adapted: cand / null / inf / v2 |
|---|---|---|---|---|---|---|---|---|---|
| petrosius2023 | MCAR | 2.742 | 0.189 | 0.532 | 0.072 | 100.000 | 0.0300 | -0.011 / 0.312 | 1.1529 / 0.9986 / 1.0000 / 1.0862 |
| petrosius2023 | MNAR | 2.566 | 0.196 | 0.532 | 0.047 | 100.000 | 0.0314 | -0.411 / -0.063 | 0.9223 / 0.9932 / 1.0003 / 0.9046 |
| leduc2022 | MCAR | 1.597 | 0.143 | 0.352 | 0.000 | 3.400 | 0.0349 | 0.002 / 0.402 | 1.3093 / 1.0518 / 1.0004 / 1.3090 |
| leduc2022 | MNAR | 1.595 | 0.142 | 0.352 | 0.000 | 3.398 | 0.0322 | -0.423 / -0.049 | 0.7780 / 0.8685 / 1.0003 / 0.7779 |
| weiss2026_liver | MCAR | 2.841 | 0.059 | 0.427 | 0.065 | 100.000 | 0.0276 | 0.073 / 0.336 | 1.2211 / 1.0018 / 1.0001 / 1.0625 |
| weiss2026_liver | MNAR | 2.671 | 0.059 | 0.427 | 0.018 | 100.000 | 0.0346 | -0.528 / -0.289 | 0.9787 / 0.9982 / 1.0002 / 0.9345 |
| gregoire2023_mixCTRL | MCAR | 1.692 | 0.191 | 0.440 | 0.001 | 92.086 | 0.0309 | 0.002 / 0.528 | 1.3360 / 0.9997 / 0.9998 / 1.3349 |
| gregoire2023_mixCTRL | MNAR | 1.656 | 0.192 | 0.440 | 0.000 | 79.133 | 0.0294 | -0.475 / 0.042 | 0.7946 / 0.9881 / 0.9998 / 0.7945 |
| px063590_ADbrain | MCAR | 1.697 | 0.096 | 0.325 | 0.000 | 27.178 | 0.0277 | 0.023 / 0.286 | 1.2151 / 1.0023 / 1.0002 / 1.2151 |
| px063590_ADbrain | MNAR | 1.695 | 0.096 | 0.325 | 0.000 | 32.946 | 0.0404 | -0.370 / -0.122 | 0.8380 / 0.9895 / 0.9999 / 0.8380 |
| sadiku2025 | MCAR | 2.125 | 0.149 | 0.342 | 0.040 | 100.000 | 0.0278 | -0.003 / 0.471 | 1.4022 / 1.0000 / 1.0003 / 1.1645 |
| sadiku2025 | MNAR | 2.099 | 0.146 | 0.342 | 0.023 | 100.000 | 0.0267 | -0.462 / -0.027 | 1.0197 / 0.9963 / 1.0004 / 0.9483 |

- **The selection curve is identified, and it is far from MAR on every target.**
  - Fitted sigma_d is 1.6–2.8 × median sd on all six targets, well inside the upper bound of 100.
  - Under the null (permuted indicators), sigma_d runs to the bound of 100 on petrosius, weiss and sadiku, reaches
    79–92 on gregoire and 27–33 on px063590, and is 3.4 on leduc. The null keeps a tilt only on leduc (0.123), which
    is consistent with leduc's null still reaching 0.8686 on MNAR: per-cell missing counts carry information there.
  - The hypothesis that TMT/carrier targets would collapse to the champion on their own is **not** supported.
    leduc2022 (TMT/carrier) gets one of the strongest fits.
- **Bias.** Under the MNAR eval mask the champion over-predicts the hidden entries by 0.37–0.53 log2 (true −
  champion). The candidate cuts this to −0.29 to +0.04. Under MCAR the champion is unbiased (−0.01 to +0.07), and the
  candidate under-predicts by 0.29–0.53. The tilt is applied to every missing entry, and the model cannot tell a
  randomly hidden entry from a detection-limit dropout.
- **v2's MAR weight pi** is 0.00 on leduc, gregoire and px063590, and only 0.04–0.07 on petrosius, weiss and
  sadiku. It helps exactly where it is non-zero: MCAR cand_raw on weiss 1.1860 → 1.0459, sadiku 1.4206 → 1.1998 and
  petrosius 1.1444 → 1.0825, with MNAR better too. Elsewhere it is identical to v1.
- **Unit check (impute_inf): FAILED as specified, passed in the harness.**
  - The spec required max |inf − champion| / sd_j < 1e-3; the observed values are 0.027–0.040 on every target and
    mask (`inf_unit_check_pass: false`).
  - The cause is the expected one, larger than anticipated. The champion ends on a V-step, so its U is not the
    U-step fixed point for the final V. Four EM U-steps with pseudo-observations equal to the predictions converge to
    the observed-only ridge U for that V.
  - The effect on accuracy is nil: harness cand_raw for inf is 0.9996 (spec: 1.000 ± 0.002), and the diag rmse vs
    the adapted projection is 0.9998–1.0004.

## Honest read
1. **The selection model finds real, value-dependent missingness.** Under the MNAR eval mask the raw candidate is
   12–14% better than the champion stack (v1 0.8774, v2 0.8587). Against the permuted-indicator null (0.9615 / 0.9674)
   and the machinery-only control (inf, 0.9967), the margin is 8–11 points. The gain
   comes from *which* entries are missing.
2. **It cannot be cashed in under the harness rules.** Under MCAR it is 17–27% worse. The champion stack chooses its
   weights on a *random* inner mask, so the candidate gets weight 0.0000–0.04 and cand_stack = 1.0000. Under the
   README criterion it is dead.
3. **Caveat: the MNAR eval mask is itself a logistic in the value below q25.** The model partly re-learns the eval mask
   generator, because the hidden 10% of observed entries become extra low-valued NaNs in the input. The null keeps
   per-cell missing counts, so the margin over it is from the protein-level pattern. Still, the size of the MNAR gain
   is specific to this mask design. The harness cannot score native missing entries, which are where a real user
   would apply the tilt.
4. **cal did not help.** The self-simulated MNAR internal mask should almost always favour a tilt. It lowered the
   MCAR loss (1.2721 → 1.1106) only by sometimes choosing rho < 1. The per-call rho choices were not persisted,
   because they were made in worker processes. Its stack weight was still 0.0104.

## Next step (for the king)
- The idea's value depends on a decision about the evaluation, not about the method. Is MNAR the realistic deployment
  case? Missing values in the data are largely native dropouts, and the diag shows they are value-dependent on every
  target.
- If it is, the one preregistrable route is a **mask-matched stack**: the inner mask drawn from the same mechanism as
  the evaluation, or a stack mixed in fixed proportions. This needs a harness change and a king decision. Under the
  current random inner mask, no tilted method can get stack weight.
- Cheaper hedges inside the method look unpromising:
  - v2's MAR component cannot explain the MCAR holes when native missingness is 47–68%; pi stays at 0–0.07.
  - A per-entry hedge would need information that separates "randomly lost" from "below detection", for example
    evidence or precursor counts. pepev's DIA data has these.
- Do not run the full pool: the main criterion failed.
