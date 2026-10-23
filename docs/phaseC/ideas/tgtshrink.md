# tgtshrink: Fay-Herriot / limma-trend shrinkage of the target's per-protein mean and log-sd

- **Queen:** innovation abduction (round-3 placeholder in `LEDGER.md`; spec from the abd queen)
- **Status (worker's recommendation):** **no trend credit.** The shrinkage machinery is the largest fast-pool effect in
  Phase C so far, but the trend covariates (abundance, observed fraction) make it *worse* than their own permuted-covariate
  null.
  - Main run (trend + psi chosen per target on an inner mask): fast cand_stack **0.9597** (all) / 0.9586 (in-complex),
    11/12 splits <= 0.99. null_stack **0.9575** / 0.9563. The margin is **-0.0022**, the wrong sign (bar: >= +0.005).
  - The null is the strongest number: shrink-to-(near-)constant with per-target psi selection, 0.9575, which beats
    nfbase_psi4's null (0.9691) by 0.0116.
  - Prediction (1) holds (both masks below the champion). Predictions (2) and (3) fail (see "Predictions").
  - Full pool (15 dev targets x 2 masks x 3 seeds; run because fast cand_stack <= 0.99, per the spec): cand_stack
    **0.9683** / 0.9684, 24/30 target-masks <= 0.99; null_stack **0.9650** / 0.9653. Margin -0.0033, again the wrong
    sign. Both the trend version and its null clear the champion-revision bar (< 1 on MCAR and MNAR, <= 0.99 overall);
    the null is better on both masks.
- **Code:** `src/scp_graphfs_mcp/phasec/tgtshrink.py`
- **Results** (in `results/phaseC/`): `tgtshrink_fast.json`, `tgtshrink_psi4_fast.json`,
  `tgtshrink_const_psi4_fast.json`, `tgtshrink_full.json`
- **Cache** (in `data/cache/phasec/tgtshrink/`, not committed): `unit_check.json`; per-call logs `log_<variant>.jsonl`
  (`fast_log_main.jsonl` / `fast_log_null.jsonl` are the fast-pool rows of the main run and its null, snapshotted before
  the full-pool run appended to `log_main.jsonl` / `log_null.jsonl`).

## Idea
The champion standardises each target protein with its raw 30-cell nanmean and nanstd and writes predictions back as
mu_j + sd_j (U V^T)_ij, so noise in these two statistics enters every imputed entry. nfbase's null showed that shrinking
them toward a calibrated constant at psi = 4 gives 0.9691. Fay-Herriot small-area estimation (Fay & Herriot 1979,
*JASA* 74:269-277) says the shrinkage target should be a covariate regression ("synthetic estimate") and the shrinkage
strength should be estimated; limma-trend (Smyth 2004, *SAGMB* 3:3) says log-variance depends on abundance; survey
non-response says the response rate (observed fraction) belongs in the location regression under MNAR.

Mapping: area -> protein; respondents -> the n_j observed cells; direct estimate -> nanmean, log nanstd; sampling
variance -> s^2/n, 1/(2(n-1)); area covariates -> abundance, observed fraction.

## Mechanism (as built)
Everything is computed from the matrix M the harness passes in (eval input or inner-mask input). No external data.
0. Stats: n_j, f_j = n_j / n_cells, xbar_j, s_j (ddof 0), base_ok = finite & n >= 2 & s > 1e-6, cal = base_ok & n >= 5,
   ybar = xbar - c_t (c_t = median xbar over n >= 2), vy = s^2/n, l = log s, vl = 1/(2(n-1)).
1. Trend designs. Covariates standardised over cal and clipped to the cal range.
   - Scale: D_l = [1, ns(m), f, f^2, f m] with m = standardised xbar.
   - Location: D_y = [1, ns(f)].
   - ns: natural cubic spline in the ESL (eq. 5.4-5.5) truncated-power form, boundary knots at the cal min / max and
     interior knots at the 10/50/90% cal quantiles. **Deviation:** with three fixed interior knots plus boundary knots a
     natural spline has 4 non-constant columns (R's `ns` with 3 interior knots), not 3. The fallbacks ([1, f, f^2] for
     location, cubic polynomial for m) apply if the knots are not strictly increasing; they were never triggered in the
     fast-pool logs (all `ns_m` and `ns_f` = "ns").
   - WLS over cal with weights n - 1.
2. Fay-Herriot, separately for location and scale: r = stat - t, tau2 = max(1e-4, var(r[cal]) - psi mean(v[cal])),
   w = tau2 / (tau2 + psi v), posterior = t + w r. Non-base_ok proteins keep the champion's xbar and s.
3. `nfbase.adapt_project_mb(M, V, warm, base, mu, sd, lam_v=100, extra=10)`, iters = 10.
4. psi selection (`impute` / `impute_null`): Msel hides 10% of the observed entries of proteins with observed fraction
   >= 0.5 in M; half uniform, half with `run.make_mask`'s MNAR weights expit((q25_j - x) / (0.5 sd_j)) computed on M;
   rng = default_rng(seed + 31337). For psi in {1, 2, 4, 8, 16}: steps 0-2 on Msel, `adapt_project_mb(..., iters=5)`
   with ctx['base'] (preregistered speed deviation). Score = mean(RMSE on the uniform part, RMSE on the MNAR part);
   argmin, ties to the smaller psi; refit on M.

Entry points: `impute` (trend + selected psi), `impute_null`, `impute_psi4` / `impute_psi4_null` (fixed psi = 4),
`impute_const_psi4` (D = [1], **unweighted** OLS, restricted to nfbase-covered proteins, psi = 4; unweighted so that it
matches nfbase's `_fuse` OLS exactly, as required by unit check (b)).

### Unit checks (`data/cache/phasec/tgtshrink/unit_check.json`)
- (a) psi = 0 vs `ctx['champion']` (`transfer.adapt_project`) on MCAR seed 100 of all six fast targets: max-abs
  difference 0.0 on all six (pass < 1e-10). Note that psi = 0 is short-circuited in `shrink` (returns xbar, s), so this
  checks the `adapt_project_mb` path and the guard logic, not the floating-point identity c_t + (xbar - c_t) = xbar.
- (b) Constant-trend posteriors vs `nfbase.baselines(Xm, c*1, c*1, psi=4)` on nfbase-covered proteins: max-abs
  difference <= 1.8e-15 for mu and <= 1.8e-15 for log sd on all six targets (pass < 1e-8). Harness reproduction:
  `tgtshrink_const_psi4` cand_stack 0.9691 vs nfbase_psi4's null_stack 0.9691, difference 0.0000 (within 0.002).

## Null control
`impute_null` keeps the full pipeline (same designs, WLS, Fay-Herriot weights, the same Msel and psi grid, the same
`adapt_project_mb`). Before the trend fit it permutes the covariate rows (m_j, f_j) jointly across cal proteins
(rng = default_rng(seed + 4242)); non-cal base_ok proteins get the covariates of a random cal protein from the same rng.
Sampling variances stay with their own proteins. The permuted trend is near-constant (mean trend R^2 0.001-0.007 in
`fast_log_null.jsonl`), i.e. shrink-to-constant with per-target psi. `impute_psi4_null` is the same at fixed psi = 4.

## Results (fast pool: 6 targets x 2 masks x 2 seeds; rmse_all mean_rel_to_champion)

| id | cand_raw | cand_solo | cand_stack (all / in-cx) | null_stack (all / in-cx) | splits <= 0.99 | weight |
|---|---|---|---|---|---|---|
| tgtshrink (trend, psi selected) | 0.9607 | 0.9640 | 0.9597 / 0.9586 | 0.9575 / 0.9563 | 11/12 | 0.677 |
| tgtshrink_psi4 (trend, psi = 4) | 0.9742 | 0.9756 | 0.9730 / 0.9726 | 0.9683 / 0.9684 | 11/12 | 0.573 |
| tgtshrink_const_psi4 (D = [1], nfbase cover) | 0.9670 | 0.9692 | 0.9691 / 0.9688 | — | 12/12 | 0.750 |

Per mask (mean over the 6 targets of the per-target-mask ratios in `per_target_mask`):

| id | mask | cand_raw | cand_solo | cand_stack | null_stack |
|---|---|---|---|---|---|
| tgtshrink | MCAR | 0.9782 | 0.9749 | 0.9708 | 0.9672 |
| tgtshrink | MNAR | 0.9433 | 0.9531 | 0.9487 | 0.9478 |
| tgtshrink_psi4 | MCAR | 0.9900 | 0.9851 | 0.9788 | 0.9754 |
| tgtshrink_psi4 | MNAR | 0.9583 | 0.9662 | 0.9671 | 0.9613 |
| tgtshrink_const_psi4 | MCAR | 0.9787 | 0.9772 | 0.9760 | — |
| tgtshrink_const_psi4 | MNAR | 0.9552 | 0.9612 | 0.9621 | — |

Per target-mask cand_stack / null_stack (tgtshrink): gregoire2023_mixCTRL MCAR 0.9744 / 0.9665, MNAR 0.9299 / 0.9329;
leduc2022 MCAR 0.9753 / 0.9784, MNAR 0.9752 / 0.9694; petrosius2023 MCAR 0.9672 / 0.9688, MNAR 0.9540 / 0.9681;
px063590_ADbrain MCAR 0.9785 / 0.9810, MNAR 0.9920 / 0.9799; sadiku2025 MCAR 0.9735 / 0.9641, MNAR 0.9507 / 0.9439;
weiss2026_liver MCAR 0.9559 / 0.9446, MNAR 0.8905 / 0.8929. The only split above 0.99 is px063590_ADbrain MNAR (0.9920).

Chosen psi (`fast_log_*.jsonl`; 48 calls = 24 splits x {eval, inner-mask input}):
- trend: psi 1: 3, 2: 5, 4: 4, 8: 7, **16: 29**. Per target: gregoire 16 x 8; leduc 4 x 3, 16 x 5; petrosius 8 x 6,
  16 x 2; px063590 1 x 3, 2 x 5; sadiku 4, 8, 16 x 6; weiss 16 x 8.
- null: psi 1: 2, 2: 9, 4: 2, 8: 5, **16: 30**.
- Trend R^2 (mean over calls, trend run): location 0.006 (gregoire), 0.010 (leduc), 0.433 (petrosius), 0.019
  (px063590), 0.346 (sadiku), 0.454 (weiss); scale 0.020-0.270.
- tau2_l hits the 1e-4 floor in most calls (mean_w_l = 0.000 on 4/6 targets in the trend run): the log-sd is replaced
  almost entirely by the trend (or, in the null, by a constant).

### Full pool (`tgtshrink_full.json`, runtime 946 s)

| mask | cand_raw | cand_solo | cand_stack | null_stack |
|---|---|---|---|---|
| all (mean_rel_to_champion) | 0.9737 | 0.9699 | 0.9683 (in-cx 0.9684) | 0.9650 (in-cx 0.9653) |
| MCAR (mean of 15 targets) | 0.9941 | 0.9799 | 0.9776 | 0.9736 |
| MNAR (mean of 15 targets) | 0.9534 | 0.9599 | 0.9589 | 0.9564 |

- Stack weight 0.614. cand_stack <= 0.99 on 24/30 target-masks; null_stack <= 0.99 on 25/30.
- cand_stack > 0.99: inns2026_CCS MCAR 0.9929, px063590_ADbrain MNAR 0.9911, schoof2021 MCAR 0.9975,
  williams2020_tmt MCAR 0.9992 / MNAR 0.9912, woo2025_PC3 MCAR 1.0002 (the only one above the champion).
- The trend beats its null on 13/30 target-masks. Largest trend wins: guise2024 MNAR 0.9515 vs 0.9636, petrosius2023
  MNAR 0.9493 vs 0.9623. Largest losses: inns2026_CCS MNAR 0.9672 vs 0.9400, weiss2026_liver MCAR 0.9382 vs 0.9155,
  amrute2026_MFS MCAR 0.9616 vs 0.9431.
- Chosen psi (180 calls each, full-pool rows of `log_main.jsonl` / `log_null.jsonl` after the first 48): trend
  1: 28, 2: 40, 4: 19, 8: 23, 16: 70; null 1: 20, 2: 31, 4: 17, 8: 20, 16: 92.

## Predictions (preregistered in the spec)
1. cand_stack <= 0.99 with both masks below the champion: **holds** (0.9597; MCAR 0.9708, MNAR 0.9487).
2. Trend beats its permuted-covariate null by >= 0.005, largest on MNAR: **fails.** Fast margin -0.0022 overall, -0.0036
   on MCAR, -0.0009 on MNAR; full-pool margin -0.0033 overall, -0.0040 MCAR, -0.0025 MNAR. The trend wins on 6/12 target-masks (petrosius both, gregoire MNAR, leduc MCAR,
   px063590 MCAR, weiss MNAR) and loses on the other 6.
3. Fixed psi = 4 with the trend beats 0.9691: **fails** (0.9730; its own null 0.9683, the constant reference 0.9691).

## Honest read
- The gain is the shrinkage machinery plus the psi selection, not the covariates. Selecting psi per target on an inner
  mask moves shrink-to-constant from 0.9691 (nfbase null, fixed psi = 4, nfbase cover only) to 0.9575 (this null,
  all base_ok proteins). The covariate trend gives that back partly (0.9597 trend vs 0.9575).
- The selection pushes psi to the grid edge (16) in about 60% of calls, and the scale tau2 is floored: the best
  per-protein sd is close to a single common sd. In other words, the champion's per-protein unit-variance scaling looks
  harmful; replacing it with a common (or trend) scale lets the ridge-regularised loadings absorb the protein-specific
  spread. This is a statement about the encoding, not about biology.
- The abundance/observed-fraction trend in log-sd is real (R^2 up to 0.27) but using it as the shrinkage target makes
  the scale vary in a way the imputer does not want; the location trend on f (R^2 0.35-0.45 on the DIA/Astral-like
  targets petrosius, sadiku, weiss) helps on petrosius MNAR (0.9540 vs 0.9681) but not on average.
- Per the spec's own fallback: no trend credit; tgtshrink is a candidate only as a constant-shrinkage revision, and the
  version to carry forward is the null pipeline (`impute_null`: near-constant target, per-target psi), not `impute`.
- I did not make a v2: nothing failed for a fixable reason (no bug, unit checks pass, the constant reference reproduces
  nfbase to 4 decimals).

## Next step
- On the full pool the champion-revision bar is met by both pipelines, but the trend has no credit, so the number to
  quote for the revision is the null pipeline's 0.9650 (full) / 0.9575 (fast), and it should be re-run as a named,
  un-permuted constant version before any preregistration.
- Promote the *constant* version explicitly: a clean `impute_const_sel` (D = [1], no permutation, all base_ok proteins,
  psi selected on the inner mask) with an extended grid {1, ..., 16, 32, 64, inf}, since 16 is the modal choice; psi =
  inf on the scale is "one common sd". Its natural null is the champion (psi = 0). Then full pool.
- For the other queens (varprior, seqprior): credit any covariate only as the margin over this selected-psi constant
  shrinkage (0.9575 fast), a stricter baseline than 0.9691.
- Adoption as a champion revision needs a full-pool win on both masks, PREREGISTRATION_v9 and Amit's approval (king).

## Data provenance
None. No external data was downloaded. `data/cache/phasec/nfbase/<target>_ref.npz` (existing) is read only to define
protein coverage for `impute_const_psi4` and unit check (b). No PROVENANCE rows.
