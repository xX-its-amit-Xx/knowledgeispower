# mnarjoint: missingness as a second channel of a joint low-rank encoding [Z_obs | w·(M − mean)]

- **Queen:** representation & mathematics (round-4 directive; adjacent pivot from `selmnar`)
- **Status (worker's recommendation):** **dead** (fast cand_stack 0.9576 = null_stack 0.9576; base 0.9575, so the
  idea is 0.0001 *worse* than its own w = 0 base). Full pool not run.
- **Code:** `src/scp_graphfs_mcp/phasec/mnarjoint.py`: `adapt_project_mj`, `impute`, `impute_null`, `impute_base`,
  `impute_w1`, `impute_w1_null`; `python -m scp_graphfs_mcp.phasec.mnarjoint unit` and `... report`.
- **Results:** `results/phaseC/mnarjoint_fast.json`, `mnarjoint_base_fast.json`, `mnarjoint_nullraw_fast.json`,
  `mnarjoint_w1_fast.json`, `mnarjoint_unit.json`, `mnarjoint_report.json` (per-mask means and chosen-w counts,
  computed from the result JSONs and the run logs `data/cache/phasec/mnarjoint/log_*.jsonl`).
- **Data provenance:** no external data; the only input is the target's own NaN pattern plus ctx. No PROVENANCE row.
  Method: Sportisse A, Boyer C, Josse J (2020) "Imputation and low-rank estimation with Missing Not At Random data",
  *Statistics and Computing* 30:1629-1643, https://doi.org/10.1007/s11222-020-09963-5 (DOI as verified in the
  directive's idea spec; not re-verified by this worker).

## Idea
selmnar showed the target's NaN pattern is value-informative (MNAR cand_raw 0.8774 vs permuted-indicator null 0.9615)
but a global probit tilt of every hole cost 17-27% on MCAR. Here the cell factors U must also explain *which*
proteins are missing (low-rank fit of [X | M], Sportisse et al. 2020), so a hole is moved only through the shared
factors, entry by entry, and the channel strength w is selected on the 50/50 MCAR/MNAR inner mask (w = 0 = base).

## Mechanism (as built, following the spec)
- **Base** (round-4 rule; cshrink not landed): tgtshrink.impute_null's pipeline, imported, not copied:
  `shrink(M, psi, "null", seed)`, psi from {1,2,4,8,16} on `tgtshrink.selection_mask` (rng seed + 31337, half
  uniform, half MNAR-weighted), score 0.5·(RMSE_MCAR-part + RMSE_MNAR-part), inner fits iters = 5, then
  `nfbase.adapt_project_mb` (lam_v 100, extra 10, lam 5, iters 10, seed 0).
- **adapt_project_mj:** adapt_project_mb plus an M block on ok columns: Mc = Mo − Mo.mean(0), Wm = w·Mc, free loadings
  Bt (Gaussian surrogate). Per iteration: (a) U[i] = solve(V_oᵀV_o + BtᵀBt + lam I, V_oᵀZ_io + Btᵀ Wm[i]);
  (b) champion ridge V-step unchanged; (c) Bt = Wmᵀ U (UᵀU + lam_b I)⁻¹, lam_b = 5. Output exactly as adapt_project_mb.
- **w selection:** w ∈ {0} ∪ {0.1, 0.3, 1, 3}, fitted on Msel with the chosen psi (iters 5), scored on the same two
  hole sets; ties → smaller w. Final fit on M, iters 10.
- **Unit check** (`mnarjoint_unit.json`, seed-100 MCAR input built as the harness builds it, all 6 fast targets):
  max |impute_base − tgtshrink.impute_null| = 0.0 and max |adapt_project_mj(w = 0) − adapt_project_mb| = 0.0 (also
  with the M-block branch forced on at w = 0). all_pass = true.

## Null and references
- `impute_null`: same machinery and w selection, but the indicators fed to the M block are permuted within each cell
  across ok proteins (rng seed + 4242, rows in order), in both the inner-selection fits and the final fit.
- `mnarjoint_base` (w forced to 0) = round-4 base run; `mnarjoint_nullraw` (cand = impute_null) = null's cand_raw.
- `mnarjoint_w1` (fixed w = 1, no selection) vs `impute_w1_null`: diagnostic of the channel at full strength.

## Results (fast pool: 6 targets × MCAR/MNAR × seeds 100, 101; rmse_all `mean_rel_to_champion`)
| run | cand_raw | cand_solo | **cand_stack** | **null_stack** | in-cx cand_stack / null | splits ≤ 0.99 | mean cand weight | runtime |
|---|---|---|---|---|---|---|---|---|
| mnarjoint (w selected) | 0.9555 | 0.9575 | **0.9576** | **0.9576** | 0.9564 / 0.9564 | 12 / 12 | 0.8021 | 424.6 s |
| mnarjoint_base (w = 0) | 0.9555 | 0.9575 | **0.9575** | 0.9575 | 0.9563 / 0.9563 | 12 / 12 | 0.8021 | 196.6 s |
| mnarjoint_nullraw | 0.9555 | 0.9575 | 0.9576 | — | 0.9564 / — | 12 / 12 | 0.8021 | 414.3 s |
| mnarjoint_w1 (w = 1 fixed) | 0.9716 | 0.9673 | 0.9654 | 0.9773 | 0.9631 / 0.9753 | 11 / 12 | 0.5521 | 459.8 s |

The base reproduces tgtshrink_fast null_stack (0.9575 / 0.9563) exactly, as required (bar: within 0.001).

**Per-mask means** (mean over the 6 targets of the JSON `per_target_mask` rows; from `mnarjoint_report.json`):

| run | MCAR raw / solo / stack / null_stack | MNAR raw / solo / stack / null_stack |
|---|---|---|
| mnarjoint | 0.9698 / 0.9684 / 0.9673 / 0.9673 | 0.9412 / 0.9465 / 0.9478 / 0.9478 |
| mnarjoint_base | 0.9697 / 0.9684 / 0.9672 / 0.9672 | 0.9412 / 0.9465 / 0.9478 / 0.9478 |
| mnarjoint_nullraw | 0.9698 / 0.9684 / 0.9673 / — | 0.9412 / 0.9465 / 0.9478 / — |
| mnarjoint_w1 | 0.9896 / 0.9798 / 0.9749 / 0.9871 | 0.9535 / 0.9548 / 0.9560 / 0.9675 |

**Per target, cand_raw MCAR / MNAR:** base vs w1 — gregoire2023_mixCTRL 0.9665 / 0.9267 vs 0.9684 / 0.9263;
leduc2022 0.9861 / 0.9693 vs 1.0051 / 0.9811; petrosius2023 0.9757 / 0.9690 vs 1.0001 / 0.9877;
px063590_ADbrain 0.9831 / 0.9796 vs 1.0171 / 1.0016; sadiku2025 0.9624 / 0.9190 vs 0.9786 / 0.9290;
weiss2026_liver 0.9446 / 0.8835 vs 0.9685 / 0.8954. mnarjoint (selected) differs from the base by ≤ 0.0002 on every
target-mask.

**Chosen w** (log rows matched to splits by NaN-pattern signature; 12 calls per cell):

| | MCAR eval | MCAR inner | MNAR eval | MNAR inner |
|---|---|---|---|---|
| impute | 0: 9, 0.1: 2, 0.3: 1 | 0: 8, 0.1: 2, 0.3: 2 | 0: 11, 0.1: 1 | 0: 9, 0.1: 1, 0.3: 2 |
| impute_null | 0: 10, 0.1: 2 | 0: 11, 0.1: 1 | 0: 12 | 0: 10, 0.1: 2 |

w = 1 and w = 3 were never chosen. Mean relative change of U on the final fit, mean|U − U_w0| / mean|U_w0|:
impute 0.0188 (MCAR eval) / 0.0035 (MNAR eval); w1 0.4098 / 0.4168; w1 null 0.4713 / 0.4986.

### Pass bar (fixed in the directive)
- cand_stack ≤ 0.99: **pass** (0.9576), but only through the base.
- ≥ 0.005 below null_stack: **fail** (0.9576 vs 0.9576, margin 0.0000).
- ≥ 0.005 below mnarjoint_base cand_stack: **fail** (0.9576 vs 0.9575, margin −0.0001).
- MCAR cand_raw no more than 0.005 worse than base: pass (0.9698 vs 0.9697).
→ Full pool not run.

## Honest read
- The inner-mask selection does what it was designed to do as a safety valve: it switches the channel off (w = 0 in
  77 of 96 selected calls; never above 0.3), so the candidate equals the base to 4 decimals. There is no gain to credit.
- At full strength (w = 1) the channel is **harmful on both masks**, including MNAR: MNAR cand_raw 0.9535 vs base
  0.9412, MCAR 0.9896 vs 0.9697. The true pattern is less harmful than the permuted one (w1 cand_stack 0.9654 vs its
  null 0.9773; MNAR raw 0.9535 vs null-stack 0.9675 is not a like-for-like raw comparison; the w1 null's raw was not
  run), so the pattern carries information, but in this encoding it never beats ignoring it.
- The selection picked w > 0 slightly more often on MCAR splits than on MNAR splits: the opposite of the hypothesis.
- Likely reason (not tested): the value loadings are anchored hard to the reference (lam_v = 100 toward Vp, extra
  columns toward 0), so U moved by the free M block cannot be matched by V; the shift of U (≈ 40% at w = 1) mostly adds
  noise to UVᵀ instead of lowering the specific holes. In Sportisse et al. both blocks' loadings are free. selmnar's
  MNAR raw gain (0.8774) came from an explicit value-dependent tilt, which the shared-factor route does not reproduce.
- No revision was made: the failure is not a bug or a scale error (unit checks exact; w grid spans 0.1-3 and the
  smallest w already degrades the inner score in most calls).

## Next step
Close the [X | M] shared-factor route on the anchored projection. If missingness is revisited, it should be as a
per-protein, inner-mask-calibrated detection-limit offset inside the shrinkage (where the lever is), with per-mask
cand_raw against a permuted-indicator null, not through U.
