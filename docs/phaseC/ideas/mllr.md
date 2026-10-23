# mllr: MLLR-style shared r0 x r0 transform of the reference loadings, on the cshrink base

- **Queen:** innovation abduction (round 5, directive in the round-4 entry of `KING_LOG.md`)
- **Status (worker's recommendation):** **dead.**
  - Fast cand_stack is 0.9552. The base (`mllr_base`, A = I, identical to cshrink) is 0.9554, so the margin is 0.0002.
  - The protein-permuted null scores 0.9538, which is **0.0014 better than the candidate**.
  - The diagonal arm is 0.9550 (0.0004 below the base).
  - No arm meets the bar of ≥ 0.005 below both the null and the base, so the full pool was not run.
  - The TMT/carrier prediction is also falsified: A barely moves on those targets.
- **Code:** `src/scp_graphfs_mcp/phasec/mllr.py`
- **Results** (in `results/phaseC/`): `mllr_fast.json`, `mllr_base_fast.json`, `mllr_diag_fast.json`, plus the
  housekeeping run `cshrink_inf_full.json`.
- **Cache** (in `data/cache/phasec/mllr/`, not committed):
  - `unit_check.json`;
  - per-call logs: `fast_mllr_log_{main,null}.jsonl`, `fast_mllr_base_log_base.jsonl`,
    `fast_mllr_diag_log_{diag,base}.jsonl`, and `log_tmtdiag.jsonl` (the TMT diagnostic);
  - smoke-test logs `smoke_log_{main,null}.jsonl` (the smoke JSON was deleted);
  - harness stdout `*.out` / `cshrink_inf_full.log`.

## Step 0 (housekeeping): `cshrink_inf` on the full pool

This was the king's request: the parameter-free ψ = ∞ arm (one constant mean and one constant sd for every protein),
with the old-grid selected-ψ arm in the null slot as in its fast run. Run:

```
.venv/bin/python benchmarks/phasec.py --cand scp_graphfs_mcp.phasec.cshrink:impute_psiinf \
    --null scp_graphfs_mcp.phasec.cshrink:impute_oldgrid --id cshrink_inf --pool full --jobs 2
```

It took 1737.3 s, with another worker sharing the CPUs.

**Overall (all / MCAR / MNAR, cand_stack vs the champion).** Per-mask values are means over the 15 per-target rows of
`per_target_mask`, computed the same way as cshrink's documented 0.9714 / 0.9504.

| run | all | MCAR | MNAR | cand_raw all / MCAR / MNAR | target-masks ≤ 0.99 | worst target-mask | stack weight |
|---|---|---|---|---|---|---|---|
| `cshrink_inf_full` (ψ = ∞) | **0.9607** | 0.9704 | 0.9509 | 0.9653 / 0.9818 / 0.9487 | 26/30 | woo2025_PC3 MCAR 1.0009 | 0.664 |
| its null slot (old-grid selected ψ) | 0.9631 | 0.9728 | 0.9535 | | | | |
| `cshrink_full` (selected ψ, for reference) | 0.9609 | 0.9714 | 0.9504 | 0.9605 / 0.9790 / 0.9420 | 26/30 | woo2025_PC3 MCAR 1.0006 | 0.683 |

In-complex cand_stack is 0.9606 for ψ = ∞ and 0.9610 for selected ψ.

**Target-masks above 0.99 for ψ = ∞:**
- woo2025_PC3 MCAR 1.0009;
- schoof2021 MCAR 0.9979;
- woo2025_PC3 MNAR 0.9969;
- williams2020_tmt MCAR 0.9946.

williams2020_tmt MNAR is 0.9900, exactly at the bar.

**TMT/carrier rows** (cand_raw / cand_solo / cand_stack / null_stack):

| target-mask | ψ = ∞ (`cshrink_inf_full`) | selected ψ (`cshrink_full`) |
|---|---|---|
| woo2025_PC3 MCAR | 1.0404 / 1.0039 / 1.0009 / 1.0006 | 1.0189 / 1.0005 / 1.0006 / 1.0015 |
| woo2025_PC3 MNAR | 1.0250 / 1.0029 / 0.9969 / 0.9898 | 0.9717 / 0.9865 / 0.9898 / 0.9954 |
| williams2020_tmt MCAR | 1.0161 / 1.0015 / 0.9946 / 0.9989 | 1.0049 / 1.0006 / 0.9989 / 0.9971 |
| williams2020_tmt MNAR | 0.9957 / 0.9883 / 0.9900 / 0.9947 | 0.9888 / 0.9875 / 0.9947 / 0.9891 |
| schoof2021 MCAR | 1.0146 / 1.0044 / 0.9979 / 0.9968 | 1.0123 / 1.0035 / 0.9972 / 1.0008 |
| schoof2021 MNAR | 0.9835 / 0.9940 / 0.9869 / 0.9862 | 0.9727 / 0.9887 / 0.9862 / 0.9881 |

**Read for the king.**
- ψ = ∞ ties selected ψ on the stacked mean (0.9607 vs 0.9609), with a slightly better MCAR (0.9704 vs 0.9714) and a
  slightly worse MNAR (0.9509 vs 0.9504).
- Its *raw* output is clearly worse (0.9653 vs 0.9605), and the TMT targets are where removing per-protein
  standardisation hurts most (woo2025_PC3 raw MCAR 1.0404). The stack rescues it.
- For a preregistration, the parameter-free arm is defensible only as a stacked component. It does not fix the
  woo2025_PC3 MCAR > 1.00 problem (1.0009).

## Idea (the analogy)

**MLLR speaker adaptation.** Leggetter & Woodland (1995), "Maximum likelihood linear regression for speaker adaptation
of continuous density hidden Markov models", *Computer Speech & Language* 9(2):171–185, doi:10.1006/csla.1995.0010.
The citation was given in the directive as verified via Crossref (api.crossref.org/works/10.1006/csla.1995.0010); this
worker did not re-query it.

A speaker-independent model with thousands of Gaussian means is adapted to a new speaker with very little data. Instead
of adapting each mean separately (MAP), MLLR moves all of them with one shared linear transform, estimated from all the
adaptation frames pooled together. Mapped onto our problem:

| speech | here |
|---|---|
| speaker-independent model | the reference loadings Vref (~21 datasets, leave-target-out) |
| Gaussian means | per-protein loading rows Vref_j |
| per-Gaussian MAP | the champion's per-protein pull toward Vref_j (lam_v = 100) |
| the shared MLLR transform | A, an r0 × r0 matrix (r0 = 20, so 400 parameters) shared by all proteins |

As in speech practice (MLLR, then MAP), each protein keeps its MAP deviation D_j, which is now pulled toward
Vref_j A instead of Vref_j.

## Mechanism (as built)

**`adapt_project_mllr`** copies `nfbase.adapt_project_mb` line for line (same Z, ok, Vp, V init, U, rng and solve
order), then adds `A = I_r0` and `Vr = where(usable, Vref, 0)`. Each of the 10 iterations runs:

1. **U-step**, unchanged.
2. **A-step**, only when κ is finite.
   - D = V1 − Vr A.
   - R = Z − U2 Veᵀ − U1 Dᵀ on obs & ok.
   - Model: R_ij ≈ vr_jᵀ A u1_i. Using the row-major vec(A):
     - H = Σ_i kron(G_i, u1_i u1_iᵀ), with G_i = Vr[o_i]ᵀ Vr[o_i];
     - b = vec(Vrᵀ Rᵀ U1);
     - solve (H + κ·n_usable·I) vec(A) = b + κ·n_usable·vec(I).
   - In the `diag` arm, the same equations are restricted to the diagonal of A.
3. **Prior update:** Vp[:, :r0] = Vr A, then the unchanged V-step.

With κ = ∞, steps 2 and 3 are skipped and the code path is exactly `adapt_project_mb`.

**`run`:**
- ψ is cshrink's own selection: `C._score` with the same score cache (tag 'c'), joint ψ on cshrink's grid, ties to the
  smaller ψ. Then mu, sd = `cshrink_stats(M, ψ, ψ)`.
- κ ∈ {0.1, 1, 10, 100, ∞} is selected on `tgtshrink.selection_mask(M, seed)` (50/50 MCAR/MNAR, rng seed + 31337):
  - mu/sd come from cshrink_stats(Msel);
  - fits use iters = 5;
  - score = 0.5 (RMSE_MCAR + RMSE_MNAR);
  - the κ = ∞ score is reused from cshrink's score (equality is asserted in the unit check);
  - take the argmin, with ties going to the larger κ.
- Refit on M with iters = 10.

**Per-call log:** psi, κ scores, chosen κ, and for the final fit ‖A − I‖_F, the singular values, the real parts of the
eigenvalues, trace(H)/r0², and the polar parts ‖Q − I‖ and ‖S − I‖. The same diagnostics are logged for each finite-κ
selection fit.

**Deviation from the spec:** the spec asked for the spectrum "from a diagnostic refit at each finite κ". To save compute,
the logged per-κ diagnostics come from the selection fits on Msel (iters = 5), not from extra refits on M.

**Entry points:**

| function | what it does |
|---|---|
| `impute` | main arm |
| `impute_null` | protein-permuted A-step |
| `impute_base` | κ = ∞ |
| `impute_diag` | diagonal A |
| `impute_cshrink` | `cshrink.run(..., "main", log=False)`, i.e. `cshrink.impute` without writing to cshrink's log files; used as `mllr_base`'s null slot |

### Identifiability (found while building; this explains the result)
- The per-cell factors U are free and their ridge (lam = 5) is isotropic. So the objective is invariant under
  A → A Q, u_i → Qᵀ u_i, D → D Q for any orthogonal Q.
- Only the ‖A − I‖ ridge breaks this, so the **rotation part of A is not identified by the data at all**. Only A Aᵀ
  (a re-weighting of the reference factors relative to the U and D ridges) can be learned.
- That is why the polar diagnostics below show ‖Q − I‖ ≈ 0 and all the movement in S.
- The synthetic full-fit check confirms it: with κ tiny, A drifts far from A0 (rel. error 2.52), while the A-step with
  known U recovers A0 exactly.

## Null control
- **`impute_null`** keeps the whole machinery: the A-step, the κ grid, inner-mask selection, the ridge and the refit.
  The only change is in the A-step, both in G_i and b and during selection: each usable protein's reference row is
  replaced by another usable protein's row, using the fixed permutation rng(seed + 4242). The residual R and the prior
  mean Vr A still use the true Vr. This destroys the protein correspondence between the target's covariance and the
  reference.
- **`mllr_base`** (κ = ∞) is the same code with the idea switched off.

## Unit checks (`data/cache/phasec/mllr/unit_check.json`)

| check | result |
|---|---|
| U1: `impute_base` vs `cshrink.impute`, MCAR seed 100 | max abs diff **0.0** on petrosius2023 and leduc2022 (pass); κ = ∞ selection score equals C._score (assert pass) |
| U2: κ = 1e12 vs κ = ∞ (petrosius2023) | max abs diff 1.07e-13 (pass) |
| U3: synthetic A-step with known U recovers a rotation-times-scaling A0 (p = 1500, r0 = 20, 30 cells, 30% missing) | rel. Frobenius error 0.0175 at noise sd 0.1 (pass). A first run at noise sd 0.3 gave 0.0526 (just over 5%; it scales with the noise, so it is estimation noise, not a bug), and the noise was lowered. |
| U3 (info only): full alternating fit, κ = 1e-3 | rel. error of A 2.52, of A Aᵀ 3.64: U absorbs rotation and scale (see Identifiability) |

## Results: fast pool (6 targets × 2 masks × 2 seeds)

All numbers are `rmse_all` `mean_rel_to_champion` from the JSONs. Per-mask values are means of the 6 per-target rows.

| run | cand_stack | in-cx | cand_solo | cand_raw | null_stack | MCAR raw / stack | MNAR raw / stack | ≤ 0.99 | weight | runtime |
|---|---|---|---|---|---|---|---|---|---|---|
| `mllr` (main; null = permuted) | **0.9552** | 0.9531 | 0.9546 | 0.9538 | **0.9538** (in-cx 0.9528) | 0.9715 / 0.9671 | 0.9362 / 0.9433 | 12/12 | 0.792 | 1673 s |
| `mllr_base` (κ = ∞; null = cshrink.impute) | 0.9554 | 0.9533 | 0.9548 | 0.9539 | 0.9554 | 0.9713 / 0.9670 | 0.9365 / 0.9437 | 12/12 | 0.792 | 214 s |
| `mllr_diag` (diagonal A; null = base) | 0.9550 | 0.9530 | 0.9545 | 0.9537 | 0.9554 | 0.9713 / 0.9669 | 0.9360 / 0.9432 | 12/12 | 0.792 | 253 s |
| reference: `cshrink_fast` | 0.9554 | 0.9533 | 0.9548 | 0.9539 | 0.9637 | 0.9713 / 0.9670 | 0.9365 / 0.9437 | 12/12 | 0.792 | 236 s |

The null's per-mask values are MCAR stack 0.9635 and MNAR stack 0.9441.

- `mllr_base` reproduces cshrink_fast to the fourth decimal on every per-target row (bar: within 0.001).
- **Bars** (≤ 0.99; ≥ 0.005 below both the null and the base; MCAR raw no worse than the base by > 0.005):
  - main: 0.9552 is 0.0002 below the base and 0.0014 *above* the null. **Fails.**
  - diag: 0.9550 is 0.0004 below the base (its null slot). **Fails.**
  - MCAR raw guard: 0.9715 vs 0.9713 for main, so the guard holds.
- **Per target, main vs null (cand_stack):**
  - The null wins clearly on gregoire2023_mixCTRL MCAR (0.9542 vs 0.9668), sadiku2025 MCAR (0.9633 vs 0.9679) and
    px063590_ADbrain MCAR (0.9765 vs 0.9811).
  - Everything else ties within 0.002.
- Full pool was **not run**: the fast bars failed.

### Diagnosis per target (fast pool, all 8 calls per target: 2 masks × 2 seeds × {eval, inner} input)
Columns are means of the final fit. The last column is the mean relative inner-mask gain of the best κ over κ = ∞.

| target | arm | chosen κ (counts) | ‖A − I‖_F | top σ(A) | bottom σ(A) | ‖Q − I‖ (rotation) | ‖S − I‖ (re-weight) | trace(H)/r0² (median) | inner gain |
|---|---|---|---|---|---|---|---|---|---|
| petrosius2023 | main | 0.1:5, 1:2, ∞:1 | 0.128 | 1.002 | 0.903 | 0.025 | 0.126 | 1340 | 0.102% |
| leduc2022 | main | ∞:3, 0.1:3, 1:2 | 0.010 | 1.002 | 0.991 | 0.001 | 0.010 | 102 | 0.003% |
| weiss2026_liver | main | 1:6, 0.1:2 | 0.121 | 1.020 | 0.890 | 0.025 | 0.118 | 3213 | 0.138% |
| gregoire2023_mixCTRL | main | ∞:3, 0.1:4, 1:1 | 0.093 | 1.001 | 0.908 | 0.014 | 0.092 | 200 | 0.034% |
| px063590_ADbrain | main | ∞:5, 0.1:2, 1:1 | 0.012 | 1.000 | 0.989 | 0.001 | 0.012 | 103 | 0.002% |
| sadiku2025 | main | 0.1:3, 1:4, ∞:1 | 0.142 | 1.032 | 0.877 | 0.025 | 0.140 | 2573 | 0.087% |
| petrosius2023 | null | ∞:5, 100:3 | 0.015 | 1.001 | 0.988 | 0.005 | 0.015 | 1065 | 0.003% |
| leduc2022 | null | 1:6, 10:1, ∞:1 | 0.274 | 1.003 | 0.785 | 0.058 | 0.269 | 129 | 0.011% |
| weiss2026_liver | null | ∞:7, 100:1 | 0.016 | 1.000 | 0.984 | 0.004 | 0.016 | 3271 | 0.000% |
| gregoire2023_mixCTRL | null | 0.1:7, 1:1 | 2.487 | 1.025 | 0.005 | 3.033 | 2.445 | 7933 | 0.862% |
| px063590_ADbrain | null | 0.1:6, 1:2 | 1.611 | 1.014 | 0.216 | 0.599 | 1.597 | 3212 | 0.520% |
| sadiku2025 | null | ∞:6, 0.1:1, 1:1 | 0.531 | 1.007 | 0.753 | 0.435 | 0.521 | 15622 | 0.140% |

**TMT/carrier targets (diagnostic only, `python -m scp_graphfs_mcp.phasec.mllr --tmtdiag`).** These are not in the
fast pool. The diagnostic runs one call per mask at seed 100 (split recipe `nfbase._split`), outside the harness, so
there are no RMSE ratios.

| target-mask | ψ | chosen κ | ‖A − I‖_F | σ range | ‖Q − I‖ / ‖S − I‖ | trace(H)/r0² | inner score κ = ∞ → best |
|---|---|---|---|---|---|---|---|
| woo2025_PC3 MCAR | 2 | ∞ | 0.000 | 1.0–1.0 | 0 / 0 | — | 0.93666 → 0.93666 |
| woo2025_PC3 MNAR | 2 | 1 | 0.001 | 0.999–1.0 | 0.000 / 0.001 | 51 | 0.93538 → 0.93538 |
| williams2020_tmt MCAR | 1 | 1 | 0.001 | 0.999–1.0 | 0.000 / 0.001 | 111 | 0.74168 → 0.74168 |
| williams2020_tmt MNAR | 1 | 1 | 0.001 | 0.999–1.0 | 0.000 / 0.001 | 109 | 0.69713 → 0.69712 |
| schoof2021 MCAR | 16 | 0.1 | 0.157 | 0.853–1.0 | 0.018 / 0.156 | 427 | 0.96725 → 0.96667 |
| schoof2021 MNAR | 8 | 0.1 | 0.073 | 0.936–1.001 | 0.006 / 0.073 | 204 | 0.95294 → 0.95277 |

## Honest read
1. **A learns almost nothing on real data.**
   - ‖A − I‖_F ≤ 0.14 on average, and σ(A) stays in [0.88, 1.03].
   - All of the movement is a mild *down-weighting* of some reference factors (the symmetric part S). There is no
     rotation: ‖Q − I‖ ≤ 0.025.
   - The rotation half of the hypothesis cannot work in this model: rotation is non-identifiable once U is free (see
     Identifiability).
   - The re-weighting half is mostly already absorbed by the per-cell U-step (hundreds to thousands of observed proteins per cell
     against a ridge of 5) and by the per-protein D. The inner-mask gain is ≤ 0.14% on every target.
2. **The TMT/carrier prediction is falsified.**
   - On woo2025_PC3 and williams2020_tmt, A stays at I to three decimals.
   - trace(H)/r0² is small there (51–111), meaning the reference factors carry little target signal (small u1), so
     there is nothing to re-weight.
   - Only schoof2021 moves (‖A − I‖ 0.16 / 0.07), like the DIA targets.
3. **The null beats the candidate, and that is informative.**
   - With a permuted A-step and small κ, A collapses some reference directions toward zero (bottom σ 0.005 on
     gregoire2023_mixCTRL, 0.216 on px063590_ADbrain).
   - This *switches off* part of the reference prior, and it improves MCAR on gregoire2023_mixCTRL (0.9542 vs 0.9668),
     sadiku2025 and px063590_ADbrain.
   - So the transferable signal here is "on some targets the reference pull is too strong in some directions". It is
     not "the reference subspace needs re-rotation".
   - That is a generic attenuation effect, already in the spirit of ebshrink (per-protein lam_v, dead). It suggests a
     per-target scalar or per-factor strength of the reference pull, selected on the inner mask, is worth one cheap
     test.
4. **No revision (v2) was run.** The failure is not a bug or a scale problem:
   - at κ = 0.1 the ridge (κ·n_usable ≈ 270) is already well below the data curvature (trace(H)/r0² 100–3,200), and
     A still stays near I;
   - the unit checks pass.

## Next step (for the king)
- **Close the shared-linear-transform route.** The diagonal arm (pure factor re-weighting) adds 0.0004.
- **Cheapest follow-up the null points at.** A per-target attenuation of the reference prior, i.e. lam_v or a scalar
  shrink c·Vref with c selected on the 50/50 inner mask, credited against cshrink. Its null would be the same selection
  on a permuted inner mask. Expect a gain concentrated on MCAR for gregoire2023_mixCTRL, sadiku2025 and
  px063590_ADbrain.
- **Cross-pollination notes (from the proposal; untested):**
  - **bio:** MLLR regression-class trees become one A per curated protein class (compartment or complex), with random
    classes as the null. Given the result above, expect a small effect unless the classes have very different factor
    strengths.
  - **rep:** feature-space (constrained) MLLR transforms the target's Z before projection. It pairs with lwkrige's
    covariance, but suffers the same identifiability issue unless the cell factors are constrained.
  - **abd, next round:** speaker-adaptive training, i.e. one A_d per reference dataset when refitting Vref so the
    reference becomes technology-neutral. The citation must be verified before use; it is not cited here.

## Data provenance
None. The idea is data-only. Internal inputs: the leave-target-out reference loadings (`data/cache/phasec/<target>_V.npz`)
and cshrink's inner-mask score cache (`data/cache/phasec/cshrink/scores/`, read and written with the same exact keys).
Nothing was downloaded, so there is no PROVENANCE row.
