# lwkrige: kriging with a locally re-estimated covariance, as an exact correction on cshrink

- **Queen:** representation & mathematics (round 5; directive in the latest `KING_LOG.md` entry, abd's kriging analogy)
- **Status (worker's recommendation):** **dead at the fast bar, but with a real, null-separated signal.**
  - v1 (as specified, pairwise-complete S_target) never moves off d = 1, so it equals cshrink exactly.
  - v2 (the one allowed revision: a PSD S_target from the base-completed block) gives fast cand_stack **0.9522**.
    That is 0.0032 below the base (0.9554) and 0.0031 below its permuted-S null (0.9553), both under the 0.005 bar.
    The full pool was not run.
  - The gain sits on the TMT targets, as predicted. The technology-matched Sigma_ref arm is dead: it loses to the base
    on every call, and it does worse than its all-reference control.
- **Code:** `src/scp_graphfs_mcp/phasec/lwkrige.py`
- **Results** (in `results/phaseC/`): `lwkrige_fast.json` (v1), `lwkrige_base_fast.json`, `lwkrige_d85_fast.json`
  (v1 diagnostic), `lwkrige_v2_fast.json`, `lwkrige_v2_d98_fast.json`, `lwkrige_tech_fast.json`
- **Cache** (in `data/cache/phasec/lwkrige/`, not committed):
  - `unit_check.json`;
  - per-call logs: `v1_log_main.jsonl` / `v1_log_null.jsonl` (v1), and `log_<variant>.jsonl` for base, d85, d85_null,
    v2, v2_null, v2_d98, v2_d98_null, tech and allraw;
  - `<target>_Vtech.npz` (technology-matched loadings);
  - `smoke/` (smoke-test logs; the smoke result JSON was deleted).

## Idea
- **Base.** cshrink imputes a cell's holes through the covariance implied by its adapted factor model,
  Sigma_ref = V_ad V_ad^T + 5 I. The conditional-Gaussian mean under this matrix is exactly the U-step ridge posterior
  (unit check U4).
- **Change.** Kriging with a locally re-estimated variogram, i.e. Ledoit-Wolf / Schafer-Strimmer linear shrinkage toward
  a structured target. It blends the target's own 30-cell covariance into Sigma_ref:
  Sigma_d = (1 - d) S_target + d c Sigma_ref.
- **Selection.** d is chosen on the 50/50 MCAR/MNAR inner mask.
- **Crediting.** The result is applied as an exact correction: Y = Yc + sd (zhat_d - zhat_1). So d = 1 is cshrink
  bit for bit.

## Mechanism (as built)
- **Step 0 (base).** psi is selected exactly as `cshrink.run` selects it: same inner mask and same on-disk score cache
  (`cshrink._score`, tag `c`). Then (mu, sd) = `cshrink_stats(M, psi, psi)`, and
  (Yc, U, Vad) = `adapt_project_uv` (a verbatim copy of `nfbase.adapt_project_mb` that also returns U and V;
  lam_v = 100, extra = 10, lam = 5, iters = 10).
- **Step 1 (block).** B = the 1,500 ok proteins with the most observed cells (stable argsort). Z = (M[:, B] - mu) / sd.
  - **v1:** S_target = pairwise-complete covariance (each pair centred on its common cells). Entries with < 8 common
    cells are filled from c Sigma_ref.
  - **v2:** the holes of Z are filled with the base's conditional mean (zhat_1, the cshrink projection), then
    S_target is the ordinary sample covariance. It is PSD with rank <= 29.
  - In both: c = sum of valid diag(S) / sum of the same diag(Sigma_ref).
- **Step 2 (blend).** Sigma_d = (1 - d) S + d c Sigma_ref. For d < 1, eigenvalues are clipped at 1e-3 × their mean.
- **Step 3 (conditional mean).** Dense Cholesky per cell: zhat_d[i, m] = Sigma_d[m, o] Sigma_d[o, o]^-1 z[i, o].
  d = 1 runs through the same path.
- **Step 4 (output).** Y = Yc, plus sd_j (zhat_d - zhat_1) on the missing entries in B. Every other entry keeps Yc.
- **Step 5 (selection).** d is selected on `tgtshrink.selection_mask(M, seed)`, rebuilding steps 0-4 on Msel with the
  psi selected on M and iters = 5. Score = 0.5 × (RMSE_MCAR-part + RMSE_MNAR-part). Ties go to the larger d. The final
  fit is on M.
  - Grid v1 = {0.5, 0.7, 0.85, 0.95, 1}.
  - Grid v2 = {0.9, 0.95, 0.98, 0.99, 1}.
- **Diagnostic only.** The Schafer-Strimmer analytic d* (toward c Sigma_ref, on the zero-filled standardised block) is
  logged and never used.
- **Pluggable for mod's `linerna`.** `krige_correct(M, ctx, Yc, mu, sd, Vad, sigma_ref_fn=None, d_grid=..., S_fn=None,
  s_est=...) -> (Ys, info)` returns 'base' and every d. `lw_run(...)` runs the full select-and-refit pipeline with the
  same plug-ins. `sigma_ref_fn(B, ctx)` returns a |B| × |B| PSD matrix.

**Unit checks** (`data/cache/phasec/lwkrige/unit_check.json`, petrosius2023 MCAR seed 100; psi = 16):

| Check | Result |
|---|---|
| U1: `adapt_project_uv` vs `adapt_project_mb` | max diff 0.0 |
| U2: Yc vs `cshrink.run` | bit-exact (true) |
| U3: `impute_base` vs cshrink | bit-exact (true) |
| U4: conditional mean under Vad Vad^T + 5 I vs explicit ridge U-step on B | max diff 2.3e-14 |

## Null controls
- **impute_null / impute_v2_null.** The same B, c, blend, PSD repair, grid, inner-mask selection and conditioning
  values. Only S_target changes: it comes from Zp, in which each protein column of Z (values and NaN pattern together)
  is permuted across cells with rng(seed + 4242 + global column index). Marginals and observed counts are kept. In v2,
  Zp is completed by the same base conditional mean before the covariance.
- **lwkrige_base.** d = 1, which equals cshrink by construction.
- **impute_allraw.** The content control for the technology arm: the same raw-loading Sigma_ref path, built from all
  reference datasets.

## Data provenance
No external data, so there is no PROVENANCE row. Inputs: the target matrix, `ctx['V']`, `ctx['warm']`,
`ctx['ref']()`, `ctx['ref_src']()` and `ttprior.META` (acquisition classes).

Citations (as supplied in the directive, verified there through Crossref on 2026-10-23; not re-verified by this
worker):
- Ledoit & Wolf 2004, *J. Multivariate Analysis* 88(2):365-411, doi:10.1016/S0047-259X(03)00096-4.
- Schäfer & Strimmer 2005, *Stat. Appl. Genet. Mol. Biol.* 4(1), doi:10.2202/1544-6115.1175.

## Results (fast pool: 6 targets × 2 masks × 2 seeds; rmse_all mean_rel_to_champion)

| run | cand_stack | cand_solo | cand_raw | null_stack | in-cx cand_stack | weight | splits ≤ 0.99 |
|---|---|---|---|---|---|---|---|
| `lwkrige_base` (d = 1) | 0.9554 | 0.9548 | 0.9539 | — | 0.9533 | 0.7917 | 12/12 |
| `lwkrige` v1 (pairwise S, grid 0.5–1) | 0.9554 | 0.9548 | 0.9539 | 0.9554 | 0.9533 | 0.7917 | 12/12 |
| `lwkrige_d85` v1, d forced 0.85 | 1.0 | 1.0772 | 1.5903 | 1.0 | 1.0 | 0.0 | 0/12 |
| **`lwkrige_v2`** (completed S, grid 0.9–1) | **0.9522** | 0.9516 | 0.952 | 0.9553 | 0.9494 | 0.7292 | 12/12 |
| `lwkrige_v2_d98` (v2, d forced 0.98) | 0.9524 | 0.9548 | 0.9541 | 0.9598 | 0.9495 | 0.7083 | 12/12 |
| `lwkrige_tech` (tech Sigma_ref, v2 S; null = allraw) | 0.9554 | 0.9548 | 0.9539 | 0.9554 | 0.9533 | 0.7917 | 12/12 |

`lwkrige_base` reproduces cshrink_fast (0.9554 / 0.9548 / 0.9539) exactly.

Per-mask means of the 6 per-target values (computed from `per_target_mask` in the JSONs):

| run | MCAR cand_raw | MCAR cand_stack | MCAR null | MNAR cand_raw | MNAR cand_stack | MNAR null |
|---|---|---|---|---|---|---|
| base | 0.9713 | 0.967 | — | 0.9365 | 0.9437 | — |
| v2 | 0.9708 | 0.964 | 0.9668 | 0.9333 | 0.9403 | 0.9438 |
| v2_d98 | 0.9741 | 0.9654 | 0.9713 | 0.9341 | 0.9395 | 0.9484 |

TMT_DDA (leduc2022, gregoire2023_mixCTRL, px063590_ADbrain) vs DIA (petrosius2023, weiss2026_liver, sadiku2025), per
`ttprior.META`. Each value is the mean of 6 target-masks:

| run | TMT cand_stack | TMT null | DIA cand_stack | DIA null |
|---|---|---|---|---|
| base | 0.9678 | — | 0.9429 | — |
| v2 | 0.9627 | 0.9676 | 0.9417 | 0.943 |
| v2_d98 | 0.9622 | 0.9738 | 0.9426 | 0.9458 |

Per target-mask, `lwkrige_v2` cand_stack (base → v2):

| target | MCAR | MNAR |
|---|---|---|
| gregoire2023_mixCTRL | 0.9668 → 0.9593 | 0.9325 → 0.9194 |
| leduc2022 | 0.9778 → 0.9747 | 0.9688 → 0.9646 |
| px063590_ADbrain | 0.9812 → 0.9757 | 0.9800 → 0.9822 (worse) |
| petrosius2023 | 0.9692 → 0.9700 (worse) | 0.9668 → 0.9655 |
| sadiku2025 | 0.9676 → 0.9667 | 0.9311 → 0.9285 |
| weiss2026_liver | 0.9397 → 0.9378 | 0.8828 → 0.8816 |

### Chosen d (per-call logs; 48 calls each = 24 splits × eval/inner input)

| variant | chosen d |
|---|---|
| v1 main | d = 1 in 48/48 |
| v1 null | d = 1 in 48/48 |
| v2 main | d = 0.99 in 40/48, d = 1 in 8/48 |
| v2 null | d = 1 in 45/48, d = 0.99 in 3/48 (px063590 2, weiss2026 1) |
| tech | base in 48/48 |
| allraw | base in 48/48 |

- **v2 main by target:**
  - gregoire2023_mixCTRL: 0.99 in 8/8.
  - leduc2022, weiss2026_liver, px063590_ADbrain, sadiku2025: 0.99 in 7/8 each.
  - petrosius2023: 4/8.
- **Inner-mask score gain at d = 0.99** (median of score(1) − score(0.99)): +0.00145 for v2 against −0.00184 for its
  null.
- **Tech vs allraw**, median inner-mask score relative to base:

  | d | tech | allraw |
  |---|---|---|
  | 0.99 | 1.0438 | 1.0250 |
  | 1.0 | 1.1700 | 1.1564 |

### Diagnostics (from the logs)

**v1 S_target is grossly indefinite.** Scratch check on Msel, seed 100:

| target / mask | negative eigenvalues | λ_min | λ_max |
|---|---|---|---|
| petrosius2023 MCAR | 39% | −88.9 | 837 |
| leduc2022 | 43% | −40 | 154 |

- After clipping, every d < 1 explodes. The v1 inner-mask score at d = 0.5 was 1.83 against 0.76 at d = 1
  (petrosius2023), and forced d = 0.85 gives cand_raw 1.5903.
- The median pairwise n is 9 in v1, and frac_filled (< 8 common cells) has a median of 0.339.

**Other diagnostics:**
- **c (median):** 0.149 in v1, 0.122 in v2.
- **Schafer-Strimmer analytic d* (median):**

  | variant | median d* | min |
  |---|---|---|
  | v2 | 0.477 | 0.184 |
  | v2 null | 0.718 | 0.665 |

  The analytic d* separates real from permuted S, but it is far from the inner-mask optimum (0.99). The zero-filled
  estimator ignores missingness, so it is not a usable selector here.
- **Technology-matched loadings:** fitted on 5 TMT_DDA or 9 DIA reference datasets. There was no fallback on any fast
  target (the dataset lists are in `<target>_Vtech.npz`).

## Honest read
- **Prediction (1), that the inner mask picks d < 1 on most splits:** v1 no (0/48). v2 yes (40/48 at 0.99), against
  3/48 for the null. The target's own covariance carries information the adapted factor model lacks, but only at a
  1% blend weight.
  - The rank-29 sample covariance from 30 cells is very large on its subspace compared with the reference structure
    (the trace scaling by c spreads it over 1,500 proteins).
  - So the useful blend is tiny. Selection never went to 0.98 or below.
- **Prediction (2), a margin of at least 0.005 below both base and null:** **fails.** The margins are 0.0032 vs base
  and 0.0031 vs null.
  - The forced-d diagnostic separates S's content from its machinery more clearly: v2_d98 0.9524 vs its null 0.9598, a
    margin of 0.0074.
  - But it still beats the base only by 0.0030.
  - MCAR cand_raw (0.9708) is not worse than the base (0.9713), so the MCAR guard holds.
- **Prediction (3), a larger gain on TMT:** **holds.**
  - TMT cand_stack moves 0.9678 → 0.9627 (−0.0051; null 0.9676).
  - DIA moves 0.9429 → 0.9417 (−0.0012; null 0.943).
  - gregoire2023_mixCTRL gains most (MNAR 0.9325 → 0.9194).
  - This matches the round-4 weak spot (TMT/carrier MCAR), but it is a fast-pool, 3-target observation.
- **Prediction (4), technology-matched Sigma_ref beats all-reference raw:** **fails.**
  - Raw-loading covariances (V V^T + diag(1 − ‖v‖²)) lose to the adapted base by 2.5–17% on the inner mask.
  - The technology-matched version is worse than the all-reference version, 1.0438 vs 1.0250 at d = 0.99.
  - Amit's technology seed is dead at the covariance level too. The adapted loadings (pulled to the target with 30
    cells) dominate any unadapted reference structure.
- **Overall:** lwkrige v2 is a real but small (~0.3%) data-only correction, concentrated on TMT targets. It sits below
  the 0.005 credit bar, so no full pool was run.

## Next step
- The signal says the in-domain covariance matters, but a raw 30-cell covariance is too noisy to blend directly. The
  natural follow-up is a *low-rank* in-domain correction: blend only the top k (2–5) eigenvectors of S_target's
  residual after the adapted model, with d on a finer grid near 1 (0.97–1.0).
  - This is close to abd's `mllr` (a shared re-rotation) and to adding free factors. Compare against `mllr` before
    spending a round on it.
- For `linerna`, use `lw_run(..., sigma_ref_fn=..., s_est="completed", d_grid=D_GRID_V2)`. Expect any external
  Sigma_ref to need d ≈ 0.99 too, and expect it to lose to the adapted Vad base if used in place of it (as tech and
  allraw did). An external covariance is better blended on the S side, alongside the adapted Sigma_ref, than used as Sigma_ref.
  That needs a small new hook, because `S_fn` currently returns the data matrix that S is computed from, not a
  covariance.
