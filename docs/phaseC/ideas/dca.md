# dca: mean-field direct-coupling analysis (shrinkage-inverse precision, Gaussian conditional mean)

- **Queen:** innovation abduction (round-1 directive `dca`)
- **Status (worker's recommendation):** **dead** at the fast screen.
  - cand_stack is 1.0004 (rmse_all) and 1.0006 (in-complex), which misses the 0.99 bar. The null_stack is 1.0.
  - The stack gave the candidate weight 0 in 23 of the 24 inner-mask fits; the mean weight was 0.0104.
  - The full pool was not run, and no v2 was built, because no bug or scale error turned up (see the honest read).
- **Code:** `src/scp_graphfs_mcp/phasec/dca.py`, with `impute`, `impute_null` and the diagnostic `couplings_diag`
  (`python -m scp_graphfs_mcp.phasec.dca`).
- **Results:**
  - `results/phaseC/dca_fast.json` is the screen.
  - `results/phaseC/dca_couplings_diag.json` is the optional CORUM diagnostic.
  - The single-target smoke run (`dca_smoke`) was deleted.

## Idea
Every model so far encodes protein co-regulation as a low-rank covariance. The Morcos et al. (PNAS 2011) mfDCA
estimator inverts a pseudocount-regularised covariance instead. In MSAs, that separates direct couplings from
transitive ones. Here the hypothesis was that a shrinkage-inverse precision matrix, fitted on the 21-dataset
leave-target-out single-cell reference, would carry direct pairwise couplings (complex subunits) that the champion's
20+10 factors drop. Imputing with E[z_m | z_o] = −Θ_mm⁻¹ Θ_mo z_o would then add information to the stack.

## Mechanism (as built, following the queen's spec)
1. **Reference statistics** (per target and variant; cached in `data/cache/phasec/dca/{target}_{real|null}.npz`;
   built under a per-target `fcntl` file lock).
   - Z = `ctx['ref']()`, the reference already z-scored within each dataset.
   - S = the 1,500 target proteins with the most reference observations, taken from the real reference.
   - C is the pairwise-complete correlation on Z[:, S], with pairs that have n < 30 set to 0, rescaled to unit
     diagonal. It is then made PSD by clipping eigenvalues at 1e-3·mean and rescaled to unit diagonal again.
2. **λ selection, inside the reference only.**
   - The largest reference dataset is held out. It was px063590_ADbrain for all targets except px063590_ADbrain
     itself, where it was schoof2021.
   - C is recomputed on the remaining datasets.
   - From the held-out dataset, 300 rows are taken and 20% of their observed S-entries are hidden (rng 0).
   - Each λ in {0.3, 0.5, 0.7, 0.85} is scored by the RMSE of the conditional mean on the hidden entries.
   - Θ = inv((1−λ)C + λI).
3. **Imputation.**
   - Xm is standardised exactly as in `transfer.adapt_project`.
   - Each cell with at least 5 observed S-proteins gets the conditional mean, computed through
     C_λ,mo C_λ,oo⁻¹ z_o or −Θ_mm⁻¹Θ_mo z_o, whichever system is smaller. Values are clipped to [−4, 4] and
     de-standardised.
   - Every other entry keeps `ctx['champion']`.
4. **Export (real variant), not used for imputation.** `{target}_couplings.npz` holds:
   - the partial correlations −Θ_ij/√(Θ_iiΘ_jj);
   - the APC-corrected F;
   - the top-30 partners per protein, with their accessions;
   - the coupling degree (the sum of the top-30 F).
   The null's couplings are written to `{target}_null_couplings.npz` for the diagnostic.

**Reference facts recorded in the cache `info`:**
- 21 reference datasets per target, with 9,366–11,889 rows.
- p = 1,878–5,140 target proteins, of which the 1,500 in S are modelled.
- Before the PSD repair, C had 418–463 negative eigenvalues (real) and 545–685 (null).

## Null control
`impute_null` is the identical pipeline: the same S, correlation, repair, in-reference λ selection, inverse and
solver. The one difference is that, before any statistic is computed, each reference dataset's protein columns are
permuted with an independent permutation (rng 12345 + the dataset's index in sorted order). Missingness moves with
the values. This destroys the cross-dataset agreement on which protein couples to which.

## Data provenance
There is no external data. The idea uses only the existing leave-target-out dev single-cell reference
(`ctx['ref']`, `ctx['ref_src']`), so no holdout is touched and no PROVENANCE row is needed. The cache in
`data/cache/phasec/dca/` is 420 MB.

## λ per target (from the cache `info` / build log; held-out-reference RMSE, z-units)
- λ = **0.85**, the upper edge of the grid, was chosen for every target and for both variants.
- The table compares the real-variant conditional mean (at 0.85) with predicting z = 0 on the same hidden entries.

| target | real λ=0.85 | real z=0 | null λ=0.85 | null z=0 |
|---|---|---|---|---|
| petrosius2023 | 0.9679 | 0.9757 | 0.9782 | 0.9769 |
| leduc2022 | 0.9629 | 0.9710 | 0.9879 | 0.9865 |
| weiss2026_liver | 0.9678 | 0.9751 | 1.0040 | 1.0029 |
| gregoire2023_mixCTRL | 0.9609 | 0.9707 | 0.9820 | 0.9809 |
| px063590_ADbrain | 0.9706 | 0.9930 | 0.9896 | 0.9900 |
| sadiku2025 | 0.9756 | 0.9870 | 0.9863 | 0.9840 |

- **Real reference:** across datasets, the conditional mean beats the column mean by only about 1–2%.
- **Null:** it is about equal to the column mean or slightly worse, as intended.

## Time per target
- **Statistics build** (`build_s` in the cache, both variants):
  - petrosius2023: 218.6 s (real) and 47.1 s (null). This was the cold run, which includes assembling the reference.
  - Other targets: 6–11 s per variant.
- **Fast screen:** `runtime_s` = 645.6 s for 24 splits (6 targets × 2 masks × 2 seeds, `--jobs 2`), which is about
  108 s per target.

## Results: fast pool (`results/phaseC/dca_fast.json`, mean_rel_to_champion)

| metric | rmse_all | rmse_in_complex |
|---|---|---|
| cand_stack | **1.0004** | **1.0006** |
| cand_solo | 1.032 | 1.0365 |
| cand_raw | 1.059 | 1.0676 |
| null_stack | 1.0 | 1.0 |
| stack_lowrank | 1.0773 | 1.0904 |
| knn | 1.1114 | 1.1312 |

- Splits with cand_stack ≤ 0.99: 0 of 12.
- mean_cand_weight_in_stack = 0.0104.
- All target × mask cells are exactly 1.0 except weiss2026_liver MNAR, at 1.0049. That is the one split where the
  inner mask gave the candidate weight, and the weight did not transfer.

## Diagnostic: do the couplings form a complex graph? (`results/phaseC/dca_couplings_diag.json`)
The metric is the fraction of each protein's top-30 partners that share a CORUM complex, among proteins that are in
CORUM (835–871 of 1,500 per target). The overall figures are means over the 6 targets.

| partners by | mean | range over targets |
|---|---|---|
| APC-corrected DCA F (real) | 0.0799 | 0.0662–0.0891 |
| APC-corrected DCA F (null) | 0.0177 | 0.0145–0.0226 |
| raw \|correlation\| (real) | 0.1063 | 0.0891–0.1264 |
| random partners | 0.0169 | 0.0143–0.0204 |

- The real couplings are about 4.5× enriched for shared complexes over random; the null is at the random level.
- **However, the DCA claim fails here.** Raw |correlation| recovers complexes *better* than the inverted,
  APC-corrected couplings, on every target.

## Honest read
- **Not a bug, and not additive.**
  - The candidate carries real signal: cand_solo is 1.032, well below stack_lowrank at 1.0773, and the real
    couplings are complex-enriched while the null's are not.
  - But that signal is strictly dominated by the champion's adapted projection. The stack weight is ~0, and
    cand_stack is about equal to null_stack, both at about 1.0.
- **The precision encoding adds nothing over the low rank, for three reasons.**
  1. **λ pins at 0.85.** The reference itself asks for heavy shrinkage toward I. Pairwise-complete, cross-lab
     single-cell covariance is noisy enough (more than 400 negative eigenvalues before repair) that the inverse
     must be mostly identity. What survives is close to the leading covariance structure, which the low rank
     already holds.
  2. **Direct coupling does not beat correlation.** Unlike MSAs, inverting did not sharpen complex recovery; raw
     correlation is better. In SCP, complex co-variation appears to be carried by the shared programmes, so it is
     not a sparse residual that the factors miss.
  3. **Only 1,500 proteins.** The candidate covers at most the 1,500 most-observed proteins. These are the ones the
     champion already predicts best.
- **Pushing λ beyond the grid would not help.** It moves the estimate toward the column mean. That is why I did not
  spend the one revision: no fixable cause was found.

## Next step
- **Retire `dca` as a stack member.** The planned low-rank + sparse follow-up (the same estimator on Z − U_ref Vᵀ)
  is the only variant that targets the residual directly. Given result 2 above (couplings ⊂ programmes), I would
  rank it low.
- **Cross-pollination value:**
  - The exported couplings are a complex-enriched, data-driven graph (8% vs 1.7% random). bio could test them as a
    graph penalty on V.
  - But raw |correlation| top-30 neighbours are even more enriched (10.6%), so that graph is the stronger candidate
    if bio wants one.
  - The per-protein coupling degree and diag(Θ) are in `{target}_couplings.npz` / `{target}_real.npz` for rep's
    per-protein λ_v,j.
- **Next queued abduction:** the Netflix-Prize baseline predictors (explicit per-reference-dataset × protein
  offsets).
