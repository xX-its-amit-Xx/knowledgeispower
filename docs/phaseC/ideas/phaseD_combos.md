# Phase D (dev only): stacking the champion and its Phase-C revision with bpca_oba and gcn_nograph

- **Status:** dev screen, finished (24/24 splits). **DEV ONLY.** No holdout is involved: nothing here loads, touches or
  scores a holdout dataset or seed. The targets are the Phase-C fast pool, and the harness asserts that no holdout
  target is in it. No verdict changes, and nothing is preregistered.
- **Code:** `benchmarks/phased_dev.py`
- **Result:** `results/phaseD/dev_combos.json`. Every number below comes from this file.
- **Cache** (not committed): `data/cache/phaseD/parts/`. It holds one file per (split, component) with the scored
  values only, so a run can resume. cshrink's psi-score cache is redirected to `data/cache/phaseD/cshrink/`, and the
  frozen `data/cache/phasec/cshrink/` is not written.

## Question
v9 found two data-only comparators that roughly match the H15 champion on dev: bpca_oba (the addendum's Oba-orientation
BPCA) and gcn_nograph (the GCN with no edges, i.e. a per-protein MLP). Does stacking them with the champion's components
(kNN, low-rank, adapt) or with the Phase-C revision (cshrink) help?

## Protocol
- **Splits:** `phasec.FAST` (6 targets) x `FAST_SEEDS` (100, 101) x {MCAR, MNAR} = 24 splits, 30 cells each. Cells
  are subsampled with `explore_phaseb.subsample` using the labels. The recipe is exactly that of
  `panel_v9.evaluate_split`: a 10% evaluation mask, plus a 10% inner mask (rng seed + 7919) that is used only to pick
  the stack weights.
- **Components:** each one is computed on both the evaluation input and the inner input.
  - kNN: `phase2.loo_knn` fill.
  - low-rank: `grmf`, lam_g = 0.
  - adapt: `transfer.adapt_project`, lam_v = 100, extra = 10, with V / warm from `phasec.target_models`.
  - cshrink: `cshrink.run(M, ctx, 'phaseD', log=False)`.
  - `baselines_addendum.bpca_oba`.
  - `baselines.gcn_nograph`.
- **Stacking:**
  - Simplex weights on a 0.25 grid for every combination of 2 to 6 components. These are the `explore_phaseb.GRID`
    tables for 2 to 4 components, and the same construction extended to 5 and 6. The **6-component stack uses the 0.25
    grid**.
  - An unconstrained non-negative least-squares stack of the same 6 components on the inner mask is reported as an
    extra row (`…__nnls`).
- **Metric:** RMSE relative to the champion = stack(kNN, low-rank, adapt). The headline is the mean, over the 12 target
  x mask cells, of the seed-averaged RMSE ratio (the `phasec.summarize` convention).
- **Reproduction check:** on all 24 splits, the champion, revised, bpca_oba+stack and gcn_nograph+stack RMSEs equal
  those in the frozen v9 dev panel (`results/v9/dev_splits/`) and in the v9 addendum dev cache, with a difference of
  exactly 0.

## Results (RMSE, all masked entries; < 1 is better)

| method (stack components) | rel. to champion | splits ≤ 0.99 / < 1 / > 1 of 24 | rel. to revised | splits < 1 / > 1 vs revised | in-complex rel. to champion | mean stack weights |
|---|---|---|---|---|---|---|
| champion (kNN, lr, adapt) | 1.0000 | – | 1.0477 | 0 / 24 | 1.0000 | kNN 0.198, lr 0.021, adapt 0.781 |
| revised (kNN, lr, adapt, cshrink) | 0.9554 | 24 / 24 / 0 | 1.0000 | – | 0.9533 | kNN 0.125, lr 0, adapt 0.083, cshrink 0.792 |
| bpca_oba+stack (kNN, lr, bpca) | 0.9916 | 8 / 15 / 9 | 1.0386 | 2 / 22 | 0.9985 | **bpca 0.906** (>0 in 24/24) |
| gcn_nograph+stack (kNN, lr, gcn) | 0.9917 | 10 / 14 / 10 | 1.0390 | 2 / 22 | 1.0004 | **gcn 0.865** (>0 in 24/24) |
| kNN, lr, adapt, bpca | 0.9809 | 14 / 22 / 2 | 1.0275 | 2 / 22 | 0.9848 | adapt 0.354, **bpca 0.594** (24/24) |
| kNN, lr, adapt, gcn | 0.9813 | 13 / 20 / 4 | 1.0280 | 2 / 22 | 0.9860 | adapt 0.396, **gcn 0.563** (24/24) |
| kNN, lr, cshrink, bpca | 0.9470 | 24 / 24 / 0 | 0.9912 | 20 / 3 | 0.9454 | cshrink 0.604, **bpca 0.385** (23/24) |
| kNN, lr, adapt, cshrink, bpca | **0.9469** | 24 / 24 / 0 | **0.9911** | 20 / 3 | **0.9451** | adapt 0.010, cshrink 0.604, **bpca 0.375** (23/24) |
| kNN, lr, adapt, cshrink, bpca, gcn (0.25 grid) | 0.9474 | 24 / 24 / 0 | 0.9920 | 18 / 6 | 0.9473 | cshrink 0.521, **bpca 0.219** (14/24), **gcn 0.250** (18/24) |
| same 6, NNLS (extra row) | 0.9464 | 24 / 24 / 0 | 0.9907 | 18 / 6 | 0.9458 | cshrink 0.583, bpca 0.186, gcn 0.216 |
| raw cshrink alone | 0.9539 | 23 / 24 / 0 | 0.9983 | 10 / 9 | 0.9508 | – |
| raw bpca_oba alone | 0.9918 | 9 / 14 / 10 | – | – | 0.9979 | – |
| raw gcn_nograph alone | 0.9929 | 10 / 13 / 11 | – | – | 1.0014 | – |

Descriptive split-level intervals: geometric mean of per-split ratios, t 95% interval over 24 splits. These are
optimistic, because the splits share 6 targets.
- (kNN, lr, adapt, cshrink, bpca) vs revised: 0.9911 [0.9869, 0.9953]. No split is ≥ 1.01 worse.
- The 6-component grid stack vs revised: 0.9917 [0.9831, 1.0005].
- bpca_oba+stack vs champion: 0.9912 [0.9787, 1.0039].
- gcn_nograph+stack vs champion: 0.9912 [0.9798, 1.0027].

### Per target x mask (rel. to champion, seed-averaged)

| target | mask | revised | bpca+stk | gcn+stk | kNN,lr,ad,bpca | kNN,lr,ad,gcn | kNN,lr,cs,bpca | kNN,lr,ad,cs,bpca | all 6 (grid) |
|---|---|---|---|---|---|---|---|---|---|
| gregoire2023_mixCTRL | MCAR | 0.9668 | 0.9414 | 0.9665 | 0.9414 | 0.9653 | 0.9390 | 0.9390 | 0.9390 |
| gregoire2023_mixCTRL | MNAR | 0.9325 | 0.9386 | 0.9548 | 0.9386 | 0.9570 | 0.9206 | 0.9206 | 0.9223 |
| leduc2022 | MCAR | 0.9778 | 1.0044 | 1.0155 | 0.9986 | 1.0025 | 0.9782 | 0.9782 | 0.9806 |
| leduc2022 | MNAR | 0.9688 | 1.0032 | 1.0034 | 0.9985 | 0.9968 | 0.9745 | 0.9745 | 0.9749 |
| petrosius2023 | MCAR | 0.9692 | 1.0062 | 1.0068 | 0.9857 | 0.9840 | 0.9591 | 0.9591 | 0.9610 |
| petrosius2023 | MNAR | 0.9668 | 1.0130 | 0.9952 | 0.9890 | 0.9869 | 0.9599 | 0.9587 | 0.9595 |
| px063590_ADbrain | MCAR | 0.9812 | 0.9925 | 0.9952 | 0.9921 | 0.9924 | 0.9738 | 0.9738 | 0.9711 |
| px063590_ADbrain | MNAR | 0.9800 | 0.9955 | 1.0027 | 0.9940 | 0.9957 | 0.9774 | 0.9774 | 0.9774 |
| sadiku2025 | MCAR | 0.9676 | 1.0085 | 0.9704 | 0.9877 | 0.9599 | 0.9618 | 0.9618 | 0.9362 |
| sadiku2025 | MNAR | 0.9311 | 0.9939 | 0.9531 | 0.9769 | 0.9426 | 0.9101 | 0.9101 | 0.9019 |
| weiss2026_liver | MCAR | 0.9397 | 1.0254 | 0.9952 | 0.9959 | 0.9794 | 0.9357 | 0.9357 | 0.9233 |
| weiss2026_liver | MNAR | 0.8828 | 0.9767 | 1.0412 | 0.9728 | 1.0130 | 0.8743 | 0.8743 | 0.9216 |

## Honest read
1. **Adding bpca_oba to the revision helps a little, consistently, on dev.**
   - (kNN, lr, adapt, cshrink, bpca) is 0.9469 vs the champion, against 0.9554 for revised, i.e. 0.9911 of revised.
   - It is better than revised on 20/24 splits and ≥ 1% better on 7/24. It is worse on 3, none of them by ≥ 1%.
   - bpca_oba takes a mean weight of 0.375 and is used on 23/24 splits.
   - The gain is about 0.9% RMSE. That is below the 1% SESOI that has been used throughout this project, and it is
     measured on the same fast pool that cshrink was selected on, so it is optimistic.
2. **adapt (the reference-based component) adds nothing once cshrink and bpca are in.** Its mean weight is 0.010, and
   (kNN, lr, cshrink, bpca) is 0.9470 vs 0.9469 with adapt. Only the reference-informed path inside cshrink is left: it
   still uses V / warm via `adapt_project_mb`.
3. **gcn_nograph adds nothing on top of bpca + cshrink.**
   - On the 0.25 grid the 6-component stack is 0.9474, slightly worse than the 5-component stack's 0.9469. Its split
     interval vs revised reaches 1.0005, and it is worse than revised on 6 splits (vs 3 for the 5-component stack).
   - NNLS on the 6 components gives 0.9464, a 0.0005 difference from the best 5-component grid stack. That is noise
     level, and NNLS is not sum-constrained.
   - The weight that gcn_nograph gets (mean 0.250 on the grid) mostly displaces bpca. bpca's weight drops from 0.375
     to 0.219.
4. **Without cshrink, bpca_oba and gcn_nograph mostly substitute for adapt rather than complement it.**
   - Each takes the majority weight with the champion's components: bpca 0.594, gcn 0.563.
   - They improve on the champion by about 1.9% (0.9809 / 0.9813), with 2 to 4 splits worse.
   - Alone in a stack they are about equal to the champion: 0.9916 / 0.9917, worse on 9 to 10 splits, with intervals
     crossing 1. On the in-complex endpoint they are about 1.000.
5. **cshrink alone (raw) is 0.9539, about the same as revised (0.9554).** The revision's gain is essentially the
   cshrink encoding, not the stack. The stack's main job is to blend in bpca.
6. **For the thesis:** none of this involves a curated complex graph. bpca_oba and gcn_nograph are graph-free, and the
   gains are data-only or reference-only. In-complex RMSE moves in step with all-entry RMSE (0.9451 vs 0.9533 for
   revised), so there is no complex-specific signal.
7. **Recommendation:** if anything is carried forward, it is a single candidate, "revised + bpca_oba"
   (kNN, lr, adapt, cshrink, bpca). Before any claim, it would need the full dev pool (FULL_SEEDS, all dev targets) and
   then a preregistration. Its expected gain over revised is below 1% and may not clear the SESOI.
