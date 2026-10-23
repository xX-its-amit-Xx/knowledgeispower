# Phase-2 exploration log (dev masks / dev CV seeds only)

Protocol: `benchmarks/EXPLORATION.md`. Every result goes here, including failures. Raw numbers are in
`results/phase2/dev/`. These are **development** numbers: they guide the choice of method and are
**not** evidence for the claim.

## 2026-10-21 · P1 + W1: panel design as sensor placement (issue #6)
Run with `python benchmarks/explore_panel.py` (CV seeds 100–102, all 215 cells). Values are mean held-out R² over targets.

| k | HVP | HVP + W1 | pivoted QR + W1 | random (10) | coverage CORUM + W1 | its null (3) | coverage data graph (matched) |
|---|---|---|---|---|---|---|---|
| 25 | 0.213 | 0.218 | 0.367 | 0.447 | 0.452 | 0.314 | 0.487 |
| 50 | 0.199 | 0.340 | 0.368 | 0.416 | 0.406 | 0.329 | 0.449 |
| 100 | 0.237 | 0.228 | 0.178 | 0.243 | 0.261 | 0.267 | 0.292 |

(Complex Portal coverage is close to CORUM's; see the JSON.)

**Read:**
- **W1 (keratin filter):** HVP puts 8 keratins in 25 slots; the filter fixes that.
- **Structure is real:** greedy coverage over curated complexes beats its own degree-preserving null by about 0.08–0.14 R² at k = 25–50.
- **Not better than data or random:** curated coverage only *ties* random panels, and coverage over a data-derived co-variation graph beats both.
- **Caveat:** the data graph is built on all cells, a mild leak in its favour.
- **Verdict on dev:** P1 does not qualify under the selection rule, because the data-derived graph wins.

## 2026-10-21 · Imputation ideas E1–E6, dev grid (issues #2–#5)
Run: `python benchmarks/explore.py --n-cells 0 60 30`. In-complex RMSE, mean over dev seeds 100–102.
Raw rows: `results/phase2/dev/dev_results.json`.

**All 215 cells.** Each cell shows MCAR / MNAR; the "data graph" column is the matched-edge data-derived co-variation graph.

| Method (curated graph) | Curated | Its null (3) | Data graph | kNN | Low-rank, no graph |
|---|---|---|---|---|---|
| E1 residual-on-kNN (Complex Portal) | 0.672 / 0.747 | 0.684 / 0.756 | 0.681 / 0.776 | 0.687 / 0.762 | 0.659 / 0.723 |
| E1+E4 gated (Complex Portal) | 0.671 / 0.746 | 0.684 / 0.754 | 0.681 / 0.776 | 〃 | 〃 |
| E2 stoich-on-kNN (Complex Portal) | 0.677 / 0.754 | 0.687 / 0.761 | 0.684 / 0.773 | 〃 | 〃 |
| E5 fused (Complex Portal) | 0.671 / 0.750 | 0.694 / 0.768 | 0.674 / 0.766 | 〃 | 〃 |
| E6 residual-on-low-rank (Complex Portal) | 0.656 / 0.720 | 0.659 / 0.723 | 0.658 / 0.722 | 〃 | 〃 |
| E3 low-rank + Laplacian (CORUM) | 0.640 / 0.720 | 0.641 / 0.721 | 0.654 / 0.723 | 〃 | 〃 |

**Read (all cells):**
- **Kinship structure matters most.** Every kNN-based variant with Complex Portal beats kNN, its null and the data graph on both masks. That is the first time the curated prior beats "data neighbours" on this benchmark.
- **But a plain low-rank model beats them all**, so they still lose to the strongest data-alone baseline.
- **E6 (curated residuals on top of low-rank) beats low-rank, its null and the data graph,** by about 0.003 RMSE.
- **E3's graph penalty helps,** but its null does almost as well (Δ ≈ 0.0005): this is mostly Laplacian regularisation, not curation.

**Low-n regimes: BUG FOUND, rerun pending.** `data_graph` required ≥ 30 co-observed cells per pair. That is impossible at 28 cells, so the data graph was empty and "data" collapsed to kNN. It is also restrictive at 60 cells. Fix: threshold = min(30, max(5, n/2)), unchanged at 215 cells. The low-n numbers from this run are **void** and will be replaced by the rerun (`dev_results_lown_fixed.json`).

**Selection (rule fixed in `EXPLORATION.md` before these runs):**
- **Eligible:** among methods whose curated version beats kNN, its null and the data graph on dev (215 cells, both masks), the lowest mean RMSE is **E3 low-rank + Laplacian with CORUM** (0.680).
- **Chosen:** it goes to confirmation (`PREREGISTRATION_v2.md`). Its dev margin over its own null is tiny, so we expect the confirmation to **fail** on the null comparison. We follow the rule regardless.

## 2026-10-21 · v2 CONFIRMATION (seeds 1000–1004, single run; `results/phase2/confirm/confirm_v2.json`)

**Primary (E3 low-rank + CORUM Laplacian): NOT CONFIRMED, 6/8 comparisons.** It beats kNN, low-rank-no-graph
and the data graph on both masks. It does **not** beat its own shuffled null:
- MCAR Δ = −0.0007 [−0.0031, +0.0016]
- MNAR Δ = +0.0001 [−0.0018, +0.0017]

The Laplacian helps through regularisation, not curation. This is the outcome predicted from dev.

**Secondary (reported, not claimed):**
- **E6, Complex Portal residuals on a low-rank base,** beats low-rank, its own null and the data graph on **both** masks, with all CIs excluding 0:

  | E6 (Complex Portal) vs. | MCAR ΔRMSE [95% CI] | MNAR ΔRMSE [95% CI] |
  |---|---|---|
  | Low-rank, no graph | −0.0033 [−0.0039, −0.0028] | −0.0030 [−0.0035, −0.0025] |
  | Own null | −0.0033 [−0.0039, −0.0028] | −0.0030 [−0.0035, −0.0025] |
  | Data-derived graph | −0.0026 [−0.0032, −0.0019] | −0.0016 [−0.0021, −0.0010] |

  The effect is small (~0.5% of RMSE). Its null gives the same RMSE as the low-rank base, so null "partners" add nothing (fitted β ≈ 0).
- **E1+E4** (gated, kNN base) beats kNN, its null and the data graph, but loses to low-rank.
- **E3 with STRING** has the lowest RMSE of all (MCAR 0.627, MNAR 0.706). It was not part of any pass test.

**What this suggests:** E6 is a candidate for an **independent** test. Next: `PREREGISTRATION_v3.md` on the
held-out Leduc 2022 dataset, chosen from Petrosius evidence only.

## 2026-10-21 · Low-data regimes, rerun with the data-graph fix (`dev_results_lown_fixed.json`)
In-complex RMSE, dev seeds 100–102, for E1+E4 (gated residual-on-kNN) with Complex Portal:

| Cells | Mask | data + CP | data + null | data alone | CP alone | kNN | low-rank |
|---|---|---|---|---|---|---|---|
| 28 | MCAR | **0.841** | 0.871 | 0.883 | 0.864 | 0.922 | 1.046 |
| 28 | MNAR | **0.861** | 0.874 | 0.904 | 0.876 | 0.932 | 1.065 |
| 60 | MCAR | **0.759** | 0.777 | 0.772 | 0.779 | 0.813 | 0.802 |
| 60 | MNAR | **0.849** | 0.856 | 0.869 | 0.854 | 0.886 | 0.889 |

**Read:**
- **Low-rank collapses with few cells,** so kNN-based stacking is the right base there.
- **Complex Portal adds information beyond the data graph:** data + CP beats data + null and data alone in every cell of the table.
- **CORUM is mixed:** at n = 28 MNAR, data + null (0.828) ≤ data + CORUM (0.835).

This motivates **v4** (low-data regime, independent Leduc dataset). v4 was preregistered *before* the v3 result was known.

## 2026-10-21 · v4 CONFIRMATION (Leduc 2022, low-data regime; `results/phase3/confirm_v4.json`)

**Primary: NOT CONFIRMED (2/6).** At n = 30 cells, data + Complex Portal:
- **vs. data + null:** beats it on MCAR (−1.05%, CI excludes 0). On MNAR the CI includes 0 (−0.35%).
- **vs. data alone:** beats it on MCAR (−0.50%). On MNAR the CI includes 0.
- **vs. kNN: loses on both masks** (+6.25% MCAR, +2.93% MNAR). On Leduc, residual correction on top of kNN *hurts*, whatever the graph.

**Secondary, worth recording (consistent across all 4 regime × mask cells, CIs excluding 0):** using **Complex Portal
partners alone** gives lower error than using **data-derived co-variation partners**. Both are in the same residual-on-kNN model, and both are still worse than kNN:

| Cells | Mask | Complex Portal partners vs. data-derived partners |
|---|---|---|
| 30 | MCAR | −5.3% RMSE |
| 30 | MNAR | −2.8% RMSE |
| 60 | MCAR | −4.7% RMSE |
| 60 | MNAR | −3.7% RMSE |

**Mixed:** CORUM gives −2.3% / −2.1% at n = 60, but at n = 30 it is not better (−0.2%, CI includes 0) and is worse on MNAR (+0.7%).

**Read:** on an independent dataset with few cells, curated complexes are a *better choice of partners* than
partners estimated from the data. However, the partner-residual model itself doesn't beat kNN here.

## 2026-10-21 · v3 CONFIRMATION (Leduc 2022, all 1543 cells; `results/phase3/confirm_v3.json`)

**Primary: passes the preregistered criterion (8/8), but the effect is practically negligible.**

| E6 + Complex Portal vs. (in-complex) | MCAR | MNAR |
|---|---|---|
| Low-rank base | −0.01% RMSE | −0.01% |
| Own null | −0.01% | −0.01% |
| Data-derived graph | −0.14% | −0.07% |
| kNN | −1.00% | −0.98% |

**Why "passes" overstates it:**
- **The correction is tiny:** the fitted correction weight is β ≈ 0.02 on every seed, so curated partners barely change the low-rank estimate.
- **The CIs are optimistic:** separation comes from pooling ~100k masked entries with correlated errors.
- **The kNN gain isn't the prior's:** the ~1% gain over kNN comes from the low-rank base.

**Process error, disclosed:** v3 had no *smallest effect size of interest*. With one, for example 1% of RMSE, the primary
would **fail** against the base and the null. We report it as: *passes the letter of v3, not its spirit*.

**Consistent signal across both datasets:** curated partners beat data-derived partners. This holds for
v3 vs. the data graph (small), v4 secondary (3–5% at 30–60 cells), Petrosius v2 secondary and Petrosius dev.
That is the defensible, narrower form of "better than data alone": *as the source of which proteins inform
a missing one*. It is not a claim that the curated prior beats the best data-only imputer.

## 2026-10-21 · v5 CONFIRMATION (SCoPE2, third dataset, SESOI 1%; `results/phase3/confirm_v5.json`)

**Primary: NOT CONFIRMED (0/4).** At n = 30, Complex Portal partners vs. data-derived partners:

| Mask | Δ RMSE | 95% CI | % of reference | Pass? |
|---|---|---|---|---|
| MCAR | −0.0043 | [−0.0104, +0.0016] | −0.78% | fails (CI includes 0 and below the 1% SESOI) |
| MNAR | −0.0063 | [−0.0121, −0.0004] | −0.95% | fails (below the 1% SESOI) |

Against its own null: −0.38% (MCAR) and −0.21% (MNAR), both with CIs including 0.

**Secondary:**
- **At n = 60,** Complex Portal beats data-derived partners by −1.45% (MCAR) and −1.48% (MNAR), with CIs excluding 0.
- **CORUM partners are worse than data-derived partners** on SCoPE2 (+3.7% at n = 30, MCAR).
- **No graph-based variant meaningfully beats kNN.**

**Where the evidence stands after three datasets:**
1. **Robust:** curated complex graphs carry real signal. They beat degree-preserving shuffles across tests and datasets (phase 1, v4, v5 n = 60), though the gain is often small.
2. **Directionally consistent but small, and not confirmed at a 1% SESOI:** Complex Portal partners beat
   data-derived partners. That's 3–5% on Leduc, ~1% on SCoPE2, and it isn't true for CORUM on SCoPE2.
3. **Not supported:** that a curated prior gives a *better imputer than data alone*. kNN or low-rank models match or beat
   every prior-based method by any practically meaningful margin, on all three datasets.

Further attempts on these three datasets would be forking paths. A new claim needs new data.

## 2026-10-21 · Amendment 4 dev: external-evidence priors (`results/phase4/dev/dev_ext.json`)

**Setup:** leave-one-dataset-out over the 5-dataset corpus. Petrosius uses 215/28 cells; Leduc and SCoPE2 use subsamples of 200/30. Dev seeds 100–102.
Each prior graph was tested in residual-on-kNN (E1) and residual-on-low-rank (E6). The priors were:
- Complex Portal (CP);
- CP gated by external co-variation;
- the external-only graph;
- external ∪ CP, vs. external ∪ CP-null;
- hu.MAP 3.0, alone and gated;
- an in-dataset data graph.

**Findings:**
- **External co-variation adds nothing.** The external-only graph never beats kNN or low-rank, gating CP by external
  evidence never beats plain CP, and ext ∪ CP is worse than CP alone.
  - Curated CP edges do have higher external correlation than random pairs (median r 0.057 vs 0.003), so the evidence is real but weak.
- **On Leduc and SCoPE2, every partner correction is ≈ 0.** RMSE equals the kNN or low-rank base to 3 decimals; β ≈ 0.
- **Only Petrosius at 28 cells shows a gain:** CP + E1 gives 0.862 / 0.878 vs. best data-only 0.884 / 0.905 (MCAR / MNAR). This is the
  dataset that motivated v4 and v5, which did not replicate.
- **hu.MAP 3.0 ≈ CP** wherever there is any signal.

**Selection rule (Amendment 4): no method qualifies.** No curated variant beats all baselines across corpus targets.
Following the rule, **v6 is not run as a win**, and khan2023 stays an untouched holdout for a future, different idea.

## 2026-10-22 · Amendment 5 dev: knowledge features in side-information MF (`results/phase5/dev/dev_imc.json`)
Overall RMSE, mean over dev seeds 100–102. "curated_all" = GO + HPA + UniProt keywords + PaxDb/length/mass + spectral
embeddings of Complex Portal and STRING.

| Target, cells | kNN | low-rank | IMC, no features | IMC, shuffled features | IMC, curated_all |
|---|---|---|---|---|---|
| Leduc 30, MCAR / MNAR | **0.546 / 0.633** | 0.835 / 0.868 | 0.725 / 0.772 | 0.682 / 0.742 | 0.664 / 0.733 |
| SCoPE2 30 | **0.551 / 0.646** | 0.912 / 0.943 | 0.765 / 0.831 | 0.710 / 0.773 | 0.692 / 0.764 |
| Petrosius 28 | 0.981 / 0.962 | 1.126 / 1.131 | 0.994 / 1.003 | 0.933 / 0.970 | **0.891 / 0.936** |
| at 200 / 215 cells | ≈ everything equal (±0.01) | | | | |

**Read:**
- **Knowledge features beat shuffled features in every 30-cell cell.** Combining sources is best, and no single source is enough.
- **But kNN dominates the low-rank family at 30 cells** on Leduc and SCoPE2.
- **Rule outcome:** IMC qualifies on 1 of 3 targets, so no v6 is run.

## 2026-10-22 · Amendment 6 dev: kNN + knowledge-IMC ensemble, w = 0.5 (`results/phase5/dev/dev_ens.json`)
Overall RMSE, dev seeds 100–102.

| Target, cells | kNN | ensemble (curated) | null ensemble (shuffled features) | kNN + IMC (no features) |
|---|---|---|---|---|
| Leduc 30, MCAR / MNAR | **0.546 / 0.633** | 0.566 / 0.648 | 0.572 / 0.653 | 0.582 / 0.657 |
| SCoPE2 30 | **0.551 / 0.646** | 0.565 / 0.653 | 0.576 / 0.665 | 0.593 / 0.681 |
| Petrosius 28 | 0.981 / 0.962 | 0.878 / 0.897 | **0.876** / 0.900 | 0.893 / 0.907 |
| Leduc 200 | 0.503 / 0.594 | 0.492 / 0.584 | 0.491 / 0.585 | **0.490 / 0.584** |
| SCoPE2 200 | 0.509 / 0.611 | 0.486 / 0.575 | 0.487 / 0.576 | **0.484 / 0.575** |

**Read:**
- **At 30 cells,** blending with knowledge-IMC makes kNN *worse* on Leduc and SCoPE2. The knowledge-feature ensemble does beat the shuffled-feature ensemble by 0.8–1.9%, so the knowledge signal is consistent but too weak to lift a strong kNN.
- **At 200 cells and on Petrosius,** ensembles beat kNN through model averaging, and the no-feature and shuffled ensembles do just as well.
- **Rule: does not qualify on any target.** No v6 is run, and khan2023 stays untouched.

## Overall status after Amendments 4–6 (2026-10-22 00:00 ET)
Nine distinct approaches have been tried, including:
- edge corrections;
- graph-regularised MF;
- residual stacking;
- curated × external gating;
- hu.MAP 3.0;
- a pooled 5-dataset corpus;
- knowledge feature representations from GO, HPA, UniProt, PaxDb, and spectral STRING/CP embeddings;
- kNN ensembles.

They were tested across 5 datasets, with 4 single-use confirmations. The pattern is the same throughout:
- **knowledge priors carry real, reproducible structure** (they beat shuffled versions in nearly every test);
- **they do not beat strong data-only imputers (kNN or low-rank) by a meaningful margin.**

## 2026-10-22 · Amendment 7 dev: protein-specific kNN, knowledge neighbourhoods, α = 0.5 (`results/phase5/dev/dev_pknn.json`)
Overall RMSE at 30 cells (Petrosius 28), mean over dev seeds; change of the curated version relative to each comparator:

| Target / mask | kNN | pkNN curated | vs kNN | vs null neighbourhoods | vs data neighbourhoods |
|---|---|---|---|---|---|
| Leduc MCAR | 0.546 | 0.552 | +1.1% | −1.4% | −1.3% |
| Leduc MNAR | 0.633 | 0.639 | +0.9% | −1.1% | −1.6% |
| Petrosius MCAR | 0.981 | **0.955** | −2.7% | −3.9% | −1.9% |
| Petrosius MNAR | 0.962 | **0.940** | −2.4% | −3.6% | −5.1% |
| SCoPE2 MCAR | 0.551 | 0.554 | +0.6% | −1.1% | −1.7% |
| SCoPE2 MNAR | 0.646 | 0.648 | +0.2% | −1.3% | +0.3% |

**Read:**
- **Knowledge-defined neighbourhoods beat shuffled ones in every case, and beat data-defined ones in 5/6.**
- **They beat plain kNN only on Petrosius.** On Leduc and SCoPE2 they are within ~1% of kNN.
- **Rule: qualifies on 1 of 3, so it does not qualify.** It is the closest approach so far.

## 2026-10-22 · Amendment 8 dev: pkNN α sweep (`results/phase5/dev/dev_pknn_sweep.json`)
**α = 0.1:**
- **Leduc:** kNN 0.5458 / 0.6331 vs. pkNN curated 0.5459 / 0.6304 (MCAR / MNAR). On MCAR a shuffled-neighbourhood null scores 0.5454.
- **SCoPE2:** curated pkNN is 0.5527 vs. kNN 0.5513 on MCAR.

**α = 0.25:** wins on Petrosius only.

**Rule: no α qualifies** (≥ 2 of 3 targets, both masks), so the sweep stops as pre-committed.
As the local weight shrinks, pkNN converges to kNN and the knowledge neighbourhood stops mattering.

## 2026-10-22 · Amendment 9 dev: exhaustive representation sweep (`results/phase6/dev/dev_reps.json`, `dev_reps_selection.json`)
**Scope:**
- 22 protein representations;
- three estimators: RPR, RPR-stack and pkNN;
- each with its own label-permutation nulls;
- baselines: kNN, low-rank, the kNN + low-rank stack and data-only neighbours;
- three dev targets, MCAR and MNAR, seeds 100–102, at 30 and 200 cells.

6,138 rows in total. The first run had a NaN bug, which is fixed and logged as deviation 1; the buggy run is kept in `run1_buggy/`.

**Selection rule:** 0 of 66 (representation, estimator) pairs qualify. **No Amendment-9 confirmation**, so khan2023 passes to
Amendment 10.

**Baselines at 30 cells** (RMSE relative to kNN): low-rank 1.389; kNN + low-rank stack **0.983**; data-only RPR 1.039;
data-only pkNN 1.016. At 200 cells: low-rank 0.979; stack **0.953**.

**RPR alone at 30 cells**, ranked: mean RMSE relative to kNN, then its own null.

| Representation | rel. kNN | own null |
|---|---|---|
| **Complex Portal** | **1.005** | 1.020 |
| **PDB co-structure** | 1.009 | 1.026 |
| **CORUM** | 1.012 | 1.028 |
| Monarch phenotypes | 1.022 | 1.024 |
| NCBI co-citation | 1.028 | 1.058 |
| Reactome | 1.028 | 1.036 |
| hu.MAP 3.0 | 1.031 | 1.053 |
| STRING | 1.032 | 1.061 |
| SCP corpus loadings | 1.032 | 1.047 |
| *data-only partners* | 1.039 | — |
| CCLE (reference coefficients) | 1.041 | 1.061 |
| BioGRID, ProtT5, GO, Pathway Commons, UniProt KW, HPA, IntAct, PaxDb/physical | 1.043–1.058 | ~1.065 |
| ProteomeHD, Pharos | 1.062 | 1.062–1.065 |

**Read:**
- **The thesis's ranking claim holds.** Among 22 representations, *curated physical complex membership* is the best
  protein prior. Complex Portal ranks first, PDB co-structure (physical, structural) second and CORUM third.
  - Each beats its own null by 1.5–1.7%, and data-only partners by 2.6–3.3%.
  - Broader functional knowledge (pathways, PPI networks, sequence embeddings, drug-target annotations) helps less.
  - Bulk co-regulation panels (CCLE, ProteomeHD) help little or not at all at the single-cell level.
- **Usefulness does not hold.** No representation beats kNN. When stacked with kNN, the inner-mask weight mostly
  returns to kNN (RPR-stack ≈ 0.999–1.000 vs. kNN). The best pkNN (concatenated features, 0.9965) is 0.35% better than
  kNN, below the 1% SESOI and worse than the generic kNN + low-rank stack.
- The knowledge is real: RPR beats its null in 19/22 representations and pkNN in 18/22. But at the single-cell
  level the protein-side signal it adds is smaller than what kNN already captures from cell similarity.

## 2026-10-22 · Amendment 10 dev: transfer from external references (`results/phase7/dev/dev_transfer.json`)
**Mean RMSE relative to kNN at 30 cells** (Leduc MCAR / MNAR · Petrosius · SCoPE2):

| Method | Leduc | Petrosius | SCoPE2 |
|---|---|---|---|
| kNN + low-rank stack (baseline) | 1.000 / 1.000 | **0.932 / 0.962** | 1.000 / 1.000 |
| `refmf_all` alone (corpus + CCLE + ProteomeHD loadings) | 1.023 / 1.021 | 0.980 / 1.028 | 1.002 / 0.990 |
| its permuted-loadings null | 1.05 / 1.05 | 1.07–1.08 / 1.15–1.16 | 1.03 / 1.01–1.02 |
| **kNN + `refmf_all` stack** | 0.995 / 0.997 | 0.954 / 0.993 | **0.985 / 0.986** |
| kNN + `refmf_corpus` stack | 0.998 / 1.001 | 0.978 / 1.012 | 0.988 / 0.982 |
| `refmf_all_kg` (KG cold start) vs its KG null | identical to within 0.1% | | |
| reference-donor kNN | 1.07–1.13 | | |

**Selection rule:**
- **0 of 10 qualify.** The best is kNN + `refmf_all` stack: it wins on SCoPE2 under both masks, but on Petrosius the
  generic kNN + low-rank stack is better.
- At 200 cells, transfer does not help: target-learned low-rank is better.

**Read:**
- **Transfer is the first mechanism that adds signal kNN does not already have.** Reference loadings beat permuted ones
  by 3–12% at 30 cells, adding bulk panels (`all`) beats the single-cell corpus alone, and the stack wins
  1.4–1.5% on SCoPE2.
- **The knowledge-graph cold start did nothing measurable** (≈ its null): the 4–26% of cold proteins are not where the
  error is.
- **The kNN + low-rank and kNN + transfer stacks help *different* targets** (Petrosius vs. SCoPE2). The reference is
  small (≈ 4,000–5,600 cells from 4 datasets, plus 669 bulk samples).

## 2026-10-22 · Amendment 11 dev: bigger reference (17 datasets, ~8,600–9,900 reference cells) + 3-way stack (`results/phase8/dev/`)
**Mean RMSE relative to kNN at 30 cells** (MCAR / MNAR):

| Method | Leduc | Petrosius | SCoPE2 |
|---|---|---|---|
| kNN + low-rank stack | 1.000 / 1.000 | 0.932 / 0.962 | 1.000 / 1.000 |
| `stack2_all` (Amendment 10 reference) | 0.995 / 0.997 | 0.954 / 0.993 | 0.985 / 0.986 |
| `refmf_big` alone | 1.020 / 1.019 | 0.898 / 0.972 | 1.002 / 0.983 |
| `stack2_big` | 0.993 / 0.997 | 0.894 / 0.959 | 0.985 / 0.983 |
| **`stack3_big`** (kNN + low-rank + `refmf_big`) | 0.993 / 0.997 | **0.860 / 0.933** | **0.985 / 0.983** |
| `stack3_big` with permuted loadings (null) | 1.003 / 1.004 | 0.924 / 0.974 | 0.999 / 0.998 |

**Selection rule:** **0/3 qualify.** `stack3_big` wins on Petrosius only (1/3 targets).
- On SCoPE2 it beats kNN, low-rank, the kNN + low-rank stack and its null by ≥ 1.4%. It fails only against the preregistered
  extra baseline `stack2_all` (the smaller Amendment-10 reference): the bigger reference adds nothing on SCoPE2.
- On Leduc nothing beats kNN by 1%.

At 200 cells, `stack3_big` ≈ the kNN + low-rank stack (the projection weight goes to 0).

**Read:**
- The bigger reference makes transfer much stronger on Petrosius: `refmf_big` alone is −10% vs kNN, where the
  Amendment-10 reference gave −2%.
- *Post-hoc observation (not a claim):* had `stack2_all` not been a preregistered bar, `stack3_big` would have passed on
  2/3 targets. We keep the rule as written.

## 2026-10-22 · v6 confirmation on the khan2023 holdout (Amendment 12, post-gate): **H12 NOT CONFIRMED (3/8 primary)** (`results/confirm_v6/confirm_v6.json`)
**Setup:**
- 30-cell subsamples, seeds 6000–6009.
- Reference: 17 datasets plus bulk panels, 10,129 samples × 4096 proteins. 84% of khan2023 proteins have reference
  loadings.

| Comparator | MCAR ΔRMSE (95% CI) | MNAR ΔRMSE (95% CI) |
|---|---|---|
| kNN (0.446 / 0.525) | **−1.40%** [CI < 0] PASS | −0.98% [CI < 0] fail (just under SESOI) |
| low-rank | −28.3% PASS | −22.9% PASS |
| kNN + low-rank stack | −0.55% [CI < 0] fail | −0.20% [CI < 0] fail |
| own permuted-loading null | −0.61% [CI < 0] fail | −0.20% [CI < 0] fail |
| *secondary:* `stack2_all` | −1.04% | −0.85% |

**Read:**
- On a completely untouched dataset (MCF-10A EMT, a different lab-year and cell type), the 3-way transfer stack is the best
  method, and **every difference is statistically resolved (CIs exclude 0)**.
- **But the margin over the strongest generic baseline (kNN + low-rank) and over the permuted reference is 0.2–0.6%,
  below the preregistered 1% SESOI.** The selected stack weights put 0.25–0.5 on the reference projection in 9/10
  subsamples, so the projection is used, but its *knowledge-specific* contribution over a permuted reference is small.
- **Conclusion: H12 is not confirmed.** This was the single use of khan2023.

### Where this leaves the thesis (Oct 22, after Amendments 4–12)
- **Curated physical complex membership is the best protein-level prior we tested** (ranked #1 of 22), and
  knowledge priors consistently beat their own nulls.
- But neither curated priors nor any of 21 other knowledge sources beats plain kNN by a meaningful margin at single-cell
  resolution.
- The only mechanism that adds new signal is **transfer of protein covariance from large external references** (≈ 10k
  cells). On the holdout it is real but small (≈ 0.2–0.6% over strong baselines; 1.0–1.4% over kNN).
- **Honest headline:** a rigorous, preregistered negative result with one small, reproducible positive (transfer).

## 2026-10-22 · Amendment 13 dev: biology-integrated model, in-complex endpoint (`results/phase9/dev/`)
**In-complex RMSE relative to kNN at 30 cells** (MCAR / MNAR):

| Method | Leduc | Petrosius | SCoPE2 |
|---|---|---|---|
| kNN + low-rank stack | 1.000 / 1.000 | 0.934 / 0.956 | 1.000 / 1.000 |
| `stack3_plain` (transfer, no biological prior) | **0.988 / 0.993** | 0.858 / 0.917 | **0.990 / 0.982** |
| `stack3_kg` (transfer with κ = 0.25 complex prior on the loadings) | 0.991 / 0.996 | **0.856 / 0.913** | 0.992 / 0.985 |
| **BIM** (+ complex RPR in a 4-way stack) | 0.991 / 0.996 | 0.856 / 0.913 | 0.992 / 0.985 |
| `bim_null` (prior labels permuted) | 0.993 / 0.997 | 0.880 / 0.934 | 0.991 / 0.988 |

**Gate: fails** (0/3 targets pass both H13a and H13b under both masks). brunner2022 stays untouched.
- **H13a (useful vs kNN and the kNN + low-rank stack by ≥ 1%)** holds on Petrosius under both masks and on SCoPE2 MNAR only.
  As in Amendment 11, this comes from *transfer*.
- **H13b (biology-specific)** fails. Adding the complex prior to the transferred loadings is ~0.2–0.3% *worse* than plain
  transfer on Leduc and SCoPE2, and only 0.2–0.4% better on Petrosius.
  - BIM beats its permuted-prior null clearly only on Petrosius (−2.2–2.7%). There the null is a *wrong* prior that
    hurts, not a correct prior that helps.
- **The inner-mask stack gave the complex RPR component weight > 0 in only 1 of 36 subsamples.**
- At 200 cells, BIM ≈ the kNN + low-rank stack.

**Read:** once protein covariance is learned from ~10k external cells, curated complex membership adds no further
information, because the reference already contains it.
- Complex partners are among the most strongly co-varying proteins in the data: Complex Portal ranked #1 in
  Amendment 9 *because* it matches real covariance.
- So the prior is **correct but redundant** given enough data. Without that data (30 cells alone), it is correct but too
  weak to beat kNN.

## 2026-10-22 · v7 confirmation on the brunner2022 holdout (Amendment 14): **H14 CONFIRMED, 6/6 primary** (`results/confirm_v7/confirm_v7.json`)
**Setup:**
- HeLa, timsTOF DIA, 434 cells; 30-cell subsamples (6 per cell-cycle group), seeds 7000–7009.
- Reference: 16 public single-cell datasets + CCLE + ProteomeHD (9,695 samples). 99.9% of brunner proteins have
  reference loadings, and 62% of proteins are complex members.

| In-complex RMSE, `stack3_plain` vs | MCAR (Δ%, 95% CI of ΔRMSE) | MNAR |
|---|---|---|
| kNN (0.457 / 0.605) | **−5.99%** [−0.036, −0.020] ✅ | **−8.88%** [−0.057, −0.050] ✅ |
| kNN + low-rank stack | **−2.44%** [−0.015, −0.007] ✅ | **−5.42%** [−0.034, −0.029] ✅ |
| permuted-reference null | **−2.51%** [−0.015, −0.008] ✅ | **−4.30%** [−0.027, −0.022] ✅ |

**Secondary:**
- Overall RMSE passes all six too: −4.9 / −7.9% vs kNN, −1.9 / −4.9% vs the stack, −1.9 / −3.7% vs the null.
- **Curated-complex-specific:** BIM vs `stack3_plain` is −0.17% / −0.06%, with CIs including 0. BIM vs `bim_null` is
  −0.54% / −0.34%, CIs including 0. **Not significant.**
- The stack weights gave the reference projection 0.5–0.75 in all 20 subsamples.

**What is and isn't shown:**
- ✅ **Empirical biological prior knowledge transfers.** Protein co-regulation learned from public single-cell and bulk
  proteomics, combined with kNN and low-rank, imputes an unseen dataset of a different technology substantially better
  than the best data-only method and than a reference with scrambled protein identities. This holds on complex members
  and overall.
- ❌ **Curated complex membership adds nothing measurable** on top of that empirical prior.
- ⚠️ **Caveats:**
  - H14 was refined after dev results (Amendment 14).
  - On the earlier holdout khan2023 (MCF-10A, Amendment 12, overall endpoint), the closely related `stack3_big` gained
    only 0.2–0.6% over the kNN + low-rank stack and did not pass. **So the effect replicates on 1 of 2 holdouts, and its
    size likely depends on dataset and reference similarity.**
  - Only 7 HeLa single cells are in the reference. We did not check the cell-line composition of ProteomeHD.

---

# Phase B (post-hackathon; `main` frozen at `076a946` as the submission)

## 2026-10-23 · Ceiling diagnostic: how much could transfer gain with an *ideal* reference? (`results/phaseB/dev/dev_ceiling.json`)
**Setup:**
- Same projection method as `stack3_plain`, but the reference is the target dataset's **own other cells** (same cell
  types, lab and instrument; disjoint from the 30-cell subsample).
- It is a diagnostic, not a candidate.
- Phase-B reference: wu2026, brunner2022 and derks2022 excluded as holdouts.

RMSE relative to the best data-only method (kNN + low-rank stack), MCAR / MNAR:

| Target | `stack3_plain` (public reference) | ideal in-domain reference |
|---|---|---|
| Leduc | 0.993 / 0.996 | **0.903 / 0.907** |
| Petrosius | 0.939 / 0.983 | **0.790 / 0.829** |
| SCoPE2 | 0.984 / 0.984 | **0.860 / 0.841** |
| **mean** | **0.980** | **0.855** |

**Read:**
- With an ideal reference, transfer would give −10% to −21% (mean −14.5%) vs the best data-only method. Today's
  public reference captures about 1/7 of that.
- 20–30% is beyond this linear ceiling on average, though not on Petrosius.
- The gap points at reference *relevance*: matching cell types and acquisition.
- Next: reference selection (choose the external cells most similar to the target, as genotype imputation chooses the
  closest haplotypes), and nonlinear models.

## 2026-10-23 · Phase B dev: reference selection and masked autoencoder (3 original dev targets, `results/phaseB/dev/`)
RMSE relative to the kNN + low-rank stack, mean over 3 targets × 2 masks:

| Method | `dev_select` | `dev_mae30` |
|---|---|---|
| `stack3_plain` (full public reference) | 0.974 | 0.970 |
| top-2000 most similar reference cells | 0.973 | — |
| 2000 random reference cells (control) | 0.981 | — |
| top-3 most similar reference datasets | 0.983 | — |
| MAE alone (30 epochs) | — | 1.025 |
| 4-way stack + MAE / + scrambled-reference MAE | — | 0.970 / 0.970 |

**Read:**
- Similarity-based selection beats a random subset of equal size by 0.8%, but not the full reference.
- The MAE is worse than kNN alone, and the stack gives it no weight.
- Neither closes the gap to the in-domain ceiling (0.855).
- Run-to-run ALS variability of `stack3_plain` is about 0.5% (0.974 vs 0.970 vs 0.980 for the same configuration
  with different row orders).

## 2026-10-23 · Phase B dev: which reference carries the signal? (`dev_bigbulk`, `dev_sc_only`, `dev_bulk_only`)
`stack3_plain` RMSE relative to the kNN + low-rank stack (3 dev targets × 2 masks):

| Reference | Score |
|---|---|
| single-cell corpus only (13 datasets) | **0.970** |
| single-cell + 6 bulk atlases (CCLE, ProteomeHD, ProCan, GTEx, CPTAC, ATLANTiC; ~4,000 samples) | 0.971 |
| bulk atlases only | 0.982 |

**Read:**
- Adding ~2,800 bulk samples (ProCan, GTEx, CPTAC, ATLANTiC) does not help.
- The transferable structure comes from *single-cell* references: same measurement physics, carrier and missingness
  patterns.
- Next: more single-cell datasets (agent in progress), and technology-matched references.

## 2026-10-23 · Phase B dev: more single-cell data, then domain-adapted projection (`dev_sc25`, `dev_adapt*`)
**More data.** A single-cell-only reference with 25 datasets (8 new PRIDE datasets, ~5,400 more cells) gives
`stack3_plain` = 0.973, the same as with 13 datasets (0.970). More generic reference data has plateaued.

**Domain adaptation.** `transfer.adapt_project` is MAP low-rank on the 30 target cells, with the reference loadings as
a Gaussian prior (strength λ_v) plus `extra` free factors.

Scores are relative to the kNN + low-rank stack, mean of 3 targets × 2 masks; single-cell reference only.

| Configuration | With reference | Scrambled reference | No reference (zero prior, same free factors) |
|---|---|---|---|
| plain projection (Phase A) | 0.973 | | |
| λ_v = 30, 3 extra | 0.954 | | |
| λ_v = 100, 0 extra | 0.970 | | |
| **λ_v = 100, 10 extra** | **0.953** | 0.975 | 0.975 |
| λ_v = 100, 20 extra | 0.953 | 0.976 | 0.976 |
| λ_v = 300 / 1000, 3 extra | 0.971 / 0.973 | | |

**Read:**
- Letting the target bend the reference loadings, and add a few factors of its own, roughly doubles the transfer gain:
  −4.75% vs −2.7%.
- The reference accounts for about 2.3 points of that (vs. scrambled or no reference).
- λ_v and `extra` were tuned on these 3 dev targets. Next: check that this generalises across the full LODO dev pool
  before anything goes near a holdout.

## 2026-10-23 · Phase B dev: adapted transfer generalises across the LODO dev pool (`results/phaseB/dev/dev_lodo16.json`)
**Setup:**
- 15 dev targets: the 3 originals plus 12 corpus datasets with ≥ 100 cells.
- woo2022 was excluded: it turned out to be **mouse**, with 0 accessions shared with the human datasets. Its reference
  rows were always empty after alignment, so it never contributed. It is now in `transfer.NON_HUMAN`.
- Single-cell reference, leave-target-out. Settings: λ_v = 100, 10 extra factors.

Mean relative to the kNN + low-rank stack, 15 targets × 2 masks:

| Method | All entries | In-complex |
|---|---|---|
| Phase-A `stack3_plain` | 0.962 | 0.953 |
| **adapted** | **0.936** | **0.927** |
| adapted, scrambled reference | 0.965 | 0.963 |
| adapted, no reference | 0.975 | 0.973 |

**Per target:**
- The adapted method beats the best data-only method by ≥ 1% on 14/15 targets.
- It beats both controls by ≥ 1% on 9/15.
- Best: weiss2026 liver (0.853 / 0.804), amrute2026 aorta (0.90 / 0.89), sadiku2025 neutrophils, Petrosius, guise2024
  neurons.
- Weakest: gregoire2023 (≈ 0.99–1.00).

**Gate:** this passes the Amendment-15 dev gate (≥ 1% mean gain over the current best, and beats its controls).
`PREREGISTRATION_v8.md` commits the single wu2026 confirmation.

## 2026-10-23 · v8 confirmation on the wu2026 holdout (Amendment 15, H15): **CONFIRMED, 16/16** (`results/confirm_v8/confirm_v8.json`)
**Setup:**
- Prenatal human brain, 2268 cells, DIA, 90% missing.
- 30-cell random subsamples, seeds 8000–8009.
- Single-cell reference: 21 datasets, no holdouts.
- Adapted projection (λ_v = 100, 10 extra factors) stacked with kNN and low-rank.

| `adapt` vs | MCAR, all | MCAR, in-complex | MNAR, all | MNAR, in-complex |
|---|---|---|---|---|
| kNN | −14.2% | −15.4% | −15.4% | −16.1% |
| **kNN + low-rank (best data-only)** | **−5.0%** | **−4.9%** | **−6.8%** | **−6.8%** |
| scrambled reference | −5.2% | −5.3% | −6.1% | −6.1% |
| no reference | −5.1% | −5.1% | −7.2% | −7.2% |

All 16 paired-bootstrap 95% CIs are below 0, and every difference is beyond the 1% SESOI.

**Read:**
- On a tissue and technology not used in development, borrowing protein co-regulation from 21 public single-cell
  datasets, and letting the target adapt it, cuts imputation error by 5–7% vs. the best data-only method and by
  14–16% vs. kNN.
- The scrambled-reference and no-reference controls show the gain comes from the *biological content* of the
  reference.
- This is about 2× the Phase-A holdout effect on Brunner (2.4–5.4%), on a harder (90%-missing) dataset.
- derks2022 remains reserved.

## 2026-10-23 · Phase C round 1 (dev, fast pool): four new information sources / encodings, all dead (`results/phaseC/*_fast.json`)
**Setup:** `benchmarks/phasec.py --pool fast` (6 dev targets × 2 masks × 2 seeds); every number is rmse_all mean
relative to the H15 champion (< 1 is better), all / in-complex. Bars for the full pool: cand_stack ≤ 0.99 and ≥ 0.005
below the null. Charter and per-idea docs: `docs/phaseC/`.

| id (queen) | idea | cand_stack | null_stack | cand_solo |
|---|---|---|---|---|
| `ttprior` (bio) | metadata-kernel (tissue + technology) dataset weights in `fit_loadings`, τ = 1 | 0.9951 / 0.9953 | 0.9991 / 1.0005 | 0.9944 |
| `ebshrink` (rep) | per-protein λ_v,j from cross-dataset loading variance | 0.9986 / 0.9983 | 1.0000 / 0.9999 | 0.9991 |
| `scrnaprior` (mod) | HPA single-cell-type RNA SVD (k = 8) as an extra factor block | 0.9973 / 0.9974 | 0.9987 / 0.9988 | 0.9989 |
| `dca` (abd) | shrinkage-inverse precision (direct couplings), Gaussian conditional mean | 1.0004 / 1.0006 | 1.0000 / 1.0000 | 1.032 |

**Read:**
- None reaches the bar; none went to the full pool. The champion stands.
- Re-shaping the prior on V is saturated: the best gain is 0.5% (ttprior) and only 0.4% beyond its null.
- ttprior's kernels separate the axes: technology-only 0.9970 (null 0.9980) vs tissue-only 0.9992 (null 0.9997). Its
  gain sits on Astral-DIA targets (sadiku2025 MCAR 0.9797); TMT-carrier targets stay at 0.997–1.001.
- dca's conditional mean is real signal on its own (cand_solo 1.032 vs stack_lowrank 1.0773) but gets stack weight
  0.0104: the adapted projection covers it. Its couplings put 0.0799 of top-30 partners in a shared CORUM complex vs
  0.1063 for raw correlation and 0.0177 for its null (`dca_couplings_diag.json`).
- Round 2 (king's directives in `docs/phaseC/KING_LOG.md`) moves to cell-level marker weighting (`ctprior`), a
  nonparanormal encoding (`copula`), peptide-level precision weights (`pepev`) and Netflix-style shrunken target
  baselines (`nfbase`).

## 2026-10-23 · Phase C round 2 (dev, fast pool): four more ideas dead; the target's standardisation is the lead (`results/phaseC/*_fast.json`)
**Setup:** as in round 1 (`benchmarks/phasec.py --pool fast`, 6 dev targets × 2 masks × 2 seeds; rmse_all mean
relative to the H15 champion, all / in-complex; full-pool bars: cand_stack ≤ 0.99 and ≥ 0.005 below the null).

| id (queen) | idea | cand_stack | null_stack | cand_solo |
|---|---|---|---|---|
| `ctprior` (bio) | HPA-marker cell-type soft cell weights in `fit_loadings` | 0.9990 / 0.9986 | 0.9990 / 0.9987 | 0.9988 |
| `copula` (rep) | nonparanormal encoding of reference and target | 1.0005 / 1.0011 | 1.0004 / 1.0004 | 1.0193 |
| `pepev` (mod) | precursor counts as U-step precision weights | 0.9996 / 0.9987 | 0.9998 / 0.9996 | 0.9988 |
| `pepev_v2` (mod) | protein-level precursor-count weights | 0.9969 / 0.9967 | 1.0003 / 1.0002 | 0.9969 |
| `nfbase` (abd), ψ = 1 | shrink the target's mean/sd toward reference tech-group baselines | 0.9868 / 0.9859 | 0.9867 / 0.9863 | 0.9864 |
| `nfbase_psi4` (abd) | the same, ψ = 4 | 0.9663 / 0.9647 | 0.9691 / 0.9688 | 0.9663 |

**Read:**
- No idea clears both bars; nothing went to the full pool; the champion stands.
- Re-weighting or re-encoding the reference is closed: ctprior's margin over its null is 0.0000; copula's
  reference-side-only variant is 0.9980 vs null 0.9991; the target-side quantile decode hurts (tgtonly cand_raw 1.0284).
- pepev's weights help MCAR and hurt MNAR (weiss2026_liver 0.9759 / 1.0121): in DIA, precursor counts track the value
  within protein (`pepev_diag.json`).
- **The lead is data-only.** nfbase's protein-permuted null at ψ = 4 (James–Stein shrinkage of the target's own
  30-cell per-protein mean and log-sd toward a calibrated constant) is 0.9691, ≤ 0.99 on 12/12 splits. Post hoc:
  location only 0.9920, scale only 0.9832, joint 0.9663. Protein-specific reference baselines add only on DIA MCAR
  (weiss2026_liver 0.9490 vs null 0.9645).
- Round 3 (directives in `docs/phaseC/KING_LOG.md`): `tgtshrink` (data-only champion revision, ψ on the inner mask,
  full pool if it holds), `selmnar` (probit selection model for the missing entries), `seqprior` (sequence
  detectability features as shrinkage targets) and `varprior` (half-life / complex-buffering prior on log-sd). The
  last two must beat shrink-to-constant at the same ψ.

## 2026-10-23 · Phase C round 3 (dev): shrinkage is the lever, not its target; selected-ψ constant shrinkage reaches 0.9650 on the full pool (`results/phaseC/*.json`)
**Setup:** as in rounds 1–2 (`benchmarks/phasec.py`, fast pool = 6 dev targets × 2 masks × 2 seeds; full pool = 15
LODO dev targets × 2 masks × 3 seeds; rmse_all mean relative to the H15 champion, all / in-complex).

| id (queen) | idea | cand_stack | null_stack | cand_solo |
|---|---|---|---|---|
| `tgtshrink` (abd), fast | shrink target mean/log-sd toward a (mean, obs-fraction) trend, ψ selected on a 50/50 inner mask | 0.9597 / 0.9586 | 0.9575 / 0.9563 | 0.9640 |
| `tgtshrink` (abd), **full** | the same | 0.9683 / 0.9684 | **0.9650 / 0.9653** | 0.9699 |
| `tgtshrink_psi4` / `_const_psi4` | fixed ψ = 4: trend / constant | 0.9730 / 0.9691 | 0.9683 / — | 0.9756 / 0.9692 |
| `selmnar` (rep) | Heckman probit selection model in an EM U-step | 1.0000 / 1.0000 | 0.9985 / 0.9986 | 1.0192 |
| `seqprior` (mod) | sequence/HPA detectability features as shrinkage targets | 0.9686 / 0.9689 | 0.9687 / 0.9694 | 0.9687 |
| `varprior` (bio) | half-life (Mathieson 2018) + CORUM as a log-sd shrinkage target | 0.9679 / 0.9686 | 0.9718 / 0.9728 | 0.9678 |

**Read:**
- No idea clears its bars; the champion stands. tgtshrink's trend is worse than its permuted null (fast −0.0022, full
  −0.0033). seqprior ties its shrink-to-constant baseline (0.9686); varprior beats its constant baseline by 0.0007 and
  its null by 0.0039 (v2 0.9682 vs null 0.9686).
- **The useful result is data-only:** tgtshrink's null (shrink the target's per-protein mean and log-sd toward a
  constant, ψ chosen per target on an inner mask) is 0.9575 fast and 0.9650 on the full pool (MCAR 0.9736, MNAR 0.9564;
  25/30 target-masks ≤ 0.99; worst woo2025_PC3 MCAR 1.0000), vs 0.9691 for fixed ψ = 4 on fast. ψ = 16, the top of the
  grid, was chosen in 30/48 null calls.
- selmnar's missingness model is strongly mask-dependent: cand_raw MNAR 0.8774 (permuted-indicator null 0.9615, no
  selection 0.9967) but MCAR 1.2721; the stack gives it weight 0.0000. The MNAR eval mask is value-logistic, so part of
  that gain may be the model re-learning the mask generator.
- Round 4 (directives in `docs/phaseC/KING_LOG.md`): `cshrink` (explicit selected-ψ constant shrinkage, extended grid,
  random-ψ null, full pool), `mnarjoint` (joint low-rank of values and missingness indicators), `subthr` (sub-threshold
  DIA precursor evidence at naturally missing entries) and `techhyper` (technology/tissue-aware hyperprior on ψ). All
  are credited only against the selected-ψ constant base.

## 2026-10-23 · Phase C round 4 (dev): the standardisation lever is used up; `cshrink` reaches 0.9609 on the full pool (`results/phaseC/*.json`)
**Setup:** as in rounds 1–3 (`benchmarks/phasec.py`; fast = 6 dev targets × 2 masks × 2 seeds; full = 15 LODO dev
targets × 2 masks × 3 seeds; rmse_all mean relative to the H15 champion, all / in-complex). From round 4 every idea is
credited against the selected-ψ shrink-to-constant base.

| id (queen) | idea | cand_stack | null_stack | base / control | cand_solo |
|---|---|---|---|---|---|
| `cshrink` (abd), fast | explicit selected-ψ shrink-to-constant of target mean/log-sd, ψ ∈ {1…64, ∞} | 0.9554 / 0.9533 | 0.9637 / 0.9641 (random ψ) | fixed ψ 16: 0.9576 | 0.9548 |
| `cshrink` (abd), **full** | the same | **0.9609 / 0.9610** | 0.9716 / 0.9728 | fixed ψ 16 (`cshrink_gsd`): 0.9624 | 0.9606 |
| `mnarjoint` (rep) | [X ∣ w·M] shared cell factors, w on the inner mask | 0.9576 / 0.9564 | 0.9576 / 0.9564 | w = 0: 0.9575 | 0.9575 |
| `subthr` (mod), 2 targets | sub-threshold PSM / PEP evidence at naturally missing entries | 0.9647 / 0.9628 | 0.9609 / 0.9598 | κ = 0: 0.9640 | 0.9631 |
| `techhyper` (bio) | metadata-kernel hyperprior on ψ | 0.9568 / 0.9544 | 0.9569 / 0.9543 | flat prior: 0.9567 | 0.9559 |

**Read:**
- cshrink beats its random-ψ null (full MCAR 0.9714, MNAR 0.9504 vs null 0.9793 / 0.9640; 26/30 target-masks ≤ 0.99)
  but misses the bar for a v9 draft (MCAR > 0.97; woo2025_PC3 MCAR 1.0006). Selection over fixed ψ = 16 adds only
  0.0022 fast / 0.0015 full, ψ = ∞ is 0.9559 on fast, and tau²(log-sd) sits at its floor: the gain is from dropping the
  champion's per-protein unit-variance scaling, not from tuning ψ. techhyper confirms the ψ curve is flat above ≈ 16.
- mnarjoint: the inner mask chose w = 0 in 37/48 calls; forced w = 1 is 0.9654 (permuted-pattern null 0.9773), with
  cand_raw MCAR 0.9896 / MNAR 0.9535 vs base 0.9697 / 0.9412.
- subthr: MCAR cand_raw 0.9771 vs base 0.9743; the sadiku MNAR gain (raw 0.9021 vs 0.9190) is reproduced better by
  the null (stacked 0.9293 vs 0.9409), i.e. it is "missing = low", not the PEP content.
- Weak spot left: TMT/carrier MCAR on the full pool (woo2025_PC3 1.0006, williams2020_tmt 0.9989, schoof2021 0.9972).
- Round 5 (directives in `docs/phaseC/KING_LOG.md`) moves to covariance in the target's own domain: `mllr` (shared
  r × r transform of the reference loadings), `lwkrige` (Ledoit–Wolf / kriging blend of target and reference
  covariance, with a technology-matched arm), `linerna` (within-cell-line scRNA co-variation as a new modality) and
  `plmemb` (protein language-model embeddings as a prior on weakly supported loadings).

## 2026-10-23 · Phase C round 5 (dev): only the target's own covariance carries signal; all four ideas dead (`results/phaseC/*.json`)
**Setup:** as in round 4 (fast = 6 dev targets × 2 masks × 2 seeds; rmse_all mean relative to the H15 champion, all /
in-complex). Every idea is credited against the `cshrink` base (fast 0.9554 / 0.9533; each `<id>_base` run reproduces
it exactly).

| id (queen) | idea | cand_stack | null_stack | base / control | cand_solo |
|---|---|---|---|---|---|
| `lwkrige` v2 (rep) | Ledoit–Wolf blend of the target's own covariance with the adapted Σ_ref, δ on the inner mask | 0.9522 / 0.9494 | 0.9553 / 0.9531 | base 0.9554; tech-matched Σ_ref 0.9554 | 0.9516 |
| `plmemb` (bio) | ProtT5 embedding kernel prior on weak / cold Vref rows | 0.9529 / 0.9520 | 0.9539 / 0.9517 | Gaussian embeddings 0.9535; base 0.9554 | 0.9524 |
| `mllr` (abd) | shared r × r MLLR transform of Vref | 0.9552 / 0.9531 | 0.9538 / 0.9528 | base 0.9554; diagonal 0.9550 | 0.9546 |
| `linerna` v2 (mod), 4 targets | within-cell-line scRNA programmes on the cshrink residual | 0.9668 / 0.9686 | 0.9661 / 0.9682 (mismatched line) | base 0.9668 / 0.9688 | 0.9668 |
| `cshrink_inf` (abd), **full** | ψ = ∞ housekeeping | 0.9607 / 0.9606 | 0.9631 (old-grid ψ) | selected ψ 0.9609 | 0.9620 |

**Read:**
- lwkrige v2 beats its cell-permuted null by 0.0031 and the base by 0.0032 (bar 0.005); δ = 0.99 was chosen in 40/48
  calls, and forced δ = 0.98 is 0.9524 vs null 0.9598. The gain is on TMT (leduc2022, gregoire2023_mixCTRL,
  px063590_ADbrain: 0.9678 → 0.9627, null 0.9676) more than DIA (0.9429 → 0.9417). The v1 pairwise-complete S was
  indefinite and never used (0.9554).
- plmemb: embeddings predict pool loadings (LOO R² 0.057–0.115 vs ≈ −0.14 permuted, `plmemb_diag.json`), but Gaussian
  embeddings reproduce most of the gain (0.9535 vs 0.9529); the real-over-null part is weiss2026_liver only.
- mllr: the permuted null wins (MCAR on gregoire2023_mixCTRL 0.9542 vs 0.9668), i.e. attenuating the reference helps
  some targets; A itself stays near I.
- linerna: the matched line loses to a mismatched line (0.9668 vs 0.9661); inner-mask correlation of the RNA-programme
  prediction with the held-out residual is ≈ 0.01.
- cshrink_inf on the full pool (MCAR 0.9704, MNAR 0.9509, 26/30 ≤ 0.99, worst woo2025_PC3 MCAR 1.0009) ties selected
  ψ; the TMT/carrier MCAR weak spot remains. The H15 champion stands; no holdout touched.
- Round 6 (directives in `docs/phaseC/KING_LOG.md`): `taperkrige` (curated-biology covariance localisation of the
  in-domain covariance), `kronmt` (Kronecker cell × protein GP), `otanalog` (optimal-transport analog reference cells)
  and `plexre` (TMT plex / batch random effects).
