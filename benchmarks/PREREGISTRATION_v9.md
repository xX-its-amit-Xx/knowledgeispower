# Preregistration v9: the H15 champion and its revision against a broad comparator panel, on new holdouts

Committed **before** any v9 holdout is imputed or scored. Amit asked for this on 2026-11-07:
- redo the earlier validation;
- add at least two more holdouts;
- compare against the methods researchers typically use, including graph-prior methods distinct from those already
  tried.

## 0. What was already done (no holdout was spent on it)

**The khan2023 test (v6, H14) was rerun with the frozen v6 code.**
- Setup: commit `411d439` with its 17-dataset reference corpus.
- Result: an exact reproduction. Every RMSE is identical, and the verdict is again **3/8, not confirmed**
  (`results/v9/v6_khan2023_reproduction.json`).
- Holdout record before v9:
  - khan2023 (H14): not confirmed;
  - brunner2022 (H14): confirmed;
  - wu2026 (H15): confirmed 16/16.
- The current champion (H15) has been tested on **one** holdout, wu2026.

**Comparator panel** (`src/scp_graphfs_mcp/baselines.py`), screened on Phase-C dev (6 targets × 2 masks × 2 seeds):

| Group | Methods |
|---|---|
| Data-only | mean; MinDet, MinProb; kNN; softImpute; iterative SVD; BPCA; MICE; missForest (reduced: 10 trees, 3 iterations); denoising autoencoder |
| Graph prior (STRING v12.5 ≥ 700) | Laplacian-regularised MF; the same with an equal-size data-derived graph (control); random walk with restart on kNN residuals; 2-layer GCN; inductive MF on spectral graph embeddings |
| Controls | GCN on a degree-preserving shuffle of STRING; GCN with no edges |

- PPCA (ML) was dropped: it diverges at 30 cells.
- Each method is scored raw and as `+stack` (simplex stack with kNN and low-rank, the same help the champion gets).
- Dev results are in `results/v9/dev_panel_dev_fast.json` and `dev_panel_devctl_fast.json`.
- The two strongest competitors, chosen by **mean dev RMSE (all entries) relative to the champion**, are frozen in
  `results/v9/v9_selection.json`:
  - best data-only = `gcn_nograph+stack` (**0.9917**: it beats the champion on dev);
  - best graph prior = `gcn+stack` (1.0143).
- **Amendment made before any holdout was run.** The first draft of this file named `mice+stack` (1.0275).
  - The GCN controls then showed that the same network with **no edges** beats it with STRING: 0.993 vs 1.042 raw,
    with no edges better on 11/12 dev splits.
  - The no-edge network is a protein-wise masked MLP across the subsample's cells. It uses no prior and no reference,
    so it is a data-only method.
  - STRING beats a degree-preserving shuffled STRING (1.042 vs 1.096), but loses to no graph.
  - Selection now includes it.

## 1. Methods under test (frozen)

- **`champion` (H15):** stack(kNN, low-rank, `adapt`), exactly as in PREREGISTRATION_v8.
  - Reference loadings: ALS rank 20, λ = 5, 15 iterations, learned on `explore_phaseb.reference_datasets()`, aligned to
    the holdout's proteins.
  - The 22 reference datasets exclude every holdout of every version, and the non-human woo2022.
  - λ_v = 100, 10 extra factors.
- **`revised` (Phase-C `cshrink`):** stack(kNN, low-rank, `adapt`, `cshrink`).
  - `cshrink` is selected-ψ constant shrinkage of the target's per-protein mean and log-sd inside the adapted
    projection (`src/scp_graphfs_mcp/phasec/cshrink.py`).
  - Dev full-pool result was 0.9609 vs the champion; it has never been scored on a holdout.

## 2. Holdouts

**Primary** (never used by any method in this repository):

| Holdout | Source | Cells × proteins | Missing |
|---|---|---|---|
| `derks2022` | plexDIA, mTRAQ | 125 × 1093 | 0.57 |
| `karagach2025_S1` | PXD058457, timsTOF Slice-PASEF | 258 × 3687 | 0.33 |
| `px061710_GSTP1KO` | PXD061710 | 78 × 3899 | |
| `px061710_MDAfix` | PXD061710 | 76 × 3325 | |
| `px061710_RKOfix` | PXD061710 | 75 × 5088 | |

- These are **3 independent projects**. The three px061710 sets share a lab and instrument, so they count as one
  project.
- Provenance and QC are in `docs/holdout_v9_datasets.md`.

**Secondary** (reported, not part of the verdict):

| Holdout | Why it is secondary |
|---|---|
| `khan2023`, `brunner2022` | Never seen by H15/`cshrink` development, but already used for the H14 tests, so we have seen results on them. |
| `karagach2025_S2` (39 cells), `px077946_H2O2` (40 cells) | Fewer than 60 cells, so 30-cell subsamples overlap heavily. |

## 3. Design

- Seeds 9000–9009.
- 30-cell **random** subsamples, no stratification: labels are technical or unverified on several holdouts.
- MCAR and MNAR 10% masks from `benchmarks/run.py::make_mask`, on proteins observed in ≥ 50% of the subsample.
- The stack uses a 10% inner mask of observed entries.
- Endpoints: RMSE on all masked entries, and RMSE on complex-member proteins.
- Script: `benchmarks/confirm_v9.py <holdout>`. It runs once per holdout, refuses uncommitted code, and reads the
  frozen selection file.

## 4. Hypotheses and pass rules

Every test is a 95% paired-bootstrap CI (1000 resamples, seed 17) of ΔRMSE. A test passes when the CI is entirely
below 0 **and** ΔRMSE ≤ −1% of the comparator's RMSE (the SESOI).

**H16 (the champion generalises and adds meaningful value).**
- Test `champion` against each of:
  - `stack_lowrank`, the best data-only stack;
  - `champion@null`, reference loadings permuted across proteins;
  - `champion@noref`, no reference;
  - `gcn_nograph+stack`;
  - `gcn+stack`.
- That is 5 comparators × 2 endpoints × 2 masks = **20 tests per holdout**.
- A holdout confirms H16 if all 20 pass.
- A project confirms if its single set confirms; for px061710, if ≥ 2 of its 3 sets confirm.
- **H16 is confirmed if ≥ 2 of the 3 primary projects confirm.**

**H18 (the revision beats the strongest typical methods).**
- Test `revised` against `stack_lowrank`, `gcn_nograph+stack` and `gcn+stack`: 3 × 2 × 2 = **12 tests per holdout**.
- The project rules are the same as for H16. **H18 is confirmed if ≥ 2 of the 3 primary projects confirm.**

**H17 (the revision adds to the champion).**
- Test `revised` against `champion`: 4 tests per holdout.
- The project rules are the same as for H16. **H17 is confirmed if ≥ 2 of the 3 primary projects confirm.**

**Reported for every holdout, without a verdict:**
- the champion against every panel method, raw and stacked, with CIs;
- the GCN graph controls: does STRING beat the shuffled graph and no graph?

**Prediction, written now:**
- **H16 is expected to fail** against `gcn_nograph+stack`, which already beats the champion on dev (0.9917).
- H16 against `gcn+stack` is marginal: the dev margin is only 1.4%.
- The claim most likely to survive is H18: the revision is 4.5% better than the champion on dev and better than
  `gcn_nograph+stack` on every dev split.
- Any failure will be reported as such, not reinterpreted.
