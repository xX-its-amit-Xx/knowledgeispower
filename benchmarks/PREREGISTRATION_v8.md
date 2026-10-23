# Preregistration v8: Phase-B holdout wu2026 (Amendment 15, H15)

Committed **before** wu2026 is evaluated. wu2026 (prenatal human brain, 2268 cells, DIA, 90% missing) was used only as
part of an external reference in Amendments 11–14. It was never imputed or scored, and it has been in no reference since
Amendment 15.

**Method `adapt`** (fixed by Phase-B dev, `docs/EXPLORATION_LOG.md`):
- Reference loadings: ALS rank 20, λ = 5, 15 iterations, ≥ 20 observations. They are learned on the single-cell corpus
  `explore_phaseb.reference_datasets()`, which excludes every holdout and the non-human woo2022. No bulk atlases.
- `transfer.adapt_project` with λ_v = 100 and 10 extra free factors.
- The result is stacked with kNN and low-rank (simplex weights, step 0.25), chosen on a 10% inner mask of observed
  entries.

**Data:** `data/cache/corpus/wu2026.npz`, rebuilt by `scripts/corpus_ext/build_wu2026.py` (log2, per-cell median-centred).

**Design:**
- 30-cell random subsamples (wu2026 has no per-cell labels), seeds 8000–8009.
- MCAR and MNAR masks from `benchmarks/run.py::make_mask`, on proteins observed in ≥ 50% of the subsample.

**Endpoints:** RMSE on all masked entries (primary) and on complex-member proteins (primary).

**Comparators:**
- kNN;
- the kNN + low-rank stack (best data-only);
- `adapt` with the reference loadings permuted across proteins (scrambled reference);
- `adapt` with no reference (zero prior, same free factors).

**Pass:** for every comparator × endpoint × mask (16 tests), the 95% paired-bootstrap CI (1000 resamples, seed 17) of
ΔRMSE is entirely < 0, and ΔRMSE ≤ −1% of the comparator's RMSE.

**Outcomes:**
- **H15 is confirmed** if all 16 pass.
- **Partial:** the result is reported test by test.

**Single use:** `confirm_v8.py` refuses to run twice or with uncommitted code. derks2022 stays reserved.
