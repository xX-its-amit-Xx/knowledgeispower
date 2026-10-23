# Phase C idea ledger

**Status values:**
- `proposed`
- `building`
- `screened-fast`
- `screened-full`
- `dead`: with the reason
- `promoted`

All numbers are mean ratios vs the H15 champion (< 1 is better), from `results/phaseC/<id>_<pool>.json`.

| id | queen | idea (one line) | status | fast cand_stack (all / in-cx) | null_stack | full | notes |
|---|---|---|---|---|---|---|---|
| sanity | — | the champion fed back as a candidate (harness check) | screened-fast | 1.000 / 1.000 | 1.000 (noise; weight 0) | | harness OK; champion = 0.928 of stack_lowrank on the fast pool; 107 s per screen |

## Seed backlog (from Amit, Oct 23 · not yet built)
- **Tissue/technology-aware priors.** Weight reference datasets, or use separate loadings, by tissue and technology
  similarity to the target.
- **Per-protein adaptation strength.** Learn λ_v per protein (empirical Bayes) instead of a single value.
- **Change the encoding space.** Hyperbolic or higher-dimensional embeddings, kernels instead of linear factors.
- **Synthetic features.** Use simple sub-models to impute features (e.g., predicted abundance or half-life), then feed
  them in.
- **New modalities.** HPA imaging and subcellular localisation, scRNA co-expression atlases, peptide-level evidence.
- **Curated complexes inside the adapted model.** Use them as a graph penalty on V.
