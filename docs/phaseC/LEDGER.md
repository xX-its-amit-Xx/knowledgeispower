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
| ttprior | bio | tissue/technology-aware reference prior: metadata-kernel dataset weights in a weighted `fit_loadings`, fed to `adapt_project` (FULL kernel, τ = 1) | dead (fast: 0.9951 > 0.99; margin over null 0.0040 < 0.005) | 0.9951 / 0.9953 | 0.9991 / 1.0005 | not run | **best of round 1.** cand_solo 0.9944; stack weight 0.510; 1/12 splits ≤ 0.99 (sadiku2025 MCAR 0.9797). Variants (`ttprior_<v>_fast.json`): full τ=2 0.9975 (null 0.9998), full τ=0.5 0.9953 (null 0.9987), tissue-only 0.9992 (null 0.9997), tech-only 0.9970 (null 0.9980), two-block 0.9965 (null 0.9984). Signal is **technology**, concentrated on Astral-DIA targets; tissue-only = its null |
| ebshrink | rep | per-protein λ_v,j by empirical Bayes from cross-dataset loading variance τ²_j | dead (fast: 0.9986; margin 0.0014) | 0.9986 / 0.9983 | 1.0000 / 0.9999 | not run | cand_solo 0.9991; weight 0.29; 1/12 ≤ 0.99 (gregoire2023 MNAR 0.9876). v2 (λ ∝ 1/τ⁴): 0.9988 / 0.9985, null 1.0002, weight 0.17. Per-protein anchoring is as saturated as the scalar λ_v |
| scrnaprior | mod | HPA single-cell-type RNA SVD (k = 8) as an extra ridge-penalised factor block in `adapt_project` (λ_s = 30) | dead (fast: 0.9973; margin 0.0014) | 0.9973 / 0.9974 | 0.9987 / 0.9988 | not run | cand_solo 0.9989; weight 0.260; 0/12 ≤ 0.99; best weiss2026_liver 0.9903 / 0.9908. λ_s = 100: 0.9976 (null 0.9989). Zero block (S = 0): 0.9998. Ordering real < permuted < zeros holds, but the effect is ~0.1–0.3% |
| dca | abd | direct-coupling (shrinkage-inverse precision) Gaussian conditional-mean imputation on the top 1,500 proteins | dead (fast: 1.0004; stack weight ≈ 0) | 1.0004 / 1.0006 | 1.0000 / 1.0000 | not run | cand_solo 1.032 (beats stack_lowrank 1.0773, so real signal, but fully covered by the adapted projection); weight 0.0104; λ = 0.85 (grid edge) on all 6 targets. CORUM top-30 partner fraction (`dca_couplings_diag.json`, mean of 6): DCA 0.0799, null 0.0177, raw corr 0.1063, random 0.0169 — inversion recovers complexes *worse* than raw correlation |
| ctprior | bio | marker-knowledge (HPA scRNA markers) cell-type-aware soft *cell* weighting of the reference in a weighted `fit_loadings` | proposed (round 2) | | | | null: size- and observation-matched random marker sets |
| copula | rep | nonparanormal (rank-Gaussian) encoding of each reference dataset and the target before fitting / projecting | proposed (round 2) | | | | null: per-protein marginal transforms reassigned to random other proteins |
| pepev | mod | peptide/precursor-level evidence as per-entry precision weights in the U-step of `adapt_project` | proposed (round 2) | | | | null: weights permuted among each cell's observed entries; fallback HPA IF compartments |
| nfbase | abd | Netflix-style hierarchical baselines: shrink the target's 30-cell protein mean/sd toward its technology group | proposed (round 2) | | | | null: protein labels permuted in reference baselines; 2nd null: tech labels permuted |

## Seed backlog (from Amit, Oct 23)
Round-1 status of the seeds: tissue/technology dataset weighting → `ttprior` (dead; technology axis carries a small
signal); per-protein adaptation strength → `ebshrink` (dead); scRNA atlas modality → `scrnaprior` (dead). The rest are
not yet built.

- **Tissue/technology-aware priors.** Weight reference datasets, or use separate loadings, by tissue and technology
  similarity to the target.
- **Per-protein adaptation strength.** Learn λ_v per protein (empirical Bayes) instead of a single value.
- **Change the encoding space.** Hyperbolic or higher-dimensional embeddings, kernels instead of linear factors.
- **Synthetic features.** Use simple sub-models to impute features (e.g., predicted abundance or half-life), then feed
  them in.
- **New modalities.** HPA imaging and subcellular localisation, scRNA co-expression atlases, peptide-level evidence.
- **Curated complexes inside the adapted model.** Use them as a graph penalty on V.
