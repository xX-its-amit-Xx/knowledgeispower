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
| ctprior | bio | marker-knowledge (HPA scRNA markers) cell-type-aware soft *cell* weighting of the reference in a weighted `fit_loadings` | dead (fast: 0.999; margin over null 0.0000) | 0.9990 / 0.9986 | 0.9990 / 0.9987 | not run | cand_solo 0.9988; weight 0.5521; 0/12 ≤ 0.99. τ = ½ median (`ctprior_half`): 0.9980 / 0.9982 (null 0.9982 / 0.9978); τ = ∞ drift control (`ctprior_inf`): 0.9998. Best target weiss2026_liver MNAR 0.9961 (null 0.9974). Weight shifts are the same on every target (generic per-dataset score scaling), not cell identity |
| copula | rep | nonparanormal (rank-Gaussian) encoding of each reference dataset and the target before fitting / projecting | dead (fast: 1.0005; above champion) | 1.0005 / 1.0011 | 1.0004 / 1.0004 | not run | cand_raw 1.0274, cand_solo 1.0193; weight 0.156; 0/12. `copula_refonly` (reference side only): 0.9980 / 0.9976, null 0.9991 / 0.9988, 1/12 (gregoire2023 MNAR 0.9881 vs null 0.9998), weight 0.323. `copula_tgtonly`: 1.0010 / 1.0014, null 1.0001, weight 0.146. Target-side empirical-quantile decode adds noise; TMT is the *most* Gaussian platform (`copula_diag.json`) |
| pepev | mod | peptide/precursor-level evidence as per-entry precision weights in the U-step of `adapt_project` | dead (fast: 0.9996; v2 0.9969, margin 0.0034) | 0.9996 / 0.9987 | 0.9998 / 0.9996 | not run | cand_solo 0.9988; weight 0.5729; 1/12. 5/6 targets aligned (px063590 skipped). α = 1: 1.0005 (null 1.0001); U+V weights: 0.9988 (null 0.9998); protein-permuted null 1.0001. Helps MCAR, hurts MNAR (weiss2026_liver 0.9759 / 1.0121): in DIA, N correlates with the value within protein (Spearman 0.21–0.44, `pepev_diag.json`). v2 (protein-level column-mean N, `pepev_v2`): 0.9969 / 0.9967, null 1.0003 / 1.0002 |
| nfbase | abd | Netflix-style hierarchical baselines: shrink the target's 30-cell protein mean/sd toward its technology group | dead as a biological prior (margin over null ≤ 0.0028); **machinery is the largest Phase-C effect** | 0.9868 / 0.9859 | 0.9867 / 0.9863 | not run | ψ = 1 main. cand_solo 0.9864; weight 0.7188; 6/12. **ψ = 4 (`nfbase_psi4`): 0.9663 / 0.9647, null 0.9691 / 0.9688, 12/12 ≤ 0.99, weight 0.7604.** technull 0.9868, flat 0.9869 (tech × protein term adds nothing). meanonly 0.9916 (null 0.9915), sdonly 0.9947 (null 0.9950); post hoc meanonly_psi4 0.9920, sdonly_psi4 0.9832 (null 0.9829). Real-over-null only on DIA MCAR (weiss2026_liver 0.9490 vs 0.9645). The null = James–Stein shrinkage of the target's own mean/log-sd toward a constant → round-3 `tgtshrink` |
| tgtshrink | abd | **champion revision (data-only):** EB/James–Stein shrinkage of the target's per-protein mean and log-sd inside `adapt_project`, ψ chosen on the inner mask, log-sd shrunk toward a limma-style trend in (mean, observed fraction) | proposed (round 3) | | | | null: trend covariates permuted across proteins (= shrink to a constant, nfbase's null); full pool if fast ≤ 0.99 |
| selmnar | rep | Heckman-style probit selection model: missing entries as informative (likely-low) in an EM U-step of `adapt_project` | proposed (round 3) | | | | null: the same EM with missing indicators permuted within each cell |
| seqprior | mod | protein-intrinsic features (UniProt sequence detectability: length, tryptic peptides, MW, GRAVY, pI; HPA subcellular location) predicting per-protein location / log-sd / detection rate, as protein-specific shrinkage targets | proposed (round 3) | | | | null: feature rows permuted across proteins |
| varprior | bio | turnover/buffering prior on cell-to-cell variability: per-protein log-sd shrinkage target from measured half-lives + complex membership | proposed (round 3) | | | | null: half-life and complex covariates permuted across proteins |

## Seed backlog (from Amit, Oct 23)
Round-1 status of the seeds: tissue/technology dataset weighting → `ttprior` (dead; technology axis carries a small
signal); per-protein adaptation strength → `ebshrink` (dead); scRNA atlas modality → `scrnaprior` (dead).
Round-2 status: tissue prior at cell-type resolution → `ctprior` (dead); peptide-level evidence → `pepev` (dead);
change of encoding → `copula` (dead). Round 3 takes the synthetic-features seed (predicted half-life, detectability) into
the target's standardisation step, where `nfbase` found the largest effect.

- **Tissue/technology-aware priors.** Weight reference datasets, or use separate loadings, by tissue and technology
  similarity to the target.
- **Per-protein adaptation strength.** Learn λ_v per protein (empirical Bayes) instead of a single value.
- **Change the encoding space.** Hyperbolic or higher-dimensional embeddings, kernels instead of linear factors.
- **Synthetic features.** Use simple sub-models to impute features (e.g., predicted abundance or half-life), then feed
  them in.
- **New modalities.** HPA imaging and subcellular localisation, scRNA co-expression atlases, peptide-level evidence.
- **Curated complexes inside the adapted model.** Use them as a graph penalty on V.
