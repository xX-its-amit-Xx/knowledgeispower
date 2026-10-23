# seqprior: sequence and localisation features as per-protein shrinkage targets

- **Queen:** new data modalities (round-3 directive in `KING_LOG.md`).
- **Status (worker's recommendation):** **dead** at the fast screen. It clears 0.99 vs the champion on all 12 splits,
  but it ties the shrink-to-constant baseline and its own permuted-feature null.

  | Run (rmse_all, fast pool) | cand_stack | null_stack | const cand_stack | Margin vs min(null, const) |
  |---|---|---|---|---|
  | `seqprior` (GBM) | 0.9686 | 0.9687 | 0.9686 | 0.0000 |
  | `seqprior_ridge` (post hoc) | 0.9683 | 0.9687 | 0.9686 | 0.0003 |

  The pass bar was a margin of at least 0.005, so the full pool was not run.
- **Code:** `src/scp_graphfs_mcp/phasec/seqprior.py`. It imports `ref_baselines`, `adapt_project_mb`, `_prior_cached`,
  `_load`, `MIN_N_CAL`, `LAM_V`, `EXTRA` and `FORBIDDEN` from `nfbase`, and `META` from `ttprior`. `baselines` is a
  local copy with a multi-regressor calibration.
- **Results** (in `results/phaseC/`): `seqprior_fast.json`, `seqprior_const_fast.json`, `seqprior_ridge_fast.json`,
  and `seqprior_diag.json` (sub-model OOF R², coverage, and fusion coefficients and tau2 on seed 100).
- **Cache** (in `data/cache/phasec/seqprior/`, not committed):
  - The two downloads, plus response headers and the HPA licence page.
  - `features_by_acc.json`.
  - Per target: `<t>_feat.npz`, `<t>_of.npz` and `<t>_pred_{true,null,ridge,ridge_null}.npz`.
  - `unit_check.json` and `log_*.txt`.

## Idea
A protein's mean abundance relative to the dataset median, its cell-to-cell log-sd and its detection rate are partly
fixed by its sequence and location. A sub-model trained on the LODO reference predicts all three for every target
protein, including proteins the reference does not cover. Its predictions replace the constant that nfbase's ψ = 4
null shrinks toward.

The expected gain was largest on DIA MNAR targets.

## Mechanism (as built, following the mod spec)
**A. Data.** See the provenance table below. The downloads are parsed as data only, with gzip, zipfile and string
splitting. Nothing from the cache directory is executed. The first parse ran under `python -I` from the repository
root.

**B. Features.** There are 29 columns per target protein: the first accession from `T._target_accs`, with any `-n`
isoform suffix stripped.
- Sequence (7):
  - log length.
  - log(1 + the number of unique fully tryptic peptides of 7–30 aa). These use the Keil rules: cleave after K/R, not
    before P, with 0 missed cleavages.
  - log average MW: the sum of average residue masses + 18.015.
  - GRAVY, the Kyte–Doolittle mean.
  - pI by bisection with the EMBOSS pKa values.
  - frac K+R and frac C.
- Localisation (20): HPA 'Main location' multi-hot over the 20 most frequent terms. It is joined on the FASTA `GN=`
  name.
- Flags: has_hpa and has_seq.

Missing values get the column median, and every column is z-scored over the target's proteins. Coverage
(`seqprior_diag.json`):

| Target | p | Sequence coverage | HPA coverage |
|---|---|---|---|
| petrosius2023 | 2813 | 1.000 | 0.834 |
| leduc2022 | 2844 | 1.000 | 0.847 |
| weiss2026_liver | 4520 | 1.000 | 0.791 |
| gregoire2023_mixCTRL | 1878 | 0.998 | 0.848 |
| px063590_ADbrain | 5140 | 0.899 | 0.726 |
| sadiku2025 | 3039 | 1.000 | 0.797 |

Sanity check: HBA1 (P69905) gives length 142 and MW 15,258 Da.

**C. Responses.** The source is `nfbase.ref_baselines(t, names)`, which is LODO and asserts that no holdout is in it.
- Responses are y (relative mean), l (log sd) and the observed fraction.
- The observed fraction is computed from `_load(d)` and is NaN where the protein is absent from d. It is cached as
  `<t>_of.npz`.
- Each response is the protein-wise nanmean over the datasets in the target's acquisition group. The observed fraction
  enters as logit(clip(of, 0.02, 0.98)).
- Each protein is weighted by the number of datasets covering it.
- Every fast target had ≥ 2 group datasets and ≥ 500 covered proteins, so the fallback to all datasets was never used:
  - 9 group datasets for the DIA targets.
  - 5 for the TMT targets.

**D. Sub-models.** One HistGradientBoostingRegressor per response, with the spec's hyperparameters.
- Covered proteins get the 5-fold out-of-fold prediction (KFold, shuffle, seed 0).
- Uncovered proteins get the full-fit prediction.
- The post-hoc ridge variant uses RidgeCV with alphas logspace(−2, 3, 12).
- Predictions are disk-cached per (target, mode) and also lru-cached.

**E. Fusion.**
- Location: OLS of ybar on [1, py, po]. Scale: OLS of lbar on [1, pl, po]. Both are fitted over the cal proteins
  (base_ok & n ≥ 5).
- Then nfbase's method-of-moments tau2 and BLUP at ψ = 4, applied to every base_ok protein.
- Then `adapt_project_mb(lam_v = 100, extra = 10)`.

**F. Entry points.**
- `impute`: real features.
- `impute_null`: feature rows permuted across the target's proteins with `np.random.default_rng(777)`, before both
  training and prediction.
- `impute_const`: design [1], the shrink-to-constant baseline.
- `impute_ridge` and `impute_ridge_null`: the post-hoc ridge variant.
- `impute_nf` and `impute_nf_null`: the optional variant with the nfbase reference prior as an extra regressor. These
  are implemented but **not run**, because sequence alone added nothing.

**G. Unit checks.**
- ψ = 0 reproduces `transfer.adapt_project` on MCAR seed 100 of all 6 fast targets. The max |diff| is 1.1e-15 to
  5.6e-15, so all pass < 1e-10 (`unit_check.json`).
- `impute_const` gives cand_stack 0.9686 (all) / 0.9692 (in-complex). nfbase_psi4's null_stack is 0.9691 / 0.9688, so
  the differences are 0.0005 and 0.0004, inside the 0.002 tolerance. 0.9686 also equals `varprior_const`.

## Null control
`impute_null` is the identical pipeline (features, GBM, OOF, 3-column calibration, ψ = 4 BLUP, adapt_project_mb), with
the feature rows permuted across proteins. The required second baseline is `impute_const`.

Credit = min(null_stack, const cand_stack) − cand_stack.

## Data provenance (for the king to merge into `data/PROVENANCE.md`)
| File | URL | Licence | Release / version | Downloaded | Bytes | SHA-256 |
|---|---|---|---|---|---|---|
| `uniprot_human_sp.fasta.gz` | `https://rest.uniprot.org/uniprotkb/stream?format=fasta&compressed=true&query=%28organism_id%3A9606%29+AND+%28reviewed%3Atrue%29` | CC BY 4.0 (UniProt) | x-uniprot-release 2026_03 (release date 02-September-2026); 20,431 entries | 2026-10-23 | 7,508,673 | dc8ef8cbf0fc1b857b3be9a559613533ec49824079fb057e31fdfcf578cb229a |
| `subcellular_location.tsv.zip` | `https://www.proteinatlas.org/download/tsv/subcellular_location.tsv.zip` | CC BY 4.0 (see note) | last-modified 2025-11-06; 13,603 gene rows (13,597 unique gene names) | 2026-10-23 | 252,331 | a9466882bd3328416bcde77c32a6bcb2a17d70183af90e98d032f6f5f781d504 |

**Licence note.** The directive said CC BY-SA 4.0 for HPA. The HPA licence page (https://www.proteinatlas.org/about/licence,
fetched 2026-10-23) states: "The Human Protein Atlas is licensed under the Creative Commons Attribution 4.0
International License for all copyrightable parts of our database." So CC BY 4.0 is recorded here. The HPA SHA-256
matches the hash the queen recorded at her fetch.

## Results (exact numbers from the JSON, rmse_all mean_rel_to_champion)

**Means.**

| Run | cand_raw | cand_solo | cand_stack | null_stack | Splits ≤ 0.99 | Stack weight |
|---|---|---|---|---|---|---|
| `seqprior` | 0.9666 | 0.9687 | 0.9686 | 0.9687 | 12/12 | 0.75 |
| `seqprior_const` (null = seqprior null) | 0.9664 | 0.9686 | 0.9686 | 0.9687 | 12/12 | 0.75 |
| `seqprior_ridge` (null = ridge null) | 0.9661 | 0.9682 | 0.9683 | 0.9687 | 12/12 | 0.75 |

In-complex cand_stack:
- `seqprior`: 0.9689 (null 0.9694).
- `seqprior_const`: 0.9692.
- `seqprior_ridge`: 0.9687 (null 0.9694).

**Per target × mask, `seqprior`.** Each cell is cand_raw / cand_solo / cand_stack / null_stack. The last column is the
`seqprior_const` cand_stack.

| Target | Mask | seqprior | const |
|---|---|---|---|
| gregoire2023_mixCTRL | MCAR | 0.9724 / 0.9724 / 0.9724 / 0.9724 | 0.9721 |
| gregoire2023_mixCTRL | MNAR | 0.9374 / 0.9442 / 0.9442 / 0.9446 | 0.9443 |
| leduc2022 | MCAR | 0.9853 / 0.9779 / 0.9779 / 0.9781 | 0.9781 |
| leduc2022 | MNAR | 0.9676 / 0.9678 / 0.9678 / 0.9683 | 0.9681 |
| petrosius2023 | MCAR | 0.9878 / 0.9878 / 0.9814 / 0.9816 | 0.9817 |
| petrosius2023 | MNAR | 0.9843 / 0.9843 / 0.9802 / 0.9781 | 0.9781 |
| px063590_ADbrain | MCAR | 0.9833 / 0.9783 / 0.9791 / 0.9792 | 0.9790 |
| px063590_ADbrain | MNAR | 0.9793 / 0.9798 / 0.9831 / 0.9831 | 0.9830 |
| sadiku2025 | MCAR | 0.9763 / 0.9784 / 0.9784 / 0.9800 | 0.9798 |
| sadiku2025 | MNAR | 0.9398 / 0.9670 / 0.9670 / 0.9668 | 0.9672 |
| weiss2026_liver | MCAR | 0.9665 / 0.9665 / 0.9665 / 0.9662 | 0.9659 |
| weiss2026_liver | MNAR | 0.9195 / 0.9195 / 0.9254 / 0.9264 | 0.9264 |

The predicted DIA-MNAR signature did not appear:
- petrosius2023 MNAR is *worse* than const: 0.9802 vs 0.9781.
- sadiku2025 MNAR ties: 0.9670 vs 0.9672.
- weiss2026_liver MNAR is 0.0010 better: 0.9254 vs 0.9264.

The ridge variant's best split is weiss2026_liver MNAR, 0.9242 vs a null of 0.9265.

**Sub-model OOF R²** (weighted, covered proteins; `seqprior_diag.json`). Each cell is y / l / o.

| Target | GBM real | GBM null | Ridge real | Ridge null |
|---|---|---|---|---|
| petrosius2023 | 0.149 / 0.188 / 0.060 | −0.103 / −0.114 / −0.110 | 0.153 / 0.159 / 0.134 | −0.018 / −0.013 / −0.005 |
| leduc2022 | −0.062 / −0.016 / 0.147 | −0.132 / −0.130 / −0.102 | 0.028 / 0.041 / 0.149 | −0.005 / −0.019 / −0.013 |
| weiss2026_liver | 0.168 / 0.151 / 0.098 | −0.078 / −0.067 / −0.073 | 0.139 / 0.123 / 0.131 | −0.009 / −0.010 / −0.012 |
| gregoire2023_mixCTRL | 0.006 / 0.087 / 0.077 | −0.156 / −0.154 / −0.153 | 0.043 / 0.065 / 0.111 | −0.015 / −0.028 / −0.014 |
| px063590_ADbrain | −0.029 / −0.006 / 0.144 | −0.126 / −0.096 / −0.124 | 0.031 / 0.030 / 0.165 | −0.019 / −0.013 / −0.012 |
| sadiku2025 | 0.149 / 0.158 / 0.120 | −0.146 / −0.088 / −0.103 | 0.140 / 0.139 / 0.157 | −0.013 / −0.025 / −0.009 |

## Honest read
**The modality carries real but modest signal.**
- Sequence and localisation explain 14–19% of the reference's per-protein relative mean and log-sd on DIA targets.
  On TMT targets they explain about 0–9%.
- They explain 6–17% of the detection logit.
- The permuted null is negative everywhere.
- In the target's own calibration the prior is used:
  - On DIA targets the location slope on py is 0.55–0.99.
  - The scale slope on pl is 0.65–0.87 on DIA and 0.13–0.47 on TMT.
  - Under the null these slopes are near 0, except gregoire2023 scale, where the slope on the permuted pl is about 0.40.

  These coefficients are in `seqprior_diag.json` → `fusion_seed100`.

**But it cannot move imputation RMSE, for two structural reasons visible in that diagnostic:**
1. **Location is not where shrinkage acts.** The weights below are in `seqprior_diag.json` → `loc_weights_seed100`
   and `fusion_seed100`.
   - On the DIA targets the between-protein tau2 for ybar is 1.77–2.45 under the real prior. The median s²/n of an
     evaluated protein is 0.022–0.057. So the mean BLUP weight on the target's own mean is 0.88–0.93, and the prior
     barely enters the location estimate. The prior lowers tau2 by 4–11% relative to const, which changes little.
   - On the TMT targets the location is shrunk almost fully (mean w 0.003 for leduc2022 and px063590, 0.75 for
     gregoire2023). But there the calibrated slope on py is about 0 (−0.04 to 0.13), because ratio-normalised TMT data
     have no between-protein location for the prior to explain.
2. **Few of the evaluated proteins are sparse.** The harness masks only proteins with ≥ 50% of the 30 cells observed;
   the median n of masked columns is 11–24. So the sparse, reference-uncovered proteins the idea targets are rarely
   evaluated. In addition, the reference already covers 51–95% of each target's proteins.
   - For scale, tau2 sits at the 1e-4 floor on most targets, so sd is set entirely by the calibrated prior. A prior
     that explains about 15% of between-protein log-sd variance then changes sd_j by a few percent, and that is
     absorbed by the low-rank fit and the stack.

The result is a 0.0000 margin (GBM) or 0.0004 margin (ridge, post hoc) over the constant. This repeats varprior's
finding: a predicted per-protein *shrinkage target* is not the lever. The lever is the shrinkage itself, toward any
constant.

No revision (v2) was run: there was no bug or scale error to fix. The ψ = 0 check passes, and the const run reproduces
nfbase_psi4's null.

## Next step
- Close the "protein-specific shrinkage target" line. nfbase (reference), varprior (turnover and complexes) and
  seqprior (sequence and localisation) have all tied shrink-to-constant at ψ = 4.
- The features may still have value where the target has *no* own estimate, that is, as a cold-start prior for
  proteins with n < 2 in the target. The fast harness does not evaluate that case. It would need a different benchmark
  (masking whole proteins) before it is worth building.
- The cached feature table (`features_by_acc.json`, 20,431 accessions × 7 sequence features plus HPA) is reusable by
  any later idea that needs protein-intrinsic covariates. One example is a detectability term in rep's `selmnar`
  probit (c_j ~ features).
