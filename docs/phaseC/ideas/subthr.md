# subthr: forced photometry on naturally missing entries (round 4)

- **Queen:** new data modalities (analogy from abd: astronomy forced photometry).
- **Status:** **dead** on the fast screen. cand_stack 0.9647, which is above the null (0.9609) and the base (0.9640). The bar is ≥ 0.005 below both.
- **Code:** `src/scp_graphfs_mcp/phasec/subthr.py`.
- **Results:**
  - `results/phaseC/subthr_fast.json` (main)
  - `subthr_base_fast.json` (κ = 0, same code)
  - `subthr_k1_fast.json` (fixed κ = 1, diagnostic)
  - `subthr_kappa.json` (ψ/κ selections)
  - `subthr_diag.json` (unit check, evidence per split)
  - `subthr_build.json` (audit/build).

## Idea
In astronomy, a source that falls below the detection threshold is still measured at its known position and kept as a low-precision observation. Here the FDR filter plays the role of the threshold. Sub-threshold search output at a (cell, protein) entry that the published matrix leaves NaN is added as a weighted pseudo-observation. These pseudo-observations feed two things:
- the target's per-protein mean and log-sd, which the base's Fay–Herriot shrinkage acts on;
- the U-step of the adapted projection.

Evidence never sits at a scored entry: masked entries were observed, and evidence exists only where X is NaN. This is asserted at build time and again at call time.

## Audit (repeated by the worker; numbers from `results/phaseC/subthr_build.json`)
- **weiss2026_liver, petrosius2023:** no sub-threshold evidence. I accepted the queen's audit and did not re-read the 1.97 GB zip. The petrosius export is pre-filtered at EG.Qvalue < 0.01 (pepev's file). I did not re-measure either.
- **leduc2022 (TMT pSCoPE):** quantitative evidence exists in the 134 per-batch PSM assays of `leduc2022_pSCoPE.Rda`.
  - There are 822,823 PSMs in total, all with dart_qval ≤ 0.1.
  - 91,102 PSMs have 0.01 < dart_qval ≤ 0.1. Of these, **48,110** remain after removing Reverse/contaminant hits; about half of the sub-threshold PSMs are decoys or contaminants.
  - 39,718 of those PSMs map to a target protein.
  - 132 batches contain target cells.
  - Result: **78,749 evidence entries** (2.72% of X's NaNs). **38,979 of them have P > 0** (1.35% of X's NaNs), because many sub-threshold PSMs have dart_PEP = 1.
  - Scale check on passing PSMs (4,746,748 PSM×cell pairs; diagnostic only, nothing fitted):
    - Spearman(within-batch PSM deviation, X's within-batch deviation) = 0.474 with the per-cell offset correction, 0.386 without it.
    - PSM deviations are wider than X's: MAD ratio 1.55; sd 0.604 vs 0.428.
- **sadiku2025 (Spectronaut, Astral DIA):** detection-score evidence only. Every entry has a run-wise PG.PEP (fraction 1.0).
  - On the corpus matrix (314 × 3039; the queen's figures were on the raw pivot):
    - Max PEP among quantified entries: 0.1416.
    - PEP quantiles at missing entries: p1 0.0103, p10 0.1607, p50 0.8141.
    - Fraction of missing entries with PEP < 0.05 / 0.15 / 0.5: 0.0326 / 0.0929 / 0.2951.
    - Within-protein Spearman(PEP, value) median: −0.211.
  - Result: 643,068 positions, of which 189,750 have P > 0 (29.5% of the NaNs).
  - Filter settings re-fetched (`setup.txt`): precursor Q-value cutoff 0.01, precursor PEP 0.15, protein Q-value 0.01 at both run and experiment level, protein PEP 0.15, imputation None.
- **gregoire2023_mixCTRL:** not included.
  - 14,492 sage PSMs have protein_fdr > 0.01.
  - Only **2** (PSM, channel) entries fall on a protein that is NaN in that batch's `proteins_<batch>_norm_log`.
  - 0 entries fall on a protein that is absent from the batch assay but present in the target names.
  - That is ≤ 0.504 evidence entries per 30-cell split, below the ≥ 1 rule.
- **Substitution, stated openly:** the directive listed weiss2026_liver / petrosius2023 / sadiku2025. The two targets that actually carry evidence are leduc2022 (TMT, not on the directive's list) and sadiku2025.
  - Fewer than 5 dev targets have evidence, so the full pool is not allowed.

Per 30-cell split (`subthr_diag.json`, mean over the 4 evaluation calls in `subthr_kappa.json`):
- leduc2022: 614 evidence entries (sum P ≈ 234) against 42,611–45,521 natural NaNs on usable columns.
- sadiku2025: 15,469 evidence entries (sum P ≈ 8,189) against ~47,000–50,000 NaNs.

## Mechanism
Implemented as in the queen's spec, with one stated deviation for leduc. The leduc PSM deviation (log2 reporter intensity minus the batch's median over target cells) also has a **per-cell offset** subtracted: the median of the same deviation over that cell's passing PSMs (dart_qval ≤ 0.01). The reason is that X is per-cell median-centred, and this correction raised the passing-PSM Spearman from 0.386 to 0.474. It uses raw data only and fits no parameter on X.

The steps:
- **Evidence values:**
  - leduc: D = median over the entry's sub-threshold PSMs of that deviation; P = 1 − min dart_PEP.
  - sadiku: P = max(0, 1 − PEP/0.5).
- **Pseudo-values** (from Xm only):
  - leduc: colmean_j + D.
  - sadiku: the protein's 10th percentile (the cell's 10th percentile if n_j < 5).
  - Weights π = κP.
- **Weighted model:**
  - Weighted target stats: n_eff = n + Σπ; weighted mean and sd; vy = s²/n_eff; vl = 1/(2 max(n_eff − 1, 1)).
  - These feed an exact copy of `tgtshrink.shrink(mode='null')`.
  - The U-step of `nfbase.adapt_project_mb` is weighted over observed ∪ evidence entries. The V-step uses observed entries only.
- **Selection** on `tgtshrink.selection_mask` (50/50 MCAR/MNAR, seed + 31337): first ψ with κ = 0, as in the base; then κ ∈ {0, 0.1, 0.3, 1} at that ψ, with ties going to the smaller κ.
- **Unit check** (`subthr_diag.json`): with κ = 0 the code reproduces `tgtshrink.impute_null` with max |diff| = **0.0** on leduc2022 and on sadiku2025 (MNAR, seed 100). The weighted projection with π = 1e-12 differs from the unweighted one by 5.4e-13 (leduc) and 1.1e-9 (sadiku).

## Null and baselines
- **Null (`impute_null`):** within each cell, the (P, D) records are permuted among that cell's evidence positions on usable columns (seed + 777 + i). The pseudo-value is then recomputed for the protein each record lands on. Positions, the weight distribution, κ selection and all machinery are kept.
- **Base (`subthr_base`):** κ forced to 0, same code; this equals `tgtshrink.impute_null`, the cshrink stand-in.

## Results (fast pool restricted to leduc2022 and sadiku2025 × MCAR/MNAR × seeds 100/101)
Values are rmse_all `mean_rel_to_champion`; < 1 is better.

| run | cand_raw | cand_solo | cand_stack | null_stack | in-complex cand_stack | weight | splits ≤ 0.99 |
|---|---|---|---|---|---|---|---|
| `subthr` (κ selected) | 0.9565 | 0.9631 | **0.9647** | **0.9609** | 0.9628 (null 0.9598) | 0.75 | 4/4 |
| `subthr_base` (κ = 0) | 0.9592 | 0.9640 | **0.9640** | — | 0.9637 | 0.75 | 4/4 |
| `subthr_k1` (κ = 1) | 0.9565 | 0.9638 | 0.9668 | 0.9637 | 0.9662 (null 0.9625) | 0.625 | 4/4 |

Per target and mask (`per_target_mask` in each JSON):

| target / mask | cand_raw | base cand_raw | cand_stack | null_stack | base cand_stack |
|---|---|---|---|---|---|
| leduc2022 MCAR | 0.9858 | 0.9861 | 0.9784 | 0.9788 | 0.9784 |
| leduc2022 MNAR | 0.9694 | 0.9693 | 0.9695 | 0.9694 | 0.9694 |
| sadiku2025 MCAR | 0.9685 | 0.9624 | 0.9700 | 0.9662 | 0.9641 |
| sadiku2025 MNAR | **0.9021** | 0.9190 | 0.9409 | 0.9293 | 0.9439 |

Mean per mask (computed from the four rows above):
- MCAR: cand_raw 0.9771 vs base 0.9743.
- MNAR: cand_raw 0.9358 vs base 0.9442.

**MCAR cand_raw is 0.0028 worse than the base on average and 0.0061 worse on sadiku2025.** The sadiku MCAR difference breaks the "no worse than base by > 0.005" rule.

Chosen κ (`subthr_kappa.json`, 8 calls per target = 4 splits × eval/inner):
- Main: leduc {0: 2, 0.3: 1, 1: 5}; sadiku {0.1: 1, 0.3: 5, 1: 2}.
- Null: leduc {0: 4, 0.3: 2, 1: 2}; sadiku {0.3: 3, 1: 5}.
- ψ: leduc mostly 2; sadiku mostly 16 (6/8).

Bar: cand_stack ≤ 0.99 holds (0.9647). But it is 0.0038 **above** the null and 0.0007 above the base, where it needed to be ≥ 0.005 below both. **Fails; no full pool** (also barred by the < 5 evidence targets rule). No v2 run: the failure is not a bug or a scale error (see below).

## Honest read
- **leduc2022: no effect (±0.0004 everywhere).** The quantitative evidence is real: passing-PSM deviations track X with Spearman 0.47. But it covers only 614 of 42,611–45,521 natural missing entries per split, about 1.4%, with mean P ≈ 0.38. A protein's mean/sd and a cell's U-scores barely move. Sub-threshold TMT PSMs are too few to matter. About half of the 91,102 PSMs in the sub-threshold window are decoys or contaminants, and many of the rest have dart_PEP = 1.
- **sadiku2025: the gain is a missing-is-low tilt, not detection content.** Adding near-floor pseudo-values on MNAR moves cand_raw from 0.9190 to 0.9021, but the protein-permuted null does as well or better once stacked (null_stack 0.9293 vs cand 0.9409). The selection also prefers larger κ for the null (κ = 1 in 5/8 calls). On MCAR the evidence hurts: cand_raw 0.9685 vs base 0.9624. This is selmnar's pattern again (MNAR help, MCAR harm, the random-inner-mask stack discounting it), weaker because κ is calibrated on the 50/50 inner mask.
- **Why the content does not help:** in the real assignment, the low-PEP ("almost detected") entries get the floor pseudo-value. Those are exactly the missing entries whose true value is likely *higher* than a random hole. The high-PEP ("truly absent") entries, where a floor value would be right, get weight 0. Spreading the same weights over random holes (the null) is a better "missing = low" prior than the real assignment. A PEP-only detection score therefore cannot be used as a floor observation. Using it would require a value model conditioned on PEP, which is a calibrated per-entry MNAR model: rep's `mnarjoint` territory, not new information.
- Both channels add real information, but neither is large enough or correctly signed to move a 30-cell estimate. The base pipeline itself (0.9640 on these two targets) stays the result to beat.

## Next step (for the queen)
- Close forced photometry on published search outputs: TMT sub-threshold PSMs are too sparse, and DIA PEP is the wrong sign as a floor. The only sources with sub-threshold *quantities* are raw-level artefacts, such as the 27.4 GB Spectronaut `.sne` or a DIA-NN rerun at a loose q-value. Both are out of the 2 GB budget.
- If the queen wants to keep the PEP signal, the right encoding is PEP as a covariate of the *missingness model*: P(missing | value, PEP), with the probability of "almost detected" rising toward the detection limit. That belongs in `mnarjoint` as an extra indicator channel (for example, the indicator block w·(M − mean) replaced by w·(1 − PEP)). It must be credited against mnarjoint's own null.

## Data provenance (rows for the king to merge into `data/PROVENANCE.md`)
| local file | URL | use / notes | date | SHA-256 |
|---|---|---|---|---|
| `data/cache/phasec/pepev/dl_leduc/leduc2022_pSCoPE.Rda` (re-used) | https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/scpdata/leduc2022_pSCoPE.Rda | Phase C `subthr`: per-batch PSM assays (eAL*/wAL*), sub-threshold PSMs 0.01 < dart_qval ≤ 0.1. Same file and terms as the existing leduc2022 row | 2026-10-23 | `805e1248e8a1ccdf9f3aa0be32a05bc837178a18826f56c36de32369cc4be932` |
| `data/cache/corpus_raw/sadiku2025/20250204_112949_gbm_sc_full_fasta_19.4_Report.tsv` (re-used) | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/11/PXD061065/20250204_112949_gbm_sc_full_fasta_19.4_Report.tsv | Phase C `subthr`: `.raw.PG.PEP (Run-Wise)` columns. CC0 (PRIDE API v3 PXD061065 licence field, checked 2026-10-23) | 2026-10-23 | `b8a1d6203c2d2feae5024ae3f6d35c950f498880496747f35ca89df348a48e8c` |
| `data/cache/phasec/subthr/dl/setup.txt` (new, 2,513,260 bytes) | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/11/PXD061065/20250204_112949_gbm_sc_full_fasta_19.4_Report.setup.txt | Phase C `subthr` audit: Spectronaut filter settings. CC0 (PRIDE PXD061065) | 2026-10-23 | `6e1dc617ecad0e2863a14e07ef0e5cc90343a29b0df57b70bc4c4f645fa5fc76` |
| `data/cache/corpus_raw/gregoire2023_mixCTRL/gregoire2023_mixCTRL.Rda` (re-used, audit only) | https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/scpdata/gregoire2023_mixCTRL.Rda | Phase C `subthr` audit: sage PSM assays (protein_fdr); not used in the model. Terms as the existing row | 2026-10-23 | `27f208c330527ec1ee42f953cfd6a7859a637f83aa4dda9ed083b8b81fe04224` |

Derived caches (rebuildable with `python -m scp_graphfs_mcp.phasec.subthr build`): `data/cache/phasec/subthr/{leduc2022,sadiku2025}_E.npz` and `leduc2022_psm_extract.npz`.
