# pepev: per-entry precursor/peptide-count precision weights in a heteroscedastic U-step of adapt_project

- **Queen:** new data modalities (round-2 directive)
- **Status (worker's recommendation):** **dead** at the fast screen, for both v1 and the one revision.
  - **v1** (`pepev`): cand_stack (rmse_all) is 0.9996 against null_stack 0.9998. It misses the 0.99 bar, and the
    margin over the null is 0.0002.
  - **v2** (`pepev_v2`, protein-level evidence weights): cand_stack is 0.9969 against null_stack 1.0003. The margin,
    0.0034, is below 0.005, and cand_stack misses 0.99.
  - The full pool was not run.
  - The abort rule did not trigger: 5 of 6 fast targets passed the build-time X-equality check, against a threshold
    of 4. The HPA-IF fallback was therefore not built.
- **Code:** `src/scp_graphfs_mcp/phasec/pepev.py`
  - `build` CLI: one builder per target.
  - `diag` CLI: diagnostics and the unit check.
  - Model functions: `adapt_project_w`, `impute` / `impute_null`, and the variants `impute_a1`, `impute_uv`,
    `impute_prot`, `impute_null_prot`, `impute_v2` / `impute_null_v2`.
- **Results** (all in `results/phaseC/`):
  - `pepev_fast.json` (the screen that counts)
  - `pepev_a1_fast.json`, `pepev_uv_fast.json`, `pepev_prot_fast.json`, `pepev_nullprot_fast.json` (preregistered
    variants)
  - `pepev_v2_fast.json` (the one revision)
  - `pepev_diag.json` (N statistics and the unit check)

## Idea
Every model so far weights each observed (cell, protein) value equally. A protein value quantified from many
precursors or peptides should be less noisy than a single-hit value. So the U-step of the adapted projection, which
gives each cell's factor scores, is turned into a weighted ridge with per-entry precision weights
w_ij ∝ N_ij^α, where N_ij is the precursor or peptide count. The principle is inverse-variance weighting, as in
limma-voom precision weights (Law et al. 2014, *Genome Biology*) and the peptide-count-dependent variance of DEqMS
(Zhu et al. 2020, *Mol Cell Proteomics*).

## Mechanism (as built)
### Step 1: count matrices aligned to `load_target` (cached `data/cache/phasec/pepev/<target>_N.npz`, with the sha of X)
Each builder replays the corpus builder on the protein-level source. It carries N through the same selection: the
`keep_cols` step, the dedup by most-observed accession, the ≥ 50-observed row filter and the per-cell median centring.
It then asserts `allclose(X_rebuilt, X_corpus, equal_nan=True)`, compared at float32 with atol 1e-5. Where X is
observed but N is 0 or missing, N is set to 1.

| target | source of N | X equality | observed entries with N > 0 in source |
|---|---|---|---|
| sadiku2025 | Spectronaut pivot, per-run `PG.NrOfStrippedSequencesUsedForQuantification` (cached file, no download) | pass (314 x 3039) | 1.0 |
| weiss2026_liver | DIA-NN `report.pr_matrix.tsv` (range-read from the PRIDE zip): non-NaN precursors per run and `Protein.Group` | pass (763 x 4520) | 1.0 |
| petrosius2023 | Spectronaut PEPQuant export: nunique `EG.PrecursorId` with `EG.Qvalue` ≤ 0.01 and `PEP.Quantity` > 0, per (`R.FileName`, `PG.ProteinAccessions`) | pass (215 x 2813); all 215 runs and 2813/2813 PGs matched | 1.0 |
| gregoire2023_mixCTRL | scpdata `peptides_<batch>_norm_log`, the direct parent assay of the `proteins_<batch>_norm_log` used by the corpus: non-NaN peptides per protein group (cached Rda, no download) | pass (119 x 1878) | 1.0 |
| leduc2022 | scpdata `peptides_log`, the parent of `proteins_norm2`: non-NaN peptides per `Leading.razor.protein.id` | pass (1543 x 2844) | 0.9841 |
| px063590_ADbrain | — | **not aligned (size)**: PSM tables 1.22 GB + 0.94 GB = 2.16 GB, over the cap | — |

Notes:
- **gregoire:** I used the batch-level `peptides_<batch>_norm_log` assays rather than the 63 per-run
  `peptide_<batch>_<run>` assays. The protein assay is aggregated from these, so they match the protein values
  exactly. Their rowData `proteins` column is the protein-group string that names the protein rows.
- **petrosius:** the X-equality check here is weaker than for the other targets. `load_target` reads the h5ad
  directly, so there is no corpus rebuild to replay. Alignment rests on the exact match of run names (`R.FileName` +
  `.raw`) and PG accessions, and on every observed entry having N > 0.

The fractions in the last column, and everything else in this doc, come from `pepev_diag.json` and the build log.

### Step 2: weights
- p_ij = N_ij^α, with α = 0.5 fixed a priori.
- Only entries o_i = obs(Xm)[i] & ok are used, where ok means the column is observed in Xm.
- Normalise per cell to mean 1 over o_i, clip to [0.25, 4], then renormalise.
- Weights are never read at entries that are NaN in the given Xm.

### Step 3: call-time alignment
For each call, the module recomputes `idx = explore_phaseb.subsample(n, labels, 30, seed)` and then asserts:
- `names == ctx['names']`;
- obs(Xm) ⊆ obs(X[idx]);
- `allclose(Xm[obs], X[idx][obs])`;
- Xm keeps at least 70% of X[idx]'s observed entries.

All assertions held on all 5 aligned targets × 2 masks × 2 seeds × {eval, inner}. On px063590_ADbrain, `ctx['champion']`
is returned unchanged, so its ratio is exactly 1.

### Step 4: model
`adapt_project_w` is a local copy of `transfer.adapt_project`. Only the U-step changes:
U[i] = solve((V_o·w)ᵀV_o + λI, (V_o·w)ᵀZ_io). The defaults are those of the harness call: lam_v = 100, lam = 5,
extra = 10, iters = 10, seed = 0.

- **Unit check:** with W = ones it matches `T.adapt_project` bit-exactly: `array_equal` True, max |diff| 0.0 on one
  MNAR seed-100 split of each of the 5 aligned targets (`pepev_diag.json`).
- **Variant `pepev_uv`:** the V-step is weighted as well.

## Null controls
- **Primary (`impute_null`):** identical pipeline, but each cell's weights are permuted among that cell's own observed
  entries (`default_rng(seed + 1000*i)`). This keeps the same weight distribution and mean of 1, and destroys the link
  between evidence and protein.
- **Second null (`impute_null_prot`, run as the null of `pepev_nullprot`):** N's protein columns are permuted within
  the dataset, so each protein takes another protein's count profile.
- **v2 null (`impute_null_v2`):** the protein-level weights (column means) are permuted across proteins.

## Results (fast pool: 6 targets × 2 masks × 2 seeds; mean ratio vs the champion, < 1 better)
- "All 6" is `mean_rel_to_champion` from the JSON. px063590 contributes exactly 1.0 to cand_stack and null_stack.
- "Aligned 5" is the mean over the 10 aligned (target, mask) records in `per_target_mask`.

| run | cand_raw | cand_solo | **cand_stack (all 6)** | null_stack (all 6) | cand_stack / null_stack (aligned 5) | in-complex cand / null | splits ≤ 0.99 | mean cand weight |
|---|---|---|---|---|---|---|---|---|
| **pepev** (α 0.5, U-step) | 0.9997 | 0.9988 | **0.9996** | 0.9998 | 0.9995 / 0.9997 | 0.9987 / 0.9996 | 1 / 12 | 0.5729 |
| pepev_a1 (α 1) | 1.0067 | 1.0041 | 1.0005 | 1.0001 | 1.0006 / 1.0002 | 0.9997 / 0.9999 | 1 / 12 | 0.4167 |
| pepev_uv (U + V weighted) | 1.0014 | 0.9997 | 0.9988 | 0.9998 | 0.9986 / 0.9998 | 0.9986 / 0.9996 | 1 / 12 | 0.4583 |
| pepev_prot (column-mean N; null = within-cell perm.) | 0.9966 | 0.9969 | 0.9969 | 0.9998 | 0.9963 / 0.9997 | 0.9967 / 0.9996 | 1 / 12 | 0.6875 |
| pepev_nullprot (cand = pepev; null = protein-permuted N) | 0.9997 | 0.9988 | 0.9996 | 1.0001 | 0.9995 / 1.0001 | 0.9987 / 1.0001 | 1 / 12 | 0.5729 |
| **pepev_v2** (column-mean N; null = column means permuted) | 0.9966 | 0.9969 | **0.9969** | 1.0003 | 0.9963 / 1.0004 | 0.9967 / 1.0002 | 1 / 12 | 0.6875 |

Runtimes were 140.2 s (pepev), 311.1, 324.9, 326.7, 335.1 and 273.8 s. The later runs shared the box with
another worker.

**Ordering check (pepev, rmse_all, all 6):**
- real 0.9996 < within-cell permuted 0.9998 < protein-permuted 1.0001; all-ones = champion = 1.000.
- The ordering holds, but the whole spread is 0.0005.

**Per target, main run `pepev` (rmse_all; cand_raw / cand_stack / null_stack):**

| target | MCAR | MNAR |
|---|---|---|
| gregoire2023_mixCTRL | 1.0085 / 1.0002 / 0.9999 | 1.0032 / 1.0 / 1.0 |
| leduc2022 | 1.0044 / 0.999 / 1.0 | 0.9996 / 0.9976 / 1.0 |
| petrosius2023 | 0.9956 / 0.9956 / 0.9993 | 1.0097 / 1.0097 / 1.0 |
| px063590_ADbrain (not aligned) | 1.003 / 1.0 / 1.0 | 1.0008 / 1.0 / 1.0 |
| sadiku2025 | 0.9934 / 1.0 / 1.0007 | 0.9907 / 1.0052 / 0.9985 |
| weiss2026_liver | 0.9759 / 0.9759 / 0.9974 | 1.0121 / 1.0121 / 1.0014 |

Over the 5 aligned targets, pepev's cand_stack is 0.9941 on MCAR and 1.0049 on MNAR. The null is 0.9995 and 1.0.

**Per target, revision `pepev_v2` (rmse_all; cand_raw / cand_stack / null_stack):**

| target | MCAR | MNAR |
|---|---|---|
| gregoire2023_mixCTRL | 1.0072 / 1.0032 / 1.0002 | 1.0075 / 0.9999 / 1.0003 |
| leduc2022 | 1.0056 / 0.999 / 1.0 | 1.0017 / 0.9985 / 1.0001 |
| petrosius2023 | 0.9961 / 0.9967 / 1.0003 | 0.9993 / 0.9993 / 1.0 |
| sadiku2025 | 0.9939 / 0.9946 / 1.001 | 0.9676 / 0.9959 / 1.003 |
| weiss2026_liver | 0.9806 / 0.9806 / 0.9977 | 0.9955 / 0.9955 / 1.001 |

## Diagnostics (`pepev_diag.json`)
| target | N p10/p25/p50/p75/p90 | frac N = 1 | Spearman(N, value): pooled / within cell (median) / within protein (median) | share of var(log N) between proteins |
|---|---|---|---|---|
| sadiku2025 (DIA, Astral) | 1 / 1 / 2 / 4 / 8 | 0.394 | 0.401 / 0.402 / 0.207 | 0.806 |
| weiss2026_liver (DIA, Astral) | 1 / 1 / 2 / 4 / 8 | 0.391 | 0.632 / 0.638 / 0.442 | 0.776 |
| petrosius2023 (DIA, Astral) | 1 / 1 / 2 / 3 / 6 | 0.423 | 0.426 / 0.452 / 0.354 | 0.790 |
| gregoire2023_mixCTRL (TMT) | 1 / 1 / 2 / 3 / 6 | 0.491 | 0.003 / −0.010 / −0.032 | 0.790 |
| leduc2022 (TMT) | 1 / 1 / 2 / 3 / 5 | 0.483 | 0.005 / 0.002 / 0.004 | 0.785 |

- On the MNAR seed-100 split, the median N of the masked entries equals that of the kept entries (2 vs 2) on every
  target.
- Per-split weights after clipping range from about 0.57 to 4.0 (p10 ≈ 0.58–0.66, p50 ≈ 0.84–0.91, p90 ≈ 1.58–1.70).
  The lower clip at 0.25 never binds.

## Honest read
1. **The data side worked.** Peptide/precursor evidence was aligned exactly to 5 of the 6 fast targets, about 620 MB
   in all. The weighted U-step reproduces the champion bit-exactly at W = 1, and every call-time alignment assertion
   held. These count matrices now exist for abd's queued extreme-deconvolution analogy.
2. **Per-entry evidence weighting (the hypothesis) does not help on average.** The result is 0.9996 against a null of
   0.9998.
   - It helps on MCAR (aligned MCAR cand_stack 0.9941; weiss2026_liver 0.9759).
   - It hurts on MNAR (aligned 1.0049; weiss 1.0121, petrosius 1.0097).
   - α = 1 amplifies both effects: weiss MCAR 0.9685, MNAR 1.0279.
   - The mechanism is visible in the diagnostics. On the DIA targets, N tracks the entry's *own* value, even within a
     protein (within-protein Spearman 0.21–0.44). Weighting by N therefore up-weights the high-value entries, which
     makes the weights informative about the outcome. The fitted cell scores then over-explain each cell's high
     values and extrapolate worse to the low-abundance entries that an MNAR mask removes.
   - On the TMT targets, N carries no value information (Spearman ≈ 0) and the effect is ≈ 0.
3. **The usable part is per-protein, not per-entry.**
   - About 78–81% of the variance of log N lies between proteins.
   - Protein-only weights (column-mean N from the dataset's other cells, so the split's own cells are excluded) give
     0.9969 and do not degrade MNAR (aligned MNAR 0.9978 in `pepev_prot`).
   - That is why the one revision (v2) used them, with the proper protein-permutation null. v2 gives 0.9969 vs
     1.0003, a margin of 0.0034 (aligned 5: 0.9963 vs 1.0004, margin 0.0041).
   - This is the cleanest separation from its null in round 2, but it is still below both bars. It also resembles
     round 1's per-protein results (ebshrink, scrnaprior): per-protein reliability is a ~0.3% effect.
4. The uv variant (0.9988) shows that weighting the V-step too softens the MNAR damage. It is still far from 0.99.
5. The stack is part of the problem on MNAR.
   - On sadiku2025 MNAR, v2's raw candidate is 0.9676 but the stack picks 0.9959.
   - The stack weights are chosen on a random inner mask, so a candidate that wins specifically under MNAR gets
     diluted.
   - This is a harness property, not something tuned here. It may matter for other ideas.

## Next step (if anyone continues this line)
- **Decorrelate the per-entry weight from the value before weighting.** For example, use the residual of log N_ij
  after regressing it on the protein's mean log N and the entry's own value. This is a different estimator and needs
  its own preregistration, so it was not done here. Expected size: small, given that the per-protein part is ~80% of
  the signal.
- **Feed N into abd's extreme-deconvolution idea** as per-point noise variances (σ²_ij ∝ 1/N_ij) inside a likelihood
  that models MNAR censoring explicitly. The censoring model is what keeps informative weights from biasing U.
- **Do not spend a full-pool run on pepev or pepev_v2.**

## Data provenance (rows for the king to merge into `data/PROVENANCE.md`)
| local file | source URL | licence | date | SHA-256 | notes |
|---|---|---|---|---|---|
| `data/cache/phasec/pepev/dl_petro/orig.tsv` | https://dataverse.uclouvain.be/api/access/datafile/24335?format=original (`20240201_130747_PEPQuant (Normal).tsv`, dataset doi:10.14428/DVN/4DSPJM) | CC BY-SA 4.0 (Dataverse API `license.uri` = https://creativecommons.org/licenses/by-sa/4.0/) | 2026-10-23 | `f629ac7a6380c25db8328e564b953f76a9c7b25a7ffbb71d74eadad6a662b0d7` | 294,665,776 B, MD5 `587e02fbdbee52857b33dc44d46b04a5` = Dataverse's recorded MD5. That MD5 belongs to the *original* upload. The default endpoint serves the ingested `.tab` (310,158,870 B with a variable-name header, 310,158,660 B without), which has a different MD5 (`dc391db8bfa692815306ac08a9451789` for the headerless copy). Both `.tab` copies were downloaded during checking and then deleted. |
| `data/cache/phasec/pepev/dl_weiss/report.pr_matrix.tsv` | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/12/PXD062231/search_engine_results.zip, member `search_engine_results/report.pr_matrix.tsv`, HTTP range reads only | "Creative Commons Public Domain (CC0)" (PRIDE API v3 project record PXD062231, `license` field, verbatim) | 2026-10-23 | `997b1abf9b4fa17780301fb86415fad8bccc89f1fa4cc32b93813884ba79337c` | 57,535,592 B raw (23,286,660 B compressed); zip CRC-32 `0xe0e4598a`, verified by `zipfile` on extraction |
| `data/cache/phasec/pepev/dl_leduc/leduc2022_pSCoPE.Rda` | https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/scpdata/leduc2022_pSCoPE.Rda | as the existing leduc2022 PROVENANCE row (same file) | 2026-10-23 | `805e1248e8a1ccdf9f3aa0be32a05bc837178a18826f56c36de32369cc4be932` | 268,763,036 B; same hash as the existing PROVENANCE row. Used for `peptides_log` / `proteins_norm2`; the extract is cached as `leduc2022_extract.npz` |
| (no download) sadiku2025 | existing `data/cache/corpus_raw/sadiku2025/20250204_112949_gbm_sc_full_fasta_19.4_Report.tsv` | existing row | — | existing row (`b8a1d620…`) | per-run `PG.NrOfStrippedSequencesUsedForQuantification` columns |
| (no download) gregoire2023_mixCTRL | existing `data/cache/corpus_raw/gregoire2023_mixCTRL/gregoire2023_mixCTRL.Rda` | existing row | — | existing row (`27f208c3…`) | `peptides_<batch>_norm_log` assays + rowRanges `proteins` |

The derived count matrices are cached in `data/cache/phasec/pepev/<target>_N.npz`. Each holds N (float32), the
SHA-256 of the target's X and the protein names; the impute functions refuse a stale cache. Total on disk is
644 MB, which is under the 2 GB cap.
