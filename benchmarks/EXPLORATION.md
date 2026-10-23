# Exploration protocol (phase 2): develop on dev masks, confirm once on held-out masks

The preregistered phase-1 result (`benchmarks/PREREGISTRATION.md`, `results/`) **stands as reported**:
the curated priors beat their shuffled nulls, but they lose to kNN. This file fixes, *before any phase-2 run*,
how we look for a method that might make the claim true without fooling ourselves.

## Why a protocol
If we try ideas until one beats kNN on the same masks, we overfit the benchmark: a forking-paths
problem. So we split the masks:

| Use | Mask seeds | Masks | Who sees the numbers |
|---|---|---|---|
| Phase 1 (done) | 0–4 | MCAR, MNAR | already public |
| **Development** | **100, 101, 102** | MCAR, MNAR | used freely to try ideas, logged in `docs/EXPLORATION_LOG.md` |
| **Confirmation** | **1000–1004** | MCAR, MNAR | used **once**, after `benchmarks/PREREGISTRATION_v2.md` is committed |

The mask generator, eligible proteins (detected in ≥ 50% of cells), mask fraction (10%) and strata are
identical to phase 1.

## Fair "data alone" baselines (stronger than phase 1)
The claim is "curated complexes are a better prior than data alone". Every graph-based idea is
therefore run on four kinds of graph:
1. **CORUM** and **Complex Portal** (curated);
2. **their degree-preserving nulls** (seeds 0–9);
3. **a data-derived co-variation graph:** for each protein, its most-correlated proteins in the
   *masked training matrix* (pairwise-complete Pearson), with the same total edge count as Complex
   Portal. This is the strongest "data alone" neighbour set;
4. **STRING** (combined and experimental).

On top of these, kNN (k = 10) stays a baseline. A phase-2 method supports the claim only if, on the
confirmation masks and in the in-complex stratum, the curated version beats **kNN, its own null, and
the data-derived graph**, with paired 95% bootstrap CIs excluding 0 on both masks.

## Ideas to try (cross-domain analogies)
| ID | Idea | Analogy |
|---|---|---|
| E1 | **Residual-on-kNN:** kNN gives each cell's state; partners explain what kNN misses. x̂ = kNN + β · (edge-weighted mean of partners' kNN residuals in the same cell) | Gradient boosting / residual learning; "baseline + correction" in sensor fusion |
| E2 | **Stoichiometric offsets:** slope-1 log-ratio model per edge, x̂_j = mean_k (x_k + δ_jk), where δ_jk is the training median of x_j − x_k | Mass balance / conserved ratios in chemistry; fixed isotope ratios |
| E3 | **Graph-regularised low-rank completion:** X ≈ U Vᵀ with a Laplacian penalty tr(Vᵀ L V) over the protein graph | Netflix-style recommender systems with a social/item graph (GRMF) |
| E4 | **Edge reliability gating:** reweight each prior edge by its training co-variation, shrunk toward the prior weight; the prior restricts *which* pairs may be used | Attention restricted by a mask; trust calibration; Kalman gain by measured precision |
| E5 | **Inverse-variance fusion of kNN and partner estimates,** per protein, from training residual variances | Meta-analysis / sensor fusion (inverse-variance weighting) |
| W1 | **Contaminant filter** (keratins etc.) before panel selection | Data hygiene (cRAP lists in proteomics) |
| P1 | **Panel as sensor placement:** greedy submodular coverage of complexes vs. data-only greedy column subset selection | Sensor placement / facility location; pivoted-QR column selection |

## Rules
- Every idea's dev result is logged, including failures, in `docs/EXPLORATION_LOG.md` and on its GitHub issue.
- Hyperparameters are tuned on dev masks only.
- **Selection rule (fixed now):** choose the *single* imputation method with the lowest dev in-complex RMSE
  (averaged over both masks), among methods whose curated version beats kNN, its null and the data graph
  on dev. If none qualifies, we report that and do not run confirmation as a "win".
- Before the confirmation run, `benchmarks/PREREGISTRATION_v2.md` is committed with the chosen method,
  hyperparameters, baselines and pass criteria. The confirmation masks are run once.
- Phase-2 outputs go to `results/phase2/` (dev runs to `results/phase2/dev/`).

## Amendment 1 (committed before the runs it covers)

Dev seed 100, MCAR, was already run (in `results/phase2/dev/dev_results_quick.json`). It showed:
- residual-on-kNN with curated graphs beats kNN and the null, but the data-derived graph and STRING beat it;
- low-rank completion (E3) beats everything, and graph regularisation adds nothing over its null.

Before running more, three dev dimensions are added:
1. **Low-data regime (prior-matters-most hypothesis; analogies: recommender cold-start, Bayesian shrinkage).**
   - Subsample n_cells ∈ {30, 60} (cell-subsample seeds = mask seed; stratified by cell type), then mask
     within the subsample.
   - The data-derived graph and kNN use only those cells. The curated prior is unchanged.
2. **Prior + data (does the prior add information beyond the data?).** Graphs: union(data-derived,
   curated) vs. union(data-derived, null) vs. data-derived alone.
3. **E6 = residual correction on a low-rank base:** x̂ = GRMF + β · partner-weighted GRMF residuals.

**Scoping rule:** if the only qualifying method wins in a regime (e.g. n_cells ≤ 60), the v2
preregistration states that regime, and the README reports the claim *as scoped*. The phase-1
full-data result stays as reported.

## Amendment 2: an independent held-out dataset (committed before any evaluation on it)

Repeated seeds on the same 215 cells can't fully remove forking-path risk, so we add an **independent
dataset**: Leduc et al. 2022 (nPOP/pSCoPE), melanoma + monocytes, the **pre-imputation** protein assay
`proteins_norm2`. That's 1543 cells and 65.9% missing; provenance is in `data/PROVENANCE.md`.

- It's built from the scpdata `.Rda` without R (`scripts/rda_extract_leduc.py`, `scripts/build_leduc2022.py`).
- **No method has been evaluated on it.** No imputation, panel or correlation analysis may be run on it
  until `benchmarks/PREREGISTRATION_v3.md` is committed. v3 fixes the method(s) using Petrosius
  evidence only (dev + v2 confirmation).
- Allowed before v3: shape, missingness and label counts (already printed by the build script).

## Amendment 3: a third held-out dataset (committed before any evaluation on it)

**Dataset:** SCoPE2 (Specht et al. 2021) from scpdata `specht2019v3.Rda`, 1490 cells, 79.5% missing. The scpdata
`proteins` assay is imputed, so proteins are re-aggregated from the **unimputed** `peptides` assay (median per
protein per cell). Provenance is in `data/PROVENANCE.md`.

**What has been looked at:** only shape, missingness and label counts.

**Purpose:** it is reserved for v5. v5 tests the one claim that held up across Petrosius and Leduc: curated
partners beat data-derived partners. This time the test includes a minimum effect size.

## Amendment 4: "more data as inductive bias", with holdouts reserved before download

**Idea (analogy: pretraining / transfer learning):** learn the partner prior from *other* datasets instead of
from the dataset being imputed. The sources are curated complexes, predicted complexes (hu.MAP 3.0), and protein
co-variation pooled across an external corpus of single-cell proteomics datasets. Each external dataset is a
noisy, independent measurement of which proteins move together.

**Reserved holdouts (declared before they were downloaded; never used to build priors or tune anything):**
- **khan2023** (scpdata EH9477; 421 MCF-10A cells, EMT): the v6 confirmation dataset.
- **derks2022** (scpdata EH7712; plexDIA): a spare for a later confirmation.

**Pretraining corpus, which is also the dev set:**
- petrosius2023_AstralAML, leduc2022_pSCoPE and specht2019v3 (already local);
- schoof2021 and woo2022_macrophage (to be added).

Dev evaluation is **leave-one-dataset-out:** a target dataset's prior is built only from the *other* corpus
datasets. The dev mask seeds are the phase-2 ones (100–102); subsample seeds for low-n are the same.

**Required baselines:** kNN, low-rank (`grmf`, λ_g = 0), the in-dataset data graph (matched edges), degree-preserving
nulls, and **an external data-only graph** (pooled co-variation, no curation). The last one is the
"more data alone" control, so a win must come from curation × external evidence, not from external data alone.

**v6 rule (fixed now):** pick the single method with the best mean dev in-complex RMSE across corpus targets
whose curated version beats all baselines on dev. Preregister it with a **SESOI of 1% RMSE** against every
comparator, then run once on khan2023.

## Amendment 5: knowledge as protein *representations* (side-information matrix completion); committed before any run

**Why a new representation:** edge-based partner corrections are ≈ 0 after kNN or low-rank on Leduc and SCoPE2 (Amendment-4 dev).
Knowledge may enter more usefully as **protein feature vectors** that shape a low-rank model's protein loadings.
- **Analogies:** recommender systems with item side information (inductive / collective matrix factorisation), and
  pretrained embeddings.
- **Where it should matter most:** at **few cells**, where plain low-rank collapses (n = 28: 1.05 vs. kNN 0.92).

**Model (IMC):** z-score X per protein, then fit X ≈ U Vᵀ on observed entries.
- **Loadings:** V = F W + Δ, where F is a p×d protein feature matrix from knowledge sources.
- **Penalties:** W is penalised lightly; Δ (protein-specific deviation) is penalised heavily, so the prior dominates when data are scarce.
- **Variant "pretrained V":** V₀ is taken from a factor model fitted on the external corpus (leave-target-out), and v_j = V₀_j + δ_j.

**Feature sources (each logged in `data/PROVENANCE.md`):**
- GO annotations (goa_human GAF);
- Human Protein Atlas (subcellular location, tissue/cell specificity);
- UniProt (length, subcellular location);
- PaxDb integrated abundance;
- spectral embeddings of the STRING (combined, physical), Complex Portal and hu.MAP graphs;
- external-corpus loadings.

Second tier, only if time allows: Pharos (target development level), Monarch (phenotypes), INDRA (mechanisms).

**Required baselines (all must be beaten):**
1. kNN;
2. plain low-rank (no features);
3. **the same IMC with feature rows permuted across proteins** (null F, 5 seeds), the analogue of a degree-preserving shuffle;
4. for the pretrained variant, the corpus factor model *without* curated features.

**Dev:** leave-one-dataset-out on the corpus targets. Petrosius uses 215 and 28 cells; Leduc and SCoPE2 use 200- and 30-cell subsamples. Dev seeds are 100–102.

**Selection rule (fixed now):** among IMC variants whose curated-feature version beats all four baselines on dev on
≥ 2 of 3 targets in the 30-cell regime, pick the one with the lowest mean 30-cell RMSE. It is then preregistered
(v6) with a **1% SESOI** against each baseline and run once on **khan2023**. If none qualifies, that is reported
and khan2023 stays untouched.

## Amendment 6: ensemble of kNN and knowledge-IMC (committed before running it)

**Motivation:** Amendment-5 dev found that knowledge features beat shuffled features in every 30-cell case, but kNN beats
the IMC family there. kNN errors come from neighbouring cells, and IMC errors come from protein knowledge, so the two errors may be
partly independent. An ensemble may beat both (analogy: Netflix-prize blending, sensor fusion).

**Method E7:** x̂ = (1 − w) · kNN + w · IMC(curated_all), with **w = 0.5 fixed** (nothing is tuned).

**Required baselines:**
- kNN;
- **the null ensemble** (1 − w) · kNN + w · IMC(shuffled curated_all features, seeds 0–1), so a gain must come from
  the knowledge, not from averaging two models;
- the kNN + IMC(no features) ensemble.

**Dev:** the same targets, regimes and seeds as Amendment 5, in the 30-cell regime plus all-cells / 200 cells.

**Rule:** E7 qualifies if it beats all three baselines on ≥ 2 of 3 targets at 30 cells, on both masks.
If it qualifies, it is preregistered as v6 (SESOI 1% vs. each baseline) and run once on **khan2023** at 30-cell
subsamples (seeds 6000–6009).

## Amendment 7: protein-specific kNN, with knowledge-defined cell similarity (committed before running)

**Idea (analogies: attention, locally weighted regression):** to impute protein j, measure cell–cell distance on j's
knowledge-defined neighbourhood S_j instead of on all proteins. S_j = the m = 50 proteins most similar to j by cosine
similarity of the curated_all feature vectors (Amendment 5).
- **Distance:** d_j = 0.5 · d_local(S_j) + 0.5 · d_global (both nan-euclidean, rescaled to mean 1). This is fixed and not tuned.
- **Estimate:** the mean of the k = 10 nearest donors that observe j.

**Required baselines:**
1. standard kNN;
2. **null:** S_j from row-shuffled features (seeds 0–1);
3. **data-only:** S_j = the 50 proteins most correlated with j in the masked training matrix.

**Rule:** it qualifies if it beats all three on ≥ 2 of 3 targets at 30 cells, on both masks. If so: v6 on khan2023 with a 1% SESOI.

## Amendment 8: tune pkNN's local/global mix on dev (committed before running)

**Basis:** Amendment-7 dev with α = 0.5 was within ~1% of kNN on Leduc and SCoPE2. Tuning on dev is allowed by this protocol;
the confirmation set (khan2023) stays untouched.

**Sweep:**
- α ∈ {0.1, 0.25} (α = 0.5 was already run) with m = 50 and k = 10;
- the same baselines (kNN; null neighbourhoods seeds 0–1, at the same α; data neighbourhoods, at the same α);
- the same 30-cell dev configurations.

**Selection:** the α whose curated pkNN beats all baselines on ≥ 2 of 3 targets, on both masks, with the lowest mean RMSE.
If one qualifies, v6 is preregistered on khan2023 (30-cell subsamples, seeds 6000–6009, 1% SESOI vs. each baseline).
If none qualifies, the sweep stops here: there will be no further α values and no other k or m.

## Amendment 9: exhaustive representation sweep (committed before any run it covers)

**Why:** Amendments 4–8 tested five knowledge sources and three estimators. Amit asked (Oct 22) for an exhaustive attempt at
"every single representation" before we accept a negative result. This amendment fixes the whole grid, the baselines,
the selection rule and the single confirmation **before** any real-data run. Representation modules are built by
agents that are not allowed to run imputation or see results.

**Representations (protein side).** Each is a graph (neighbours with weights), a feature block (cosine neighbours), or a
bulk reference panel (correlations and ridge coefficients learned *in the reference*):
- **Already in the repo:** Complex Portal, CORUM, STRING, hu.MAP 3.0 (graphs); GO, HPA, UniProt keywords, PaxDb/physical,
  external single-cell corpus loadings (features).
- **New:**
  - Sequence and structure: ProtT5 sequence embeddings (UniProt); PDB co-structure (SIFTS: chains in the same entry).
  - Drug-target annotations: Pharos/TCRD.
  - Pathways and interaction databases: Reactome pathways; BioGRID physical; IntAct.
  - Phenotype and literature: Monarch phenotypes (and diseases); NCBI gene2pubmed co-citation.
  - Mechanistic graph: INDRA, or Pathway Commons if INDRA is not openly downloadable.
  - Bulk co-regulation panels: CCLE bulk proteomics (375 lines) and ProteomeHD.
- **Combined:**
  - `union_curated`: the union of all knowledge graphs (excluding bulk panels), weight = number of sources.
  - `concat_features`: all feature blocks concatenated.
- A source that fails to download or covers < 20% of a target's proteins is reported as such and dropped. Dropping
  never depends on imputation results.

**Estimators** (each run on every representation, with fixed hyper-parameters):
1. **RPR** (reference-panel regression, `scp_graphfs_mcp.refpanel`). It uses m = 20 neighbours and λ = 10. Bulk panels use
   reference ridge coefficients (λ_ref = 1, `signed=True`); graphs and features use similarity-weighted averages.
2. **RPR-stack:** `stack_with_base` mixes RPR with kNN. The weight w ∈ {0, .25, .5, .75, 1} is chosen on a 10% inner
   mask of *observed* eligible entries.
3. **pkNN** (Amendment 7, α = 0.5, m = 50, k = 10). Graphs enter through their spectral embedding; features directly.
   It runs at 30 cells only, as before.

**Baselines and nulls** (all must be beaten):
- kNN; low-rank (`grmf(..., lam_g=0)`);
- the generic ensemble `stack_with_base(kNN, low-rank)`;
- the data-only control: RPR and RPR-stack using neighbours from the target data's own correlations (`neighbours_from_data`);
- a **null per representation:** RPR and RPR-stack with the neighbour lists relabelled at random (`permute_neighbours`,
  seeds 0 and 1); pkNN with `permuted_features` (seeds 0 and 1).

**Grid (dev only):**
- targets: Petrosius 2023, Leduc 2022, SCoPE2;
- MCAR and MNAR masks;
- seeds 100–102;
- regimes: 30 cells (**primary**) and 200 cells (all 215 for Petrosius).

The data representation stays log2 as loaded. Cell-median centring is applied only in the robustness check below.

**Selection rule:**
- A (representation, estimator) pair *qualifies* if, at 30 cells, its mean RMSE (over seeds) is ≥ 1% lower than **every**
  baseline and its own nulls on ≥ 2 of 3 targets, under both masks.
- Among qualifiers, the one with the lowest mean RMSE relative to kNN is selected. Only that pair goes to confirmation,
  so there is a single confirmatory test despite the large dev grid.
- Every cell of the grid is published (`results/phase6/dev/`) whatever the outcome.

**Robustness check (selected pair only, dev):** repeat with cell-median-centred data (`--normalise median`) and at
200 cells. This is reported but is not a gate.

**Confirmation (v6, only if something qualifies):**
- The holdout is khan2023 (reserved in Amendment 4; never inspected beyond shape and missingness).
- It is run on 30-cell subsamples, seeds 6000–6009, under MCAR and MNAR.
- Success: the selected pair beats each baseline and its nulls by ≥ 1% mean RMSE (SESOI). The 95% paired bootstrap CI of
  the relative difference must exclude 0 in favour of the method, under both masks.
- `PREREGISTRATION_v6.md` will be committed before the holdout is loaded for evaluation.

**If nothing qualifies:** we report the full grid as an exhaustive negative result. There is no second sweep.

**Amendment 9: deviation 1 (implementation bug, fixed before any selection was acted on).** In the first full run,
`refpanel.rpr` returned NaN for any protein whose neighbour list contained a protein never observed in that subsample. The
kNN fill leaves such proteins NaN. As a result, 3,948 of 6,138 RPR / RPR-stack rows were NaN. The data-only control never
selects empty columns, so it was unaffected, as was pkNN. Fix: those proteins contribute z = 0 (test
`test_rpr_neighbour_never_observed_does_not_poison`). The whole grid is re-run unchanged. The first run's files are kept as
`results/phase6/dev/run1_buggy/` for transparency. The selection rule, grid and hyper-parameters are unchanged.

## Amendment 10: transfer from external references for tiny targets (committed before any run it covers)

**Why (dev-derived motivation, disclosed):**
- In Amendment 9's first run, unaffected parts only, low-rank completion was ~39% *worse* than kNN at 30 cells but
  ~2% *better* at 200, and the kNN + low-rank stack was the strongest baseline.
- Our reading: at 30 cells the bottleneck is learning *which proteins move together*, not knowing the prior.
- Fields with tiny target samples import that structure from elsewhere:
  - ancestry PCA projects new genomes onto reference PCs;
  - recommender systems handle cold-start items with side information;
  - genotype imputation copies from reference haplotypes;
  - NLP uses a masked model pretrained on a large corpus. Our version is linear, since no deep-learning stack is
    installed.

Amit asked (Oct 22) for "a wildly new representation or strategy" drawn from other domains.

**References** (`scp_graphfs_mcp.transfer`; the target is always excluded and holdouts are never included):
- `corpus`: the other single-cell datasets in `external.CORPUS`, z-scored per protein within each dataset;
- `all`: `corpus` plus CCLE (375 lines) and ProteomeHD (294 experiments) bulk samples.

**Strategies** (hyper-parameters fixed: rank 20, λ = 5, 15 ALS iterations, ≥ 20 reference observations to count a
protein as warm; donor k = 10 from 2000 reference cells, ≥ 20 co-observed proteins):
1. `refmf_<ref>`: protein loadings learned on the reference by masked ALS. Each target cell fits its own scores on its
   observed proteins (ridge), and the missing entries are predicted from them. Proteins the reference never saw keep
   the kNN fill.
2. `refmf_<ref>_kg`: as 1, but cold proteins take the weighted mean loading of their *warm* knowledge neighbours. The
   neighbours come from the Amendment-9 `union_curated` graph (top 20). This is where the curated prior is tested
   directly.
3. `donor_corpus`: kNN whose donors are external corpus cells, mapped to target scale with the target's per-protein
   mean and SD.
4. `stack(kNN, s)` for each s in 1–3: the same inner-mask weight rule as Amendment 9.

**Nulls and baselines:**
- Baselines: kNN, low-rank, and `stack(kNN, low-rank)`.
- Null for `refmf`: loadings permuted among warm proteins (`permute_rows`, seeds 0 and 1).
- Null for `_kg`: the same strategy with the knowledge neighbour lists relabelled (`refpanel.permute_neighbours`,
  seed 0).
- Null for `donor`: reference columns permuted among proteins (seed 0).
- **Secondary (thesis-specific, not a gate):** `refmf_<ref>_kg` vs `refmf_<ref>`.

**Grid:** as in Amendment 9 (three dev targets; MCAR and MNAR; seeds 100–102; 30 cells primary; 200 cells / all of
Petrosius secondary).

**Selection:** the Amendment-9 rule. A pair qualifies if it is ≥ 1% better than every baseline and its own null(s) on
≥ 2/3 targets under both masks at 30 cells. The single best qualifier is confirmed.

**Holdout:**
- If Amendment 9 has no qualifier, Amendment 10's single confirmation uses khan2023, as specified there: seeds
  6000–6009, 30-cell subsamples, ≥ 1% SESOI, and the paired-bootstrap CI must exclude 0 under both masks.
- If Amendment 9 has a qualifier, khan2023 goes to it, and Amendment 10 waits for derks2022 to be built under the same
  rule.
- Every result is published either way.

## Amendment 11: a bigger reference, plus a three-way stack (committed before any run it covers)

**Why (from Amendment 10 dev):**
- Transfer beats its permuted null by 3–12% at 30 cells, and a larger reference (`all` > `corpus`) was better.
- The kNN + low-rank and kNN + transfer stacks win on different targets.
- Amit asked for more data as the inductive bias.

**Bigger corpus:**
- Candidate scpdata datasets, fixed now: williams2020_lfq, williams2020_tmt, cong2020AC, liang2020_hela, brunner2022,
  woo2022_lung, gregoire2023_mixCTRL, guise2024, krull2024, hu2023_K562, bubis2025, leduc2025, wu2026, ai2025a.
- Each is added as single cells only. The rule is human only (no ortholog mapping). The protein assay must be
  unimputed; if only an imputed protein assay exists, proteins are rebuilt as the median of unimputed peptides, and if
  that is not possible the dataset is dropped. Values are log2 and per-cell median-centred, as in `scripts/build_corpus.py`.
- **Excluded to protect the holdouts:** khan2023 and derks2022 themselves, plus leduc2022_plexDIA, which may share
  plexDIA runs with derks2022.
- **Also excluded:** specht2019v2, a subset of the SCoPE2 dev target.
- Inclusion depends only on whether the build succeeds, never on imputation results. The final list is recorded in
  `data/PROVENANCE.md` before the run.

**Methods:**
- `refmf_big` uses the Amendment-10 projection with the reference = the big corpus (leave target out) + CCLE + ProteomeHD,
  with unchanged hyper-parameters.
- `stack2_big` is the kNN + `refmf_big` stack, as in Amendment 10.
- `stack3_big`: w_knn·kNN + w_lr·low-rank + w_t·`refmf_big`. The weights are on the simplex grid with step 0.25 and are
  chosen on the same 10% inner mask of observed entries.

**Nulls and baselines:**
- Baselines: kNN, low-rank, the kNN + low-rank stack, and `stack2_all` from Amendment 10 (the smaller reference).
- Null: `refmf_big` with loadings permuted among warm proteins (seeds 0 and 1), passed through the same stack.

**Grid and selection:** the Amendment-10 grid and rule. Each candidate must be ≥ 1% better than every baseline and its
null on ≥ 2/3 targets under both masks at 30 cells.

**Holdout:** khan2023 (Amendment 9 had no qualifier, and Amendment 10 had none). The single confirmation goes to the best
Amendment-11 qualifier, using seeds 6000–6009, 30-cell subsamples, a 1% SESOI, and a paired-bootstrap CI that excludes 0
under both masks. The reference for khan2023 is the full big corpus + bulk panels; khan2023 is never in it.

## Amendment 12: spend the khan2023 holdout on the best dev candidate (post-gate confirmatory test; committed before the holdout is evaluated)

**Status:** no candidate in Amendments 9–11 passed the dev gate.

**Decision:**
- Amit has asked us to push for a positive result while staying honest. We therefore spend the single khan2023
  confirmation on the strongest dev candidate, **`stack3_big`** (kNN + low-rank + reference projection with the
  Amendment-11 reference).
- It is labelled a **post-gate** test: it did not pass our own dev qualification rule, and any write-up must say so.
- It is the only use of khan2023. If it fails, the project reports a negative result. derks2022, if ever rebuilt from
  unimputed precursors, remains for future work.

**Hypothesis (H12):** on 30-cell subsamples of khan2023, `stack3_big` has lower overall RMSE than:
- kNN;
- low-rank;
- the kNN + low-rank stack;
- its own permuted-loading null (mean of seeds 0 and 1).

The margin must be ≥ 1% (SESOI), and the 95% paired-bootstrap CI must lie entirely below 0. This must hold under **both**
MCAR and MNAR. `stack2_all` is reported as a **secondary** comparator (whether the bigger reference adds over the smaller).

**Details:** `benchmarks/PREREGISTRATION_v6.md` and `benchmarks/confirm_v6.py`, run once, at the commit recorded in the output.

## Amendment 13: a biology-integrated model (BIM) with a slightly refined hypothesis (committed before any run it covers)

**Why:**
- Amendments 9–12 show that curated *physical* complex membership is the best protein prior (A9 ranking: Complex Portal,
  PDB co-structure, CORUM).
- They also show that external references add real signal, but the prior was never *integrated* into the transferred
  structure. It was only tried as a separate imputer or as a cold-start patch.
- Amit (Oct 22): "integrate as much inductive bias from biological priors as possible; tweak the hypothesis slightly if
  needed."

**Refined hypothesis.** The thesis is about complexes, so the endpoint is RMSE on **complex-member proteins**
(`in_any_complex`, the phase-1 in-complex endpoint), which is where a complex prior can act. Overall RMSE is secondary.
- **H13a (useful):** BIM beats kNN and the kNN + low-rank stack by ≥ 1% in-complex RMSE.
- **H13b (biology-specific):** BIM beats (i) the same pipeline without the biological prior (`stack3_plain`) and
  (ii) the same pipeline with the prior's protein labels permuted (`bim_null`), by ≥ 0.25% in-complex RMSE.

**BIM (fixed before running):**
1. **Physical prior graph `phys`:** the union of Complex Portal, CORUM and PDB co-structure. Edge weight = number of
   sources; each protein keeps its top 20 neighbours.
2. **Knowledge-regularised reference factorisation.** ALS loadings on the external reference (the Amendment-11 corpus
   minus the target, **minus brunner2022**, plus CCLE and ProteomeHD). Each warm protein's loading is shrunk toward the
   mean of its `phys` partners with strength κ·n_obs, κ = 0.25 (`transfer.fit_loadings(kg=..., kappa=0.25)`). Rank 20,
   λ = 5, 15 iterations, and ≥ 20 observations are as before. Cold proteins get the `phys` cold start (Amendment 10).
3. **Within-target complex regression:** `refpanel.rpr` with the `phys` neighbours (m = 20, λ = 10).
4. **Four-way stack** of kNN, low-rank, step-2 projection and step-3 RPR. The weights are on the simplex grid with step
   0.25, chosen on the 10% inner mask of observed entries (overall squared error, as before).

**Comparators:**
- kNN;
- the kNN + low-rank stack;
- `stack3_plain`: kNN + low-rank + projection with κ = 0 and no RPR (the A11 method, re-fit on the same reference);
- `bim_null`: BIM with `phys` neighbour lists relabelled at random (`refpanel.permute_neighbours`, seed 0), used in both
  the factorisation and the RPR.
- Also reported: `stack3_kg` (no RPR component).

**Dev:**
- The Amendment-9 grid: 3 targets, MCAR and MNAR, seeds 100–102, 30 cells primary; 200 cells / all of Petrosius reported.
- **Gate:** H13a and H13b both hold (point estimates, mean over seeds) on ≥ 2/3 targets under both masks at 30 cells.

**New holdout: brunner2022.**
- 434 HeLa single cells in 5 cell-cycle groups (timsTOF DIA, a different technology).
- It has only ever been part of the external *reference*. It was never imputed or scored, and no decision used its accuracy.
- From now on it is removed from every reference. It joins `EXCLUDED_FROM_CORPUS` for Amendment 13.
- **If the dev gate passes:** confirmation on brunner2022 uses 30-cell subsamples stratified by `CellCycleStage` (6
  per stage), seeds 7000–7009, MCAR and MNAR, with in-complex RMSE.
  - H13a passes if the paired-bootstrap CI is < 0 and the difference is ≤ −1% vs kNN and vs the kNN + low-rank stack.
  - H13b passes if the CI is < 0 and the difference is ≤ −0.25% vs `stack3_plain` and vs `bim_null`.
  - All of this must hold under both masks.
- **If the dev gate fails:** brunner2022 stays untouched and we report the result.

## Amendment 14: refined hypothesis H14, tested once on the brunner2022 holdout (committed before the holdout is evaluated)

**Honest status:** the curated-complex-specific hypothesis (H13b) has no dev support. In dev, the prior is correct but
redundant once ~10k external cells are used.

**Refinement.** Amit asked to tweak the hypothesis slightly and to integrate as much biological prior knowledge as
possible. We broaden "biological prior" from *curated* knowledge to *empirical* prior knowledge: protein co-regulation
learned from public data (16 single-cell datasets, excluding all holdouts, plus the CCLE and ProteomeHD bulk atlases).

**H14:** on 30-cell subsamples of an unseen dataset, `stack3_plain` has lower **in-complex RMSE** than:
- kNN;
- the kNN + low-rank stack;
- its own permuted-reference null (`transfer.permute_rows`, seeds 0 and 1, mean);

by ≥ 1% with a 95% paired-bootstrap CI below 0, under both MCAR and MNAR (6 primary comparisons). `stack3_plain` is
kNN + low-rank + reference projection, with simplex weights chosen on an inner mask of observed entries.

**Disclosure:**
- H14 was formulated *after* seeing dev results. In Amendment 13 dev it met the criteria vs kNN and the kNN + low-rank
  stack on Petrosius and SCoPE2, under both masks.
- On khan2023 a closely related test (overall RMSE, Amendment 12) did not pass vs the kNN + low-rank stack.
- This is therefore a *replication attempt with a refined endpoint* on a new, untouched dataset of a different technology.
  It is the single use of brunner2022.

**Secondary (reported, not part of H14):**
- overall RMSE;
- BIM vs `stack3_plain` (curated-complex-specific);
- BIM vs `bim_null`.

**Details:** `benchmarks/PREREGISTRATION_v7.md` and `benchmarks/confirm_v7.py`.

---

# Phase B: post-hackathon scaling program (from Oct 23, 2026)

**Submission snapshot:** the OSI Hackathon entry, as judged, is `main` at `076a946`. All later work happens on
`claude/happy-hamilton-kdexid` and never changes the numbers of the submitted claims.

## Amendment 15: goal, data, holdouts and protocol for scaling (committed before any Phase-B data is downloaded or any run)

**Goal (Amit, Oct 23):**
- Push the gain of prior-informed imputation over the best data-only method from 2–5% toward 20–30%.
- Do this with much more biological data, new representations and bigger models.
- The goal is aspirational. We report what we get.

**Primary metric:**
- Relative RMSE vs `stack_lowrank` (kNN + low-rank, inner-mask weights), on 30-cell subsamples, MCAR and MNAR, all
  proteins and in-complex.
- Secondary: 10- and 100-cell regimes.

**Holdouts (reserved now, before any Phase-B work):**
- **H-A: `wu2026`** (prenatal human brain, 2268 cells, DIA). It is removed from every reference from now on. It was
  reference-only in Amendments 11–14 and was never imputed or scored.
- **H-B: `derks2022`**, to be rebuilt from *unimputed* precursor data (scpdata's protein assay is imputed). Only shape
  and missingness may be inspected while building.
- Each holdout is used **once**, for the single best Phase-B method, with a preregistration committed beforehand. The
  pass rule is the Amendment-14 rule (CI < 0 and a 1% SESOI vs kNN, the kNN + low-rank stack, and the method's own
  scrambled-reference null).

**Dev protocol:**
- **Targets:** leave-one-dataset-out over all non-holdout single-cell datasets with ≥ 100 cells and cell labels.
  That is the original three dev targets plus any reference dataset that qualifies, a much larger dev set than in
  Phase A.
- **Seeds:** dev seeds 100–102.
- **Rule:** a new idea moves forward only if its *mean* gain over the current best method is ≥ 1% across dev targets,
  and it beats its own scrambled control.
- **Logging:** every idea is listed in `docs/EXPLORATION_LOG.md` with its result, including failures.

**Workstreams:**
1. **More reference data (bulk):**
   - Sanger ProCan cell-line proteomes;
   - GTEx tissue proteomes;
   - CPTAC pan-cancer proteomes;
   - any other public protein × sample matrix with ≥ 100 samples.

   Each needs SHA-256 provenance and a `reference()` module like `reps/ccle.py`.
2. **More reference data (single-cell):** additional public SCP datasets (PRIDE/MassIVE/Zenodo), built with the
   Amendment-11 rules.
3. **Reference selection:** weight reference samples or datasets by similarity to the target, analogous to selecting
   the most similar haplotypes.
4. **Nonlinear transfer:** a masked autoencoder pretrained on the full reference, with only the cell embedding fitted on
   the target. This is the nonlinear version of reference projection.
5. **Peptide-level priors:** impute peptides, then aggregate, where complex stoichiometry constrains ratios.
6. **Head-to-heads:** GIRAFI, Lupine and omicsGMF, if their code is available, as extra baselines.
