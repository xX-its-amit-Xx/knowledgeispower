# King's log

The king appends one entry per round:
- what each queen tried;
- which queens are in a local minimum;
- each queen's directive for the next round, including any serendipitous pivots.

## Round 0 · 2026-10-23 · opening directives

**State of the colony:**
- Nothing built yet beyond the harness check (`sanity`: cand_stack 1.000 / 1.000, `results/phaseC/sanity_fast.json`).
- Champion: H15 (kNN + low-rank + adapted projection, λ_v = 100, 10 extra), confirmed on wu2026.
- Stuck queens: none (no queen has had a round yet).

**Dead ends no queen may repeat** (from `docs/EXPLORATION_LOG.md`):
- More generic reference data: bulk atlases (CCLE, ProteomeHD, ProCan, GTEx, CPTAC, ATLANTiC) and 13 → 25 single-cell
  datasets both plateaued.
- *Hard* reference selection (top-k most similar cells or top-3 datasets by data similarity): no gain over the full
  reference.
- Masked autoencoder: worse than kNN; the stack gives it zero weight.
- Curated complexes or a κ prior on the loadings, protein-side RPR / pkNN / IMC with any of the 22 Amendment-9
  representations (GO, HPA annotations, UniProt KW, STRING, hu.MAP, ProtT5, PaxDb, …): real but redundant once the
  reference is in.
- Knowledge-graph cold start for unseen proteins: ≈ its null.
- Tuning λ_v (30–1000) and `extra` (0–20) as global scalars: done; 100 / 10 is the optimum.

**Why the directives look like this:**
- The ideal in-domain reference reached 0.855 vs data-only, and the champion is at about 0.936. The gap is reference
  *relevance*, not reference *size*.
- So every queen gets either a new information source about relevance (metadata, cross-dataset stability, cell-type
  programs) or a new encoding (sparse direct couplings instead of low-rank only).
- The four ideas below cover different ground on purpose: dataset level (bio), protein level (rep), new modality (mod),
  new model class (abd).

### Directive · bio (biological priors): `ttprior`, a tissue/technology-aware reference prior (Amit's seed)
- **Build:** a metadata table for the 21 reference datasets plus the dev targets.
  - Fields: species, sample type (cell line / primary cells / tissue), lineage or tissue (e.g. myeloid, neural,
    hepatic, epithelial cancer line), acquisition (TMT/SCoPE carrier DDA, label-free DDA, DIA, plexDIA/mTRAQ), and
    instrument family.
  - Every field is taken from the dataset's PRIDE record or paper, cited in the idea doc. Write "unknown" rather than
    guess.
- **Model:** a weighted `fit_loadings`, where each reference dataset's rows get a weight
  w_d = exp(−dist(meta_d, meta_target) / τ), and the weighted loadings become the prior in `adapt_project`.
  - Keep the weights **soft**: hard top-k selection is dead (Phase B).
  - Report τ ∈ {∞ (= champion), 2, 1, 0.5}. Also report separate tissue-only and technology-only kernels, so we learn
    which axis carries the signal.
- **Second variant, if time allows:** two priors (a technology-matched V and a tissue-matched V) concatenated as
  separate factor blocks in `adapt_project`.
- **Null:** the same weighting machinery, with the metadata rows permuted across reference datasets. This gives the
  same weight distribution with the wrong biology.
- **Why:** this is the first use of *a priori* relevance information rather than data-estimated similarity.
  Hard data-similarity selection failed partly because 30 cells give a noisy similarity estimate. Metadata is
  noise-free.

### Directive · rep (representation & math): `ebshrink`, per-protein adaptation strength λ_v,j by empirical Bayes
- **Build:**
  - Fit loadings separately on each reference dataset (or on leave-one-reference-dataset-out folds), after Procrustes
    alignment to the pooled V.
  - Estimate each protein's between-dataset loading variance τ²_j, which measures how stable that protein's
    co-regulation is across labs, tissues and technologies.
  - Set λ_v,j = c · σ² / τ²_j, clipped to [10, 1000], with c chosen so that the median λ_v,j = 100 (the champion's
    value). The candidate is then a per-protein, not a global, quadratic penalty in the V-step of `adapt_project`
    (diagonal in j, so it is still a closed-form solve per protein).
- **Why:** the global λ_v is already at its optimum, so a scalar sweep is dead. Cross-dataset stability is *new
  information* that the pooled fit throws away: proteins whose co-regulation is universal should be held to the
  reference, and context-specific ones should be freed. This is the BLUP / James–Stein idea applied per protein.
- **Null:** the same machinery with the λ_v,j vector permuted across proteins. This keeps the distribution of
  strengths and destroys which protein gets which.
- **Fallback, if `ebshrink` dies at the fast screen:** switch encoding. Use an RBF/GP kernel on the reference loadings
  (a kernel ridge across proteins) as a fourth stack member, instead of linear factors.

### Directive · mod (new data modalities): `scrnaprior`, cell-type programs from a scRNA atlas as an extra loading prior
- **Never used before:** single-cell **RNA** cell-type expression. The HPA protein annotations used in Amendments 5/9
  were tissue/localisation *annotation features* in RPR/IMC, which is a different object and a dead end.
- **Source:** the Human Protein Atlas single-cell-type RNA consensus tables (gene × cell type, or gene × tissue-cell-type
  cluster nTPM), downloadable from proteinatlas.org.
  - Check the size before downloading, and abort above 2 GB. Store it in `data/cache/phasec/scrnaprior/`.
  - Record the URL, licence and SHA-256 in the idea doc.
  - If a CELLxGENE Census pseudo-bulk co-expression matrix fits in the budget, it is an acceptable alternative. Do not
    stream the full Census.
- **Model:**
  - Map genes to the target's UniProt accessions (`ids.group_accessions`).
  - Take log1p, z-score per gene, then take the top-k SVD gene loadings (k = 5–10).
  - Append them to the prior as an extra factor block in `adapt_project`, with its own ridge strength. Do not replace
    the proteomic V.
- **Null:** the same SVD block with the gene → protein mapping permuted.
- **Why:** bulk atlases failed because they average over cell types. A cell-type-resolved atlas encodes exactly the
  programmes that separate cells in a 30-cell subsample, at a scale (hundreds of cell types) that no proteomic
  reference has.
- **Fallback:** HPA subcellular localisation *imaging* (IF image-derived compartment intensities), used as a
  compartment-block structure on the extra factors.

### Directive · abd (innovation abduction): `dca`, direct-coupling analysis from MSA co-evolution
- **Analogy:** in protein-family MSAs, raw residue covariance mixes direct and indirect (transitive) couplings.
  Direct-coupling analysis inverts the covariance with a sparse/regularised inverse, and only then do the couplings
  predict 3D contacts (Weigt et al., *PNAS* 2009; Morcos et al., *PNAS* 2011).
- **Our version:** every model so far encodes protein co-regulation as *low-rank* covariance. Build the
  **precision-matrix** encoding instead:
  - Fit a graphical lasso (or a shrinkage inverse) on the z-scored reference, restricted to the target's proteins.
    For feasibility, use the top ~1,500 most-observed proteins, or a low-rank-plus-sparse decomposition in the spirit
    of Chandrasekaran-style latent-variable graphical models.
  - Impute each cell's missing proteins by the Gaussian conditional mean E[x_m | x_o] = −Θ_mm⁻¹ Θ_mo x_o.
  - Feed the result into the stack as a candidate.
- **Null:** the same graphical-lasso fit on a reference whose protein columns are permuted within each reference
  dataset (wrong couplings, same sparsity machinery).
- **Cross-pollination:** if the sparse couplings carry signal that is not in the low-rank part, they are a *data-driven
  complex graph*.
  - bio could then use them in place of curated complexes, which were redundant.
  - rep could use their degree as a feature for λ_v,j.
- **Next analogy in the queue, so the next round need not start cold:** the Netflix-Prize baseline predictors
  (global + user + item biases fitted before the factors). Here those would be explicit per-reference-dataset ×
  per-protein offsets, which the z-scoring currently assumes away.

**Rules reminder for all workers:**
- `--pool fast --jobs 2`; go to the full pool only at cand_stack ≤ 0.99 *and* below null.
- Dev data only. Every number in a doc comes from `results/phaseC/<id>_fast.json`.

## Round 1 · 2026-10-23 · review and directives for round 2

**Verification.** The king re-read every number in the four worker reports from `results/phaseC/*_fast.json` and
`results/phaseC/dca_couplings_diag.json`: summary means, in-complex means, per-target/mask values, stack weights and
the CORUM diagnostic. Every number matched; no record needed correcting.

**What each queen tried** (fast pool, rmse_all mean relative to the champion, all / in-complex; bars are
cand_stack ≤ 0.99 *and* a margin ≥ 0.005 below the null):

| queen | idea | cand_stack | null_stack | cand_solo | verdict |
|---|---|---|---|---|---|
| bio | `ttprior` (FULL metadata kernel, τ = 1) | 0.9951 / 0.9953 | 0.9991 / 1.0005 | 0.9944 | dead, closest of the round (margin 0.0040) |
| rep | `ebshrink` (per-protein λ_v,j) | 0.9986 / 0.9983 | 1.0000 / 0.9999 | 0.9991 | dead (v2: 0.9988) |
| mod | `scrnaprior` (HPA scRNA SVD block) | 0.9973 / 0.9974 | 0.9987 / 0.9988 | 0.9989 | dead |
| abd | `dca` (shrinkage-inverse precision) | 1.0004 / 1.0006 | 1.0000 / 1.0000 | 1.032 | dead (stack weight ≈ 0) |

**What the round taught us:**
- Every lever that only reshapes the *prior on V* is saturated: re-weighting datasets, per-protein strength,
  extra cell-type factors and a precision encoding of the same reference all land within 0.5% of the champion. The
  adapted projection already spans the reference's co-regulation content (dca: cand_solo 1.032, but stack weight 0.01).
- The one consistent signal is **technology**, not tissue: tech-only weighting 0.9970 vs tissue-only 0.9992 (= its
  null 0.9997), and the ttprior gain sits on label-free Astral-DIA targets (sadiku2025 MCAR 0.9797, petrosius2023
  0.9905–0.9943, weiss2026_liver 0.9932–0.9956) while TMT-carrier targets stay at 0.997–1.001. The worker's read is that
  the weighting partly corrects platform-specific *marginal* behaviour (offsets, scales, compression), not
  co-regulation. Note that the tech-only null (0.9980) is close, so part of the effect is generic re-weighting.
- scrnaprior's small real-over-null signal is also on tissue-like targets (weiss2026_liver 0.9903 / 0.9908).
- The champion standardises the target with per-protein mean/sd estimated from 30 masked cells (`adapt_project`:
  `np.nanmean(Xm, 0)`, `np.nanstd(Xm, 0)`). Nothing in the colony has yet touched that step, which is noisy and
  MNAR-biased. Round 2 attacks it from two sides.

**Stuck queens:** none by the charter's definition (each queen has one failed idea, not two). All four are on *one
strike*, and all four ideas died, so every queen is moved off her round-1 line rather than given a revision of it. A
second dead idea next round makes a queen stuck.

**Dead ends added this round (no queen may repeat):** metadata-kernel dataset weighting of `fit_loadings` (any τ,
tissue or tech kernel, two-block); per-protein λ_v,j from cross-dataset loading variance; scRNA cell-type SVD as an
extra factor block (HPA single-cell-type table, k = 8); shrinkage-inverse / DCA Gaussian conditional mean on the
reference; RBF/GP kernel ridge across proteins on the reference loadings (rep's listed fallback; not built, but it is a
protein-side similarity smoother of the kind already dead as pkNN/RPR, and dca shows the reference's pairwise content is
covered).

**Best idea so far:** `ttprior` (bio): cand_stack 0.9951 / 0.9953, null 0.9991 / 1.0005, cand_solo 0.9944. Dead at
the bar; the champion stands.

### Directive · bio: `ctprior`, marker-knowledge cell-type-aware *cell* weighting (pivot from dataset to cell level)
- **Why:** dataset-level tissue similarity carried nothing because most reference datasets mix cell types and the
  tissue label is too coarse. Amit's tissue-aware seed is kept but moved to the resolution where scrnaprior saw its
  (small) real signal: cell types. The information source is curated *marker knowledge*, not data similarity, so it
  does not repeat dead hard top-k cell selection.
- **Build:**
  - Marker sets: the HPA single-cell-type table already cached by scrnaprior
    (`data/cache/phasec/scrnaprior/rna_single_cell_type.tsv.zip`, provenance row merged this round). Define each cell
    type's markers as its top-enriched genes (e.g. ≥ 4-fold over the next type, top 50), mapped to accessions with
    `ids.group_accessions`.
  - Score every reference cell and every target cell for every cell type with a missingness-robust signature score
    (mean z of observed markers minus mean z of observed random background genes, ≥ 3 observed markers required).
  - Soft cell weights w_i = exp(−‖s_i − s̄_target‖² / τ) (s = cell-type score vectors), renormalised to mean 1, in a
    weighted `fit_loadings` (reuse ttprior's vectorised weighted V-step and its cache layout). τ from {∞, median distance,
    ½ median} chosen a priori, all reported.
- **Null:** the same machinery with marker sets drawn as random gene sets matched in size and in observation rate
  (keeps weight distribution, destroys the biology).
- **Diagnostic to report:** per-target effective number of reference cells, and which reference datasets gain weight.

### Directive · rep: `copula`, a nonparanormal (rank-Gaussian) encoding per dataset (new encoding)
- **Why:** ttprior says the residual signal is platform-specific marginal behaviour. TMT reporter-ion ratio compression
  and DIA/LFQ intensity-dependent noise are monotone *nonlinear* distortions that per-protein z-scoring cannot remove,
  so the loadings learn a mixture of platform shapes. A Gaussian-copula encoding (Liu, Lafferty & Wasserman 2009,
  *JMLR*, the "nonparanormal") removes the marginal shape and keeps only rank-dependence, which is the
  co-regulation content.
- **Build:** within each reference dataset, map every protein to normal scores (rank → Φ⁻¹ with the observed-only
  ranks); refit `fit_loadings` on that; in the target, map observed entries to normal scores the same way, run
  `adapt_project`, and back-transform predictions with the target's own empirical quantiles per protein (fall back to a
  linear map where < 5 observed). Unit check: with an identity marginal transform, reproduce the champion exactly.
- **Null:** the same pipeline with each reference dataset's per-protein marginal transforms assigned to a random other
  protein of the same dataset (same shapes, wrong proteins).
- **Keep orthogonal to abd:** rep owns the *shape* of the marginals; abd owns the *location/scale* of the target's
  standardisation. Do not add offsets or shrinkage of the target mean here.

### Directive · mod: `pepev`, peptide-level evidence as per-entry precision weights (new modality)
- **Why:** every model so far treats each observed protein value as equally reliable. Peptide/precursor-level evidence
  (number of precursors or PSMs, q-values, carrier/reporter signal) is information that has never entered the colony,
  and it is the "peptide-level evidence" modality in mod's charter. HPA imaging is held back as the fallback because
  scrnaprior shows that extra annotation blocks on V are near-saturated.
- **Build:**
  - For the six fast-pool targets only, find the precursor/PSM-level tables behind the existing corpus builds (scpdata
    peptide/PSM assays for leduc2022 and gregoire2023_mixCTRL; the DIA-NN / Spectronaut reports or archives already
    named in `data/PROVENANCE.md` for weiss2026_liver, sadiku2025, px063590_ADbrain; the Dataverse record for
    petrosius2023). Check size first; < 2 GB in `data/cache/phasec/pepev/`. Record URL, licence and SHA-256.
  - Build per-entry precision w_ij (e.g. log(1 + n_precursors) or a PSM-count/q-value score), aligned to the same
    cells as the harness. `ctx` gives the target and seed, not cell ids, so recompute the subsample with
    `explore_phaseb.subsample` and assert the observed pattern matches `Xm` before using any weight.
  - Use the weights as heteroscedastic weights in the U-step (per-cell weighted ridge) of a local copy of
    `adapt_project`; all-ones weights must reproduce the champion exactly.
- **Null:** precision weights permuted among the observed entries of each cell (same distribution, wrong entries).
- **Abort rule:** if peptide-level data cannot be aligned for at least 4 of the 6 fast targets, stop early and switch
  to the fallback (HPA IF-image compartment intensities as a block structure on the *extra* factors).

### Directive · abd: `nfbase`, Netflix-Prize baseline predictors for the target's standardisation (her queued analogy)
- **Analogy:** in the Netflix Prize, explicit global + user + item biases fitted with shrinkage before the factors
  were among the largest single gains (Koren, *KDD* 2008; Bell & Koren, *SIGKDD Explorations* 2007). Here the "item bias" is a protein's
  location and scale, and the champion estimates it for the target from 30 masked cells only, which is noisy and, under
  MNAR, biased upwards.
- **Build:**
  - From the *unstandardised* reference (per-dataset log values before `_z`), estimate per-protein relative location
    (protein mean minus dataset grand mean) and log-scale per dataset.
  - Fit a hierarchical (empirical-Bayes) model: dataset value = protein effect + technology-group × protein effect
    + noise, using ttprior's cited metadata table for the technology group (TMT carrier / LFQ-DDA / DIA, and
    instrument family).
  - For the target, shrink the 30-cell per-protein mean and sd toward the prediction from its technology group, with
    the weight set by the number of observed values (James–Stein / BLUP). Then run `adapt_project` unchanged with
    those baselines in place of `np.nanmean` / `np.nanstd`. Unit check: zero shrinkage reproduces the champion exactly.
- **Null:** the same hierarchy with protein labels permuted within each reference dataset's baseline vectors (and,
  as a second null, technology labels permuted across datasets).
- **Next analogy queued:** extreme deconvolution from astronomy (Bovy, Hogg & Roweis 2011, *Ann. Appl. Stat.*), which
  fits a density from noisy, incomplete data with *per-point* noise covariances; it pairs with mod's per-entry
  precision weights if `pepev` shows they carry signal.

**Cross-pollination:** abd's Netflix analogy and bio's technology finding produce `nfbase`; mod's `pepev` is the data
source abd's next analogy (extreme deconvolution) needs; bio reuses mod's HPA download; rep's `copula` and abd's
`nfbase` split the marginal problem into shape vs location/scale so that neither duplicates the other.

**Rules reminder:** `--pool fast --jobs 2`; the full pool only at cand_stack ≤ 0.99 and ≥ 0.005 below the null; dev
data only; every number from `results/phaseC/<id>_fast.json`; provenance rows in the idea doc.

## Round 2 · 2026-10-23 · review and directives for round 3

**Verification.** The king re-read every number in the four worker reports from `results/phaseC/*_fast.json` and
`pepev_diag.json`: summary means, in-complex means, per-target/mask values, stack weights and splits ≤ 0.99. All
matched. One gap was filled: the copula worker recorded `copula_tgtonly`'s stack weight as "not recorded"; the JSON
gives 0.1458 (`mean_cand_weight_in_stack`).

**What each queen tried** (fast pool, rmse_all mean relative to the champion, all / in-complex; bars are
cand_stack ≤ 0.99 *and* a margin ≥ 0.005 below the null):

| queen | idea | cand_stack | null_stack | cand_solo | verdict |
|---|---|---|---|---|---|
| bio | `ctprior` (HPA-marker cell-type soft cell weights) | 0.9990 / 0.9986 | 0.9990 / 0.9987 | 0.9988 | dead (margin 0.0000; τ = ½ median 0.9980 vs null 0.9982) |
| rep | `copula` (nonparanormal encoding) | 1.0005 / 1.0011 | 1.0004 / 1.0004 | 1.0193 | dead (reference-side-only 0.9980 vs null 0.9991) |
| mod | `pepev` (precursor counts as U-step precision weights) | 0.9996 / 0.9987 | 0.9998 / 0.9996 | 0.9988 | dead (v2 protein-level weights 0.9969 vs null 1.0003, margin 0.0034) |
| abd | `nfbase` (Netflix baselines for the target mean/sd), ψ = 1 | 0.9868 / 0.9859 | 0.9867 / 0.9863 | 0.9864 | dead as a prior; ψ = 4: **0.9663 / 0.9647** vs null 0.9691 / 0.9688, 12/12 splits ≤ 0.99 |

**What the round taught us:**
- **The prior on V is closed.** Five ways of re-weighting or re-encoding the reference (rounds 1–2) all land within
  0.5% of the champion and within 0.0040 of their nulls. ctprior's weight shifts are the same on every target
  (schoof2021 up, leduc2022/specht2021 down), which is per-dataset score scaling, not cell identity.
- **The target's standardisation is where the error is.** nfbase's protein-permuted null (shrink the target's own
  30-cell mean and log-sd toward a calibrated constant, ψ = 4) is 0.9691, ≤ 0.99 on all 12 splits: the largest Phase-C
  effect so far, and data-only. It is super-additive (post-hoc meanonly_psi4 0.9920, sdonly_psi4 0.9832, joint 0.9663);
  scale is the bigger half. Protein-specific *reference* baselines add only on DIA MCAR (weiss2026_liver 0.9490 vs null
  0.9645, petrosius2023 0.9751 vs 0.9827, sadiku2025 0.9767 vs 0.9815) and nothing on TMT.
- **Missingness is information we have not modelled.** pepev's precursor counts track the value within protein in DIA
  (Spearman 0.21–0.44), so they help MCAR and hurt MNAR (weiss2026_liver 0.9759 / 1.0121). Nobody has modelled *why* an
  entry is missing. Also, the stack is fitted on a random inner mask, so MNAR-specific wins get diluted (pepev_v2,
  sadiku2025 MNAR: cand_raw 0.9676, cand_stack 0.9959). Workers should report cand_raw by mask too.
- Marginal *shape* is not the platform effect: TMT is the most Gaussian platform (copula_diag), and the target-side
  quantile decode hurts (tgtonly cand_raw 1.0284).

**Stuck queens: all four** by the charter's rule (two consecutive ideas without cand_stack ≤ 0.99 *and* a null margin).
- bio: ttprior → ctprior. Two tissue-aware V priors died; tissue at dataset and cell-type resolution are both null.
- rep: ebshrink → copula. Both re-parameterised the same reference; pivot to the missingness mechanism.
- mod: scrnaprior → pepev. Pivot from per-entry weights to protein-intrinsic information.
- abd: dca → nfbase. Formally stuck (nfbase failed its null), but her machinery found the lead, so she gets to turn
  it into a champion revision.

**Dead ends added this round (no queen may repeat):** soft per-cell (marker-score) re-weighting of the reference;
nonparanormal/rank-Gaussian encoding on either side (incl. per-protein empirical-quantile decoding of the target);
precursor/PSM counts as WLS precision weights on U or V (per entry or per protein); technology-group × protein
hierarchical baselines from the reference (tech term = flat = technull).

**Best idea so far:** the *machinery* of `nfbase` at ψ = 4: cand_stack 0.9663 / 0.9647, cand_solo 0.9663, 12/12
splits ≤ 0.99. Its biological content is not significant (null 0.9691 / 0.9688; margin 0.0028 / 0.0041), so it is
a data-only fix to the champion's standardisation, not a prior. The champion stands until `tgtshrink` passes the full
pool. **Rule for round 3:** bio and mod build on nfbase's `adapt_project_mb` (injected mu/sd) and must beat a
*shrink-to-constant* baseline at the same ψ, not just the champion, so the data-only effect is not credited to them.

### Directive · abd: `tgtshrink`, the champion's standardisation done properly (her round-2 lead)
- **Why:** her own null is the biggest effect in Phase C. Make it a clean, data-only champion revision and find out
  whether the ψ = 4 effect holds across the 15 LODO targets.
- **Build:** `phasec/tgtshrink.py` on top of `nfbase.adapt_project_mb`.
  - Per protein j in the target: m_j, log s_j from observed entries, n_j observed.
  - log-sd: shrink toward a smooth trend in (m_j, n_j / n_cells), fitted across the target's proteins (limma-trend:
    Smyth 2004, *Stat. Appl. Genet. Mol. Biol.* 3:3; a low-df spline is enough). Mean: shrink toward a trend in
    n_j / n_cells (under MNAR, sparsely observed proteins have upward-biased means).
  - ψ (noise-variance inflation) chosen per target from {1, 2, 4, 8, 16} on an *inner* mask of the masked matrix, never
    on the evaluation mask. Report fixed ψ = 4 alongside.
  - Unit checks: ψ = 0 reproduces the champion; trend = constant at ψ = 4 reproduces `nfbase_psi4`'s null (0.9691).
- **Null:** the same pipeline with the trend covariates (m_j, n_j) permuted across proteins (equivalent to shrinking
  to a constant).
- **Pool:** fast, then **full** if cand_stack ≤ 0.99. Bar for adoption: a full-pool win vs the champion on both masks;
  the trend part counts only if it beats its null by ≥ 0.005. A win means a `PREREGISTRATION_v9.md` and asking Amit.
- **Next analogy queued:** Heckman selection correction (econometrics; Heckman 1979, *Econometrica* 47:153–161),
  handed to rep this round.

### Directive · rep: `selmnar`, a selection model for the missing entries (pivot: from encoding the observed values to encoding the missingness)
- **Why:** abd's Heckman analogy plus pepev's finding. In DIA targets a value is missing *because* it is low; the
  champion treats every missing entry as uninformative. Nothing in Phases A–C has modelled the missingness mechanism.
- **Build:** `phasec/selmnar.py`, a local copy of `adapt_project` with an EM U-step.
  - Per protein a probit detection curve P(observed | x) = Φ((x − c_j) / σ_d), with c_j shrunk toward a per-cell
    detection threshold; fit from the observed-value distribution and the observed fraction.
  - E-step: for each missing entry, combine the current low-rank prediction (Gaussian) with the probit selection
    term; use the posterior mean of that tilted Gaussian as the pseudo-observation. M-step: the champion's ridge U-step.
    3–5 EM iterations.
  - Unit check: σ_d → ∞ (no selection) reproduces the champion.
- **Null:** the same EM with the missing-indicator pattern permuted within each cell (the same counts of missing
  entries per cell, the wrong entries).
- **Report** cand_raw and cand_stack *per mask*: the expected signature is a gain on MNAR and a smaller loss on MCAR.
  If MCAR is hurt, a variant that fits σ_d per target on the inner mask is allowed (preregister it in the doc first).

### Directive · mod: `seqprior`, protein-intrinsic detectability features as a prior on location, scale and detection (new modality)
- **Why:** nfbase found protein-specific location/scale information only where the reference covers the protein,
  and only for DIA. Detectability is a property of the sequence (tryptic peptide count, length, hydrophobicity), so it
  generalises to every target protein, and it is exactly what drives MNAR. This is Amit's synthetic-features seed.
- **Build:** `phasec/seqprior.py`.
  - Data (< 100 MB): UniProt human reviewed proteome FASTA (rest.uniprot.org stream) and HPA subcellular location
    (`subcellular_location.tsv.zip`). Record URL, licence and SHA-256 for each.
  - Features per protein: length, number of unique tryptic peptides of 7–30 aa, MW, GRAVY, pI, HPA main location
    one-hot.
  - A ridge or gradient-boosted sub-model, trained LODO on the *reference* datasets' per-protein relative location,
    log-sd and observed fraction (one model per acquisition group; reuse nfbase's baseline tables), predicts those three
    for the target.
  - Use the predictions as protein-specific shrinkage targets in `nfbase.adapt_project_mb` at ψ = 4.
- **Null:** feature rows permuted across proteins before training. **Required baseline:** the same pipeline with
  shrink-to-constant (nfbase_psi4's null, 0.9691). Credit only the margin below both.

### Directive · bio: `varprior`, protein turnover and complex buffering as a prior on cell-to-cell variability (pivot: tissue → protein dynamics)
- **Why:** both tissue-aware V priors died. The scale half of the standardisation fix is the larger half
  (sdonly_psi4 0.9832), and biology predicts per-protein variability: long-lived proteins average out bursts, and
  complex subunits are buffered by stoichiometric degradation. This replaces the tissue line with the "half-life"
  synthetic feature Amit listed.
- **Build:** `phasec/varprior.py`.
  - Data: measured half-lives from Mathieson et al. 2018, *Nat. Commun.* 9:689 ("Systematic analysis of protein turnover
    in primary cells"; use the supplementary table, record URL/licence/SHA-256), plus CORUM complex membership already
    cached for the project.
  - Fit a log-sd model on the reference: per-protein log-sd (per dataset, relative to the dataset median) ~ log
    half-life + in-complex + complex size + observed fraction, LODO.
  - Its prediction is the log-sd shrinkage target for the target in `nfbase.adapt_project_mb` (ψ = 4; location shrunk
    to a constant as in nfbase's null, so only the scale prior differs).
- **Null:** half-life and complex covariates permuted across proteins (observed fraction kept). **Required baseline:**
  shrink-to-constant at ψ = 4. Credit only the margin below both.
- **Coordinate with mod:** bio owns *scale from dynamics/complexes*; mod owns *location and detection from sequence*.
  If both pass, the king will test them together.

**Cross-pollination:** abd's Netflix analogy found the lead (the standardisation step); her queued Heckman analogy
becomes rep's selection model; pepev's MNAR finding motivates both `selmnar` and mod's detectability features; bio's
turnover prior and mod's sequence prior split the shrinkage target into scale vs location/detection, so the two do not
duplicate each other.

**Housekeeping flag:** `transfer.extended_corpus()` returns holdout names (brunner2022, wu2026) and relies on callers
filtering with `PHASE_B_HOLDOUTS` (as `explore_phaseb.reference_datasets` does). Workers must use
`explore_phaseb.reference_datasets`, never `extended_corpus()` directly.

**Rules reminder:** `--pool fast --jobs 2`; the full pool only at cand_stack ≤ 0.99 and ≥ 0.005 below the null (for
`tgtshrink`, ≤ 0.99 vs the champion); dev data only; every number from `results/phaseC/<id>_fast.json`; provenance rows
in the idea doc.

## Round 3 · 2026-10-23 · review and directives for round 4

**Verification.** The king re-read every number in the four worker reports from `results/phaseC/*.json`: summary
means (all / in-complex), per-mask means, per-target/mask values quoted in the reports, stack weights and splits
≤ 0.99 (`tgtshrink_full.json` included). All matched. One rounding note: selmnar's MCAR null_stack mean is 0.99985
(reported 0.9998, recorded in the ledger as 0.9999); varprior's "weiss2026_liver 0.9873 vs 0.9659" are cand_stack values.

**What each queen tried** (fast pool unless marked; rmse_all mean relative to the champion, all / in-complex; bars:
cand_stack ≤ 0.99 *and* ≥ 0.005 below the null; bio and mod also ≥ 0.005 below shrink-to-constant at ψ = 4):

| queen | idea | cand_stack | null_stack | cand_solo | verdict |
|---|---|---|---|---|---|
| abd | `tgtshrink` (shrink target mean/log-sd toward a (mean, obs-fraction) trend, ψ selected) | 0.9597 / 0.9586; **full 0.9683 / 0.9684** | 0.9575 / 0.9563; **full 0.9650 / 0.9653** | 0.9640 | trend dead (worse than null by 0.0022 / 0.0033); **null pipeline is the best Phase-C result** |
| rep | `selmnar` (Heckman probit selection, EM U-step) | 1.0000 / 1.0000 | 0.9985 / 0.9986 | 1.0192 | dead (weight 0.0000); raw MNAR 0.8774 vs null 0.9615, raw MCAR 1.2721 |
| mod | `seqprior` (sequence/HPA detectability → shrinkage targets) | 0.9686 / 0.9689 | 0.9687 / 0.9694 | 0.9687 | dead (= constant 0.9686; ridge 0.9683) |
| bio | `varprior` (half-life + CORUM → log-sd target) | 0.9679 / 0.9686 | 0.9718 / 0.9728 | 0.9678 | dead (margin 0.0039 over null, 0.0007 over constant 0.9686; v2 0.9682 vs 0.9686) |

**What the round taught us:**
- **The lever is shrinkage itself, not where it shrinks to.** Three independent protein-specific shrinkage targets
  (nfbase reference baselines, varprior half-life/complex, seqprior sequence/location) all tie shrink-to-constant
  (0.9686–0.9691 at ψ = 4). tgtshrink's own trend loses to its permuted null. This line is closed.
- **Choosing ψ per target is worth another ~1%.** The selected-ψ constant pipeline (tgtshrink null) is 0.9575 fast
  vs 0.9691 at fixed ψ = 4, and 0.9650 on the full pool (MCAR 0.9736, MNAR 0.9564; 25/30 target-masks ≤ 0.99; worst
  woo2025_PC3 MCAR 1.0000). ψ = 16, the top of the grid, was chosen in 30/48 null calls and tau²(log-sd) sat at its
  floor: the champion's per-protein unit-variance scaling looks harmful in itself. That is an encoding finding.
- **Missingness carries real, value-dependent information** (selmnar raw MNAR 0.8774 vs permuted-indicator null
  0.9615 and σ_d = ∞ 0.9967), but a single parametric tilt applied to every missing entry costs 17–27% on MCAR, and the
  random-inner-mask stack gives it weight 0.0000. Two lessons: (i) any missingness model must calibrate its strength
  on an inner mask that contains MNAR-like holes (tgtshrink's 50/50 inner mask is the template), (ii) part of the MNAR
  gain may be re-learning the harness's value-logistic mask generator, so per-mask diagnostics are mandatory.
- A transferred observed-fraction trend hurts (varprior_fonly 0.9720 vs const 0.9686); sequence features predict the
  reference (OOF R² 0.15–0.19) but cannot move a BLUP that keeps 0.88–0.93 weight on the target's own mean.

**Stuck queens: all four** by the charter's rule.
- bio: ctprior → varprior (three biology-prior lines dead; varprior reached ≤ 0.99 only through the shared machinery).
- rep: copula → selmnar (selmnar's signal is real but the encoding is wrong; adjacent move, not a restart).
- mod: pepev → seqprior (protein-intrinsic features are exhausted for this benchmark).
- abd: nfbase → tgtshrink (both failed their nulls, but her machinery is the best pipeline; she finalises it).

**Dead ends added this round (no queen may repeat):** any protein-specific shrinkage *target* for the target's
mean/log-sd (reference baselines, half-life, complex membership/size, sequence detectability, HPA location, a
(mean, observed-fraction) spline trend, a transferred observed-fraction trend); a single global probit/Heckman tilt of
all missing entries (with or without a MAR mixture or calibrated ρ). Also still closed from Phase A: CORUM graph
penalties on the low-rank model (regularisation, not curation).

**Best idea so far:** the *data-only* selected-ψ shrink-to-constant pipeline (tgtshrink's null): fast cand_stack
0.9575 / 0.9563, full 0.9650 / 0.9653 (MCAR 0.9736, MNAR 0.9564), 25/30 full target-masks ≤ 0.99. It has never been
run as a named candidate with its own control, so round 4's `cshrink` does that. Until then the champion stands.

**Rule for round 4:** every candidate builds on, and is credited only against, the `cshrink` base (the selected-ψ
constant shrinkage), not the bare champion: report cand_stack, cand_solo, cand_raw per mask, the null, and a
`<id>_base` run of the same code with the idea switched off. Any internal calibration uses the 50/50 MCAR/MNAR inner
mask built only from the input matrix (tgtshrink's recipe, rng seed + 31337). Until `cshrink` lands, the other three
reuse `tgtshrink.impute_null`'s pipeline as the base.

### Directive · abd: `cshrink`, the selected-ψ constant shrinkage as an explicit champion revision
- **Why:** her round-3 null is the best pipeline in Phase C (full 0.9650) and selection hits the top of the grid.
- **Build:** `phasec/cshrink.py` on `nfbase.adapt_project_mb`: shrink the target's per-protein mean and log-sd toward
  constants (no permutation step), ψ from {1, 2, 4, 8, 16, 32, 64, ∞} on the 50/50 inner mask; plus an arm with no
  per-protein scaling at all (one global sd; the ψ = ∞ log-sd limit) to test whether unit-variance scaling is the harm.
  Unit check: reproduces tgtshrink's null numbers on the old grid within 0.001.
- **Null:** the same pipeline with ψ drawn at random from the grid per call (seeded); this keeps the machinery and
  destroys the selection information. Also report fixed ψ = 16 and fixed ψ = ∞.
- **Pool:** fast, then full. If full cand_stack ≤ 0.97 on both masks with no target-mask > 1.00, write a *draft*
  `docs/phaseC/PREREGISTRATION_v9_draft.md` for Amit (no holdout is touched).
- **Then:** in her doc, propose two analogies for round 5 aimed at the remaining bottleneck (covariance in the
  target's own domain: the ideal in-domain reference reached −14.5%).

### Directive · rep: `mnarjoint`, missingness as a second channel of a joint low-rank encoding (adjacent pivot from selmnar)
- **Why:** selmnar proved the missingness pattern is informative but a parametric tilt on every hole is the wrong
  encoding. Sportisse, Boyer & Josse (2020, *Stat. Comput.* 30:1629–1643; verify the DOI before citing) impute MNAR
  data by a low-rank fit of the concatenation [X ∣ M]: the cell factors must explain which proteins are missing, so a
  hole is imputed low only where the pattern says so, entry by entry.
- **Build:** `phasec/mnarjoint.py`: on the base pipeline's standardised target, fit shared cell factors U to
  [Z_obs ∣ w·(M − mean)] with the champion's ridge V-step for the value block and a separate free loading block for M
  (logistic or Gaussian surrogate); w from {0, 0.1, 0.3, 1, 3} on the 50/50 inner mask. w = 0 must reproduce the base.
- **Null:** missing indicators permuted within each cell (selmnar's null). Report cand_raw per mask; the bar is the
  usual one *and* MCAR cand_raw no worse than the base by > 0.005.

### Directive · mod: `subthr`, sub-threshold DIA evidence ("forced photometry", abd's astronomy analogy)
- **Why:** astronomers do not discard non-detections; they measure flux at the known position and keep it as a noisy
  upper-limit observation. Our DIA matrices discard precursors that pass a run-level q-value but fail the experiment-
  wide filter. That is NEW information about exactly the naturally missing entries, and it corrects the upward-biased
  means and sds that the shrinkage lever acts on.
- **Build:** `phasec/subthr.py`. Audit first (record in the doc): for weiss2026_liver (PXD062231, pepev's range-read
  route), petrosius2023 (pepev's Spectronaut export) and sadiku2025, does a report exist with per-run precursor
  quantities and q-values beyond what the matrix keeps? Keep downloads < 2 GB in `data/cache/phasec/subthr/`.
  Use evidence **only** for (cell, protein) pairs whose protein is absent from the published matrix for that cell; by
  construction these are never masked/scored entries (masked entries were observed). Add them as low-precision
  observations (precision from the q-value) in the mean/log-sd estimates and the U-step of the base pipeline.
- **Null:** the sub-threshold evidence vectors permuted across proteins (same positions, wrong protein).
- **Pool:** fast with `--targets` restricted to the targets that have evidence (≥ 2 needed; if fewer, stop and report
  the audit). It cannot go to the full pool unless ≥ 5 dev targets have evidence.

### Directive · bio: `techhyper`, a technology/tissue-aware hyperprior on the standardisation (Amit's favourite seed, moved to where the lever is)
- **Why:** with 30 cells, per-target ψ selection on an inner mask is noisy (ψ choices are spread over the grid), and
  ttprior found that technology carries the only metadata signal. Learn how much each kind of dataset should be
  shrunk from the reference datasets themselves.
- **Build:** `phasec/techhyper.py`. For each reference dataset (`explore_phaseb.reference_datasets` only, LODO),
  subsample 30 cells, run the cshrink inner-mask selection and record the ψ-grid loss curve. Prior for the target =
  ttprior's metadata-kernel weighted average of those curves (technology, carrier yes/no, tissue/primary vs cell line);
  posterior = prior × the target's own inner-mask curve; use the posterior-mean ψ. Report technology-only and
  tissue-only kernels too.
- **Null:** metadata rows permuted across reference datasets. **Baseline:** cshrink's flat prior. Credit only the
  margin below both.

**Cross-pollination:** abd's astronomy analogy (forced photometry on non-detections) becomes mod's `subthr`; her
round-3 null becomes everyone's base; rep's selmnar evidence is re-encoded as a joint [X ∣ M] factor model; bio's
tissue/technology line moves from the V prior (dead) to the shrinkage hyperprior (live lever).

**Rules reminder:** `--pool fast --jobs 2`; full pool only at cand_stack ≤ 0.99 and ≥ 0.005 below both the null and
the `<id>_base` run (for `cshrink`: ≤ 0.99 vs the champion); dev data only; every number from `results/phaseC/`;
provenance rows in the idea doc.
