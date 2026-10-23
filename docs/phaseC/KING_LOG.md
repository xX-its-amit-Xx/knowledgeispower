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
