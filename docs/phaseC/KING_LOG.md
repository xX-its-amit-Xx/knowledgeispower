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
