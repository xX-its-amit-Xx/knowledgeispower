# scp-graphfs-mcp

## 1. Claim

> Curated protein-complex membership (CORUM, EBI Complex Portal) is a better prior for **reconstructing a cell's protein state from incomplete single-cell proteomics measurements** than data alone, than generic PPI networks (STRING), and than a degree-preserving shuffled version of the same graph.

**Why it could be true:** complex subunits are made in fixed ratios, and subunits that fail to assemble are degraded. So co-complex proteins co-vary at the protein level, more than their mRNAs do (Gonçalves et al., *Cell Systems* 2017, [10.1016/j.cels.2017.08.013](https://doi.org/10.1016/j.cels.2017.08.013)). That makes a missing subunit partly inferable from its measured partners.

## v9 validation (Nov 2026)

> **The latest result; it supersedes the Oct 22 update below.** It is a multi-holdout test of the reference-transfer
> champion (H15) against a 17-method panel, preregistered in `benchmarks/PREREGISTRATION_v9.md`. There were 3 new
> projects (5 holdout sets), 30-cell subsamples, MCAR and MNAR masks, and a 1% SESOI.
>
> | Hypothesis | Projects | Verdict | Strict reading (post hoc) |
> |---|---|---|---|
> | **H16:** the champion beats 5 comparators | 1 of 3 (karagach2025) | **NOT CONFIRMED** | 0 of 3 |
> | **H17:** the revision (`cshrink` target standardisation) beats the champion | 3 of 3 | **CONFIRMED** | 2 of 3 |
> | **H18:** the revision beats the strongest competitors | 1 of 3 (karagach2025) | **NOT CONFIRMED** | 1 of 3 |
>
> **What the tests show:**
> - **Strongest competitor: a data-only one.** It is the no-edge "GCN" (really a protein-wise masked MLP), stacked with
>   kNN and low-rank (`gcn_nograph+stack`).
>   - It ties the champion on average (RMSE ratio 0.998).
>   - The champion is significantly worse than it in 5 of 10 set × mask cells.
>   - The best graph-prior competitor, `gcn+stack` (STRING), averages 1.033.
> - **What the champion does beat** (every set and mask, all entries; every CI excludes 0):
>   - the kNN + low-rank stack, by 3.0–14.0%;
>   - its own permuted-reference null, by 2.7–4.1%;
>   - its no-reference ablation, by 2.6–8.3%.
> - **The revision** improves on the champion by 1.1–3.4% (MCAR) and 1.9–7.0% (MNAR), across both endpoints. It is a
>   data-only preprocessing change, not a biological prior.
>
> **BPCA (post hoc):** the frozen panel's BPCA was broken (it treated cells as observations). On derks2022 it equalled
> per-protein mean imputation exactly, as it did on 20 of 24 dev splits.
> - **derks2022:** a correctly oriented BPCA (`bpca_oba`) beats the champion (0.978 MCAR / 0.967 MNAR) and narrowly
>   beats (MCAR) or ties (MNAR) the revision, by less than the 1% SESOI. The pcaMethods BPCA, on 4 splits without CIs, gives 0.947 /
>   0.934 relative to the champion.
> - **karagach2025_S1** (newer files, see `docs/HANDOFF_v9.md` §3): the revision beats `bpca_oba+stack` by 4.1% /
>   5.4% (MCAR / MNAR). The champion only roughly ties it (−0.89% / +0.70%), and raw `bpca_oba` beats the champion
>   under MNAR (champion +1.43% worse, CI excluding 0).
> - The px061710 addenda and all four secondary holdouts are still running or queued (`docs/HANDOFF_v9.md` §2).
>
> **Graph priors did not help.**
> - In the GCN, STRING beat a degree-preserving shuffle but lost to no edges in all 10 holdout × mask cells.
> - The other STRING methods did no better than a generic data-only third stack component.
> - No v9 method used CORUM or Complex Portal, so v9 itself does not test the thesis.
> - On the whole record (v2–v9), the thesis in §1 should be reported as **not supported**. What has replicated is
>   transfer of protein-specific structure from public single-cell proteomics data, which is empirical, not curated.
>
> **Read more:**
> - [`docs/v9_validation.md`](docs/v9_validation.md): the authoritative report;
> - [`docs/v9_results.md`](docs/v9_results.md): the generated verdict tables;
> - [`docs/v9_addendum_results.md`](docs/v9_addendum_results.md): the post-hoc addenda;
> - [`docs/v9_review.md`](docs/v9_review.md): the independent review;
> - [`docs/HANDOFF_v9.md`](docs/HANDOFF_v9.md): what is still running and how to resume.

## 2. Headline result

> **Update, Oct 22 (final).** Our preregistered refined hypothesis **H14 was confirmed (6/6) on an untouched holdout**
> (Brunner 2022, HeLa, timsTOF DIA, 30-cell subsamples).
> - **What H14 tests:** an *empirical biological prior* (protein co-regulation learned from 16 public single-cell datasets
>   + CCLE + ProteomeHD) stacked with kNN and low-rank.
> - **Result on complex-member proteins:** 2.4–5.4% lower RMSE than the best data-only method (kNN + low-rank), 6–9% lower
>   than kNN, and 2.5–4.3% lower than the same pipeline with a scrambled reference.
> - **Caveats:**
>   - Curated complexes are the best of 22 knowledge priors, but they add nothing measurable on top of the empirical prior.
>   - On a second holdout (Khan 2023) the gain was below the 1% bar, so the result replicates on 1 of 2 holdouts.
> - Full record: `docs/EXPLORATION_LOG.md` (Amendments 9–14), `results/confirm_v7/`.
>
> The original phase-1 headline follows.

<!-- HEADLINE:START -->
Dataset: 215 cells × 2813 protein groups, 40.7% missing; 1545 proteins are in ≥1 curated complex. 5 seeds × 2 masks, 10% of observed values masked; 10 shuffled nulls per prior; 95% paired-bootstrap CIs. Run: `7626455`, 567 s.

**Imputation, proteins in a curated complex (RMSE, lower is better; Δ = prior − reference, 95% CI):**

| Prior | Mask | RMSE | kNN | Δ vs kNN | own null (mean) | Δ vs own null |
|---|---|---|---|---|---|---|
| Complex Portal | MCAR | 0.730 | 0.687 | +0.043 [+0.037, +0.048] | 0.807 | -0.076 [-0.082, -0.071] |
| CORUM | MCAR | 0.771 | 0.687 | +0.083 [+0.077, +0.089] | 0.843 | -0.073 [-0.078, -0.068] |
| STRING (combined ≥700) | MCAR | 0.761 | 0.687 | +0.074 [+0.068, +0.079] | 0.895 | -0.134 [-0.140, -0.128] |
| STRING (experimental ≥700) | MCAR | 0.784 | 0.687 | +0.097 [+0.091, +0.103] | 0.897 | -0.113 [-0.119, -0.108] |
| Complex Portal | MNAR | 0.816 | 0.772 | +0.044 [+0.040, +0.048] | 0.901 | -0.085 [-0.089, -0.081] |
| CORUM | MNAR | 0.865 | 0.772 | +0.093 [+0.088, +0.097] | 0.940 | -0.075 [-0.080, -0.071] |
| STRING (combined ≥700) | MNAR | 0.878 | 0.772 | +0.106 [+0.101, +0.110] | 1.002 | -0.124 [-0.128, -0.121] |
| STRING (experimental ≥700) | MNAR | 0.905 | 0.772 | +0.133 [+0.128, +0.137] | 1.004 | -0.099 [-0.104, -0.095] |

Not in a complex, the curated priors fall back to kNN, so they are identical to it (MCAR 0.796, MNAR 0.850). STRING has partners for these proteins and is worse there (MCAR 0.927, MNAR 1.034).

**Preregistered verdicts:**

- **Test 1 (headline): complex prior beats kNN *and* its own null, in-complex, on both masks:** **fail**.
  - Every prior beats its own shuffled null, so the curated graph structure carries real signal.
  - No prior beats kNN imputation, which uses whole-cell similarity, so "better than data alone" is not supported.
- **Test 2: panel design beats HVP and the null:** **fail**. Complex-propagated panels won 0/3 (CORUM) and 0/3 (Complex Portal) k values.
  - At k=50, mean held-out R² was 0.193 for Complex Portal, 0.181 for CORUM, 0.189 for HVP and 0.419 for **random panels**.
  - Random beat every scored panel at k=25 and k=50. The HVP seed scores favour high-variance proteins, including likely contaminants.
- **Test 3: stability (k=50, Jaccard):** **pass**, by a negligible margin.
  - Complex Portal 0.761 vs HVP 0.751 vs its null 0.760.
  - CORUM 0.742 (fail).
- **Test 4: curated ≥ STRING-exp ≥ STRING on in-complex RMSE:** **fail**. Complex Portal is best of the four priors on both masks, but STRING-combined beats STRING-experimental.
- **Sanity (classification, k=50, macro-F1, not evidence):** HVP 0.787, random 0.767, Complex Portal 0.768, all proteins 0.865. Ceiling flag: False.

**Caveat found after the run (not preregistered):** for well-measured complex proteins, the weighted mean of *shuffled-null* partners already correlates with the protein at median r = 0.63 (IQR 0.38–0.79), against median r = 0.75 for real curated partners (`results/discordant_subunits.csv`). Much of the co-variation is shared broadly across proteins, not specific to complexes.

**Exploratory, not preregistered: per-cell median normalisation** (`results/exploratory/normalise_median/`). The verdicts are unchanged (test 1 **fail**). In-complex MCAR RMSE: Complex Portal 0.725, kNN 0.669, own null 0.797. Null-partner median r stays 0.63, so per-cell scale does not explain kNN's advantage.
<!-- HEADLINE:END -->

The design was preregistered in [`benchmarks/PREREGISTRATION.md`](benchmarks/PREREGISTRATION.md), committed before any benchmark ran. Every row, including losses, is in `results/*.json` and on the dashboard.

## 2b. Phase 2–3: looking for a version of the claim that survives held-out tests

<!-- PHASE23:START -->
*Exploration protocol:* [`benchmarks/EXPLORATION.md`](benchmarks/EXPLORATION.md). Every dev result is in [`docs/EXPLORATION_LOG.md`](docs/EXPLORATION_LOG.md). Each confirmation was preregistered before its held-out data was touched, and run once.

**v2: Petrosius held-out masks.** Primary E3 (low-rank + CORUM Laplacian): **not confirmed** (6/8). It does not beat its own shuffled null:
- MCAR: -0.0007 [-0.0031, +0.0016]
- MNAR: +0.0001 [-0.0018, +0.0017]

**v3: independent dataset (Leduc 2022, 1543 cells); E6, Complex Portal residuals on low-rank.** Passes 8/8 as written, but the gain is practically zero (no minimum effect size was preregistered; disclosed):
- MCAR vs lowrank_nograph: -0.0000 [-0.0001, -0.0000] (-0.01%)
- MCAR vs own_null_mean: -0.0000 [-0.0001, -0.0000] (-0.01%)
- MCAR vs data_graph: -0.0007 [-0.0008, -0.0005] (-0.14%)
- MCAR vs knn: -0.0047 [-0.0053, -0.0039] (-1.00%)
- MNAR vs lowrank_nograph: -0.0001 [-0.0001, -0.0001] (-0.01%)
- MNAR vs own_null_mean: -0.0001 [-0.0001, -0.0001] (-0.01%)
- MNAR vs data_graph: -0.0004 [-0.0005, -0.0002] (-0.07%)
- MNAR vs knn: -0.0055 [-0.0062, -0.0048] (-0.98%)

**v4: Leduc, 30-cell subsamples; data graph + Complex Portal.** **Not confirmed** (2/6). It loses to kNN. Secondary: Complex Portal partners vs. data-derived partners:
- n=30 MCAR: -0.0298 [-0.0340, -0.0258] (-5.33%)
- n=30 MNAR: -0.0177 [-0.0219, -0.0138] (-2.77%)
- n=60 MCAR: -0.0253 [-0.0284, -0.0225] (-4.73%)
- n=60 MNAR: -0.0230 [-0.0257, -0.0202] (-3.72%)

**v5: third dataset (SCoPE2), 30-cell subsamples, minimum effect size 1%.** The claim is narrow: *curated Complex Portal partners beat data-derived partners and a shuffled null.* Result: **not confirmed** (0/4).
- n=30 MCAR vs data_alone: -0.0043 [-0.0104, +0.0016] (-0.78%) ✗
- n=30 MCAR vs own_null_mean: -0.0021 [-0.0050, +0.0009] (-0.38%) ✗
- n=30 MCAR vs knn (secondary): +0.0057 [+0.0029, +0.0085] (+1.04%) ✗
- n=30 MNAR vs data_alone: -0.0063 [-0.0121, -0.0004] (-0.95%) ✗
- n=30 MNAR vs own_null_mean: -0.0014 [-0.0041, +0.0013] (-0.21%) ✗
- n=30 MNAR vs knn (secondary): -0.0002 [-0.0028, +0.0026] (-0.02%) ✗
- n=60 MCAR vs data_alone (secondary): -0.0078 [-0.0119, -0.0038] (-1.45%) ✓
- n=60 MCAR vs own_null_mean (secondary): -0.0032 [-0.0046, -0.0018] (-0.60%) ✗
- n=60 MCAR vs knn (secondary): +0.0010 [-0.0004, +0.0024] (+0.19%) ✗
- n=60 MNAR vs data_alone (secondary): -0.0094 [-0.0136, -0.0049] (-1.48%) ✓
- n=60 MNAR vs own_null_mean (secondary): -0.0029 [-0.0046, -0.0012] (-0.46%) ✗
- n=60 MNAR vs knn (secondary): -0.0026 [-0.0043, -0.0009] (-0.42%) ✗

**Bottom line after three datasets (no hand-typed numbers above; see the JSON):**
- **Robust:** curated complexes beat degree-preserving shuffles of themselves.
- **Small, and not confirmed at a 1% minimum effect:** Complex Portal partners beat data-derived partners.
- **Not supported:** "better than data alone" as an imputer. kNN or low-rank models match or beat every prior-based method by any meaningful margin.
<!-- PHASE23:END -->

## 3. Quickstart

```bash
git clone https://github.com/xX-its-amit-Xx/knowledgeispower && cd knowledgeispower
python -m venv .venv && source .venv/bin/activate && pip install -e ".[dev]"
pytest -q -m "not network"
python benchmarks/run.py                  # ~10 min on a laptop -> results/
scp-graphfs-mcp                           # MCP server over stdio
claude mcp add scp-graphfs -- "$PWD/.venv/bin/scp-graphfs-mcp"   # connect Claude Code
```

The MCP tools are `get_complexes`, `explain_protein`, `build_prior`, `propagate_scores`, `impute_with_prior` and `select_panel`. See a real session in [`demo/agent_session.md`](demo/agent_session.md).

## 4. Prior art and differentiation

Every DOI and arXiv ID resolved, and its title, authors, year and venue matched (re-checked by an independent verifier). Some specific claims about papers could not be checked from the full text, and those entries say so in `unverified_note`. GIRAFI's use of a STRING v11.5 prior could not be re-checked (bioRxiv was rate-limited). The full list (35 entries) and query log in [`docs/prior-art.json`](docs/prior-art.json) and [`docs/prior-art-search-log.md`](docs/prior-art-search-log.md).

| Work | Class | How we differ |
|---|---|---|
| **GIRAFI** (bioRxiv 2026, [10.64898/2026.05.18.726081](https://doi.org/10.64898/2026.05.18.726081)) | **Direct overlap** | Graph-prior imputation for SCP **already exists**. It uses a STRING prior, and in the text we fetched we found no shuffled-network control. We test *which* prior (curated complexes vs. STRING vs. a degree-preserving null) |
| Petrosius et al., *MCP* 2025 ([10.1016/j.mcpro.2025.100982](https://doi.org/10.1016/j.mcpro.2025.100982)) | Direct overlap (mechanism) | Shows known complexes co-vary in Astral SCP data, descriptively. We use that co-variation as a prior and benchmark it |
| Fuzzy-FishNET, Goh, *BMC Med Genomics* 2016 ([10.1186/s12920-016-0228-z](https://doi.org/10.1186/s12920-016-0228-z)) | Direct overlap | Complex-based feature selection, but on bulk proteomics for class discrimination |
| Protrec2, *Brief Bioinform* 2025 ([10.1093/bib/bbaf692](https://doi.org/10.1093/bib/bbaf692)) | Direct overlap | Complex membership to recover missing proteins, but in bulk data and predicting detection, not single-cell abundance |
| NIFty (bioRxiv 2026, [10.64898/2026.03.06.710179](https://doi.org/10.64898/2026.03.06.710179)) | Direct overlap (task) | Missingness-tolerant SCP feature selection with no prior. Candidate baseline (slot S4) |
| GEARS, *Nat Biotechnol* 2024 ([10.1038/s41587-023-01905-6](https://doi.org/10.1038/s41587-023-01905-6)) | Adjacent | Transcriptomic perturbation prediction with GO and co-expression graphs; no complex prior |
| Wefers et al., *Nat Methods* 2026 ([10.1038/s41592-026-03142-6](https://doi.org/10.1038/s41592-026-03142-6)) | Adjacent | A localisation-prediction benchmark where generic PPI (STRING/BioPlex) GraphSAGE gave no substantial gain. Not SCP |
| AROMA (arXiv [2604.20263](https://arxiv.org/abs/2604.20263)) | Adjacent | Merges STRING and CORUM into one graph for perturbation reasoning. We keep the sources separate and test each against a null |
| PRINCE, *PLoS Comput Biol* 2010 ([10.1371/journal.pcbi.1000641](https://doi.org/10.1371/journal.pcbi.1000641)) | Building block | Network propagation; our RWR panel scoring follows it |
| STRING MCP ([string-db.org/help/mcp](https://string-db.org/help/mcp/)), ToolUniverse (arXiv [2509.23426](https://arxiv.org/abs/2509.23426)) | Building block | Agent access to STRING and to Complex Portal lookups already exists. A bare complex-lookup MCP server is not novel |
| MCPmed (arXiv [2507.08055](https://arxiv.org/abs/2507.08055)) | Collaborator | MCP wrappers for GEO, STRING and the UCSC Cell Browser. Our server could join |
| CORUM ([10.1093/nar/gkae1033](https://doi.org/10.1093/nar/gkae1033)), Complex Portal ([10.1093/nar/gkae1085](https://doi.org/10.1093/nar/gkae1085)), STRING ([10.1093/nar/gkae1113](https://doi.org/10.1093/nar/gkae1113)) | Data sources | Snapshots, versions and SHA-256 in [`data/PROVENANCE.md`](data/PROVENANCE.md) |

**Corrected from our starting notes:**
- Wefers 2026 is not an SCP paper.
- The Leduc 2022 nPOP paper does not analyse complex co-variation.
- The complex-attenuation paper's first author is Gonçalves.
- AROMA and MAP-KG are unrelated works.

## 5. What's novel, narrowly

Four things together:
- curated complex priors (CORUM and Complex Portal, each tested separately against STRING);
- the single-cell proteomics reconstruction tasks: imputation under MCAR and MNAR masks, and panel design;
- a degree-preserving shuffled-graph null for every prior;
- an MCP interface so an LLM agent can run it.

None of the pieces is new alone.

## 6. Limitations

- **Coverage:** curated databases have limited coverage and favour well-studied proteins.
  - Only 1317 of 2813 protein groups here have a curated-complex partner that is also measured (CORUM 1147, Complex Portal 847).
  - CORUM is mammal-centric.
- **One dataset:** 215 cells of one AML model (Petrosius et al. 2023 preprint, [10.1101/2023.06.06.543943](https://doi.org/10.1101/2023.06.06.543943)).
- **Masking undercounts the hardest cases:** masked *observed* values under-represent the hardest true dropouts. Our MNAR mask only mimics them.
- **`layers['raw']` is reconstructed** as `2**X`, because the source ships log2 values only.
- **Contaminants:** the HVP panel scores pick up keratins (likely sample contaminants), as the agent transcript shows. No contaminant filter was preregistered.
- **The CORUM null keeps ~30% of real edges.** Its large cliques restrict how far degree-preserving swaps can move, so that null is conservative.
- **Optimistic CIs:** the bootstrap resamples masked entries independently and pools overlapping masks from 5 seeds, so the CIs are narrower than they should be.
- **Hackathon scale:** built in about one day.

## 7. Contributors

- Thesis, priors, benchmark, core MCP: **Amit Shenoy**
- Hard-filter arm: **Collaborator A** *(in progress; credited here when merged)*
- Variant→complex bridge (slot S1): *in progress; credited here once merged*
- One line per merged slot: see [`CONTRIBUTORS.md`](CONTRIBUTORS.md)

## 8. Where you could plug in

Open slots are listed in [`docs/SLOTS.md`](docs/SLOTS.md); start with [`docs/ONBOARDING.md`](docs/ONBOARDING.md) (15 minutes). Coding agents: read [`AGENTS.md`](AGENTS.md).

| Slot | Task |
|---|---|
| S2 | Second dataset loader + cross-context test |
| S3 | Extra prior source (hu.MAP, BioPlex, ontology) |
| S4 | NIFty baseline wrapper |
| S5 | MCP client recipes |
| S6 | Glossary + dashboard explainers |
| S7 | Test hardening |

## 9. Citation and license

See [`CITATION.cff`](CITATION.cff); also cite the data sources in [`data/PROVENANCE.md`](data/PROVENANCE.md). Code is under the Apache-2.0 license. The `scp-prep` dependency is GPL-3.0-or-later, so anyone redistributing a combined build is bound by GPL terms.
