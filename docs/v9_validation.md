# v9 validation: what was tested, what held, and what it means

*Written 2026-11-08 for Amit and for a skeptical outside reader.*

**Status of the inputs:**
- The five primary holdout results are final (`results/confirm_v9/*.json`).
- No secondary holdout has finished. The run log `results/v9/holdout_runs.log` shows khan2023 starting at
  2026-11-08T06:16:41Z, and no secondary result file exists yet.
- The post-hoc addendum has finished for derks2022 and six dev targets only (`results/v9/addendum/`).
- The second post-hoc addendum (Addendum 2, 4 seeds, no CIs) has finished for derks2022 only
  (`results/v9/addendum2/derks2022.json`, written 2026-11-08T06:51:04Z per `results/v9/addendum2_runs.log`);
  karagach2025_S1 started at the same time.

**How to read the numbers:**
- Every number below is copied from the file named next to it, or the nearest named file above it.
- Percentages are ΔRMSE as a percentage of the comparator's RMSE. A negative value means the first-named method has
  *lower* error.
- "✗" marks a test that fails the preregistered rule.

---

## 1. Summary

**What was asked** (Amit, 2026-11-07, recorded in `benchmarks/PREREGISTRATION_v9.md` §0):
1. Redo the earlier validation.
2. Add at least two more holdouts.
3. Compare against the methods researchers typically use, including graph-prior methods not tried before.
4. Answer the underlying question: does our method add meaningful value?

**What was done:**
- The v6 khan2023 test was rerun with the frozen v6 code. It reproduced exactly (§2).
- A comparator panel was built: 10 data-only methods, 4 STRING graph-prior methods (`grmf_string`, `rwr_resid`, `gcn`,
  `imc_embed`), a data-graph control (`grmf_datagraph`) and 2 GCN graph controls (`results/v9/v9_selection.json`),
  each scored raw and stacked. The panel was screened on dev, and the two strongest competitors were frozen before any
  holdout was touched (§4).
- Three hypotheses were preregistered and tested on 3 new primary projects (5 holdout sets):
  - **H16:** the current champion beats 5 comparators.
  - **H17:** a revision (`revised`, the Phase-C `cshrink`) beats the champion.
  - **H18:** the revision beats the strongest competitors.
- A 39-agent independent review, a 23-agent expansion with four red-team reports, and an independent verdict audit
  followed (§6).

**Preregistered verdicts** (rule in §4 of the preregistration; counts from `docs/v9_results.md`, generated from
`results/confirm_v9/*.json`):

| Hypothesis, as preregistered | Projects confirming | Verdict | Strict-SESOI reading (post hoc) |
|---|---|---|---|
| **H16:** the champion beats `stack_lowrank`, `champion@null`, `champion@noref`, `gcn_nograph+stack` and `gcn+stack` on all 20 tests per holdout | 1 of 3 (karagach2025 only) | **NOT CONFIRMED** | 0 of 3, not confirmed |
| **H17:** the revision beats the champion on all 4 tests per holdout | 3 of 3 | **CONFIRMED** | 2 of 3 (derks2022 drops out), still confirmed |
| **H18:** the revision beats `stack_lowrank`, `gcn_nograph+stack` and `gcn+stack` on all 12 tests per holdout | 1 of 3 (karagach2025 only) | **NOT CONFIRMED** | 1 of 3, not confirmed |

The preregistration applies the 1% SESOI to the point estimate. The strict reading requires the CI's upper bound itself
to be ≤ −1%. The strict project-level verdicts come from the independent audit in `results/v9/expansion_raw.json`.

**Answer to "does our method add meaningful value?"**

- **By the preregistered test, no.** H16 is not confirmed: 1 of 3 projects, and 0 of 3 under the strict reading. The
  preregistration defines H16 as "the champion generalises and adds meaningful value" and says "Any failure will be
  reported as such, not reinterpreted" (`benchmarks/PREREGISTRATION_v9.md` §4).
- **What does hold is a subset of H16's tests** (a post-hoc reading: the comparators of v8's H15, "not the v9 verdict",
  in the words of red-team report 2, entry 20 of the list in `results/v9/expansion_raw.json`). On all 5 sets, under
  both masks, the H15 champion (kNN + low-rank + reference-adapted projection) had lower all-entry RMSE than:
  - the kNN + low-rank stack, the strongest data-only baseline of v6–v8: by 3.0–14.0%;
  - the same model with the reference loadings permuted across proteins: by 2.7–4.1%;
  - the same model with no reference: by 2.6–8.3%.

  Every one of those CIs excludes 0 (`results/confirm_v9/*.json`). The last two are internal ablations, not comparisons
  with independent data-only methods. `champion@noref` also reuses the champion's reference-tuned settings (λ_v = 100,
  10 extra factors, per-protein unit-variance scaling), so it is a weak stand-in for a data-only method (red-team
  report 4, entry 22).
- **The champion does not beat the strongest data-only method in our own panel.** That method is a no-edge "GCN", which
  is really a protein-wise masked MLP, stacked with kNN and low-rank (`gcn_nograph+stack`). Averaged over masks and sets,
  its RMSE relative to the champion is 0.998 (`docs/v9_results.md`). On all entries the champion was significantly
  *worse* than it in 5 of 10 set × mask cells, including all three px061710 sets under MCAR (§5.2).
- **Post hoc, a correctly oriented BPCA beats the champion on derks2022 under both masks, and narrowly beats (MCAR) or
  ties (MNAR) the revision.** Its edge over the revision is below the 1% SESOI in every mask × endpoint cell (§7). The
  frozen panel's BPCA was broken (§6). BPCA has so far been scored on no other holdout. A second post-hoc run with
  the pcaMethods reference BPCA (Addendum 2: 4 of the 10 splits, no CIs, revision not scored) points the same way
  against the champion: 0.947 (MCAR) and 0.934 (MNAR) relative to it (`results/v9/addendum2/derks2022.json`).
- **The revision improves on the champion in every set and mask (H17 confirmed):** by 1.1–3.4% under MCAR and 1.9–7.0%
  under MNAR, across both endpoints. The direction is consistent on all 5 sets, but the size of the gain in a new project
  is not pinned down (project-level MCAR 95% CI [−4.75, +0.20], §5.4). It does not reliably beat the no-edge MLP under
  MCAR. It is a data-only change to how the target is standardised, not a biological prior.
- **Graph priors did not help.** In the 2-layer GCN, STRING edges beat a degree-preserving shuffle but lost to no edges
  in all 10 holdout × mask cells.
- **The project's original thesis was not tested in v9.** No v9 method uses CORUM or Complex Portal. Earlier holdouts
  showed curated complex partners beating data-derived partners in some tests, but a curated prior never beat a strong
  data-only or reference-transfer imputer by a meaningful margin (§8).

---

## 2. Redo of the earlier validation

### 2.1 Exact reproduction of the v6 khan2023 test

The v6 confirmation was rerun at its original commit `411d439`, with its 17-dataset reference
(10,129 × 4096 reference matrix, 84% warm proteins) and seeds 6000–6009.

The result is `results/v9/v6_khan2023_reproduction.json`; the original is `results/confirm_v6/confirm_v6.json`. The two
files differ in a single field, `runtime_seconds` (87.9 vs 98.6). Every RMSE, every CI, every stack weight and the
verdict are identical: **3/8 passes, not confirmed.**

| `stack3_big` vs | MCAR | MNAR |
|---|---|---|
| kNN | −1.40% (pass) | −0.98% (fail: SESOI) |
| low-rank | −28.28% (pass) | −22.91% (pass) |
| kNN + low-rank stack | −0.55% (fail: SESOI) | −0.20% (fail: SESOI) |
| own permuted-loading null (mean of 2) | −0.61% (fail: SESOI) | −0.20% (fail: SESOI) |

Every CI in this table lies below 0. The failures are all failures of the 1% SESOI, not of significance.

**Naming note.** `PREREGISTRATION_v9.md` calls the khan2023 test "H14". In `benchmarks/PREREGISTRATION_v6.md` and
`docs/EXPLORATION_LOG.md` it is **H12**. H14 is the brunner2022 test (v7). The method and the result are the same under
either label.

### 2.2 What the prior holdout record means

All entries are from `docs/EXPLORATION_LOG.md`.

| Test | Holdout | Method | Verdict | Size of the reference effect |
|---|---|---|---|---|
| v6, H12 | khan2023 (MCF-10A) | `stack3_big` (plain projection, 17-dataset reference plus bulk) | not confirmed, 3/8 | 0.2–0.6% over the stack and the null |
| v7, H14 | brunner2022 (HeLa, timsTOF DIA) | `stack3_plain` | confirmed, 6/6 | in-complex −2.44% / −5.42% vs the stack, −2.51% / −4.30% vs the null (MCAR / MNAR) |
| v8, H15 | wu2026 (prenatal brain, DIA, 90% missing) | `adapt` champion | confirmed, 16/16 | −5.0% / −6.8% vs kNN + low-rank, −5.2% / −6.1% vs the scrambled reference (all entries, MCAR / MNAR) |

**What the reproduction shows:**
- The khan2023 failure is not a pipeline or nondeterminism artefact.
- It was a real, small effect that fell below the SESOI.

**What the record shows overall:**
- The three tests used three successively refined methods.
- H14 was refined after seeing dev results, and this was disclosed in `PREREGISTRATION_v7.md`.
- So before v9, the current champion (H15) had been tested on exactly one holdout.
- Reference transfer had replicated on 2 of 3 holdouts, with effect sizes that varied by dataset.

v9 is therefore the first multi-holdout test of the H15 champion.

---

## 3. Holdouts

### 3.1 Primary holdouts

All rows are from `docs/holdout_v9_datasets.md`, except derks2022, which is from `data/PROVENANCE.md`.

| Set | Project | Acquisition | Cells × proteins | Missing | Biological content |
|---|---|---|---|---|---|
| derks2022 | scpdata `derks2022` (Derks et al. 2022) | plexDIA, mTRAQ, Q-Exactive | 125 × 1093 | 0.566 | PDAC 45, melanoma 43, U-937 37 |
| karagach2025_S1 | PXD058457 (Karagach et al. 2025) | timsTOF SCP, Slice-PASEF DIA, DIA-NN 1.9.1 | 258 × 3687 | 0.334 | HeLa only; 19 sample-prep variants (technical) |
| px061710_GSTP1KO | PXD061710 | timsTOF SCP, diaPASEF, Spectronaut | 78 × 3899 | 0.425 | RKO sgCTRL 39, sgGSTP1 39 |
| px061710_MDAfix | PXD061710 | same | 76 × 3325 | 0.377 | MDA fresh, plus three PFA levels (19 each) |
| px061710_RKOfix | PXD061710 | same | 75 × 5088 | 0.411 | RKO fresh 18, plus three PFA levels (19 each) |

**Provenance:**
- The PRIDE projects were confirmed through the PRIDE API on 2026-11-07. All are CC0. File sizes match PRIDE, and
  SHA-256 was computed locally because PRIDE gives no checksums.
- derks2022 was rebuilt from the unimputed precursor assay (`qe_prec_extracted`). The scpdata `proteins` assay is
  imputed, so it was not used.
- derks2022 had been reserved as a holdout since Phase B.
- No QC beyond counts was computed before the preregistration.

**Why px061710 counts as one project.** The three px061710 sets come from one PRIDE project, one lab (PI A. von
Kriegsheim) and one instrument. GSTP1KO and RKOfix are both RKO cells. The preregistration therefore treats them as one
project, which confirms a hypothesis if ≥ 2 of its 3 sets confirm. So the verdicts rest on 3 independent units, not 5.

The red team notes a side effect of this rule (`results/v9/expansion_raw.json`). Suppose each set passes independently
with probability q. Then "≥ 2 of 3" passes with probability 3q² − 2q³. That is higher than q when q > 0.5 and lower when
q < 0.5, so px061710 is not treated exactly like a single-set project.

### 3.2 Secondary holdouts

These are reported only, and none has finished.

| Set | Why it is secondary |
|---|---|
| khan2023, brunner2022 | Already used for the H12 and H14 tests. |
| karagach2025_S2 (39 cells) | Fewer than 60 cells, and the same lab and cell line as S1. |
| px077946_H2O2 (40 cells) | Fewer than 60 cells. |

### 3.3 Design (preregistration §3)

- Seeds 9000–9009.
- 30-cell random subsamples.
- 10% MCAR and MNAR masks, on proteins observed in ≥ 50% of the subsample.
- Two endpoints: all masked entries, and complex-member proteins.
- 95% paired entry bootstrap, 1000 resamples, seed 17.
- Pass rule: CI entirely below 0 **and** ΔRMSE ≤ −1% of the comparator's RMSE.

**Reference.** The champion's reference is 22 public single-cell datasets: a 12,008-row reference matrix aligned to each
holdout. The fraction of holdout proteins with reference loadings is:

| Set | Warm fraction |
|---|---|
| derks2022 | 0.994 |
| karagach2025_S1 | 0.983 |
| px061710_GSTP1KO | 0.981 |
| px061710_MDAfix | 0.978 |
| px061710_RKOfix | 0.942 |

**Leakage.** The independent review found no path from holdout values into the reference, the tuning or the selection
(`docs/v9_review.md`).

Two overlaps are known, and the check is incomplete:
- derks2022's reference contains leduc2022 and specht2021. They come from the same lab and share cell lines, though
  they are different runs.
- karagach2025_S1 (HeLa) shares its cell line with liang2020_hela and cong2020AC in the reference. karagach2025 is also
  the only project that confirms H16; red-team report 2 (`results/v9/expansion_raw.json`, entry 20) flags this as a
  possible confound (an unverified hypothesis that three projects cannot resolve).
- Whether the px061710 cell lines (RKO, MDA-MB-231) appear in the reference "was not checked" (red-team report 3,
  entry 21).

---

## 4. Methods compared

### 4.1 The frozen panel

From `PREREGISTRATION_v9.md` §0 and `src/scp_graphfs_mcp/baselines.py`:

| Group | Methods |
|---|---|
| Data-only | mean; MinDet; MinProb; kNN; softImpute; iterative SVD; BPCA; MICE; missForest (10 trees, 3 iterations); denoising autoencoder |
| Graph prior (STRING v12.5, score ≥ 700) | Laplacian-regularised MF (`grmf_string`); the same with an equal-size data-derived graph (`grmf_datagraph`, control); random walk with restart on kNN residuals (`rwr_resid`); 2-layer GCN (`gcn`); inductive MF on spectral graph embeddings (`imc_embed`) |
| GCN graph controls | `gcn_null` (degree-preserving shuffle of STRING); `gcn_nograph` (no edges) |

- Each method was scored raw and as `+stack`. The stack is a simplex stack with kNN and low-rank, which is the same help
  the champion gets.
- PPCA was dropped because it diverges at 30 cells.

### 4.2 The GCN graph controls

All three GCN variants have the same architecture and differ only in the graph.

**Dev screen** (mean RMSE relative to the champion; `results/v9/v9_selection.json`):

| Variant | Raw | Stacked |
|---|---|---|
| `gcn` (STRING) | 1.0415 | 1.0143 |
| `gcn_null` (shuffled STRING) | 1.0959 | 1.0487 |
| `gcn_nograph` (no edges) | 0.9929 | 0.9917 |

- No edges beat STRING on 11 of 12 dev target × mask cells (`PREREGISTRATION_v9.md` says "splits"). Per split, raw
  `gcn_nograph` beats raw `gcn` on 20 of 24 and the stacked versions on 21 of 24 (computed from
  `results/v9/dev_panel_devctl_fast.json`).
- The no-edge network is a protein-wise masked MLP across the subsample's cells. It uses no prior and no reference, so
  it is a data-only method.

### 4.3 The dev selection and its amendment

**The rule** (frozen in `results/v9/v9_selection.json`): pick the lowest mean dev RMSE (all entries) relative to the
champion, separately among graph-prior and data-only methods, on the Phase-C fast pool (6 targets × 2 masks × 2 seeds).

**Timeline** (`git log` of `benchmarks/PREREGISTRATION_v9.md`, and `results/v9/holdout_runs.log`):

| Time (UTC) | Commit | Event |
|---|---|---|
| 2026-11-07 23:06:58 | `a794dab` | Draft. Selects `mice+stack` (1.0275) as best data-only and `gcn+stack` (1.0143) as best graph prior. Only H16 and H17 exist. |
| 23:23:33 | `ee28470` | **Amendment.** The GCN controls showed that no edges beat STRING, so `gcn_nograph+stack` (0.9917) replaces `mice+stack`, and H18 is added. |
| 23:34:37 | `0062e02` | Final preregistration, with the dev dry-run mechanics check. |
| 23:34:53 | (run) | First holdout (derks2022) starts. |

**The amendment made the tests harder, not easier.**
- It replaced a comparator that the champion beats on dev (1.0275) with one that beats the champion on dev (0.9917).
- The preregistration then predicted in writing that H16 would fail against `gcn_nograph+stack`. That prediction held.
  It also predicted that "the claim most likely to survive is H18". H18 failed (1 of 3 projects).
- The draft had no H18. Commit `a794dab` lists "`revised` against `mice+stack` and `gcn+stack`" only under "Reported …
  without a verdict" (`git show a794dab:benchmarks/PREREGISTRATION_v9.md`; `docs/v9_review.md` makes the same slip).
- On derks2022, the revision's RMSE ratio to the draft comparator `mice+stack` is 0.975 (MCAR) and 0.957 (MNAR),
  computed from `results/confirm_v9/derks2022.json`. Had H18 been written with the draft comparator `mice+stack`,
  derks2022 would very likely have passed it. Across all 5 sets × 2 masks × 2 endpoints the ratio is ≤ 0.976 (computed
  from `rmse` and `rmse_in_complex` in `results/confirm_v9/*.json`), and the revision passes every test against
  `stack_lowrank` and `gcn+stack`. So such an H18 would most likely have been confirmed; no CI against `mice+stack` was
  computed for the revision, so this is not certain.

**Three textual inconsistencies,** none affecting any number (`results/v9/review_raw.json`, verdicts 22 and 23):
- The table still files the no-edge GCN under "Controls".
- §4 of the preregistration calls `stack_lowrank` "the best data-only stack", a label carried over from v8. Under v9's
  own selection it is not: 1.0773, against 1.0275 for `mice+stack` and 0.9917 for `gcn_nograph+stack`
  (`results/v9/v9_selection.json`).
- The prediction that the revision beats `gcn_nograph+stack` "on every dev split" is overstated. The true count is 22 of
  24 splits, or 11 of 12 target × mask cells (`docs/v9_review.md`).

### 4.4 Champion and revision

**`champion` (H15):** stack(kNN, low-rank, `adapt`).
- `adapt` is MAP low-rank on the 30 target cells. It uses reference loadings as a Gaussian prior (λ_v = 100) plus 10
  free factors.
- The loadings come from ALS (rank 20, λ = 5, 15 iterations) on the 22-dataset single-cell reference.
- The stack weights are chosen on a 10% inner mask of observed entries.

**Champion controls:**
- `champion@null`: the loading rows permuted across proteins. One permutation (seed 0) is reused for all splits.
- `champion@noref`: no reference, with the same free factors.

**`revised`:** stack(kNN, low-rank, `adapt`, `cshrink`). `cshrink` shrinks the target's per-protein mean and log-sd
toward a constant inside the adapted projection, with ψ selected on an inner mask (`src/scp_graphfs_mcp/phasec/cshrink.py`).

**Dev history of `cshrink`** (`docs/EXPLORATION_LOG.md`, Phase C rounds 4–5, RMSE relative to the champion):

| Pool | Setting | Score |
|---|---|---|
| Full pool | selected ψ | 0.9609 (MCAR 0.9714, MNAR 0.9504) |
| Full pool | random-ψ null | 0.9716 |
| Full pool | fixed ψ = ∞ | 0.9607 |
| Fast pool | fixed ψ = 16 | 0.9576 |

- The log attributes the gain to dropping the champion's per-protein unit-variance scaling, not to tuning ψ.
- `cshrink` **missed its own Phase-C promotion bar**: full-pool MCAR was 0.9714 against a bar of 0.97, and woo2025_PC3
  MCAR was 1.0006 against a limit of 1.00. The v9 preregistration does not say this (§6).
- The revision gets 4 stack slots, while every comparator gets 3.

---

## 5. Results (preregistered; `results/confirm_v9/*.json`)

The tests and verdicts in §5.1 and §5.2 are the preregistered ones. Several readings in this section are post hoc and are
labelled where they appear: the mask-specific verdicts (§5.3), the SESOI sweep and the project-level intervals (§5.4),
and the CI-widening check (§5.1, §6). The panel CIs reported without a verdict (§5.2, §5.5) are unadjusted for
multiplicity: about 780 panel CIs, plus the GCN-control and addendum CIs, with no error control (red-team report 2,
`results/v9/expansion_raw.json`, entry 20). Read them as estimates, not as confirmed findings.

### 5.1 Tests passed per holdout

| Set | H16 (of 20) | H17 (of 4) | H18 (of 12) | H16 strict | H17 strict | H18 strict |
|---|---|---|---|---|---|---|
| derks2022 | 12 | 4 | 11 | 12 | 2 | 8 |
| karagach2025_S1 | **20** | 4 | **12** | 19 | 4 | 12 |
| px061710_GSTP1KO | 18 | 4 | 12 | 18 | 4 | 11 |
| px061710_MDAfix | 16 | 4 | 11 | 16 | 4 | 10 |
| px061710_RKOfix | 17 | 4 | 10 | 17 | 4 | 10 |

Project verdicts:
- **H16:** derks2022 ✗, karagach2025 ✓, px061710 ✗ (0 of 3 sets).
- **H17:** ✓ on all three projects.
- **H18:** derks2022 ✗, karagach2025 ✓, px061710 ✗ (1 of 3 sets: GSTP1KO). GSTP1KO's set-level pass rests on one
  test at −1.033% with CI upper bound −0.318% (MCAR, all entries); it flips if the CIs are widened 1.5× (post hoc, §6).

An independent audit recomputed all 330 reported values and found 0 mismatches with `docs/v9_results.md`
(`results/v9/expansion_raw.json`).

**The counts are not graded evidence.** Each H16 or H18 holdout verdict is the minimum over its comparators. Every H16
and H18 failure is against `gcn_nograph+stack`, plus `gcn+stack` in derks2022's H16. Every test against
`stack_lowrank`, `champion@null` and `champion@noref` passes, by at least 2.6% (smallest: derks2022 MCAR in-complex vs
`@noref`, −2.600%). The passing tests against `gcn+stack` have smaller floors: in H16, at least 1.98% (RKOfix MNAR, all
entries, −1.977%); in H18, at least 1.25% (derks2022 MCAR, all entries, −1.245%; in-complex −1.296%), with the largest
margin −11.38% (`results/confirm_v9/*.json`).

### 5.2 Key contrasts, all masked entries

Negative means the first-named method has lower RMSE.

**MCAR:**

| Set | champ vs `stack_lowrank` | champ vs `@null` | champ vs `@noref` | champ vs `gcn_nograph+stack` | champ vs `gcn+stack` | revised vs champ | revised vs `stack_lowrank` | revised vs `gcn_nograph+stack` | revised vs `gcn+stack` |
|---|---|---|---|---|---|---|---|---|---|
| derks2022 | −3.02 | −2.73 | −2.64 | +0.25 ✗ | −0.01 ✗ | −1.23 | −4.21 | −0.98 ✗ | −1.24 |
| karagach2025_S1 | −7.66 | −4.00 | −4.79 | −1.16 | −3.12 | −3.22 | −10.63 | −4.34 | −6.24 |
| px061710_GSTP1KO | −9.87 | −3.62 | −4.60 | +1.13 ✗ | −4.53 | −2.14 | −11.80 | −1.03 | −6.58 |
| px061710_MDAfix | −7.27 | −3.32 | −2.80 | +2.62 ✗ | −4.15 | −3.07 | −10.12 | −0.54 ✗ | −7.09 |
| px061710_RKOfix | −8.61 | −2.80 | −3.83 | +1.70 ✗ | −2.52 | −1.90 | −10.34 | −0.23 ✗ | −4.37 |

**MNAR:**

| Set | champ vs `stack_lowrank` | champ vs `@null` | champ vs `@noref` | champ vs `gcn_nograph+stack` | champ vs `gcn+stack` | revised vs champ | revised vs `stack_lowrank` | revised vs `gcn_nograph+stack` | revised vs `gcn+stack` |
|---|---|---|---|---|---|---|---|---|---|
| derks2022 | −4.43 | −4.09 | −4.34 | +0.63 ✗ | −0.50 ✗ | −3.03 | −7.33 | −2.42 | −3.52 |
| karagach2025_S1 | −11.10 | −4.12 | −7.20 | −2.41 | −4.05 | −6.02 | −16.46 | −8.28 | −9.83 |
| px061710_GSTP1KO | −13.96 | −4.06 | −8.30 | −1.55 | −5.17 | −4.70 | −18.01 | −6.18 | −9.63 |
| px061710_MDAfix | −12.66 | −3.80 | −6.86 | +2.30 ✗ | −5.44 | −6.07 | −17.96 | −3.91 | −11.18 |
| px061710_RKOfix | −7.32 | −2.93 | −5.13 | −0.78 ✗ | −1.98 | −3.21 | −10.30 | −3.97 | −5.13 |

On all entries, the champion's CI against `gcn_nograph+stack` lies entirely above 0 (champion significantly worse) in 5
of 10 set × mask cells: GSTP1KO, MDAfix and RKOfix under MCAR, and derks2022 (+0.63%, CI [+0.005, +1.33]) and MDAfix
(+2.30%, CI [+1.81, +2.86]) under MNAR (CIs expressed as % of the comparator's RMSE).

**In-complex endpoint.** The pattern is the same. The in-complex tests that change the count:
- These fail: derks2022 champion vs `gcn_nograph+stack` (MCAR −0.12, MNAR +0.59) and vs `gcn+stack` (MCAR −0.19,
  MNAR −0.41); GSTP1KO MCAR vs `gcn_nograph+stack` +0.19; MDAfix MCAR +2.22 and MNAR +0.63; RKOfix MCAR +0.55; and
  revised vs `gcn_nograph+stack` on RKOfix MCAR, −0.86.
- These pass: RKOfix MNAR champion vs `gcn_nograph+stack` (−3.23), and revised vs `gcn_nograph+stack` in-complex on
  derks2022 MCAR (−1.22) and MDAfix MCAR (−1.05).

**Absolute scale** (log2 RMSE; MCAR / MNAR):

| Set | champion | revised | `gcn_nograph+stack` | `stack_lowrank` |
|---|---|---|---|---|
| derks2022 | 0.922 / 1.142 | 0.910 / 1.107 | 0.920 / 1.135 | 0.951 / 1.195 |
| karagach2025_S1 | 0.493 / 0.627 | 0.477 / 0.589 | 0.499 / 0.642 | 0.534 / 0.705 |
| px061710_GSTP1KO | 0.384 / 0.478 | 0.375 / 0.456 | 0.379 / 0.486 | 0.426 / 0.556 |
| px061710_MDAfix | 0.342 / 0.424 | 0.331 / 0.399 | 0.333 / 0.415 | 0.368 / 0.486 |
| px061710_RKOfix | 0.399 / 0.464 | 0.391 / 0.449 | 0.392 / 0.467 | 0.436 / 0.500 |

Among the leading methods, the differences are a few hundredths of a log2 unit.

**Raw (unstacked) no-edge network.** From `champion_vs_every_method`, the champion was significantly worse than raw
`gcn_nograph` (CI lower bound above 0) in these cells:

| Set | Mask | Champion vs raw `gcn_nograph` |
|---|---|---|
| GSTP1KO | MCAR | +1.13% |
| MDAfix | MCAR | +2.68% |
| RKOfix | MCAR | +1.70% |
| MDAfix | MNAR | +2.30% |
| derks2022 | MNAR | +1.74% |

**Mean over masks and sets** (relative to the champion; `docs/v9_results.md` full-panel table):

| Method | Mean RMSE ratio |
|---|---|
| revised | 0.965 |
| `gcn_nograph+stack` | 0.998 |
| `gcn_nograph` | 0.998 |
| champion | 1.000 |
| `gcn+stack` | 1.033 |
| `champion@null` | 1.037 |
| `mice+stack` | 1.048 |
| `champion@noref` | 1.054 |
| `stack_lowrank` | 1.095 |

### 5.3 MCAR vs MNAR

The revision's margins are larger under MNAR than under MCAR on every set (tables above):
- vs the champion, about twice as large (1.7–2.5× per set): MCAR −1.23 to −3.22%, MNAR −3.03 to −6.07%;
- vs `gcn_nograph+stack`, 2–17× larger (per set: derks2022 2.5×, karagach2025_S1 1.9×, GSTP1KO 6.0×, MDAfix 7.3×,
  RKOfix 17×): MCAR −0.23 to −4.34%, MNAR −2.42 to −8.28%.

**Verdicts by mask** (post hoc; red-team recomputation in `results/v9/expansion_raw.json`, checked against the per-test
JSON):

| Reading | H16 | H17 | H18 |
|---|---|---|---|
| MCAR tests only | 1/3 | 3/3 | 1/3 |
| MNAR tests only | 1/3 | 3/3 | **3/3** |
| Both masks (preregistered) | 1/3 | 3/3 | 1/3 |

**Why MCAR is the cleaner test.**
- The MNAR evaluation mask is a synthetic, value-logistic generator: `expit((q25 − x)/(0.5·sd))` in
  `benchmarks/run.py::make_mask`.
- `cshrink` picks ψ on an inner mask that is half MCAR and half built with the same functional form.
- The champion and every `+stack` comparator choose their weights on a uniform (MCAR) inner mask.

**How much the ψ tuning explains is disputed.**
- Both skeptics in `results/v9/review_raw.json` argued that ψ selection explains little, because ψ = ∞ (no selection) gives
  dev MNAR 0.9509 against 0.9504 with selection (`docs/EXPLORATION_LOG.md`, round 5). On that view the MNAR lean comes
  from the shrinkage itself, not from tuning ψ to the generator.
- The other side: `docs/v9_review.md` still says the MNAR margins "should therefore be read as partly tuned to the
  simulator", and all four red-team reports keep that caveat (`results/v9/expansion_raw.json`, entries 19–22).
- **A second asymmetry, which the ψ = ∞ ablation does not address:** every `+stack` comparator chooses its weights on a
  uniform inner mask, and this hurts comparators that are strong under MNAR. On derks2022 MNAR (raw vs stacked RMSE,
  `results/v9/addendum/derks2022.json`; red-team report 4, entry 22): `gcn_nograph` 1.1223 vs 1.1347, `knn_protein`
  1.1067 vs 1.1234, `bpca_oba` 1.1038 vs 1.1056. The preregistered comparators are the stacked versions. This bears on
  the MNAR-only H18 reading of 3/3.
- Either way, the MNAR numbers measure performance on one synthetic missingness model.

### 5.4 Fragile passes and near misses

From the audit in `results/v9/expansion_raw.json`; "pct" is ΔRMSE as % of the comparator.

**Passes within 0.5 points of the SESOI:**

| Set | Test | pct | CI (%) |
|---|---|---|---|
| derks2022 | H17 MCAR vs champion, all | −1.231 | [−1.693, −0.798] |
| derks2022 | H17 MCAR vs champion, in-complex | −1.106 | [−1.573, −0.655] |
| derks2022 | H18 MCAR vs `gcn_nograph+stack`, in-complex | −1.222 | [−2.135, −0.296] |
| derks2022 | H18 MCAR vs `gcn+stack`, all | −1.245 | [−1.951, −0.527] |
| derks2022 | H18 MCAR vs `gcn+stack`, in-complex | −1.296 | [−2.111, −0.583] |
| karagach2025_S1 | H16 MCAR vs `gcn_nograph+stack`, all | −1.159 | [−1.758, −0.641] |
| px061710_GSTP1KO | H18 MCAR vs `gcn_nograph+stack`, all | −1.033 | [−1.760, −0.318] |
| px061710_MDAfix | H18 MCAR vs `gcn_nograph+stack`, in-complex | −1.051 | [−1.867, −0.293] |
| px061710_RKOfix | H17 MCAR vs champion, in-complex | −1.403 | [−1.779, −1.037] |

**Failures within 0.5 points:**

| Set | Test | pct | CI (%) |
|---|---|---|---|
| derks2022 | H16 MNAR vs `gcn+stack`, all | −0.503 | [−1.047, +0.013] |
| derks2022 | H18 MCAR vs `gcn_nograph+stack`, all | −0.980 | [−1.781, −0.133] |
| px061710_MDAfix | H18 MCAR vs `gcn_nograph+stack`, all | −0.535 | [−1.123, +0.066] |
| px061710_RKOfix | H16 MNAR vs `gcn_nograph+stack`, all | −0.781 | [−1.196, −0.363] |
| px061710_RKOfix | H18 MCAR vs `gcn_nograph+stack`, in-complex | −0.857 | [−1.401, −0.315] |

**What this means:**
- **H16's only project pass (karagach2025) rests on one test at −1.16%,** with a CI upper bound of −0.64%.
- **derks2022's H18 failure rests on a −0.98% result against a −1.00% threshold.**
- **The verdicts move with the SESOI** (post hoc red-team sweep, `results/v9/expansion_raw.json`):

  | Rule | H16 | H17 | H18 |
  |---|---|---|---|
  | Preregistered, 0.5% | 1/3 | 3/3 | 2/3 (confirmed) |
  | Preregistered, 1% | 1/3 | 3/3 | 1/3 |
  | Preregistered, 1.5% | 0/3 | 2/3 | 1/3 |
  | Strict, 1% | 0/3 | 2/3 | 1/3 |
  | Strict, 2% | 0/3 | 1/3 | 1/3 |

- **The SESOI itself is a carried-forward convenience number.** It first appears after v3 as "for example 1% of RMSE".
  In absolute terms it is about 0.009–0.012 log2 units. Because the rule applies it to the point estimate, a true effect
  sitting exactly at the SESOI passes about half the time (red-team report 2, `results/v9/expansion_raw.json`, entry 20).
- **With projects as the unit, the effect against the hardest comparator is not pinned down.** The red team ran a post
  hoc t interval on 2 df over the 3 project means (same entry):

  | Contrast (MCAR) | Mean | 95% CI | Prediction interval for a new project |
  |---|---|---|---|
  | champion vs `gcn_nograph+stack` (decides H16) | +0.30% | [−3.4, +4.0] | [−7.1, +7.7] |
  | revised vs `gcn_nograph+stack` (decides H18) | −1.97% | [−7.1, +3.1] | — |
  | H17, revised vs champion | −2.27% | [−4.75, +0.20] | — |

- **Project counts carry little information on their own.** Exact binomial 95% intervals for the per-project success
  rate are [0.008, 0.906] for 1/3, [0.094, 0.992] for 2/3 and [0.292, 1.000] for 3/3 (same entry).

### 5.5 What the graph-prior methods show

**GCN with STRING vs its controls** (all entries; raw networks):

| Set | Mask | STRING vs shuffled STRING | STRING vs no graph |
|---|---|---|---|
| derks2022 | MCAR | −7.14% | +6.07% ✗ |
| derks2022 | MNAR | −6.00% | +5.09% ✗ |
| karagach2025_S1 | MCAR | −5.20% | +7.23% ✗ |
| karagach2025_S1 | MNAR | −4.75% | +3.50% ✗ |
| px061710_GSTP1KO | MCAR | −9.91% | +9.37% ✗ |
| px061710_GSTP1KO | MNAR | −8.92% | +5.48% ✗ |
| px061710_MDAfix | MCAR | −10.42% | +14.86% ✗ |
| px061710_MDAfix | MNAR | −11.27% | +11.83% ✗ |
| px061710_RKOfix | MCAR | −8.23% | +6.99% ✗ |
| px061710_RKOfix | MNAR | −8.40% | +3.00% ✗ |

All 40 CIs (both contrasts × 10 cells × 2 endpoints) lie entirely on the side shown (`results/confirm_v9/*.json`).

**What this shows:**
- STRING is less harmful than a degree-matched random graph, so it carries some real structure.
- It is more harmful than no graph, in every cell.
- The red team found the STRING-vs-shuffled gain is about the same for complex members and non-members. For example,
  derks2022 MCAR is −7.1% vs −7.2%, and RKOfix MCAR is −7.8% vs −8.7% (`results/v9/expansion_raw.json`). If co-complex
  stoichiometry were the mechanism, the gain should concentrate in complex members. It does not.
- **Architecture caveat:** the GCN has symmetric-normalised A + I and no skip connection, so a protein with d neighbours
  keeps self-weight 1/(d+1) at each hop. These results are a statement about this architecture. They do not show that
  STRING carries no usable information.

**Stacked versions** (mean RMSE ratio to the champion; `docs/v9_results.md`):

| Set | `gcn+stack` | `gcn_nograph+stack` |
|---|---|---|
| derks2022 | 1.003 | 0.996 |
| karagach2025_S1 | 1.037 | 1.018 |
| px061710_GSTP1KO | 1.051 | 1.002 |
| px061710_MDAfix | 1.050 | 0.976 |
| px061710_RKOfix | 1.023 | 0.996 |

**The other STRING methods** (means over the 5 sets):

| Method | Mean RMSE ratio | Comparison |
|---|---|---|
| `rwr_resid+stack` | 1.045 | data-only `mice+stack`: 1.048 |
| `imc_embed+stack` | 1.047 | data-only `mice+stack`: 1.048 |
| `grmf_string+stack` | 1.072 | data-graph control `grmf_datagraph+stack`: 1.089 |
| `grmf_string` (raw) | 1.173 | `grmf_datagraph` (raw): 1.284 |

- These three methods had no shuffled-STRING null in v9.
- Their stacked gains are about the same as those of a generic data-only third component.

---

## 6. Independent review and red team

**Sources:**
- `docs/v9_review.md` and `results/v9/review_raw.json`: 39 agents; four reviewers, each finding given to two skeptics;
  a literature check.
- `results/v9/expansion_raw.json`: 23 agents; new comparators each checked by two adversarial reviewers; four red-team
  reports; an independent verdict audit.

Nothing below changed the frozen v9 code or any preregistered verdict.

| Finding | Status after the skeptics | What was done |
|---|---|---|
| **Frozen BPCA was broken.** It treated the ~30 cells as observations, so ARD pruned every component. It equals per-protein mean imputation on 20 of 24 dev splits. On derks2022, `bpca` = `mean` exactly (MCAR RMSE 1.1079025) and `bpca+stack` = `stack_lowrank` (`results/confirm_v9/derks2022.json`). The selection never saw a working BPCA. | Confirmed by both skeptics: critical in the review; in the skeptics' corrected impact, "moderate, not critical" (verdict 28) and "high severity … rather than critical" (verdict 29, `results/v9/review_raw.json`). Both judged it unlikely that a working BPCA would displace `gcn_nograph+stack`; the post-hoc dev result below narrowly contradicts that | Post hoc `bpca_oba` (proteins as observations, as in pcaMethods/MSnbase). On dev (post hoc), `bpca_oba+stack` scores 0.9914 against 0.9918 for `gcn_nograph+stack` (0.9917 in `v9_selection.json`; 0.9918 is the addendum's recomputation). So under the frozen rule a working BPCA would narrowly have been selected (`docs/v9_addendum_results.md`, generated from `results/v9/addendum/dev_*.json`). Holdout results: §7. |
| **Entry bootstrap ignores clustering** by protein, cell and split. The SE is understated by about 1.5–2×. In the review's dev proxy, a cell-cluster bootstrap inflated the SE 3.42×, a seed-cluster bootstrap 1.89× and a protein-cluster bootstrap 1.27× (`results/v9/review_raw.json`, verdict 7). The skeptics' CI-widening check (upper bound moved to point + k × (upper − point)) was run when only derks2022 and karagach2025_S1 had finished. Recomputed here on all five sets (post hoc, from `results/confirm_v9/*.json`): every project-level verdict holds at 2×; at 2.5× and 3×, karagach2025 H16 and derks2022 H17 fail, and no other project-level verdict changes. GSTP1KO's set-level H18 pass, px061710's only H18 set, already flips at 1.5×. | Confirmed (rated major, then moderate) | Post hoc, the addendum reruns the splits, keeps per-entry errors with coordinates, and recomputes every test with a protein-cluster bootstrap and a per-split t interval. The per-split interval replaced a degenerate split bootstrap. On derks2022 the protein-cluster CI is a median 1.42× wider than the entry CI (range 1.05–2.23), and the per-split CI a median 1.25× (`results/v9/expansion_raw.json`). The per-split interval is itself anticonservative, because the 10 subsamples overlap (red-team report 2, entry 20). No cell-cluster bootstrap has been run. Only derks2022 is finished; the px061710 sets, with the most cell overlap, have no robustness check yet. |
| **SESOI wording is ambiguous.** The code applies 1% to the point estimate, as v8 did. A strict reader would require the CI bound to clear −1%. | Both skeptics held that the text names the point estimate, as the code applies it; the review kept the finding because a reader can still take the strict view | Every test re-scored under the strict reading (`docs/v9_results.md`; §1, §5.4). The decisive reading should be fixed *before* the remaining addendum and secondary results arrive. |
| **`cshrink`'s MNAR margin may be tuned to the simulator.** Its ψ-selection mask shares the MNAR generator's form. | Both skeptics refuted the material impact: ψ = ∞ gives the same dev MNAR (§5.3). The review and all four red-team reports keep the caveat, and the uniform stack mask for comparators is a second, unrefuted asymmetry (§5.3) | Kept as a caveat. MCAR is treated as the clean test, and the mask-specific verdicts are reported (§5.3). |
| **`cshrink` missed its own Phase-C promotion bar** (full MCAR 0.9714 > 0.97; woo2025_PC3 MCAR 1.0006 > 1.00). The preregistration does not disclose this. | Confirmed (minor; not a protocol breach, since v9 was requested separately) | Disclosed here and in `docs/v9_review.md`. H17 and H18 should be read with it. |
| **Panel coverage gaps.** Missing: LLS, QRILC, Perseus down-shift, kNN over proteins, SeqKNN, impSeq, EM, msImpute, a working BPCA, and every dedicated proteomics/SCP tool (PIMMS, omicsGMF, Lupine, DreamAI, GIRAFI, SoftHybrid). Both H18 comparators are in-house GCN variants. | Confirmed by the literature check | **Addendum 1:** `bpca_oba`, LLS, `knn_protein`, QRILC, Perseus, PIMMS CF and VAE. **Addendum 2:** an R bridge (pcaMethods BPCA/LLS, impute.knn, imputeLCMD QRILC, missForest defaults, msImpute) and Python ports of SeqKNN, trKNN, impSeq/impSeqRob, RegEM and simple single-value imputers. Each was checked by two adversarial reviewers (commit `e23cbcf`; `results/v9/expansion_raw.json`). The commit message says "no confirmed bugs"; the reviewers' finding is narrower: no confirmed bug that corrupts results (the wording of entry 13). They did confirm two edge-case crashes (`rem_protein` on a rank-0 pattern, entry 17; `r_knn` when a cell is more than 80% missing among the kept proteins, entry 13), one edge-case mismatch with the reference (`r_knn` vs `impute.knn` when exactly one protein exceeds rowmax, entry 13), an upstream pcaMethods `llsImpute` indexing defect that affects the `r_lls` dev numbers (entry 12), and wrong citation attributions in `simple.py` (entry 6). None of these is fixed in the code as of the current commit. Addendum 2's holdout run (§7) scores a subset of these methods; msImpute is bridged but not in that run. **Not added:** omicsGMF, Lupine, GIRAFI, DreamAI, SoftHybrid. **Wording:** "strongest typical methods" in H18 should read "strongest methods in our panel". |
| **The mask design cannot test left-censored methods or real missingness.** Only *observed* values are masked, only in proteins seen in ≥ 50% of the subsample. On derks2022, 538 of 1093 proteins carry masked entries under MCAR and 537 under MNAR (`n_clusters` in `results/v9/addendum/derks2022.json`). Entries are masked independently, with no per-cell, per-run or per-plex structure. No primary holdout is DDA or TMT. QRILC (RMSE 5.37) and Perseus (3.60) against the champion's 0.92 on derks2022 MCAR reflect this design, not bugs. | Confirmed (red team, objection 1) | Not addressed in v9. Proposed: pattern-transplant masks, run or plex block masks, errors by detection rate and intensity decile, and the PXD061710 multi-cell runs as a check on truly missing values (§9). |
| **No downstream metric.** RMSE can reward over-smoothing, and `cshrink` explicitly shrinks per-protein scale. | Confirmed (literature check, red team) | Not addressed. Proposed on dev: calibration slope, variance ratio, fold-change recovery, DE-rank overlap, CORUM-pair correlation recovery, clustering ARI (§9). |
| **Tuning effort was lopsided** (red-team report 4, `results/v9/expansion_raw.json`, entry 22, ranked within its most serious objection) | Raised by the red team; not a finding of the 39-agent review, so no skeptic rated it | Not addressed. The champion was tuned on dev over Phases B–C (λ_v = 100 with 10 extra factors; ALS rank 20, λ = 5). Every comparator ran one fixed configuration, including the decisive no-edge network (400 epochs, lr 1e-2, no validation or early stopping). `champion@noref` reuses the champion's reference-tuned settings and unit-variance scaling, so it is "not a fair stand-in for a data-only method" (same entry). Proposed: equal inner-mask tuning for every comparator (§9). |
| **Other points** | Minor | The DAE overfits (no early stopping), so it is a weak stand-in for deep models. missForest and MICE were cut back further than disclosed (50 nearest features; MICE is single imputation). `champion@null` uses one permutation; spread between draws was ≤ about 0.25% in v6/v7 (`results/v9/review_raw.json`). The restart log was gitignored (no restart occurred). The dry-run JSON records an earlier commit hash. The in-complex endpoint is 61–82% of masked entries (61.2% on RKOfix MCAR to 81.7% on derks2022 MNAR, from `n_entries` in `results/confirm_v9/*.json`), so it is not an independent test. |

**Red-team "claims that survive"** (`results/v9/expansion_raw.json`, entries 19–22). Not every report states every item
(report 3, entry 21, for example, makes no H18 claim), but taken together the four reports support this list:
1. Within the champion's own architecture, the protein identity of the external reference matters.
2. The champion beats the kNN + low-rank stack.
3. H16 is not confirmed.
4. `cshrink` improves the champion under these masks.
5. H18 is not confirmed, and the panel's strongest comparators are in-house.
6. In this GCN, STRING is less harmful than a shuffle and more harmful than no graph.
7. No v9 method tests curated complexes.

---

## 7. Post-hoc addendum (exploratory)

**This section is post hoc.** Its results are being filled in from `results/v9/addendum/` and `results/v9/addendum2/`.
The addendum has finished for derks2022 and six dev targets. karagach2025_S1 started at 2026-11-08T05:05:53Z
(`results/v9/addendum_runs.log`), and runs on px061710_GSTP1KO and px061710_MDAfix have started
(`results/v9/addendum_px061710_*.log`). Addendum 2 has finished for derks2022 only (`results/v9/addendum2_runs.log`).

**Known so far for derks2022** (`results/v9/addendum/derks2022.json`):

- **Exact reproduction.** The max |RMSE difference| from the frozen run is 0.0 under both masks.
- **The verdicts are unchanged under protein-cluster and per-split resampling.** Only H18's count moves, from 11/12 to
  10/12 under the per-split interval. These are the only two resamplings run, on derks2022 only; no cell-cluster
  bootstrap exists, and the per-split interval is itself anticonservative because the subsamples overlap (§6):

  | Interval | H16 | H17 | H18 |
  |---|---|---|---|
  | Entry bootstrap | 12/20 | 4/4 | 11/12 |
  | Protein-cluster bootstrap | 12/20 | 4/4 | 11/12 |
  | Per-split t interval | 12/20 | 4/4 | 10/12 |

  Under the strict reading, H17 is 2/4 and H18 is 8/12 under all three intervals.
- **A correctly oriented BPCA beats the champion.** Relative to the champion, `bpca_oba` is 0.978 (MCAR) and 0.967 (MNAR).
  The champion is worse by:
  - MCAR: +2.23%, entry CI [+1.36, +3.08], protein-cluster [+1.27, +3.18];
  - MNAR: +3.44%, entry CI [+2.72, +4.14], protein-cluster [+2.63, +4.21].
- **Against the revision it narrowly wins under MCAR and ties under MNAR.**
  - MCAR, all entries: `revised` is worse by +0.97%, entry CI [+0.14, +1.80], protein-cluster [+0.10, +1.83]. The
    per-split interval [−0.16, +2.13] includes 0.
  - MNAR, all entries: +0.31%, CI [−0.37, +0.97], a tie.
  - BPCA's point estimate is better in all four mask × endpoint cells (+0.97%, +0.89%, +0.31%, +0.98%), but every one
    is below the 1% SESOI, so under the preregistered rule BPCA would not "beat" the revision in any cell.
- **Counterfactual: had a working BPCA been the frozen data-only comparator,** every derks2022 test against
  `bpca_oba+stack` fails, for both the champion and the revision (8 of 8; `docs/v9_addendum_results.md` §3).
- **Dev corroboration is partial.** On dev, `revised` beats `bpca_oba+stack` in 11 of 12 target × mask cells; the
  exception is gregoire2023_mixCTRL (TMT) MCAR, 0.941 for `bpca_oba+stack` against 0.967 for `revised`
  (`docs/v9_addendum_results.md` §4). These are the same dev data on which `cshrink` was selected.
- **Other added comparators:**
  - Raw kNN over proteins (`knn_protein`) beats the champion under MNAR (champion +3.18%) but not under MCAR (champion
    −5.10%).
  - Stacked PIMMS-VAE (1.014 / 1.030) and stacked LLS (1.016 / 1.035) do not beat the champion (MCAR / MNAR ratios).

**Addendum 2 on derks2022** (`results/v9/addendum2/derks2022.json`; post hoc, seeds 9000–9003 only, RMSE relative to
the champion on those 4 splits, no CIs, and `revised` not scored):

| Method | MCAR | MNAR |
|---|---|---|
| pcaMethods BPCA (`r_bpca`) | 0.947 | 0.934 |
| missForest, package defaults (`r_missforest`) | 0.985 | 0.947 |
| RegEM, protein-space dual form (`rem_protein`) | 0.997 | 0.963 |

- `r_bpca` is below the champion on all 4 splits under both masks (`per_split_rel_to_champion`).
- These ratios use 4 of the 10 preregistered splits, so they are not directly comparable with the 10-split `bpca_oba`
  ratios above.

These results bear on H16 and H18 but do not change them. They are exploratory: a holdout already spent, and a
comparator added after the results were known.

---

## 8. What this means for the original thesis and for the project

**The thesis** (`AGENTS.md`): curated complex membership (CORUM, Complex Portal) is a better prior for imputing
single-cell protein values than data alone, than STRING, and than a degree-preserving shuffle of the same graph.

**v9 does not test it.**
- No v9 method uses CORUM or Complex Portal.
- They enter only through the in-complex endpoint, which covers 61–82% of the masked entries.
- The one knowledge graph v9 tested, STRING in a GCN, beat its shuffle but lost to no graph in every cell (§5.5).

**The earlier holdouts that did test curated priors** (`docs/EXPLORATION_LOG.md`):

| Test | Result |
|---|---|
| v2 (CORUM Laplacian) | 6/8, not confirmed: no better than its own shuffled null |
| v3 | 8/8 by the letter, but the effect was about −0.01% (β ≈ 0.02); it would fail a 1% SESOI |
| v4 | 2/6 |
| v5 | 0/4 |
| v7 secondary | BIM vs plain transfer −0.17% / −0.06%, CIs including 0 |

Related dev result, not a holdout: in Amendment 13 (`results/phase9/dev/`), the complex component got stack weight > 0
in 1 of 36 subsamples.

- On dev, curated priors usually beat shuffled versions of themselves. On holdouts the record is mixed: v2's primary was
  no better than its own null (v2's secondary E6 beat its null by about 0.5% of RMSE, CIs excluding 0); v3 beat it by
  about 0.01%; v4 beat it under MCAR (−1.05%) but not under MNAR (CI including 0); v5 (−0.38% / −0.21%) and v7's BIM
  test (−0.54% / −0.34%) had CIs including 0. Only v4 MCAR cleared 1%.
- On holdouts, Complex Portal partners did beat data-derived partners in some tests (v4 secondary: −5.3% / −2.8% at
  30 cells; v5 secondary: −1.45% / −1.48% at n = 60, CIs excluding 0).
- Curated priors never beat a strong data-only or reference-transfer imputer by a meaningful margin.
- The log's own reading after Amendment 13 still fits: the prior is correct but redundant once protein covariance can be
  learned from about 10k external cells, and too weak on its own at 30 cells.

**Claims we can defend, worded narrowly:**
1. **Protein-specific information in public single-cell proteomics data transfers to new datasets.**
   - Setting: 30-cell random subsamples of 5 sets from 3 DIA projects, scored on simulated masks of detected values.
   - Borrowing loadings from 22 public SCP datasets lowered RMSE by 2.7–4.1% against a protein-permuted reference and by
     2.6–8.3% against no reference, on every set and both masks.
   - Counting brunner2022 (v7) and wu2026 (v8), reference transfer has beaten its scrambled-reference control by at
     least the 1% SESOI on 5 of 6 holdout projects. On khan2023 it was significantly better, but by only 0.20–0.61%
     (`results/confirm_v6/confirm_v6.json`). Three different transfer methods are involved.
   - Both controls are internal ablations with limits: one permutation and one ALS fit per holdout, and a
     no-reference arm that reuses the champion's reference-tuned settings (§6).
   - The permutation null cannot tell biological co-variation apart from per-protein scale.
   - The reference is empirical, not curated.
2. **The champion clearly beats the kNN + low-rank stack** (3.0–14.0%), the strongest data-only baseline of v6–v8.
3. **It does not beat the strongest data-only methods we now know of:**
   - it ties a no-edge masked MLP on average (0.998);
   - it is significantly worse than that MLP in 5 of 10 set × mask cells (all entries), including all three px061710
     sets under MCAR;
   - post hoc, it loses to a working BPCA on derks2022 (§7).

   H16 is not confirmed.
4. **A data-only change to the target's standardisation (`cshrink`) improves the champion** (H17 confirmed).
   - MCAR: 1.2–3.2%; MNAR: 3.0–6.1% (all entries).
   - It passes the strict reading on 2 of 3 projects.
   - It is a preprocessing choice found by search on dev, it missed its Phase-C bar, and whether it needs the reference is
     untested.
5. **"Beats the strongest typical methods" is not supported** (H18: 1 of 3 projects). Post hoc, a working BPCA ties or
   narrowly beats the revision on derks2022, by less than the 1% SESOI in every cell (§7).

**The honest conclusion.**
- After nine preregistered phases, the project has one positive result that has replicated across phases: transfer of
  protein-specific structure from public single-cell proteomics data. v9 adds a second, so far confirmed only in v9:
  the `cshrink` revision improves the champion (H17, 3 of 3 projects).
- These results come out of a long search. Phases A–C tried more than 25 ideas, and the claim under test moved from
  curated complexes, to reference transfer, to target standardisation (red-team report 3, `results/v9/expansion_raw.json`,
  entry 21). `cshrink` came from about 20 Phase-C ideas on the same 6-target dev pool and is the explicit form of the
  failed `tgtshrink` idea's null pipeline; a random-ψ control already gets about 73% of its dev gain (red-team report
  4, entry 22).
- The transfer result is empirical and data-derived. Its advantage over the best data-only imputers is unresolved and
  looks dataset-dependent. derks2022 is the only v9 primary holdout with distinct cell types, and the only holdout on
  which `bpca_oba` has been scored so far. On it the gain over the kNN + low-rank stack is smallest (3.0–4.4%), the
  champion's stack weight on the projection drops to about 0.5 (0.88–1.00 on the other four sets, `stack_weights` in
  `results/confirm_v9/*.json`), and, post hoc, a working BPCA beats the champion (§7; Addendum 2's pcaMethods BPCA,
  on 4 splits without CIs, agrees).
- Dev data partly corroborate derks2022's pattern (red-team reports 1 and 2, entries 19–20): on gregoire2023_mixCTRL
  (TMT), `bpca_oba+stack` beats `revised` under MCAR (0.941 vs 0.967); Phase C showed a TMT/carrier MCAR weak spot
  (woo2025_PC3 1.0006, williams2020_tmt 0.9989, schoof2021 0.9972, `docs/EXPLORATION_LOG.md`); and raw `knn_protein`
  beats the champion under MNAR on dev sadiku2025 (0.959) and weiss2026_liver (0.955), not only on derks2022. For
  balance, `revised` beats `bpca_oba+stack` on the other 11 of 12 dev cells (`docs/v9_addendum_results.md` §4), which
  are the dev data used to select `cshrink`.
- No test so far shows that curated biological knowledge, and complexes in particular, improves single-cell proteomics
  imputation by a meaningful margin.
- The thesis as stated in `AGENTS.md` should be reported as **not supported**.
- The project's useful output is a rigorous negative result on curated priors, plus a moderate reference transfer
  method that is competitive with, but not shown to be better than, the strongest data-only imputers we have tried. Its
  controls have known limits: one permutation and one ALS fit per holdout, a permutation null that cannot separate
  co-variation from per-protein scale, and a weak no-reference arm.

---

## 9. Limitations and next steps

### Limitations of v9

- **Few independent units.** There are 3 projects, and the px061710 sets are correlated (shared lab, and RKO twice). The
  30-cell subsamples overlap heavily on the 75–78-cell px061710 sets.
- **Narrow biology.** Four of the five sets are a single cultured line with technical contrasts. None is primary tissue,
  and none is DDA or TMT.
- **Simulated missingness only.** No natural missing value, no value below the detection limit, and no protein seen in
  fewer than 50% of cells was scored.
- **Only n = 30 was evaluated.** Earlier phases found the reference gains shrink with more cells.
- **No downstream analysis:** no differential expression, clustering or correlation recovery.
- **Unequal tuning** favours the champion and the revision: they were tuned over Phases B–C, every comparator ran one
  fixed configuration (§6). **Unequal stack slots** favour the revision (4 slots against 3 for every other method).
- **The stack mask is uniform for every comparator,** which handicaps comparators that are strong under MNAR (§5.3).
- **One reference fit and one null permutation per holdout.** Phase B measured ALS run-to-run spread at about 0.5%.
- **Within-dataset CIs understate uncertainty.** The robust intervals (protein-cluster and per-split only; no
  cell-cluster bootstrap) are available for derks2022 only.
- **Weak inference at the project level.** With 3 projects, exact binomial 95% intervals on the per-project success rate
  are wide (§5.4). The 1% SESOI is a convenience number applied to the point estimate. The roughly 780 panel CIs and the
  addendum CIs have no multiplicity control.
- **Circularity for downstream use.** Imputing through reference loadings injects external co-variation, so any
  downstream analysis of complex co-variation on imputed data would be partly manufactured (red-team report 1,
  `results/v9/expansion_raw.json`, entry 19).
- **No noise floor.** The "truth" values carry measurement noise and there is no replicate-based noise floor, so the
  fraction of reducible error captured is unknown (same entry).
- **Only squared errors were saved,** so signed bias and calibration cannot be checked without reruns (same entry).
- **Same-lab and same-cell-line reference datasets** for derks2022 and karagach2025; overlap with the px061710 cell
  lines was not checked.

### Next steps, in priority order

1. **Finish and report the addendum** on karagach2025_S1 and the px061710 sets, then the secondary holdouts. Keep
   everything labelled post hoc.
2. **Fix the decisive SESOI reading now,** before those results arrive.
3. **Settle the data-only comparison on dev before any new holdout:**
   - stack(kNN, low-rank, `adapt`, best data-only method) against stack(kNN, low-rank, best data-only method);
   - `revised@null` and `revised@noref`;
   - every comparator given the same inner-mask tuning and the same number of stack slots.
4. **Test what the reference actually contributes:**
   - a norm-preserving null (keep each loading row's norm and randomise its direction);
   - a diagonal-only reference-variance prior;
   - references that exclude the same lab and cell line.
5. **Give curated complexes a fair, direct test** inside the champion, against degree-preserving and no-graph controls,
   on dev first. Options:
   - a κ-weighted CORUM/Complex Portal prior on the loadings;
   - a GCN with a skip connection or a learned self-weight.
6. **Make the benchmark closer to real use:**
   - pattern-transplant and run/plex block masks;
   - errors by detection rate and intensity;
   - n = 30, 60, 100 and all cells;
   - at least one downstream metric (ARI, fold-change or DE recovery).
7. **Add a dedicated SCP imputer** (omicsGMF, Lupine) where installable, and score msImpute on the holdouts (it is
   bridged in Addendum 2 but not in that addendum's holdout run).
8. **Choose the next holdouts deliberately:**
   - DDA or TMT;
   - primary tissue;
   - several real cell types;
   - no lab or cell-line overlap with the reference.
9. **Report a project-level random-effects estimate with a prediction interval** alongside any vote count.
