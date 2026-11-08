# Independent review of the v9 validation (2026-11-08)

**How it was run:**
- A read-only `Workflow` ran 39 agents while the preregistered holdout runs executed.
- Four reviewers each took one dimension: holdout leakage, statistics, baseline fidelity, preregistration compliance.
- Every finding was then given to two skeptics, each told to refute it.
- A literature agent checked comparator coverage against published benchmarks.

The raw output is in the session. This file records each finding and what was done about it.

Nothing below changed the frozen v9 code or any preregistered verdict. Fixes and robustness checks are **post hoc**:
- `src/scp_graphfs_mcp/baselines_addendum.py`;
- `benchmarks/addendum_v9.py`;
- results in `results/v9/addendum/`.

## Leakage: clean

The reviewer checked six things and found no path from any holdout value into a reference, a tuning choice, a cache used
by another evaluation, or the comparator selection:
- reference exclusions;
- the corpus directory contents;
- the cshrink score-cache keys;
- the inner and selection masks, which use only observed entries;
- the graph caches;
- the selection timeline: the final preregistration was committed at 23:34:37 and the first holdout started at
  23:34:53.

**Minor items (no effect on any number):**
- **cshrink score cache:** the key does not hash V, warm or base. It is correct under frozen code, but a restart after
  a code fix should first delete that holdout's cache entries.
- **Inner mask:** its eligibility comes from the pre-mask subsample. This is shared by every stacked method, so it
  largely cancels in paired tests.
- **brunner2022:** this secondary holdout sits in the corpus directory, protected only by the `reference_datasets()`
  filter.
- **Same lab, same cell lines:** derks2022's reference contains leduc2022 and specht2021, which come from the same lab
  and share cell lines with it. They are different runs, but this should be stated.

## Statistics

**1. The entry bootstrap understates uncertainty** (major; confirmed by both skeptics).
- Masked entries cluster by protein, cell and split. The reviewer estimates the standard error is understated by about
  1.5–2×.
- The finished verdicts hold if every CI is widened 2×. At 3×, karagach2025_S1 H16 and derks2022 H17 would fail.
- **Done:** the addendum rerun keeps per-entry errors with coordinates. It recomputes every H16/H17/H18 test with the
  preregistered entry bootstrap, a protein-cluster bootstrap and a split-level bootstrap.

**2. The SESOI wording is ambiguous** (major).
- The code applies the 1% minimum to the point estimate, as v8 did.
- A strict reader could require the CI's upper bound to be at or below −1%.
- **Done:** every test is also re-scored under the strict reading, in the addendum and in `docs/v9_results.md`.

**3. Minor statistical points:**
- Each holdout's scrambled-reference control uses one permutation.
- `revised` tunes ψ on an inner mask built with the same MNAR generator as the evaluation mask. Its MNAR margins should
  therefore be read as partly tuned to the simulator; its MCAR margins are the clean test.
- The preregistration's verdict logic is implemented exactly as written.

## Baseline fidelity

**1. BPCA ran on the wrong axis** (critical; confirmed by both skeptics, and checked by hand).
- `baselines.bpca` used the ~30 cells as observations, so ARD pruned every component.
- On 20 of 24 dev splits it is numerically identical to per-protein mean imputation.
- Oba's BPCAfill and MSnbase use proteins as observations.
- The frozen comparator selection therefore never saw a working BPCA.
- **Done:** `baselines_addendum.bpca_oba`, correctly oriented, is scored post hoc on dev and on every holdout. Early
  dev check: 1.002 and 1.067 × the adapted projection on two splits, against 1.057 and 1.299 for the broken version.
- The write-up reports whether it would have displaced `gcn_nograph+stack` as the strongest data-only method.

**2. The transductive DAE overfits** (no early stopping), so it is a weak stand-in for the deep-learning family.
PIMMS's own paper advises against autoencoders below 50 samples.

**3. missForest and MICE settings:**
- missForest was cut to 10 trees, 3 iterations and the 50 nearest features. This is a disclosed handicap, but the
  50-feature limit was not disclosed.
- MICE is single imputation.
- `missforest_full` (100 trees, 10 iterations) exists in the addendum for a dev-only check.

**4. Checked and correct:**
- `gcn_nograph` has no edges.
- `gcn_null` is degree-preserving.
- Observed entries are never overwritten.

## Preregistration compliance

- Every hypothesis, comparator, seed range, holdout and pass rule is implemented as written.
- No code changed after the first holdout result.
- **The mice+stack → gcn_nograph+stack amendment:**
  - It was committed before any holdout ran.
  - It made H16 and H18 *harder*: under the draft comparator, derks2022 H18 would very likely have passed.
  - Filing the no-edge GCN under "controls" while selecting it as the best data-only method is an internal
    inconsistency in the text.
- **Not disclosed in the preregistration:** `cshrink` missed its own Phase-C promotion bar (full MCAR 0.9714 > 0.97;
  woo2025_PC3 MCAR 1.0006 > 1.00) before Amit asked for v9. H17 and H18 should be read with that caveat.
- **Overstated prediction:** the preregistration says the revision beats `gcn_nograph+stack` "on every dev split". The
  true count is 22 of 24 splits, or 11 of 12 target×mask cells.
- **Restart log:** the log §5 requires is gitignored. It will be committed with the final results.
- **Guard and dry-run commit:** the integrity guard does not cover untracked files or holdout data, and the dry-run JSON
  records an earlier commit hash.

## Literature coverage

The literature check found that the frozen panel leaves out several of the most commonly reported proteomics imputers:
- LLS;
- QRILC;
- Perseus down-shift;
- kNN over proteins (the `impute.knn` orientation);
- SeqKNN and impSeq;
- EM;
- msImpute;
- a working BPCA.

It also leaves out every dedicated proteomics or SCP imputation tool: PIMMS CF/VAE, omicsGMF, Lupine, DreamAI, GIRAFI and
SoftHybrid.

So "the strongest typical methods" in H18 overstates what was tested: both H18 comparators are our own GCN variants. The
accurate phrase is "the strongest methods in our panel".

**Added post hoc:** a working BPCA, LLS, kNN over proteins, QRILC and Perseus down-shift.

**Not added in this window, and why:**

| Tool | Reason |
|---|---|
| omicsGMF, msImpute | R/Bioconductor, and there is no R here |
| Lupine | GitHub install and retraining on our corpus |
| GIRAFI | code release not confirmed |
| PIMMS | attempted in an isolated environment (see `docs/v9_results.md`) |

**Benchmark limitation:** masking *observed* values cannot fairly score left-censored methods (QRILC, MinProb, Perseus),
because they assume missing values sit below the detection limit. Their large errors here reflect that, not a bug.

The literature also asks for at least one downstream metric, such as clustering agreement (ARI) or differential
expression, beyond RMSE. v9 has none.

**Verified citations from the literature check:**

| Source | DOI |
|---|---|
| NAguideR | doi:10.1093/nar/gkaa498 |
| Jin et al. 2021 | doi:10.1038/s41598-021-81279-4 |
| Harris et al. 2023 | doi:10.1021/acs.jproteome.3c00205 |
| PIMMS | doi:10.1038/s41467-024-48711-5 |
| omicsGMF | doi:10.1038/s41467-026-73402-8 |
| Lupine | doi:10.1016/j.mcpro.2025.101461 |
| SoftHybrid | doi:10.1021/acs.jproteome.5c01221 |
| msImpute | doi:10.1016/j.mcpro.2023.100558 |
| Vanderaa & Gatto | doi:10.1021/acs.jproteome.3c00227 |

## Erratum: hypothesis labels (found by the write-up agents)

- `benchmarks/PREREGISTRATION_v9.md` §0 calls the khan2023 v6 test "H14". In `PREREGISTRATION_v6.md` and
  `docs/EXPLORATION_LOG.md` that hypothesis is **H12**; H14 was the brunner2022 (v7) test.
- The same mislabel appeared in the session's messages to Amit.
- The preregistration is frozen and has been left as committed. The correct record is:
  - H12 on khan2023 (v6): not confirmed, 3/8;
  - H14 on brunner2022 (v7): confirmed;
  - H15 on wu2026 (v8): confirmed 16/16.
- No test, number or verdict depends on the label.
