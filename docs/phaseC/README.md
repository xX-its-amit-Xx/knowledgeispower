# Phase C: the colony (king → queens → workers)

**Goal:** keep climbing past the H15 champion, which was confirmed on wu2026 at −5.0 to −6.8% RMSE vs the best
data-only method. Amit's target is another 20–30%. The method is to bake in as much biological inductive bias as
possible, from as many data sources and representations as we can.

**Honest constraint (logged before we start):**
- The Phase-B ceiling diagnostic used an ideal in-domain reference (other cells of the *same* dataset).
- That ideal reference reached −14.5% vs data-only on dev.
- So a 20–30% mean gain over the champion probably needs information that no reference-projection method has.
- That is why this phase searches *new information sources and encodings*, not just tuning.
- Any number we report is a dev number until a preregistered holdout says otherwise.

## Structure

Amit asked for this structure on Oct 23. It runs as a `Workflow` script; the script itself lives in the session, and
this file is the durable spec.

| Role | Count | Job |
|---|---|---|
| **King** | 1 | Reads the ledger after every round. Detects a queen stuck in a local minimum (two rounds with no idea at ≤ 0.99 vs the champion) and pushes her to an adjacent or brand-new topic. Updates `LEDGER.md`, `KING_LOG.md`, then commits and pushes. |
| **Queen: biological priors** | 1 | Curated knowledge as inductive bias: tissue- and technology-aware priors, complexes inside the adapted model, pathways, localisation, half-lives, abundance. |
| **Queen: representation & math** | 1 | How the information is encoded: embedding geometry (hyperbolic, higher-dimensional), kernels/GPs, per-protein adaptation strength, synthetic features from simple sub-models, optimal transport between domains. |
| **Queen: new data modalities** | 1 | Data we haven't used yet: HPA imaging / subcellular localisation, scRNA co-expression atlases, peptide-level evidence, structure, interactomes at scale. |
| **Queen: innovation abduction** | 1 | Reasons from historical innovations (Netflix-prize collaborative filtering, kriging, compressed sensing, BLUP/genomic prediction, word2vec, MSA co-evolution, astronomy missing data, …) to *new* ideas for this problem. Proposes her own ideas and cross-pollinates the other queens through the king's directive. |
| **Workers** | 1 per queen per round | Implement one idea as `src/scp_graphfs_mcp/phasec/<id>.py`, screen it with `benchmarks/phasec.py`, write `docs/phaseC/ideas/<id>.md`. |

## Rules every agent follows

- **Dev data only.**
  - Never load, reference or evaluate on `derks2022`, `wu2026`, `brunner2022` or `khan2023`.
  - Never run a `confirm_v*.py`.
  - The harness asserts this.
- **Never fabricate.**
  - Every number in a doc comes from a JSON file in `results/phaseC/`.
  - Every external data source gets a row in `data/PROVENANCE.md` (URL, licence, SHA-256). The king merges these rows.
  - Workers write the row inside their idea doc.
- **A null control for every idea.**
  - The `--null` function keeps the idea's machinery but destroys its biological content: permuted proteins,
    scrambled labels, etc.
  - Without one, a gain does not count.
- **Screen on the `fast` pool first:** 6 dev targets × 2 masks × 2 seeds.
  - Go to `--pool full` (15 LODO dev targets × 3 seeds) only if `cand_stack` ≤ 0.99 on the fast pool and beats its
    null.
- **File ownership:** workers write only their own idea files. Only the king edits `LEDGER.md`, `KING_LOG.md` and
  `PROVENANCE.md`, and only the king commits.
- **Compute:** the box has 4 CPUs / 15 GB and two workers run at once. Use `--jobs 2`. Keep downloads under about 2 GB
  per idea, in `data/cache/phasec/<id>/` (git-ignored).
- **Promotion:** a holdout is spent only after a full-pool win, and only with a new `PREREGISTRATION_v9.md`. Amit is
  asked before derks2022 (the last reserved holdout) is used.

## Harness

```
python benchmarks/phasec.py --cand scp_graphfs_mcp.phasec.<id>:impute --null scp_graphfs_mcp.phasec.<id>:impute_null \
    --id <id> --pool fast --jobs 2
```

- **Metrics:** `cand_stack` (does the idea add to the champion?), `cand_solo` (could it replace the adapted
  projection?) and `null_stack`. All are relative to the champion; < 1 is better.
- **Candidate signature:** `impute(Xm, ctx) -> full matrix`. The `ctx` keys are documented at the top of
  `benchmarks/phasec.py`.

## Picking up in a new session

1. Read `LEDGER.md` (every idea, its status and its numbers).
2. Read the latest entry in `KING_LOG.md` (current directives per queen).
3. Read `ideas/*.md`.
4. Continue the round loop, or run one queen's next idea by hand with the harness above.

## Run history
- **Rounds 1–3:** 28 agents, Oct 23. Every biological-prior idea failed its null control. The lead is a data-only
  improvement: selected-ψ shrinkage of the target's per-protein mean and log-sd, at full-pool 0.9650 vs the champion.
  Details are in `KING_LOG.md`.
- **Rounds 4–6:** launched from the king's round-3 directives (`techhyper`, `mnarjoint`, `subthr`, `cshrink`). From
  round 4 on, every idea is credited against the `cshrink` base, not the bare champion.
- **To continue in a new session:** take the latest `KING_LOG.md` directives as each queen's next task and follow the
  same queen → worker → king loop.
