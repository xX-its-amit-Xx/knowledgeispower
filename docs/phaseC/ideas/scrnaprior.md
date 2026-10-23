# scrnaprior: HPA single-cell-type RNA programmes as an extra prior factor block in adapt_project

- **Queen:** new data modalities (round-1 directive)
- **Status (worker's recommendation):** **dead** at the fast screen.
  - cand_stack (rmse_all) is 0.9973, which misses the 0.99 bar.
  - The margin over the null (0.9987) is 0.0014, below the 0.005 bar.
  - The full pool was not run. No v2 revision was made: there was no bug or scale fault to fix (see the checks below).
- **Code:** `src/scp_graphfs_mcp/phasec/scrnaprior.py`, with `impute`, `impute_null` and the sanity check `impute_zero`.
- **Results:**
  - `results/phaseC/scrnaprior_fast.json` (the screen that counts)
  - `results/phaseC/scrnaprior_zero_fast.json` (S = 0 sanity check)
  - `results/phaseC/scrnaprior_lams100_fast.json` (single a-priori diagnostic, lam_s = 100; not cherry-picked)

## Idea
The champion's prior holds only proteomic co-regulation from the SCP reference corpus. The HPA single-cell-type RNA
consensus covers 154 cell types. Its leading gene programmes (SVD of z-scored log1p nCPM, gene x cell type) say which
genes rise and fall together across cell identities. We add them as k = 8 extra prior-anchored factors in the adapted
projection, with their own ridge.

## Mechanism (as built, following the spec)
1. **Atlas.** `rna_single_cell_type.tsv` was pivoted to genes x cell types (20,162 genes x 154; missing -> 0) and cached
   as `pivot.npz`. log1p(nCPM) was taken, genes with max nCPM < 1 were dropped (18,915 kept), and each gene was z-scored
   across the 154 cell types. A full SVD followed, keeping k = 8: L = A_8 Sigma_8 / sqrt(154). The 8 components hold
   0.472 of the variance; the leading singular values are 647.0, 500.8, 415.2, 387.5, 340.3, 317.2, 302.7, 274.7.
   Cached as `svd_k8.npz`. Genome-wide, so it does not depend on the target.
2. **Mapping.** Each protein name -> `ids.group_accessions` -> Ensembl genes. Primary map: `proteinatlas.tsv.zip`
   (Ensembl, Uniprot). Fallback: HGNC `ensembl_gene_id` -> `uniprot_ids`. The row is the mean of L over the matched
   genes. Coverage:

   | target | proteins | mapped | fraction |
   |---|---|---|---|
   | petrosius2023 | 2813 | 2793 | 0.993 |
   | leduc2022 | 2844 | 2813 | 0.989 |
   | weiss2026_liver | 4520 | 4483 | 0.992 |
   | gregoire2023_mixCTRL | 1878 | 1842 | 0.981 |
   | px063590_ADbrain | 5140 | 4603 | 0.896 |
   | sadiku2025 | 3039 | 3005 | 0.989 |

3. **Scale.** S is multiplied by c = RMS(nonzero rows of V[warm]) / RMS(nonzero rows of S). The RMS values afterwards
   match exactly (for example, 0.473 vs 0.473 on petrosius2023). Per-names S is cached as `S_<sha1>.npz`.
4. **Adapted projection.** A local copy of `transfer.adapt_project` (`transfer.py` is not edited) with
   Vp = [V_warm (20) | S (8) | 0 (10)]. Noise goes on the free block only, from the same rng call as the original.
   The V-step ridge is D = diag([100]*20 + [30]*8 + [100]*10); lam = 5, iters = 10, seed = 0 (the champion's default).
   - **Equivalence check:** with an empty S block (k = 0), the copy matches `T.adapt_project` exactly. Setup:
     petrosius2023, seed-100 subsample, 20% random mask; max |diff| = 0.0. This was a run-time check, not saved to JSON.
   - Note: the spec said seed = ctx['seed']. I used the default seed = 0, as in the harness's champion and ebshrink. The
     seed only sets the 0.01-scale noise on the free block.
5. k = 8 and lam_s = 30 were fixed a priori and not tuned. `SCRNAPRIOR_LAM_S` is an env override used only for the
   single lam_s = 100 diagnostic.

## Null control
`impute_null` uses the same S, scale, k, lam_s and smask, but the rows of S are permuted among the mapped proteins
(`default_rng(seed + 7)`).
- **Kept:** the spectrum, the scale, the coverage pattern and the extra degrees of freedom.
- **Destroyed:** which protein gets which cell-type programme.

**Secondary check** (`impute_zero`): S = 0 with lam_s = 30, i.e. 8 more ridge-to-zero factors.

## Results (fast pool: 6 targets x 2 masks x 2 seeds; mean ratio vs champion, < 1 better)
| run | metric | cand_raw | cand_solo | cand_stack | null_stack | splits <= 0.99 | mean cand weight | runtime |
|---|---|---|---|---|---|---|---|---|
| scrnaprior (lam_s 30) | rmse_all | 1.0034 | 0.9989 | **0.9973** | 0.9987 | 0 / 12 | 0.2604 | 381.8 s |
| scrnaprior (lam_s 30) | rmse_in_complex | 1.004 | 0.999 | 0.9974 | 0.9988 | | | |
| scrnaprior_zero (S = 0) | rmse_all | 1.0002 | 1.0007 | 0.9998 | — | 0 / 12 | | 738.7 s |
| scrnaprior_zero (S = 0) | rmse_in_complex | 1.0002 | 1.0003 | 0.9997 | — | | | |
| scrnaprior_lams100 (diagnostic) | rmse_all | 0.9998 | 0.9978 | 0.9976 | 0.9989 | 0 / 12 | 0.3854 | 829.6 s |
| scrnaprior_lams100 (diagnostic) | rmse_in_complex | 1.0004 | 0.9979 | 0.9978 | 0.999 | | | |

**Per target, main run (rmse_all; cand_raw / cand_stack / null_stack):**

| target | mask | cand_raw | cand_stack | null_stack |
|---|---|---|---|---|
| gregoire2023_mixCTRL | MCAR | 1.0133 | 1.0004 | 1.0 |
| gregoire2023_mixCTRL | MNAR | 1.0144 | 1.0002 | 1.0003 |
| leduc2022 | MCAR | 1.0156 | 1.0 | 1.0 |
| leduc2022 | MNAR | 1.008 | 1.0001 | 0.9998 |
| petrosius2023 | MCAR | 0.9965 | 0.9947 | 0.9976 |
| petrosius2023 | MNAR | 0.9966 | 0.9958 | 0.9976 |
| px063590_ADbrain | MCAR | 1.018 | 0.9998 | 0.9998 |
| px063590_ADbrain | MNAR | 1.0132 | 0.9991 | 0.9996 |
| sadiku2025 | MCAR | 1.0059 | 1.0001 | 0.9984 |
| sadiku2025 | MNAR | 0.9749 | 0.996 | 1.0011 |
| weiss2026_liver | MCAR | 0.9915 | 0.9903 | 0.9966 |
| weiss2026_liver | MNAR | 0.9932 | 0.9908 | 0.9935 |

Runtimes are inflated by the parallel worker sharing the CPUs. The smoke test (petrosius2023 only, id
`scrnaprior_smoke`) was deleted.

## Honest read
- **The ordering is consistent but tiny.** Real RNA (0.9973) < permuted RNA (0.9987) < zeros (0.9998). The same
  ordering holds at lam_s = 100 (0.9976 vs 0.9989). The real mapping carries a little signal beyond "more factors", but
  the effect is about 0.1–0.3% RMSE, well below both bars. No split reaches 0.99.
- **Where it helps.** Only weiss2026_liver (0.990–0.991) and petrosius2023 (about 0.995), with some help on
  sadiku2025 MNAR. On gregoire2023, leduc2022 and px063590 it does nothing. Its raw predictions are worse than the
  champion (1.01–1.02), so the stack mostly turns it off there.
- **Why it is small (my reading, not tested):**
  - The 20 proteomic factors fitted on about 21 SCP datasets probably already span most of the cell-identity axes
    present in a 30-cell subsample.
  - RNA-protein concordance across cell types is moderate, so the RNA axes are a noisy version of what the proteomic
    prior already holds.
  - The targets with some gain (liver tissue, petrosius) may contain cell-type contrasts that the HPA atlas resolves
    well.
- **Not a bug or a scale fault:**
  - The copy is exact at k = 0.
  - The scale matches the proteomic block.
  - Coverage is 90–99%.
  - A stiffer ridge (lam_s = 100) gives the same result.
  - So the one allowed revision was not used.

## Data provenance (row for the king to merge into data/PROVENANCE.md)
| file | URL | licence | downloaded | SHA-256 | notes |
|---|---|---|---|---|---|
| `data/cache/phasec/scrnaprior/rna_single_cell_type.tsv.zip` | https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip | CC BY-SA 4.0 (Human Protein Atlas) | 2026-10-23 | `2ddd3d90fb050a8efea508792ba7741c510fb5aa3d89f352aab99794ec2c9d64` | 16,346,850 bytes; inner TSV dated 2025-12-12; matches the pre-checked hash |

Derived caches (git-ignored, regenerable): `pivot.npz`, `svd_k8.npz`, `S_<sha1 of names>.npz`. ID mapping uses the
existing `data/cache/features/proteinatlas.tsv.zip` (already in PROVENANCE) and `data/raw/hgnc_complete_set.txt`.

## Next step
- Do not take this to the full pool.
- The signal that survives the null is concentrated in tissue-like targets. A cheaper and more targeted variant could
  restrict the atlas to cell types plausibly present in the target, e.g. the tissue-cluster table
  `rna_single_cell_cluster.tsv.zip`, about 197 MB. But the ceiling looks low.
- The queen's listed fallback (HPA subcellular-localisation imaging as compartment blocks) is a different kind of
  information and is the better next bet for this modality.
