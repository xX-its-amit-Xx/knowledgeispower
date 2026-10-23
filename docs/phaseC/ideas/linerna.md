# linerna: within-line scRNA programmes as a second-stage factor block on the cshrink residual

- **Queen:** new data modalities (round 5, directive `linerna`)
- **Status:** **dead** (fast, 4 targets). v1 chose gamma = 0 in 32/32 calls, so it is identical to the base (kill rule).
  v2 (one revision, a column-scale fix) has cand_stack 0.9668 vs null 0.9661, perm 0.9669 and base 0.9668: no margin
  over anything, and the mismatched-line null is slightly *better*.
- **Code:** `src/scp_graphfs_mcp/phasec/linerna.py`. Results: `results/phaseC/linerna_fast.json`,
  `linerna_perm_fast.json`, `linerna_base_fast.json`, `linerna_v2_fast.json`, `linerna_v2_perm_fast.json`.
  Per-call logs: `data/cache/phasec/linerna/log_*.jsonl`. Unit checks: `data/cache/phasec/linerna/unit_check.json`.
  Mapping: `data/cache/phasec/linerna/mapping_report.json`.

## Idea

After cshrink, part of each cell's leftover residual may come from programmes that vary from cell to cell inside one
cell line: cell cycle, stress, and monocyte-to-macrophage state. Single-cell RNA of the same line records the gene-gene
co-variation of those programmes. A 20-factor model from that line's scRNA, mapped to proteins, should then predict the
missing part of the residual from its observed part. It should do this better than the same model from a different
line (PC9), better than the same model with permuted protein rows, and better than no second stage (gamma = 0, which is
exactly cshrink).

## Audit: which dev targets have a usable matched line

Steps 1 to 3 are the queen's audit, copied here as given. I re-downloaded only the files in the provenance table and
re-checked their SHA-256. I did not re-check the Kinker or SCP542 claims beyond downloading the pool-composition file.

1. **The Kinker atlas is verified as a source:** GEO GSE157220 (PMID 33128048). Kinker GS, Greenwald AC, Tal R et al.,
   "Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity", *Nat Genet* 2020,
   52(11):1208-1218, doi:10.1038/s41588-020-00726-6. It covers 198 lines, 53,513 cells and 22,722 genes.
2. **The Kinker atlas cannot be used this round:**
   - Its pool composition lists PC3_PROSTATE, the only dev-target line it contains (woo2025_PC3).
   - It lists no haematopoietic line, no U-2 OS and no WM989.
   - The CPM file (1,778,171,720 B) has cell barcodes only.
   - The per-cell line labels (SCP542 `Metadata.txt`) return HTTP 401 without a Single Cell Portal login.
   - So zero dev targets are usable from Kinker.
3. **Replacement single-line sources:**
   - U-2 OS: GSE146773 (Mahdessian et al., *Nature* 2021, 590:649-654, doi:10.1038/s41586-021-03232-9).
   - U-937 and PC9: GSM3972674 in GSE134842 (Aissa et al., *Nat Commun* 2021, 12:1628,
     doi:10.1038/s41467-021-21884-z).
4. **Target-to-line map** (my checks in the "verification" column):

| target | line(s) | match | verification |
|---|---|---|---|
| krull2024 | U-2 OS | exact | `scripts/corpus_ext/build_krull2024.py` header and the ttprior metadata (Cellosaurus CVCL_0042) |
| specht2021 | U-937 monocytes + PMA macrophages | exact | Not named in the scpdata `specht2019v3` manual (`man/specht2019v3.Rd` only says "macrophages and monocytes"), nor in the PubMed abstract. **It is named in the full text** (Europe PMC PMC7839219 fullTextXML): "indicate that U-937-derived macrophages are more heterogeneous than their precursor monocytes", and the benchmark compares "embryonic kidney cells (HEK-293) or … monocytes (U-937 cells)" |
| gregoire2023_mixCTRL | THP-1 / U-937 | partial | scpdata `data.R` (devel): "two monocytes cell lines (THP1 and U937) as well as controled mixtures of both and macrophage-like cells produced upon PMA treatment" |
| leduc2022 | melanoma + U-937 | partial | scpdata `data.R` says "melanoma cells and monocytes"; this target's labels are Monocyte 775 / Melanoma 768. That the monocytes are U-937 is the queen's audit claim, which I did not re-verify |
| woo2025_PC3 | PC-3 | none usable | Kinker has it but with no labels (see step 2) |
| petrosius2023, schoof2021 | OCI-AML8227 | none found | |
| williams2020_tmt | MOLM-14 / K562 / CMK | none | |
| the rest | primary tissue | n/a | |

There are 4 screenable targets: enough for a fast screen (at least 2), not enough for the full pool (at least 5).

5. **Kill check on mapping.** All four targets have at least 30% of proteins mapped
   (`mapping_report.json`), so this kill rule did not fire:

| target | proteins | matched line, mapped | PC9 (null), mapped |
|---|---|---|---|
| krull2024 | 1411 | U2OS 1250 (0.8859) | 862 (0.6109) |
| specht2021 | 3042 | U937 1180 (0.3879) | 1326 (0.4359) |
| gregoire2023_mixCTRL | 1878 | U937 1170 (0.6230) | 1242 (0.6613) |
| leduc2022 | 2844 | U937 1499 (0.5271) | 1672 (0.5879) |

## Mechanism (as built)

**A. RNA factor loadings per line**, cached in `data/cache/phasec/linerna/<line>_pcs.npz`.

- **U2OS:** `GSE146773_Counts.csv.gz` gives 1,152 cells. All 42,728 columns are Ensembl ids, so no non-gene column
  was dropped.
- **Drop-seq split (U937 and PC9):** each cell is assigned by its myeloid marker sum Mk vs its epithelial marker sum Ek.
  - U937 = 658 cells; PC9 = 1,445 cells; 531 cells dropped.
  - The audit's "457 ambiguous" does not add up: 2,634 − 658 − 1,445 = 531.
  - The marker genes are removed before the PCA. Only 17 were removed because CLDN4 is not in the file.
- **PCA per line:** log1p(1e4 · count / total), keep genes detected in at least 5% of the line's cells, z-score each
  gene, truncated SVD with k = 20. L = v·s/√(n−1); h = clip(1 − ΣL², 0.05, 1).

| line | cells | genes kept | variance explained by 20 PCs | mean h | PC1, PC2, PC3 |
|---|---|---|---|---|---|
| U2OS | 1152 | 17616 | 0.1513 | 0.8487 | 0.0339, 0.0236, 0.0108 |
| U937 | 658 | 4403 | 0.1089 | 0.8911 | 0.0380, 0.0128, 0.0057 |
| PC9 | 1445 | 5348 | 0.1023 | 0.8977 | 0.0415, 0.0239, 0.0047 |

- **Gene-to-protein map.** U2OS uses HPA Ensembl → UniProt with an HGNC `ensembl_gene_id` fallback. Symbol lines use
  HPA `Gene` → UniProt with an HGNC `symbol` fallback.
- **Protein rows:** a protein's row is the mean of its mapped genes' rows. Unmapped proteins get Lp = 0, h = 1.

**B. `_fit`** is a copy of `nfbase.adapt_project_mb` that also returns Z and P = UVᵀ.

**C. `impute`:**
- ψ is selected exactly as cshrink does (`TS.selection_mask`, `cshrink._score` over `cshrink.GRID`, tag `'c'`).
- On the same 50/50 inner mask:
  - E = Z − P on observed entries, divided per column by c_j.
  - Per cell: f_i = (Lp_oᵀH_o⁻¹Lp_o + I)⁻¹ Lp_oᵀH_o⁻¹E_io, with the cell skipped if |o| < 20; then ehat = Lp_m f_i.
  - Y = out + γ·sd·c·ehat, with γ in {0, 0.25, 0.5, 1, 2} and ties going to the smaller γ.
- The final fit on the input uses iters = 10.
- Entry points:
  - `impute`: matched line
  - `impute_null`: PC9
  - `impute_perm`: matched line, Lp/h rows permuted among mapped proteins, rng seed + 7
  - `impute_base`: γ = 0
  - `impute_v2`, `impute_null_v2`, `impute_perm_v2`: the same three arms with the v2 scale fix

**Unit checks** (`unit_check.json`, krull2024 MCAR seed 100):
- `_fit` equals `nfbase.adapt_project_mb` bitwise: true.
- `impute_base` vs `cshrink.impute`: max |diff| on missing entries = 0.0. Pass.

**v1 → v2 (the one allowed revision, a scale bug).**
- **The bug:** in v1, c_j = the observed residual sd with a 1e-3 floor. Columns with 1-2 observed residuals get
  c_j ≈ 1e-3, so E/c explodes and corrupts f_i.
- **Smoke diagnostic** (krull2024 inner mask, psi 2): |ehat| was at 95th percentile 10.6 and 99.9th 85.0. Raising γ
  from 0 to 0.25 doubled the inner RMSE.
- **v2 fix:** c_j = the observed sd only if the column has at least 5 observed residuals, otherwise 1; floor 0.1. This
  rule was fixed from the diagnosis, not tuned on scores.
- **After the fix:** the same split gives |ehat| at 95th percentile 0.43 and 99.9th 0.92.
- **The base is shared:** γ = 0 never touches c, so `linerna_base` is the base for both versions.

## Null control

- **`impute_null` (mismatched line):** the same machinery (cshrink ψ, inner-mask γ, k = 20, the same map) with PC9
  covariance. PC9 is lung adenocarcinoma and matches no dev target. For the U-937 targets it comes from the same
  Drop-seq run, so only the line changes.
- **`impute_perm` (permuted protein rows):** keeps the spectrum and sparsity but scrambles which protein loads on which
  programme.
- **`impute_base` (γ = 0):** equals cshrink.
- **Confound:** for krull2024 the matched line is Smart-seq2 but PC9 is Drop-seq, so the null also changes platform
  there. The perm arm is the clean control for that target.

## Results

All values are rmse_all `mean_rel_to_champion` (in-complex after the slash), fast pool restricted to
`krull2024,specht2021,gregoire2023_mixCTRL,leduc2022`, 2 masks × 2 seeds = 8 target-masks. Per-mask, exact-match
(krull2024, specht2021) and partial-match (gregoire2023_mixCTRL, leduc2022) values are means of the JSON
`per_target_mask` rows.

| run | cand_stack | cand_solo | cand_raw | null_stack | stack weight | target-masks ≤ 0.99 | runtime (s) |
|---|---|---|---|---|---|---|---|
| `linerna` (v1) | 0.9668 / 0.9688 | 0.9668 / 0.9696 | 0.9648 / 0.9668 | 0.9668 / 0.9688 (PC9) | 0.7188 | 8/8 | 206.6 |
| `linerna_perm` (v1) | 0.9668 / 0.9688 | 0.9668 / 0.9696 | 0.9648 / 0.9668 | 0.9668 / 0.9688 (PC9) | 0.7188 | 8/8 | 200.8 |
| `linerna_base` (γ = 0) | 0.9668 / 0.9688 | 0.9668 / 0.9696 | 0.9648 / 0.9668 | 0.9716 / 0.9747 (cshrink random ψ) | 0.7188 | 8/8 | 184.1 |
| `linerna_v2` | **0.9668 / 0.9686** | 0.9668 / 0.9693 | 0.9651 / 0.9667 | **0.9661 / 0.9682** (PC9) | 0.7188 | 8/8 | 118.3 |
| `linerna_v2_perm` | 0.9669 / 0.9694 | 0.9666 / 0.9695 | 0.9647 / 0.9668 | 0.9661 / 0.9682 (PC9) | 0.7344 | 8/8 | 122.1 |

Per mask and per match group (rmse_all):

| run | cand_stack MCAR / MNAR | cand_raw MCAR / MNAR | null_stack MCAR / MNAR | cand_stack exact / partial | null_stack exact / partial |
|---|---|---|---|---|---|
| base | 0.9749 / 0.9586 | 0.9780 / 0.9516 | (cshrink null) 0.9777 / 0.9654 | 0.9720 / 0.9615 | 0.9707 / 0.9724 |
| v1 main | 0.9749 / 0.9586 | 0.9780 / 0.9516 | 0.9749 / 0.9586 | 0.9720 / 0.9615 | 0.9720 / 0.9615 |
| v2 main | 0.9756 / 0.9580 | 0.9792 / 0.9511 | 0.9753 / 0.9569 | 0.9716 / 0.9620 | 0.9706 / 0.9616 |
| v2 perm | 0.9754 / 0.9584 | 0.9781 / 0.9514 | 0.9753 / 0.9569 | 0.9722 / 0.9615 | 0.9706 / 0.9616 |

Per target-mask, v2 main cand_stack (null_stack PC9):

| target | MCAR | MNAR |
|---|---|---|
| krull2024 | 0.9872 (0.9837) | 0.9830 (0.9828) |
| specht2021 | 0.9686 (0.9718) | 0.9475 (0.9440) |
| gregoire2023_mixCTRL | 0.9686 (0.9677) | 0.9326 (0.9324) |
| leduc2022 | 0.9781 (0.9779) | 0.9689 (0.9685) |

**γ histograms** (from the per-call logs; each run calls the candidate twice per split, so 32 calls per arm):

| arm | γ = 0 | 0.25 | 0.5 | 1 | 2 |
|---|---|---|---|---|---|
| v1 main | **32** | 0 | 0 | 0 | 0 |
| v1 perm | 32 | 0 | 0 | 0 | 0 |
| v1 null (both runs, 64 calls) | 64 | 0 | 0 | 0 | 0 |
| v2 main | 9 | 13 | 7 | 3 | 0 |
| v2 perm | 21 | 6 | 4 | 1 | 0 |
| v2 null (both runs, 64 calls) | 20 | 8 | 18 | 16 | 2 |

- **Mean inner-mask γ score curve, v2,** relative to γ = 0 at γ = 0 / 0.25 / 0.5 / 1 / 2:
  - main: 1 / 0.9998 / 1.0004 / 1.0045 / 1.0234
  - null: 1 / 0.9992 / 0.9990 / 1.0002 / 1.0092
  - perm: 1 / 1.0002 / 1.0013 / 1.0061 / 1.0264
- **Inner-mask Pearson(ehat, held-out residual), v2, by target (krull / specht / gregoire / leduc):**
  - main: 0.0109 / 0.0100 / 0.0101 / 0.0087
  - null (PC9): 0.0526 / 0.0322 / 0.0234 / 0.0050
  - perm: 0.0080 / 0.0001 / −0.0215 / −0.0152
- **Mean |o|** (observed mapped proteins per cell, inner mask, v2 main): krull 367.1, specht 328.7, gregoire 446.4,
  leduc 560.8.

**The bars (v2):**
- cand_stack 0.9668 ≤ 0.99: yes.
- Margin below the null (0.9661): −0.0007 (the candidate is worse). Below perm (0.9669): 0.0001. Below the base
  (0.9668): 0.0000. All are far short of 0.005.
- MCAR cand_raw 0.9792 vs base 0.9780: 0.0012 worse, within 0.005.
- **Fails.** v1 also hits the kill rule (γ = 0 in 100% of calls, at least 75% kills). The full pool is not run: the
  fast bars fail, and there are fewer than 5 targets.

## Honest read

- **v1 never switched on.** Its scale bug made every γ > 0 catastrophic on the inner mask, so selection returned
  cshrink exactly (`linerna` = `linerna_perm` = `linerna_base` to 4 decimals).
- **v2 switches on in 72% of calls but buys nothing.**
  - The best inner-mask γ improves the score by about 0.02% on average.
  - The matched-line ehat correlates with the held-out residual at about 0.01.
  - The **mismatched PC9 covariance correlates better** (0.05 on krull2024, the target with a platform confound) and
    stacks marginally better (0.9661 vs 0.9668).
- **What the residual shows:** whatever residual structure is left after cshrink is not the within-line programme
  structure seen in scRNA. The small part that any external covariance captures is generic (housekeeping /
  depth-like axes shared across lines), not line-specific.
- **Likely reasons:**
  1. 20 PCs explain only 10-15% of gene variance, so h ≈ 0.85-0.90 and the factor block is weak by construction.
  2. Within-line mRNA-protein co-variation is weak.
  3. The cshrink projection's 20 reference + 10 free factors already absorb the dominant shared axes of these 30-cell
     targets.
  4. Shallow Drop-seq depth (U937 has only 4,403 genes at ≥ 5% detection) leaves specht2021 with 39% of proteins
     mapped.
- **The exact-match targets do not stand out:** krull2024 and specht2021 average 0.9716 vs the null's 0.9706.

## Next step

Do not pursue within-line scRNA covariance as a second-stage residual prior. Getting labelled PC-3 cells (Kinker SCP542)
to reach 5 targets is **not** justified by this screen. If within-line covariance is retried, it should be only as a
Σ_ref swap inside rep's `lwkrige` (an import, not a residual block), and only if lwkrige's own in-domain covariance
shows a gain over its null first.

## Data provenance (rows for the king to merge into `data/PROVENANCE.md`)

Licence for all GEO files: "GEO public record; licence not stated in the series". See
https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html (HTTP 200 on 2026-10-23; I did not read or quote its text).

| local file | URL | content | date | SHA-256 | bytes |
|---|---|---|---|---|---|
| `data/cache/phasec/linerna/GSE146773_Counts.csv.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146773/suppl/GSE146773_Counts.csv.gz | FUCCI U-2 OS Smart-seq2 counts, 1,152 cells × 42,728 Ensembl genes (Mahdessian et al. 2021, *Nature* 590:649-654, doi:10.1038/s41586-021-03232-9, PMID 33627808) | 2026-10-23 | `e2f2e79dba3ae39f3e9528a96f12caf4d240e42f3ccde73a47fe12a9dbc0c563` | 27,824,191 |
| `data/cache/phasec/linerna/GSM3972674_PC9U937.dge.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/GSM3972674/suppl/GSM3972674_PC9U937.dge.txt.gz | Drop-seq PC9 + U-937 mix, 25,132 genes × 2,634 cells (GSE134842; Aissa et al. 2021, *Nat Commun* 12:1628, doi:10.1038/s41467-021-21884-z, PMID 33712615) | 2026-10-23 | `e05217a5feff79bb6dee79ea11b9fc53e3f350b22198539af65640e122b3344e` | 3,195,381 |
| `data/cache/phasec/linerna/GSE157220_Pool_composition.xlsx` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157220/suppl/GSE157220_Pool_composition.xlsx | Kinker et al. 2020 pool composition (*Nat Genet* 52:1208-1218, doi:10.1038/s41588-020-00726-6, PMID 33128048); audit only, not used for modelling | 2026-10-23 | `7f71f3d7223a9d92ed51c810fafc0e150dd5be0b362e0e3f49f973b949ee62b0` | 15,219 |

**Verification-only fetches (not cached in the repo):**
- scpdata `man/specht2019v3.Rd` (https://raw.githubusercontent.com/UCLouvain-CBIO/scpdata/master/man/specht2019v3.Rd).
- scpdata `R/data.R` (https://raw.githubusercontent.com/UCLouvain-CBIO/scpdata/devel/R/data.R).
- Europe PMC full text of PMC7839219 (https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7839219/fullTextXML).

**Gene-protein mapping** uses the existing `data/cache/features/proteinatlas.tsv.zip` (already in PROVENANCE.md) and
the HGNC table via `scp_graphfs_mcp.ids`. Total cache in `data/cache/phasec/linerna/`: 36 MB.
