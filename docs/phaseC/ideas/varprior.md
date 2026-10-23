# varprior: protein turnover and complex buffering as the log-sd shrinkage target

- **Queen:** biological priors (round-3 directive in `KING_LOG.md`, pivot from tissue to protein dynamics).
- **Status (worker's recommendation):** **dead** at the fast screen. It clears 0.99 vs the champion, but it does not
  clear the required 0.005 margin over its null or over the shrink-to-constant baseline.
  - Main run (`varprior`): cand_stack 0.9679 (all) / 0.9686 (in-complex). The null is 0.9718 / 0.9728, a margin of
    0.0039. The const baseline is 0.9686 / 0.9692, a margin of 0.0007.
  - Revision (`varprior_v2`, no observed-fraction trend): cand_stack 0.9682 / 0.9689. Its null is 0.9686 / 0.9691, a
    margin of 0.0004. The const baseline is 0.9686, also a margin of 0.0004.
  - The full pool was not run.
- **Code:** `src/scp_graphfs_mcp/phasec/varprior.py`. It imports `ref_baselines`, `baselines`, `adapt_project_mb`,
  `LAM_V`, `EXTRA`, `_load` and `FORBIDDEN` from `nfbase`; nothing is copied.
- **Results** (in `results/phaseC/`): `varprior_fast.json`, `varprior_const_fast.json`, `varprior_fonly_fast.json`,
  `varprior_v2_fast.json`.
- **Cache** (in `data/cache/phasec/varprior/`, not committed):
  - `mathieson2018_MOESM5.xlsx`, `halflife_by_symbol.json`, `halflife_meta.json`, `mapping_meta.json`
  - `<target>_cov.npz`, `ncells.json`
  - `diagnostics.json`, `unit_check.json`, `log_*.txt`

## Idea
The champion standardises each target protein with its 30-cell nanstd, and that sd_j multiplies every imputed value of
protein j. nfbase showed that shrinking log-sd at psi = 4 is a large lever: post-hoc sdonly_psi4 reached 0.9832.

The hypothesis was that cell-to-cell variability is partly a protein property:
- Long-lived proteins average out synthesis bursts.
- Orphan complex subunits are degraded, so complex members are buffered.

If that holds, a shrinkage *target* predicted from measured half-life and CORUM membership should beat shrinking to a
constant.

## Mechanism (as built, following the bio spec)
**0. Data.**
- Mathieson 2018 Supplementary Data 3 (MOESM5), sheet 1 ('protein half lives high qual'). It is parsed with
  zipfile + xml.etree; openpyxl is not installed.
- The 8 human `half_life` columns are: B cells, NK cells, hepatocytes and monocytes, replicates 1 and 2 each. Mouse
  columns are excluded. Cells with t='e' are NaN.
- hl = nanmedian per gene. This gives 6,855 genes, median 72.539 h (`halflife_meta.json`).
- `ids.uniprot_for_symbol` resolves 6,454 of the 6,855 symbols (94%) to 6,473 accessions (`mapping_meta.json`). Every
  accession of a symbol is used.
- CORUM comes from `complexes.load_complexes(("corum",))`.

**1. Covariates per target protein (first accession via `T._target_accs`).**
- h = log2 hl. Where hl is missing, h is set to the median of the covered h, and mh = 1.
- cx = in any CORUM complex.
- sz = log of the median member count of the protein's complexes.
- h and sz are z-scored over the target's p proteins.
- Coverage of h (`diagnostics.json`):

  | Target | Coverage of h |
  |---|---|
  | px063590_ADbrain | 0.597 |
  | leduc2022 | 0.719 |
  | sadiku2025 | 0.723 |
  | weiss2026_liver | 0.758 |
  | gregoire2023_mixCTRL | 0.786 |
  | petrosius2023 | 0.797 |

- CORUM coverage runs from 0.370 to 0.522.

**2. Reference model.**
- The reference is `nfbase.ref_baselines(target, names)`. It is LODO and excludes the target and all holdouts; the
  code also asserts this.
- Response: r_dj = l_dj − median_j l_dj.
- f_dj = (1 + 1/(2 vl_dj)) / N_d. N_d comes from `nfbase._load(d)`, with a holdout assert, cached in `ncells.json`.
- Design: x = [1, h, mh, cx, sz, f, f²].
- Fit: pooled ridge, λ = 1, intercept unpenalised. It is fitted once per target (lru_cache).

**3. `impute`.**
- f_j is the target's own observed fraction.
- pl = x·β.
- py = 0, so location is shrunk to a constant.
- `nfbase.baselines(Xm, py, pl, psi=4)`, then `adapt_project_mb(lam_v=100, extra=10)`.

**4. Variants.**
- `impute_const`: pl = 0, the required baseline.
- `impute_fonly`: x = [1, f, f²], attribution only.
- `impute_v2`: the one allowed revision, x = [1, h, mh, cx, sz], with no f trend.

**Unit check (`unit_check.json`).** At psi = 0, `impute` reproduces `transfer.adapt_project` on MCAR seed 100 of all 6
fast targets. The max-abs difference is 1.1e-15 to 5.6e-15, so all six pass < 1e-10.

## Null control
`impute_null` permutes the covariate rows (h, mh, cx, sz) jointly across the target's p proteins.
- rng: `default_rng(4242 + index of the target in sorted(FAST))`.
- The same permutation is applied in the reference fit and in the target prediction, and β is refitted.
- f (reference f_dj and target f_j) is kept.
- `impute_v2_null` does the same for the v2 design.
- `varprior_const` and `varprior_fonly` were run with `--null impute_null`, so their null_stack equals the main run's.

## Diagnostics before the harness (`diagnostics.json`, dev reference only)
Pooled ridge β for the full model on the 6 fast targets, ordered [1, h, mh, cx, sz, f, f²]:

| Coefficient | Range of β | Typical SE |
|---|---|---|
| h | −0.0026 to 0.0066 | ≈ 0.002 |
| mh | 0.033 to 0.155 | ≈ 0.006 |
| cx | −0.0011 to 0.0205 | ≈ 0.006 |
| sz | −0.018 to −0.037 | ≈ 0.0025 |
| f | 0.41 to 0.58 | not reported |
| f² | −0.39 to −0.56 | not reported |

- **Expected signs, h < 0 and cx < 0: not met.**
  - The half-life coefficient is ≈ 0 or slightly positive.
  - The partial Spearman of h with r given f is −0.006 to 0.034 over all proteins, and 0.011 to 0.042 over proteins
    with a measured hl. The null's values are −0.011 to 0.011.
  - The cx marginal partial Spearman is negative (−0.046 to −0.082), but in the joint fit that signal sits in sz
    (larger complexes → lower log-sd), not in cx.
  - The largest bio term is **mh**: proteins *without* a Mathieson half-life have higher log-sd. That is most
    plausibly a low-abundance / detectability proxy, not turnover.
- **LODO R² within the reference.** The full model reaches 0.0037 to 0.0363; f-only reaches −0.0048 to 0.0058; the
  null reaches −0.0051 to 0.0059. The incremental R² from (h, mh, cx, sz) is 0.0031 (leduc2022) to 0.0321
  (sadiku2025). The bio covariates explain at most about 3% of the cross-protein variance in relative log-sd.

## Results (fast pool, rmse_all mean_rel_to_champion; in-complex in brackets)
| id | cand_stack | null_stack | cand_solo | cand_raw | splits ≤ 0.99 | stack weight |
|---|---|---|---|---|---|---|
| `varprior` (full model) | **0.9679** (0.9686) | 0.9718 (0.9728) | 0.9678 | 0.9656 | 12/12 | 0.75 |
| `varprior_const` (pl = 0, required baseline) | 0.9686 (0.9692) | 0.9718 (0.9728) | 0.9686 | 0.9664 | 12/12 | 0.75 |
| `varprior_fonly` (f, f² only) | 0.9720 (0.9731) | 0.9718 (0.9728) | 0.9722 | 0.9705 | 12/12 | 0.71875 |
| `varprior_v2` (bio covariates, no f) | 0.9682 (0.9689) | 0.9686 (0.9691) | 0.9682 | 0.9660 | 12/12 | 0.75 |

For reference, nfbase_psi4's null was 0.9691. The const baseline here is 0.9686, 0.0005 better; the difference comes
from the wider coverage: every base_ok protein is shrunk, not only those the reference covers.

**Pass bar.**
- Main run: cand_stack ≤ 0.99 passes. 0.9679 ≤ const − 0.005 = 0.9636 **fails**. 0.9679 ≤ null − 0.005 = 0.9668
  **fails**.
- v2: fails both margins (0.0004 each).

**Per target (cand_stack, from the JSONs).**
- The main run's real-over-null margin is concentrated on weiss2026_liver (MCAR 0.9651 vs null 0.9859; MNAR 0.9264 vs
  0.9457) and sadiku2025 (MCAR 0.9787 vs 0.9816; MNAR 0.9658 vs 0.9715). On weiss the main run only *matches* const
  (0.9659 / 0.9264).
- The f-only trend hurts weiss2026_liver badly (MCAR 0.9873 vs const 0.9659).
- So the main run's null margin is mostly the null *losing* ground: it keeps the harmful f-trend but has no bio content
  to counter it. That is not the bio prior adding information over constant shrinkage.

## Honest read
- Measured protein half-life does not predict single-cell log-sd in this corpus once observed fraction is accounted for
  (β_h ≈ 0, wrong sign).
- Complex size has a small, consistent negative effect, as buffering predicts. The largest "bio" term is "has no
  Mathieson half-life", which is a detectability proxy.
- Together these explain ≤ 3% of the reference's relative log-sd variance. Fed into the psi = 4 EB fusion, they change
  cand_stack by 0.0004 to 0.0007 relative to shrinking to a constant.
- The scale lever nfbase found is a *shrinkage* effect (variance reduction of a 30-cell sd), not a *target* effect.
  Where the shrinkage target points barely matters, unless it points somewhere wrong (the f-trend on weiss2026_liver).
- Biological dynamics priors therefore add nothing measurable on top of `tgtshrink`-style constant shrinkage.
- A useful by-product for the king: shrink-to-constant over every base_ok protein (`varprior_const`, location + log-sd,
  psi = 4) reaches 0.9686, with 12/12 splits ≤ 0.99 and no external data. This is the data-only `tgtshrink` lever.
  Note that a reference-trained observed-fraction trend makes it *worse* (`varprior_fonly` 0.9720), which is relevant
  to tgtshrink's planned (mean, observed fraction) trend.

## Next step
- Close the half-life line.
- If the bio queen stays on scale, the only covariate with a consistent sign is complex size. It is too weak (|β| ≈
  0.03 on a log-sd scale) to clear 0.005.
- Better levers are probably:
  - the shrinkage *strength* (psi chosen per target on the inner mask, which is tgtshrink's design);
  - a target-side trend, rather than one transferred from the reference.

## Data provenance (for the king to merge into `data/PROVENANCE.md`)
| Source | URL | Licence | SHA-256 | Retrieved | Notes |
|---|---|---|---|---|---|
| Mathieson et al. 2018, *Nat Commun* 9:689, doi:10.1038/s41467-018-03106-1, Supplementary Data (MOESM5) | https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-03106-1/MediaObjects/41467_2018_3106_MOESM5_ESM.xlsx | CC BY 4.0 | e0503f49230680cbf86fbc86db145ba439c03382e584fde8f0d83241c9826515 (3,338,308 bytes; matches the queen's check) | 2026-10-23 | Sheet 'protein half lives high qual' (sheet1.xml); the 8 human half_life columns (B, NK, hepatocyte, monocyte × 2 replicates); per-gene nanmedian; 6,855 genes, median 72.539 h |
| CORUM human complexes | already cached: data/raw/corum_humanComplexes.txt | already in PROVENANCE.md | 5e40556eb59bb767f3396ba8c080f64c777903ffca9102345fdc35662a2e0282 | (existing) | membership and size covariates |
| HGNC complete set | already cached: data/raw/hgnc_complete_set.txt | already in PROVENANCE.md | (existing) | (existing) | symbol → UniProt via `ids.uniprot_for_symbol` |
