# Holdout v9 datasets (built 2026-11-07, not yet used by any method)

These six matrices are **holdouts**. Do not run any imputation, kNN, low-rank, projection, benchmark or confirm
script on them until a v9 preregistration is committed. They are written only to `data/cache/holdout/<name>.npz`
(keys `X`, `accs`, `cells`, `labels`) and `data/datasets/<name>_holdout.h5ad` (`obs['cell_type']` = label). Nothing
was written to `data/cache/corpus/`, the reference corpus. The shared writer `scripts/holdout_ext/_holdout_common.py` refuses any path there.

Processing follows the corpus `save()` in `scripts/corpus_ext/build_woo2025_PC3.py`, applied in this order:
1. Linear intensity becomes log2, and values <= 0 or missing become NaN.
2. Each protein group is identified by its first canonical UniProt accession (`ids.group_accessions`). Groups with no accession are dropped. When two groups share an accession, the most-observed one is kept.
3. Cells with fewer than 50 observed proteins are dropped.
4. Each cell is median-centred.
5. Proteins never observed in a kept cell are dropped.

Only the QC counts below were computed. No correlation or structure statistic was computed.

## Provenance

All three projects were confirmed through `https://www.ebi.ac.uk/pride/ws/archive/v3/projects/<PXD>` on 2026-11-07. Each is Homo sapiens (two also list Mus musculus), label-free DIA single-cell proteomics, with licence **CC0** ("Creative Commons Public Domain (CC0)"). File URLs come from the PRIDE files API (`/projects/<PXD>/files`). Local file sizes equal the PRIDE `fileSizeBytes`. PRIDE gives no checksums, so the SHA-256 values below were computed locally.

| Dataset | PXD | PRIDE title | Publication | Raw file URL | SHA-256 | Builder | Outputs |
|---|---|---|---|---|---|---|---|
| karagach2025_S1 | PXD058457 | High throughput single-cell proteomics of in-vivo cells (lab PI T. Geiger; published 2025-07-09) | Karagach et al. 2025, Mol Cell Proteomics, doi:10.1016/j.mcpro.2025.101018 (linked in PRIDE; Crossref-confirmed) | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD058457/S1_rehyd_cover_384.pg_matrix.tsv | `60715ab23c38a24bd0746fdd94dfb1d9c71082a42622b652048dcffb07641b1a` | `scripts/holdout_ext/build_karagach2025_S1.py` | `data/cache/holdout/karagach2025_S1.npz`, `data/datasets/karagach2025_S1_holdout.h5ad` |
| karagach2025_S2 | PXD058457 | same | same | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD058457/S2_Automation.pg_matrix.tsv | `4cdbf497294d88762ebf3b971bca088801c2bd553be858971a3b5377043ab5be` | `scripts/holdout_ext/build_karagach2025_S2.py` | `data/cache/holdout/karagach2025_S2.npz`, `data/datasets/karagach2025_S2_holdout.h5ad` |
| (SDRF, both above) | PXD058457 | | | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/07/PXD058457/Karagach_Smollich_Resubmision.sdrf.tsv | `5db8a6af13ef9e5dfd41baec5ee3468c8ef087909e60eb6073416495ee28678c` | | |
| px061710_GSTP1KO | PXD061710 | Overcoming preservation challenges to enable single-cell proteomics of fixed cell and tissue samples with retained proteome integrity (lab PI A. von Kriegsheim; published 2025-05-23) | **Not linked in PRIDE.** A Crossref title match gives Makar et al. 2025, J Proteome Res, doi:10.1021/acs.jproteome.5c00268, and the bioRxiv preprint doi:10.1101/2025.03.10.642380. | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD061710/20230613_GSTP1KO_Figure1_search_engine_results.tsv | `740e3c618cf5897c3099b71cffb08995806aa32a7c8874da0c34c49f9214abe7` | `scripts/holdout_ext/build_px061710_GSTP1KO.py` | `data/cache/holdout/px061710_GSTP1KO.npz`, `data/datasets/px061710_GSTP1KO_holdout.h5ad` |
| px061710_MDAfix | PXD061710 | same | same | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD061710/20240709MDA_fresh_fixed_24hrs_Figure3_search_engine_results.tsv | `72d750b7847a0af2978019f37e1002251e23d3f9fa2efd4e1c31cfc3d434e5fb` | `scripts/holdout_ext/build_px061710_MDAfix.py` | `data/cache/holdout/px061710_MDAfix.npz`, `data/datasets/px061710_MDAfix_holdout.h5ad` |
| px061710_RKOfix | PXD061710 | same | same | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD061710/20240809RKOFreshFixed20min_Figure4_search_engine_results.tsv | `2bb0dd829c3137c2fa54c7d1d843169229aeff9ea9146d5764b989759b36a03d` | `scripts/holdout_ext/build_px061710_RKOfix.py` | `data/cache/holdout/px061710_RKOfix.npz`, `data/datasets/px061710_RKOfix_holdout.h5ad` |
| (SDRF, all three above) | PXD061710 | | | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD061710/PXD061710_community_annotated.sdrf.tsv | `160cd973fe524b587654d3738b63798b4fada60c894b47eab1ca5da2792f6402` | | |
| px077946_H2O2 | PXD077946 | Refining liquid chromatography conditions across platforms to democratize single-cell proteomics (lab PI H. Yu; published 2026-09-14) | **Not linked in PRIDE.** A Crossref title match gives Wang et al., bioRxiv preprint doi:10.64898/2025.12.16.694739 (posted 2025-12-19). | https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD077946/diann_singlecell.zip | zip `0bd4c3b9bbd3ad8e94acf7e3ef9dba4330dd1c25667630f177539b5bc5714a05`. Members: `report.unique_genes_matrix.tsv` `aeb2df163cedc56543f02ada9b99d6bfcda5edf3f36c7fbe418c72351c33d197`, `report.pr_matrix.tsv` `674317e5472390f1a9485ff8bd4bdcd46f5126eb3f3524f7d177126ab55c244d`, `report.stats.tsv` `43dfa61592e2417af8a8bac3527db6aa31f147cfbe4734cf59eca42a677ee9e2`. The extracted files are byte-identical to the zip members. | `scripts/holdout_ext/build_px077946_H2O2.py` | `data/cache/holdout/px077946_H2O2.npz`, `data/datasets/px077946_H2O2_holdout.h5ad` |

Other files downloaded under `data/cache/corpus_raw/`, not built, with their SHA-256:
- `karagach2025/S3_1536.pg_matrix.tsv`: `2ca716913cabb7666faa2b681d91fc3a542898dc24ed0b36f84a5a7ea912d9be`
- `karagach2025/S4_fixation.pg_matrix.tsv`: `da59f91b3c09142d8bdee42ff2456bc33da82e80e90629843910714ebac0cf58`
- `karagach2025/S6_Macrophages.pg_matrix.tsv` (mouse, unusable): `d3c64e1c819ea5477461f93ffa5057742da9920140fe0ae1fce7dc94a12c69d9`
- `px061710/20241022RKO_Fresh_frozen_Figure2_search_engine_results.tsv`: `e279a3b1fcf512db4c305a38843097e021d2c56f69a50c2ff065b07fcad64512`
- `px061710/20241221_PDAC_Figure5_search_engine_results.tsv` (mouse KPC, unusable): `8f7061d7ed98982a1a965b7a954158daeb5f3ed98817d5d7f058908a8a743ca2`

All of these are at the same PRIDE directory as their siblings.

## QC (counts only)

Corpus overlap is the number of accessions shared with the union of `accs` over all 25 `data/cache/corpus/*.npz` files. That union holds 15,154 accessions.

| Dataset | Instrument / acquisition / software | Cells | Proteins | Missing | Proteins in >= 50% cells | Corpus overlap | Labels (count) | Flags |
|---|---|---|---|---|---|---|---|---|
| karagach2025_S1 | timsTOF SCP, Slice-PASEF DIA, Evosep ZOOM40, CellenONE; DIA-NN 1.9.1 | 258 (283 runs; 25 with < 50 proteins dropped) | 3687 | 0.334 | 2632 | 3671 / 3687 | 19 sample-prep variants (technical, all HeLa): MM50_rehyd25 34, MM200_rehyd100 20, MM100_rehyd50 19, MM1000_open 14, MM500_open 14, MM200_open 13, MM1000_cover 12, MM1000_rehyd200 12, MM50_cover 12, MM50_rehyd100 12, MM50_rehyd50 12, MM100_open 12, MM100_cover 11, MM500_cover 11, MM500_rehyd200 11, MM500_rehyd250 11, MM50_open 11, MM200_cover 9, MM1000_rehyd500 8 | HeLa only; the labels are prep conditions, not biology |
| karagach2025_S2 | as above (automated Bravo prep) | 39 (40 runs; 1 dropped) | 3704 | 0.143 | 3327 | 3683 / 3704 | Slice-PASEF frames: 1F 10, 2F 10, 2F_Demichev 7, 4F 12 | **< 60 cells**; same lab and cell line as S1 |
| px061710_GSTP1KO | timsTOF SCP, diaPASEF, CellenONE; Spectronaut (`PG.MS2Quantity`) | 78 | 3899 | 0.425 | 2231 | 3885 / 3899 | RKO_sgCTRL 39, RKO_sgGSTP1 39 | |
| px061710_MDAfix | timsTOF SCP, diaPASEF; Spectronaut (`PG.Quantity`) | 76 | 3325 | 0.377 | 2091 | 3300 / 3325 | MDA_fresh 19, MDA_PFA0.1_24h 19, MDA_PFA0.3_24h 19, MDA_PFA1_24h 19 | |
| px061710_RKOfix | timsTOF SCP, diaPASEF; Spectronaut (`PG.Quantity`) | 75 | 5088 | 0.411 | 2904 | 4988 / 5088 | RKO_fresh 18, RKO_PFA0.1_20min 19, RKO_PFA0.3_20min 19, RKO_PFA1_20min 19 | same lab and cell line (RKO) as GSTP1KO |
| px077946_H2O2 | timsTOF HT, dia-PASEF, Evosep ("WhisperZoom"); DIA-NN (version not checked) | 40 | 3268 | 0.450 | 1763 | 3259 / 3268 | 293T_ctrl 7, 293T_H2O2 6, HCT116_ctrl 5, HCT116_H2O2 8, K562_ctrl 7, K562_H2O2 7 | **< 60 cells** |

## How non-single-cell runs were excluded

- **karagach2025.** Every run in S1 and S2 is named `HeLa_SC_` or `Hela_SC_` and is annotated as HeLa in the SDRF. The matrices contain no blanks or bulk runs.
- **px061710.** Columns were joined to the SDRF. The builders keep runs with `cells per well` = 1, which are the `_SC_` runs, and assert that count. They drop blanks and 2-, 5-, 10- and 50-cell runs (18 in GSTP1KO, 20 each in MDAfix and RKOfix).
- **px077946.** All 40 runs come from the archive named `diann_singlecell.zip`, from run directory `20250117h2o2sc`. The 3 ng runs are in a separate archive, `diann_3ng.zip`, which was not used. There is no SDRF, so single-cell status rests on the archive and directory names.
- **px077946 protein IDs.** The archive has no pg_matrix. Values come from DIA-NN's `report.unique_genes_matrix.tsv`. Each gene is mapped to the first accession of the Protein.Group of its proteotypic precursors in `report.pr_matrix.tsv`. TMPO is the one gene with two groups; it takes the most frequent one.

## Duplicate and re-analysis check

- PXD058457, PXD061710 and PXD077946 do not appear in `data/PROVENANCE.md`.
- No file in the repo referenced them before this build.
- None of them is a re-analysis of a corpus dataset.
- Some cell lines also appear in the corpus, from other labs and projects:
  - HeLa in brunner2022, cong2020AC and liang2020_hela
  - K562 in hu2023_K562 and williams2020_tmt

  These are not duplicates.
- **Within-holdout dependence:** the two karagach sets share a lab, a cell line and an instrument. So do the three px061710 sets; GSTP1KO and RKOfix are both RKO. Treat them as clustered by project, not as six independent datasets.
