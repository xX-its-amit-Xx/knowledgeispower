#!/usr/bin/env bash
# Download the protein-feature sources used by scp_graphfs_mcp.features (Amendment 5) into data/cache/features/.
# Checksums of the snapshots used are in data/PROVENANCE.md ("current" releases may move on).
set -euo pipefail
cd "$(dirname "$0")/../data/cache" && mkdir -p features && cd features
curl -sS -L -O https://current.geneontology.org/annotations/goa_human.gaf.gz
curl -sS -L -O https://www.proteinatlas.org/download/proteinatlas.tsv.zip
curl -sS -L -O https://pax-db.org/downloads/latest/datasets/9606/9606-WHOLE_ORGANISM-integrated.txt
curl -sS -L -O https://stringdb-downloads.org/download/protein.physical.links.v12.5/9606.protein.physical.links.v12.5.txt.gz
curl -sS -L -o uniprot_human.tsv "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:9606+AND+reviewed:true&fields=accession,length,mass,cc_subcellular_location,keyword&format=tsv"
sha256sum *

# --- Amendment 9 representations (agent: prott5/pdb/pharos)
# Raw files -> data/cache/reps/<name>/; derived caches (pca128.npz, pairs.npz, targets.tsv.gz) are rebuilt on first use.
(
# cwd is still data/cache/features from the section above ($0 may be relative, so don't re-resolve it)
if [ "$(basename "$PWD")" = features ]; then cd ..; else cd "$(dirname "$0")/../data/cache"; fi
mkdir -p reps/prott5 reps/pdb_costructure reps/pharos && cd reps
UA='scp-graphfs-mcp/0.1 (research)'
(cd prott5 && curl -sS -L -A "$UA" -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/UP000005640_9606/per-protein.h5 \
  && curl -sS -L -A "$UA" -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/embeddings/UP000005640_9606/RELEASE.metalink)
(cd pdb_costructure && curl -sS -L -A "$UA" -O https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz)
# Pharos/TCRD: one GraphQL request returns all ~20k targets (~20 s).
(cd pharos && curl -sS -m 900 -A "$UA" -H 'Content-Type: application/json' -X POST https://pharos-api.ncats.io/graphql -o pharos_targets.json \
  -d '{"query":"{ dbVersion targets { count targets(top:25000, skip:0) { uniprot sym tdl fam novelty ligandCounts { name value } } } }"}')
sha256sum prott5/per-protein.h5 prott5/RELEASE.metalink pdb_costructure/pdb_chain_uniprot.tsv.gz pharos/pharos_targets.json
)

# --- Amendment 9 representations (agent: pathways/ppi/phenotype/citation)
# Raw files -> data/cache/reps/{reactome,biogrid,intact,monarch,cocitation}/; derived caches (derived_*.npz/json)
# are rebuilt on first use by scp_graphfs_mcp.reps.<name>. IntAct micluster is ~850 MB, gene2pubmed ~290 MB.
(
if [ "$(basename "$PWD")" = features ] || [ "$(basename "$PWD")" = reps ]; then cd "$(dirname "$PWD")"; fi
case "$(basename "$PWD")" in cache) ;; *) cd "$(dirname "$0")/../data/cache" ;; esac
mkdir -p reps/reactome reps/biogrid reps/intact reps/monarch reps/cocitation && cd reps
UA='scp-graphfs-mcp/0.1 (research)'
MONARCH=https://data.monarchinitiative.org/monarch-kg/latest
(cd reactome && curl -sS -L -R -A "$UA" -O https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt \
  && curl -sS -L -R -A "$UA" -O https://reactome.org/download/current/UniProt2Reactome.txt \
  && curl -sS -A "$UA" -o reactome_version.txt https://reactome.org/ContentService/data/database/version)
(cd biogrid && curl -sS -L -R -A "$UA" -O https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip)
(cd intact && curl -sS -L -R -A "$UA" -O https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact-micluster.zip)
(cd monarch && for f in gene_to_phenotypic_feature_association causal_gene_to_disease_association correlated_gene_to_disease_association; do
   curl -sS -L -R -A "$UA" -O "$MONARCH/tsv/all_associations/$f.all.tsv.gz"; done \
  && curl -sS -L -R -A "$UA" -O "$MONARCH/metadata.yaml")
(cd cocitation && curl -sS -L -R -A "$UA" -O https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz)
sha256sum reactome/UniProt2Reactome_All_Levels.txt reactome/UniProt2Reactome.txt biogrid/BIOGRID-ORGANISM-LATEST.tab3.zip \
  intact/intact-micluster.zip monarch/*.tsv.gz cocitation/gene2pubmed.gz
)

# --- Amendment 9 representations (agent: ccle/proteomehd/indra)
# Bulk co-regulation reference panels + mechanistic graph -> data/cache/reps/{ccle,proteomehd,pathwaycommons}/;
# derived caches (ccle_matrix.npz, proteomehd_matrix.npz, proteomehd_coreg.npz, pc14_edges.npz) are rebuilt on
# first use by scp_graphfs_mcp.reps.<name>. INDRA DB dumps were not publicly reachable (S3 403), so the
# mechanistic graph uses Pathway Commons v14 (see reps/pathwaycommons.py). coregulation_scores.zip is ~194 MB.
(
if [ "$(basename "$PWD")" = features ] || [ "$(basename "$PWD")" = reps ]; then cd "$(dirname "$PWD")"; fi
case "$(basename "$PWD")" in cache) ;; *) cd "$(dirname "$0")/../data/cache" ;; esac
mkdir -p reps/ccle reps/proteomehd reps/pathwaycommons && cd reps
UA='scp-graphfs-mcp/0.1 (research)'
PC=https://download.baderlab.org/PathwayCommons/PC2/v14
(cd ccle && curl -sS -L -R -A "$UA" -O https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz)
# ProteomeHD_v1_1.7z was obtained via `git clone --depth 1 https://github.com/Rappsilber-Laboratory/ProteomeHD`
# (Data/ProteomeHD_v1_1.7z); the raw URL below is the equivalent direct link.
(cd proteomehd && curl -sS -L -R -A "$UA" -O https://github.com/Rappsilber-Laboratory/ProteomeHD/raw/2ee6d87110b9d4932af0d106927eb289dfbce321/Data/ProteomeHD_v1_1.7z \
  && curl -sS -L -R -A "$UA" -O https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/11/PXD008888/coregulation_scores.zip \
  && if command -v 7z >/dev/null; then 7z x -y ProteomeHD_v1_1.7z >/dev/null; \
     else python3 -I -c 'import sys, py7zr; py7zr.SevenZipFile(sys.argv[1]).extractall(sys.argv[2])' ProteomeHD_v1_1.7z . \
       || echo "extract ProteomeHD_v1_1.7z manually (needs 7z or py7zr) -> ProteomeHD_v1_1.csv" >&2; fi)
(cd pathwaycommons && curl -sS -L -R -A "$UA" -O "$PC/pc-hgnc.txt.gz" && curl -sS -L -R -A "$UA" -O "$PC/datasources.txt")
sha256sum ccle/protein_quant_current_normalized.csv.gz proteomehd/ProteomeHD_v1_1.7z proteomehd/ProteomeHD_v1_1.csv \
  proteomehd/coregulation_scores.zip pathwaycommons/pc-hgnc.txt.gz pathwaycommons/datasources.txt
)

# --- Phase B bulk atlases
# Bulk proteomics reference panels -> data/cache/reps/{procan,gtex_proteome,cptac,atlantic}/; derived caches
# (procan_matrix.npz, gtex_proteome_matrix.npz, cptac_matrix.npz, atlantic_matrix.npz) are rebuilt on first use by
# scp_graphfs_mcp.reps.<name>. Checksums: data/PROVENANCE.md. ~420 MB in total; no file > 1 GB.
(
case "$(basename "$PWD")" in cache) ;; *) cd "$(dirname "$0")/../data/cache" ;; esac
mkdir -p reps/procan reps/gtex_proteome reps/cptac reps/atlantic && cd reps
UA='scp-graphfs-mcp/0.1 (research)'
# ProCan-DepMapSanger (Goncalves et al. 2022), figshare 10.6084/m9.figshare.19345397 v1 (figshare MD5s match).
(cd procan && curl -sS -L -R -A "$UA" -o ProCan-DepMapSanger_protein_matrix_8498_averaged.txt https://ndownloader.figshare.com/files/34411172 \
  && curl -sS -L -R -A "$UA" -o ProCan-DepMapSanger_mapping_file_averaged.txt https://ndownloader.figshare.com/files/34411133 \
  && curl -sS -L -R -A "$UA" -o README.md https://ndownloader.figshare.com/files/36155439)
# GTEx proteome (Jiang et al. 2020, Cell) Table S2 (Elsevier CDN; the PMC copy is behind a captcha).
(cd gtex_proteome && curl -sS -L -R -A "$UA" -O https://ars.els-cdn.com/content/image/1-s2.0-S0092867420310783-mmc3.xlsx)
# CPTAC BCM harmonised gene-level proteome tables, Zenodo record 8394329 (data store of the `cptac` package).
(cd cptac && for f in \
    brca:BRCA:Tumor ccrcc:CCRCC:Normal ccrcc:CCRCC:Tumor coad:COAD:Normal coad:COAD:Tumor gbm:GBM:Tumor \
    hnscc:HNSCC:Normal hnscc:HNSCC:Tumor lscc:LSCC:Normal lscc:LSCC:Tumor luad:LUAD:Normal luad:LUAD:Tumor \
    ov:OV:Normal ov:OV:Tumor pdac:PDAC:Normal pdac:PDAC:Tumor ucec:UCEC:Normal ucec:UCEC:Tumor; do
  IFS=: read -r c C k <<< "$f"
  n="bcm-$c-proteomics-${C}_proteomics_gene_abundance_log2_reference_intensity_normalized_$k.txt.gz"
  curl -sS -L -R -A "$UA" -o "$n" "https://zenodo.org/api/records/8394329/files/$n/content"
done && n=bcm-brca-mapping-gencode.v34.basic.annotation-mapping.txt.gz \
  && curl -sS -L -R -A "$UA" -o "$n" "https://zenodo.org/api/records/8394329/files/$n/content")
# ATLANTiC NCI60 + CRC65 proteomes (Frejno et al. 2020, Nat Commun), Supplementary Data 3.
(cd atlantic && curl -sS -L -R -A "$UA" -O "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-020-17336-9/MediaObjects/41467_2020_17336_MOESM6_ESM.xlsx")
sha256sum procan/ProCan-DepMapSanger_protein_matrix_8498_averaged.txt procan/ProCan-DepMapSanger_mapping_file_averaged.txt \
  procan/README.md gtex_proteome/1-s2.0-S0092867420310783-mmc3.xlsx cptac/bcm-*.txt.gz atlantic/41467_2020_17336_MOESM6_ESM.xlsx
)
