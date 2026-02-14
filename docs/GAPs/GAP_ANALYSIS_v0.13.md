# SATORI v0.13.0 スキル拡張ギャップ分析

> 作成日: 2026-02-12  
> 対象バージョン: v0.13.0 (86 スキル / 26 カテゴリ A–Z, 50 TU-linked)  
> 前回分析: GAP_ANALYSIS_v0.12.md (76 → 86 へ実装済み)

---

## 1. 調査概要

### 1.1 背景

v0.12.0 で 86 スキル / 50 ToolUniverse 連携を達成。Phase 4 (v0.12) で実装した 10 スキルはパスウェイ富化、文献検索、PPI ネットワーク、バリアント効果予測、がんゲノミクス、メタボロミクス DB、分子ドッキング、系統的レビュー、臨床レポーティング、タンパク質ドメインファミリーをカバーした。

### 1.2 調査ソース

| リポジトリ | 規模 | 最終確認日 |
|-----------|------|-----------|
| **mims-harvard/ToolUniverse** | `default_config.py`: **188 カテゴリ**, 1229+ ツール | 2026-02-12 |
| **K-Dense-AI/claude-scientific-skills** | **140 スキル** (databases: 28+, packages: 55+, integrations: 15+, analysis: 30+) | 2026-02-12 |
| **SATORI v0.12.0** | **86 スキル**, **50 TU-linked**, 26 カテゴリ (A–Z) | 現在 |

---

## 2. ToolUniverse カバレッジ全体像

### 2.1 SATORI が既にカバー済みの TU カテゴリ (72 categories)

#### v0.10 以前 (35 categories)
| # | Config Key | SATORI スキル |
|---|-----------|--------------|
| 1 | `fda_drug_adverse_event` | pharmacovigilance |
| 2 | `fda_drug_adverse_event_detail` | pharmacovigilance |
| 3 | `fda_drug_label` | regulatory-science |
| 4 | `dailymed` | regulatory-science |
| 5 | `oncokb` | precision-oncology |
| 6 | `civic` | precision-oncology |
| 7 | `opentarget` | drug-target-profiling |
| 8 | `monarch` | disease-research |
| 9 | `uniprot` | bioinformatics |
| 10 | `ChEMBL` | cheminformatics |
| 11 | `dgidb` | drug-target-profiling |
| 12 | `clinvar` | variant-interpretation |
| 13 | `gnomad` | variant-interpretation |
| 14 | `clingen` | variant-interpretation |
| 15 | `admetai` | admet-pharmacokinetics |
| 16 | `pubchem` | cheminformatics |
| 17 | `rcsb_pdb` | protein-structure-analysis |
| 18 | `rcsb_search` | protein-structure-analysis |
| 19 | `alphafold` | protein-structure-analysis |
| 20 | `geo` | gene-expression-transcriptomics |
| 21 | `gtex` | gene-expression-transcriptomics |
| 22 | `gtex_v2` | gene-expression-transcriptomics |
| 23 | `expression_atlas` | gene-expression-transcriptomics |
| 24 | `arrayexpress` | gene-expression-transcriptomics |
| 25 | `mgnify` | microbiome-metagenomics |
| 26 | `obis` | environmental-ecology |
| 27 | `gbif` | environmental-ecology |
| 28 | `iedb_tools` | immunoinformatics |
| 29 | `imgt` | immunoinformatics |
| 30 | `sabdab` | immunoinformatics |
| 31 | `cdc` | epidemiology-public-health |
| 32 | `who_gho` | epidemiology-public-health |
| 33 | `euhealth` | epidemiology-public-health |
| 34 | `cellxgene_census` | single-cell-genomics |
| 35 | `hca_tools` | single-cell-genomics |

#### v0.11 (11 categories)
| # | Config Key | SATORI スキル |
|---|-----------|--------------|
| 36 | `chipatlas` | epigenomics-chromatin |
| 37 | `jaspar` | epigenomics-chromatin |
| 38 | `screen` | epigenomics-chromatin |
| 39 | `encode` | epigenomics-chromatin |
| 40 | `pride` | proteomics-mass-spectrometry |
| 41 | `clinical_trials` | clinical-trials-analytics |
| 42 | `clinical_trials_tools` | clinical-trials-analytics |
| 43 | `fda_orange_book` | regulatory-science |
| 44 | `pharmgkb` | pharmacogenomics |
| 45 | `fda_pharmacogenomic_biomarkers` | pharmacogenomics |
| 46 | `gwas` | population-genetics |

#### v0.12 (25 categories)
| # | Config Key | SATORI スキル |
|---|-----------|--------------|
| 47 | `kegg` | pathway-enrichment |
| 48 | `reactome` | pathway-enrichment |
| 49 | `go` | pathway-enrichment |
| 50 | `wikipathways` | pathway-enrichment |
| 51 | `pathway_commons_tools` | pathway-enrichment |
| 52 | `pubmed` | literature-search |
| 53 | `EuropePMC` | literature-search |
| 54 | `semantic_scholar` | literature-search |
| 55 | `OpenAlex` | literature-search |
| 56 | `crossref` | literature-search |
| 57 | `intact` | protein-interaction-network |
| 58 | `ppi` | protein-interaction-network |
| 59 | `stitch` | protein-interaction-network |
| 60 | `HumanBase` | protein-interaction-network |
| 61 | `alphamissense` | variant-effect-prediction |
| 62 | `cadd` | variant-effect-prediction |
| 63 | `spliceai` | variant-effect-prediction |
| 64 | `cosmic` | cancer-genomics |
| 65 | `cbioportal` | cancer-genomics |
| 66 | `depmap` | cancer-genomics |
| 67 | `hmdb` | metabolomics-databases |
| 68 | `metacyc` | metabolomics-databases |
| 69 | `metabolomics_workbench` | metabolomics-databases |
| 70 | `interpro` | protein-domain-family |
| 71 | `interproscan` | protein-domain-family |

---

### 2.2 残存未カバー TU カテゴリ (科学的に重要なもの: 74 categories)

#### A. 遺伝性疾患・希少疾患 ★★★ (15+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| OMIM | `omim` | `OMIM_get_entry`, `OMIM_search_entries` | 2+ |
| Orphanet | `orphanet` | `Orphanet_get_disease`, `Orphanet_get_genes`, `Orphanet_get_classification`, `Orphanet_search_diseases`, `Orphanet_search_by_name` | 5 |
| DisGeNET | `disgenet` | `DisGeNET_get_gene_diseases`, `DisGeNET_search_associations` | 2+ |
| IMPC | `impc` | `IMPC_search_genes`, `IMPC_get_gene_summary`, `IMPC_get_phenotypes_by_gene`, `IMPC_get_gene_phenotype_hits` | 4 |

#### B. Human Protein Atlas ★★★ (15+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| HPA | `hpa` | `HPA_search_genes_by_query`, `HPA_get_subcellular_location`, `HPA_get_rna_expression_by_source`, `HPA_get_rna_expression_in_specific_tissues`, `HPA_get_biological_processes_by_gene`, `HPA_get_cancer_prognostics_by_gene`, `HPA_get_protein_interactions_by_gene`, `HPA_get_comprehensive_gene_details_by_ensembl_id`, `HPA_get_contextual_biological_process_analysis` + more | 15+ |

#### C. 薬理学・ターゲット ★★★ (18+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| BindingDB | `bindingdb` | `BindingDB_get_ligands_by_pdb`, `BindingDB_get_ligands_by_uniprot`, `BindingDB_get_ligands_by_uniprots`, `BindingDB_get_targets_by_compound` | 4 |
| GPCRdb | `gpcrdb` | `GPCRdb_get_ligands`, `GPCRdb_get_receptor_info` | 2+ |
| GtoPdb | `gtopdb` | `GtoPdb_get_targets`, `GtoPdb_list_diseases`, `GtoPdb_list_ligands`, `GtoPdb_search_interactions` | 4 |
| BRENDA | `brenda` | `BRENDA_get_inhibitors`, `BRENDA_get_km`, `BRENDA_get_kcat`, `BRENDA_get_enzyme_info` | 4 |
| Pharos/TCRD | `pharos` | `Pharos_get_target`, `Pharos_get_disease_targets`, `Pharos_get_tdl_summary`, `Pharos_search_targets` | 4 |

#### D. ゲノム・配列アノテーション ★★★ (12+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| Ensembl | `ensembl` | `ensembl_lookup_gene`, Ensembl REST tools | 2+ |
| dbSNP | `dbsnp` | `dbSNP_get_variant`, `dbSNP_search` | 2+ |
| BLAST | `blast` | `BLAST_nucleotide_search`, `BLAST_protein_search` | 2 |
| NCBI Nucleotide | `ncbi_nucleotide` | `NCBI_search_nucleotide`, `NCBI_get_sequence`, `NCBI_fetch_accessions` | 3 |
| GDC | `gdc` | GDC cancer genomics tools | 2+ |

#### E. クロスDB アノテーション ★★★ (13+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| BioThings | `biothings` | `MyGene_get_gene_annotation`, `MyGene_query_genes`, `MyGene_batch_query`, `MyChem_get_chemical_annotation`, `MyChem_query_chemicals`, `MyVariant_get_variant_annotation`, `MyVariant_query_variants` | 7 |
| ID Mapping | `idmap` | UniProt/PDB/ChEMBL/Gene cross-database ID mapping tools | 6+ |

#### F. RNA 生物学 ★★☆ (11 ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| Rfam | `rfam` | `Rfam_get_family`, `Rfam_search_sequence`, `Rfam_get_alignment`, `Rfam_get_covariance_model`, `Rfam_get_sequence_regions`, `Rfam_get_structure_mapping`, `Rfam_get_tree_data`, `Rfam_id_to_accession`, `Rfam_accession_to_id` | 9 |
| RNAcentral | `rnacentral` | `RNAcentral_search`, `RNAcentral_get_by_accession` | 2 |

#### G. 構造生物学（拡張） ★★☆ (21+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| EMDB | `emdb` | `EMDB_get_structure`, `EMDB_get_map_info`, `EMDB_get_sample_info`, `EMDB_get_publications`, `EMDB_get_imaging_info`, `EMDB_get_validation`, `EMDB_search_structures` | 7 |
| PDBe API | `pdbe_api` | `PDBe_get_summary`, PDBe tools | 2+ |
| Proteins API | `proteins_api` | `proteins_api_search`, `proteins_api_get_protein`, `proteins_api_get_variants`, `proteins_api_get_features`, `proteins_api_get_xrefs` + 3 more | 8+ |
| Complex Portal | `complex_portal` | `ComplexPortal_get_complex`, `ComplexPortal_search_complexes` | 2 |
| DeepGO | `deepgo` | `DeepGO_predict_function` | 1 |
| EVE | `eve` | EVE evolutionary variant predictions | 1+ |

#### H. 化合物ライブラリ ★★☆ (3+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| ZINC | `zinc` | `ZINCTool` (search 1.4B+ compounds) | 1+ |
| Enamine | `enamine` | Enamine compound catalog search | 1+ |
| eMolecules | `emolecules` | eMolecules vendor aggregator | 1+ |

#### I. 治療用抗体 ★★☆ (3 ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| Thera-SAbDab | `therasabdab` | `TheraSAbDab_get_all_therapeutics`, `TheraSAbDab_search_therapeutics`, `TheraSAbDab_search_by_target` | 3 |

#### J. ゲノムスケール代謝モデル ★★☆ (8+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| BiGG Models | `bigg_models` | `BiGG_get_model`, `BiGG_get_metabolite`, `BiGG_get_reaction`, `BiGG_search`, `BiGG_list_models`, `BiGG_get_model_reactions`, `BiGG_get_database_version` | 7 |
| BioModels | `biomodels_tools` | BioModels repository tools | 1+ |

#### K. テキストマイニング ★★☆ (2 ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| PubTator | `pubtator` | `PubTator3_LiteratureSearch`, `PubTator3_EntityAutocomplete` | 2 |

#### L. プレプリント・OAアーカイブ ★★☆ (15+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| bioRxiv | `biorxiv` | `BioRxiv_search_preprints`, `BioRxiv_get_preprint` | 2 |
| medRxiv | `medrxiv` | `MedRxiv_get_preprint` | 1 |
| arXiv | `arxiv` | `ArXiv_search_papers`, `ArXiv_get_pdf_snippets` | 2 |
| PMC | `pmc` | `PMC_search_papers` | 1 |
| DOAJ | `doaj` | `DOAJ_search_articles` | 1 |
| Unpaywall | `unpaywall` | `Unpaywall_check_oa_status` | 1 |
| HAL | `hal` | `HAL_search_archive` | 1 |
| CORE | `core` | `CORE_search_papers`, `CORE_get_fulltext_snippets` | 2 |
| Zenodo | `zenodo` | Zenodo dataset/publication search | 1+ |
| OpenAIRE | `openaire` | `OpenAIRE_search_publications` | 1 |
| OSF Preprints | `osf_preprints` | `OSF_search_preprints` | 1 |
| Fatcat | `fatcat` | `Fatcat_search_scholar` | 1 |
| DBLP | `dblp` | `DBLP_search_publications` | 1 |

#### M. 公衆衛生・患者情報 ★☆☆ (12+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| NHANES | `nhanes` | NHANES health data tools | 1+ |
| Health Disparities | `health_disparities` | `health_disparities_get_county_rankings_info`, `health_disparities_get_svi_info` | 2 |
| MedlinePlus | `medlineplus` | `MedlinePlus_search_topics_by_keyword`, `MedlinePlus_connect_lookup_by_code`, `MedlinePlus_get_genetics_condition_by_name`, `MedlinePlus_get_genetics_gene_by_name`, `MedlinePlus_get_genetics_index` | 5 |
| ODPHP | `odphp` | ODPHP health prevention tools | 1+ |
| RxNorm | `rxnorm` | `RxNorm_get_drug_names` | 1 |
| Guidelines | `guidelines` | Unified guideline tools | 2+ |

#### N. オントロジー・エンリッチメント ★☆☆ (6+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| EFO | `EFO` | `OSL_get_efo_id_by_disease_name` | 1 |
| OLS | `ols` | OLS ontology lookup tools | 2+ |
| Enrichr | `Enrichr` | Enrichr gene set enrichment | 1+ |
| UMLS | `umls` | `umls_search_concepts`, `umls_get_concept_details` | 2 |

#### O. EBI データベース群 ★☆☆ (5+ ツール)

| TU Category | Config Key | 主要ツール関数名 | ツール数 |
|---|---|---|---|
| EBI Search | `ebi_search` | EBI search tools | 1+ |
| ENA Browser | `ena_browser` | ENA nucleotide archive tools | 1+ |
| BioStudies | `biostudies` | BioStudies tools | 1+ |
| dbfetch | `dbfetch` | dbfetch data retrieval | 1+ |
| MetaboLights | `metabolights` | MetaboLights study tools | 2+ |

#### P. 追加ツール ★☆☆

| TU Category | Config Key | ツール数 | 備考 |
|---|---|---|---|
| Cellosaurus | `cellosaurus` | 3 | 細胞株データベース |
| Disease Target Score | `disease_target_score` | 5+ | ターゲットエビデンススコア |
| RegulomeDB | `regulomedb` | 1+ | 調節バリアント (epigenomics 拡張) |
| ReMap | `remap` | 1+ | 調節マップ (epigenomics 拡張) |
| 4DN | `fourdn` | tools | 3Dゲノム |
| Wikidata SPARQL | `wikidata_sparql` | 1 | ナレッジグラフ |
| DBpedia | `dbpedia` | 1 | ナレッジグラフ |
| Wikipedia | `wikipedia` | 1 | 一般知識 |
| SIMBAD | `simbad` | 2 | 天文学 |
| WoRMS | `worms` | 1+ | 海洋分類学 |
| Paleobiology | `paleobiology` | 1 | 古生物学 |
| MPD | `mpd` | 1 | マウス表現型 |
| NVIDIA NIM | `nvidia_nim` | tools | GPU推論API |

#### Q. インフラ・内部ツール (スキップ推奨: ~42 categories)

| TU Category | Config Key | 理由 |
|---|---|---|
| special_tools | `special_tools` | 内部制御 |
| tool_finder | `tool_finder` | TU内部ツール検索 |
| agents | `agents` | Agentic tools |
| smolagents | `smolagents` | Smolagents wrapper |
| tool_discovery_agents | `tool_discovery_agents` | ツール発見 |
| web_search_tools | `web_search_tools` | Web検索 |
| package_discovery_tools | `package_discovery_tools` | パッケージ発見 |
| pypi_package_inspector_tools | `pypi_package_inspector_tools` | PyPI検査 |
| drug_discovery_agents | `drug_discovery_agents` | 創薬AI |
| dataset | `dataset` | データセットアクセス |
| mcp_auto_loader_txagent | `mcp_auto_loader_txagent` | MCP loader |
| mcp_auto_loader_expert_feedback | `mcp_auto_loader_expert_feedback` | MCP loader |
| mcp_auto_loader_boltz | `mcp_auto_loader_boltz` | MCP loader |
| mcp_auto_loader_uspto_downloader | `mcp_auto_loader_uspto_downloader` | MCP loader |
| compose | `compose` | ツール合成 |
| python_executor | `python_executor` | Python実行 |
| tool_composition | `tool_composition` | ツール合成 |
| embedding | `embedding` | 埋め込み |
| output_summarization | `output_summarization` | 出力要約 |
| optimizer | `optimizer` | 最適化 |
| compact_mode | `compact_mode` | コンパクトモード |
| xml | `xml` | XML処理 |
| url | `url` | URL取得 |
| file_download | `file_download` | ファイルDL |
| markitdown | `markitdown` | ドキュメント変換 |
| software_bioinformatics | `software_bioinformatics` | パッケージ情報 |
| software_genomics | `software_genomics` | パッケージ情報 |
| software_single_cell | `software_single_cell` | パッケージ情報 |
| software_structural_biology | `software_structural_biology` | パッケージ情報 |
| software_cheminformatics | `software_cheminformatics` | パッケージ情報 |
| software_machine_learning | `software_machine_learning` | パッケージ情報 |
| software_visualization | `software_visualization` | パッケージ情報 |
| software_scientific_computing | `software_scientific_computing` | パッケージ情報 |
| software_physics_astronomy | `software_physics_astronomy` | パッケージ情報 |
| software_earth_sciences | `software_earth_sciences` | パッケージ情報 |
| software_image_processing | `software_image_processing` | パッケージ情報 |
| software_neuroscience | `software_neuroscience` | パッケージ情報 |
| visualization_protein_3d | `visualization_protein_3d` | 可視化ユーティリティ |
| visualization_molecule_2d | `visualization_molecule_2d` | 可視化ユーティリティ |
| visualization_molecule_3d | `visualization_molecule_3d` | 可視化ユーティリティ |
| adverse_event | `adverse_event` | faers と重複 |
| faers_analytics | `faers_analytics` | faers の拡張版 |

---

## 3. K-Dense-AI/claude-scientific-skills ギャップ分析

### 3.1 K-Dense 全 140 スキルリスト

K-Dense リポジトリの `scientific-skills/` ディレクトリに存在する全スキル:

adaptyv, aeon, alphafold-database, anndata, arboreto, astropy, benchling-integration, biopython, biorxiv-database, bioservices, brenda-database, cellxgene-census, chembl-database, cirq, citation-management, clinical-decision-support, clinical-reports, clinicaltrials-database, clinpgx-database, clinvar-database, cobrapy, cosmic-database, dask, datacommons-client, datamol, deepchem, deeptools, denario, diffdock, dnanexus-integration, document-skills, drugbank-database, ena-database, ensembl-database, esm, etetoolkit, exploratory-data-analysis, fda-database, flowio, fluidsim, fred-economic-data, gene-database, generate-image, geniml, geo-database, geopandas, get-available-resources, gget, gtars, gwas-database, histolab, hmdb-database, hypogenic, hypothesis-generation, imaging-data-commons, infographics, iso-13485-certification, kegg-database, labarchive-integration, lamindb, latchbio-integration, latex-posters, literature-review, market-research-reports, markitdown, matchms, matlab, matplotlib, medchem, metabolomics-workbench-database, modal, molfeat, networkx, neurokit2, neuropixels-analysis, offer-k-dense-web, omero-integration, openalex-database, opentargets-database, opentrons-integration, paper-2-web, pathml, pdb-database, peer-review, pennylane, perplexity-search, plotly, polars, pptx-posters, protocolsio-integration, pubchem-database, pubmed-database, pufferlib, pydeseq2, pydicom, pyhealth, pylabrobot, pymatgen, pymc, pymoo, pyopenms, pysam, pytdc, pytorch-lightning, qiskit, qutip, rdkit, reactome-database, research-grants, research-lookup, rowan, scanpy, scholar-evaluation, scientific-brainstorming, scientific-critical-thinking, scientific-schematics, scientific-slides, scientific-visualization, scientific-writing, scikit-bio, scikit-learn, scikit-survival, scvi-tools, seaborn, shap, simpy, stable-baselines3, statistical-analysis, statsmodels, string-database, sympy, torch_geometric, torchdrug, transformers, treatment-plans, umap-learn, uniprot-database, uspto-database, vaex, venue-templates, zarr-python, zinc-database

### 3.2 SATORI がカバー済みの K-Dense スキル (113/140)

以下は SATORI の 86 スキルが既にカバーしている K-Dense スキル:

| K-Dense スキル | カバーする SATORI スキル | カバー度 |
|---|---|---|
| alphafold-database | protein-structure-analysis | ◎ |
| anndata | single-cell-genomics | ◎ |
| arboreto | gene-expression-transcriptomics | ○ |
| benchling-integration | lab-data-management | ○ |
| biopython | bioinformatics, sequence-analysis | ◎ |
| biorxiv-database | literature-search | ◎ |
| bioservices | bioinformatics | ○ |
| brenda-database | metabolomics-databases | ○ |
| cellxgene-census | single-cell-genomics | ◎ |
| chembl-database | cheminformatics, drug-target-profiling | ◎ |
| cirq | quantum-computing | ◎ |
| citation-management | citation-checker | ◎ |
| clinical-decision-support | clinical-decision-support | ◎ |
| clinical-reports | clinical-reporting | ◎ |
| clinicaltrials-database | clinical-trials-analytics | ◎ |
| clinpgx-database | pharmacogenomics | ◎ |
| clinvar-database | variant-interpretation | ◎ |
| cosmic-database | cancer-genomics | ◎ |
| dask | data-preprocessing | ○ |
| datamol | cheminformatics | ○ |
| deepchem | cheminformatics, deep-learning | ○ |
| deeptools | epigenomics-chromatin | ○ |
| denario | multi-omics | △ |
| diffdock | molecular-docking | ◎ |
| dnanexus-integration | lab-data-management | ○ |
| drugbank-database | drug-target-profiling | ○ |
| ensembl-database | bioinformatics | ○ |
| esm | protein-design | ◎ |
| exploratory-data-analysis | eda-correlation | ◎ |
| fda-database | regulatory-science, pharmacovigilance | ◎ |
| gene-database | bioinformatics | ○ |
| geniml | epigenomics-chromatin | ○ |
| geo-database | gene-expression-transcriptomics | ◎ |
| geopandas | environmental-ecology | ○ |
| get-available-resources | pipeline-scaffold | ○ |
| gget | bioinformatics | ○ |
| gtars | epigenomics-chromatin | ○ |
| gwas-database | population-genetics | ◎ |
| histolab | medical-imaging | ○ |
| hmdb-database | metabolomics-databases | ◎ |
| hypogenic | hypothesis-pipeline | ○ |
| hypothesis-generation | hypothesis-pipeline | ◎ |
| iso-13485-certification | regulatory-science | ○ |
| kegg-database | pathway-enrichment | ◎ |
| labarchive-integration | lab-data-management | ○ |
| lamindb | lab-data-management | △ |
| latchbio-integration | lab-data-management | ○ |
| latex-posters | presentation-design | ○ |
| literature-review | systematic-review, deep-research | ◎ |
| matchms | proteomics-mass-spectrometry | ◎ |
| matplotlib | publication-figures | ◎ |
| medchem | cheminformatics | ○ |
| metabolomics-workbench-database | metabolomics-databases | ◎ |
| molfeat | cheminformatics | ○ |
| networkx | network-analysis | ◎ |
| neurokit2 | neuroscience-electrophysiology | ◎ |
| neuropixels-analysis | neuroscience-electrophysiology | ◎ |
| openalex-database | literature-search | ◎ |
| opentargets-database | drug-target-profiling | ◎ |
| opentrons-integration | lab-automation | ◎ |
| pathml | medical-imaging | ○ |
| pdb-database | protein-structure-analysis | ◎ |
| peer-review | peer-review-response | ◎ |
| pennylane | quantum-computing | ◎ |
| perplexity-search | deep-research | ○ |
| plotly | publication-figures | ○ |
| pptx-posters | presentation-design | ○ |
| protocolsio-integration | lab-automation | ○ |
| pubchem-database | cheminformatics | ◎ |
| pubmed-database | literature-search | ◎ |
| pydeseq2 | gene-expression-transcriptomics | ◎ |
| pydicom | medical-imaging | ○ |
| pylabrobot | lab-automation | ◎ |
| pymatgen | computational-materials | ◎ |
| pymc | bayesian-statistics | ◎ |
| pymoo | process-optimization | ○ |
| pyopenms | proteomics-mass-spectrometry | ◎ |
| pysam | bioinformatics | ○ |
| pytdc | drug-target-profiling | ○ |
| pytorch-lightning | deep-learning | ◎ |
| qiskit | quantum-computing | ◎ |
| qutip | quantum-computing | ◎ |
| rdkit | cheminformatics | ◎ |
| reactome-database | pathway-enrichment | ◎ |
| research-grants | grant-writing | ◎ |
| research-lookup | deep-research | ○ |
| scanpy | single-cell-genomics | ◎ |
| scholar-evaluation | critical-review | ○ |
| scientific-brainstorming | hypothesis-pipeline | ○ |
| scientific-critical-thinking | critical-review | ◎ |
| scientific-schematics | scientific-schematics | ◎ |
| scientific-slides | presentation-design | ◎ |
| scientific-visualization | publication-figures | ○ |
| scientific-writing | academic-writing | ◎ |
| scikit-bio | bioinformatics | ○ |
| scikit-learn | ml-classification, ml-regression | ◎ |
| scikit-survival | survival-clinical | ◎ |
| scvi-tools | single-cell-genomics | ◎ |
| seaborn | publication-figures | ○ |
| shap | explainable-ai | ◎ |
| statistical-analysis | statistical-testing | ◎ |
| statsmodels | statistical-testing | ○ |
| string-database | protein-interaction-network | ◎ |
| torch_geometric | graph-neural-networks | ◎ |
| torchdrug | drug-repurposing | ○ |
| transformers | deep-learning | ○ |
| treatment-plans | clinical-reporting | ◎ |
| umap-learn | pca-tsne | ○ |
| uniprot-database | bioinformatics | ◎ |
| uspto-database | regulatory-science | ○ |
| zinc-database | drug-repurposing | ○ |
| cobrapy | systems-biology | ○ |
| omero-integration | lab-data-management | ○ |
| venue-templates | academic-writing | ○ |

### 3.3 SATORI 未カバーの K-Dense スキル (27/140)

| # | K-Dense スキル | ドメイン | 重要度 | 理由 |
|---|---|---|---|---|
| 1 | **adaptyv** | タンパク質工学 | ★★☆ | クラウド実験プラットフォーム（タンパク質テスト） |
| 2 | **aeon** | 時系列ML | ★☆☆ | 既存 time-series と重複大 |
| 3 | **astropy** | 天文学 | ★☆☆ | 天文学ドメイン |
| 4 | **datacommons-client** | 公共データ | ★☆☆ | 統計データアクセス |
| 5 | **document-skills** | ドキュメント処理 | ★☆☆ | XLSX処理 |
| 6 | **ena-database** | ヌクレオチドDB | ★★☆ | European Nucleotide Archive |
| 7 | **etetoolkit** | 系統樹解析 | ★★☆ | 系統発生学 |
| 8 | **flowio** | フローサイトメトリー | ★★☆ | FCSファイル解析 |
| 9 | **fluidsim** | CFD | ★☆☆ | 流体シミュレーション |
| 10 | **fred-economic-data** | 経済データ | ★☆☆ | FRED API (NEW) |
| 11 | **generate-image** | AI画像生成 | ★☆☆ | FLUX/Gemini 画像生成 |
| 12 | **imaging-data-commons** | 医療画像DB | ★★☆ | NCI画像データ (NEW) |
| 13 | **infographics** | インフォグラフィック | ★☆☆ | 視覚的データ表現 (NEW) |
| 14 | **market-research-reports** | 市場調査 | ★☆☆ | ビジネスレポート |
| 15 | **markitdown** | 文書変換 | ★☆☆ | ファイル→Markdown |
| 16 | **matlab** | 数値計算 | ★★☆ | MATLAB/Octave |
| 17 | **modal** | クラウドコンピュート | ★☆☆ | サーバーレスGPU |
| 18 | **paper-2-web** | 学術出版 | ★☆☆ | 論文→Web/動画 |
| 19 | **polars** | 高速データ処理 | ★☆☆ | 既存 data-preprocessing と重複 |
| 20 | **pufferlib** | 強化学習 | ★☆☆ | 高性能RL |
| 21 | **pyhealth** | ヘルスケアAI | ★★★ | 臨床ML (EHR, ICU, 薬剤推薦) |
| 22 | **rowan** | 計算化学 | ★★☆ | クラウド量子化学 (NEW) |
| 23 | **simpy** | 離散イベントSIM | ★☆☆ | 離散事象シミュレーション |
| 24 | **stable-baselines3** | 強化学習 | ★☆☆ | RLアルゴリズム |
| 25 | **sympy** | 記号計算 | ★☆☆ | 記号数学 |
| 26 | **vaex** | 大規模データ処理 | ★☆☆ | 既存 data-preprocessing と重複 |
| 27 | **zarr-python** | 配列ストレージ | ★☆☆ | N次元配列フォーマット |

---

## 4. v0.13.0 推奨 TOP 10 新スキル

### 4.1 優先度スコアリング基準

| 基準 | 重み |
|------|------|
| 未使用 TU カテゴリのツール数 | 30% |
| K-Dense スキルカバレッジ | 20% |
| 既存パイプラインフロー完成度 | 25% |
| 科学的インパクト | 25% |

---

### 4.2 新スキル提案

---

### #1: `scientific-rare-disease-genetics` (Q-4)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | Q. 腫瘍学・疾患研究 |
| **優先度** | **CRITICAL** — 10/10 |
| **ToolUniverse カテゴリ** | `omim` (2+), `orphanet` (5), `disgenet` (2+), `impc` (4) |
| **ToolUniverse ツール関数** | `OMIM_get_entry`, `OMIM_search_entries`, `Orphanet_get_disease`, `Orphanet_get_genes`, `Orphanet_get_classification`, `Orphanet_search_diseases`, `Orphanet_search_by_name`, `DisGeNET_get_gene_diseases`, `DisGeNET_search_associations`, `IMPC_search_genes`, `IMPC_get_gene_summary`, `IMPC_get_phenotypes_by_gene`, `IMPC_get_gene_phenotype_hits` |
| **合計 TU ツール** | **13+** |
| **K-Dense カバー** | gene-database, ensembl-database (間接的) |
| **説明** | OMIM/Orphanet/DisGeNET/IMPCを統合した遺伝性疾患・希少疾患研究パイプライン。メンデル遺伝病の遺伝子-疾患関連マッピング、希少疾患分類・遺伝子候補同定、遺伝子-疾患アソシエーションスコアリング、マウスKO表現型からのヒト疾患類推。ACMG/AMPバリアント分類と直結する精密医療パイプラインの基盤。 |
| **正当化** | ① 4つの主要TUカテゴリが完全未活用 ② disease-research スキルとの相補性高 ③ 精密医療 (variant-interpretation → rare-disease-genetics) の自然な拡張 ④ 6,000+希少疾患は患者3.5億人に影響 |

---

### #2: `scientific-human-protein-atlas` (F-10)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | F. 生命科学・オミクス |
| **優先度** | **CRITICAL** — 9.5/10 |
| **ToolUniverse カテゴリ** | `hpa` (15+) |
| **ToolUniverse ツール関数** | `HPA_search_genes_by_query`, `HPA_get_subcellular_location`, `HPA_get_rna_expression_by_source`, `HPA_get_rna_expression_in_specific_tissues`, `HPA_get_biological_processes_by_gene`, `HPA_get_cancer_prognostics_by_gene`, `HPA_get_protein_interactions_by_gene`, `HPA_get_comprehensive_gene_details_by_ensembl_id`, `HPA_get_contextual_biological_process_analysis` + 6 more |
| **合計 TU ツール** | **15+** |
| **K-Dense カバー** | (直接対応なし — bioservices で間接的) |
| **説明** | Human Protein Atlas API を活用したタンパク質発現アトラス解析パイプライン。組織レベルタンパク質発現プロファイリング、細胞内局在解析、がん予後バイオマーカー探索、RNA発現ソース比較 (HPA/GTEx/FANTOM5)、遺伝子機能コンテキスト解析。ドラッグターゲット検証における組織特異性評価に必須。 |
| **正当化** | ① 単一カテゴリで15+ツール — SATORI最大の未活用TU リソース ② 組織/細胞レベル発現データは創薬・精密医療の基盤 ③ gene-expression-transcriptomics, drug-target-profiling との統合でパイプライン完成 |

---

### #3: `scientific-pharmacology-targets` (J-5)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | J. 創薬・ファーマコロジー |
| **優先度** | **HIGH** — 9.0/10 |
| **ToolUniverse カテゴリ** | `bindingdb` (4), `gpcrdb` (2+), `gtopdb` (4), `brenda` (4), `pharos` (4) |
| **ToolUniverse ツール関数** | `BindingDB_get_ligands_by_pdb`, `BindingDB_get_ligands_by_uniprot`, `BindingDB_get_ligands_by_uniprots`, `BindingDB_get_targets_by_compound`, `GPCRdb_get_ligands`, `GtoPdb_get_targets`, `GtoPdb_list_diseases`, `GtoPdb_list_ligands`, `GtoPdb_search_interactions`, `BRENDA_get_inhibitors`, `BRENDA_get_km`, `BRENDA_get_kcat`, `BRENDA_get_enzyme_info`, `Pharos_get_target`, `Pharos_get_disease_targets`, `Pharos_get_tdl_summary`, `Pharos_search_targets` |
| **合計 TU ツール** | **18+** |
| **K-Dense カバー** | brenda-database, opentargets-database (部分的) |
| **説明** | 包括的薬理学ターゲットプロファイリングパイプライン。BindingDB (2.8M+ アフィニティ測定)、GPCRdb (GPCR特化)、GtoPdb (IUPHAR公式ガイド)、BRENDA (酵素キネティクスKm/kcat/Vmax)、Pharos (IDGアンダースタディド・ターゲット) を統合。ターゲットドラッガビリティ評価、結合親和性プロファイル、酵素阻害剤スクリーニング。 |
| **正当化** | ① 5つのTUカテゴリ統合で18+ツール ② 創薬パイプライン (drug-target-profiling → pharmacology-targets → molecular-docking) を完成 ③ BRENDA は酵素研究の世界最大DB、BINDINGDBは結合アフィニティの権威 |

---

### #4: `scientific-genome-sequence-tools` (F-11)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | F. 生命科学・オミクス |
| **優先度** | **HIGH** — 8.5/10 |
| **ToolUniverse カテゴリ** | `ensembl` (2+), `dbsnp` (2+), `blast` (2), `ncbi_nucleotide` (3), `gdc` (2+) |
| **ToolUniverse ツール関数** | `ensembl_lookup_gene`, `dbSNP_get_variant`, `dbSNP_search`, `BLAST_nucleotide_search`, `BLAST_protein_search`, `NCBI_search_nucleotide`, `NCBI_get_sequence`, `NCBI_fetch_accessions`, GDC tools |
| **合計 TU ツール** | **12+** |
| **K-Dense カバー** | ensembl-database, ena-database, gene-database, gget, pysam, scikit-bio |
| **説明** | Ensembl/dbSNP/BLAST/NCBI/GDC を統合したゲノム配列・アノテーションパイプライン。遺伝子座情報取得、バリアントアノテーション、配列アライメント (BLAST)、ヌクレオチド配列検索・取得、NCI Genomic Data Commons (GDC) がんゲノムデータアクセス。bioinformatics スキルの API レベル拡張。 |
| **正当化** | ① Ensembl は脊椎動物ゲノムのリファレンス ② BLAST は配列比較の事実上の標準 ③ bioinformatics, variant-interpretation, cancer-genomics の下流パイプラインに配列レベルデータを供給 |

---

### #5: `scientific-biothings-idmapping` (F-12)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | F. 生命科学・オミクス |
| **優先度** | **HIGH** — 8.5/10 |
| **ToolUniverse カテゴリ** | `biothings` (7), `idmap` (6+) |
| **ToolUniverse ツール関数** | `MyGene_get_gene_annotation`, `MyGene_query_genes`, `MyGene_batch_query`, `MyChem_get_chemical_annotation`, `MyChem_query_chemicals`, `MyVariant_get_variant_annotation`, `MyVariant_query_variants`, UniProt/PDB/ChEMBL ID mapping tools |
| **合計 TU ツール** | **13+** |
| **K-Dense カバー** | bioservices (間接的) |
| **説明** | BioThings (MyGene/MyChem/MyVariant) + IDMap を統合したクロスデータベースアノテーション・ID変換パイプライン。遺伝子 → タンパク質 → 化合物 → バリアント間の ID マッピング、バッチアノテーション、マルチ DB 統合クエリ。全てのオミクス解析パイプラインのハブとしてデータ統合を実現。 |
| **正当化** | ① 13+ツールで最もユーティリティ性が高い ② 他の全スキルの入出力を橋渡し ③ BioThings は cross-database ID mapping の事実上の標準 |

---

### #6: `scientific-noncoding-rna` (F-13)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | F. 生命科学・オミクス |
| **優先度** | **HIGH** — 8.0/10 |
| **ToolUniverse カテゴリ** | `rfam` (9), `rnacentral` (2) |
| **ToolUniverse ツール関数** | `Rfam_get_family`, `Rfam_search_sequence`, `Rfam_get_alignment`, `Rfam_get_covariance_model`, `Rfam_get_sequence_regions`, `Rfam_get_structure_mapping`, `Rfam_get_tree_data`, `Rfam_id_to_accession`, `Rfam_accession_to_id`, `RNAcentral_search`, `RNAcentral_get_by_accession` |
| **合計 TU ツール** | **11** |
| **K-Dense カバー** | (直接対応なし) |
| **説明** | Rfam/RNAcentral を活用したノンコーディングRNA解析パイプライン。RNA ファミリー検索・分類、配列からの RNA ファミリー同定、二次構造マッピング、共分散モデル解析、RNA 配列領域解析。miRNA/lncRNA/rRNA/tRNA/ribozyme の包括的アノテーション。 |
| **正当化** | ① 11ツールで充実、特に Rfam は 9ツール ② ノンコーディング RNA は近年急成長中の分野 ③ K-Dense にも直接対応スキルなし — SATORI 独自の差別化要素 |

---

### #7: `scientific-structural-proteomics` (K-5)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | K. 構造生物学・タンパク質工学 |
| **優先度** | **HIGH** — 8.0/10 |
| **ToolUniverse カテゴリ** | `emdb` (7), `pdbe_api` (2+), `proteins_api` (8+), `complex_portal` (2), `deepgo` (1), `eve` (1+) |
| **ToolUniverse ツール関数** | `EMDB_get_structure`, `EMDB_get_map_info`, `EMDB_get_sample_info`, `EMDB_get_publications`, `EMDB_get_imaging_info`, `EMDB_get_validation`, `EMDB_search_structures`, PDBe tools, `proteins_api_search`, `proteins_api_get_protein`, `proteins_api_get_variants`, `proteins_api_get_features`, `proteins_api_get_xrefs`, `ComplexPortal_get_complex`, `ComplexPortal_search_complexes`, `DeepGO_predict_function`, EVE tools |
| **合計 TU ツール** | **21+** |
| **K-Dense カバー** | pdb-database (部分的) |
| **説明** | EMDB/PDBe/Proteins API/Complex Portal/DeepGO/EVE を統合した拡張構造プロテオミクスパイプライン。クライオ電顕構造解析 (EMDB)、PDB ヨーロッパ API 経由の高度構造検索、タンパク質機能予測 (DeepGO)、進化的バリアント効果 (EVE)、タンパク質複合体構造解析。既存 protein-structure-analysis を超える構造データ統合。 |
| **正当化** | ① 21+ツールで最大の未活用 TU リソース群 ② EMDB (クライオ電顕) は急成長分野 ③ protein-structure-analysis, protein-domain-family と合わせて構造生物学トライアドを完成 |

---

### #8: `scientific-compound-screening` (J-6)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | J. 創薬・ファーマコロジー |
| **優先度** | **MEDIUM-HIGH** — 7.5/10 |
| **ToolUniverse カテゴリ** | `zinc` (1+), `enamine` (1+), `emolecules` (1+), `therasabdab` (3) |
| **ToolUniverse ツール関数** | `ZINCTool`, Enamine search, eMolecules search, `TheraSAbDab_get_all_therapeutics`, `TheraSAbDab_search_therapeutics`, `TheraSAbDab_search_by_target` |
| **合計 TU ツール** | **6+** |
| **K-Dense カバー** | zinc-database, drugbank-database |
| **説明** | ZINC (14億化合物)/Enamine/eMolecules/TheraSAbDab を統合したバーチャルスクリーニング・化合物ライブラリパイプライン。市販化合物検索、Make-On-Demand 化合物探索、治療用抗体データベース統合、ドッキング用リガンドセット準備。molecular-docking の上流パイプライン。 |
| **正当化** | ① molecular-docking の上流 (リガンド準備) が欠落 ② ZINC は世界最大のバーチャルスクリーニングライブラリ ③ 治療用抗体 (TheraSAbDab) はバイオ医薬品研究に必須 |

---

### #9: `scientific-healthcare-ai` (H-6)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | H. 臨床・疫学・メタ科学 |
| **優先度** | **MEDIUM-HIGH** — 7.5/10 |
| **ToolUniverse カテゴリ** | — (K-Dense 独自) |
| **ToolUniverse ツール関数** | — |
| **合計 TU ツール** | 0 (K-Dense専用) |
| **K-Dense カバー** | **pyhealth** (★★★), flowio |
| **説明** | PyHealth/FlowIO を中心としたヘルスケアAI・臨床機械学習パイプライン。EHR データ (MIMIC-III/IV, eICU) からの予測モデル構築、ICU死亡率予測、再入院リスク評価、安全な薬剤推薦 (DDI制約付き)、フローサイトメトリーデータ解析 (FCSファイル処理)、医療コード変換 (ICD-9/10, NDC, RxNorm)。33+モデルアーキテクチャ対応。 |
| **正当化** | ① pyhealth は K-Dense で★★★級の未カバースキル ② 臨床ML は survival-clinical, clinical-decision-support との自然な拡張 ③ フローサイトメトリーは免疫学・腫瘍学の基盤技術 |

---

### #10: `scientific-metabolic-modeling` (W-3)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | W. システム生物学 |
| **優先度** | **MEDIUM** — 7.0/10 |
| **ToolUniverse カテゴリ** | `bigg_models` (7), `biomodels_tools` (1+) |
| **ToolUniverse ツール関数** | `BiGG_get_model`, `BiGG_get_metabolite`, `BiGG_get_reaction`, `BiGG_search`, `BiGG_list_models`, `BiGG_get_model_reactions`, `BiGG_get_database_version`, BioModels tools |
| **合計 TU ツール** | **8+** |
| **K-Dense カバー** | cobrapy |
| **説明** | BiGG Models/BioModels + COBRApy を統合したゲノムスケール代謝モデリングパイプライン。FBA (Flux Balance Analysis)/FVA、遺伝子ノックアウトシミュレーション、代謝フラックス予測、バイオプロダクション最適化、代謝モデル検索・比較。systems-biology スキルの代謝特化拡張。 |
| **正当化** | ① BiGG は 7ツールで充実 ② COBRApy はゲノムスケール代謝モデルの標準ツール ③ systems-biology, metabolomics-databases との統合でシステム代謝学パイプラインを完成 |

---

## 5. 実装優先度マトリクス

| 順位 | スキル名 | TU ツール数 | TU カテゴリ数 | KD カバー | 総合スコア |
|------|---------|-----------|------------|----------|-----------|
| **1** | rare-disease-genetics | 13+ | 4 | ○ | **10.0/10** |
| **2** | human-protein-atlas | 15+ | 1 | △ | **9.5/10** |
| **3** | pharmacology-targets | 18+ | 5 | ○ | **9.0/10** |
| **4** | genome-sequence-tools | 12+ | 5 | ◎ | **8.5/10** |
| **5** | biothings-idmapping | 13+ | 2 | △ | **8.5/10** |
| **6** | noncoding-rna | 11 | 2 | ✗ | **8.0/10** |
| **7** | structural-proteomics | 21+ | 6 | ○ | **8.0/10** |
| **8** | compound-screening | 6+ | 4 | ○ | **7.5/10** |
| **9** | healthcare-ai | 0 | 0 | ◎ | **7.5/10** |
| **10** | metabolic-modeling | 8+ | 2 | ◎ | **7.0/10** |

---

## 6. カバレッジ予測

| 指標 | v0.12.0 | v0.13.0 後 |
|------|---------|-----------|
| スキル数 | 86 | **96** |
| ToolUniverse 連携カテゴリ数 | 50 (71 config keys) | **71 (97 config keys)** |
| 新規 TU ツール参照数 | — | **130+** |
| ToolUniverse 科学カテゴリカバー率 | ~49% (71/146) | **~66% (97/146)** |
| K-Dense カバー率 | ~81% (113/140) | **~84% (118/140)** |

---

## 7. v0.14.0+ ロードマップ候補 (Phase 6 以降)

v0.13.0 後も残る主要ギャップ:

### Tier 1 (v0.14.0 候補)
| # | 候補スキル | TU カテゴリ | ツール数 |
|---|-----------|-----------|---------|
| 1 | preprint-archive-search | biorxiv, medrxiv, arxiv, pmc, doaj, unpaywall + 6 more | 15+ |
| 2 | public-health-surveillance | nhanes, health_disparities, medlineplus, odphp, guidelines, rxnorm | 12+ |
| 3 | medical-ontology-mapping | umls, EFO, ols, Enrichr | 6+ |
| 4 | cell-line-database | cellosaurus | 3 |
| 5 | text-mining-biomedical | pubtator | 2 |

### Tier 2 (v0.15.0 候補)
| # | 候補スキル | TU/KD 出典 |
|---|-----------|-----------|
| 6 | regulatory-genomics-extended | regulomedb, remap, fourdn |
| 7 | computational-chemistry-cloud | K-Dense: rowan |
| 8 | reinforcement-learning | K-Dense: stable-baselines3, pufferlib |
| 9 | symbolic-mathematics | K-Dense: sympy |
| 10 | phylogenetics | K-Dense: etetoolkit |

---

## 8. サマリ統計

| 区分 | カテゴリ/スキル数 |
|------|-----------------|
| TU 全カテゴリ | ~188 |
| TU 科学的カテゴリ | ~146 |
| TU インフラ/内部 (スキップ) | ~42 |
| SATORI 既カバー TU (v0.12) | 71 config keys |
| 残存未カバー TU (科学的) | 75 config keys |
| v0.13.0 で新規カバー予定 | 26 config keys |
| K-Dense 全スキル | 140 |
| SATORI 既カバー K-Dense | 113 |
| 残存未カバー K-Dense | 27 |
| v0.13.0 で新規カバー予定 | 5 |

---

*このドキュメントは `mims-harvard/ToolUniverse` の `default_config.py` (188 categories, 1229+ tools) および `K-Dense-AI/claude-scientific-skills` (140 skills) の 2026-02-12 時点のスナップショット分析に基づいています。*
