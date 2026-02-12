# ToolUniverse 未使用カテゴリ分析 — SATORI v0.11.1

> **生成日**: 2025-07  
> **ToolUniverse**: mims-harvard/ToolUniverse (default_config.py, ~188 categories, 1229+ tools)  
> **SATORI v0.11.1**: 42 existing skills mapped  
> **未使用カテゴリ数**: ~100+ (うち科学的に重要なもの ~65)

---

## 凡例

- **重要度**: ★★★ = 必須級 / ★★☆ = 高優先 / ★☆☆ = 有用だが低優先 / ○ = インフラ/内部ツール
- **ツール数**: TU docs/tools_config_index.rst に基づく概算
- **提案スキル名**: SATORI skill として推奨されるマッピング名

---

## A. バリアント効果予測・スプライシング ★★★

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `alphamissense` | `AlphaMissense_get_variant_score`, `AlphaMissense_get_residue_scores`, `AlphaMissense_get_protein_scores` | 3 | `variant-effect-prediction` |
| `cadd` | `CADD_get_variant_score`, `CADD_get_position_scores`, `CADD_get_range_scores` | 3 | `variant-effect-prediction` |
| `spliceai` | `SpliceAI_get_max_delta`, `SpliceAI_predict_pangolin`, `SpliceAI_predict_splice` | 3 | `splice-prediction` |
| `eve` | EVE evolutionary variant effect predictions | 1+ | `variant-effect-prediction` |

**科学的根拠**: ACMG/AMPバリアント分類の計算エビデンス(PP3/BP4)に必須。AlphaMissense, CADD, SpliceAIは臨床遺伝学の標準ツール。

---

## B. パスウェイ・機能アノテーション ★★★

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `kegg` | `kegg_find_genes`, `kegg_get_gene_info`, `kegg_search_pathway`, `kegg_get_pathway_info`, `kegg_list_organisms` | 5 | `pathway-analysis` |
| `reactome` | `Reactome_get_pathway_reactions`, `Reactome_map_uniprot_to_pathways`, `Reactome_get_pathway`, `Reactome_list_top_pathways`, + 15 more | 18+ | `pathway-analysis` |
| `go` | `GO_get_annotations_for_gene`, `GO_get_genes_for_term`, `GO_get_term_by_id`, `GO_get_term_details`, `GO_search_terms` | 5 | `gene-ontology` |
| `wikipathways` | `WikiPathways_get_pathway`, `WikiPathways_search` | 2 | `pathway-analysis` |
| `pathway_commons_tools` | `pc_search_pathways`, pathway interactions | 2+ | `pathway-analysis` |

**科学的根拠**: パスウェイ解析は drug target validation、メカニズム研究、precision medicine の基盤。KEGG + Reactome + GO は生命科学研究の3大リソース。

---

## C. 文献検索・ナレッジベース ★★★

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `pubmed` | `PubMed_search_articles`, `PubMed_get_article`, `PubMed_get_cited_by`, `PubMed_get_related`, `PubMed_get_links` | 5 | `literature-search` |
| `semantic_scholar` | `SemanticScholar_search_papers`, `SemanticScholar_get_pdf_snippets` | 2 | `literature-search` |
| `EuropePMC` | `EuropePMC_search_articles` | 1 | `literature-search` |
| `OpenAlex` | `openalex_literature_search` | 1 | `literature-search` |
| `pubtator` | `PubTator3_LiteratureSearch`, `PubTator3_EntityAutocomplete` | 2 | `text-mining` |
| `crossref` | `Crossref_search_works`, `Crossref_get_work`, `Crossref_get_journal`, `Crossref_get_funder`, `Crossref_list_funders`, `Crossref_list_types` | 6 | `literature-search` |
| `biorxiv` | `BioRxiv_search_preprints`, `BioRxiv_get_preprint` | 2 | `preprint-search` |
| `medrxiv` | `MedRxiv_get_preprint` | 1 | `preprint-search` |
| `arxiv` | `ArXiv_search_papers`, `ArXiv_get_pdf_snippets` | 2 | `preprint-search` |
| `pmc` | `PMC_search_papers` | 1 | `literature-search` |
| `literature_search` | `IntentAnalyzerAgent`, `KeywordExtractorAgent`, `LiteratureSearchTool`, `LiteratureSynthesisAgent`, `MultiAgentLiteratureSearch`, `OverallSummaryAgent` | 6 | `literature-deep-research` |

**科学的根拠**: PubMed単独で3600万論文。Semantic Scholarの引用ネットワーク解析は研究インパクト評価に必須。文献エビデンスはすべてのバイオメディカル判断の基盤。

---

## D. タンパク質ドメイン・構造解析（拡張） ★★★

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `interpro` | `InterPro_get_domain_details`, `InterPro_get_protein_domains`, `InterPro_search_domains` | 3 | `protein-domain-analysis` |
| `interproscan` | `InterProScan_scan_sequence`, `InterProScan_get_job_status`, `InterProScan_get_job_results` | 3 | `protein-domain-analysis` |
| `pdbe_api` | PDBe API tools | 2+ | `protein-structure-analysis` |
| `rcsb_search` | `PDB_search_similar_structures` | 1+ | `protein-structure-analysis` |
| `emdb` | `EMDB_get_structure`, `EMDB_get_map_info`, `EMDB_get_sample_info`, `EMDB_get_publications`, `EMDB_get_imaging_info`, `EMDB_get_validation`, `EMDB_search_structures` | 7 | `cryo-em-analysis` |
| `deepgo` | `DeepGO_predict_function` | 1 | `protein-function-prediction` |
| `complex_portal` | `ComplexPortal_get_complex`, `ComplexPortal_search_complexes` | 2 | `protein-complex-analysis` |

**科学的根拠**: InterProは 100M+ タンパク質のドメイン分類データベース。InterProScan でカスタム配列解析可能。EMDBはクライオ電顕構造のゴールドスタンダード。

---

## E. 分子相互作用ネットワーク ★★★

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `ppi` | `STRING_get_protein_interactions`, `BioGRID_get_interactions` | 2 | `protein-interaction-network` |
| `stitch` | `STITCH_get_chemical_protein_interactions`, `STITCH_get_interaction_partners`, `STITCH_resolve_identifier` | 3 | `chemical-protein-interactions` |
| `intact` | `intact_get_interactions`, `intact_search_interactions`, `intact_get_interaction_network`, `intact_get_complex_details`, `intact_get_interactor`, + 3 more | 8 | `molecular-interactions` |
| `HumanBase` | `humanbase_ppi_analysis` (tissue-specific NetworkX graph) | 1 | `tissue-specific-networks` |

**科学的根拠**: STRING は 67M タンパク質。PPI ネットワーク解析は drug target discovery、genetic pathway 分析の中核。STITCHは化合物-タンパク質相互作用で drug repurposing に必須。

---

## F. 遺伝性疾患・希少疾患 ★★★

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `omim` | OMIM disease-gene relationships | 2+ | `mendelian-disease` |
| `orphanet` | `Orphanet_get_disease`, `Orphanet_get_genes`, `Orphanet_get_classification`, `Orphanet_search_diseases`, `Orphanet_search_by_name` | 5 | `rare-disease-research` |
| `impc` | `IMPC_search_genes`, `IMPC_get_gene_summary`, `IMPC_get_phenotypes_by_gene`, `IMPC_get_gene_phenotype_hits` | 4 | `mouse-phenotype-analysis` |

**科学的根拠**: OMIM は遺伝性疾患の権威データベース。Orphanet は 6000+ 希少疾患を網羅。IMPCは 9000+ ノックアウトマウスの表現型データ。

---

## G. ゲノム・配列解析 ★★☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `ensembl` | `ensembl_lookup_gene` + Ensembl REST tools | 2+ | `genome-annotation` |
| `blast` | `BLAST_nucleotide_search`, `BLAST_protein_search` | 2 | `sequence-alignment` |
| `dbsnp` | dbSNP variant lookup tools | 2+ | `variant-annotation` |
| `encode` | `ENCODE_get_experiment`, `ENCODE_get_biosample`, `ENCODE_get_file`, `ENCODE_list_files`, `ENCODE_search_biosamples`, `ENCODE_search_experiments` | 6+ | `functional-genomics` |
| `ncbi_nucleotide` | `NCBI_search_nucleotide`, `NCBI_get_sequence`, `NCBI_fetch_accessions` | 3 | `sequence-retrieval` |
| `rnacentral` | `RNAcentral_search`, `RNAcentral_get_by_accession` | 2 | `non-coding-rna` |
| `rfam` | `Rfam_get_family`, `Rfam_search_sequence`, `Rfam_get_alignment`, `Rfam_get_covariance_model`, `Rfam_get_sequence_regions`, `Rfam_get_structure_mapping`, `Rfam_get_tree_data`, `Rfam_id_to_accession`, `Rfam_accession_to_id` | 9 | `rna-families` |

**科学的根拠**: Ensembl は vertebrate genome のリファレンスアノテーション。BLASTは配列比較の事実上の標準。ENCODEは機能ゲノミクスの大規模コンソーシアム。

---

## H. 創薬・化合物ライブラリ ★★☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `bindingdb` | `BindingDB_get_ligands_by_pdb`, `BindingDB_get_ligands_by_uniprot`, `BindingDB_get_ligands_by_uniprots`, `BindingDB_get_targets_by_compound` | 4 | `binding-affinity-data` |
| `zinc` | `ZINCTool` (search 1.4B+ purchasable compounds) | 1+ | `compound-library` |
| `enamine` | Enamine compound catalog search | 1+ | `compound-library` |
| `emolecules` | eMolecules compound search | 1+ | `compound-library` |
| `pharos` | `Pharos_get_target`, `Pharos_get_disease_targets`, `Pharos_get_tdl_summary`, `Pharos_search_targets` | 4 | `target-druggability` |
| `gpcrdb` | `GPCRdb_get_ligands`, GPCRdb receptor tools | 2+ | `gpcr-pharmacology` |
| `brenda` | `BRENDA_get_inhibitors`, `get_km`, `get_kcat`, `get_enzyme_info` | 4 | `enzyme-kinetics` |
| `therasabdab` | `TheraSAbDab_get_all_therapeutics`, `TheraSAbDab_search_therapeutics`, `TheraSAbDab_search_by_target` | 3 | `therapeutic-antibodies` |
| `sabdab` | `SAbDab_get_structure`, `SAbDab_get_summary`, `SAbDab_search_structures` | 3 | `antibody-structures` |

**科学的根拠**: ZINCは14億化合物のvirtual screening ライブラリ。BindingDB は 2.8M+ のaffinity測定値。BRENDAは酵素反応速度パラメータの世界最大データベース。

---

## I. 免疫・抗体 ★★☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `imgt` | `IMGT_get_gene_info`, `IMGT_get_sequence`, `IMGT_search_genes` | 3 | `immunogenetics` |
| `therasabdab` | (上記参照) | 3 | `therapeutic-antibodies` |
| `sabdab` | (上記参照) | 3 | `antibody-structures` |

**科学的根拠**: IMGTは免疫グロブリン/TCR遺伝子の国際標準。バイオ医薬品(抗体医薬)開発に必須。

---

## J. メタボロミクス・代謝パスウェイ ★★☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `hmdb` | `HMDB_get_metabolite`, `HMDB_get_diseases`, `HMDB_search` | 3 | `metabolomics` |
| `metacyc` | `MetaCyc_get_pathway`, `MetaCyc_get_compound`, `MetaCyc_get_reaction`, `MetaCyc_search_pathways` | 4 | `metabolic-pathways` |
| `metabolomics_workbench` | `MetabolomicsWorkbench_get_study`, `MetabolomicsWorkbench_search_compound_by_name`, `MetabolomicsWorkbench_get_refmet_info`, `MetabolomicsWorkbench_search_by_exact_mass`, `MetabolomicsWorkbench_search_by_mz`, `MetabolomicsWorkbench_get_compound_by_pubchem_cid` | 6 | `metabolomics` |
| `bigg_models` | `BiGG_get_model`, `BiGG_get_metabolite`, `BiGG_get_reaction`, `BiGG_search`, `BiGG_list_models`, `BiGG_get_model_reactions`, `BiGG_get_database_version` | 7 | `genome-scale-models` |
| `metabolights` | MetaboLights study/compound tools | 2+ | `metabolomics` |

**科学的根拠**: HMDBは 220,000+ ヒト代謝物。代謝フェノタイピング・biomarker discoveryに必須。BiGGはゲノムスケール代謝モデルのリファレンス。

---

## K. 公衆衛生・医薬品情報 ★★☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `dailymed` | `DailyMed_search_spls`, `DailyMed_get_spl_by_setid`, `DailyMed_parse_adverse_reactions`, `DailyMed_parse_dosing`, `DailyMed_parse_drug_interactions`, `DailyMed_parse_contraindications`, `DailyMed_parse_clinical_pharmacology` | 7 | `drug-labeling` |
| `rxnorm` | `RxNorm_get_drug_names` | 1 | `drug-terminology` |
| `medlineplus` | `MedlinePlus_search_topics_by_keyword`, `MedlinePlus_connect_lookup_by_code`, `MedlinePlus_get_genetics_condition_by_name`, `MedlinePlus_get_genetics_gene_by_name`, `MedlinePlus_get_genetics_index` | 5 | `patient-health-info` |
| `umls` | `umls_search_concepts`, `umls_get_concept_details` | 2 | `medical-ontology` |
| `cdc` | `cdc_data_get_dataset`, CDC REST tools | 2 | `public-health-surveillance` |
| `nhanes` | NHANES health/nutrition data | 1+ | `public-health-surveillance` |
| `health_disparities` | `health_disparities_get_county_rankings_info`, `health_disparities_get_svi_info` | 2 | `health-equity-analysis` |
| `who_gho` | `who_gho_get_data`, `who_gho_query_health_data` | 2 | `global-health-data` |
| `odphp` | ODPHP health prevention tools | 1+ | `health-prevention` |
| `euhealth` | EU health info tools (vaccination, surveillance, mortality, primary care workforce) | 5+ | `eu-health-data` |

**科学的根拠**: DailyMedはFDA承認ラベル全文。UMLSは統一医学言語システムで 200+ 語彙統合。CDC/NHANESは米国公衆衛生サーベイランスの基盤。

---

## L. 単一細胞・遺伝子発現（拡張） ★★☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `hpa` | `HPA_search_genes_by_query`, `HPA_get_subcellular_location`, `HPA_get_rna_expression_by_source`, `HPA_get_rna_expression_in_specific_tissues`, `HPA_get_biological_processes_by_gene`, `HPA_get_cancer_prognostics_by_gene`, `HPA_get_protein_interactions_by_gene`, `HPA_get_comprehensive_gene_details_by_ensembl_id`, `HPA_get_contextual_biological_process_analysis`, + more | 15+ | `human-protein-atlas` |
| `cellxgene_census` | `CELLxGENE_get_expression_data`, `CELLxGENE_get_cell_metadata`, `CELLxGENE_get_gene_metadata`, `CELLxGENE_get_embeddings`, `CELLxGENE_get_presence_matrix`, `CELLxGENE_get_census_versions`, `CELLxGENE_download_h5ad` | 7 | `single-cell-atlas` |
| `hca_tools` | `hca_search_projects`, `hca_get_file_manifest` | 2 | `human-cell-atlas` |
| `cellosaurus` | `CellosaurusSearchTool`, `CellosaurusGetCellLineInfoTool`, `CellosaurusQueryConverterTool` | 3 | `cell-line-database` |
| `arrayexpress` | ArrayExpress experiment search | 1+ | `microarray-data` |

**科学的根拠**: HPAは タンパク質発現の組織・細胞レベルアトラス。CELLxGENE Census は 50M+ 単一細胞データ。single-cell は最も成長著しい分野。

---

## M. がん（追加） ★★☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `cbioportal` | `cBioPortal_get_cancer_studies`, `cBioPortal_get_mutations` | 2 | `cancer-genomics-portal` |
| `depmap` | `DepMap_get_gene_dependencies`, `DepMap_search_genes`, `DepMap_get_cell_line`, `DepMap_get_cell_lines`, `DepMap_search_cell_lines` | 5 | `cancer-dependency-map` |
| `disease_target_score` | `reactome_disease_target_score`, `europepmc_disease_target_score`, + other scores | 5+ | `target-evidence-scoring` |

**科学的根拠**: cBioPortalは 100+ がんゲノムスタディ。DepMapは 1000+ がん細胞株の遺伝子依存性データで synthetic lethality 研究に必須。

---

## N. 医学用語・オントロジー ★☆☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `EFO` | `OSL_get_efo_id_by_disease_name` | 1 | `ontology-mapping` |
| `ols` | `TermSearchResponse`, OLS ontology lookup | 2+ | `ontology-lookup` |
| `Enrichr` | Enrichr gene set enrichment | 1+ | `gene-set-enrichment` |
| `idmap` | ID mapping tools (MyGene, MyChem, MyVariant) | 6+ | `id-mapping` |
| `biothings` | `MyGene_get_gene_annotation`, `MyGene_query_genes`, `MyGene_batch_query`, `MyChem_get_chemical_annotation`, `MyChem_query_chemicals`, `MyVariant_get_variant_annotation`, `MyVariant_query_variants` | 7 | `biothings-annotation` |

**科学的根拠**: BioThings (MyGene/MyChem/MyVariant) は cross-database ID マッピングの事実上の標準。Enrichr は遺伝子セットエンリッチメント解析の主要ツール。

---

## O. 知識グラフ・セマンティックウェブ ★☆☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `wikidata_sparql` | Wikidata SPARQL query | 1 | `knowledge-graph-query` |
| `wikipedia` | Wikipedia article retrieval | 1 | (utility) |
| `dbpedia` | `DBpedia_SPARQL_query` | 1 | `knowledge-graph-query` |

---

## P. 特許・知財 ★☆☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `uspto` | `get_abstract_from_patent_app_number`, `get_claims_from_patent_app_number` | 2+ | `patent-search` |

**科学的根拠**: 創薬における IP landscape 分析、FTO (Freedom to Operate) 評価。

---

## Q. 文献アーカイブ（補完） ★☆☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `doaj` | `DOAJ_search_articles` | 1 | `literature-search` |
| `unpaywall` | `Unpaywall_check_oa_status` | 1 | `open-access-check` |
| `hal` | `HAL_search_archive` | 1 | `literature-search` |
| `core` | `CORE_search_papers`, `CORE_get_fulltext_snippets` | 2 | `literature-search` |
| `zenodo` | Zenodo dataset/publication search | 1+ | `research-data` |
| `openaire` | `OpenAIRE_search_publications` | 1 | `literature-search` |
| `osf_preprints` | `OSF_search_preprints` | 1 | `preprint-search` |
| `fatcat` | `Fatcat_search_scholar` | 1 | `literature-search` |
| `dblp` | `DBLP_search_publications` | 1 | `cs-literature` |
| `simbad` | `SIMBAD_query_object`, `SIMBAD_advanced_query` | 2 | `astronomy` |

---

## R. 環境・生態系（追加） ★☆☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `mgnify` | `MGnify_search_studies`, `MGnify_list_analyses` | 2 | `metagenomics` |
| `worms` | WoRMS marine species tools | 1+ | `marine-taxonomy` |
| `paleobiology` | `Paleobiology_get_fossils` | 1 | `paleobiology` |
| `mpd` | `MPD_get_phenotype_data` | 1 | `mouse-phenotype-data` |

---

## S. 薬理学・ガイドライン（追加） ★☆☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `gtopdb` | `GtoPdb_get_targets`, `GtoPdb_list_diseases`, `GtoPdb_list_ligands`, `GtoPdb_search_interactions` | 4 | `pharmacology-guide` |
| `clinical_trials_tools` | Additional CT extraction tools (`extract_clinical_trial_outcomes`, `extract_clinical_trial_adverse_events`) | 2+ | `clinical-trials-analytics` (既存拡張) |

---

## T. 可視化 ★☆☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `visualization_protein_3d` | `visualize_protein_structure_3d` | 1 | `molecular-visualization` |
| `visualization_molecule_2d` | `visualize_molecule_2d` | 1 | `molecular-visualization` |
| `visualization_molecule_3d` | `visualize_molecule_3d` | 1 | `molecular-visualization` |

---

## U. EBI データベース群 ★☆☆

| TU Category | 主要ツール | ツール数 | 提案スキル名 |
|---|---|---|---|
| `ebi_search` | EBI search tools | 1+ | `ebi-data-access` |
| `ena_browser` | ENA nucleotide archive browser | 1+ | `nucleotide-archive` |
| `biostudies` | BioStudies tools | 1+ | `biostudies` |
| `dbfetch` | dbfetch data retrieval | 1+ | `ebi-data-access` |
| `proteins_api` | `proteins_api_search`, `proteins_api_get_protein`, `proteins_api_get_variants`, `proteins_api_get_features`, `proteins_api_get_xrefs`, + 3 more | 8+ | `protein-api` |

---

## V. インフラ・ユーティリティ ○

以下は内部ツール/インフラで、SATORI skill としては不要:

| TU Category | 説明 |
|---|---|
| `special_tools` | CallAgent, Finish (内部制御) |
| `tool_finder` | Tool_Finder, Tool_Finder_Keyword, Tool_Finder_LLM (ツール検索) |
| `agents` / `smolagents` / `tool_discovery_agents` | Agentic tools (TU内部AI) |
| `drug_discovery_agents` | Drug discovery AI agents |
| `web_search_tools` | `web_search`, `web_api_documentation_search` |
| `package_discovery_tools` / `pypi_package_inspector_tools` | Package inspection |
| `compose` / `tool_composition` | Tool composition workflows |
| `python_executor` | Python code execution |
| `embedding` | Embedding tools |
| `output_summarization` | Output summarization |
| `optimizer` / `compact_mode` | Performance optimization |
| `dataset` | Dataset access |
| `xml` / `url` / `file_download` / `markitdown` | Data processing utilities |
| `mcp_auto_loader_*` | MCP auto-loader tools |
| `nvidia_nim` | NVIDIA NIM inference |
| `software_*` (8 categories) | Package info tools (bioinformatics, genomics, etc.) |
| `gtex` | Older GTEx (gtex_v2 already used) |

---

## 推奨: 新スキル優先実装リスト

### Phase 1 (★★★ — 即時追加推奨)

| 優先順 | 提案スキル名 | TU Categories | 推定ツール数 |
|---|---|---|---|
| 1 | `pathway-analysis` | kegg, reactome, go, wikipathways, pathway_commons_tools | 32+ |
| 2 | `literature-search` | pubmed, semantic_scholar, EuropePMC, OpenAlex, crossref, pmc, biorxiv, medrxiv, arxiv | 20+ |
| 3 | `variant-effect-prediction` | alphamissense, cadd, eve | 7+ |
| 4 | `splice-prediction` | spliceai | 3 |
| 5 | `protein-domain-analysis` | interpro, interproscan | 6 |
| 6 | `protein-interaction-network` | ppi, stitch, intact, HumanBase | 14+ |
| 7 | `rare-disease-research` | omim, orphanet | 7+ |

### Phase 2 (★★☆ — 高優先)

| 優先順 | 提案スキル名 | TU Categories | 推定ツール数 |
|---|---|---|---|
| 8 | `human-protein-atlas` | hpa | 15+ |
| 9 | `single-cell-atlas` | cellxgene_census, hca_tools | 9 |
| 10 | `metabolomics` | hmdb, metabolomics_workbench, metabolights | 11+ |
| 11 | `drug-labeling` | dailymed, rxnorm | 8 |
| 12 | `cancer-dependency-map` | depmap, cbioportal | 7 |
| 13 | `binding-affinity-data` | bindingdb | 4 |
| 14 | `compound-library` | zinc, enamine, emolecules | 3+ |
| 15 | `biothings-annotation` | biothings, idmap | 13+ |
| 16 | `genome-annotation` | ensembl, dbsnp, encode | 10+ |
| 17 | `therapeutic-antibodies` | therasabdab, sabdab, imgt | 9 |
| 18 | `mouse-phenotype-analysis` | impc, mpd | 5 |

### Phase 3 (★☆☆ — 有用)

| 優先順 | 提案スキル名 | TU Categories | 推定ツール数 |
|---|---|---|---|
| 19 | `public-health-surveillance` | cdc, nhanes, who_gho, health_disparities, odphp, euhealth | 12+ |
| 20 | `medical-ontology` | umls, EFO, ols | 5+ |
| 21 | `gene-set-enrichment` | Enrichr | 1+ |
| 22 | `metabolic-pathways` | metacyc, bigg_models | 11 |
| 23 | `non-coding-rna` | rnacentral, rfam | 11 |
| 24 | `sequence-alignment` | blast, ncbi_nucleotide | 5 |
| 25 | `target-druggability` | pharos, gpcrdb, gtopdb | 10+ |
| 26 | `enzyme-kinetics` | brenda | 4 |
| 27 | `cell-line-database` | cellosaurus | 3 |
| 28 | `patent-search` | uspto | 2+ |
| 29 | `molecular-visualization` | visualization_* (3 categories) | 3 |
| 30 | `text-mining` | pubtator | 2 |
| 31 | `metagenomics` | mgnify | 2 |
| 32 | `patient-health-info` | medlineplus | 5 |

---

## サマリ統計

| 区分 | カテゴリ数 | 推定ツール数 |
|---|---|---|
| SATORI 既使用 | ~40 | ~300+ |
| 科学的に重要な未使用 | ~65 | ~350+ |
| インフラ/内部 (スキップ推奨) | ~35 | ~100+ |
| **合計 TU カテゴリ** | **~140** | **1229** |

---

*このドキュメントは `mims-harvard/ToolUniverse` の `default_config.py` および `tools/__init__.py` (1229 tools) の分析に基づいています。*
