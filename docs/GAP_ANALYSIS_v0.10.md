# SATORI v0.10.0 スキル拡張ギャップ分析

> 作成日: 2025-07  
> 対象バージョン: v0.10.0 (56 スキル / 19 カテゴリ A–S)  
> 前回分析: GAP_ANALYSIS.md (47 → 56 へ Phase 1 実装済み)

---

## 1. 調査概要

### 1.1 比較対象リポジトリ

| リポジトリ | 規模 | 特徴 |
|-----------|------|------|
| **mims-harvard/ToolUniverse** | 1000+ ツール, 350+ SMCP 経由 | 実 API ラッパー群。IEDB, IMGT, MGnify, OBIS, CDC, CELLxGENE, ChIP-Atlas, BioModels 等の専門 DB を直接呼び出し可能 |
| **K-Dense-AI/claude-scientific-skills** | 140 スキル, ~17 カテゴリ | Claude Code 向け SKILL.md 形式。Scanpy, scvi-tools, AnnData, CellxGENE Census, COBRApy, Arboreto 等の解析パッケージスキルが充実 |

### 1.2 調査対象ドメイン (13 領域)

1. 環境科学・生態学
2. 疫学・公衆衛生
3. シングルセルゲノミクス
4. 空間トランスクリプトミクス
5. 科学テキストマイニング / NLP
6. マイクロバイオーム・メタゲノミクス
7. システム生物学
8. 集団遺伝学
9. 免疫情報学・ワクチンデザイン
10. ラジオミクス・デジタル病理 (既存 medical-imaging 拡張)
11. 食品科学・栄養学
12. 農業・植物科学
13. 臨床 NLP

---

## 2. ドメイン別ギャップマッピング

### 凡例
- ◎ = 強力なツール/スキル支援あり
- ○ = 中程度の支援あり
- △ = 限定的な支援
- ✗ = 支援なし

| # | ドメイン | ToolUniverse | K-Dense | SATORI 現状 | ギャップ度 |
|---|---------|-------------|---------|------------|----------|
| 1 | シングルセルゲノミクス | ◎ CELLxGENE (4ツール), HCA (2ツール) | ◎ Scanpy, scvi-tools (25+モデル), AnnData, Census, geniml | △ multi-omics が部分的 | **Critical** |
| 2 | 空間トランスクリプトミクス | ○ CELLxGENE 一部 | ◎ scvi-tools (DestVI, Stereoscope, Tangram, scVIVA) | ✗ 未カバー | **Critical** |
| 3 | 免疫情報学・ワクチンデザイン | ◎ IEDB (7+ツール), IMGT (3ツール), TheraSAbDab | △ なし | ✗ 未カバー | **Critical** |
| 4 | マイクロバイオーム・メタゲノミクス | ○ MGnify (2ツール) | △ なし | ✗ 未カバー | **High** |
| 5 | 疫学・公衆衛生 | ○ CDC, health_disparities (3ツール) | △ なし | △ survival-clinical が部分的 | **High** |
| 6 | 科学テキストマイニング / NLP | ○ PubMed, EuropePMC, LiteratureSynthesis, MedicalTermNormalizer | ○ OpenAlex, research-lookup, citation-management | △ deep-research が部分的 | **High** |
| 7 | システム生物学 | ○ BioModels, Reactome, KEGG, WikiPathways, Pathway Commons | ◎ COBRApy, Denario, HypoGeniC, LaminDB, Arboreto | △ network-analysis, multi-omics が部分的 | **Moderate-High** |
| 8 | 集団遺伝学 | ◎ gnomAD (5+ツール), GWAS (13+ツール), Ensembl (gene trees) | ○ gget, scikit-bio | △ variant-interpretation が部分的 | **Moderate** |
| 9 | 環境科学・生態学 | △ OBIS (2ツール) | ✗ なし | ✗ 未カバー | **Moderate** |
| 10 | 臨床 NLP | △ MedicalTermNormalizer, ICD, LOINC | △ clinical-decision-support | △ clinical-decision-support が部分的 | **Moderate** |
| 11 | 農業・植物科学 | △ IMPC (マウス表現型のみ) | ✗ なし | ✗ 未カバー | **Low** |
| 12 | 食品科学・栄養学 | ✗ なし | ✗ なし | ✗ 未カバー | **Low** |
| 13 | ラジオミクス・デジタル病理 | ○ 一部 | ○ imaging-data-commons, OMERO | ○ medical-imaging 既存 | **Already Covered** |

---

## 3. 新スキル提案 (10 スキル / 6 新カテゴリ)

### カテゴリ構成

| カテゴリ | 名称 (日本語) | スキル数 |
|---------|-------------|---------|
| **T** | シングルセル・空間オミクス | 2 |
| **U** | 免疫インフォマティクス | 2 |
| **V** | マイクロバイオーム・環境科学 | 2 |
| **W** | 疫学・公衆衛生 | 1 |
| **X** | 科学テキストマイニング | 1 |
| **Y** | システム生物学・集団遺伝学 | 2 |

---

### T-1. `scientific-single-cell-genomics`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | T. シングルセル・空間オミクス |
| **説明** | シングルセル RNA-seq データの前処理・クラスタリング・細胞型アノテーション・軌跡推定・差次発現解析の統合パイプライン。Scanpy/AnnData エコシステムをベースに、CELLxGENE Census (6100万+細胞) からの参照データ取得、scvi-tools による確率的モデリングを統合。 |
| **主要テクニック** | QC (doublet除去, Solo)、正規化 (scran)、HVG 選択、PCA/Harmony バッチ補正、Leiden クラスタリング、UMAP 可視化、CellAssign/scANVI 細胞型アノテーション、RNA velocity (scVelo/VeloVI)、擬時間解析、DGE (Wilcoxon/DESeq2-like)、GRN推定 (Arboreto) |
| **ToolUniverse SMCP ツール** | `CELLxGENE_get_expression_data`, `CELLxGENE_get_embeddings`, `CELLxGENE_get_gene_metadata`, `CELLxGENE_get_presence_matrix`, `hca_search_projects`, `hca_get_file_manifest` |
| **K-Dense 対応スキル** | Scanpy, AnnData, scvi-tools, CellxGENE Census, geniml/scEmbed, Arboreto |
| **正当性** | シングルセル解析は現代ゲノミクスの最重要技術。SATORI の既存 `multi-omics` ではバルク解析中心で、シングルセル特有のワークフロー (sparse matrix, 細胞型解像度, trajectory) が未カバー。ToolUniverse と K-Dense の両方で強力なサポートがあり、即座に統合可能。 |

---

### T-2. `scientific-spatial-transcriptomics`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | T. シングルセル・空間オミクス |
| **説明** | 空間トランスクリプトミクスデータ (10x Visium, MERFISH, Slide-seq, STARmap) の解析パイプライン。空間的遺伝子発現パターン、細胞タイプデコンボリューション、空間的変動遺伝子検出、リガンド-受容体相互作用マッピングを支援。 |
| **主要テクニック** | 空間的自己相関 (Moran's I, Geary's C)、空間的変動遺伝子検出 (SpatialDE, SPARK)、細胞タイプデコンボリューション (DestVI, Stereoscope, Tangram, Cell2location)、空間的近傍解析、リガンド-受容体ペア検出 (CellChat, NicheNet)、組織領域セグメンテーション、空間ドメイン検出 (BayesSpace)、scVIVA 可視化 |
| **ToolUniverse SMCP ツール** | `CELLxGENE_get_expression_data` (空間データ含む), `hca_search_projects` |
| **K-Dense 対応スキル** | scvi-tools (DestVI, Stereoscope, Tangram, scVIVA), Scanpy (squidpy 連携) |
| **正当性** | 空間トランスクリプトミクスは 2020年代最注目の技術 (Nature Method of the Year 2020)。組織内の遺伝子発現を空間的に解像でき、がん微小環境、神経回路、発生生物学に革命的インパクト。SATORI では完全に未カバーであり、K-Dense の scvi-tools 空間モデル群が即座に参照可能。 |

---

### U-1. `scientific-immunoinformatics`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | U. 免疫インフォマティクス |
| **説明** | 免疫エピトープ予測、MHC 結合親和性解析、B 細胞/T 細胞エピトープマッピング、免疫グロブリン/TCR レパトア解析、ワクチン候補抗原スクリーニングの統合パイプライン。IEDB・IMGT を中心データソースとし、免疫レパトアシーケンシング (AIRRseq) 解析も支援。 |
| **主要テクニック** | MHC-I/II 結合予測 (NetMHCpan)、B 細胞エピトープ予測 (BepiPred)、T 細胞エピトープマッピング、抗体 CDR 解析、VDJ リコンビネーション解析、免疫原性スコアリング、クロスリアクティビティ評価、人口カバレッジ算出 (HLA 多型)、抗体-抗原ドッキング |
| **ToolUniverse SMCP ツール** | `iedb_search_epitopes`, `iedb_search_bcell`, `iedb_search_mhc`, `iedb_search_antigens`, `iedb_get_epitope_mhc`, `iedb_get_epitope_antigens`, `iedb_get_epitope_references`, `iedb_search_references`, `IMGT_get_gene_info`, `IMGT_get_sequence`, `IMGT_search_genes`, `TheraSAbDab_*` (治療用抗体 DB) |
| **K-Dense 対応スキル** | なし (ギャップ — K-Dense でも未カバー) |
| **正当性** | COVID-19 以降、ワクチンデザイン・免疫療法への需要が急増。ToolUniverse に IEDB (7+ツール)、IMGT (3ツール)、TheraSAbDab という強力な免疫学 DB アクセスが存在し、SATORI が統合すれば市場で唯一のワクチンデザイン対応スキルセットとなる。K-Dense にも存在しない差別化ポイント。 |

---

### U-2. `scientific-infectious-disease`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | U. 免疫インフォマティクス |
| **説明** | 感染症アウトブレイク監視、病原体ゲノムサーベイランス、薬剤耐性予測、疫学モデリング (SIR/SEIR)、ドラッグリパーパシング統合パイプライン。CDC データ、WHO 報告、ゲノム配列データベースを横断的に統合。 |
| **主要テクニック** | 系統解析 (Nextstrain style)、AMR 遺伝子検出、変異モニタリング、疫学曲線フィッティング (SIR/SEIR/SIS)、基本再生産数 (R₀) 推定、薬剤感受性予測、パンデミック早期警告指標、接触追跡解析 |
| **ToolUniverse SMCP ツール** | `CDCRESTTool` (CDC 公衆衛生データ), `iedb_search_epitopes` (交差免疫), `NCBI_search_nucleotide`, `BLAST_nucleotide_search`, `kegg_search_pathway` (感染経路), `ChEMBL_search_drugs` (リパーパシング), `DrugBank_search` |
| **K-Dense 対応スキル** | なし (ToolUniverse に `tooluniverse-infectious-disease` スキルが存在) |
| **正当性** | ToolUniverse が既に `tooluniverse-infectious-disease` スキルを保有し、完全なフォールバックチェーン (Drug Search, Pathway Analysis 等) が文書化済み。パンデミック対応は公衆衛生の最重要課題であり、CDC ツールとの統合は SATORI にとって高インパクト。 |

---

### V-1. `scientific-microbiome-metagenomics`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | V. マイクロバイオーム・環境科学 |
| **説明** | 腸内細菌叢・皮膚/口腔マイクロバイオーム・環境メタゲノムの解析パイプライン。16S rRNA アンプリコンシーケンシングおよびショットガンメタゲノミクスのワークフローを統合。MGnify (EBI) を中心データソースとし、多様性解析・機能アノテーション・比較メタゲノミクスを支援。 |
| **主要テクニック** | OTU/ASV クラスタリング (DADA2, QIIME 2)、α/β 多様性解析 (Shannon, Simpson, UniFrac)、分類学的プロファイリング (MetaPhlAn, Kraken2)、機能アノテーション (HUMAnN3)、メタゲノムアセンブリ (MEGAHIT, metaSPAdes)、MAG (Metagenome-Assembled Genome) 回収、比較メタゲノミクス (LEfSe, ANCOM)、パスウェイ濃縮解析 |
| **ToolUniverse SMCP ツール** | `MGnify_search_studies`, `MGnify_list_analyses`, `kegg_search_pathway`, `kegg_get_pathway_info`, `Reactome_search_pathway`, `NCBI_search_nucleotide`, `BLAST_nucleotide_search` |
| **K-Dense 対応スキル** | なし (K-Dense でも未カバー — 差別化ポイント) |
| **正当性** | マイクロバイオーム研究は NIH Human Microbiome Project を皮切りに急拡大。腸脳軸、がん免疫療法応答、自己免疫疾患との関連が次々と報告。MGnify (EBI メタゲノミクス DB) がToolUniverse にあり、SATORI も K-Dense もカバーしていない完全な空白領域。 |

---

### V-2. `scientific-environmental-ecology`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | V. マイクロバイオーム・環境科学 |
| **説明** | 海洋・淡水・陸域の生物多様性データ解析、種分布モデリング、生態系サービス評価、気候変動影響予測の統合パイプライン。OBIS (海洋生物多様性情報システム) を中心データソースとし、GBIF 連携、生態学的ネットワーク解析、保全優先度評価を支援。 |
| **主要テクニック** | 種分布モデリング (MaxEnt, BRT)、生物多様性指標 (種の豊富さ、均等度、β多様性)、群集生態学解析 (NMDS, CCA/RDA)、系統的多様性 (PD, MPD)、食物網モデリング、eDNA メタバーコーディング、気候変動シナリオ (RCP/SSP) 応答予測、保全優先度評価 (Zonation)、ランドスケープ生態学指標 (FRAGSTATS) |
| **ToolUniverse SMCP ツール** | `OBIS_search_occurrences`, `OBIS_search_taxa`, `NCBI_search_nucleotide` (eDNA), `BLAST_nucleotide_search` (種同定), `kegg_list_organisms` |
| **K-Dense 対応スキル** | GeoPandas (空間データ), scikit-bio (多様性指標), ETE Toolkit (系統樹) |
| **正当性** | 生物多様性損失は国連 SDGs の重点課題。OBIS 海洋生物多様性データがToolUniverse にあり、K-Dense の地理空間スキル(GeoPandas)と組み合わせ可能。環境科学は両リポとも専用スキルが存在しない完全な空白領域。 |

---

### W-1. `scientific-epidemiology-public-health`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | W. 疫学・公衆衛生 |
| **説明** | 疫学調査デザイン・疾病サーベイランス・健康格差分析・公衆衛生介入評価の統合パイプライン。CDC 公衆衛生データ、健康格差指標 (SVI)、ICD/LOINC 医療用語標準化を統合し、エビデンスに基づく公衆衛生政策立案を支援。 |
| **主要テクニック** | 疫学研究デザイン (コホート、ケースコントロール、断面)、罹患率/死亡率/致死率算出、年齢調整率、生命表解析、空間疫学 (クラスター検出, Kulldorff SaTScan)、因果推論 (差の差法, 操作変数, 回帰不連続)、DAG ベース交絡調整、スクリーニング評価 (感度/特異度/LR)、費用対効果分析 (ICER/QALY) |
| **ToolUniverse SMCP ツール** | `CDCRESTTool`, `health_disparities_get_county_rankings_info`, `health_disparities_get_svi_info`, `icd_search_codes`, `loinc_search_codes`, `MedicalTermNormalizer`, `search_clinical_trials`, `MedlinePlus_search_topics_by_keyword` |
| **K-Dense 対応スキル** | statsmodels (回帰), scikit-survival (生存解析), GeoPandas (空間疫学) |
| **正当性** | SATORI の既存 `survival-clinical` はイベント時間解析に特化し、疫学固有のメトリクス (罹患率、年齢調整率、空間クラスター) やサーベイランスパイプラインが欠落。COVID-19 以降、公衆衛生データ解析の需要は恒久的に増大。ToolUniverse の CDC・健康格差ツールが即座に統合可能。 |

---

### X-1. `scientific-text-mining-nlp`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | X. 科学テキストマイニング |
| **説明** | 科学論文のテキストマイニング、バイオメディカル固有表現認識 (Bio-NER)、関係抽出、自動要約、知識グラフ構築の統合パイプライン。PubMed/EuropePMC の大規模文献コーパスから構造化知識を自動抽出し、仮説生成・文献レビュー自動化を支援。 |
| **主要テクニック** | バイオメディカル NER (遺伝子、疾患、薬物、変異)、関係抽出 (PPI、DDI、Gene-Disease)、知識グラフ構築 (異種ネットワーク)、トピックモデリング (BERTopic, LDA)、引用ネットワーク解析、エビデンスマッピング (PICO フレームワーク)、自動系統的レビュー支援、アブストラクトスクリーニング (ML ベース) |
| **ToolUniverse SMCP ツール** | `PubMed_search`, `EuropePMC_search`, `LiteratureSynthesisAgent`, `MedicalLiteratureReviewer`, `MedicalTermNormalizer`, `LiteratureSearchTool`, `MultiAgentLiteratureSearch`, `OpenAIRE_search_publications`, `ArXiv_search_papers` |
| **K-Dense 対応スキル** | OpenAlex, citation-management, research-lookup, literature-review |
| **正当性** | PubMed に年間 100 万件超の論文が追加される中、テキストマイニングは情報爆発への唯一の対処法。SATORI の `deep-research` は一般的な文献調査にとどまり、NER・関係抽出・知識グラフ等の構造化抽出が未実装。ToolUniverse の文献エージェント群と K-Dense の OpenAlex/文献管理が統合基盤に。 |

---

### Y-1. `scientific-systems-biology`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | Y. システム生物学・集団遺伝学 |
| **説明** | ゲノムスケール代謝モデリング (GEM)、遺伝子制御ネットワーク (GRN) 推定、シグナル伝達パスウェイモデリング、フラックスバランス解析 (FBA)、マルチオミクスデータ統合の統合パイプライン。BioModels・Reactome・KEGG を横断的に統合し、細胞レベルのシステム挙動を予測。 |
| **主要テクニック** | フラックスバランス解析 (FBA/pFBA)、フラックス変動解析 (FVA)、遺伝子ノックアウトシミュレーション、GRN 推定 (GENIE3, GRNBoost2/Arboreto)、常微分方程式 (ODE) モデリング、感度解析、パラメータ推定、SBML モデル操作、パスウェイ濃縮解析 (GSEA/ORA)、マルチオミクス因子分解 (MOFA+) |
| **ToolUniverse SMCP ツール** | `biomodels_search`, `biomodels_get_model`, `Reactome_map_uniprot_to_pathways`, `Reactome_search_pathway`, `kegg_search_pathway`, `kegg_get_pathway_info`, `kegg_find_genes`, `WikiPathways_search`, `pathway_commons_search`, `pathway_commons_get_interactions`, `intact_search_interactions`, `intact_get_interaction_network` |
| **K-Dense 対応スキル** | COBRApy (FBA), Denario (研究パイプライン), HypoGeniC (仮説生成), LaminDB (データ管理), Arboreto (GRN 推定) |
| **正当性** | SATORI の `network-analysis` はグラフ理論的解析に特化し、生物学的パスウェイモデリングや FBA は未実装。`multi-omics` もデータ統合に主眼でモデル駆動型解析が欠落。ToolUniverse の BioModels/Reactome/KEGG/IntAct ツール群と K-Dense の COBRApy が直接的な基盤。 |

---

### Y-2. `scientific-population-genetics`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | Y. システム生物学・集団遺伝学 |
| **説明** | 集団遺伝学解析パイプライン。アレル頻度推定、Hardy-Weinberg 平衡検定、LD 解析、自然選択検出、集団構造推定、系統地理学解析を統合。gnomAD (14万+ゲノム) と GWAS Catalog の大規模データを活用し、遺伝的多様性と進化の理解を支援。 |
| **主要テクニック** | アレル頻度推定・集団分化 (F_ST, D_ST)、LD 構造解析 (r², D')、自然選択統計量 (Tajima's D, iHS, nSL, XP-EHH)、集団構造推定 (PCA, ADMIXTURE)、祖先推定 (ローカルアダプテーション)、人口統計推定 (PSMC, MSMC)、系統樹推定 (最尤法, ベイズ)、GWAS メタアナリシス、多遺伝子リスクスコア (PRS) 算出、メンデルランダム化 (MR) |
| **ToolUniverse SMCP ツール** | `gnomad_get_gene`, `gnomad_get_gene_constraints`, `gnomad_get_region`, `gnomad_get_transcript`, `gnomad_get_variant`, `gnomad_search_genes`, `gnomad_search_variants`, `gwas_search_studies`, `gwas_search_associations`, `gwas_get_associations_for_snp`, `gwas_get_snps_for_gene`, `gwas_get_variants_for_trait`, `ensembl_get_genetree`, `ensembl_get_variants` |
| **K-Dense 対応スキル** | gget (Ensembl ラッパー), scikit-bio (多様性指標), ETE Toolkit (系統樹) |
| **正当性** | SATORI の `variant-interpretation` は臨床的バリアント解釈に特化し、集団レベルの進化遺伝学 (F_ST, 選択圧, 人口統計) が未カバー。ToolUniverse の gnomAD (5+ツール) と GWAS (13+ツール) は集団遺伝学に最適な API 群。PRS 算出は個別化医療と直結し臨床需要も高い。 |

---

## 4. ToolUniverse SMCP ツール詳細マッピング

新スキルに必要な ToolUniverse 専用ツールの完全一覧:

### 免疫学ドメイン (U-1, U-2)
```
# IEDB (Immune Epitope Database) — 7+ ツール
iedb_search_epitopes          # エピトープ検索 (配列, 構造タイプ)
iedb_search_bcell             # B 細胞エピトープ検索
iedb_search_mhc               # MHC 拘束性アッセイ検索
iedb_search_antigens           # 抗原検索
iedb_get_epitope_mhc           # エピトープ-MHC 関連データ取得
iedb_get_epitope_antigens      # エピトープ-抗原関連データ取得
iedb_get_epitope_references    # エピトープ-文献関連データ取得
iedb_search_references         # 免疫学文献検索

# IMGT (Immunogenetics) — 3 ツール
IMGT_get_gene_info             # 免疫グロブリン/TCR 遺伝子情報
IMGT_get_sequence              # 配列取得 (FASTA/EMBL)
IMGT_search_genes              # 遺伝子検索 (IGHV, TRAV 等)

# CDC — 公衆衛生データ
CDCRESTTool                    # CDC 公衆衛生統計 REST API
```

### シングルセルドメイン (T-1, T-2)
```
# CELLxGENE Census — 4 ツール
CELLxGENE_get_expression_data  # 発現データ取得
CELLxGENE_get_embeddings       # 細胞埋め込み取得
CELLxGENE_get_gene_metadata    # 遺伝子メタデータ取得
CELLxGENE_get_presence_matrix  # 存在マトリクス取得

# Human Cell Atlas — 2 ツール
hca_search_projects            # HCA プロジェクト検索
hca_get_file_manifest          # ファイルマニフェスト取得
```

### マイクロバイオームドメイン (V-1)
```
# MGnify (EBI Metagenomics) — 2 ツール
MGnify_search_studies          # メタゲノミクス研究検索 (バイオーム/キーワード)
MGnify_list_analyses           # 研究のアナリシス一覧取得
```

### 生態学ドメイン (V-2)
```
# OBIS (Ocean Biodiversity) — 2 ツール
OBIS_search_occurrences        # 海洋生物出現記録検索
OBIS_search_taxa               # 分類群検索
```

### 集団遺伝学ドメイン (Y-2)
```
# gnomAD — 7 ツール
gnomad_get_gene                # 遺伝子アレル頻度
gnomad_get_gene_constraints    # 制約メトリクス (pLI, Z スコア)
gnomad_get_region              # 領域バリアント
gnomad_get_transcript          # トランスクリプト
gnomad_get_variant             # バリアント詳細
gnomad_search_genes            # 遺伝子検索
gnomad_search_variants         # バリアント検索

# GWAS Catalog — 13+ ツール
gwas_search_studies            # GWAS 研究検索
gwas_search_associations       # 関連検索
gwas_get_associations_for_snp  # SNP 関連取得
gwas_get_associations_for_study # 研究別関連
gwas_get_associations_for_trait # 形質別関連
gwas_get_snp_by_id             # SNP 詳細
gwas_get_snps_for_gene         # 遺伝子 SNP
gwas_get_studies_for_trait     # 形質別研究
gwas_get_study_by_id           # 研究詳細
gwas_get_variants_for_trait    # 形質バリアント
gwas_search_snps               # SNP 検索
gwas_get_association_by_id     # 関連詳細
```

### システム生物学ドメイン (Y-1)
```
# BioModels — 2+ ツール
biomodels_search               # SBML モデル検索
biomodels_get_model            # モデル詳細取得

# Pathway DB — 5+ ツール
Reactome_map_uniprot_to_pathways
Reactome_search_pathway
kegg_search_pathway
kegg_get_pathway_info
WikiPathways_search
pathway_commons_search
pathway_commons_get_interactions

# IntAct (分子相互作用) — 8 ツール
intact_search_interactions
intact_get_interactions
intact_get_interaction_details
intact_get_interaction_network
intact_get_interactions_by_complex
intact_get_interactions_by_organism
intact_get_complex_details
intact_get_interactor
```

---

## 5. 実装優先度マトリクス

| 順位 | スキル | ギャップ度 | ツール支援度 | 科学的インパクト | 差別化度 | 総合スコア |
|------|--------|-----------|------------|---------------|---------|-----------|
| **1** | T-1 single-cell-genomics | Critical | ◎◎ (TU+KD) | ◎ | ○ | **9.5/10** |
| **2** | U-1 immunoinformatics | Critical | ◎ (TU 12+ツール) | ◎ | ◎ | **9.5/10** |
| **3** | T-2 spatial-transcriptomics | Critical | ○◎ (KD 豊富) | ◎ | ◎ | **9.0/10** |
| **4** | V-1 microbiome-metagenomics | High | ○ (TU MGnify) | ◎ | ◎ | **8.5/10** |
| **5** | Y-1 systems-biology | Moderate-High | ◎◎ (TU+KD) | ◎ | ○ | **8.0/10** |
| **6** | W-1 epidemiology-public-health | High | ○ (TU CDC) | ◎ | ○ | **8.0/10** |
| **7** | X-1 text-mining-nlp | High | ○◎ (TU+KD) | ○ | ○ | **7.5/10** |
| **8** | U-2 infectious-disease | High | ○ (TU 豊富) | ◎ | ○ | **7.5/10** |
| **9** | Y-2 population-genetics | Moderate | ◎ (TU 20+ツール) | ○ | ○ | **7.0/10** |
| **10** | V-2 environmental-ecology | Moderate | △ (TU OBIS) | ○ | ◎ | **6.5/10** |

### 推奨実装フェーズ

| フェーズ | スキル | 目標バージョン |
|---------|--------|--------------|
| Phase 2A (即時) | T-1, U-1, T-2 | v0.11.0 |
| Phase 2B (次期) | V-1, Y-1, W-1 | v0.12.0 |
| Phase 2C (後期) | X-1, U-2, Y-2, V-2 | v0.13.0 |

---

## 6. カバレッジ予測

| 指標 | 現在 (v0.10.0) | Phase 2A 後 | Phase 2C 後 |
|------|---------------|-------------|-------------|
| スキル数 | 56 | 59 | 66 |
| カテゴリ数 | 19 (A–S) | 21 (A–U) | 25 (A–Y) |
| ToolUniverse カバー率 | ~60% | ~75% | ~85% |
| K-Dense カバー率 | ~45% | ~60% | ~70% |
| 未カバー主要ドメイン | 13 | 7 | 3 |

### Phase 2C 後も残る主要ギャップ

| ドメイン | 理由 |
|---------|------|
| 農業・植物科学 | 両リポにツール/スキル支援がほぼ皆無 |
| 食品科学・栄養学 | 両リポにツール/スキル支援が皆無 |
| 臨床 NLP | 専用ツールが限定的 (MedicalTermNormalizer のみ) |

---

## 7. 前回分析 (Tier 2/3) との関係

前回の GAP_ANALYSIS.md で Tier 2 として挙げた以下のギャップは引き続き有効だが、本分析では上記 10 スキルを優先:

| 前回 Tier 2 項目 | 状態 | 関連性 |
|-----------------|------|--------|
| proteomics-ms | 未実装 | 独立ギャップ (次回分析で再評価) |
| metabolic-modeling | 未実装 | → Y-1 systems-biology に統合可能 |
| neuroscience | 未実装 | 独立ギャップ |
| multi-objective-optimization | 未実装 | 独立ギャップ |
| discrete-event-simulation | 未実装 | 独立ギャップ |
| large-scale-data | 未実装 | 独立ギャップ |
| schematics-generation | 未実装 | 独立ギャップ |
| ELN/LIMS integration | 未実装 | 独立ギャップ |
| cloud-infrastructure | 未実装 | 独立ギャップ |

---

## 8. 結論

SATORI v0.10.0 は 56 スキルで生命科学・創薬・臨床を良好にカバーするが、**シングルセルゲノミクス**、**免疫情報学**、**空間トランスクリプトミクス**の 3 領域が Critical なギャップとして特定された。これらは:

1. **科学的インパクト**: 2020 年代の生命科学で最も急成長している技術分野
2. **ツール支援**: ToolUniverse (CELLxGENE, IEDB, IMGT) と K-Dense (Scanpy, scvi-tools) の両方で強力なサポートが存在
3. **差別化**: K-Dense にも存在しない免疫情報学スキルは市場唯一の差別化ポイント

10 スキルの段階的実装により、v0.13.0 で 66 スキル / 25 カテゴリ (A–Y) に到達し、科学 AI ツールとしての包括性を大幅に向上できる。
