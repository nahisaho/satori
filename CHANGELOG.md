# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.11.0] - 2026-02-16

### Added
- **10 新スキル (既存 9 カテゴリ拡張)** を追加 — ToolUniverse & K-Dense-AI/claude-scientific-skills ギャップ分析に基づく Phase 3 ドメイン拡張 (66→76 スキル)
- **E. 信号・スペクトル・時系列（3→4 種）**
  - **scientific-neuroscience-electrophysiology** スキル (#67): SpikeInterface/Kilosort4 スパイクソート・Allen Institute 品質基準・MNE-Python EEG (前処理/ERP/マイクロステート)・NeuroKit2 ECG HRV/EDA 解析・脳機能結合 (wPLI/グラフ理論)
- **F. 生命科学・オミクス（5→7 種）**
  - **scientific-proteomics-mass-spectrometry** スキル (#68): pyOpenMS LC-MS/MS 前処理・ペプチド ID (FDR)・タンパク質定量 (LFQ/TMT/SILAC/iBAQ)・PTM マッピング・matchms スペクトル類似度・GNPS 分子ネットワーク・PRIDE/UniProt/KEGG/Reactome ツール連携
  - **scientific-gene-expression-transcriptomics** スキル (#69): GEO データセット取得 (GEOparse)・PyDESeq2 差次発現解析・GTEx 組織発現/eQTL 照会・GSEA/ORA 遺伝子セット濃縮解析・geo/GTEx/ExpressionAtlas/ArrayExpress (24 ツール) 連携
- **G. 化学・材料・イメージング（3→4 種）**
  - **scientific-computational-materials** スキル (#70): pymatgen 結晶構造操作・対称性解析・Materials Project API 照会・相図凸包解析・電子バンド構造/DOS 可視化・VASP/QE 入出力
- **H. 臨床・疫学・メタ科学（3→4 種）**
  - **scientific-clinical-trials-analytics** スキル (#71): ClinicalTrials.gov API v2 多基準検索・試験詳細取得・競合ランドスケープ解析・AE/アウトカム抽出・バルクエクスポート・ClinicalTrials (7)/FDA (2) ツール連携
- **M. 実験室自動化・データ管理（1→2 種）**
  - **scientific-lab-data-management** スキル (#72): Benchling ELN/DNA 設計/レジストリ・DNAnexus ゲノミクス PaaS・OMERO バイオイメージング管理・Protocols.io プロトコル共有
- **N. 科学プレゼンテーション・図式（1→2 種）**
  - **scientific-scientific-schematics** スキル (#73): CONSORT フロー図・NN アーキテクチャ図・分子パスウェイ図 (Mermaid)・TikZ 出版品質ベクター図
- **O. 研究計画・グラント・規制（2→3 種）**
  - **scientific-regulatory-science** スキル (#74): FDA Orange Book 承認履歴/特許/排他性・510(k) 医療機器クリアランス・ISO 13485 設計管理/CAPA・USPTO 特許検索・FDA_OrangeBook (6)/FAERS ツール連携
- **P. ファーマコビジランス・薬理ゲノミクス（1→2 種）**
  - **scientific-pharmacogenomics** スキル (#75): PharmGKB 遺伝子-薬物相互作用・CPIC ガイドライン・Star アレルアノテーション・代謝型判定 (PM/IM/NM/RM/UM)・FDA PGx バイオマーカー・PharmGKB (7)/FDA (3)/OpenTargets (1) ツール連携
- **T. シングルセル・空間・エピゲノミクス（2→3 種）**
  - **scientific-epigenomics-chromatin** スキル (#76): ChIP-seq MACS2/MACS3 ピークコール・ATAC-seq NFR 検出・WGBS/RRBS DMR 検出・ChromHMM 15 状態モデル・Hi-C TAD/A-B コンパートメント・DiffBind 差次的結合・ChIPAtlas (4)/JASPAR (3)/SCREEN (1)/ENCODE (3) ツール連携

### Changed
- **README.md**: スキル数を 66→76 に更新、カテゴリ数 26 (変更なし, 既存カテゴリ拡張)、ToolUniverse 連携スキル数を 32→42 に更新、各カテゴリのスキル数・説明・ディレクトリ構造・ファイル I/O テーブルを更新
- **docs/GAP_ANALYSIS_v0.11.md**: v0.11.0 ギャップ分析ドキュメント新規追加

## [0.10.0] - 2026-02-15

### Added
- **7 新カテゴリ (T-Z)・10 新スキル** を追加 — ToolUniverse ギャップ分析に基づくドメイン拡張
- **T. シングルセル・空間オミクス（2 種）**
  - **scientific-single-cell-genomics** スキル (#57): scRNA-seq QC (Scanpy)・Leiden クラスタリング・DEG (Wilcoxon)・RNA velocity (scVelo)・CellChat 細胞間通信・CELLxGENE/HCA ツール連携
  - **scientific-spatial-transcriptomics** スキル (#58): Visium/MERFISH 前処理 (Squidpy)・SVG 検出 (Moran's I)・空間ドメイン検出・cell2location デコンボリューション・空間 L-R 解析
- **U. 免疫・感染症（2 種）**
  - **scientific-immunoinformatics** スキル (#59): MHC-I/II 結合予測 (MHCflurry)・B 細胞エピトープ・TCR/BCR レパトア多様性・抗体 CDR 解析 (ANARCI)・ワクチン候補ランキング・IEDB/IMGT/SAbDab ツール連携
  - **scientific-infectious-disease** スキル (#60): 病原体 WGS QC・AMR 遺伝子検出 (ResFinder/RGI)・MLST/cgMLST 型別・系統解析 (IQ-TREE)・SIR/SEIR 数理モデル・CDC/EUHealthInfo ツール連携
- **V. マイクロバイオーム・環境（2 種）**
  - **scientific-microbiome-metagenomics** スキル (#61): DADA2 ASV パイプライン・MetaPhlAn/Kraken2 分類プロファイリング・α/β 多様性 (Shannon/Bray-Curtis/PERMANOVA)・ANCOM-BC 差次的存在量・HUMAnN 機能プロファイリング・MGnify/KEGG ツール連携
  - **scientific-environmental-ecology** スキル (#62): SDM (MaxEnt/RF/GBM)・生物多様性指数 (α/β/γ)・群集序列化 (NMDS/CCA/RDA)・保全優先順位ランキング・OBIS/GBIF ツール連携
- **W. システム生物学（1 種）**
  - **scientific-systems-biology** スキル (#63): SBML/RoadRunner 動的シミュレーション・定常状態/ヤコビアン解析・FBA/pFBA (cobrapy)・GRN 推定 (GENIE3)・Sobol 感度解析 (SALib)・BioModels/Reactome/BiGG ツール連携
- **X. 疫学・公衆衛生（1 種）**
  - **scientific-epidemiology-public-health** スキル (#64): RR/OR/RD/NNT/AF リスク指標・直接/間接年齢標準化 (SMR)・LISA/Getis-Ord 空間クラスタリング・DAG バックドア基準交絡分析・WHO/CDC/EUHealthInfo ツール連携
- **Y. 集団遺伝学（1 種）**
  - **scientific-population-genetics** スキル (#65): PLINK2 遺伝子型 QC・HWE 検定・PCA/ADMIXTURE 集団構造推定・Weir-Cockerham Fst・iHS/Tajima's D 選択スキャン・gnomAD/GWAS Catalog ツール連携
- **Z. 科学テキストマイニング（1 種）**
  - **scientific-text-mining-nlp** スキル (#66): BioBERT/SciSpaCy NER・関係抽出 (PPI/DDI/GDA)・知識グラフ構築 (Louvain)・BERTopic トピックモデリング・引用ネットワーク分析 (PageRank/HITS)・PubTator/PubMed/EuropePMC ツール連携

### Changed
- **README.md**: スキル数を 56→66 に更新、カテゴリ数を 19→26 に拡張 (T-Z 追加)、ToolUniverse 連携スキル数を 22→32 に更新、次世代オミクス・疫学パイプラインフロー図を追加、ディレクトリ構造に 7 新カテゴリを追加

## [0.9.0] - 2026-02-14

### Added
- **ToolUniverse MCP ツール連携**: 22 スキル（HIGH 13 + MEDIUM 9）の SKILL.md に `### 利用可能ツール` セクションを追加 — [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で 1,200 以上の外部科学データベースツール（ClinVar, gnomAD, OncoKB, CIViC, FAERS, UniProt, ChEMBL, PubMed 等）への参照を記載
  - **HIGH（13 スキル）**: variant-interpretation, precision-oncology, disease-research, drug-target-profiling, drug-repurposing, pharmacovigilance, clinical-decision-support, protein-structure-analysis, admet-pharmacokinetics, protein-design, bioinformatics, multi-omics, metabolomics
  - **MEDIUM（9 スキル）**: deep-research, cheminformatics, sequence-analysis, citation-checker, meta-analysis, network-analysis, graph-neural-networks, survival-clinical, grant-writing

### Changed
- **README.md**: Overview に ToolUniverse MCP 連携の説明を追加、ファイル入出力セクションに「ToolUniverse MCP ツール連携」サブセクションを追加（アーキテクチャ図付き）

## [0.8.0] - 2026-02-13

### Added
- **4 新カテゴリ (P-S)・9 新スキル** を追加 — ToolUniverse (mims-harvard) と claude-scientific-skills (K-Dense-AI) のギャップ分析 Tier-1 に基づく MECE 設計
- **P. ファーマコビジランス（1 種）**
  - **scientific-pharmacovigilance** スキル (#48): FAERS 不均衡分析 (PRR/ROR/IC/EBGM)・MedDRA 階層構造 (LLT→PT→HLT→HLGT→SOC)・時系列トレンド・人口統計層別化・Naranjo 因果評価スケール・安全性シグナルレポート生成
- **Q. 腫瘍学・疾患研究（2 種）**
  - **scientific-precision-oncology** スキル (#49): CIViC GraphQL API・OncoKB Annotation API・cBioPortal TCGA 変異頻度・TMB/MSI 判定・AMP/ASCO/CAP Tiering (I-IV)・分子腫瘍ボード (MTB) レポート生成
  - **scientific-disease-research** スキル (#50): GWAS Catalog (EBI) 検索・DisGeNET GDA スコア・Orphanet 希少疾患カタログ・HPO 表現型マッチング・Polygenic Risk Score (PRS) 算出
- **R. 量子・先端計算（5 種）**
  - **scientific-quantum-computing** スキル (#51): Qiskit VQE (H₂ 基底エネルギー)・PennyLane 量子 ML (VQC)・Cirq ノイズモデリング・QAOA MaxCut・QuTiP Jaynes-Cummings ダイナミクス
  - **scientific-graph-neural-networks** スキル (#52): PyG GCN/GAT/GIN 分子特性予測・TorchDrug MoleculeNet ベンチマーク・知識グラフ推論 (TransE/RotatE)・Scaffold Split・GNNExplainer
  - **scientific-bayesian-statistics** スキル (#53): PyMC/Stan 階層ベイズ回帰・NUTS MCMC サンプリング・事後予測チェック (PPC)・WAIC/LOO-CV モデル比較・Gaussian Process ベイズ最適化
  - **scientific-explainable-ai** スキル (#54): TreeSHAP/KernelSHAP/DeepSHAP 特徴量寄与分解・LIME ローカル説明・Captum Integrated Gradients・反実仮想説明・公平性監査
  - **scientific-deep-learning** スキル (#55): PyTorch Lightning 構造化パイプライン・timm 事前学習 CNN/ViT・HuggingFace Transformer Fine-tune・Optuna HPO (Hyperband Pruning)・ONNX/TorchScript エクスポート
- **S. 医用イメージング（1 種）**
  - **scientific-medical-imaging** スキル (#56): DICOM/NIfTI 読み込み・HU ウィンドウ処理・Isotropic リサンプリング・MONAI U-Net/SwinUNETR 3D セグメンテーション・WSI パッチ抽出 (openslide)・PyRadiomics 特徴量抽出

### Changed
- **README.md**: スキル数を 47→56 に更新、カテゴリ数を 15→19 に拡張 (P-S 追加)、先端計算・医用パイプラインフロー図を追加、ディレクトリ構造に 4 新カテゴリを追加

## [0.7.1] - 2026-02-13

### Fixed
- J-O カテゴリ 11 スキル全てに `## References` セクション（Output Files + 参照スキル テーブル）を追加 — パイプライン連携の欠落を修正

## [0.7.0] - 2026-02-13

### Added
- **6 新カテゴリ (J-O)・11 新スキル** を追加 — ToolUniverse (mims-harvard) と claude-scientific-skills (K-Dense-AI) のギャップ分析に基づく MECE 設計
- **J. 創薬・ファーマコロジー（3 種）**
  - **scientific-drug-target-profiling** スキル (#37): 9-path 標的プロファイリング (Identity→Safety)・TDL 分類 (Tclin/Tchem/Tbio/Tdark)・ドラッガビリティマトリクス (SM/Ab/PROTAC/ASO/遺伝子治療)・競合ランドスケープ
  - **scientific-admet-pharmacokinetics** スキル (#38): 5 段階 ADMET パイプライン (吸収→分布→代謝→排泄→毒性)・Lipinski/Veber/Ghose/QED ドラッグライクネス・PAINS/Brenk 構造アラート・CYP 阻害/基質予測・1-コンパートメント PK モデル
  - **scientific-drug-repurposing** スキル (#39): 7 戦略リポジショニング (Target/Compound/Disease/Mechanism/Network/Phenotype/Structure)・ネットワーク近接解析 (Guney et al.)・6 因子重み付き多基準スコアリング
- **K. 構造生物学・タンパク質工学（2 種）**
  - **scientific-protein-structure-analysis** スキル (#40): PDB Search API v2・AlphaFold DB 取得・品質評価 (Resolution/R-factor/pLDDT)・結合サイト検出 (fpocket/PDBe)
  - **scientific-protein-design** スキル (#41): ESM-2 変異スキャン (ΔLL)・de novo 設計パイプライン (RFdiffusion→ProteinMPNN→ESMFold)・バインダー/スキャフォールド/酵素設計・検証基準 (pLDDT>70, pTM>0.5)
- **L. 精密医療・臨床意思決定（2 種）**
  - **scientific-variant-interpretation** スキル (#42): ACMG/AMP 28 基準分類アルゴリズム・薬理ゲノミクス (CYP2D6/2C19/2C9/DPYD/TPMT/HLA-B, CPIC ガイドライン)・体細胞変異 OncoKB レベル・in silico 予測 (CADD/REVEL/AlphaMissense/SpliceAI)
  - **scientific-clinical-decision-support** スキル (#43): GRADE エビデンス枠組 (1a-5)・推奨グレード (A-D)・精密腫瘍学ワークフロー (OncoKB)・ClinicalTrials.gov API マッチング・リスク-ベネフィット分析
- **M. 実験室自動化（1 種）**
  - **scientific-lab-automation** スキル (#44): PyLabRobot ユニバーサル API・Opentrons OT-2 プロトコル・SOP テンプレート・Protocols.io API・ELN/LIMS 連携・QC 検証ワークフロー
- **N. 科学プレゼンテーション（1 種）**
  - **scientific-presentation-design** スキル (#45): 15 スライド口頭発表テンプレート (時間配分付)・tikzposter テンプレート・matplotlib ワークフロー図・3m ルール・アクセシビリティ (コントラスト比 ≥4.5)
- **O. 研究計画・グラント（2 種）**
  - **scientific-grant-writing** スキル (#46): NIH Specific Aims 1 ページテンプレート・Research Strategy (Significance/Innovation/Approach 12 頁)・JSPS 科研費フォーマット・予算計画・Budget Justification
  - **scientific-research-methodology** スキル (#47): SCAMPER/Cross-Domain/TRIZ ブレインストーミング・仮説型分類・研究デザインマトリクス・FINER 基準・バイアス分類 (6 種)・IRB/倫理チェックリスト

### Changed
- **README.md**: スキル数を 36→47 に更新、カテゴリ数を 9→15 に拡張 (J-O 追加)、ディレクトリ構造に 6 新カテゴリを追加

## [0.6.0] - 2026-02-12

### Added
- **scientific-deep-research** スキル (#36): SHIKIGAMI の WebResearcher パラダイムを科学研究に適応した深層リサーチスキル
  - Think→Search→Evaluate→Synthesize 反復サイクル（最大 15 ラウンド）
  - 学術データベース検索統合（PubMed, Google Scholar, arXiv, Semantic Scholar, CiNii, J-STAGE 等）
  - エビデンス階層評価（Level 1a〜5 + プレプリント）
  - ソース信頼性スコアリング（IF, h-index, サンプルサイズ, 統計手法, 再現性, 出版年）
  - ハルシネーション防止マーキング（✅ 検証済 / 📎 単一ソース / ⚠️ 未査読 / ❓ AI推定 / 🔄 古いデータ / ⚡ 矛盾）
  - 交差検証（数値乖離率判定・内容矛盾処理）
  - 品質ゲート（Phase 1→2, 2→3, 3→4, 4 完了）
  - PICO/PECO/SPIDER 構造化・検索戦略設計
  - PRISMA 2020 フローテンプレート（Systematic Review 用）
  - 批判的評価チェックリスト（Cochrane RoB / NOS / AMSTAR-2 準拠）
  - 日英並列検索必須・学術ドメイン優先
  - 分野別検索戦略テンプレート（生命科学・材料科学・CS/AI）
  - 出力ファイル: research_report.md, evidence_table.json, search_log.md, source_registry.json, prisma_flow.md

### Changed
- **README.md**: スキル数を 35→36 に更新、カテゴリ I (Deep Research) を追加、ディレクトリ構造に [I] Deep Research を追加

## [0.5.2] - 2026-02-12

### Fixed
- **scientific-academic-writing** スキル: SATORI バージョン参照ルールを追加 — 論文原稿内でバージョン番号をハードコードせず `package.json` から動的取得するよう明記
- **scientific-academic-writing** スキル: AI 使用開示 (AI Usage Disclosure) セクションを追加 — ジャーナル別記載場所・テンプレート・Authors セクション禁止ルールを整備

## [0.5.0] - 2026-02-12

### Added
- **scientific-peer-review-response** スキル (#11): 査読コメント構造化 (Major/Minor/Editorial 自動分類)・ポイントバイポイント回答レター生成・改訂版カバーレター・コメント-改訂マッピング
- **scientific-revision-tracker** スキル (#12): 改訂履歴追跡・バージョンスナップショット・差分検出・変更マークアップ (Markdown/LaTeX/プレーンテキスト)・トレーサビリティ検証
- **scientific-paper-quality** スキル (#13): 可読性スコア (Flesch-Kincaid/Gunning Fog)・IMRAD バランス分析・語彙品質 (冗長表現/過剰主張検出)・ジャーナル適合性チェック・再現可能性チェック・総合品質スコアカード

### Changed
- **README.md**: スキル数を 32→35 に更新、カテゴリ A を 10→13 種に拡張、パイプラインフローに査読対応・改訂追跡・品質評価フェーズを追加
- スキル番号を再付番 (B-H カテゴリ: #11-32 → #14-35)

## [0.4.0] - 2025-06-13

### Added
- **scientific-supplementary-generator** スキル (#8): Supplementary Information 自動生成・SI 図表整理・本文-SI 相互参照検証・ジャーナル別 SI フォーマット対応 (Nature/Science/ACS/IEEE/Elsevier)
- **scientific-latex-formatter** スキル (#9): Markdown→LaTeX 変換・ジャーナルテンプレート適用 (Nature/Science/ACS/IEEE/Elsevier/RevTeX)・数式/図表/引用の LaTeX 構文変換・BibTeX 生成
- **scientific-citation-checker** スキル (#10): 引用文献の自動検索・網羅性チェック・整合性検証 (孤立参考文献/未解決引用/重複検出)・引用カバレッジ分析

### Changed
- **README.md**: スキル数を 29→32 に更新、カテゴリ A を 7→10 種に拡張、パイプラインフローテーブルに引用検証・SI 生成・LaTeX 変換フェーズを追加
- スキル番号を再付番 (B-H カテゴリ: #8-29 → #11-32)

## [0.3.0] - 2026-02-11

### Added
- **scientific-critical-review** スキル (#29): 草稿の批判的レビュー・5パスレビューワークフロー・7層考察深化フレームワーク・修正案生成
- **scientific-hypothesis-pipeline**: 仮説定義の永続化 (`docs/hypothesis.{md,json}`)、ワークフロー設計の永続化 (`docs/workflow_design.{md,json}`)、`load_hypothesis()` / `load_workflow_design()` 関数
- **scientific-critical-review**: レビュー結果の永続化 (`manuscript/review_report.{md,json}`, `manuscript/manuscript_revised.md`, `manuscript/revision_diff.md`)、`run_review_pipeline()` 統合関数
- **scientific-pipeline-scaffold**: Step 0 (仮説・ワークフロー保存) / Step 8 (論文執筆・レビュー) を追加、参照スキル表を追加
- **scientific-academic-writing**: 参照スキル表・パイプラインフロー図を追加、`docs/hypothesis.md` / `results/analysis_summary.json` からの自動参照を明記
- 4スキル間の双方向相互参照を完備（hypothesis-pipeline ↔ critical-review ↔ academic-writing ↔ pipeline-scaffold）

### Fixed
- **scientific-critical-review**: 全 save 関数のパス方式を `Path()` (CWD相対) から `BASE_DIR` ベースに統一
- **README.md**: パイプラインフロー図追加、インストール方法追加、ディレクトリ構造に hypothesis-pipeline / critical-review を追加

## [0.2.0] - 2025-06-12

### Added
- **scientific-academic-writing**: 図の埋め込みワークフロー（セクション 8）を追加
- 全 6 ジャーナルテンプレートに `![Figure N](figures/...)` 構文を追加
  - IMRaD 標準、Nature、Science、ACS、IEEE、Elsevier

## [0.1.0] - 2025-06-12

### Added
- 初回リリース: 27 スキルを含む Agent Skills コレクション
- CLI ツール (`npx @nahisaho/satori init`) による `.github/skills/` のインストール
- `--force` / `--dry-run` オプション対応
- 8 カテゴリ (A-H) のスキル分類体系
- 7 種のジャーナルテンプレート (IMRaD, Nature, Science, ACS, IEEE, Elsevier, Qiita)

[0.8.0]: https://github.com/nahisaho/satori/compare/v0.7.1...v0.8.0
[0.7.1]: https://github.com/nahisaho/satori/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/nahisaho/satori/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/nahisaho/satori/compare/v0.5.2...v0.6.0
[0.5.2]: https://github.com/nahisaho/satori/compare/v0.5.0...v0.5.2
[0.5.0]: https://github.com/nahisaho/satori/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/nahisaho/satori/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/nahisaho/satori/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/nahisaho/satori/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/nahisaho/satori/releases/tag/v0.1.0
