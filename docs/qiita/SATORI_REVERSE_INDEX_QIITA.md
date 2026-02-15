---
title: 【SATORI v0.26.0】190スキル逆引き辞書 完全索引
tags:
  - バイオインフォマティクス
  - 機械学習
  - パイプライン
  - 創薬
  - データサイエンス
private: false
updated_at: ''
id: null
organization_url_name: null
slide: false
ignorePublish: false
---

# はじめに

**[SATORI](https://github.com/nahisaho/satori)** は **GitHub Copilot** 上で動作する、190 の専門スキルを組み合わせて仮説構築から論文出版まで、あらゆる科学研究ワークフローを自動化するパイプラインフレームワークです。

本記事では **190 スキルの逆引き辞書（Reverse Index）** を提供します。ToolUniverse API キー、成果物パス、参照関係、トリガーフレーズから目的のスキルを素早く検索できます。

## SATORI とは

**SATORI（悟り）** は禅仏教における「直観的な覚醒（awakening）」を意味する。段階的な学習（漸悟）ではなく、一瞬にしてすべてを理解する体験だ。

この名前には二つの意図がある。

**第一に、AI Agent への知識注入のメタファー。** 人間の科学者は何年もかけてドメイン知識を蓄積するが、`npm install @nahisaho/satori` の一行で、AI Agent は 190 のスキルを一瞬で「悟る」。教科書を 10 年読む代わりに、構造化されたスキルファイルが Agent に「悟り」を与える。

**第二に、バクロニム（逆頭字語）としての設計。**

***S**cientific **A**nalysis **T**oolkit for **O**rganized **R**esearch **I**ntelligence*

「科学分析のための、体系化された研究知能ツールキット」。名前が機能を説明し、機能が名前を体現する。

## なぜ「悟り」なのか

[AI Agent が代替できる業務は 88%](https://qiita.com/hisaho/items/c3f3ea9b2be822571ce8) だ。残る 12% ——研究の問いの設定、実験の意思決定、査読対応の戦略、研究倫理の判断——これらは **人間にしかできない**。

禅の「悟り」は、知識の蓄積ではない。**本質を見抜く力** だ。

SATORI は Agent に知識（88%）を与える。人間は本質（12%）を見抜く。この役割分担こそが「悟り」の構造であり、**Co-scientist の設計原理** だ。

# 逆引き辞書の概要

:::note info
**SATORI v0.26.0** — 2026-02-15 生成
:::

| 指標 | 件数 |
|------|------|
| スキル数 | **190** |
| ToolUniverse 連携スキル数 | **190** |
| ToolUniverse キー数（ユニーク） | **93** |
| 成果物パス数（ユニーク） | **648** |

# 使い方

:::note
**4 つの逆引きルート**

1. **ToolUniverse の key からスキルを探す** → [ToolUniverse key → スキル](#tooluniverse-key--スキル)
2. **手元の `results/*.csv` から生成元スキルを探す** → [成果物パス → スキル](#成果物パスresultsfigures-スキル)
3. **あるスキルがどのスキルから参照されているか** → [参照関係（逆引き）](#参照関係逆引き)
4. **トリガーフレーズからスキルを探す** → [トリガーフレーズ → スキル](#トリガーフレーズ説明文の-スキル)
:::

@[toc]

# スキル一覧（アルファベット順）

全 190 スキルをアルファベット順に掲載。各エントリの末尾に ToolUniverse 連携キーを記載（該当スキルのみ）。

<details><summary>A–C（37 スキル）</summary>

- `scientific-academic-writing` — 科学技術・学術論文の執筆スキル。IMRaD 標準、Nature/Science 系、ACS 系、IEEE 系、 Elsevier 系のジャーナル形式に対応した論文構成・セクション設計・文章パターンを提供。 「論文を書いて」「Abstr...
- `scientific-active-learning` — アクティブラーニング (能動学習) スキル。不確実性サンプリング・ Query-by-Committee・期待モデル変化・プール型/ストリーム型・ バッチアクティブラーニング・停止基準判定・ モデル改善パイプライン。
- `scientific-adaptive-experiments` — 適応的実験計画スキル。多腕バンディット (Thompson Sampling/UCB)・ ベイズ適応設計・逐次検定 (SPRT)・ Response-Adaptive Randomization・早期停止規則。
- `scientific-admet-pharmacokinetics` — ADMET 予測・薬物動態モデリングスキル。吸収(A)・分布(D)・代謝(M)・排泄(E)・毒性(T)の 包括的予測パイプライン。DeepChem/ADMET-AI/PyTDC を活用した分子特性予測、 PK/PD モデリング、ドラッ... （TU: `pubchem`）
- `scientific-advanced-imaging` — 高度バイオイメージング解析スキル。CellProfiler によるモフォロジカル プロファイリング・Cell Painting 解析、Cellpose による深層学習 セルセグメンテーション、napari によるインタラクティブ 3D...
- `scientific-advanced-visualization` — 科学データ高度可視化スキル。Plotly インタラクティブ 3D ・ Altair 宣言的可視化・Seaborn 統計プロット・ アニメーション・Parallel Coordinates・出版品質図。
- `scientific-alphafold-structures` — AlphaFold 構造予測スキル。AlphaFold Protein Structure Database REST API による予測構造取得・pLDDT 信頼度解析・PAE 残基間 距離予測・構造カバレッジ分析。ToolUniv... （TU: `alphafold`）
- `scientific-anomaly-detection` — 異常検知・外れ値検出スキル。Isolation Forest・LOF・ One-Class SVM・Autoencoder 異常検知・統計的工程管理 (SPC)・ 多変量異常検知・異常スコアリング・閾値最適化。
- `scientific-arrayexpress-expression` — ArrayExpress 発現アーカイブスキル。BioStudies/ArrayExpress REST API によるマイクロアレイ・RNA-seq 発現実験検索・メタ データ取得・データ再解析。ToolUniverse 連携: a... （TU: `arrayexpress`）
- `scientific-automl` — AutoML パイプラインスキル。Optuna ハイパーパラメータ最適化・ FLAML 高速 AutoML・Auto-sklearn モデル選択・ NAS (Neural Architecture Search)・ 特徴量エンジニアリ...
- `scientific-bayesian-statistics` — ベイズ統計スキル。PyMC・Stan・ArviZ を活用し、ベイズ回帰・階層モデル・ MCMC サンプリング・ベイズ最適化・事後予測チェック・モデル比較を支援。 「ベイズ回帰して」「MCMC で推定して」「事後分布を求めて」で発火。
- `scientific-biobank-cohort` — バイオバンク・大規模コホートデータ解析スキル。UK Biobank / BBJ / All of Us 等の大規模コホートデータに対するフェノタイプ 辞書検索・GWAS サマリー統計処理・PheWAS パイプライン。
- `scientific-bioinformatics` — バイオインフォマティクス解析パイプラインのスキル。scRNA-seq（Scanpy）、ゲノム配列解析 （BioPython）、PPI ネットワーク解析（NetworkX）、メタボロミクスの前処理を行う際に使用。 Scientific ...
- `scientific-biomedical-pubtator` — バイオメディカルテキストマイニングスキル。PubTator3 API による 遺伝子・疾患・化合物・変異・種のエンティティ認識、関係抽出、 バイオメディカル文献アノテーション自動化パイプライン。
- `scientific-biosignal-processing` — 生体信号処理スキル。ECG（R波検出・HRV時間/周波数ドメイン・Poincaréプロット）、 EEG（マルチチャネル・バンドパワーδ/θ/α/β/γ・スペクトログラム・ERP）、 EMG（バースト検出・包絡線）、呼吸信号（RSA）の...
- `scientific-biothings-idmapping` — BioThings API (MyGene.info, MyVariant.info, MyChem.info) を活用した 遺伝子・変異・化合物の横断的 ID マッピングおよびアノテーション統合スキル。 （TU: `biothings`）
- `scientific-cancer-genomics` — がんゲノミクスポータル統合スキル。COSMIC (体細胞変異カタログ)、 cBioPortal (がんゲノミクスデータ解析)、DepMap (がん細胞依存性) の 3 大がんゲノミクスデータベースを統合した変異プロファイリング、 変異... （TU: `cosmic`, `cbioportal`）
- `scientific-causal-inference` — 因果推論スキル。傾向スコアマッチング（PSM）、逆確率重み付け（IPW / IPTW）、 操作変数法（2SLS）、差分の差分法（DID）、回帰不連続デザイン（RDD）、 DAG ベースの共変量選択（backdoor criterion...
- `scientific-causal-ml` — 因果機械学習スキル。DoWhy 因果モデル・EconML CATE 推定・ Double/Debiased ML・Causal Forest・メタラーナー (S/T/X)・ 異質的処置効果 (HTE)・因果特徴量選択。
- `scientific-cell-line-resources` — 細胞株リソーススキル。Cellosaurus 細胞株データベース検索、 STR プロファイルマッチング、コンタミネーション検出、 細胞株メタデータ (由来組織・疾患・種) 取得パイプライン。 （TU: `cellosaurus`）
- `scientific-cellxgene-census` — CELLxGENE Census 大規模シングルセルアトラススキル。 CZ CELLxGENE Census API によるヒト/マウス全アトラスの メタデータ検索・遺伝子発現クエリ・セルタイプ分布解析・ データセット横断統合パイプラ... （TU: `cellxgene_census`）
- `scientific-chembl-assay-mining` — ChEMBL アッセイ・活性データマイニングスキル。ChEMBL REST API による アッセイ検索・バイオアクティビティデータ取得・IC50/Ki/EC50 SAR 解析・ ターゲット-化合物マッピング・選択性プロファイリング・... （TU: `chembl`）
- `scientific-cheminformatics` — ケモインフォマティクス解析のスキル。RDKit を用いた分子記述子計算、Morgan フィンガープリント、 Tanimoto 類似度、構造アラート検出、Lipinski Rule of 5 評価を行う際に使用。 Scientific ...
- `scientific-citation-checker` — 原稿中の引用文献の自動検索・網羅性チェックを行うスキル。 参照リスト抽出、DOI/タイトルベース自動検索、引用カバレッジ分析、 フォーマット一貫性検証、重複検出を実行する。 「引用をチェックして」「参考文献を検索」「citation ...
- `scientific-civic-evidence` — CIViC 臨床エビデンススキル。CIViC (Clinical Interpretation of Variants in Cancer) REST API を用いたバリアント臨床解釈・ エビデンスアイテム・分子プロファイル・アサー... （TU: `civic`）
- `scientific-clingen-curation` — ClinGen 臨床ゲノム資源キュレーションスキル。ClinGen API に よる遺伝子-疾患バリディティ、臨床アクショナビリティ、 投与量感受性、バリアントレベルエビデンス評価パイプライン。 ToolUniverse 連携: cl... （TU: `clingen`）
- `scientific-clinical-decision-support` — 臨床意思決定支援スキル。エビデンスに基づく治療推奨、臨床パスウェイアルゴリズム、 患者コホート解析、バイオマーカー分類、精密腫瘍学（ToolUniverse の Precision Oncology パラダイム + claude-sc...
- `scientific-clinical-nlp` — 臨床自然言語処理スキル。MedSpaCy / cTAKES / scispaCy による臨床テキスト NER、セクション検出、否定文検出、 ICD-10/SNOMED-CT エンティティリンキング、 匿名化 (De-identific...
- `scientific-clinical-pharmacology` — 臨床薬理学モデリングスキル。PopPK (NLME 混合効果モデル)・ PBPK シミュレーション・TDM 投与量最適化・ Emax/Sigmoid PD モデリング・薬物間相互作用予測・ 臨床薬理パイプライン。 TU 外スキル (P...
- `scientific-clinical-reporting` — 臨床レポート自動生成スキル。検査結果サマリー (SOAP ノート)、バイオマーカー プロファイルレポート、薬理ゲノミクスレポート、臨床試験要約を構造化テンプレート (PDF/LaTeX/HTML) で出力。HL7 FHIR Diagn...
- `scientific-clinical-standards` — 臨床標準用語・コードマッピングスキル。LOINC 臨床検査コード・ ICD-10/ICD-11 疾病分類・FHIR R4 リソースマッピング・ SNOMED CT 用語変換・臨床用語相互運用パイプライン。 （TU: `loinc`, `icd`）
- `scientific-clinical-trials-analytics` — 臨床試験レジストリ解析スキル。ClinicalTrials.gov API v2 経由の 多基準試験検索 (疾患×介入×地域×フェーズ×ステータス)、試験詳細取得 (適格基準・アウトカム・施設・有害事象)、競合ランドスケープ解析、 A...
- `scientific-compound-screening` — 化合物スクリーニングスキル。ZINC データベースを活用した購入可能化合物検索、 SMILES/名前ベースの類似性検索、カタログフィルタリング、 バーチャルスクリーニング前処理パイプライン。 （TU: `zinc`）
- `scientific-computational-materials` — 計算材料科学スキル。pymatgen による結晶構造操作・対称性解析、 Materials Project API による材料データベース照会、 相図計算 (凸包解析)、電子バンド構造・状態密度 (DOS) 可視化、 VASP/Qua...
- `scientific-crispr-design` — CRISPR gRNA 設計スキル。Cas9/Cas12a PAM 配列検索・ オフターゲットスコアリング (CFD/MIT)・ CRISPRscan/Rule Set 2 活性予測・検証プライマー設計・ sgRNA スクリーニングラ...
- `scientific-critical-review` — 学術論文の草稿に対する批判的レビュー・修正スキル。論理構成、考察の深さ、 データ解釈の妥当性、先行研究との整合性、統計的主張の正確性を多角的に検証し、 具体的な修正案を生成する。「論文をレビューして」「考察を深めて」「草稿を改善して」...
- `scientific-crossref-metadata` — CrossRef メタデータスキル。CrossRef REST API による DOI 解決・論文メタデータ・引用数・ジャーナル情報・ 助成金情報検索。ToolUniverse 連携: crossref。 （TU: `crossref`）

</details>

<details><summary>D–F（26 スキル）</summary>

- `scientific-data-preprocessing` — 科学データの前処理パイプラインスキル。欠損値補完（KNNImputer/SimpleImputer）、 エンコーディング（LabelEncoder/OneHot/ダミー変数）、スケーリング（Standard/MinMax/Robust...
- `scientific-data-profiling` — データプロファイリング・品質スキル。ydata-profiling 自動 EDA ・ Great Expectations データバリデーション・データ品質スコア・ 型推論・相関検出・外れ値フラグ・データカタログ生成。
- `scientific-data-simulation` — 物理・化学・生物学に基づく合成データ生成のスキル。実験データが不足する場合に、 ドメイン知識を反映したシミュレーションデータを生成する際に使用。 Scientific Skills Exp-06, 07, 08, 09, 12, 13...
- `scientific-data-submission` — 科学データ登録・アーカイブスキル。GenBank/SRA 配列登録・ ENA 配列アーカイブ・GEO 発現データ登録・BioProject/BioSample メタデータ管理・FAIR 原則準拠データ共有。
- `scientific-deep-chemistry` — 深層学習分子特性予測スキル。DeepChem による GCN/MPNN/AttentiveFP 分子特性予測・MoleculeNet ベンチマーク・ChemBERTa/GROVER 事前学習モデル・分子フィンガープリントフィーチャライザ。
- `scientific-deep-learning` — 深層学習スキル。PyTorch Lightning・Hugging Face Transformers・timm を活用し、 NN アーキテクチャ設計・転移学習・分散トレーニング・ハイパーパラメータ最適化・ モデルデプロイを支援。 「...
- `scientific-deep-research` — 科学文献の深層リサーチスキル。SHIKIGAMI の WebResearcher パラダイム （Think→Report→Action 反復サイクル）を科学研究に特化させた実装。 学術データベース検索、エビデンス階層評価、ソース追跡、...
- `scientific-depmap-dependencies` — DepMap 依存性スキル。Cancer Dependency Map (DepMap) Portal API によるがん細胞株 CRISPR/RNAi 依存性スコア・薬剤 感受性データ・遺伝子効果取得。ToolUniverse 連携... （TU: `depmap`）
- `scientific-disease-research` — 疾患研究スキル。GWAS Catalog・Orphanet・OMIM・HPO・DisGeNET を統合し、 疾患-遺伝子関連解析・希少疾患診断支援・表現型-遺伝型マッピング・ 疫学的特性評価を支援。 「疾患と遺伝子の関連を調べて」「希... （TU: `disgenet`）
- `scientific-doe` — 実験計画法（DOE）スキル。直交配列表（L9/L16/L27）、中心複合計画（CCD）、 Box-Behnken 設計、D-最適計画、応答曲面法（RSM）、交互作用解析、 ベイズ最適化（Gaussian Process）、効果プロット...
- `scientific-drug-repurposing` — ドラッグリポジショニング（既存薬再創薬）スキル。ToolUniverse の Drug Repurposing パラダイムに準拠し、7 つの戦略（ターゲット型、化合物型、疾患駆動型、メカニズム型、 ネットワーク型、表現型、構造型）で候... （TU: `pharos`）
- `scientific-drug-target-profiling` — 創薬ターゲットプロファイリングスキル。ToolUniverse / Open Targets / ChEMBL / UniProt を活用したドラッグターゲットインテリジェンス。ドラッガビリティ評価、安全性プロファイリング、 ターゲッ... （TU: `dgidb`）
- `scientific-drugbank-resources` — DrugBank リソーススキル。DrugBank API を用いた薬剤記述・ 薬理情報・標的タンパク質・薬物相互作用検索。 ToolUniverse 連携: drugbank。 （TU: `drugbank`）
- `scientific-ebi-databases` — EBI データベース群統合アクセススキル。EBI Search 横断検索、ENA Browser ヌクレオチドアーカイブ、BioStudies 研究データ、dbfetch エントリ取得、 MetaboLights メタボロミクスリポジ...
- `scientific-eda-correlation` — 探索的データ解析（EDA）と相関分析のスキル。データの分布可視化、相関ヒートマップ、 散布図行列の作成を行う際に使用。Scientific Skills Exp-02, 12, 13 で確立したパターン。
- `scientific-encode-screen` — ENCODE / ChIP-Atlas エピゲノムアトラススキル。ENCODE REST API 実験・ファイル・バイオサンプル検索、SCREEN cis 制御エレメント、 ChIP-Atlas エンリッチメント解析、エピゲノムアノテ... （TU: `encode`, `chipatlas`）
- `scientific-ensembl-genomics` — Ensembl REST API ゲノミクススキル。遺伝子ルックアップ・配列取得・ VEP (Variant Effect Predictor) バリアントアノテーション・ クロスリファレンス・制御要素・系統樹・相同性検索・分類学統合...
- `scientific-ensemble-methods` — アンサンブル学習スキル。Stacking/Blending 多段積層・ Boosting (XGBoost/LightGBM/CatBoost) 勾配ブースティング・ Bagging/Random Subspace・Voting 分類...
- `scientific-environmental-ecology` — 環境科学・生態学解析スキル。種分布モデリング（SDM / MaxEnt）・ 生物多様性指標（α/β/γ 多様性）・群集構造解析（NMDS/CCA/RDA）・ 生態学的ニッチモデリング・保全優先順位評価・OBIS/GBIF データ統合パ... （TU: `gbif`）
- `scientific-environmental-geodata` — 環境地理空間データスキル。SoilGrids REST API による土壌特性 取得、WorldClim/CHELSA 気候データ、生物多様性-環境モデリング 統合。直接 REST API 連携 (TU 外)。
- `scientific-epidemiology-public-health` — 疫学・公衆衛生解析スキル。観察研究デザイン（コホート/症例対照/横断）・ リスク指標（RR/OR/HR/NNT）・標準化死亡比（SMR）・年齢調整率・ 空間疫学（GIS / 空間クラスタリング）・因果推論ダイアグラム（DAG）・ WH... （TU: `who_gho`）
- `scientific-epigenomics-chromatin` — エピゲノミクス・クロマチン生物学解析スキル。ChIP-seq ピーク呼び出し (MACS2/MACS3)、 ATAC-seq ヌクレオソームフリー領域検出、DNA メチル化パターン解析 (WGBS/RRBS)、 ヒストン修飾クロマチン... （TU: `chipatlas`）
- `scientific-explainable-ai` — 説明可能 AI (XAI) スキル。SHAP・LIME・Captum・InterpretML を活用し、 モデル予測の根拠説明・特徴量寄与分解・反実仮想説明・公平性監査を支援。 「モデルの予測を説明して」「SHAP 値を計算して」「L...
- `scientific-expression-comparison` — Expression Atlas / GTEx / HPA 統合発現比較スキル。EBI Expression Atlas ベースライン/差次的発現検索、実験アクセション取得、組織間・条件間 発現比較、マルチソース統合発現プロファイリン...
- `scientific-feature-importance` — 特徴量重要度分析のスキル。Tree-based Feature Importance と Permutation Importance を 用いて予測モデルの説明可能性を向上させる際に使用。 Scientific Skills Exp...
- `scientific-federated-learning` — 連合学習スキル。Flower フレームワークによる FL パイプライン・ FedAvg/FedProx/FedOpt 集約戦略・差分プライバシー (DP-SGD)・ 非 IID データ分割・通信効率化。

</details>

<details><summary>G–I（22 スキル）</summary>

- `scientific-gdc-portal` — NCI Genomic Data Commons ポータルスキル。GDC REST API を用いたがんゲノムプロジェクト横断検索・ケースメタデータ・ 体細胞変異 (SSM)・遺伝子発現・ファイル取得。 ToolUniverse 連携... （TU: `gdc`）
- `scientific-gene-expression-transcriptomics` — 遺伝子発現・トランスクリプトミクス解析スキル。GEO (Gene Expression Omnibus) からの 公開データセット取得・前処理、DESeq2 (PyDESeq2) による差次発現解析、 GTEx 組織発現参照・eQTL...
- `scientific-genome-sequence-tools` — ゲノム配列解析総合スキル。Ensembl ゲノムブラウザ、dbSNP 変異データ、 BLAST 相同性検索、NCBI Nucleotide 配列取得、GDC がんゲノミクスデータの 統合パイプライン。
- `scientific-geo-expression` — GEO (Gene Expression Omnibus) 発現プロファイルスキル。GEO REST API データセット検索・サンプル情報・発現マトリクス取得・バルク RNA-seq/マイクロアレイ差次的発現解析。ToolUnive... （TU: `geo`）
- `scientific-geospatial-analysis` — 地理空間データ解析スキル。GeoPandas ベクターデータ処理・ Rasterio ラスター解析・Folium/Kepler.gl インタラクティブ地図・ 空間自己相関 (Moran's I)・クリギング補間・CRS 変換。
- `scientific-glycomics` — 糖鎖構造解析スキル。GlyConnect / GlyGen / GlyCosmos 糖鎖データベース統合検索・糖鎖構造描画・糖タンパク質 グリコシル化部位予測・レクチンバインディング・ 糖鎖マスフラグメンテーション解析パイプライン。 ...
- `scientific-gnomad-variants` — gnomAD バリアントスキル。gnomAD (Genome Aggregation Database) GraphQL API を用いた集団アレル頻度・遺伝子制約スコア (pLI/LOEUF)・リージョンクエリ・トランスクリプトレベ... （TU: `gnomad`）
- `scientific-gpu-singlecell` — GPU アクセラレーション シングルセル解析スキル。 rapids-singlecell / cuML / cuGraph による GPU 並列処理。 大規模 (>1M cells) データの高速前処理・クラスタリング・ 次元削減。K...
- `scientific-grant-writing` — 研究グラント（助成金申請書）執筆スキル。NIH R01/R21、NSF、JSPS 科研費、 ERC 等のフォーマットに対応。Specific Aims、Research Strategy、予算計画、 Biosketch の構造化作成を...
- `scientific-graph-neural-networks` — グラフニューラルネットワーク (GNN) スキル。PyTorch Geometric・TorchDrug・ DeepChem を活用し、分子特性予測・タンパク質モデリング・知識グラフ推論・ ノード/グラフ分類を支援。 「GNN で分子...
- `scientific-gtex-tissue-expression` — GTEx 組織発現スキル。GTEx Portal REST API v2 による 組織特異的遺伝子発現パターン解析・eQTL ルックアップ・ 多組織比較。ToolUniverse 連携: gtex_v2。 （TU: `gtex_v2`）
- `scientific-gwas-catalog` — GWAS カタログスキル。NHGRI-EBI GWAS Catalog REST API によるゲノム ワイド関連研究メタデータ・関連シグナル・形質・遺伝子座検索。 ToolUniverse 連携: gwas。 （TU: `gwas`）
- `scientific-healthcare-ai` — ヘルスケア AI スキル。PyHealth 臨床 ML パイプライン、 フローサイトメトリー (FlowIO) 解析、電子健康記録 (EHR) 処理、 臨床予測モデル構築のガイダンス。
- `scientific-hgnc-nomenclature` — HGNC 遺伝子命名法スキル。HUGO Gene Nomenclature Committee REST API による公式遺伝子シンボル検索・エイリアス解決・ 遺伝子ファミリー/グループクエリ・ID クロスリファレンス パイプライン...
- `scientific-human-cell-atlas` — Human Cell Atlas (HCA) データポータルスキル。HCA Data Portal API プロジェクト検索・ファイルダウンロード・CELLxGENE Census 統合・ 細胞型アノテーション・アトラス構築。Tool... （TU: `hca_tools`, `cellxgene_census`）
- `scientific-human-protein-atlas` — Human Protein Atlas (HPA) 統合スキル。組織/細胞タンパク質発現、 がん予後バイオマーカー、RNA 発現プロファイル、細胞内局在、 タンパク質相互作用の包括的検索・解析パイプライン。 （TU: `hpa`）
- `scientific-hypothesis-pipeline` — ユーザーのプロンプト（研究テーマ・データ記述）から仮説を立案し、 検証用の解析パイプラインを自動生成するスキル。PICO/PECO フレームワークによる 仮説構造化、適切な統計検定の選択、パイプラインコード生成を行う。 「仮説を立てて...
- `scientific-icgc-cancer-data` — ICGC がんゲノムデータスキル。ICGC ARGO DCC API および レガシー API による国際がんゲノムデータ検索・ドナー/ 検体/変異解析。直接 API (ToolUniverse 非連携)。
- `scientific-image-analysis` — 科学画像解析スキル。顕微鏡画像のセグメンテーション（Otsu/Watershed/Felzenszwalb）、 粒径分布解析、形態計測（面積・周囲長・真円度・アスペクト比）、テクスチャ解析 （GLCM/LBP）、強度プロファイル、マル...
- `scientific-immunoinformatics` — 免疫情報学スキル。エピトープ予測（MHC-I/II バインディング）・ T 細胞/B 細胞エピトープマッピング・抗体構造解析（CDR ループ）・ 免疫レパトア解析（TCR/BCR クロノタイプ）・ワクチン候補設計・ IEDB/IMGT... （TU: `iedb`, `imgt`, `sabdab`, `therasabdab`）
- `scientific-infectious-disease` — 感染症ゲノミクス・疫学スキル。病原体ゲノム解析（SNP/系統樹）・ AMR（薬剤耐性）遺伝子検出・分子疫学（MLST/cgMLST）・ アウトブレイク調査トレーシング・疫学的 SIR/SEIR コンパートメントモデル・ 伝播ネットワー...
- `scientific-interactive-dashboard` — インタラクティブダッシュボードスキル。 Streamlit / Dash / Panel / Voilà による 科学データダッシュボード構築・リアルタイムパラメータ探索 UI ・ ウィジェット連動・データアップロード・解析パイプライ...

</details>

<details><summary>L–M（27 スキル）</summary>

- `scientific-lab-automation` — 実験室自動化・プロトコル管理スキル。PyLabRobot（液体ハンドリング）、 Protocols.io（プロトコル共有）、Benchling/LabArchives（ELN/LIMS 統合）、 Opentrons（ロボティクス）によ...
- `scientific-lab-data-management` — ラボデータ管理スキル。Benchling (ELN/DNA 設計/レジストリ)、 DNAnexus (ゲノミクス PaaS)、LatchBio (ワークフロー)、 OMERO (バイオイメージング)、Protocols.io (プロト...
- `scientific-latex-formatter` — Markdown 原稿を LaTeX 形式に変換し、ジャーナル指定のテンプレート（.cls/.sty）に 適合するフォーマッティングを行うスキル。数式・図表・引用・相互参照の LaTeX 構文変換、 ジャーナル別スタイル適用、コンパイ...
- `scientific-lipidomics` — リピドミクス解析スキル。LipidMAPS / SwissLipids / LION 脂質データベース統合検索・脂質サブクラス分類・ 脂質 MS/MS スペクトル同定・脂質パスウェイエンリッチメント・ 脂質プロファイリングパイプライン...
- `scientific-literature-search` — 学術文献検索・取得スキル。PubMed E-utilities、Semantic Scholar、 OpenAlex、EuropePMC、CrossRef の 5 大学術データベース API を統合した 文献検索パイプライン。MeSH...
- `scientific-marine-ecology` — 海洋生態学統合スキル。OBIS 海洋生物分布・WoRMS 海洋分類体系・ GBIF 生物多様性レコード・FishBase 魚類データ。ToolUniverse 連携: obis, worms, gbif。 （TU: `obis`, `worms`, `gbif`）
- `scientific-materials-characterization` — 薄膜・材料キャラクタリゼーション解析のスキル。Thornton-Anders 構造ゾーンモデル（SZM）、 XRD 結晶子サイズ解析（Scherrer 方程式）、Williamson-Hall プロット、多技法融合データ解析、 PSP...
- `scientific-md-simulation` — 分子動力学シミュレーション解析スキル。MDAnalysis によるトラジェクトリ解析・ RMSD/RMSF/Rg 時系列指標・水素結合解析・二次構造変化追跡・ OpenFF Toolkit力場パラメータ化・溶媒和自由エネルギー推定パイ...
- `scientific-medical-imaging` — 医用イメージングスキル。DICOM/NIfTI 処理・WSI (Whole Slide Image) 解析・ PathML・MONAI・3D Slicer 連携・放射線画像解析・病理組織画像解析を支援。 「DICOM を解析して」「W...
- `scientific-meta-analysis` — メタ解析スキル。固定効果・ランダム効果モデル（DerSimonian-Laird）、Forest プロット、 異質性評価（I²/Q 検定/τ²）、出版バイアス検出（Funnel プロット/Egger/Begg 検定）、 サブグループ解...
- `scientific-metabolic-atlas` — 代謝アトラススキル。Metabolic Atlas / Human-GEM REST API による 代謝反応・代謝産物・コンパートメント検索、フラックス解析統合、 代謝ネットワーク可視化。K-Dense 連携: metabolic-...
- `scientific-metabolic-flux` — 代謝フラックス解析スキル。13C/15N 安定同位体トレーサー データを用いた代謝フラックス推定・EMU モデリング・ フラックスバランス制約統合パイプライン。
- `scientific-metabolic-modeling` — 代謝モデリングスキル。BiGG Models ゲノムスケール代謝モデル、 BioModels SBML リポジトリを統合した代謝ネットワーク解析・ モデル検索パイプライン。 （TU: `biomodels`）
- `scientific-metabolomics` — メタボロミクス解析スキル。Pareto スケーリング、PLS-DA + VIP スコア、置換検定（Q²）、 代謝パスウェイ濃縮解析（Fisher exact test）、代謝物相関ネットワーク、 Volcano プロット/箱ひげ図によ... （TU: `hmdb`, `metabolomics_workbench`）
- `scientific-metabolomics-databases` — メタボロミクスデータベース統合スキル。HMDB (Human Metabolome Database、 220,000+ 代謝物)、MetaCyc (代謝パスウェイ)、Metabolomics Workbench (NIH メタボロミ... （TU: `metacyc`）
- `scientific-metabolomics-network` — 代謝物ネットワーク構築スキル。KEGG/Reactome 代謝パスウェイ グラフ抽出・代謝物相関ネットワーク構築 (GGM/WGCNA)・ ハブ代謝物同定・MetaboAnalyst 統合エンリッチメント パイプライン。 TU 外スキ...
- `scientific-metagenome-assembled-genomes` — メタゲノムアセンブルゲノム (MAG) 解析スキル。 MetaBAT2 / CONCOCT / MaxBin2 ビニング・CheckM2 品質評価・ GTDB-Tk 分類学的分類・dRep 脱重複・Prokka アノテーション・ MA...
- `scientific-microbiome-metagenomics` — マイクロバイオーム・メタゲノミクス解析スキル。16S rRNA アンプリコン解析（DADA2）・ ショットガンメタゲノム解析（MetaPhlAn / HUMAnN）・α/β 多様性・ 差次存在量解析（DESeq2 / ANCOM-BC... （TU: `mgnify`）
- `scientific-missing-data-analysis` — 欠損データ解析スキル。欠損パターン診断 (MCAR/MAR/MNAR) ・ Little's MCAR テスト・多重代入法 (MICE) ・KNN 補完・ MissForest・VAE/GAIN 補完・欠損パターン可視化・Rubin'...
- `scientific-ml-classification` — 機械学習分類パイプラインのスキル。複数の分類モデル（Logistic Regression, Random Forest, SVM, XGBoost）を StratifiedKFold 交差検証で比較し、ROC 曲線・混同行列で評価す...
- `scientific-ml-regression` — マルチターゲット回帰モデルの学習・評価・比較スキル。複数の回帰モデル（Ridge, Lasso, Random Forest, Gradient Boosting, Extra Trees）を KFold 交差検証で比較する際に使用。...
- `scientific-model-monitoring` — MLOps モデル監視スキル。データドリフト検出 (Evidently/NannyML)・ モデル性能劣化検出・特徴量ドリフト・コンセプトドリフト・ A/B テスト統計・モデルレジストリ・再学習トリガー。
- `scientific-model-organism-db` — モデル生物データベース統合スキル。FlyBase (ショウジョウバエ)、 WormBase (線虫)、ZFIN (ゼブラフィッシュ)、RGD (ラット)、 MGI (マウス) の REST API を統合した モデル生物遺伝子・表現型... （TU: `impc`, `mpd`）
- `scientific-molecular-docking` — 構造ベース分子ドッキングスキル。DiffDock (拡散生成モデル)、 AutoDock Vina (スコアリング関数)、GNINA (CNN ベーススコアリング) を統合した タンパク質-リガンド結合ポーズ予測、バーチャルスクリーニ...
- `scientific-monarch-ontology` — Monarch Initiative 疾患-表現型オントロジースキル。 Monarch Initiative API を用いた疾患-遺伝子-表現型 アソシエーション・HPO フェノタイピング・ 遺伝子-疾患推定・オントロジーセマンティ... （TU: `monarch`）
- `scientific-multi-omics` — マルチオミクス統合解析スキル。ゲノム・トランスクリプトーム・プロテオーム・メタボローム データの統合手法（MOFA/SNF/DIABLO）、オミクス間相関解析、CCA/PLS 統合、 パスウェイレベル統合、ネットワーク統合のテンプレー...
- `scientific-multi-task-learning` — マルチタスク学習スキル。Hard/Soft Parameter Sharing・ GradNorm 勾配正規化・PCGrad 勾配投影・ タスクバランシング・補助タスク設計。

</details>

<details><summary>N–P（35 スキル）</summary>

- `scientific-nci60-screening` — NCI-60 がん細胞株薬剤応答スキル。CellMiner API 薬剤感受性・ NCI-60 GI50/LC50 データ・DepMap cancer dependency 統合・ 薬剤-分子マーカー相関・細胞株パネル比較解析。
- `scientific-network-analysis` — ネットワーク解析・相関ネットワーク構築のスキル。NetworkX を用いたグラフ構築、 中心性解析、コミュニティ検出、ネットワーク可視化を行う際に使用。 Scientific Skills Exp-04, 07 で確立したパターン。P...
- `scientific-network-visualization` — ネットワーク解析・可視化スキル。NetworkX グラフ構築・ コミュニティ検出 (Louvain/Leiden)・中心性指標・ PyVis インタラクティブ・ネットワーク統計量・動的ネットワーク。
- `scientific-neural-architecture-search` — ニューラルアーキテクチャ探索 (NAS) スキル。DARTS 微分可能 NAS・ Optuna NAS 統合・効率的ネットワーク設計・探索空間定義・ Pareto 最適化 (精度 vs FLOPS)。
- `scientific-neuroscience-electrophysiology` — 神経科学・電気生理学解析スキル。Neuropixels/マルチ電極アレイの スパイクソーティング (SpikeInterface + Kilosort4)、品質指標 (SNR/ISI 違反率)、 EEG マイクロステート・事象関連電位...
- `scientific-noncoding-rna` — 非コード RNA (ncRNA) 解析スキル。Rfam RNA ファミリー検索、 RNAcentral 統合 ncRNA データベース、共分散モデル、構造マッピング、 系統樹解析パイプライン。
- `scientific-ontology-enrichment` — オントロジー・エンリッチメント解析スキル。EFO 実験ファクターオントロジー、 OLS オントロジー検索サービス、Enrichr 遺伝子セット濃縮解析、 UMLS メタシソーラス統一医学言語体系の統合パイプライン。
- `scientific-opentargets-genetics` — Open Targets Platform 遺伝学スキル。Open Targets Platform GraphQL API を用いた標的-疾患アソシエーション・薬剤 エビデンス・L2G 遺伝的関連・ファーマコゲノミクス検索。 Too... （TU: `opentarget`）
- `scientific-paleobiology` — 古生物学データベーススキル。Paleobiology Database (PBDB) REST API による化石産出記録・分類群・コレクション検索、地質年代 多様性曲線・古地理解析。ToolUniverse 連携: paleobio... （TU: `paleobiology`）
- `scientific-paper-quality` — 論文品質の定量的評価スキル。可読性スコア、セクションバランス分析、 語彙多様性、学術語使用率、冗長表現検出、ジャーナル要件適合チェック、 再現可能性チェックを実行する。 「論文の品質をチェックして」「可読性スコアを出して」「投稿前チェ...
- `scientific-parasite-genomics` — 寄生虫ゲノミクススキル。PlasmoDB/VectorBase/ToxoDB REST API による寄生虫ゲノム検索・遺伝子情報・薬剤標的同定・比較 ゲノミクス。直接 REST API 連携 (TU 外)。
- `scientific-pathway-enrichment` — パスウェイ・Gene Ontology 富化解析スキル。KEGG パスウェイ検索・マッピング、 Reactome パスウェイ階層解析、Gene Ontology (BP/MF/CC) アノテーション、 WikiPathways コミュ...
- `scientific-pca-tsne` — PCA・t-SNE・UMAP による次元削減と空間マッピングのスキル。化学空間・特徴量空間・ 多技法融合空間の可視化を行う際に使用。Scientific Skills Exp-02, 03, 05, 07, 11, 13 で汎用的に使...
- `scientific-peer-review-response` — 査読コメントへの体系的対応とリバッタルレター生成スキル。 査読コメントの構造化解析（Major/Minor/Editorial 分類）、 ポイント・バイ・ポイント回答生成、改訂箇所マッピング、 複数ラウンド対応を行う。 「査読に回答し...
- `scientific-perturbation-analysis` — シングルセル摂動解析スキル。pertpy による CRISPR スクリーン解析・ 薬剤応答分析・scGen 摂動予測・Augur 摂動応答性スコアリング・ scIB 統合ベンチマーク・差次的摂動応答パイプライン。
- `scientific-pharmacogenomics` — ファーマコゲノミクス (薬理ゲノム学) 解析スキル。PharmGKB/ClinPGx による 遺伝子-薬物相互作用照会、CPIC ガイドライン取得・解釈、Star アレル分類、 代謝酵素表現型判定 (PM/IM/NM/RM/UM)、F... （TU: `fda_pharmacogenomic_biomarkers`）
- `scientific-pharmacology-targets` — 薬理学的ターゲットプロファイリングスキル。BindingDB 結合親和性、 GPCRdb GPCR 構造-活性、GtoPdb 薬理学、BRENDA 酵素動態、 Pharos 未解明ターゲット(TDL)の統合解析パイプライン。 （TU: `bindingdb`, `gtopdb`, `brenda`）
- `scientific-pharmacovigilance` — ファーマコビジランス（医薬品安全性監視）スキル。FAERS/FDA 有害事象報告データベースを活用し、 不均衡分析（PRR/ROR/IC）、MedDRA 階層構造、時系列トレンド、人口統計層別化を実施。 市販後安全性シグナル検出と定量...
- `scientific-pharmgkb-pgx` — PharmGKB 薬理ゲノミクススキル。PharmGKB REST API による 臨床アノテーション・薬物遺伝子関連・投与量ガイドライン・ スターアレル解析。ToolUniverse 連携: pharmgkb。 （TU: `pharmgkb`）
- `scientific-pharos-targets` — Pharos/TCRD ターゲットプロファイリングスキル。Illuminating the Druggable Genome (IDG) Pharos GraphQL API による ターゲット開発レベル (TDL) 分類・疾患関連・... （TU: `pharos`）
- `scientific-phylogenetics` — 系統解析スキル。ete3/ETE Toolkit による系統樹構築・可視化、 scikit-bio 系統的多様性、配列アライメントベース進化解析、 分子時計・分岐年代推定、祖先配列再構成パイプライン。
- `scientific-pipeline-scaffold` — 科学データ解析パイプラインの基盤スキル。ディレクトリ構造の自動構築、再現性のためのシード管理、 進捗ログ出力、実行時間計測、JSON サマリー生成、ダッシュボード総括図の作成を行う際に使用。 全 13 実験に共通する足場パターンを統合。
- `scientific-plant-biology` — 植物バイオロジー統合スキル。Plant Reactome 代謝パスウェイ・ TAIR Arabidopsis ゲノム情報・Phytozome 比較ゲノミクス・ Ensembl Plants 種間オーソログ解析。
- `scientific-population-genetics` — 集団遺伝学解析スキル。アレル頻度解析・Hardy-Weinberg 平衡検定・ 集団構造解析（PCA / ADMIXTURE）・Fst 分化指標・選択圧検出（iHS / XP-EHH）・ 連鎖不平衡（LD）解析・GWAS Catalo...
- `scientific-precision-oncology` — 精密腫瘍学スキル。CIViC・OncoKB・cBioPortal・COSMIC・GDC/TCGA を統合し、 腫瘍ゲノムプロファイリング・分子標的選定・バイオマーカー評価・治療推奨を支援。 「がんゲノム解析して」「腫瘍プロファイリング... （TU: `oncokb`）
- `scientific-preprint-archive` — プレプリント・オープンアクセスアーカイブ検索スキル。bioRxiv/medRxiv プレプリント検索、arXiv 論文取得、PMC フルテキスト、DOAJ OA ジャーナル、 Unpaywall OA リンク、CORE/HAL/Zen...
- `scientific-presentation-design` — 科学プレゼンテーション・ポスター・模式図設計スキル。学会発表スライド、 LaTeX/PPTX ポスター、科学模式図（ワークフロー図・メカニズム図）、 ビジュアルアブストラクトの作成を支援。claude-scientific-skill...
- `scientific-process-optimization` — 応答曲面法（ML-RSM）とパレート多目的最適化のスキル。プロセスパラメータの最適条件探索、 コンターマップ、プロセスウィンドウ可視化を行う際に使用。 Scientific Skills Exp-12, 13 で確立したパターン。
- `scientific-protein-design` — タンパク質設計スキル。ESM タンパク質言語モデル、de novo 設計、指向性進化の 計算的ガイド、安定性予測をカバー。ToolUniverse の Protein Therapeutic Design パラダイムと claude-...
- `scientific-protein-domain-family` — タンパク質ドメイン・ファミリー解析スキル。InterPro アノテーション検索、 InterProScan によるシーケンスベースドメイン予測、Pfam/SMART/CDD ドメイン分類、ドメインアーキテクチャ可視化、ファミリー系統樹構築。
- `scientific-protein-interaction-network` — タンパク質-タンパク質相互作用 (PPI) ネットワーク解析スキル。STRING、IntAct、 BioGRID、STITCH (化学-タンパク質) 相互作用データベースを統合した ネットワーク構築・解析パイプライン。GO/KEGG ... （TU: `intact`）
- `scientific-protein-structure-analysis` — タンパク質構造解析スキル。PDB / AlphaFold DB / PDBe を活用した 3D 構造解析。 構造アラインメント、結合部位検出、分子ドッキング準備、構造品質評価。 ToolUniverse の Protein Struc... （TU: `proteinsplus`）
- `scientific-proteomics-mass-spectrometry` — プロテオミクス・質量分析解析スキル。LC-MS/MS データ前処理、ペプチド同定 (PSM/FDR 制御)、 蛋白質定量 (LFQ/TMT/SILAC/iBAQ)、翻訳後修飾 (PTM) マッピング、 スペクトル類似度スコアリング (...
- `scientific-public-health-data` — 公衆衛生データアクセススキル。NHANES 疫学調査データ、MedlinePlus 一般向け 健康情報、RxNorm 薬剤標準語彙、ODPHP 健康目標・ガイドライン、 Health Disparities 健康格差データ統合パイプラ... （TU: `nhanes`, `medlineplus`, `odphp`）
- `scientific-publication-figures` — 論文品質（Nature/Science/Cell レベル）の科学図表を作成するスキル。matplotlib rcParams 設定、 DPI 300、spines 制御、カラーパレット選択、マルチパネル構成を行う際に使用。 Scien...

</details>

<details><summary>Q–V（43 スキル）</summary>

- `scientific-quantum-computing` — 量子計算スキル。Qiskit・Cirq・PennyLane・QuTiP を活用し、 量子回路設計・シミュレーション・変分量子アルゴリズム（VQE/QAOA）・ 量子化学計算・量子機械学習を支援。 「量子回路を設計して」「VQE で基底...
- `scientific-radiology-ai` — 放射線診断支援 AI スキル。CADe/CADx パイプライン・ CT/MRI 分類・セグメンテーション・Grad-CAM 説明可能性・ 構造化レポート・AI-RADS グレーディング。 ※ scientific-medical-im...
- `scientific-rare-disease-genetics` — 希少疾患遺伝学スキル。OMIM 遺伝子-疾患マッピング、Orphanet 希少疾患 分類・遺伝子照会、DisGeNET 疾患-遺伝子関連スコア、IMPC マウス表現型 参照、遺伝子-表現型統合解析パイプライン。 （TU: `orphanet`）
- `scientific-rcsb-pdb-search` — RCSB PDB 構造検索スキル。RCSB PDB Search API および Data API によるタンパク質立体構造検索・メタデータ取得・ リガンド情報・解像度フィルタリング。ToolUniverse 連携: rcsb_pdb... （TU: `rcsb_pdb`, `rcsb_search`）
- `scientific-reactome-pathways` — Reactome パスウェイスキル。Reactome Content Service REST API によるパスウェイ検索・階層取得・UniProt マッピング・ パスウェイ図データ取得。ToolUniverse 連携: react... （TU: `reactome`）
- `scientific-regulatory-genomics` — レギュラトリーゲノミクススキル。RegulomeDB バリアント制御機能スコア、 ReMap 転写因子結合マッピング、4D Nucleome (4DN) 三次元ゲノム構造 解析の統合パイプライン。
- `scientific-regulatory-science` — 規制科学パイプラインスキル。FDA (医薬品/医療機器/食品)・EMA・PMDA 規制データベース横断照会、 Orange Book 承認履歴・特許・排他性情報、510(k) デバイスクリアランス、 ISO 13485 品質管理システ...
- `scientific-reinforcement-learning` — 強化学習スキル。Stable-Baselines3 による RL エージェント訓練、 Gymnasium 環境構築、PufferLib 大規模マルチエージェント、 科学応用 (分子生成・実験最適化・ロボット制御) パイプライン。
- `scientific-reproducible-reporting` — 再現可能レポーティングスキル。Quarto 科学文書・ Jupyter Book 多章構成・Papermill パラメトリック実行・ nbconvert 自動変換・Sphinx-Gallery コード例ドキュメント。
- `scientific-research-methodology` — 研究方法論・研究デザインスキル。体系的な研究計画策定、ブレインストーミング フレームワーク、批判的思考法、研究倫理・IRB、先行研究評価、 クロスドメイン着想法を含む研究者のメタスキル群。 「研究計画を立てて」「ブレインストーミングし...
- `scientific-revision-tracker` — 論文改訂の変更履歴追跡・差分管理スキル。原稿バージョン間の diff 生成、 変更箇所のハイライト（赤字削除/青字追加）、改訂サマリー自動生成、 査読コメントと改訂箇所のトレーサビリティ管理を行う。 「改訂をトラッキングして」「変更履...
- `scientific-rrna-taxonomy` — rRNA リファレンス・分類学スキル。SILVA SSU/LSU rRNA データベース・ Greengenes2 系統分類・MGnify メタゲノム解析・QIIME2 分類器・ scikit-bio 配列解析・系統分類パイプライン。
- `scientific-scatac-signac` — scATAC-seq 解析スキル (Signac/SnapATAC2/episcanpy)。 ピークコーリング・モチーフ解析・Gene Activity スコア・ RNA+ATAC マルチモーダル統合 (WNN)。K-Dense: s...
- `scientific-scientific-schematics` — 科学図式・作図スキル。CONSORT フロー図 (臨床試験)、実験プロトコルフロー、 ニューラルネットワークアーキテクチャ図、分子パスウェイ図、 TikZ/SVG/Mermaid ベースの出版品質ベクター図の生成、 AI レビューによ...
- `scientific-scvi-integration` — scvi-tools シングルセル統合スキル。scVI 変分オートエンコーダ統合・ scANVI 半教師有りアノテーション・totalVI CITE-seq RNA+タンパク質結合解析・SOLO ダブレット検出・潜在空間解析。
- `scientific-semantic-scholar` — Semantic Scholar 学術グラフスキル。Semantic Scholar Academic Graph API による論文検索・著者プロファイル・引用グラフ・ 推薦・TLDR 要約。ToolUniverse 連携: sem... （TU: `semantic_scholar`）
- `scientific-semi-supervised-learning` — 半教師あり学習スキル。Self-Training・Label Propagation・ MixMatch/FixMatch・Pseudo-Labeling・ラベル効率評価。
- `scientific-sequence-analysis` — ゲノム配列解析スキル。コドン使用頻度（RSCU/CAI）、ペアワイズアラインメント （Needleman-Wunsch/Smith-Waterman）、系統解析（Jukes-Cantor/UPGMA/ブートストラップ）、 ORF 探索...
- `scientific-single-cell-genomics` — シングルセルゲノミクス解析スキル。scRNA-seq データの品質管理・正規化・ 次元削減（PCA/UMAP）・クラスタリング（Leiden）・差次発現遺伝子（DEG）同定・ セルタイプアノテーション・RNA velocity・細胞間...
- `scientific-spatial-multiomics` — 空間マルチオミクス統合スキル。MERFISH/Visium 等の空間 トランスクリプトームと空間プロテオミクスのマルチモーダル 統合・空間共検出解析・セル近傍グラフ構築パイプライン。
- `scientific-spatial-transcriptomics` — 空間トランスクリプトミクス解析スキル。10x Visium / MERFISH / Slide-seq データの 前処理・空間的遺伝子発現パターン検出（Moran's I / SpatialDE）・ 空間ドメイン同定（BayesSpa...
- `scientific-spectral-signal` — 分光スペクトルおよび生体信号の前処理・解析スキル。ベースライン補正、フィルタリング、 ピーク検出、帯域パワー解析を行う際に使用。 Scientific Skills Exp-08（ECG/EEG）、Exp-11（ラマン分光）で確立した...
- `scientific-squidpy-advanced` — 高度 Squidpy 空間解析スキル。空間自己相関・共起解析・空間 近傍・リガンド受容体空間マッピング・ニッチ同定。 K-Dense 連携: squidpy-advanced。
- `scientific-statistical-simulation` — 統計シミュレーションスキル。Monte Carlo 法・Bootstrap 推論・ Permutation Test・統計的検出力分析・確率的リスク評価。
- `scientific-statistical-testing` — 統計検定・多重比較・エンリッチメント解析のスキル。t検定、カイ二乗検定、ANOVA、 Bonferroni/BH 補正、Fisher 正確検定、ベイズ推論を行う際に使用。 Scientific Skills Exp-03, 04, 0...
- `scientific-stitch-chemical-network` — STITCH 化学-タンパク質相互作用ネットワークスキル。STITCH REST API を用いた化学物質-タンパク質インタラクション検索・ 信頼度スコアリング・ネットワーク薬理学・ポリファーマコロジー解析。 ToolUniverse... （TU: `stitch`）
- `scientific-streaming-analytics` — ストリーミング解析スキル。River オンライン学習・ リアルタイム異常検知・ストリーミング統計・ 増分データ可視化・概念ドリフト検出。
- `scientific-string-network-api` — STRING/BioGRID/STITCH ネットワーク解析スキル。STRING タンパク質相互作用 ネットワーク直接 API、BioGRID 実験的 PPI、STITCH 化学-タンパク質ネットワーク、 ネットワークトポロジー解析・... （TU: `ppi`）
- `scientific-structural-proteomics` — 構造プロテオミクス統合スキル。EMDB クライオ EM、PDBe 構造データ、 Proteins API (UniProt)、Complex Portal 複合体、DeepGO 機能予測、 EVE 変異影響評価の統合パイプライン。
- `scientific-supplementary-generator` — 学術論文の Supplementary Information (SI) を自動生成するスキル。 本文から溢れた図表・手法詳細・追加データを構造化し、ジャーナル規定に準拠した SI ドキュメントを生成する。「SIを作って」「補足資料を...
- `scientific-survival-clinical` — 生存解析と臨床統計のスキル。Kaplan-Meier 曲線、Cox 比例ハザードモデル、Log-rank 検定、 検出力分析、NNT/NNH 算出を行う際に使用。 Scientific Skills Exp-03, 06 で確立したパ...
- `scientific-symbolic-mathematics` — 記号数学スキル。SymPy による解析的微積分・線形代数・微分方程式求解、 記号式の LaTeX 変換、数値計算との統合、科学モデリング用 記号計算パイプライン。
- `scientific-systematic-review` — PRISMA 2020 準拠系統的レビュースキル。マルチ DB 検索戦略立案 (PubMed/Embase/Cochrane/Web of Science)、スクリーニングワークフロー (タイトル/抄録→全文)、品質評価 (RoB 2...
- `scientific-systems-biology` — システム生物学解析スキル。動的モデリング（ODE / SBML）・ 代謝フラックス解析（FBA / pFBA）・遺伝子制御ネットワーク推定（GRN）・ シグナル伝達経路モデリング・パラメータ推定・感度解析・ BioModels/Rea... （TU: `bigg_models`, `complex_portal`, `wikipathways`）
- `scientific-text-mining-nlp` — 科学テキストマイニング・NLP スキル。生物医学 NER（遺伝子/疾患/薬物/化合物）・ 関係抽出（PPI / DDI / GDA）・文献ベースナレッジグラフ構築・ エビデンス要約・トピックモデリング・引用ネットワーク解析パイプライン...
- `scientific-time-series` — 時系列解析・予測スキル。ARIMA/SARIMA/Prophet モデリング、変化点検出（PELT/Bayesian）、 周期解析（FFT/ウェーブレット）、季節分解（STL）、異常検出、Granger 因果性検定の テンプレートを提...
- `scientific-time-series-forecasting` — ML 時系列予測スキル。Prophet/NeuralProphet・N-BEATS・ Temporal Fusion Transformer (TFT)・時系列特徴量エンジニアリング・ バックテスト・多段階予測・アンサンブル予測。
- `scientific-toxicology-env` — 毒性学・環境衛生スキル。CTD (Comparative Toxicogenomics Database) 化学-遺伝子-疾患関連・ToxCast/Tox21 高スループット毒性スクリーニング・ IRIS ヒトリスク評価・T3DB 食...
- `scientific-transfer-learning` — 転移学習・ドメイン適応スキル。事前学習モデルファインチューニング・ Few-shot / Zero-shot 学習・ドメイン適応 (DA)・ 知識蒸留・マルチタスク学習・科学ドメイン特化モデル転移。
- `scientific-uncertainty-quantification` — 不確実性定量化スキル。Conformal Prediction・MC Dropout・ 深層アンサンブル・アレアトリック / エピステミック分離・ Calibration Curve・予測区間推定・Expected Calibrati...
- `scientific-uniprot-proteome` — UniProt プロテオームスキル。UniProt REST API による タンパク質検索・ID マッピング・配列取得・機能アノテーション・ UniRef/UniParc 横断検索。ToolUniverse 連携: uniprot。 （TU: `uniprot`）
- `scientific-variant-effect-prediction` — 計算バリアント効果予測スキル。AlphaMissense (タンパク質構造ベース病原性予測)、 CADD (統合アノテーションスコア)、SpliceAI (スプライシング影響予測) の 3 大予測ツールを統合したコンセンサス病原性評価... （TU: `spliceai`, `cadd`）
- `scientific-variant-interpretation` — 遺伝子バリアント臨床解釈スキル。ClinVar / gnomAD / COSMIC / ACMG ガイドラインに 基づく病原性評価、薬理ゲノミクス（PharmGKB/ClinPGx）、バリアント-表現型相関の エビデンスグレーディング... （TU: `clinvar`）

</details>

# ToolUniverse key → スキル

ToolUniverse API キーから対応スキルを検索できます。1 キーが複数スキルにまたがる場合もあります。
v0.26.0 で全 190 スキルに TU 連携を追加。汎用キー `biotools` / `openml` / `papers_with_code` / `open_alex` はドメイン横断で多数のスキルが共有する。

| TU key | 対応スキル |
|--------|-----------|
| `alphafold` | `scientific-alphafold-structures` |
| `arrayexpress` | `scientific-arrayexpress-expression` |
| `bigg` | `scientific-metabolic-flux` |
| `bigg_models` | `scientific-systems-biology` |
| `bindingdb` | `scientific-pharmacology-targets` |
| `biomodels` | `scientific-metabolic-modeling` |
| `biothings` | `scientific-biothings-idmapping` |
| `biotools` | `scientific-adaptive-experiments`, `scientific-advanced-imaging`, `scientific-advanced-visualization`, `scientific-bayesian-statistics`, `scientific-biosignal-processing`, `scientific-data-preprocessing`, `scientific-data-profiling`, `scientific-data-simulation`, `scientific-doe`, `scientific-eda-correlation`, `scientific-image-analysis`, `scientific-interactive-dashboard`, `scientific-lab-automation`, `scientific-missing-data-analysis`, `scientific-neuroscience-electrophysiology`, `scientific-pca-tsne`, `scientific-pipeline-scaffold`, `scientific-presentation-design`, `scientific-process-optimization`, `scientific-publication-figures`, `scientific-reproducible-reporting`, `scientific-spectral-signal`, `scientific-statistical-simulation`, `scientific-statistical-testing`, `scientific-streaming-analytics`, `scientific-symbolic-mathematics`, `scientific-time-series`, `scientific-time-series-forecasting` |
| `brenda` | `scientific-pharmacology-targets` |
| `cadd` | `scientific-variant-effect-prediction` |
| `cbioportal` | `scientific-cancer-genomics` |
| `cellosaurus` | `scientific-cell-line-resources` |
| `cellxgene` | `scientific-gpu-singlecell`, `scientific-perturbation-analysis`, `scientific-scvi-integration`, `scientific-spatial-multiomics`, `scientific-squidpy-advanced` |
| `cellxgene_census` | `scientific-cellxgene-census`, `scientific-human-cell-atlas` |
| `chembl` | `scientific-chembl-assay-mining`, `scientific-deep-chemistry` |
| `chipatlas` | `scientific-encode-screen`, `scientific-epigenomics-chromatin` |
| `civic` | `scientific-civic-evidence` |
| `clingen` | `scientific-clingen-curation` |
| `clinvar` | `scientific-biobank-cohort`, `scientific-variant-interpretation` |
| `complex_portal` | `scientific-systems-biology` |
| `cosmic` | `scientific-cancer-genomics` |
| `crossref` | `scientific-academic-writing`, `scientific-critical-review`, `scientific-crossref-metadata`, `scientific-latex-formatter`, `scientific-paper-quality`, `scientific-peer-review-response`, `scientific-revision-tracker`, `scientific-supplementary-generator` |
| `ctd` | `scientific-toxicology-env` |
| `depmap` | `scientific-depmap-dependencies` |
| `dgidb` | `scientific-drug-target-profiling` |
| `disgenet` | `scientific-disease-research` |
| `drugbank` | `scientific-clinical-pharmacology`, `scientific-drugbank-resources` |
| `ena` | `scientific-data-submission` |
| `encode` | `scientific-encode-screen`, `scientific-scatac-signac` |
| `ensembl` | `scientific-crispr-design` |
| `fda_pharmacogenomic_biomarkers` | `scientific-pharmacogenomics` |
| `gbif` | `scientific-environmental-ecology`, `scientific-geospatial-analysis`, `scientific-marine-ecology` |
| `gdc` | `scientific-gdc-portal` |
| `geo` | `scientific-geo-expression` |
| `glygen` | `scientific-glycomics` |
| `gnomad` | `scientific-gnomad-variants` |
| `gtex_v2` | `scientific-gtex-tissue-expression` |
| `gtopdb` | `scientific-pharmacology-targets` |
| `gwas` | `scientific-gwas-catalog` |
| `hca_tools` | `scientific-human-cell-atlas` |
| `hgnc` | `scientific-hgnc-nomenclature` |
| `hmdb` | `scientific-metabolomics`, `scientific-metabolomics-network` |
| `hpa` | `scientific-human-protein-atlas` |
| `icd` | `scientific-clinical-standards` |
| `iedb` | `scientific-immunoinformatics` |
| `imgt` | `scientific-immunoinformatics` |
| `impc` | `scientific-model-organism-db` |
| `intact` | `scientific-protein-interaction-network` |
| `lipidmaps` | `scientific-lipidomics` |
| `loinc` | `scientific-clinical-standards` |
| `materials_project` | `scientific-materials-characterization` |
| `medlineplus` | `scientific-public-health-data` |
| `metabolic_atlas` | `scientific-metabolic-atlas` |
| `metabolomics_workbench` | `scientific-metabolomics` |
| `metacyc` | `scientific-metabolomics-databases` |
| `mgnify` | `scientific-metagenome-assembled-genomes`, `scientific-microbiome-metagenomics` |
| `monarch` | `scientific-monarch-ontology` |
| `mpd` | `scientific-model-organism-db` |
| `ncbi_taxonomy` | `scientific-phylogenetics` |
| `nci60` | `scientific-nci60-screening` |
| `nhanes` | `scientific-public-health-data` |
| `obis` | `scientific-marine-ecology` |
| `odphp` | `scientific-public-health-data` |
| `oncokb` | `scientific-precision-oncology` |
| `open_alex` | `scientific-causal-inference`, `scientific-hypothesis-pipeline`, `scientific-research-methodology` |
| `openml` | `scientific-active-learning`, `scientific-anomaly-detection`, `scientific-automl`, `scientific-causal-ml`, `scientific-ensemble-methods`, `scientific-feature-importance`, `scientific-ml-classification`, `scientific-ml-regression`, `scientific-model-monitoring`, `scientific-multi-task-learning`, `scientific-semi-supervised-learning` |
| `opentarget` | `scientific-opentargets-genetics` |
| `orphanet` | `scientific-rare-disease-genetics` |
| `paleobiology` | `scientific-paleobiology` |
| `papers_with_code` | `scientific-deep-learning`, `scientific-explainable-ai`, `scientific-federated-learning`, `scientific-neural-architecture-search`, `scientific-quantum-computing`, `scientific-reinforcement-learning`, `scientific-transfer-learning`, `scientific-uncertainty-quantification` |
| `pdb` | `scientific-md-simulation` |
| `pharmgkb` | `scientific-pharmgkb-pgx` |
| `pharos` | `scientific-drug-repurposing`, `scientific-pharos-targets` |
| `ppi` | `scientific-string-network-api` |
| `proteinsplus` | `scientific-protein-structure-analysis` |
| `pubchem` | `scientific-admet-pharmacokinetics` |
| `rcsb_pdb` | `scientific-rcsb-pdb-search` |
| `rcsb_search` | `scientific-rcsb-pdb-search` |
| `reactome` | `scientific-reactome-pathways` |
| `sabdab` | `scientific-immunoinformatics` |
| `semantic_scholar` | `scientific-semantic-scholar` |
| `spliceai` | `scientific-variant-effect-prediction` |
| `stitch` | `scientific-stitch-chemical-network` |
| `string` | `scientific-network-visualization` |
| `tair` | `scientific-plant-biology` |
| `tcia` | `scientific-medical-imaging`, `scientific-radiology-ai` |
| `therasabdab` | `scientific-immunoinformatics` |
| `umls` | `scientific-clinical-nlp` |
| `uniprot` | `scientific-uniprot-proteome` |
| `who_gho` | `scientific-epidemiology-public-health` |
| `wikipathways` | `scientific-systems-biology` |
| `worms` | `scientific-marine-ecology` |
| `zinc` | `scientific-compound-screening` |

# 成果物パス（results/figures/...）→ スキル

各スキルの SKILL.md 内で参照されている出力パスです。末尾 `_` で終わるエントリはプレフィックスパターンです。

<details><summary>data/ パス（8 件）</summary>

- `data/bbbp`: `scientific-graph-neural-networks`
- `data/dataset_processed.csv`: `scientific-hypothesis-pipeline`
- `data/dataset.csv`: `scientific-hypothesis-pipeline`
- `data/geo`: `scientific-gene-expression-transcriptomics`
- `data/hetionet`: `scientific-graph-neural-networks`
- `data/participants`: `scientific-reactome-pathways`
- `data/pathway`: `scientific-plant-biology`, `scientific-reactome-pathways`
- `data/pathways/low/entity`: `scientific-reactome-pathways`

</details>

<details><summary>figures/ パス（155 件）</summary>

- `figures/actual_vs_predicted.png`: `scientific-ml-regression`
- `figures/admixture_barplot.png`: `scientific-population-genetics`
- `figures/alpha_boxplot.png`: `scientific-microbiome-metagenomics`
- `figures/band_structure.png`: `scientific-computational-materials`
- `figures/barplot_taxonomy.png`: `scientific-microbiome-metagenomics`
- `figures/bayesian_convergence.png`: `scientific-doe`
- `figures/bayesian_model_comparison.png`: `scientific-bayesian-statistics`
- `figures/bayesian_ppc.png`: `scientific-bayesian-statistics`
- `figures/bayesian_trace.png`: `scientific-bayesian-statistics`
- `figures/bayesian_update.png`: `scientific-survival-clinical`
- `figures/causal_dag.png`: `scientific-causal-inference`
- `figures/cca_scores.png`: `scientific-multi-omics`
- `figures/changepoints.png`: `scientific-time-series`
- `figures/chemical_space_pca.png`: `scientific-cheminformatics`
- `figures/chromatin_state_heatmap.png`: `scientific-epigenomics-chromatin`
- `figures/citation_network.png`: `scientific-literature-search`, `scientific-text-mining-nlp`
- `figures/clustering_analysis.png`: `scientific-pca-tsne`
- `figures/codon_usage_heatmap.png`: `scientific-sequence-analysis`
- `figures/confusion_matrices.png`: `scientific-ml-classification`
- `figures/consort_flow.svg`: `scientific-scientific-schematics`
- `figures/correlation_heatmap.png`: `scientific-eda-correlation`
- `figures/cox_ph_forest.png`: `scientific-survival-clinical`
- `figures/cumulative_meta.png`: `scientific-meta-analysis`
- `figures/dag_diagram.png`: `scientific-epidemiology-public-health`
- `figures/deg_dotplot.png`: `scientific-single-cell-genomics`
- `figures/disease_map.png`: `scientific-epidemiology-public-health`
- `figures/distribution_boxplots.png`: `scientific-eda-correlation`
- `figures/diversity_comparison.png`: `scientific-environmental-ecology`
- `figures/dl_learning_curve.png`: `scientific-deep-learning`
- `figures/docking_scores.png`: `scientific-molecular-docking`
- `figures/domain_architecture.png`: `scientific-protein-domain-family`
- `figures/dos.png`: `scientific-computational-materials`
- `figures/eda_`: `scientific-hypothesis-pipeline`
- `figures/eeg_band_powers.png`: `scientific-spectral-signal`
- `figures/enrichment_dotplot.png`: `scientific-statistical-testing`
- `figures/enrichment_heatmap.png`: `scientific-pathway-enrichment`
- `figures/epidemic_curves.png`: `scientific-infectious-disease`
- `figures/epitope_map.png`: `scientific-immunoinformatics`
- `figures/erp_waveforms.png`: `scientific-neuroscience-electrophysiology`
- `figures/expression_heatmap.png`: `scientific-expression-comparison`
- `figures/feature_importance_`: `scientific-feature-importance`
- `figures/feature_importance_panel.png`: `scientific-feature-importance`
- `figures/fft_spectrum.png`: `scientific-time-series`
- `figures/Fig1_description.png`: `scientific-academic-writing`
- `figures/fig1_overview.png`: `scientific-academic-writing`
- `figures/fig1.png`: `scientific-academic-writing`
- `figures/fig2_composite.png`: `scientific-academic-writing`
- `figures/figure1.png`: `scientific-academic-writing`
- `figures/filename.png`: `scientific-academic-writing`
- `figures/fluorescence_merged.png`: `scientific-image-analysis`
- `figures/flux_map.png`: `scientific-systems-biology`
- `figures/forest_plot.png`: `scientific-epidemiology-public-health`, `scientific-meta-analysis`
- `figures/funnel_plot.png`: `scientific-meta-analysis`
- `figures/gnn_explanation.png`: `scientific-graph-neural-networks`
- `figures/gnn_training_curve.png`: `scientific-graph-neural-networks`
- `figures/grn_graph.png`: `scientific-systems-biology`
- `figures/gsea_dotplot.png`: `scientific-gene-expression-transcriptomics`
- `figures/gsea_running_sum.png`: `scientific-pathway-enrichment`
- `figures/hic_contact_map.png`: `scientific-epigenomics-chromatin`
- `figures/hrv_poincare.png`: `scientific-neuroscience-electrophysiology`
- `figures/hydrophobicity_profile.png`: `scientific-sequence-analysis`
- `figures/interaction_plot.png`: `scientific-doe`
- `figures/kaplan_meier.png`: `scientific-statistical-testing`, `scientific-survival-clinical`
- `figures/kg_visualization.png`: `scientific-text-mining-nlp`
- `figures/ma_plot.png`: `scientific-gene-expression-transcriptomics`
- `figures/main_effects_plot.png`: `scientific-doe`
- `figures/manhattan_fst.png`: `scientific-population-genetics`
- `figures/metabolite_class_distribution.png`: `scientific-metabolomics-databases`
- `figures/metabolite_network.png`: `scientific-metabolomics`
- `figures/model_comparison_r2.png`: `scientific-ml-regression`
- `figures/molecular_network.png`: `scientific-proteomics-mass-spectrometry`
- `figures/morans_i_dist.png`: `scientific-spatial-transcriptomics`
- `figures/multiomics_umap.png`: `scientific-multi-omics`
- `figures/mutation_spectrum.png`: `scientific-cancer-genomics`
- `figures/network_visualization.png`: `scientific-bioinformatics`, `scientific-network-analysis`
- `figures/nmds_plot.png`: `scientific-environmental-ecology`
- `figures/nn_architecture.svg`: `scientific-scientific-schematics`
- `figures/optuna_history.png`: `scientific-deep-learning`
- `figures/optuna_importance.png`: `scientific-deep-learning`
- `figures/pareto_front.png`: `scientific-process-optimization`
- `figures/partial_dependence.png`: `scientific-ml-classification`
- `figures/particle_size_distribution.png`: `scientific-image-analysis`
- `figures/patent_timeline.png`: `scientific-regulatory-science`
- `figures/pathway.md`: `scientific-scientific-schematics`
- `figures/pca_loadings.png`: `scientific-pca-tsne`
- `figures/pca_populations.png`: `scientific-population-genetics`
- `figures/pca_screeplot.png`: `scientific-pca-tsne`
- `figures/pca_tsne_panel.png`: `scientific-pca-tsne`
- `figures/pcoa_plot.png`: `scientific-microbiome-metagenomics`
- `figures/pdp_`: `scientific-feature-importance`
- `figures/peak_detection.png`: `scientific-spectral-signal`
- `figures/permutation_importance_`: `scientific-feature-importance`
- `figures/pgx_phenotype_summary.png`: `scientific-pharmacogenomics`
- `figures/phase_diagram.png`: `scientific-computational-materials`
- `figures/phylogenetic_tree.png`: `scientific-infectious-disease`, `scientific-phylogenetics`, `scientific-sequence-analysis`
- `figures/plsda_scores.png`: `scientific-metabolomics`
- `figures/poincare_plot.png`: `scientific-biosignal-processing`
- `figures/ppi_network.png`: `scientific-protein-interaction-network`
- `figures/precision_recall_curves.png`: `scientific-ml-classification`
- `figures/prisma_flow.mmd`: `scientific-systematic-review`
- `figures/prisma_flow.svg`: `scientific-systematic-review`
- `figures/process_window.png`: `scientific-process-optimization`
- `figures/propensity_distribution.png`: `scientific-causal-inference`
- `figures/psp_path_diagram.png`: `scientific-network-analysis`
- `figures/pubtator_dashboard.png`: `scientific-biomedical-pubtator`
- `figures/pv_demographics.png`: `scientific-pharmacovigilance`
- `figures/pv_temporal_trend.png`: `scientific-pharmacovigilance`
- `figures/quantum_convergence.png`: `scientific-quantum-computing`
- `figures/quantum_rabi.png`: `scientific-quantum-computing`
- `figures/radar_`: `scientific-ml-regression`
- `figures/rdd_plot.png`: `scientific-causal-inference`
- `figures/repertoire_clonality.png`: `scientific-immunoinformatics`
- `figures/response_surface_`: `scientific-process-optimization`
- `figures/rl_reward_curve.png`: `scientific-reinforcement-learning`
- `figures/roc_curves.png`: `scientific-ml-classification`
- `figures/scatter_matrix.png`: `scientific-eda-correlation`
- `figures/sdm_map.png`: `scientific-environmental-ecology`
- `figures/segmentation_overlay.png`: `scientific-medical-imaging`
- `figures/segmentation_result.png`: `scientific-image-analysis`
- `figures/shap_importance.png`: `scientific-explainable-ai`
- `figures/shap_interaction.png`: `scientific-explainable-ai`
- `figures/shap_summary.png`: `scientific-explainable-ai`
- `figures/shap_waterfall.png`: `scientific-explainable-ai`
- `figures/signature_profiles.png`: `scientific-cancer-genomics`
- `figures/similarity_heatmap.png`: `scientific-cheminformatics`
- `figures/snf_heatmap.png`: `scientific-multi-omics`
- `figures/spatial_domains.png`: `scientific-spatial-transcriptomics`
- `figures/spatial_svgs.png`: `scientific-spatial-transcriptomics`
- `figures/spectrogram_`: `scientific-biosignal-processing`
- `figures/spectrum_processed.png`: `scientific-spectral-signal`
- `figures/spike_rasters.png`: `scientific-neuroscience-electrophysiology`
- `figures/stl_decomposition.png`: `scientific-time-series`
- `figures/structure_zone_model.png`: `scientific-materials-characterization`
- `figures/symbolic_plot.png`: `scientific-symbolic-mathematics`
- `figures/tauc_plot.png`: `scientific-materials-characterization`
- `figures/tikz_figure.pdf`: `scientific-scientific-schematics`
- `figures/tikz_figure.tex`: `scientific-scientific-schematics`
- `figures/time_series_forecast.png`: `scientific-time-series`
- `figures/timecourse_plot.png`: `scientific-systems-biology`
- `figures/topic_distribution.png`: `scientific-text-mining-nlp`
- `figures/transmission_network.png`: `scientific-infectious-disease`
- `figures/trial_phase_distribution.png`: `scientific-clinical-trials-analytics`
- `figures/trial_timeline.png`: `scientific-clinical-trials-analytics`
- `figures/tsne_projection.png`: `scientific-pca-tsne`
- `figures/umap_celltypes.png`: `scientific-single-cell-genomics`
- `figures/umap_clusters.png`: `scientific-bioinformatics`
- `figures/umap_leiden.png`: `scientific-single-cell-genomics`
- `figures/variant_score_distribution.png`: `scientific-variant-effect-prediction`
- `figures/velocity_stream.png`: `scientific-single-cell-genomics`
- `figures/vip_barplot.png`: `scientific-metabolomics`
- `figures/volcano_plot.png`: `scientific-ml-classification`
- `figures/volcano_proteomics.png`: `scientific-proteomics-mass-spectrometry`
- `figures/volcano_rnaseq.png`: `scientific-gene-expression-transcriptomics`
- `figures/workflow_schematic.png`: `scientific-presentation-design`
- `figures/xrd_analysis.png`: `scientific-materials-characterization`

</details>

<details><summary>results/ パス A–L（259 件）</summary>

- `results/4dn_contacts.json`: `scientific-regulatory-genomics`
- `results/510k_clearances.csv`: `scientific-regulatory-science`
- `results/adata_spatial.h5ad`: `scientific-squidpy-advanced`
- `results/admet_profile.json`: `scientific-admet-pharmacokinetics`
- `results/admet_report.md`: `scientific-admet-pharmacokinetics`
- `results/admixture_Q.csv`: `scientific-population-genetics`
- `results/AF`: `scientific-alphafold-structures`
- `results/age_standardized_rates.csv`: `scientific-epidemiology-public-health`
- `results/alignment.csv`: `scientific-spatial-multiomics`
- `results/allele_frequencies.csv`: `scientific-population-genetics`
- `results/alpha_diversity.csv`: `scientific-microbiome-metagenomics`
- `results/alphamissense_scores.csv`: `scientific-variant-effect-prediction`
- `results/amr_genes.csv`: `scientific-infectious-disease`
- `results/analysis_summary.json`: `scientific-academic-writing`, `scientific-hypothesis-pipeline`
- `results/ancestral_sequences.fasta`: `scientific-phylogenetics`
- `results/annotations.csv`: `scientific-parasite-genomics`, `scientific-pharmgkb-pgx`, `scientific-scvi-integration`
- `results/anova_factor_effects.csv`: `scientific-doe`
- `results/antibody_structure.json`: `scientific-immunoinformatics`
- `results/arxiv_papers.csv`: `scientific-preprint-archive`
- `results/asv_table.csv`: `scientific-microbiome-metagenomics`
- `results/atacseq`: `scientific-epigenomics-chromatin`
- `results/atacseq/fragment_size_dist.csv`: `scientific-epigenomics-chromatin`
- `results/augur_scores.csv`: `scientific-perturbation-analysis`
- `results/baseline_expression.csv`: `scientific-expression-comparison`
- `results/bayesian_optimization_history.csv`: `scientific-doe`
- `results/bayesian_summary.json`: `scientific-bayesian-statistics`
- `results/bcell_epitopes.csv`: `scientific-immunoinformatics`
- `results/benchling_registry.json`: `scientific-lab-data-management`
- `results/benchling_sequences.json`: `scientific-lab-data-management`
- `results/benchmark.json`: `scientific-deep-chemistry`, `scientific-gpu-singlecell`
- `results/beta_distance_matrix.csv`: `scientific-microbiome-metagenomics`
- `results/bigg_model.json`: `scientific-metabolic-modeling`
- `results/bigg_reaction.json`: `scientific-metabolic-modeling`
- `results/bigg_search.csv`: `scientific-metabolic-modeling`
- `results/binding_sites.json`: `scientific-protein-structure-analysis`
- `results/bindingdb_ligands.csv`: `scientific-pharmacology-targets`
- `results/biodiversity_indices.csv`: `scientific-environmental-ecology`
- `results/biomodels_model.json`: `scientific-metabolic-modeling`
- `results/biomodels_search.csv`: `scientific-metabolic-modeling`
- `results/biostudies_metadata.json`: `scientific-ebi-databases`
- `results/blast_results.json`: `scientific-genome-sequence-tools`
- `results/cadd_scores.csv`: `scientific-variant-effect-prediction`
- `results/cancer_stats.csv`: `scientific-icgc-cancer-data`
- `results/canonical_correlations.csv`: `scientific-multi-omics`
- `results/capa_record.json`: `scientific-regulatory-science`
- `results/causal_estimates.csv`: `scientific-causal-inference`
- `results/cbioportal_mutations.csv`: `scientific-cancer-genomics`
- `results/cell_lines.csv`: `scientific-cell-line-resources`
- `results/cell_painting.csv`: `scientific-advanced-imaging`
- `results/cell_type_composition.csv`: `scientific-human-cell-atlas`
- `results/celltype_distribution.csv`: `scientific-cellxgene-census`
- `results/census_datasets.csv`: `scientific-cellxgene-census`
- `results/census_expression.h5ad`: `scientific-cellxgene-census`
- `results/centrality_measures.csv`: `scientific-bioinformatics`, `scientific-network-analysis`
- `results/centrality_scores.csv`: `scientific-squidpy-advanced`
- `results/changepoints.csv`: `scientific-time-series`
- `results/chembl_activities.csv`: `scientific-chembl-assay-mining`
- `results/chipatlas_enrichment.csv`: `scientific-encode-screen`
- `results/chipseq`: `scientific-epigenomics-chromatin`
- `results/chromhmm/emissions_`: `scientific-epigenomics-chromatin`
- `results/citation_metrics.csv`: `scientific-text-mining-nlp`
- `results/citation_network.graphml`: `scientific-literature-search`
- `results/citation_stats.csv`: `scientific-crossref-metadata`
- `results/citations.csv`: `scientific-semantic-scholar`
- `results/civic_assertions.csv`: `scientific-civic-evidence`
- `results/civic_evidence.csv`: `scientific-civic-evidence`
- `results/civic_gene.csv`: `scientific-civic-evidence`
- `results/civic_variants.csv`: `scientific-civic-evidence`
- `results/classification_metrics.csv`: `scientific-ml-classification`
- `results/classification.csv`: `scientific-marine-ecology`
- `results/clingen_actionability.csv`: `scientific-clingen-curation`
- `results/clingen_dosage.csv`: `scientific-clingen-curation`
- `results/clingen_validity.csv`: `scientific-clingen-curation`
- `results/clinical_decision_report.md`: `scientific-clinical-decision-support`
- `results/clinical_metrics.json`: `scientific-healthcare-ai`
- `results/clinical_ner.csv`: `scientific-clinical-nlp`
- `results/clinical_predictions.csv`: `scientific-healthcare-ai`
- `results/clinical_recommendation.json`: `scientific-clinical-decision-support`
- `results/clinical_sections.csv`: `scientific-clinical-nlp`
- `results/clinical_trials_export.csv`: `scientific-clinical-trials-analytics`
- `results/clinical_trials_search.csv`: `scientific-clinical-trials-analytics`
- `results/co_occurrence_matrix.csv`: `scientific-squidpy-advanced`
- `results/code_mapping.json`: `scientific-healthcare-ai`
- `results/codetection.csv`: `scientific-spatial-multiomics`
- `results/codon_usage.csv`: `scientific-bioinformatics`
- `results/communities.csv`: `scientific-spatial-multiomics`
- `results/community_pathway_mapping.csv`: `scientific-bioinformatics`
- `results/competitive_landscape.json`: `scientific-clinical-trials-analytics`
- `results/complex_portal.json`: `scientific-structural-proteomics`
- `results/connectivity/conn_matrix.csv`: `scientific-neuroscience-electrophysiology`
- `results/consensus_pathogenicity.csv`: `scientific-variant-effect-prediction`
- `results/consensus.csv`: `scientific-rrna-taxonomy`
- `results/conservation_priority.csv`: `scientific-environmental-ecology`
- `results/contamination_report.json`: `scientific-cell-line-resources`
- `results/cosmic_mutations.csv`: `scientific-cancer-genomics`
- `results/covariate_balance.csv`: `scientific-causal-inference`
- `results/cox_ph_results.csv`: `scientific-survival-clinical`
- `results/cpg_islands.csv`: `scientific-sequence-analysis`
- `results/cross_omics_correlation.csv`: `scientific-multi-omics`
- `results/cross_species.json`: `scientific-model-organism-db`
- `results/crossref_results.csv`: `scientific-literature-search`
- `results/ctd_disease_associations.csv`: `scientific-toxicology-env`
- `results/ctd_gene_interactions.csv`: `scientific-toxicology-env`
- `results/da_results.csv`: `scientific-microbiome-metagenomics`
- `results/dag_analysis.json`: `scientific-epidemiology-public-health`
- `results/data_extraction.csv`: `scientific-systematic-review`
- `results/data_quality.json`: `scientific-hypothesis-pipeline`
- `results/dbsnp_frequencies.csv`: `scientific-genome-sequence-tools`
- `results/dbsnp_variant.json`: `scientific-genome-sequence-tools`
- `results/de_results.csv`: `scientific-scvi-integration`
- `results/deconvolution_proportions.csv`: `scientific-spatial-transcriptomics`
- `results/deepgo_predictions.json`: `scientific-structural-proteomics`
- `results/deg_results.csv`: `scientific-geo-expression`
- `results/depmap_dependencies.csv`: `scientific-cancer-genomics`
- `results/depmap_dependency.csv`: `scientific-depmap-dependencies`, `scientific-nci60-screening`
- `results/depmap_drugs.csv`: `scientific-depmap-dependencies`
- `results/depmap_top_dependent.csv`: `scientific-depmap-dependencies`
- `results/descriptive_statistics.csv`: `scientific-eda-correlation`
- `results/deseq2_results.csv`: `scientific-gene-expression-transcriptomics`
- `results/design_candidates.json`: `scientific-protein-design`
- `results/design_report.md`: `scientific-protein-design`
- `results/device_classification.csv`: `scientific-regulatory-science`
- `results/diffbind_results.csv`: `scientific-epigenomics-chromatin`
- `results/diffdock`: `scientific-molecular-docking`
- `results/diffdock/rank`: `scientific-molecular-docking`
- `results/differential_expression.csv`: `scientific-expression-comparison`
- `results/differential_proteins.csv`: `scientific-proteomics-mass-spectrometry`
- `results/disease_research_report.json`: `scientific-disease-research`
- `results/disease_research_report.md`: `scientific-disease-research`
- `results/disgenet_gda.csv`: `scientific-rare-disease-genetics`
- `results/divergence_times.json`: `scientific-phylogenetics`
- `results/diversity.csv`: `scientific-paleobiology`
- `results/dl_training_log.json`: `scientific-deep-learning`
- `results/dnanexus_workflow_output.json`: `scientific-lab-data-management`
- `results/docking_results.csv`: `scientific-molecular-docking`
- `results/doi_resolved.csv`: `scientific-crossref-metadata`
- `results/domain_comparison.json`: `scientific-protein-domain-family`
- `results/dominant_periods.csv`: `scientific-time-series`
- `results/dosing_recommendations.csv`: `scientific-pharmacogenomics`
- `results/drug_activity.csv`: `scientific-nci60-screening`
- `results/drug_mapping.json`: `scientific-public-health-data`
- `results/drug_targets.csv`: `scientific-parasite-genomics`
- `results/drugbank_ddi.csv`: `scientific-drugbank-resources`
- `results/drugbank_detail.csv`: `scientific-drugbank-resources`
- `results/drugbank_search.csv`: `scientific-drugbank-resources`
- `results/drugbank_targets.csv`: `scientific-drugbank-resources`
- `results/druggability_matrix.json`: `scientific-drug-target-profiling`
- `results/drugs.csv`: `scientific-pharmgkb-pgx`
- `results/ebi_search.csv`: `scientific-ebi-databases`
- `results/ecg/hrv_results.csv`: `scientific-neuroscience-electrophysiology`
- `results/eda/scr_peaks.csv`: `scientific-neuroscience-electrophysiology`
- `results/edge_list.csv`: `scientific-network-analysis`
- `results/eeg/erp_evokeds.fif`: `scientific-neuroscience-electrophysiology`
- `results/eeg/microstates.csv`: `scientific-neuroscience-electrophysiology`
- `results/effect_sizes.csv`: `scientific-meta-analysis`
- `results/efo_terms.csv`: `scientific-ontology-enrichment`
- `results/embeddings.npy`: `scientific-deep-chemistry`
- `results/emdb_structure.json`: `scientific-structural-proteomics`
- `results/ena_sequences.fasta`: `scientific-ebi-databases`
- `results/encode_experiments.csv`: `scientific-encode-screen`
- `results/enrichr_results`: `scientific-ontology-enrichment`
- `results/ensembl_gene_info.json`: `scientific-ensembl-genomics`
- `results/entity_linking.csv`: `scientific-clinical-nlp`
- `results/entity_network.graphml`: `scientific-biomedical-pubtator`
- `results/entity_relations.csv`: `scientific-biomedical-pubtator`
- `results/env_stack.csv`: `scientific-environmental-geodata`
- `results/env_summary.csv`: `scientific-environmental-geodata`
- `results/epidemic_simulation.csv`: `scientific-infectious-disease`
- `results/eqtl_results.csv`: `scientific-gtex-tissue-expression`
- `results/esm_scores.json`: `scientific-protein-design`
- `results/europepmc_results.csv`: `scientific-literature-search`
- `results/eve_scores.json`: `scientific-structural-proteomics`
- `results/experiment_files.csv`: `scientific-arrayexpress-expression`
- `results/experimental_design.csv`: `scientific-doe`
- `results/experiments.csv`: `scientific-arrayexpress-expression`
- `results/expression_`: `scientific-gtex-tissue-expression`
- `results/expression_matrix.csv`: `scientific-expression-comparison`, `scientific-gtex-tissue-expression`
- `results/fba_fluxes.csv`: `scientific-systems-biology`
- `results/fcs_processed.csv`: `scientific-healthcare-ai`
- `results/fda_orange_book.json`: `scientific-regulatory-science`
- `results/feature_importance.csv`: `scientific-feature-importance`
- `results/features_detected.csv`: `scientific-proteomics-mass-spectrometry`
- `results/features.csv`: `scientific-advanced-imaging`
- `results/fluxes.csv`: `scientific-metabolic-flux`
- `results/forecast_results.csv`: `scientific-time-series`
- `results/fst_per_snp.csv`: `scientific-population-genetics`
- `results/fulltext_corpus`: `scientific-preprint-archive`
- `results/gbif_occurrences.csv`: `scientific-marine-ecology`
- `results/gdc_cases.csv`: `scientific-gdc-portal`
- `results/gdc_mutations.csv`: `scientific-genome-sequence-tools`
- `results/gdc_projects.csv`: `scientific-gdc-portal`
- `results/gdc_ssm.csv`: `scientific-gdc-portal`
- `results/gene_annotations.csv`: `scientific-pharmgkb-pgx`
- `results/gene_drug_interactions.csv`: `scientific-pharmacogenomics`
- `results/genefamilies.tsv`: `scientific-microbiome-metagenomics`
- `results/genes.csv`: `scientific-parasite-genomics`
- `results/geo_expression_matrix.csv`: `scientific-gene-expression-transcriptomics`
- `results/geo_summary.csv`: `scientific-paleobiology`
- `results/glcm_texture_features.csv`: `scientific-image-analysis`
- `results/glycan_details.csv`: `scientific-glycomics`
- `results/glycosites.csv`: `scientific-glycomics`
- `results/gnn_benchmark.json`: `scientific-graph-neural-networks`
- `results/gnn_predictions.json`: `scientific-graph-neural-networks`
- `results/gnomad_constraint.csv`: `scientific-gnomad-variants`
- `results/gnomad_rare.csv`: `scientific-gnomad-variants`
- `results/gnomad_region.csv`: `scientific-gnomad-variants`
- `results/go_enrichment.csv`: `scientific-pathway-enrichment`
- `results/gpcrdb_profile.json`: `scientific-pharmacology-targets`
- `results/gpu_singlecell.h5ad`: `scientific-gpu-singlecell`
- `results/granger_causality.csv`: `scientific-time-series`
- `results/grn_network.csv`: `scientific-systems-biology`
- `results/gsea`: `scientific-gene-expression-transcriptomics`
- `results/gsea_results.csv`: `scientific-pathway-enrichment`
- `results/gtopdb_interactions.csv`: `scientific-pharmacology-targets`
- `results/guidelines.csv`: `scientific-pharmgkb-pgx`
- `results/gwas_associations.csv`: `scientific-gwas-catalog`
- `results/gwas_significant_loci.json`: `scientific-disease-research`
- `results/gwas_significant.csv`: `scientific-biobank-cohort`
- `results/gwas_studies.csv`: `scientific-gwas-catalog`
- `results/hbond_analysis.csv`: `scientific-md-simulation`
- `results/hca_atlas.h5ad`: `scientific-human-cell-atlas`
- `results/hca_projects.csv`: `scientific-human-cell-atlas`
- `results/health_disparities.csv`: `scientific-public-health-data`
- `results/health_guidelines.json`: `scientific-public-health-data`
- `results/hgnc_alias_resolved.csv`: `scientific-hgnc-nomenclature`
- `results/hgnc_details.csv`: `scientific-hgnc-nomenclature`
- `results/hgnc_xref.csv`: `scientific-hgnc-nomenclature`
- `results/hic/compartments.csv`: `scientific-epigenomics-chromatin`
- `results/hic/tad_boundaries.bed`: `scientific-epigenomics-chromatin`
- `results/hmdb_metabolites.csv`: `scientific-metabolomics-databases`
- `results/homology_table.csv`: `scientific-ensembl-genomics`
- `results/hpa_cancer_prognostics.csv`: `scientific-human-protein-atlas`
- `results/hpa_gene_info.json`: `scientific-human-protein-atlas`
- `results/hpa_interactions.csv`: `scientific-human-protein-atlas`
- `results/hpa_subcellular.csv`: `scientific-human-protein-atlas`
- `results/hpa_tissue_expression.csv`: `scientific-human-protein-atlas`
- `results/hrv_metrics.csv`: `scientific-spectral-signal`
- `results/hub_metabolites.csv`: `scientific-metabolic-atlas`, `scientific-metabolomics-network`
- `results/hypothesis_tests.csv`: `scientific-hypothesis-pipeline`
- `results/hypothesis_verdict.json`: `scientific-hypothesis-pipeline`
- `results/id_mapping.csv`: `scientific-biothings-idmapping`
- `results/imaging_report.json`: `scientific-medical-imaging`
- `results/imaging_report.md`: `scientific-medical-imaging`
- `results/impc_phenotypes.csv`: `scientific-rare-disease-genetics`
- `results/intact_interactions.csv`: `scientific-protein-interaction-network`
- `results/integrated_target_profile.json`: `scientific-pharmacology-targets`
- `results/integrated.h5ad`: `scientific-scvi-integration`
- `results/interpro_search.csv`: `scientific-protein-domain-family`
- `results/interproscan_results.csv`: `scientific-protein-domain-family`
- `results/iso13485_checklist.json`: `scientific-regulatory-science`
- `results/kegg_pathways.csv`: `scientific-pathway-enrichment`
- `results/knowledge_graph.json`: `scientific-text-mining-nlp`
- `results/ligands.csv`: `scientific-rcsb-pdb-search`
- `results/ligrec_results.json`: `scientific-spatial-transcriptomics`
- `results/lipid_annotations.csv`: `scientific-lipidomics`
- `results/lipid_da.csv`: `scientific-lipidomics`
- `results/PSP_process_property_corr.csv`: `scientific-eda-correlation`
- `results/PSP_process_structure_corr.csv`: `scientific-eda-correlation`
- `results/PSP_structure_property_corr.csv`: `scientific-eda-correlation`

</details>

<details><summary>results/ パス M–Z（226 件）</summary>

- `results/manhattan_data.csv`: `scientific-biobank-cohort`
- `results/marker_correlations.csv`: `scientific-nci60-screening`
- `results/markers.csv`: `scientific-gpu-singlecell`
- `results/masks`: `scientific-advanced-imaging`
- `results/mass_id_results.csv`: `scientific-metabolomics-databases`
- `results/materials_query.csv`: `scientific-computational-materials`
- `results/md_summary.json`: `scientific-md-simulation`
- `results/meta_analysis_summary.csv`: `scientific-meta-analysis`
- `results/metabolic_network.graphml`: `scientific-metabolic-atlas`
- `results/metabolights_study.json`: `scientific-ebi-databases`
- `results/metabolite_network.graphml`: `scientific-metabolomics-network`
- `results/metabolites.csv`: `scientific-metabolic-atlas`
- `results/metacyc_pathways.json`: `scientific-metabolomics-databases`
- `results/methylation`: `scientific-epigenomics-chromatin`
- `results/methylation/dmr_results.csv`: `scientific-epigenomics-chromatin`
- `results/mgi_phenotypes.csv`: `scientific-model-organism-db`
- `results/mgnify_taxonomy.csv`: `scientific-rrna-taxonomy`
- `results/mhc_binding_predictions.csv`: `scientific-immunoinformatics`
- `results/mid_corrected.csv`: `scientific-metabolic-flux`
- `results/mlst_typing.json`: `scientific-infectious-disease`
- `results/model`: `scientific-advanced-imaging`, `scientific-deep-chemistry`
- `results/model_metrics.csv`: `scientific-ml-regression`
- `results/model_orthologs.csv`: `scientific-model-organism-db`
- `results/molecular_network.graphml`: `scientific-proteomics-mass-spectrometry`
- `results/molecular_properties.csv`: `scientific-cheminformatics`
- `results/monarch_diseases.csv`: `scientific-monarch-ontology`
- `results/monarch_genes.csv`: `scientific-monarch-ontology`
- `results/monarch_phenotypes.csv`: `scientific-monarch-ontology`
- `results/morphometric_data.csv`: `scientific-image-analysis`
- `results/motif_enrichment.csv`: `scientific-scatac-signac`
- `results/mtb_report.json`: `scientific-precision-oncology`
- `results/mtb_report.md`: `scientific-precision-oncology`
- `results/multimodal_wnn.h5mu`: `scientific-scatac-signac`
- `results/multiomics_clusters.csv`: `scientific-multi-omics`
- `results/mutation_signatures.csv`: `scientific-cancer-genomics`
- `results/mutations.csv`: `scientific-icgc-cancer-data`
- `results/mwb_compounds.csv`: `scientific-metabolomics-databases`
- `results/mychem_annotation.json`: `scientific-biothings-idmapping`
- `results/mygene_annotation.json`: `scientific-biothings-idmapping`
- `results/myvariant_annotation.json`: `scientific-biothings-idmapping`
- `results/ncbi_sequence.fasta`: `scientific-genome-sequence-tools`
- `results/ner_entities.csv`: `scientific-text-mining-nlp`
- `results/network_proximity.json`: `scientific-drug-repurposing`
- `results/neuroscience`: `scientific-neuroscience-electrophysiology`
- `results/nhanes_data.csv`: `scientific-public-health-data`
- `results/niche_composition.csv`: `scientific-squidpy-advanced`
- `results/node_attributes.csv`: `scientific-network-analysis`
- `results/oa_availability.json`: `scientific-preprint-archive`
- `results/obis_occurrences.csv`: `scientific-marine-ecology`
- `results/occurrences.csv`: `scientific-paleobiology`
- `results/ode_solutions.json`: `scientific-symbolic-mathematics`
- `results/omero_image_metadata.json`: `scientific-lab-data-management`
- `results/omim_search.csv`: `scientific-rare-disease-genetics`
- `results/ontology_hierarchy.json`: `scientific-ontology-enrichment`
- `results/openalex_results.csv`: `scientific-literature-search`
- `results/ora`: `scientific-gene-expression-transcriptomics`
- `results/ora_results.csv`: `scientific-pathway-enrichment`
- `results/ordination_scores.csv`: `scientific-environmental-ecology`
- `results/orf_predictions.csv`: `scientific-sequence-analysis`
- `results/orphanet_diseases.csv`: `scientific-rare-disease-genetics`
- `results/orthologs.csv`: `scientific-plant-biology`
- `results/ot_associations.csv`: `scientific-opentargets-genetics`
- `results/ot_drugs.csv`: `scientific-opentargets-genetics`
- `results/papers.csv`: `scientific-semantic-scholar`
- `results/parameter_estimates.json`: `scientific-systems-biology`
- `results/pareto_optimal.csv`: `scientific-process-optimization`
- `results/particle_size_stats.csv`: `scientific-image-analysis`
- `results/patent_search.csv`: `scientific-regulatory-science`
- `results/pathway_abundance.tsv`: `scientific-microbiome-metagenomics`
- `results/pathway_activity_scores.csv`: `scientific-multi-omics`
- `results/pathway_enrichment.csv`: `scientific-metabolomics`, `scientific-metabolomics-network`, `scientific-statistical-testing`
- `results/pca_eigenvec.csv`: `scientific-population-genetics`
- `results/pca_tsne_coordinates.csv`: `scientific-pca-tsne`
- `results/pdb_entries.csv`: `scientific-rcsb-pdb-search`
- `results/pdbe_entry.json`: `scientific-structural-proteomics`
- `results/peak_detection_results.csv`: `scientific-spectral-signal`
- `results/perturbation_de.json`: `scientific-perturbation-analysis`
- `results/perturbation_signatures.json`: `scientific-perturbation-analysis`
- `results/pgx_report.json`: `scientific-pharmacogenomics`, `scientific-variant-interpretation`
- `results/pharos_diseases.csv`: `scientific-pharos-targets`
- `results/pharos_ligands.csv`: `scientific-pharos-targets`
- `results/pharos_targets.csv`: `scientific-pharmacology-targets`, `scientific-pharos-targets`
- `results/phenotype_dict.csv`: `scientific-biobank-cohort`
- `results/phewas.csv`: `scientific-gwas-catalog`
- `results/phylo_diversity.json`: `scientific-phylogenetics`
- `results/phylogenetic_tree.nwk`: `scientific-phylogenetics`
- `results/pk_model.json`: `scientific-admet-pharmacokinetics`
- `results/plant_pathways.csv`: `scientific-plant-biology`
- `results/plddt_profiles.csv`: `scientific-alphafold-structures`
- `results/polypharmacology.csv`: `scientific-stitch-chemical-network`
- `results/power_analysis.csv`: `scientific-hypothesis-pipeline`
- `results/ppi_centrality.csv`: `scientific-protein-interaction-network`
- `results/ppi_communities.csv`: `scientific-string-network-api`
- `results/ppi_network.graphml`: `scientific-protein-interaction-network`
- `results/ppi_topology.csv`: `scientific-string-network-api`
- `results/predictions.csv`: `scientific-alphafold-structures`, `scientific-deep-chemistry`
- `results/preprint_search.csv`: `scientific-preprint-archive`
- `results/projects.csv`: `scientific-icgc-cancer-data`
- `results/protein_domains.csv`: `scientific-protein-domain-family`
- `results/protein_features.csv`: `scientific-uniprot-proteome`
- `results/protein_features.json`: `scientific-structural-proteomics`
- `results/protein_quant.csv`: `scientific-proteomics-mass-spectrometry`
- `results/protocol.json`: `scientific-lab-data-management`
- `results/psm_results.csv`: `scientific-proteomics-mass-spectrometry`
- `results/ptm_sites.csv`: `scientific-proteomics-mass-spectrometry`
- `results/publication_bias_tests.csv`: `scientific-meta-analysis`
- `results/pubmed_search.csv`: `scientific-literature-search`
- `results/pubtator_annotations.csv`: `scientific-biomedical-pubtator`
- `results/pv_signal_report.json`: `scientific-pharmacovigilance`
- `results/pv_signal_report.md`: `scientific-pharmacovigilance`
- `results/qc_report.json`: `scientific-lab-automation`
- `results/quantum_result.json`: `scientific-quantum-computing`
- `results/radiomics_features.json`: `scientific-medical-imaging`
- `results/rare_disease_profile.json`: `scientific-rare-disease-genetics`
- `results/reactions.csv`: `scientific-metabolic-atlas`
- `results/reactome_enrichment.csv`: `scientific-pathway-enrichment`
- `results/reactome_participants.csv`: `scientific-reactome-pathways`
- `results/reactome_pathways.csv`: `scientific-reactome-pathways`
- `results/references.csv`: `scientific-semantic-scholar`
- `results/refmet_standardized.csv`: `scientific-metabolomics-databases`
- `results/refs`: `scientific-rrna-taxonomy`
- `results/regulatory_features.csv`: `scientific-ensembl-genomics`
- `results/regulome_scores.csv`: `scientific-regulatory-genomics`
- `results/relations.csv`: `scientific-text-mining-nlp`
- `results/remap_binding.csv`: `scientific-regulatory-genomics`
- `results/repertoire_diversity.json`: `scientific-immunoinformatics`
- `results/repurposing_candidates.json`: `scientific-drug-repurposing`
- `results/repurposing_report.md`: `scientific-drug-repurposing`
- `results/rfam_cmscan_hits.json`: `scientific-noncoding-rna`
- `results/rfam_family.json`: `scientific-noncoding-rna`
- `results/rfam_structures.json`: `scientific-noncoding-rna`
- `results/risk_measures.json`: `scientific-epidemiology-public-health`
- `results/risk_of_bias.csv`: `scientific-systematic-review`
- `results/rl_optimization.json`: `scientific-reinforcement-learning`
- `results/rl_training_log.json`: `scientific-reinforcement-learning`
- `results/rmsd_timeseries.csv`: `scientific-md-simulation`
- `results/rmsf_per_residue.csv`: `scientific-md-simulation`
- `results/rnacentral_search.csv`: `scientific-noncoding-rna`
- `results/rscu_analysis.csv`: `scientific-sequence-analysis`
- `results/safety_analysis.csv`: `scientific-survival-clinical`
- `results/samples.csv`: `scientific-geo-expression`
- `results/sar_summary.json`: `scientific-chembl-assay-mining`
- `results/sc_cellchat_interactions.csv`: `scientific-single-cell-genomics`
- `results/sc_celltype_annotations.json`: `scientific-single-cell-genomics`
- `results/sc_deg_results.csv`: `scientific-single-cell-genomics`
- `results/sc_qc_summary.json`: `scientific-single-cell-genomics`
- `results/sc_velocity_summary.json`: `scientific-single-cell-genomics`
- `results/scatac_clustered.h5ad`: `scientific-scatac-signac`
- `results/scib_benchmark.json`: `scientific-perturbation-analysis`
- `results/screen_ccres.csv`: `scientific-encode-screen`
- `results/screening_records.csv`: `scientific-systematic-review`
- `results/scvi_model`: `scientific-scvi-integration`
- `results/sdm_predictions.tif`: `scientific-environmental-ecology`
- `results/sdrf.csv`: `scientific-arrayexpress-expression`
- `results/search_strategy.json`: `scientific-systematic-review`
- `results/selection_scan.csv`: `scientific-population-genetics`
- `results/selectivity_profile.csv`: `scientific-chembl-assay-mining`
- `results/semantic_scholar_results.csv`: `scientific-literature-search`
- `results/sensitivity_analysis.csv`: `scientific-causal-inference`, `scientific-systems-biology`
- `results/sequence_composition.csv`: `scientific-sequence-analysis`
- `results/signature_exposures.csv`: `scientific-cancer-genomics`
- `results/simulation_timecourse.csv`: `scientific-systems-biology`
- `results/snp_matrix.csv`: `scientific-infectious-disease`
- `results/spatial_autocorrelation.csv`: `scientific-squidpy-advanced`
- `results/spatial_clusters.geojson`: `scientific-epidemiology-public-health`
- `results/spatial_domains.json`: `scientific-spatial-transcriptomics`
- `results/spectral_similarity.csv`: `scientific-spectral-signal`
- `results/spike_sorting`: `scientific-neuroscience-electrophysiology`
- `results/spike_sorting/quality_metrics.csv`: `scientific-neuroscience-electrophysiology`
- `results/spike_sorting/sorting_results.npz`: `scientific-neuroscience-electrophysiology`
- `results/spliceai_scores.csv`: `scientific-variant-effect-prediction`
- `results/star_allele_calls.csv`: `scientific-pharmacogenomics`
- `results/stationarity_test.csv`: `scientific-time-series`
- `results/statistical_tests.csv`: `scientific-statistical-testing`
- `results/stitch_interactions.csv`: `scientific-protein-interaction-network`, `scientific-stitch-chemical-network`
- `results/stitch_network.csv`: `scientific-stitch-chemical-network`
- `results/str_verification.json`: `scientific-cell-line-resources`
- `results/string_enrichment.csv`: `scientific-string-network-api`
- `results/string_interactions.csv`: `scientific-protein-interaction-network`
- `results/string_network.csv`: `scientific-string-network-api`
- `results/structural_alerts.csv`: `scientific-cheminformatics`
- `results/structural_alerts.json`: `scientific-chembl-assay-mining`
- `results/structure_analysis.json`: `scientific-protein-structure-analysis`
- `results/structure_report.md`: `scientific-protein-structure-analysis`
- `results/structure_zone_statistics.csv`: `scientific-materials-characterization`
- `results/structure.cif`: `scientific-computational-materials`
- `results/subgroup_analysis.csv`: `scientific-meta-analysis`
- `results/svg_results.csv`: `scientific-spatial-transcriptomics`
- `results/symbolic_solutions.json`: `scientific-symbolic-mathematics`
- `results/tair_genes.csv`: `scientific-plant-biology`
- `results/tanimoto_similarity.csv`: `scientific-cheminformatics`
- `results/target_profile_report.md`: `scientific-drug-target-profiling`
- `results/target_profile.json`: `scientific-drug-target-profiling`
- `results/taxa.csv`: `scientific-paleobiology`
- `results/taxonomy.csv`: `scientific-microbiome-metagenomics`, `scientific-rrna-taxonomy`
- `results/test_results.json`: `scientific-hypothesis-pipeline`
- `results/tissue_patterns.csv`: `scientific-nci60-screening`
- `results/topic_model_info.csv`: `scientific-text-mining-nlp`
- `results/tox21_assays.csv`: `scientific-toxicology-env`
- `results/toxicity_pathways.json`: `scientific-toxicology-env`
- `results/transformer_ft`: `scientific-deep-learning`
- `results/transmission_network.json`: `scientific-infectious-disease`
- `results/trial_adverse_events.csv`: `scientific-clinical-trials-analytics`
- `results/trial_details.json`: `scientific-clinical-trials-analytics`
- `results/trial_matches.json`: `scientific-clinical-decision-support`
- `results/umls_mapping.json`: `scientific-ontology-enrichment`
- `results/uniprot_entries.csv`: `scientific-uniprot-proteome`
- `results/univariate_results.csv`: `scientific-metabolomics`
- `results/vaccine_candidates_ranked.csv`: `scientific-immunoinformatics`
- `results/variant_actionability.json`: `scientific-precision-oncology`
- `results/variant_classification.json`: `scientific-variant-interpretation`
- `results/variant_report.md`: `scientific-variant-interpretation`
- `results/vep_consequences.csv`: `scientific-ensembl-genomics`
- `results/vip_scores.csv`: `scientific-metabolomics`
- `results/virtual_screening.csv`: `scientific-molecular-docking`
- `results/vs_library.csv`: `scientific-compound-screening`
- `results/williamson_hall.csv`: `scientific-materials-characterization`
- `results/works.csv`: `scientific-crossref-metadata`
- `results/worms_taxonomy.csv`: `scientific-marine-ecology`
- `results/xai_report.json`: `scientific-explainable-ai`
- `results/xrd_analysis.csv`: `scientific-materials-characterization`
- `results/yearly_trend.csv`: `scientific-semantic-scholar`
- `results/zinc_catalogs.csv`: `scientific-compound-screening`
- `results/zinc_search.csv`: `scientific-compound-screening`
- `results/zinc_similar.csv`: `scientific-compound-screening`
- `results/zinc_substance.json`: `scientific-compound-screening`

</details>

# 参照関係（逆引き）

各スキルの SKILL.md 本文から他スキルへのリンクを解析し、**「どのスキルから参照されているか」** を逆引きできるようにしたものです。

> 凡例: `スキルA` ← `スキルB`, `スキルC` … スキルA はスキルB・スキルC の SKILL.md から参照されている

- `scientific-bayesian-statistics` ← `scientific-epidemiology-public-health`, `scientific-infectious-disease`, `scientific-spatial-transcriptomics`, `scientific-systems-biology`
- `scientific-bioinformatics` ← `scientific-infectious-disease`, `scientific-population-genetics`, `scientific-single-cell-genomics`, `scientific-spatial-transcriptomics`
- `scientific-causal-inference` ← `scientific-epidemiology-public-health`, `scientific-microbiome-metagenomics`
- `scientific-citation-checker` ← `scientific-text-mining-nlp`
- `scientific-clinical-trials-analytics` ← `scientific-epidemiology-public-health`
- `scientific-deep-learning` ← `scientific-single-cell-genomics`
- `scientific-deep-research` ← `scientific-text-mining-nlp`
- `scientific-disease-research` ← `scientific-population-genetics`
- `scientific-doe` ← `scientific-systems-biology`
- `scientific-epigenomics-chromatin` ← `scientific-population-genetics`, `scientific-single-cell-genomics`
- `scientific-gene-expression-transcriptomics` ← `scientific-single-cell-genomics`
- `scientific-graph-neural-networks` ← `scientific-text-mining-nlp`
- `scientific-image-analysis` ← `scientific-environmental-ecology`, `scientific-spatial-transcriptomics`
- `scientific-infectious-disease` ← `scientific-epidemiology-public-health`
- `scientific-meta-analysis` ← `scientific-epidemiology-public-health`, `scientific-text-mining-nlp`
- `scientific-metabolomics` ← `scientific-microbiome-metagenomics`, `scientific-systems-biology`
- `scientific-ml-classification` ← `scientific-environmental-ecology`
- `scientific-multi-omics` ← `scientific-microbiome-metagenomics`, `scientific-single-cell-genomics`, `scientific-systems-biology`
- `scientific-network-analysis` ← `scientific-infectious-disease`, `scientific-microbiome-metagenomics`, `scientific-single-cell-genomics`, `scientific-spatial-transcriptomics`, `scientific-systems-biology`, `scientific-text-mining-nlp`
- `scientific-pca-tsne` ← `scientific-environmental-ecology`, `scientific-population-genetics`, `scientific-single-cell-genomics`
- `scientific-pharmacogenomics` ← `scientific-population-genetics`
- `scientific-protein-design` ← `scientific-immunoinformatics`
- `scientific-protein-structure-analysis` ← `scientific-immunoinformatics`
- `scientific-sequence-analysis` ← `scientific-immunoinformatics`, `scientific-infectious-disease`
- `scientific-single-cell-genomics` ← `scientific-immunoinformatics`, `scientific-spatial-transcriptomics`
- `scientific-statistical-testing` ← `scientific-environmental-ecology`, `scientific-microbiome-metagenomics`, `scientific-population-genetics`
- `scientific-survival-clinical` ← `scientific-epidemiology-public-health`, `scientific-infectious-disease`
- `scientific-time-series` ← `scientific-environmental-ecology`
- `scientific-variant-interpretation` ← `scientific-immunoinformatics`, `scientific-population-genetics`

# トリガーフレーズ → スキル

各スキルの frontmatter `description` に含まれる全角カギ括弧「…」を抽出したものです。Claude へのプロンプト文言からスキルを探す用途に使えます。

<details><summary>英語フレーズ（27 件）</summary>

- 「Abstract を作成して」: `scientific-academic-writing`
- 「ADMET 予測して」: `scientific-admet-pharmacokinetics`
- 「citation check」: `scientific-citation-checker`
- 「DICOM を解析して」: `scientific-medical-imaging`
- 「diff を出して」: `scientific-revision-tracker`
- 「druggability 分析して」: `scientific-drug-target-profiling`
- 「ESM で評価して」: `scientific-protein-design`
- 「FAERS データを解析して」: `scientific-pharmacovigilance`
- 「GNN で分子特性を予測して」: `scientific-graph-neural-networks`
- 「GWAS 結果を解析して」: `scientific-disease-research`
- 「LaTeX に変換して」: `scientific-latex-formatter`
- 「lead optimization して」: `scientific-admet-pharmacokinetics`
- 「LIME で説明して」: `scientific-explainable-ai`
- 「MCMC で推定して」: `scientific-bayesian-statistics`
- 「Methods セクションを書いて」: `scientific-academic-writing`
- 「OncoKB で検索して」: `scientific-precision-oncology`
- 「PDB 構造を調べて」: `scientific-protein-structure-analysis`
- 「pharmacogenomics 解析して」: `scientific-variant-interpretation`
- 「reviewer response を作成して」: `scientific-peer-review-response`
- 「SHAP 値を計算して」: `scientific-explainable-ai`
- 「SIを作って」: `scientific-supplementary-generator`
- 「Specific Aims を書いて」: `scientific-grant-writing`
- 「Supplementary を作成して」: `scientific-supplementary-generator`
- 「systematic review して」: `scientific-deep-research`
- 「Transformer を Fine-tune して」: `scientific-deep-learning`
- 「VQE で基底エネルギーを求めて」: `scientific-quantum-computing`
- 「WSI を処理して」: `scientific-medical-imaging`

</details>

<details><summary>日本語フレーズ（59 件）</summary>

- 「がんゲノム解析して」: `scientific-precision-oncology`
- 「グラフ分類して」: `scientific-graph-neural-networks`
- 「グラント申請書を書いて」: `scientific-grant-writing`
- 「このデータで何がわかる？」: `scientific-hypothesis-pipeline`
- 「ジャーナルフォーマットにして」: `scientific-latex-formatter`
- 「ターゲット評価して」: `scientific-drug-target-profiling`
- 「タンパク質の構造を解析して」: `scientific-protein-structure-analysis`
- 「タンパク質を設計して」: `scientific-protein-design`
- 「ドッキング準備して」: `scientific-protein-structure-analysis`
- 「ドラッグリポジショニングして」: `scientific-drug-repurposing`
- 「ナレッジグラフ推論して」: `scientific-graph-neural-networks`
- 「ニューラルネットで学習して」: `scientific-deep-learning`
- 「バリアントの病原性を評価して」: `scientific-variant-interpretation`
- 「ブレインストーミングして」: `scientific-research-methodology`
- 「ベイズ回帰して」: `scientific-bayesian-statistics`
- 「ポスターのレイアウトを設計して」: `scientific-presentation-design`
- 「モデルの予測を説明して」: `scientific-explainable-ai`
- 「リバッタルを書いて」: `scientific-peer-review-response`
- 「事後分布を求めて」: `scientific-bayesian-statistics`
- 「仮説を立てて」: `scientific-hypothesis-pipeline`
- 「先行研究を調べて」: `scientific-deep-research`
- 「医用画像をセグメンテーションして」: `scientific-medical-imaging`
- 「参考文献を検索」: `scientific-citation-checker`
- 「可読性スコアを出して」: `scientific-paper-quality`
- 「変更履歴を作って」: `scientific-revision-tracker`
- 「学会スライドを作成して」: `scientific-presentation-design`
- 「安全性シグナルを検出して」: `scientific-pharmacovigilance`
- 「安定性を予測して」: `scientific-protein-design`
- 「実験プロトコルを作成して」: `scientific-lab-automation`
- 「希少疾患を診断して」: `scientific-disease-research`
- 「引用をチェックして」: `scientific-citation-checker`
- 「投稿前チェック」: `scientific-paper-quality`
- 「投稿用 TeX を作って」: `scientific-latex-formatter`
- 「改訂をトラッキングして」: `scientific-revision-tracker`
- 「文献調査して」: `scientific-deep-research`
- 「既存薬の新規適応を探して」: `scientific-drug-repurposing`
- 「有害事象を分析して」: `scientific-pharmacovigilance`
- 「査読に回答して」: `scientific-peer-review-response`
- 「標的タンパク質を調べて」: `scientific-drug-target-profiling`
- 「治療推奨を作成して」: `scientific-clinical-decision-support`
- 「液体ハンドリングを自動化して」: `scientific-lab-automation`
- 「深層学習モデルを構築して」: `scientific-deep-learning`
- 「疾患と遺伝子の関連を調べて」: `scientific-disease-research`
- 「研究デザインを設計して」: `scientific-research-methodology`
- 「研究計画を立てて」: `scientific-research-methodology`
- 「科研費を作成して」: `scientific-grant-writing`
- 「精密医療の解析して」: `scientific-clinical-decision-support`
- 「考察を深めて」: `scientific-critical-review`
- 「腫瘍プロファイリングして」: `scientific-precision-oncology`
- 「臨床パスウェイを設計して」: `scientific-clinical-decision-support`
- 「草稿を改善して」: `scientific-critical-review`
- 「薬物動態を評価して」: `scientific-admet-pharmacokinetics`
- 「補足資料を生成して」: `scientific-supplementary-generator`
- 「解析パイプラインを作って」: `scientific-hypothesis-pipeline`
- 「論文の品質をチェックして」: `scientific-paper-quality`
- 「論文をレビューして」: `scientific-critical-review`
- 「論文を書いて」: `scientific-academic-writing`
- 「量子シミュレーションして」: `scientific-quantum-computing`
- 「量子回路を設計して」: `scientific-quantum-computing`

</details>

