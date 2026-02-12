# SATORI v0.11.0 スキル拡張ギャップ分析

> 作成日: 2025-07  
> 対象バージョン: v0.11.0 (66 スキル / 26 カテゴリ A–Z)  
> 前回分析: GAP_ANALYSIS_v0.10.md (56 → 66 へ Phase 2 実装済み)

---

## 1. 調査概要

### 1.1 比較対象リポジトリ

| リポジトリ | 規模 | 最新特徴 |
|-----------|------|---------|
| **mims-harvard/ToolUniverse** | 350+ SMCP ツール, 90+ カテゴリ | `default_config.py` に 90+ カテゴリ登録。chipatlas (エピゲノミクス), fourdn (3D ゲノム), regulomedb/jaspar/remap/screen (調節ゲノミクス), pride (プロテオミクス), pharmgkb (ファーマコゲノミクス), spliceai (スプライシング), expression_atlas (遺伝子発現), geo (GEO), guidelines (ガイドライン) 等の専門 DB カテゴリが充実 |
| **K-Dense-AI/claude-scientific-skills** | 140+ スキル, ~17 カテゴリ | pyOpenMS/matchms (プロテオミクス), neuropixels-analysis/neurokit2 (神経科学), pymatgen (計算材料科学), clinpgx-database (ファーマコゲノミクス), fda-database (規制科学), geo-database/pydeseq2 (トランスクリプトミクス), scientific-schematics (科学図版), clinicaltrials-database (臨床試験), iso-13485-certification/uspto-database (規制・知財) |

### 1.2 現在の SATORI カバレッジ

- **66 スキル** / **26 カテゴリ (A–Z)** — 全カテゴリ使用済み
- 新スキルは既存カテゴリの**拡張**として追加
- Phase 1 (GAP_ANALYSIS.md): 47 → 56 (9 スキル追加)
- Phase 2 (GAP_ANALYSIS_v0.10.md): 56 → 66 (10 スキル追加)

### 1.3 調査対象ドメイン (10 領域)

Phase 2 完了後も残存するギャップ領域として以下を特定:

1. エピゲノミクス・クロマチン生物学
2. プロテオミクス・質量分析
3. 神経科学・電気生理学
4. ファーマコゲノミクス
5. 臨床試験解析
6. 規制科学・コンプライアンス
7. 遺伝子発現・トランスクリプトミクス
8. 計算材料科学
9. 科学図版生成 (スキマティクス)
10. ラボデータ管理・クラウド基盤

---

## 2. ドメイン別ギャップマッピング

### 凡例
- ◎ = 強力なツール/スキル支援あり (5+ ツール)
- ○ = 中程度の支援あり (2–4 ツール)
- △ = 限定的な支援
- ✗ = 支援なし

| # | ドメイン | ToolUniverse | K-Dense | SATORI 現状 | ギャップ度 |
|---|---------|-------------|---------|------------|----------|
| 1 | エピゲノミクス・クロマチン | ◎ chipatlas (4), fourdn, regulomedb, jaspar, remap, screen (20+ツール) | ○ geniml, deepTools | ✗ 未カバー | **Critical** |
| 2 | プロテオミクス・質量分析 | ○ pride | ◎ matchms, pyOpenMS | ✗ 未カバー | **Critical** |
| 3 | 神経科学・電気生理学 | ✗ | ◎ neuropixels-analysis, neurokit2 (MNE統合) | △ biosignal-processing が部分的 | **High** |
| 4 | ファーマコゲノミクス | ○ pharmgkb, fda_pharmacogenomic | ◎ clinpgx-database (9+機能) | ✗ 未カバー | **High** |
| 5 | 臨床試験解析 | ◎ clinical_trials (15+ツール), ClinicalTrialDesignAgent | ◎ clinicaltrials-database | △ survival-clinical が部分的 | **High** |
| 6 | 規制科学・コンプライアンス | ○ fda_drug_label, adverse_events, guidelines | ◎ fda-database, iso-13485, uspto-database | ✗ 未カバー | **Moderate-High** |
| 7 | 遺伝子発現・トランスクリプトミクス | ◎ geo (3), gtex_v2, expression_atlas | ◎ geo-database, pydeseq2 | △ bioinformatics が部分的 | **Moderate-High** |
| 8 | 計算材料科学 | ✗ | ◎ pymatgen (結晶構造, Materials Project) | △ materials-characterization が部分的 | **Moderate** |
| 9 | 科学図版生成 | ✗ | ◎ scientific-schematics (CONSORT, TikZ/SVG) | △ publication-figures が部分的 | **Moderate** |
| 10 | ラボデータ管理・クラウド基盤 | ✗ | ◎ benchling, dnanexus, latchbio, omero, protocolsio | △ lab-automation が部分的 | **Moderate** |

---

## 3. 新スキル提案 (10 スキル)

### カテゴリ別配置

| カテゴリ | 名称 | 追加スキル数 | 追加後合計 |
|---------|------|------------|-----------|
| **E** | 信号処理 | +1 | 4 |
| **F** | 生命科学・オミクス | +2 | 7 |
| **G** | 化学・材料・画像 | +1 | 4 |
| **H** | 臨床・疫学統計 | +1 | 4 |
| **M** | ラボ・実験 | +1 | 2 |
| **N** | プレゼンテーション・図版 | +1 | 2 |
| **O** | 研究支援・方法論 | +1 | 3 |
| **P** | 薬理学 | +1 | 2 |
| **T** | シングルセル・ゲノミクス | +1 | 3 |

---

### T-3. `scientific-epigenomics-chromatin`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | T. シングルセル・空間オミクス → ゲノミクス全般に拡張 |
| **優先度** | **HIGH** |
| **説明** | エピゲノミクス解析パイプライン。ChIP-seq、ATAC-seq、バイサルファイトシーケンシング (DNA メチル化) データの解析、ヒストン修飾マッピング、クロマチンアクセシビリティ解析、転写因子結合部位予測、3D ゲノム構造解析を統合。ChIP-Atlas (43万+実験) を中心バックエンドとし、計算エピゲノミクスの包括的ワークフローを提供。 |
| **既存との区別** | `single-cell-genomics` は scRNA-seq。`bioinformatics` は配列解析一般。エピゲノミクス特有のピーク呼び出し (MACS2)、差次結合解析、クロマチン状態モデリング (ChromHMM)、TAD/ループ検出は完全に未カバー。 |
| **主要テクニック** | ChIP-seq ピーク呼び出し (MACS2/MACS3)、差次結合解析 (DiffBind)、転写因子フットプリンティング、ヒストン修飾クロマチン状態 (ChromHMM/Segway)、ATAC-seq ヌクレオソームフリー領域検出、DNA メチル化パターン解析 (WGBS/RRBS)、TAD 検出・Hi-C 接触マップ解析、エンハンサー-プロモーター相互作用予測、コンセンサスピーク解析、モチーフ濃縮解析 |
| **ToolUniverse SMCP ツール** | |

```
# ChIP-Atlas (433,000+ 実験) — 4 ツール
ChIPAtlas_enrichment_analysis      # ChIP-seq 濃縮解析
ChIPAtlas_get_experiments          # 実験メタデータ取得
ChIPAtlas_get_peak_data            # ピークデータ取得
ChIPAtlas_search_datasets          # データセット検索

# 4DN Data Portal (3D ゲノム) — 2+ ツール
fourdn_search                      # 3D ゲノムデータ検索
fourdn_get_experiment              # Hi-C/ChIA-PET 実験データ

# 調節ゲノミクスツール群
regulomedb_search                  # 調節バリアントスコアリング
jaspar_search_profiles             # 転写因子結合モチーフ検索
jaspar_get_matrix                  # PWM (Position Weight Matrix) 取得
remap_search                       # 調節領域マッピング
screen_search_ccres                # cCRE (候補シス調節エレメント) 検索
screen_get_ccre_details            # cCRE 詳細取得
```

| **K-Dense 対応スキル** | geniml (コンセンサスピーク生成、ゲノム領域埋め込み、Region2Vec)、deepTools (ChIP-seq QC、correlation heatmap、fingerprint、computeMatrix/plotHeatmap) |
| **正当性** | ToolUniverse に **20+ の専用エピゲノミクスツール**が存在するにもかかわらず、SATORI は完全にゼロカバー。エピゲノミクスは ENCODE/Roadmap Epigenomics プロジェクトで巨大データが公開済みであり、がん、発生生物学、老化研究の基盤技術。ChIP-Atlas 一つで 43 万実験にアクセス可能な点は即座の価値。 |

---

### F-6. `scientific-proteomics-mass-spectrometry`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | F. 生命科学・オミクス |
| **優先度** | **HIGH** |
| **説明** | プロテオミクスおよび質量分析データ解析パイプライン。LC-MS/MS ワークフロー、ペプチド同定・定量、蛋白質定量 (LFQ/TMT/SILAC)、翻訳後修飾 (PTM) 解析、メタボロミクス質量分析、スペクトルライブラリ検索、化合物同定を統合。PRIDE データベース経由の公開プロテオミクスデータアクセスも支援。 |
| **既存との区別** | `spectral-signal` は一般的なスペクトル処理 (FFT 等)。`metabolomics` は代謝物プロファイリング。MS 特有のペプチド同定パイプライン (検索エンジン, FDR 制御, PSM) は未カバー。 |
| **主要テクニック** | LC-MS/MS データ前処理 (ピーク検出、ノイズ除去、ベースライン補正)、データベース検索 (PSM, FDR 制御)、de novo シーケンシング、蛋白質定量 (LFQ, iBAQ, TMT/iTRAQ, SILAC)、翻訳後修飾 (リン酸化, ユビキチン化, グリコシル化) マッピング、スペクトル類似度スコアリング (コサイン、修正コサイン)、分子ネットワーキング (GNPS)、化合物アノテーション (HMDB, MassBank 参照) |
| **ToolUniverse SMCP ツール** | |

```
# PRIDE (Proteomics Identifications Database) — 3+ ツール
pride_search_projects              # プロテオミクスプロジェクト検索
pride_get_project_details          # プロジェクト詳細取得
pride_get_project_files            # ファイル一覧取得
```

| **K-Dense 対応スキル** | matchms (マススペクトルマッチング、コサイン類似度、スペクトルライブラリ検索、化合物同定)、pyOpenMS (LC-MS/MS 完全ワークフロー — ペプチド同定、特徴検出、蛋白質定量、メタボロミクス解析) |
| **正当性** | プロテオミクスはゲノミクス・トランスクリプトミクスに次ぐ第3の主要オミクス技術であり、SATORI に完全な空白。K-Dense に pyOpenMS (完全な MS ワークフロー) と matchms (スペクトルマッチング) の 2 つの充実スキルがあり、即座に参照可能。PRIDE DB 連携でデータアクセスも確保。 |

---

### E-4. `scientific-neuroscience-electrophysiology`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | E. 信号処理 |
| **優先度** | **HIGH** |
| **説明** | 神経科学データ解析パイプライン。細胞外電気生理学 (Neuropixels/マルチ電極アレイ) のスパイクソーティング、EEG/EMG/ECG/EDA 等の生理信号解析、自律神経系評価、脳結合性解析を統合。SpikeInterface ベースの大規模神経記録解析と NeuroKit2/MNE による非侵襲生理信号解析を包括。 |
| **既存との区別** | `biosignal-processing` は一般的な生体信号処理 (フィルタリング、特徴抽出)。Neuropixels 特有のスパイクソーティングパイプライン (Kilosort4, 品質指標, AI キュレーション) や EEG マイクロステート解析は未カバー。 |
| **主要テクニック** | Neuropixels データロード (SpikeGLX/OpenEphys)、前処理 (位相シフト補正, バンドパスフィルタ, CMR)、モーション補正 (nonrigid)、スパイクソーティング (Kilosort4/SpykingCircus2/MountainSort5)、品質指標 (SNR, ISI 違反率, 存在比率)、Allen Institute/IBL 基準キュレーション、ECG (HRV 解析, R ピーク検出)、EDA (皮膚コンダクタンス応答、交感神経活動)、EEG (マイクロステート、事象関連電位、MNE 統合)、EMG (筋活動解析) |
| **ToolUniverse SMCP ツール** | なし (K-Dense 独自の強み) |
| **K-Dense 対応スキル** | neuropixels-analysis (SpikeInterface ベース完全パイプライン — ロード→前処理→ドリフト補正→スパイクソーティング→品質指標→AI キュレーション、Neuropixels 1.0/2.0 対応)、neurokit2 (ECG/EDA/EMG/EEG マルチモーダル生理信号解析、MNE-Python 統合、自律神経系評価) |
| **正当性** | 神経科学は NIH BRAIN Initiative 等で巨額投資が続く成長分野。K-Dense に neuropixels-analysis (完全な SpikeInterface パイプライン + AI キュレーション) と neurokit2 (MNE 統合 EEG + マルチモーダル生理信号) の 2 つの充実スキルが存在。SATORI の `biosignal-processing` は信号処理の汎用テクニックであり、神経科学固有のワークフローとは次元が異なる。 |

---

### P-2. `scientific-pharmacogenomics`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | P. 薬理学 |
| **優先度** | **HIGH** |
| **説明** | ファーマコゲノミクス (薬理ゲノム学) データベース照会・解釈パイプライン。遺伝子-薬物相互作用、CPIC (Clinical Pharmacogenetics Implementation Consortium) ガイドライン、アレル機能分類、薬物ラベル注釈、バリアントレベルのファーマコゲノミクスアノテーションを統合。ClinPGx (PharmGKB/CPIC/PharmCAT 統合後継) を中心データソースとし、個別化薬物療法を支援。 |
| **既存との区別** | `pharmacovigilance` は市販後安全性 (FAERS/有害事象)。`variant-interpretation` は臨床的バリアント解釈 (ClinVar)。ファーマコゲノミクスは**遺伝子型に基づく薬物選択・用量調整**という独自ドメインで、CPIC ガイドライン、アレル機能表、投与量レコメンデーション等が別系統。 |
| **主要テクニック** | 遺伝子-薬物ペア照会、CPIC ガイドライン取得 (レベル A/B/C/D)、Star アレルアノテーション、代謝酵素表現型分類 (PM/IM/NM/RM/UM)、薬物ラベル注釈 (FDA/EMA/PMDA)、バリアントアノテーション (臨床的有意性)、パスウェイ可視化 (PK/PD パスウェイ)、多遺伝子スコア算出、バッチ遺伝子照会 |
| **ToolUniverse SMCP ツール** | |

```
# PharmGKB — 3+ ツール
PharmGKB_search_drugs              # 薬物検索 (遺伝子関連)
PharmGKB_get_gene_drug_pairs       # 遺伝子-薬物ペア取得
PharmGKB_get_guidelines            # CPIC ガイドライン取得

# FDA ファーマコゲノミクス
fda_pharmacogenomic_biomarkers     # FDA 承認ファーマコゲノミクスバイオマーカー

# OpenTargets
OpenTargets_drug_pharmacogenomics_data  # 薬物ファーマコゲノミクスデータ
```

| **K-Dense 対応スキル** | clinpgx-database (ClinPGx API — 9 機能: 遺伝子照会、薬物/化学物質照会、遺伝子-薬物ペア照会、CPIC ガイドライン、アレル/バリアント情報、バリアントアノテーション、臨床アノテーション、薬物ラベル、パスウェイ。バッチ遺伝子照会、アクショナブルペア検索ユーティリティ付き) |
| **正当性** | ファーマコゲノミクスは FDA が Drug Labeling で 300+ の遺伝子-薬物ペアを公式推奨する実用段階の技術。ClinPGx (PharmGKB 後継) は K-Dense に包括的スキルが存在し、ToolUniverse にも PharmGKB ツールが登録済み。SATORI では `pharmacovigilance` が市販後安全性をカバーするが、**処方前の遺伝子型ベース薬物選択**は完全な空白。 |

---

### H-4. `scientific-clinical-trials-analytics`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | H. 臨床・疫学統計 |
| **優先度** | **HIGH** |
| **説明** | 臨床試験レジストリ (ClinicalTrials.gov) の包括的データアクセス・解析パイプライン。試験検索 (疾患・薬剤・地域・フェーズ・ステータス)、詳細情報取得 (適格基準・アウトカム・副作用・施設)、競合ランドスケープ解析、試験デザイン支援、バルクデータ取得・CSV エクスポートを統合。 |
| **既存との区別** | `survival-clinical` は統計的生存解析 (KM 曲線, Cox 回帰)。`meta-analysis` は発表済み結果の統合。臨床試験**レジストリの網羅的検索・構造化データ抽出・競合分析**は未カバー。ClinicalTrials.gov API v2 に直接対応する専用スキルが必要。 |
| **主要テクニック** | 多基準試験検索 (疾患 × 介入 × 地域 × ステータス × フェーズ)、NCT ID ベースの詳細取得 (プロトコル、適格基準、アウトカム測定、接触先、施設)、競合ランドスケープ解析 (同一適応のフェーズ分布、スポンサー分析)、試験デザイン諮問 (ClinicalTrialDesignAgent)、ページネーション/バルクデータ取得、CSV/JSON エクスポート、有害事象抽出、試験結果アウトカム抽出 |
| **ToolUniverse SMCP ツール** | |

```
# ClinicalTrials.gov — 15+ ツール
clinical_trials_search             # 臨床試験検索
clinical_trials_get_details        # 試験詳細取得
clinical_trials_get_conditions     # 疾患条件取得
clinical_trials_get_interventions  # 介入情報取得
clinical_trials_get_eligibility    # 適格基準取得
clinical_trials_get_locations      # 施設・地域取得
clinical_trials_get_outcome_measures # アウトカム測定取得
clinical_trials_get_references     # 参考文献取得
extract_clinical_trial_adverse_events  # 有害事象抽出
extract_clinical_trial_outcomes    # 試験結果抽出

# 臨床試験デザインエージェント
ClinicalTrialDesignAgent           # AI 支援試験デザイン
```

| **K-Dense 対応スキル** | clinicaltrials-database (ClinicalTrials.gov API v2 完全対応 — 疾患・介入・地域・スポンサー・ステータス検索、詳細取得、ページネーション、CSV エクスポート、データ構造ナビゲーションガイド) |
| **正当性** | ToolUniverse に **15+ の専用臨床試験ツール** (検索、詳細取得、条件・介入・適格基準・施設・アウトカム・有害事象・references 等) + AI 駆動の ClinicalTrialDesignAgent が存在し、K-Dense にも完全な API v2 スキルがある。臨床研究者にとって「どの試験が存在するか」「競合状況はどうか」は日常的な情報ニーズであり、両リポの最も充実したツール群の一つ。 |

---

### O-3. `scientific-regulatory-science`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | O. 研究支援・方法論 |
| **優先度** | **MEDIUM-HIGH** |
| **説明** | FDA・EMA・PMDA 等の規制当局データベース横断照会、医薬品/医療機器/食品の規制データアクセス、品質管理システム (ISO 13485)、特許検索 (USPTO) を統合した規制科学パイプライン。市販後サーベイランス、510(k) 510(k) クリアランス、薬事申請情報、リコール情報の体系的取得を支援。 |
| **既存との区別** | `pharmacovigilance` は市販後安全性 (FAERS シグナル検出)。規制科学は**承認プロセス・デバイス規制・食品安全・品質システム・知的財産**まで含む横断的ドメイン。FDA 6+ エンドポイント (医薬品/デバイス/食品/動物薬) を統合的にカバー。 |
| **主要テクニック** | FDA 医薬品データ (有害事象、ラベリング、NDC、リコール、承認履歴、供給不足)、FDA 医療機器データ (有害事象、510(k)、分類、PMA、UDI、リコール)、FDA 食品データ (有害事象、リコール)、ISO 13485 品質管理システム (設計管理、リスク管理、CAPA)、USPTO 特許検索 (PEDS クライアント)、規制当局間比較分析 |
| **ToolUniverse SMCP ツール** | |

```
# FDA 医薬品ツール — 既存カテゴリ活用
fda_drug_label_*                   # 添付文書ツール群
adverse_events_*                   # 有害事象ツール群
guidelines_*                       # ガイドラインツール群
```

| **K-Dense 対応スキル** | fda-database (FDA 6 カテゴリ: 医薬品 — 有害事象/ラベリング/NDC/リコール/承認/供給不足、医療機器 — 有害事象/510(k)/分類/PMA/UDI/リコール、食品 — 有害事象/リコール、動物薬 — 有害事象。Python FDAQuery クラス付き)、iso-13485-certification (医療機器品質管理システム — 設計管理/リスク管理/CAPA/プロセス検証)、uspto-database (特許検索 — PEDS クライアントによる米国特許検索) |
| **正当性** | K-Dense に FDA DB 完全スキル + ISO 13485 + USPTO の 3 本柱があり、合計で**医薬品・医療機器・食品・知的財産をカバーする包括的規制科学スキル**が構成可能。ToolUniverse にも FDA ツールとガイドラインカテゴリが存在。規制科学的観点は論文だけでなく**実用化・事業化**に直結する実務ニーズ。 |

---

### F-7. `scientific-gene-expression-transcriptomics`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | F. 生命科学・オミクス |
| **優先度** | **MEDIUM-HIGH** |
| **説明** | バルク RNA-seq / マイクロアレイの遺伝子発現データ解析パイプライン。GEO (Gene Expression Omnibus) からの公開データセット取得・前処理、差次発現解析 (DESeq2)、GTEx 組織発現参照、Expression Atlas (EBI GXA) ベースライン/差次発現統合照会を支援。 |
| **既存との区別** | `bioinformatics` は配列解析 (FASTA/VCF/BAM)。`single-cell-genomics` は scRNA-seq。**バルクトランスクリプトミクス** — GEO データセット検索・取得、DESeq2 差次発現解析、GTEx 組織発現 — は独立したワークフロー。 |
| **主要テクニック** | GEO データセット検索・メタデータ取得 (GEOparse)、発現マトリクス抽出/前処理、サンプルメタデータフィルタリング、差次発現解析 (PyDESeq2 — GLM フィッティング、Wald 検定、LFC 収縮、FDR 補正)、Volcano/MA プロット生成、GTEx 組織発現プロファイル照会、eQTL 解析、Expression Atlas ベースライン/差次発現照会、GSEA (遺伝子セット濃縮解析)、バッチ効果補正 |
| **ToolUniverse SMCP ツール** | |

```
# GEO (Gene Expression Omnibus) — 3 ツール
geo_search_datasets                # GEO データセット検索
geo_get_dataset_info               # データセット情報取得
geo_get_sample_info                # サンプル情報取得

# GTEx v2 — 3+ ツール
gtex_search_gene_expression        # 組織別遺伝子発現検索
gtex_get_eqtl                      # eQTL データ取得
gtex_get_tissue_expression         # 組織発現プロファイル

# Expression Atlas (EBI GXA)
expression_atlas_search            # 発現解析結果検索
expression_atlas_get_experiment    # 実験詳細取得
```

| **K-Dense 対応スキル** | geo-database (GEOparse ベース完全スキル — データセット検索、発現データ取得、バッチ処理、差次発現解析)、pydeseq2 (完全 DESeq2 ワークフロー — カウントデータ入力、正規化、統計検定、LFC 収縮、Volcano/MA プロット、マルチファクターデザイン対応) |
| **正当性** | ToolUniverse に GEO (3 ツール) + GTEx (3+ ツール) + Expression Atlas が存在し、K-Dense にも GEO + PyDESeq2 の充実スキルがある。バルクトランスクリプトミクスは**最も普及した実験手法**であり、SATORI では `bioinformatics` で一般的な配列操作をカバーするのみ。公開 GEO データセット (400 万+ サンプル) への体系的アクセスは高頻度の研究ニーズ。 |

---

### G-4. `scientific-computational-materials`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | G. 化学・材料・画像 |
| **優先度** | **MEDIUM** |
| **説明** | 計算材料科学パイプライン。結晶構造解析・操作、相図計算、電子構造計算、Materials Project API 経由のデータベース照会、VASP/Quantum ESPRESSO 入出力操作、ハイスループットスクリーニングワークフローを支援。Pymatgen ベースの包括的材料シミュレーション。 |
| **既存との区別** | `materials-characterization` は分析化学 (XRD パターン解析、IR/Raman スペクトル解釈)。計算材料科学は**第一原理計算 (DFT)、相図計算 (CALPHAD)、結晶構造設計**というシミュレーション駆動のワークフロー。 |
| **主要テクニック** | 結晶構造作成・操作 (スペースグループ, 格子定数, 原子位置)、欠陥/表面スラブ/ヘテロ構造生成、相図計算 (凸包安定性解析)、電子バンド構造/DOS 計算、Materials Project API 照会 (材料検索, 熱力学データ)、VASP/QE 入力ファイル生成・結果解析、ハイスループットスクリーニング (候補材料フィルタリング)、対称性解析 (spglib 統合)、弾性定数/誘電定数/磁性特性計算 |
| **ToolUniverse SMCP ツール** | なし (K-Dense 独自の強み) |
| **K-Dense 対応スキル** | pymatgen (包括的計算材料科学スキル — 結晶構造, 相図, 電子構造, Materials Project API, VASP 統合, ハイスループットスクリーニング。分子カスタム構造体パターン、GW近似、トポロジカル解析等を含む高度な機能群) |
| **正当性** | 計算材料科学は電池、触媒、半導体、超伝導体等の材料開発に不可欠。K-Dense に Pymatgen の包括的スキルが存在し、Materials Project (15 万+ 材料) へのAPI アクセスも統合済み。SATORI の `materials-characterization` は実験データ解釈であり、**in silico 材料設計**は完全に未カバー。 |

---

### N-2. `scientific-scientific-schematics`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | N. プレゼンテーション・図版 |
| **優先度** | **MEDIUM** |
| **説明** | 科学論文用スキマティック図版 (概念図・フローチャート・アーキテクチャ図) の生成パイプライン。CONSORT フロー図、神経ネットワークアーキテクチャ図、生物学的パスウェイ図、システム構成図等を TikZ/SVG/PNG で生成。AI レビュー付きイテレーション機能を統合。 |
| **既存との区別** | `publication-figures` は matplotlib/seaborn によるデータ可視化 (散布図、ヒートマップ等)。`presentation-design` はスライドレイアウト。**概念図・フローチャート・アーキテクチャ図**はデータプロットとは異なるカテゴリ。 |
| **主要テクニック** | CONSORT 参加者フロー図生成、臨床試験デザインフロー、神経ネットワーク層構成図、生物学的パスウェイダイアグラム、システムアーキテクチャ図、TikZ コード生成 (LaTeX 統合)、SVG/PNG エクスポート、プロンプト駆動図版生成、AI レビューによる反復改善 (3 イテレーション)、テンプレートベースのカスタマイズ |
| **ToolUniverse SMCP ツール** | なし (K-Dense 独自の強み) |
| **K-Dense 対応スキル** | scientific-schematics (完全スキル — CONSORT フローチャート、神経ネットワーク図、生物学的パスウェイ図、システムアーキテクチャ図。TikZ/SVG 生成、AI 生成モード (FLUX.2 Pro / Gemini 3 Pro)、レビューログ付き反復改善、generate_schematic.py スクリプト) |
| **正当性** | 科学論文の Figure 1 は多くの場合「概念図」「ワークフロー図」であり、データプロットではない。K-Dense に scientific-schematics の完全スキル (CONSORT 対応、TikZ 生成、AI レビュー) が存在。SATORI は `publication-figures` でデータ可視化をカバーするが、**概念図生成**は完全に未カバーであり、論文作成ワークフローの重要な欠落。 |

---

### M-2. `scientific-lab-data-management`

| 項目 | 内容 |
|------|------|
| **カテゴリ** | M. ラボ・実験 |
| **優先度** | **MEDIUM** |
| **説明** | 電子実験ノート (ELN)・LIMS・クラウド科学計算基盤の統合パイプライン。Benchling (DNA 配列/レジストリ)、DNAnexus (ゲノミクス PaaS)、LatchBio (バイオインフォ Verified Workflows)、OMERO (バイオイメージングデータ管理)、Protocols.io (プロトコル共有・バージョン管理) を統合し、ウェットラボ → ドライラボのシームレスなデータフローを支援。 |
| **既存との区別** | `lab-automation` はロボティクス操作 (PyLabRobot)。ELN/LIMS (データ記録・管理)、クラウド計算基盤 (DNAnexus, LatchBio)、イメージングデータ管理 (OMERO) は別系統のワークフロー。 |
| **主要テクニック** | Benchling API (DNA 配列管理, オリゴ設計, プラスミドマップ, サンプルレジストリ)、DNAnexus (dxapp.json アプリ定義, ワークフロー実行, データ管理, ジョブモニタリング)、LatchBio (Verified Workflows — bulk RNA-seq, DESeq2, AlphaFold, バリアントコーリング)、OMERO (画像データインポート/エクスポート, ROI 管理, メタデータ付与)、Protocols.io (プロトコル作成, バージョン管理, DOI 発行, フォーク/共有) |
| **ToolUniverse SMCP ツール** | なし (K-Dense 独自の強み) |
| **K-Dense 対応スキル** | benchling-integration (Benchling ELN API)、dnanexus-integration (DNAnexus PaaS — dxapp.json, ワークフロー定義, ジョブ管理)、latchbio-integration (LatchBio — Verified Workflows, バイオインフォパイプライン)、omero-integration (OMERO — バイオイメージデータ管理)、protocolsio-integration (Protocols.io — プロトコル共有) |
| **正当性** | 現代の wet lab はデジタルデータ管理が必須。K-Dense に 5 つのクラウド/ELN 統合スキルが存在し、実験室 → 計算パイプライン → 結果のデータフロー全体をカバー。SATORI の `lab-automation` はロボット操作に特化しており、**データ管理・クラウド実行・プロトコル共有**は完全な空白。研究の再現性向上に直結。 |

---

## 4. ToolUniverse 未活用カテゴリ完全一覧

`default_config.py` から抽出した 90+ カテゴリのうち、SATORI v0.10.0 で**未活用**の主要カテゴリ:

### 本提案でカバーされるカテゴリ
```
chipatlas          → T-3 epigenomics-chromatin
fourdn             → T-3 epigenomics-chromatin
regulomedb         → T-3 epigenomics-chromatin
jaspar             → T-3 epigenomics-chromatin
remap              → T-3 epigenomics-chromatin
screen             → T-3 epigenomics-chromatin
pride              → F-6 proteomics-mass-spectrometry
pharmgkb           → P-2 pharmacogenomics
clinical_trials *  → H-4 clinical-trials-analytics (* 部分的に survival-clinical で使用)
geo                → F-7 gene-expression-transcriptomics
gtex_v2            → F-7 gene-expression-transcriptomics
expression_atlas   → F-7 gene-expression-transcriptomics
fda_drug_label *   → O-3 regulatory-science (* 部分的に pharmacovigilance で使用)
guidelines         → O-3 regulatory-science
```

### 今後の分析で検討すべきカテゴリ
```
spliceai            # スプライシング予測 → RNA 生物学スキル候補
complex_portal      # タンパク質複合体 → protein-structure 拡張候補
emdb                # 電子顕微鏡データ → protein-structure 拡張候補
impc                # マウス表現型 → ニッチ (モデル生物)
gtopdb              # 薬理学ガイド → drug-target 拡張候補
mpd                 # マウス表現型 → ニッチ
who_gho             # WHO グローバル健康 → epidemiology 拡張候補
umls                # 医学用語 → clinical-decision-support 拡張候補
euhealth            # EU 健康データ → epidemiology 拡張候補
```

---

## 5. 実装優先度マトリクス

| 順位 | スキル | ギャップ度 | TU ツール支援 | KD スキル支援 | 科学的インパクト | 総合スコア |
|------|--------|-----------|-------------|-------------|---------------|-----------|
| **1** | T-3 epigenomics-chromatin | Critical | ◎ (20+ ツール) | ○ (2 スキル) | ◎ | **9.5/10** |
| **2** | F-6 proteomics-ms | Critical | ○ (3+ ツール) | ◎ (2 スキル) | ◎ | **9.0/10** |
| **3** | E-4 neuroscience-electrophysiology | High | ✗ | ◎ (2 スキル) | ◎ | **8.5/10** |
| **4** | P-2 pharmacogenomics | High | ○ (5+ ツール) | ◎ (1 スキル) | ◎ | **8.5/10** |
| **5** | H-4 clinical-trials-analytics | High | ◎ (15+ ツール) | ◎ (1 スキル) | ◎ | **8.5/10** |
| **6** | O-3 regulatory-science | Moderate-High | ○ (既存) | ◎ (3 スキル) | ○ | **7.5/10** |
| **7** | F-7 gene-expression-transcriptomics | Moderate-High | ◎ (9+ ツール) | ◎ (2 スキル) | ○ | **7.5/10** |
| **8** | G-4 computational-materials | Moderate | ✗ | ◎ (1 スキル) | ○ | **7.0/10** |
| **9** | N-2 scientific-schematics | Moderate | ✗ | ◎ (1 スキル) | ○ | **6.5/10** |
| **10** | M-2 lab-data-management | Moderate | ✗ | ◎ (5 スキル) | ○ | **6.5/10** |

### 推奨実装フェーズ

| フェーズ | スキル | スキル数 | 目標バージョン |
|---------|--------|---------|--------------|
| Phase 3A (即時) | T-3, F-6, E-4, P-2, H-4 | 5 | v0.11.0 |
| Phase 3B (次期) | O-3, F-7, G-4, N-2, M-2 | 5 | v0.12.0 |

---

## 6. カバレッジ予測

| 指標 | 現在 (v0.11.0) | Phase 3A 後 | Phase 3B 後 |
|------|---------------|-------------|-------------|
| スキル数 | 66 | 71 | 76 |
| カテゴリ数 | 26 (A–Z) | 26 (A–Z) | 26 (A–Z) |
| ToolUniverse カバー率 | ~85% | ~92% | ~95% |
| K-Dense カバー率 | ~70% | ~82% | ~90% |
| 残存主要ギャップ | 10 | 5 | ~3 |

### Phase 3 完了後のカバレッジマトリクス

```
                              SATORI    ToolUniverse    K-Dense
──────────────────────────────────────────────────────────────────
エピゲノミクス/クロマチン     △→✅      ✅✅            ✅
プロテオミクス/MS            ❌→✅      ✅              ✅✅
神経科学/電気生理             △→✅      ❌              ✅✅
ファーマコゲノミクス          ❌→✅      ✅              ✅✅
臨床試験レジストリ解析        △→✅      ✅✅            ✅✅
規制科学 (FDA/ISO/USPTO)      ❌→✅      ✅              ✅✅
遺伝子発現/トランスクリプトミクス  △→✅  ✅✅            ✅✅
計算材料科学                  △→✅      ❌              ✅✅
科学スキマティクス            △→✅      ❌              ✅✅
ラボデータ管理/クラウド       △→✅      ❌              ✅✅
──────────────────────────────────────────────────────────────────
✅✅ = 非常に充実  ✅ = カバー  △ = 部分的  ❌ = 未カバー
```

---

## 7. 前回分析 (Tier 2) との対応関係

| 前回 Tier 2 項目 (GAP_ANALYSIS.md) | 本分析での対応 | 状態 |
|-----------------------------------|---------------|------|
| #10 proteomics-ms | → **F-6** そのまま昇格 | ✅ 提案 |
| #11 metabolic-modeling (COBRApy) | → Y-1 systems-biology に統合済み | ✅ 解決済み |
| #12 neuroscience | → **E-4** そのまま昇格 | ✅ 提案 |
| #13 multi-objective-optimization | 独立ギャップとして残存 | 🔜 次回検討 |
| #14 discrete-event-simulation | 独立ギャップとして残存 | 🔜 次回検討 |
| #15 large-scale-data (Dask/Polars) | 独立ギャップとして残存 | 🔜 次回検討 |
| #16 schematics-generation | → **N-2** そのまま昇格 | ✅ 提案 |
| #17 ELN/LIMS integration | → **M-2** 拡張して統合 | ✅ 提案 |
| #18 cloud-infrastructure | → **M-2** に統合 | ✅ 提案 |

**新規発見ギャップ (本分析独自)**:
- T-3 epigenomics-chromatin — ToolUniverse の chipatlas/fourdn/regulomedb/jaspar/remap/screen カテゴリの発見
- P-2 pharmacogenomics — ClinPGx (PharmGKB 後継) の K-Dense スキル発見
- H-4 clinical-trials-analytics — ToolUniverse の 15+ 臨床試験ツールの詳細解析
- O-3 regulatory-science — K-Dense の FDA/ISO/USPTO 3 スキル統合
- F-7 gene-expression-transcriptomics — ToolUniverse の GEO/GTEx/Expression Atlas カテゴリ統合

---

## 8. 結論

SATORI v0.11.0 は 66 スキル / 26 カテゴリで科学 AI スキルの広範なカバレッジを達成しているが、以下の 3 つの系統的ギャップが残存する:

### 1. オミクス技術の深度不足
- **エピゲノミクス**: ToolUniverse に 20+ の ChIP-Atlas/4DN/調節ゲノミクスツールが存在するが完全未カバー
- **プロテオミクス**: 第3の主要オミクス技術が欠落
- **トランスクリプトミクス**: GEO 400 万+ サンプルへの体系的アクセスが未整備

### 2. 臨床/トランスレーショナルの精密化
- **ファーマコゲノミクス**: 処方前の遺伝子型ベース薬物選択 (CPIC ガイドライン) が未カバー
- **臨床試験レジストリ**: ClinicalTrials.gov の 15+ ツールが未活用
- **規制科学**: FDA/ISO/USPTO の横断的データアクセスが未整備

### 3. 実験室-計算インフラの連携
- **神経科学**: SpikeInterface/NeuroKit2 による専門ワークフローが欠落
- **計算材料科学**: Pymatgen/Materials Project による in silico 材料設計が未カバー
- **ラボデータ管理**: ELN/LIMS/クラウド基盤の統合が未実装
- **科学スキマティクス**: 概念図生成が未対応

Phase 3A (5 スキル) の即時実装で 66 → 71 スキルとなり、特に ToolUniverse の最も充実したカテゴリ群 (chipatlas 20+ ツール, clinical_trials 15+ ツール, pharmgkb 5+ ツール) を初めて活用可能にする。Phase 3B (5 スキル) で 76 スキルに到達し、両リポとのカバレッジギャップを 95%/90% まで解消できる。
