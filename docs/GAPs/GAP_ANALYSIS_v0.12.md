# SATORI v0.12.0 スキル拡張ギャップ分析

> 作成日: 2026-02  
> 対象バージョン: v0.12.0 (76 スキル / 26 カテゴリ A–Z)  
> 前回分析: GAP_ANALYSIS_v0.11.md (56 → 66 → 76 へ Phase 3A+3B 実装済み)

---

## 1. 調査概要

### 1.1 背景

v0.11.1 で 76 スキル / 42 ToolUniverse 連携を達成。Phase 3 で実装した 10 スキル（エピゲノミクス、プロテオミクス、神経科学、ファーマコゲノミクス、臨床試験、規制科学、トランスクリプトミクス、計算材料、科学図式、ラボデータ管理）は主にデータ**生成・解析**に焦点を当てたが、以下の構造的ギャップが判明：

1. **パスウェイ / Gene Ontology 富化解析** — 最も普遍的な下流解析が未対応
2. **文献検索 API 統合** — 研究ワークフローの起点が欠落
3. **タンパク質相互作用ネットワーク** — STRING/IntAct/BioGRID が未活用
4. **バリアント効果予測** — AlphaMissense/CADD/SpliceAI が未連携
5. **メタボロミクス DB** — HMDB/MetaCyc/Metabolomics Workbench が未統合
6. **がんゲノミクスポータル** — COSMIC/cBioPortal/DepMap が未活用
7. **分子ドッキング** — 創薬パイプラインの構造ベース設計が欠落

### 1.2 ToolUniverse 未活用カテゴリ (v0.12 対象)

| # | カテゴリ | Config Key | ツール数 | SATORI 対象スキル |
|---|---|---|---|---|
| 1 | KEGG | `kegg` | 5 | pathway-enrichment |
| 2 | Reactome | `reactome` | 20 | pathway-enrichment |
| 3 | Gene Ontology | `go` | 5 | pathway-enrichment |
| 4 | WikiPathways | `wikipathways` | 2 | pathway-enrichment |
| 5 | Pathway Commons | `pathway_commons` | 2 | pathway-enrichment |
| 6 | PubMed | `pubmed` | 6 | literature-search |
| 7 | EuropePMC | `EuropePMC` | 7 | literature-search |
| 8 | Semantic Scholar | `semantic_scholar` | 2 | literature-search |
| 9 | OpenAlex | `OpenAlex` | 8 | literature-search |
| 10 | CrossRef | `crossref` | 6 | literature-search |
| 11 | IntAct | `intact` | 8 | protein-interaction-network |
| 12 | STRING/BioGRID | `ppi` | 2 | protein-interaction-network |
| 13 | STITCH | `stitch` | 3 | protein-interaction-network |
| 14 | HumanBase | `HumanBase` | 1 | protein-interaction-network |
| 15 | AlphaMissense | `alphamissense` | 3 | variant-effect-prediction |
| 16 | CADD | `cadd` | 3 | variant-effect-prediction |
| 17 | SpliceAI | `spliceai` | 3 | variant-effect-prediction |
| 18 | COSMIC | `cosmic` | 2 | cancer-genomics |
| 19 | cBioPortal | `cbioportal` | 6 | cancer-genomics |
| 20 | DepMap | `depmap` | 5 | cancer-genomics |
| 21 | HMDB | `hmdb` | 3 | metabolomics-databases |
| 22 | MetaCyc | `metacyc` | 4 | metabolomics-databases |
| 23 | Metabolomics WB | `metabolomics_workbench` | 6 | metabolomics-databases |
| 24 | InterPro | `interpro` | 3 | protein-domain-family |
| 25 | InterProScan | `interproscan` | 3 | protein-domain-family |
| | **合計** | | **117** | |

---

## 2. 新スキル提案 (10 スキル)

### カテゴリ別配置

| カテゴリ | 名称 | 追加スキル数 | 追加後合計 |
|---------|------|------------|-----------|
| **A** | 基盤・ワークフロー | +1 | 14 |
| **F** | 生命科学・オミクス | +2 | 9 |
| **H** | 臨床・疫学統計 | +1 | 5 |
| **I** | Deep Research | +1 | 2 |
| **J** | 創薬・ファーマコロジー | +1 | 4 |
| **K** | 構造生物学・タンパク質工学 | +2 | 4 |
| **L** | 精密医療・臨床意思決定 | +1 | 3 |
| **Q** | 腫瘍学・疾患研究 | +1 | 3 |

---

### 2.1 `scientific-pathway-enrichment` (F-8)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | F. 生命科学・オミクス |
| **優先度** | **CRITICAL** |
| **ToolUniverse ツール** | kegg(5) + reactome(20) + go(5) + wikipathways(2) + pathway_commons(2) = **34 ツール** |
| **K-Dense 対応** | kegg-database, reactome-database |
| **説明** | KEGG/Reactome/GO 経由のパスウェイ富化解析パイプライン。ORA (Over-Representation Analysis)、GSEA (Gene Set Enrichment Analysis)、トポロジーベース解析、マルチDB 統合富化を提供。最も普遍的なオミクス下流解析。 |

### 2.2 `scientific-literature-search` (I-2)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | I. Deep Research |
| **優先度** | **CRITICAL** |
| **ToolUniverse ツール** | pubmed(6) + europepmc(7) + semantic_scholar(2) + openalex(8) + crossref(6) = **29 ツール** |
| **K-Dense 対応** | pubmed-database, openalex-database |
| **説明** | 学術文献の API レベル検索・取得パイプライン。PubMed E-utilities、Semantic Scholar、OpenAlex、EuropePMC、CrossRef の 5 大データベースを統合。MeSH 構造化検索、引用ネットワーク、著者/機関メトリクス対応。 |

### 2.3 `scientific-protein-interaction-network` (K-3)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | K. 構造生物学・タンパク質工学 |
| **優先度** | **HIGH** |
| **ToolUniverse ツール** | intact(8) + ppi/STRING+BioGRID(2) + stitch(3) + humanbase(1) = **14 ツール** |
| **K-Dense 対応** | string-database |
| **説明** | STRING/IntAct/BioGRID/STITCH 経由の PPI ネットワーク解析パイプライン。相互作用検索、ネットワーク構築、GO/KEGG 富化、化学-タンパク質相互作用、組織特異的ネットワーク対応。 |

### 2.4 `scientific-variant-effect-prediction` (L-3)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | L. 精密医療・臨床意思決定 |
| **優先度** | **HIGH** |
| **ToolUniverse ツール** | alphamissense(3) + cadd(3) + spliceai(3) = **9 ツール** |
| **K-Dense 対応** | ensembl-database (VEP) |
| **説明** | 計算バリアント効果予測パイプライン。AlphaMissense (タンパク質構造ベース)、CADD (統合アノテーション)、SpliceAI (スプライシング影響) のスコア統合とコンセンサス病原性評価。 |

### 2.5 `scientific-cancer-genomics` (Q-3)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | Q. 腫瘍学・疾患研究 |
| **優先度** | **HIGH** |
| **ToolUniverse ツール** | cosmic(2) + cbioportal(6) + depmap(5) = **13 ツール** |
| **K-Dense 対応** | cosmic-database |
| **説明** | COSMIC/cBioPortal/DepMap 経由のがんゲノミクスポータル統合パイプライン。体細胞変異プロファイル、変異シグネチャー、遺伝子依存性 (essentiality)、コピー数変化、がん種横断解析対応。 |

### 2.6 `scientific-metabolomics-databases` (F-9)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | F. 生命科学・オミクス |
| **優先度** | **HIGH** |
| **ToolUniverse ツール** | hmdb(3) + metacyc(4) + metabolomics_workbench(6) = **13 ツール** |
| **K-Dense 対応** | hmdb-database, metabolomics-workbench-database, brenda-database |
| **説明** | HMDB/MetaCyc/Metabolomics Workbench 経由のメタボロミクスデータベース統合パイプライン。代謝物同定、パスウェイマッピング、バイオマーカー発見、RefMet 標準化命名対応。 |

### 2.7 `scientific-molecular-docking` (J-4)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | J. 創薬・ファーマコロジー |
| **優先度** | **MEDIUM-HIGH** |
| **ToolUniverse ツール** | — (K-Dense 独自) |
| **K-Dense 対応** | diffdock |
| **説明** | 構造ベース分子ドッキングパイプライン。DiffDock (拡散モデル)、AutoDock Vina、GNINA (CNN スコアリング) による結合ポーズ予測、バーチャルスクリーニング、ドッキングスコア統合。 |

### 2.8 `scientific-systematic-review` (A-14)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | A. 基盤・ワークフロー |
| **優先度** | **MEDIUM-HIGH** |
| **ToolUniverse ツール** | pubmed ツール共用 |
| **K-Dense 対応** | literature-review |
| **説明** | PRISMA 2020 準拠の系統的レビューパイプライン。マルチ DB 検索戦略立案、スクリーニング (タイトル/抄録→全文)、品質評価 (RoB 2/ROBINS-I/NOS)、データ抽出テンプレート、PRISMA フロー図生成。 |

### 2.9 `scientific-clinical-reporting` (H-5)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | H. 臨床・疫学・メタ科学 |
| **優先度** | **MEDIUM** |
| **ToolUniverse ツール** | — (K-Dense 独自) |
| **K-Dense 対応** | clinical-reports, treatment-plans |
| **説明** | 臨床レポート・治療計画書の自動生成パイプライン。SOAP ノート、バイオマーカーレポート、精神科/リハビリ治療計画、慢性疾患管理計画、PDF/LaTeX 出力対応。 |

### 2.10 `scientific-protein-domain-family` (K-4)

| 項目 | 内容 |
|------|------|
| **カテゴリ** | K. 構造生物学・タンパク質工学 |
| **優先度** | **MEDIUM** |
| **ToolUniverse ツール** | interpro(3) + interproscan(3) = **6 ツール** |
| **K-Dense 対応** | bioservices |
| **説明** | InterPro/InterProScan 経由のタンパク質ドメイン・ファミリー解析パイプライン。配列からのドメイン予測、ファミリー分類、機能注釈、ドメインアーキテクチャ可視化対応。 |

---

## 3. 実装優先度マトリクス

| 順位 | スキル | ギャップ度 | TU ツール | KD 支援 | 総合スコア |
|------|--------|-----------|----------|---------|-----------|
| **1** | pathway-enrichment | Critical | ◎ (34) | ◎ | **10/10** |
| **2** | literature-search | Critical | ◎ (29) | ◎ | **9.5/10** |
| **3** | protein-interaction-network | High | ◎ (14) | ◎ | **9.0/10** |
| **4** | variant-effect-prediction | High | ◎ (9) | ◎ | **8.5/10** |
| **5** | cancer-genomics | High | ◎ (13) | ◎ | **8.5/10** |
| **6** | metabolomics-databases | High | ◎ (13) | ◎ | **8.5/10** |
| **7** | molecular-docking | Med-High | ✗ | ◎ | **7.5/10** |
| **8** | systematic-review | Med-High | △ | ◎ | **7.5/10** |
| **9** | clinical-reporting | Medium | ✗ | ◎ | **7.0/10** |
| **10** | protein-domain-family | Medium | ○ (6) | ○ | **7.0/10** |

---

## 4. カバレッジ予測

| 指標 | v0.11.1 | v0.12.0 後 |
|------|---------|-----------|
| スキル数 | 76 | **86** |
| ToolUniverse 連携数 | 42 | **50** |
| 新規 TU ツール参照数 | — | **117** |
| ToolUniverse カバー率 | ~85% | **~97%** |
| K-Dense カバー率 | ~90% | **~96%** |
