# SATORI v0.14.0 スキル拡張ギャップ分析

> 作成日: 2026-02-19  
> 対象バージョン: v0.14.0 (106 スキル / 26 カテゴリ A–Z, 66 TU-linked)  
> 前回分析: GAP_ANALYSIS_v0.13.md (96 → 106 へ実装済み)

---

## 1. 調査概要

### 1.1 背景

v0.13.0 で 96 スキル / 59 ToolUniverse 連携を達成。Phase 6 (v0.14) では ToolUniverse Tier 1 候補（プレプリントアーカイブ、公衆衛生、オントロジー、EBI データベース、細胞株）と K-Dense ギャップ（系統解析、強化学習、記号数学）を中心に 10 スキルを追加、パイプライン統合アーキテクチャを README に明示化した。

### 1.2 調査ソース

| リポジトリ | 規模 | 最終確認日 |
|-----------|------|-----------|
| **mims-harvard/ToolUniverse** | `default_config.py`: **188 カテゴリ**, 1229+ ツール | 2026-02-19 |
| **K-Dense-AI/claude-scientific-skills** | **140 スキル** (databases: 28+, packages: 55+, integrations: 15+, analysis: 30+) | 2026-02-19 |
| **SATORI v0.14.0** | **106 スキル**, **66 TU-linked**, 26 カテゴリ (A–Z) | 現在 |

---

## 2. v0.14.0 で追加した TU カテゴリ (~33 config keys)

| # | Config Key | SATORI スキル |
|---|-----------|--------------|
| 98 | `biorxiv` | preprint-archive |
| 99 | `medrxiv` | preprint-archive |
| 100 | `arxiv` | preprint-archive |
| 101 | `pmc` | preprint-archive |
| 102 | `doaj` | preprint-archive |
| 103 | `unpaywall` | preprint-archive |
| 104 | `hal` | preprint-archive |
| 105 | `core` | preprint-archive |
| 106 | `zenodo` | preprint-archive |
| 107 | `openaire` | preprint-archive |
| 108 | `osf_preprints` | preprint-archive |
| 109 | `fatcat` | preprint-archive |
| 110 | `dblp` | preprint-archive |
| 111 | `nhanes` | public-health-data |
| 112 | `health_disparities` | public-health-data |
| 113 | `medlineplus` | public-health-data |
| 114 | `odphp` | public-health-data |
| 115 | `rxnorm` | public-health-data |
| 116 | `guidelines_tools` | public-health-data |
| 117 | `efo` | ontology-enrichment |
| 118 | `ols` | ontology-enrichment |
| 119 | `enrichr` | ontology-enrichment |
| 120 | `umls` | ontology-enrichment |
| 121 | `ebi_search` | ebi-databases |
| 122 | `ena_browser` | ebi-databases |
| 123 | `biostudies` | ebi-databases |
| 124 | `dbfetch` | ebi-databases |
| 125 | `metabolights` | ebi-databases |
| 126 | `cellosaurus` | cell-line-resources |
| 127 | `regulomedb` | regulatory-genomics |
| 128 | `remap` | regulatory-genomics |
| 129 | `fourdn_portal` | regulatory-genomics |
| 130 | `pubtator` | biomedical-pubtator |

### 累積 TU カバレッジ: ~130/188 categories (69.1%)

---

## 3. v0.14.0 で追加した K-Dense カバレッジ (+5)

| # | K-Dense Key | SATORI スキル | 種別 |
|---|------------|--------------|------|
| 119 | `ena-database` | ebi-databases | database |
| 120 | `etetoolkit` | phylogenetics | package |
| 121 | `scikit-bio` | phylogenetics | package |
| 122 | `stable-baselines3` | reinforcement-learning | package |
| 123 | `pufferlib` | reinforcement-learning | package |
| 124 | `sympy` | symbolic-mathematics | package |

### 累積 K-Dense カバレッジ: ~124/140 (88.6%)

---

## 4. 残存ギャップ: ToolUniverse 未カバーカテゴリ (~58 categories)

### 4.1 Tier 1 候補 (v0.15.0 推奨、高頻度 API + 明確な科学ユースケース)

| Config Key(s) | 説明 | 推奨スキル名 | カテゴリ |
|--------------|------|-------------|---------|
| `chembl_assay`, `chembl_activity` | ChEMBL アッセイ/活性詳細 | scientific-chembl-assay-mining | G |
| `string_network`, `string_enrichment` | STRING ネットワーク直接 API | scientific-string-network-api | F |
| `ensembl_rest`, `ensembl_vep` | Ensembl REST 変異アノテーション | scientific-ensembl-genomics | F |
| `pfam`, `cath` | タンパク質ドメイン分類 DB | scientific-protein-classification | K |
| `bgee`, `expression_comparison` | 発現比較 (Bgee) | scientific-expression-comparison | F |

### 4.2 Tier 2 候補 (v0.16.0、専門性高)

| Config Key(s) | 説明 | 推奨スキル名 | カテゴリ |
|--------------|------|-------------|---------|
| `flybase`, `wormbase`, `zfin`, `rgd`, `mgi` | モデル生物 DB 統合 | scientific-model-organism-db | F |
| `toxnet`, `tox21`, `ctd` | 毒性/環境衛生 | scientific-toxicology-env | X |
| `silva`, `greengenes2` | rRNA リファレンス DB | scientific-rrna-taxonomy | V |
| `genbank_submission`, `sra_tools` | データ登録/SRA | scientific-data-submission | A |
| `cellminer`, `nci60` | NCI-60 薬剤応答 | scientific-nci60-screening | J |

### 4.3 Tier 3 候補 (v0.17.0+、ニッチ / 低優先)

| Config Key(s) | 説明 |
|--------------|------|
| `dictybase`, `plasmodb`, `vectorbase` | 稀少モデル生物 |
| `plant_reactome`, `tair` | 植物バイオ |
| `fish_base`, `ocean_biogeographic` | 海洋バイオ |
| `soil_grids`, `worldclim` | 環境地理 |
| `icgc_argo`, `beacon_network` | がんデータ共有 |
| その他 ~30 カテゴリ | 極めて専門的 |

---

## 5. 残存ギャップ: K-Dense 未カバー (~16 skills)

| # | K-Dense Key | 種別 | 優先度 | 推奨対応 |
|---|------------|------|--------|---------|
| 1 | `metabolic-atlas` | database | 中 | metabolic-modeling 拡張 |
| 2 | `bio2bel` | package | 低 | ontology-enrichment 拡張 |
| 3 | `cellprofiler` | package | 中 | image-analysis 拡張 |
| 4 | `deepchem` | package | 中 | 新スキル or cheminformatics 拡張 |
| 5 | `openff` | package | 中 | computational-materials 拡張 |
| 6 | `mdanalysis` | package | 中 | 新スキル: MD 解析 |
| 7 | `signac` | package | 低 | single-cell-genomics 拡張 (scATAC) |
| 8 | `augur` | package | 低 | single-cell-genomics 拡張 |
| 9 | `scvi-tools` | package | 中 | 新スキル: scVI 統合 |
| 10 | `squidpy-advanced` | analysis | 低 | spatial-transcriptomics 拡張 |
| 11 | `pertpy` | package | 中 | 新スキル: 摂動解析 |
| 12 | `scib` | package | 低 | single-cell-genomics 拡張 |
| 13 | `rapids-singlecell` | package | 低 | single-cell-genomics 拡張 (GPU) |
| 14 | `napari` | package | 中 | image-analysis 拡張 |
| 15 | `cellpose` | package | 中 | image-analysis 拡張 |
| 16 | `starfysh` | package | 低 | spatial-transcriptomics 拡張 |

---

## 6. カテゴリ別分布 (v0.14.0 時点)

| カテゴリ | 名称 | スキル数 | TU-linked |
|---------|------|---------|-----------|
| A | 基盤・ワークフロー | 15 | 2 |
| B | 統計・探索的解析 | 4 | 0 |
| C | 機械学習・モデリング | 3 | 0 |
| D | 実験計画・プロセス最適化 | 2 | 0 |
| E | 信号・スペクトル・時系列 | 4 | 0 |
| F | 生命科学・オミクス | 14 | 11 |
| G | 化学・材料・イメージング | 4 | 1 |
| H | 臨床・疫学・メタ科学 | 5 | 2 |
| I | Deep Research・文献検索 | 3 | 2 |
| J | 創薬・ファーマコロジー | 6 | 5 |
| K | 構造生物学・タンパク質工学 | 5 | 5 |
| L | 精密医療・臨床意思決定 | 3 | 3 |
| M | 実験室自動化・データ管理 | 2 | 0 |
| N | 科学プレゼンテーション | 2 | 0 |
| O | 研究計画・グラント・規制 | 3 | 1 |
| P | ファーマコビジランス | 2 | 2 |
| Q | 腫瘍学・疾患研究 | 5 | 4 |
| R | 量子・先端計算 | 7 | 0 |
| S | 医用イメージング | 1 | 0 |
| T | シングルセル・空間・エピゲノミクス | 4 | 3 |
| U | 免疫・感染症 | 2 | 2 |
| V | マイクロバイオーム・環境 | 3 | 1 |
| W | システム生物学 | 2 | 2 |
| X | 疫学・公衆衛生 | 2 | 2 |
| Y | 集団遺伝学 | 1 | 1 |
| Z | 科学テキストマイニング | 2 | 1 |
| **合計** | | **106** | **66** (~50 unique TU-linked) |

---

## 7. v0.15.0+ ロードマップ候補

### Phase 7 (v0.15.0): TU Tier 1 + K-Dense 中優先 — 目標 106→116 スキル

| # | 推奨スキル名 | カテゴリ | TU keys | K-Dense |
|---|-------------|---------|---------|----------|
| 1 | scientific-chembl-assay-mining | G | chembl_assay, chembl_activity | — |
| 2 | scientific-ensembl-genomics | F | ensembl_rest, ensembl_vep | — |
| 3 | scientific-string-network-api | F | string_network, string_enrichment | — |
| 4 | scientific-protein-classification | K | pfam, cath | — |
| 5 | scientific-expression-comparison | F | bgee, expression_comparison | — |
| 6 | scientific-model-organism-db | F | flybase, wormbase, zfin, rgd, mgi | — |
| 7 | scientific-md-simulation | G | — | mdanalysis, openff |
| 8 | scientific-perturbation-analysis | T | — | pertpy, scib |
| 9 | scientific-advanced-imaging | G | — | cellprofiler, napari, cellpose |
| 10 | scientific-deep-chemistry | G | — | deepchem |

### 推定達成指標 (v0.15.0 後)
- スキル: ~116
- TU カバレッジ: ~140/188 (74.5%)
- K-Dense カバレッジ: ~132/140 (94.3%)

---

## 8. まとめ

v0.14.0 では以下を達成した:

1. **10 新スキル** (96→106): プレプリント横断検索、公衆衛生データ、オントロジー濃縮、EBI 統合 DB、細胞株リソース、制御ゲノミクス、系統解析、強化学習、記号数学、バイオメディカルテキストマイニング
2. **TU 連携 59→66** (+7 スキル、~33 新規 config keys、~48 新規ツール)
3. **K-Dense 118→124** (+6: ena-database, etetoolkit, scikit-bio, stable-baselines3, pufferlib, sympy)
4. **パイプライン統合** — README にデータアクセスフロー図を新設、各 SKILL.md にパイプライン統合セクション追加

残り TU ~58 カテゴリ / K-Dense ~16 スキルのギャップは v0.15.0 以降で段階的に解消予定。
