# SATORI v0.17.0 スキル拡張ギャップ分析

> 作成日: 2026-02-22  
> 対象バージョン: v0.17.0 (132 スキル / 26 カテゴリ A–Z, 79 TU-linked)  
> 前回分析: GAP_ANALYSIS_v0.16.md (124 → 132 へ実装済み)

---

## 1. 調査概要

### 1.1 背景

v0.16.0 で 124 スキル / 74 ToolUniverse 連携を達成。Phase 9 (v0.17) では ToolUniverse Tier 3 候補（ENCODE エピゲノムアトラス・ChIP-Atlas・HCA Human Cell Atlas・GEO 発現プロファイル・PBDB 古生物学）と K-Dense 残存ギャップ（squidpy-advanced 高度空間解析・metabolic-atlas 代謝アトラス）に加え、直接 API スキル（環境地理空間データ・寄生虫ゲノミクス）を追加した。

### 1.2 調査ソース

| リポジトリ | 規模 | 最終確認日 |
|-----------|------|-----------|
| **mims-harvard/ToolUniverse** | `default_config.py`: **188 カテゴリ**, 1229+ ツール | 2026-02-22 |
| **K-Dense-AI/claude-scientific-skills** | **140 スキル** (databases: 28+, packages: 55+, integrations: 15+, analysis: 30+) | 2026-02-22 |
| **SATORI v0.17.0** | **132 スキル**, **79 TU-linked**, 26 カテゴリ (A–Z) | 現在 |

### 1.3 TU キー名の調査結果 (Phase 9)

GAP_ANALYSIS_v0.16 の Tier 3 で列挙された TU キーの実在性を subagent で検証:

| 想定キー | 実際の TU キー | 状態 |
|---------|--------------|------|
| `encode` | `encode` ✅ (L225, 5 tools) | TU 連携済み |
| `chipatlas` | `chipatlas` ✅ (L268, 1 tool) | TU 連携済み |
| `hca_tools` | `hca_tools` ✅ (L235, 1 tool) | TU 連携済み |
| `geo_profiles` | `geo` ✅ (L217, 3 tools) ※キー名異なる | TU 連携済み |
| `paleobiology` | `paleobiology` ✅ (L156, tools) | TU 連携済み |
| `soil_grids`, `worldclim` | **TU 不在** | 直接 API で実装 |
| `dictybase`, `plasmodb`, `vectorbase` | **TU 不在** | 直接 API で実装 |
| `beacon_network` | **TU 不在** | 未実装 |
| `icgc_argo` | 設定キー不在 (ツールは存在) | 未実装 |

→ Tier 3 から 5 TU キー (encode, chipatlas, hca_tools, geo, paleobiology) を発見・連携。残りは TU に不在のため直接 API で実装。

---

## 2. v0.17.0 で追加した TU カテゴリ

| # | Config Key | SATORI スキル | ツール数 |
|---|-----------|--------------|---------|
| 140 | `encode` | encode-screen | ~5 |
| 141 | `chipatlas` | encode-screen | ~1 |
| 142 | `hca_tools` | human-cell-atlas | ~1 |
| 143 | `geo` | geo-expression | ~3 |
| 144 | `paleobiology` | paleobiology | ~1 |

### 累積 TU カバレッジ: ~144/188 categories (76.6%)

---

## 3. v0.17.0 で追加した K-Dense カバレッジ (+2)

| # | K-Dense Key | SATORI スキル | 種別 |
|---|------------|--------------|------|
| 136 | `squidpy-advanced` | squidpy-advanced | analysis |
| 137 | `metabolic-atlas` | metabolic-atlas | database |

### 累積 K-Dense カバレッジ: ~137/140 (97.9%)

---

## 4. 残存ギャップ: ToolUniverse 未カバーカテゴリ (~44 categories)

### 4.1 Tier 4 候補 (v0.18.0+、極めて専門的)

| Config Key(s) | 説明 | 難易度 |
|--------------|------|--------|
| `icgc_argo` | ICGC ARGO がんデータ共有 | 低 |
| `beacon_network` | GA4GH Beacon v2 | 低 |
| `arrayexpress` | ArrayExpress マイクロアレイ | 低 (EBI DB で部分カバー) |
| `gwas_catalog` | GWAS Catalog REST API | 低 (disease-research で部分カバー) |
| `intact` | IntAct PPI | 低 (protein-interaction で部分カバー) |
| `uniprot_proteomes` | UniProt Proteomes API | 低 |
| `ensembl_regulation` | Ensembl Regulatory Build | 低 (regulatory-genomics で部分カバー) |
| その他 ~37 カテゴリ | 極めて専門的 (ニッチ分類群DB 等) | — |

---

## 5. 残存ギャップ: K-Dense 未カバー (~3 skills)

| # | K-Dense Key | 種別 | 優先度 | 推奨対応 |
|---|------------|------|--------|---------|
| 1 | `bio2bel` | package | 低 | ontology-enrichment 拡張 |
| 2 | `starfysh` | package | 低 | spatial-transcriptomics 拡張 |
| 3 | `augur` | package | ✅ | perturbation-analysis に統合済み |

> 注: `augur` は pertpy 経由で既存スキルに統合済み。実質残り ~2 スキル。

---

## 6. カテゴリ別分布 (v0.17.0 時点)

| カテゴリ | 名称 | スキル数 | TU-linked | v0.16→v0.17 差分 |
|---------|------|---------|-----------|-----------------|
| A | 基盤・ワークフロー | 16 | 2 | — |
| B | 統計・探索的解析 | 4 | 0 | — |
| C | 機械学習・モデリング | 3 | 0 | — |
| D | 実験計画・プロセス最適化 | 2 | 0 | — |
| E | 信号・スペクトル・時系列 | 4 | 0 | — |
| F | 生命科学・オミクス | **20** | **16** | **+2** |
| G | 化学・材料・イメージング | 8 | 2 | — |
| H | 臨床・疫学・メタ科学 | 5 | 2 | — |
| I | Deep Research・文献検索 | 3 | 2 | — |
| J | 創薬・ファーマコロジー | 7 | 5 | — |
| K | 構造生物学・タンパク質工学 | 5 | 5 | — |
| L | 精密医療・臨床意思決定 | 3 | 3 | — |
| M | 実験室自動化・データ管理 | 2 | 0 | — |
| N | 科学プレゼンテーション | 2 | 0 | — |
| O | 研究計画・グラント・規制 | 3 | 1 | — |
| P | ファーマコビジランス | 2 | 2 | — |
| Q | 腫瘍学・疾患研究 | 5 | 4 | — |
| R | 量子・先端計算 | 7 | 0 | — |
| S | 医用イメージング | 1 | 0 | — |
| T | シングルセル・空間・エピゲノミクス | **11** | **5** | **+3** |
| U | 免疫・感染症 | 2 | 2 | — |
| V | マイクロバイオーム・環境 | **8** | **3** | **+2** |
| W | システム生物学 | **3** | **2** | **+1** |
| X | 疫学・公衆衛生 | 3 | 2 | — |
| Y | 集団遺伝学 | 1 | 1 | — |
| Z | 科学テキストマイニング | 2 | 1 | — |
| **合計** | | **132** | **79** | **+8** |

---

## 7. v0.18.0+ ロードマップ候補

### Phase 10 (v0.18.0): TU Tier 4 + K-Dense 最終 — 目標 132→138+ スキル

| # | 推奨スキル名 | カテゴリ | TU keys | K-Dense |
|---|-------------|---------|---------|----------|
| 1 | scientific-icgc-cancer-data | Q | icgc_argo | — |
| 2 | scientific-beacon-network | F | beacon_network | — |
| 3 | scientific-starfysh-deconv | T | — | starfysh |
| 4 | scientific-bio2bel-integration | F | — | bio2bel |
| 5 | scientific-gwas-catalog | Y | gwas_catalog | — |
| 6 | scientific-arrayexpress | F | arrayexpress | — |

### 推定達成指標 (v0.18.0 後)
- スキル: ~138
- TU カバレッジ: ~148/188 (78.7%)
- K-Dense カバレッジ: ~139/140 (99.3%)

---

## 8. まとめ

v0.17.0 では以下を達成した:

1. **8 新スキル** (124→132): ENCODE/SCREEN 制御ゲノム、Human Cell Atlas、GEO 発現プロファイル、環境地理空間データ、古生物学、Metabolic Atlas、高度 Squidpy 空間解析、寄生虫ゲノミクス
2. **TU 連携 74→79** (+5 キー: encode, chipatlas, hca_tools, geo, paleobiology)
3. **K-Dense 135→137** (+2: squidpy-advanced, metabolic-atlas)
4. **TU キー実在性の精査** — GAP_ANALYSIS_v0.16 の Tier 3 想定キーのうち、5 キーが TU に実在し連携完了。残り (soil_grids/worldclim/dictybase/plasmodb 等) は TU に不在のため直接 API で独立実装

残り TU ~44 カテゴリ / K-Dense ~2 スキルのギャップは v0.18.0 以降で段階的に解消予定。
