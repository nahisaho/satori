# SATORI v0.15.0 スキル拡張ギャップ分析

> 作成日: 2026-02-20  
> 対象バージョン: v0.15.0 (116 スキル / 26 カテゴリ A–Z, 70 TU-linked)  
> 前回分析: GAP_ANALYSIS_v0.14.md (106 → 116 へ実装済み)

---

## 1. 調査概要

### 1.1 背景

v0.14.0 で 106 スキル / 66 ToolUniverse 連携を達成。Phase 7 (v0.15) では ToolUniverse Tier 1 候補（ChEMBL アッセイ、Ensembl ゲノミクス、STRING/STITCH PPI、Expression Atlas 発現比較）と K-Dense ギャップ（MD シミュレーション、摂動解析、高度イメージング、深層化学、scVI 統合）を中心に 10 スキルを追加した。

### 1.2 調査ソース

| リポジトリ | 規模 | 最終確認日 |
|-----------|------|-----------|
| **mims-harvard/ToolUniverse** | `default_config.py`: **188 カテゴリ**, 1229+ ツール | 2026-02-20 |
| **K-Dense-AI/claude-scientific-skills** | **140 スキル** (databases: 28+, packages: 55+, integrations: 15+, analysis: 30+) | 2026-02-20 |
| **SATORI v0.15.0** | **116 スキル**, **70 TU-linked**, 26 カテゴリ (A–Z) | 現在 |

### 1.3 TU キー名の調査結果

GAP_ANALYSIS_v0.14 の Tier 1 で列挙された `chembl_assay`, `string_network`, `ensembl_rest`, `pfam`, `bgee` 等の独立キーは、ToolUniverse の実際の構造では異なるキー体系であった:

| 想定キー | 実際の TU キー | ツール数 |
|---------|--------------|---------|
| `chembl_assay`, `chembl_activity` | `ChEMBL` (統一) | 29 |
| `ensembl_rest`, `ensembl_vep` | `ensembl` (統一) | 19 |
| `string_network`, `string_enrichment` | `ppi` (2) + `stitch` (3) | 5 |
| `pfam`, `cath` | `interpro` (既存) | — |
| `bgee`, `expression_comparison` | `expression_atlas` | 4 |
| `flybase`, `wormbase`, `zfin`, `rgd`, `mgi` | **TU 不在** | 0 |

→ `pfam`/`cath` は既存 interpro でカバー済み、`protein-classification` は不要と判断。代わりに `scientific-scvi-integration` を追加。

---

## 2. v0.15.0 で追加した TU カテゴリ

| # | Config Key | SATORI スキル | ツール数 |
|---|-----------|--------------|---------|
| 131 | `ChEMBL` (拡張) | chembl-assay-mining | 19 |
| 132 | `ensembl` (拡張) | ensembl-genomics | 16 |
| 133 | `ppi` (拡張) | string-network-api | 2 |
| 134 | `stitch` (拡張) | string-network-api | 3 |
| 135 | `expression_atlas` | expression-comparison | 4 |

### 累積 TU カバレッジ: ~135/188 categories (71.8%)

> 注: ChEMBL, ensembl, ppi/stitch は既存スキル (drug-target, genome-sequence-tools, protein-interaction) でも部分的に参照されていたが、v0.15.0 で専用スキルとして深掘り。

---

## 3. v0.15.0 で追加した K-Dense カバレッジ (+9)

| # | K-Dense Key | SATORI スキル | 種別 |
|---|------------|--------------|------|
| 125 | `mdanalysis` | md-simulation | package |
| 126 | `openff` | md-simulation | package |
| 127 | `cellprofiler` | advanced-imaging | package |
| 128 | `napari` | advanced-imaging | package |
| 129 | `cellpose` | advanced-imaging | package |
| 130 | `deepchem` | deep-chemistry | package |
| 131 | `pertpy` | perturbation-analysis | package |
| 132 | `scib` | perturbation-analysis | package |
| 133 | `scvi-tools` | scvi-integration | package |

### 累積 K-Dense カバレッジ: ~133/140 (95.0%)

---

## 4. 残存ギャップ: ToolUniverse 未カバーカテゴリ (~53 categories)

### 4.1 Tier 2 候補 (v0.16.0 推奨)

| Config Key(s) | 説明 | 推奨スキル名 | カテゴリ |
|--------------|------|-------------|---------|
| `toxnet`, `tox21`, `ctd` | 毒性/環境衛生 DB | scientific-toxicology-env | X |
| `silva`, `greengenes2` | rRNA リファレンス DB | scientific-rrna-taxonomy | V |
| `genbank_submission`, `sra_tools` | データ登録/SRA | scientific-data-submission | A |
| `cellminer`, `nci60` | NCI-60 薬剤応答 | scientific-nci60-screening | J |
| `flybase`, `wormbase`, `zfin`, `rgd`, `mgi` | (TU 不在 — 直接 REST で対応済み) | — | — |

### 4.2 Tier 3 候補 (v0.17.0+、ニッチ / 低優先)

| Config Key(s) | 説明 |
|--------------|------|
| `dictybase`, `plasmodb`, `vectorbase` | 稀少モデル生物 |
| `plant_reactome`, `tair` | 植物バイオ |
| `fish_base`, `ocean_biogeographic` | 海洋バイオ |
| `soil_grids`, `worldclim` | 環境地理 |
| `icgc_argo`, `beacon_network` | がんデータ共有 |
| その他 ~30 カテゴリ | 極めて専門的 |

---

## 5. 残存ギャップ: K-Dense 未カバー (~7 skills)

| # | K-Dense Key | 種別 | 優先度 | 推奨対応 |
|---|------------|------|--------|---------|
| 1 | `metabolic-atlas` | database | 中 | metabolic-modeling 拡張 |
| 2 | `bio2bel` | package | 低 | ontology-enrichment 拡張 |
| 3 | `signac` | package | 低 | single-cell-genomics 拡張 (scATAC) |
| 4 | `squidpy-advanced` | analysis | 低 | spatial-transcriptomics 拡張 |
| 5 | `rapids-singlecell` | package | 低 | single-cell-genomics 拡張 (GPU) |
| 6 | `starfysh` | package | 低 | spatial-transcriptomics 拡張 |
| 7 | `augur` | package | ✅ | perturbation-analysis に統合済み (pertpy 経由) |

> 注: `augur` は pertpy ライブラリ経由で `scientific-perturbation-analysis` に統合済み。実質残り ~6 スキル。

---

## 6. カテゴリ別分布 (v0.15.0 時点)

| カテゴリ | 名称 | スキル数 | TU-linked | v0.14→v0.15 差分 |
|---------|------|---------|-----------|-----------------|
| A | 基盤・ワークフロー | 15 | 2 | — |
| B | 統計・探索的解析 | 4 | 0 | — |
| C | 機械学習・モデリング | 3 | 0 | — |
| D | 実験計画・プロセス最適化 | 2 | 0 | — |
| E | 信号・スペクトル・時系列 | 4 | 0 | — |
| F | 生命科学・オミクス | **18** | **15** | **+4** |
| G | 化学・材料・イメージング | **8** | **2** | **+4** |
| H | 臨床・疫学・メタ科学 | 5 | 2 | — |
| I | Deep Research・文献検索 | 3 | 2 | — |
| J | 創薬・ファーマコロジー | 6 | 5 | — |
| K | 構造生物学・タンパク質工学 | 5 | 5 | — |
| L | 精密医療・臨床意思決定 | 3 | 3 | — |
| M | 実験室自動化・データ管理 | 2 | 0 | — |
| N | 科学プレゼンテーション | 2 | 0 | — |
| O | 研究計画・グラント・規制 | 3 | 1 | — |
| P | ファーマコビジランス | 2 | 2 | — |
| Q | 腫瘍学・疾患研究 | 5 | 4 | — |
| R | 量子・先端計算 | 7 | 0 | — |
| S | 医用イメージング | 1 | 0 | — |
| T | シングルセル・空間・エピゲノミクス | **6** | **3** | **+2** |
| U | 免疫・感染症 | 2 | 2 | — |
| V | マイクロバイオーム・環境 | 3 | 1 | — |
| W | システム生物学 | 2 | 2 | — |
| X | 疫学・公衆衛生 | 2 | 2 | — |
| Y | 集団遺伝学 | 1 | 1 | — |
| Z | 科学テキストマイニング | 2 | 1 | — |
| **合計** | | **116** | **70** | **+10** |

---

## 7. v0.16.0+ ロードマップ候補

### Phase 8 (v0.16.0): TU Tier 2 + K-Dense 残存 — 目標 116→124 スキル

| # | 推奨スキル名 | カテゴリ | TU keys | K-Dense |
|---|-------------|---------|---------|----------|
| 1 | scientific-toxicology-env | X | toxnet, tox21, ctd | — |
| 2 | scientific-rrna-taxonomy | V | silva, greengenes2 | — |
| 3 | scientific-data-submission | A | genbank_submission, sra_tools | — |
| 4 | scientific-nci60-screening | J | cellminer, nci60 | — |
| 5 | scientific-plant-biology | V | plant_reactome, tair | — |
| 6 | scientific-marine-ecology | V | fish_base, ocean_biogeographic | — |
| 7 | scientific-scatac-signac | T | — | signac |
| 8 | scientific-gpu-singlecell | T | — | rapids-singlecell |

### 推定達成指標 (v0.16.0 後)
- スキル: ~124
- TU カバレッジ: ~147/188 (78.2%)
- K-Dense カバレッジ: ~137/140 (97.9%)

---

## 8. まとめ

v0.15.0 では以下を達成した:

1. **10 新スキル** (106→116): ChEMBL アッセイマイニング、Ensembl ゲノミクス、STRING/STITCH PPI、発現比較、モデル生物 DB、MD シミュレーション、摂動解析、高度イメージング、深層化学、scVI 統合
2. **TU 連携 66→70** (+4 スキル、~44 新規ツール統合)
3. **K-Dense 124→133** (+9: mdanalysis, openff, cellprofiler, napari, cellpose, deepchem, pertpy, scib, scvi-tools)
4. **TU キー体系の精査** — 想定キーと実際の TU 構造のマッピングを確認し、重複排除

残り TU ~53 カテゴリ / K-Dense ~6 スキルのギャップは v0.16.0 以降で段階的に解消予定。
