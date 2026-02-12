# SATORI v0.16.0 スキル拡張ギャップ分析

> 作成日: 2026-02-21  
> 対象バージョン: v0.16.0 (124 スキル / 26 カテゴリ A–Z, 74 TU-linked)  
> 前回分析: GAP_ANALYSIS_v0.15.md (116 → 124 へ実装済み)

---

## 1. 調査概要

### 1.1 背景

v0.15.0 で 116 スキル / 70 ToolUniverse 連携を達成。Phase 8 (v0.16) では ToolUniverse Tier 2 候補（OBIS 海洋生物分布・WoRMS 海洋分類・GBIF 生物多様性・MGnify メタゲノミクス）と K-Dense 残存ギャップ（signac scATAC-seq・rapids-singlecell GPU シングルセル）に加え、直接 API スキル（環境毒性学・rRNA 分類学・NCI-60 スクリーニング・植物バイオロジー・データ投稿）を追加した。

### 1.2 調査ソース

| リポジトリ | 規模 | 最終確認日 |
|-----------|------|-----------|
| **mims-harvard/ToolUniverse** | `default_config.py`: **188 カテゴリ**, 1229+ ツール | 2026-02-21 |
| **K-Dense-AI/claude-scientific-skills** | **140 スキル** (databases: 28+, packages: 55+, integrations: 15+, analysis: 30+) | 2026-02-21 |
| **SATORI v0.16.0** | **124 スキル**, **74 TU-linked**, 26 カテゴリ (A–Z) | 現在 |

### 1.3 TU キー名の調査結果 (Phase 8)

GAP_ANALYSIS_v0.15 の Tier 2 で列挙された TU キーの実在性を subagent で検証:

| 想定キー | 実際の TU キー | 状態 |
|---------|--------------|------|
| `toxnet`, `tox21`, `ctd` | **TU 不在** | 直接 API で実装 |
| `silva`, `greengenes2` | **TU 不在** → `mgnify` ✅ | MGnify のみ TU 連携 |
| `genbank_submission`, `sra_tools` | **TU 不在** | 直接 API で実装 |
| `cellminer`, `nci60` | **TU 不在** | 直接 API で実装 |
| `plant_reactome`, `tair` | **TU 不在** | 直接 API で実装 |
| `fish_base`, `ocean_biogeographic` | **TU 不在** → `obis` ✅, `worms` ✅, `gbif` ✅ | 3 キー TU 連携 |

→ 想定 TU キーの大部分は TU に存在しなかったが、代替として `obis`/`worms`/`gbif`/`mgnify` の 4 新規 TU キーを発見・連携。

---

## 2. v0.16.0 で追加した TU カテゴリ

| # | Config Key | SATORI スキル | ツール数 |
|---|-----------|--------------|---------|
| 136 | `mgnify` | rrna-taxonomy | ~5 |
| 137 | `obis` | marine-ecology | ~3 |
| 138 | `worms` | marine-ecology | ~3 |
| 139 | `gbif` | marine-ecology | ~3 |

### 累積 TU カバレッジ: ~139/188 categories (73.9%)

---

## 3. v0.16.0 で追加した K-Dense カバレッジ (+2)

| # | K-Dense Key | SATORI スキル | 種別 |
|---|------------|--------------|------|
| 134 | `signac` | scatac-signac | package |
| 135 | `rapids-singlecell` | gpu-singlecell | package |

### 累積 K-Dense カバレッジ: ~135/140 (96.4%)

---

## 4. 残存ギャップ: ToolUniverse 未カバーカテゴリ (~49 categories)

### 4.1 Tier 3 候補 (v0.17.0+、ニッチ / 低優先)

| Config Key(s) | 説明 | 難易度 |
|--------------|------|--------|
| `dictybase`, `plasmodb`, `vectorbase` | 稀少モデル生物 (寄生虫/アメーバ) | 低 |
| `soil_grids`, `worldclim` | 環境地理・気候データ | 中 |
| `icgc_argo`, `beacon_network` | がんデータ共有 | 低 |
| `geo_profiles` | GEO 発現プロファイル | 低 (EBI DB で部分カバー) |
| `encode`, `chipatlas` | エピゲノムアトラス | 中 |
| `hca_tools` | Human Cell Atlas | 中 |
| `paleobiology` | 古生物学データ | 低 |
| その他 ~30 カテゴリ | 極めて専門的 | — |

---

## 5. 残存ギャップ: K-Dense 未カバー (~5 skills)

| # | K-Dense Key | 種別 | 優先度 | 推奨対応 |
|---|------------|------|--------|---------|
| 1 | `metabolic-atlas` | database | 中 | metabolic-modeling 拡張 |
| 2 | `bio2bel` | package | 低 | ontology-enrichment 拡張 |
| 3 | `squidpy-advanced` | analysis | 低 | spatial-transcriptomics 拡張 |
| 4 | `starfysh` | package | 低 | spatial-transcriptomics 拡張 |
| 5 | `augur` | package | ✅ | perturbation-analysis に統合済み |

> 注: `augur` は pertpy 経由で既存スキルに統合済み。実質残り ~4 スキル。

---

## 6. カテゴリ別分布 (v0.16.0 時点)

| カテゴリ | 名称 | スキル数 | TU-linked | v0.15→v0.16 差分 |
|---------|------|---------|-----------|-----------------|
| A | 基盤・ワークフロー | **16** | 2 | **+1** |
| B | 統計・探索的解析 | 4 | 0 | — |
| C | 機械学習・モデリング | 3 | 0 | — |
| D | 実験計画・プロセス最適化 | 2 | 0 | — |
| E | 信号・スペクトル・時系列 | 4 | 0 | — |
| F | 生命科学・オミクス | 18 | 15 | — |
| G | 化学・材料・イメージング | 8 | 2 | — |
| H | 臨床・疫学・メタ科学 | 5 | 2 | — |
| I | Deep Research・文献検索 | 3 | 2 | — |
| J | 創薬・ファーマコロジー | **7** | 5 | **+1** |
| K | 構造生物学・タンパク質工学 | 5 | 5 | — |
| L | 精密医療・臨床意思決定 | 3 | 3 | — |
| M | 実験室自動化・データ管理 | 2 | 0 | — |
| N | 科学プレゼンテーション | 2 | 0 | — |
| O | 研究計画・グラント・規制 | 3 | 1 | — |
| P | ファーマコビジランス | 2 | 2 | — |
| Q | 腫瘍学・疾患研究 | 5 | 4 | — |
| R | 量子・先端計算 | 7 | 0 | — |
| S | 医用イメージング | 1 | 0 | — |
| T | シングルセル・空間・エピゲノミクス | **8** | 3 | **+2** |
| U | 免疫・感染症 | 2 | 2 | — |
| V | マイクロバイオーム・環境 | **6** | **2** | **+3** |
| W | システム生物学 | 2 | 2 | — |
| X | 疫学・公衆衛生 | **3** | 2 | **+1** |
| Y | 集団遺伝学 | 1 | 1 | — |
| Z | 科学テキストマイニング | 2 | 1 | — |
| **合計** | | **124** | **74** | **+8** |

---

## 7. v0.17.0+ ロードマップ候補

### Phase 9 (v0.17.0): TU Tier 3 + K-Dense 残存 — 目標 124→130+ スキル

| # | 推奨スキル名 | カテゴリ | TU keys | K-Dense |
|---|-------------|---------|---------|----------|
| 1 | scientific-encode-chipatlas | T | encode, chipatlas | — |
| 2 | scientific-human-cell-atlas | T | hca_tools | — |
| 3 | scientific-environmental-geodata | V | soil_grids, worldclim | — |
| 4 | scientific-metabolic-atlas | W | — | metabolic-atlas |
| 5 | scientific-squidpy-advanced | T | — | squidpy-advanced |
| 6 | scientific-parasite-genomics | F | dictybase, plasmodb | — |

### 推定達成指標 (v0.17.0 後)
- スキル: ~130
- TU カバレッジ: ~145/188 (77.1%)
- K-Dense カバレッジ: ~137/140 (97.9%)

---

## 8. まとめ

v0.16.0 では以下を達成した:

1. **8 新スキル** (116→124): 環境毒性学、rRNA 分類学、データ投稿、NCI-60 スクリーニング、植物バイオロジー、海洋生態学、scATAC-seq/Signac、GPU シングルセル
2. **TU 連携 70→74** (+4 スキル: mgnify, obis, worms, gbif)
3. **K-Dense 133→135** (+2: signac, rapids-singlecell)
4. **TU キー実在性の精査** — GAP_ANALYSIS_v0.15 の Tier 2 想定キー 12 個中 8 個が TU に不在。代わりに obis/worms/gbif/mgnify の実在キーを発見して連携

残り TU ~49 カテゴリ / K-Dense ~4 スキルのギャップは v0.17.0 以降で段階的に解消予定。
