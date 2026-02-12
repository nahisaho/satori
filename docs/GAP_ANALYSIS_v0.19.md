# SATORI v0.19.0 GAP ANALYSIS

## 現状サマリー

| 指標 | v0.18.0 | v0.19.0 | 差分 |
|------|---------|---------|------|
| 総スキル数 | 140 | 148 | +8 |
| TU 連携スキル | 85 | 99 | +14 (8 new + 6 existing key additions) |
| K-Dense スキル | 140 | 148 | +8 |
| カテゴリ数 (A-Z) | 26 | 26 | ±0 |

## Phase 11 実施内容

### Track A: 既存スキル TU key 追加 (6 件)

| スキル | 追加 TU Key | ツール数 (概算) |
|--------|------------|----------------|
| scientific-rare-disease-genetics | orphanet | ~5 |
| scientific-disease-research | disgenet | ~5 |
| scientific-protein-interaction-network | intact | ~3+ |
| scientific-metabolomics-databases | metacyc | ~4 |
| scientific-compound-screening | zinc | ~4 |
| scientific-variant-interpretation | clinvar | ~3 |

### Track B: 新規 TU 連携スキル (8 件)

| # | スキル名 | Cat | TU Key | ツール数 |
|---|----------|-----|--------|---------|
| 141 | scientific-uniprot-proteome | F | uniprot | ~18 |
| 142 | scientific-rcsb-pdb-search | K | rcsb_pdb, rcsb_search | ~10+ |
| 143 | scientific-opentargets-genetics | Q | opentarget | ~55 |
| 144 | scientific-reactome-pathways | F | reactome | ~3+ |
| 145 | scientific-depmap-dependencies | Q | depmap | ~3+ |
| 146 | scientific-drugbank-resources | J | drugbank | ~3+ (PARTIAL) |
| 147 | scientific-civic-evidence | L | civic | ~12 |
| 148 | scientific-gnomad-variants | L | gnomad | ~7 |

## カテゴリ別スキル数 (v0.19.0)

| Cat | 名前 | v0.18 | v0.19 |
|-----|------|-------|-------|
| A | 基盤・ワークフロー | 17 | 17 |
| B | 統計・探索的解析 | 4 | 4 |
| C | 機械学習・モデリング | 3 | 3 |
| D | 実験計画・プロセス最適化 | 2 | 2 |
| E | 信号・スペクトル・時系列 | 4 | 4 |
| F | 生命科学・オミクス | 22 | **24** (+2) |
| G | 化学・材料・イメージング | 8 | 8 |
| H | 臨床・疫学・メタ科学 | 5 | 5 |
| I | Deep Research・文献検索 | 4 | 4 |
| J | 創薬・ファーマコロジー | 7 | **8** (+1) |
| K | 構造生物学・タンパク質工学 | 6 | **7** (+1) |
| L | 精密医療・臨床意思決定 | 3 | **5** (+2) |
| M | 実験室自動化・データ管理 | 2 | 2 |
| N | 科学プレゼンテーション・図式 | 2 | 2 |
| O | 研究計画・グラント・規制 | 3 | 3 |
| P | ファーマコビジランス・薬理ゲノミクス | 3 | 3 |
| Q | 腫瘍学・疾患研究 | 6 | **8** (+2) |
| R | 量子・先端計算 | 7 | 7 |
| S | 医用イメージング | 1 | 1 |
| T | シングルセル・空間・エピゲノミクス | 11 | 11 |
| U | 免疫・感染症 | 2 | 2 |
| V | マイクロバイオーム・環境 | 8 | 8 |
| W | システム生物学 | 3 | 3 |
| X | 疫学・公衆衛生 | 3 | 3 |
| Y | 集団遺伝学 | 2 | 2 |
| Z | 科学テキストマイニング | 2 | 2 |
| **合計** | | **140** | **148** |

## TU Key カバレッジ

### v0.19.0 で追加された TU Keys (14 件)

**新規スキル経由 (8):** uniprot, rcsb_pdb, rcsb_search, opentarget, reactome, depmap, drugbank, civic, gnomad

**既存スキル追加 (6):** orphanet, disgenet, intact, metacyc, zinc, clinvar

### K-Dense カバレッジ

- K-Dense 全 140 基盤スキル → v0.19.0 で 148 にスケール
- カバレッジ: 100%（bio2bel/starfysh は K-Dense 未収録のため除外）

## Phase 12 候補 (v0.20.0 以降)

### Priority A: 追加 TU key 候補 (既存スキルへの追加)

| スキル | 追加候補 TU Key | 理由 |
|--------|----------------|------|
| scientific-cancer-genomics | cosmic, cbioportal | COSMIC/cBioPortal 直接 API ツール |
| scientific-precision-oncology | oncokb | OncoKB アノテーション |
| scientific-immunoinformatics | iedb | IEDB エピトープ予測 |
| scientific-infectious-disease | card | AMR データベース |
| scientific-microbiome-metagenomics | mgnify | MGnify メタゲノム |

### Priority B: 新規スキル候補

| 候補名 | Cat | TU Key | 説明 |
|--------|-----|--------|------|
| scientific-biobank-data | Q | N/A | UK Biobank/BBJ 大規模コホートデータ解析 |
| scientific-spatial-multiomics | T | N/A | 空間マルチオミクス統合 (MERFISH+proteomics) |
| scientific-metabolic-flux | W | N/A | 13C/15N フラックス解析・パスウェイフロー |
| scientific-prosite-motifs | F | interpro | PROSITE/InterPro モチーフ検索 |
| scientific-monarch-ontology | Q | monarch | Monarch Initiative 疾患-表現型オントロジー |
| scientific-hgnc-nomenclature | F | ensembl/biothings | HGNC 遺伝子命名法・ID 統合 |

### Priority C: パイプライン強化

- CIViC → gnomAD → Open Targets 3段パイプライン統合テスト
- UniProt → RCSB PDB → AlphaFold 構造パイプライン連結
- DepMap + Open Targets 薬剤標的ダブルスコアリング

## 結論

v0.19.0 で Phase 11 を完了。ToolUniverse 連携は 85→99 へ大幅拡大し、
特に大規模 API (Open Targets ~55 tools, UniProt ~18 tools, CIViC ~12 tools,
RCSB PDB ~10+ tools) の統合により外部ツールカバレッジが飛躍的に向上した。

次バージョン v0.20.0 では、既存スキルへの追加 TU key (COSMIC, OncoKB, IEDB 等)
と新規スキル候補 (Biobank, Spatial-Multiomics, Metabolic-Flux 等) を検討する。
