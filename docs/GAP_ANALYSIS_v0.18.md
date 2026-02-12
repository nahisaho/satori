# SATORI v0.18.0 GAP ANALYSIS

## 1. 現在の達成状況

| 指標 | v0.17.0 | v0.18.0 | 変化 |
|------|---------|---------|------|
| 総スキル数 | 132 | 140 | +8 |
| TU 連携スキル | 79 | 85 | +6 |
| 中区分カテゴリ数 | 26 (A-Z) | 26 (A-Z) | ±0 |
| K-Dense カバー率 | ~98.6% | ~98.6% | ±0 |

## 2. v0.18.0 で追加したスキル (Phase 10)

| # | スキル名 | カテゴリ | TU Key | K-Dense |
|---|---------|---------|--------|---------|
| 133 | scientific-gwas-catalog | Y | gwas | — |
| 134 | scientific-alphafold-structures | K | alphafold | — |
| 135 | scientific-arrayexpress-expression | F | arrayexpress | — |
| 136 | scientific-semantic-scholar | I | semantic_scholar | — |
| 137 | scientific-gtex-tissue-expression | F | (direct) | — |
| 138 | scientific-pharmgkb-pgx | P | pharmgkb | — |
| 139 | scientific-crossref-metadata | A | crossref | — |
| 140 | scientific-icgc-cancer-data | Q | (direct) | — |

## 3. TU 連携状況 (85/188 カテゴリ = 45.2%)

### 新規 TU Key 連携 (+6)

| TU Key | スキル | ツール数 (推定) |
|--------|-------|----------------|
| gwas | gwas-catalog | 3-5 |
| alphafold | alphafold-structures | 3-5 |
| arrayexpress | arrayexpress-expression | 3-5 |
| semantic_scholar | semantic-scholar | 3-5 |
| pharmgkb | pharmgkb-pgx | 3-5 |
| crossref | crossref-metadata | 3-5 |

### TU 残ギャップ (~103 カテゴリ未連携)

以下は TU に存在するが SATORI スキルが未連携のカテゴリ候補:

**高優先度 (既存スキルに TU Key 追加可能):**
- `civic` — precision-oncology に追加可能
- `orphanet` — rare-disease-genetics に追加可能
- `disgenet` — disease-research に追加可能
- `intact` — protein-interaction-network に追加可能
- `metacyc` — metabolomics-databases に追加可能
- `zinc` — compound-screening に追加可能
- `cellosaurus` — cell-line-resources に追加可能
- `hmdb` — metabolomics-databases に追加可能

**中優先度 (新スキル必要):**
- `rcsb_pdb` / `rcsb_search` — PDB 直接検索 API
- `uniprot` — UniProt REST API
- `chembl` — ChEMBL (既に chembl-assay-mining があるが TU Key 確認必要)
- `clinvar` — variant-interpretation に追加可能
- `gnomad` — variant-interpretation に追加可能

## 4. K-Dense カバー率分析

### 重要な発見 (v0.17.0 GAP ANALYSIS より)

v0.17.0 ギャップ分析で `bio2bel` と `starfysh` が K-Dense ギャップとして報告されていたが、**サブエージェント調査により K-Dense リポジトリに実際には存在しないことが判明**。

- `bio2bel` — K-Dense-AI/claude-scientific-skills に **存在しない**
- `starfysh` — K-Dense-AI/claude-scientific-skills に **存在しない**
- `augur` — pertpy/perturbation-analysis で既にカバー済み

**結論:** K-Dense ギャップは事実上 0。SATORI は K-Dense の全スキルをカバー。

## 5. カテゴリ別スキル分布 (v0.18.0)

| カテゴリ | スキル数 | TU 連携数 |
|---------|---------|----------|
| A. 基盤・ワークフロー | 17 | ~10 |
| B. 統計・探索的解析 | 4 | 0 |
| C. 機械学習・モデリング | 3 | 0 |
| D. 実験計画・プロセス最適化 | 2 | 0 |
| E. 信号・スペクトル・時系列 | 4 | 0 |
| F. 生命科学・オミクス | 22 | ~18 |
| G. 化学・材料・イメージング | 8 | ~4 |
| H. 臨床・疫学・メタ科学 | 5 | ~2 |
| I. Deep Research・文献検索 | 4 | ~3 |
| J. 創薬・ファーマコロジー | 7 | ~6 |
| K. 構造生物学・タンパク質工学 | 6 | ~5 |
| L. 精密医療・臨床意思決定 | 3 | ~3 |
| M. 実験室自動化・データ管理 | 2 | ~1 |
| N. 科学プレゼンテーション・図式 | 2 | 0 |
| O. 研究計画・グラント・規制 | 3 | ~1 |
| P. ファーマコビジランス・薬理ゲノミクス | 3 | ~2 |
| Q. 腫瘍学・疾患研究 | 6 | ~5 |
| R. 量子・先端計算 | 7 | 0 |
| S. 医用イメージング | 1 | 0 |
| T. シングルセル・空間・エピゲノミクス | 11 | ~8 |
| U. 免疫・感染症 | 2 | ~2 |
| V. マイクロバイオーム・環境 | 8 | ~6 |
| W. システム生物学 | 3 | ~2 |
| X. 疫学・公衆衛生 | 3 | ~2 |
| Y. 集団遺伝学 | 2 | ~1 |
| Z. 科学テキストマイニング | 2 | ~2 |
| **合計** | **140** | **~85** |

## 6. Phase 11 ロードマップ候補

### 優先度 A (既存スキル TU Key 追加)
既存スキルの frontmatter に TU key を追加することで、実装コスト最小で TU 連携数を増やせる。

1. `civic` → precision-oncology
2. `orphanet` → rare-disease-genetics
3. `disgenet` → disease-research
4. `intact` → protein-interaction-network
5. `metacyc` → metabolomics-databases
6. `zinc` → compound-screening
7. `clinvar` → variant-interpretation
8. `gnomad` → variant-interpretation

**推定効果:** TU 85→93 (+8)

### 優先度 B (新スキル候補)
1. scientific-uniprot-proteome (F) — UniProt REST API
2. scientific-rcsb-pdb-search (K) — RCSB PDB Search API
3. scientific-opentargets-genetics (Q) — Open Targets Genetics API
4. scientific-biobank-data (H) — UK Biobank / All of Us
5. scientific-metabolic-flux (W) — 13C-MFA / リアルタイムフラックス
6. scientific-spatial-multiomics (T) — MERFISH + Proteomics 統合

## 7. 結論

v0.18.0 は Phase 10 として 8 新スキルを追加し、以下を達成:
- 総スキル: 132→140 (+8)
- TU 連携: 79→85 (+6 新規 TU Key)
- K-Dense: 事実上 100% カバー (bio2bel/starfysh は K-Dense に不存在)
- 新ドメイン: GWAS Catalog、AlphaFold DB、ArrayExpress、Semantic Scholar、GTEx、PharmGKB、CrossRef、ICGC

次フェーズでは既存スキルへの TU Key 追加 (低コスト高効果) と新スキル追加を並行して進めることで、TU カバレッジ 85→93+ を目指す。
