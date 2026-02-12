# SATORI v0.20.0 GAP ANALYSIS

## 現状サマリー

| 指標 | v0.19.0 | v0.20.0 | 差分 |
|------|---------|---------|------|
| 総スキル数 | 148 | 154 | +6 |
| TU 連携スキル | 99 | 114 | +15 (3 new + 12 existing key additions) |
| K-Dense スキル | 148 | 154 | +6 |
| カテゴリ数 (A-Z) | 26 | 26 | ±0 |

## Phase 12 実施内容

### Track A: 既存スキル TU key 追加 (12 件)

| スキル | 追加 TU Key | 説明 |
|--------|------------|------|
| scientific-cancer-genomics | cosmic, cbioportal | がん体細胞変異カタログ、がんゲノミクスポータル |
| scientific-precision-oncology | oncokb | 精密腫瘍学アノテーション |
| scientific-immunoinformatics | iedb | 免疫エピトープデータベース |
| scientific-microbiome-metagenomics | mgnify | EBI メタゲノミクス解析プラットフォーム |
| scientific-cell-line-resources | cellosaurus | 細胞株データベース (ExPASy) |
| scientific-string-network-api | ppi | STRING/BioGRID PPI ネットワーク |
| scientific-metabolic-modeling | biomodels | SBML モデルリポジトリ (EBI) |
| scientific-chembl-assay-mining | chembl | 創薬生理活性データベース (EBI) |
| scientific-pharmacology-targets | bindingdb, gtopdb, brenda | 結合親和性・GtoPdb・酵素動態 |
| scientific-drug-target-profiling | dgidb | 薬物-遺伝子相互作用データベース |
| scientific-admet-pharmacokinetics | pubchem | 化合物・物質・生理活性アッセイ |
| scientific-metabolomics | hmdb | ヒトメタボロームデータベース |

### Track B: 新規スキル (6 件、うち 3 件 TU 連携)

| # | スキル名 | Cat | TU Key | 説明 |
|---|----------|-----|--------|------|
| 149 | scientific-monarch-ontology | Q | monarch | Monarch Initiative 疾患-遺伝子-表現型 |
| 150 | scientific-gdc-portal | Q | gdc | NCI Genomic Data Commons REST API |
| 151 | scientific-biobank-cohort | H | — | UK Biobank/BBJ/AoU 大規模コホート |
| 152 | scientific-spatial-multiomics | T | — | MERFISH/CODEX 空間マルチオミクス統合 |
| 153 | scientific-metabolic-flux | W | — | 13C/15N 安定同位体フラックス解析 |
| 154 | scientific-stitch-chemical-network | G | stitch | STITCH 化学物質-タンパク質相互作用 |

## カテゴリ別スキル数 (v0.20.0)

| Cat | 名前 | v0.19 | v0.20 |
|-----|------|-------|-------|
| A | 基盤・ワークフロー | 17 | 17 |
| B | 統計・探索的解析 | 4 | 4 |
| C | 機械学習・モデリング | 3 | 3 |
| D | 実験計画・プロセス最適化 | 2 | 2 |
| E | 信号・スペクトル・時系列 | 4 | 4 |
| F | 生命科学・オミクス | 24 | 24 |
| G | 化学・材料・イメージング | 8 | **9** (+1) |
| H | 臨床・疫学・メタ科学 | 5 | **6** (+1) |
| I | Deep Research・文献検索 | 4 | 4 |
| J | 創薬・ファーマコロジー | 8 | 8 |
| K | 構造生物学・タンパク質工学 | 7 | 7 |
| L | 精密医療・臨床意思決定 | 5 | 5 |
| M | 実験室自動化・データ管理 | 2 | 2 |
| N | 科学プレゼンテーション・図式 | 2 | 2 |
| O | 研究計画・グラント・規制 | 3 | 3 |
| P | ファーマコビジランス・薬理ゲノミクス | 3 | 3 |
| Q | 腫瘍学・疾患研究 | 8 | **10** (+2) |
| R | 量子・先端計算 | 7 | 7 |
| S | 医用イメージング | 1 | 1 |
| T | シングルセル・空間・エピゲノミクス | 11 | **12** (+1) |
| U | 免疫・感染症 | 2 | 2 |
| V | マイクロバイオーム・環境 | 8 | 8 |
| W | システム生物学 | 3 | **4** (+1) |
| X | 疫学・公衆衛生 | 3 | 3 |
| Y | 集団遺伝学 | 2 | 2 |
| Z | 科学テキストマイニング | 2 | 2 |
| **合計** | | **148** | **154** |

## TU Key カバレッジ

### v0.20.0 で追加された TU Keys (15 件)

**新規スキル経由 (3):** monarch, gdc, stitch

**既存スキル追加 (12/17 unique keys):** cosmic, cbioportal, oncokb, iedb, mgnify, cellosaurus, ppi, biomodels, chembl, bindingdb, gtopdb, brenda, dgidb, pubchem, hmdb

### TU Key 非対応 (ToolUniverse 未収録)
- CARD (AMR データベース): ToolUniverse に登録なし → 将来候補
- HGNC (遺伝子命名法): ToolUniverse に登録なし → ensembl/biothings で代替

### K-Dense カバレッジ

- K-Dense 全 148 基盤スキル → v0.20.0 で 154 にスケール
- カバレッジ: 100%

## Phase 13 候補 (v0.21.0 以降)

### Priority A: 追加 TU key 候補 (既存スキルへの追加)

| スキル | 追加候補 TU Key | 理由 |
|--------|----------------|------|
| scientific-human-cell-atlas | cellxgene_census | CELLxGENE Census ~7 tools 追加 |
| scientific-drug-repurposing | pharos | Pharos/TCRD IDG ターゲットナレッジベース |
| scientific-pharmacogenomics | fda_pharmacogenomic_biomarkers | FDA PGx バイオマーカー |
| scientific-infectious-disease | card (if added to TU) | AMR 耐性遺伝子データベース |

### Priority B: 新規スキル候補

| 候補名 | Cat | TU Key | 説明 |
|--------|-----|--------|------|
| scientific-cellxgene-census | T | cellxgene_census | CELLxGENE Census 大規模セルアトラス API |
| scientific-pharos-targets | J | pharos | Pharos/TCRD IDG ターゲットプロファイリング |
| scientific-prosite-motifs | F | interpro | PROSITE/InterPro モチーフ検索拡張 |
| scientific-hgnc-nomenclature | F | — | HGNC 遺伝子命名法・ID バリデーション |
| scientific-metabolomics-network | F | — | 代謝物ネットワーク構築・MetaboAnalyst 統合 |
| scientific-clinical-nlp | Z | — | 臨床テキスト NER (MedSpaCy/cTAKES) |

### Priority C: パイプライン強化

- Monarch → GDC → precision-oncology 疾患-がんゲノムパイプライン
- STITCH → STRING → drug-target-profiling 化学-PPI 統合ネットワーク
- Biobank-cohort → population-genetics → GWAS Catalog PheWAS 連携
- Cancer-genomics (COSMIC/cBioPortal) → GDC → CIViC 変異エビデンスチェーン

## 結論

v0.20.0 で Phase 12 を完了。ハイブリッド戦略 (Track A: 既存スキル TU key 大量追加 +
Track B: 新規スキル) により、TU 連携は 99→114 へ拡大。特に Track A で 12 の既存スキルに
17 の TU key を一括追加したことで、既存スキルの外部ツール活用能力が大幅に向上した。

新規 6 スキルは、疾患オントロジー (Monarch)、がんゲノムデータポータル (GDC)、
大規模コホート (Biobank)、空間マルチオミクス、代謝フラックス、化学-タンパク質
ネットワーク (STITCH) と、多様な領域をカバー。

次バージョン v0.21.0 では、CELLxGENE Census (~7 tools)、Pharos/TCRD、
FDA PGx バイオマーカー等の新規 TU key 対応と、臨床 NLP・代謝物ネットワーク等の
新規スキル候補を検討する。
