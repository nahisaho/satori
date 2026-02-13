# GAP_ANALYSIS_v0.22.md — SATORI v0.22.0 ギャップ分析

## 現状サマリ (v0.22.0 リリース後)

| 指標 | 値 |
|---|---|
| 総スキル数 | 166 |
| TU 連携スキル | 131 |
| TU カバレッジ | 78.9% |
| カテゴリ数 | 26 (A-Z) |
| Phase | 14 完了 |

## Phase 14 実績 (v0.21.0 → v0.22.0)

### Track A: 既存スキル TU key 追加 (8 件)

| スキル | 追加 TU key | 新規 TU 連携 |
|---|---|---|
| epigenomics-chromatin | chipatlas | ★ |
| metabolomics | metabolomics_workbench | (既存 hmdb) |
| systems-biology | bigg_models, complex_portal, wikipathways | ★ |
| immunoinformatics | imgt, sabdab, therasabdab | (既存 iedb) |
| public-health-data | nhanes, medlineplus, odphp | ★ |
| epidemiology-public-health | who_gho | ★ |
| model-organism-db | impc, mpd | ★ |
| environmental-ecology | gbif | ★ |

### Track B: 新規スキル (6 件)

| # | スキル | カテゴリ | TU 連携 |
|---|---|---|---|
| 161 | scientific-glycomics | F | — |
| 162 | scientific-lipidomics | F | — |
| 163 | scientific-metagenome-assembled-genomes | V | — |
| 164 | scientific-crispr-design | M | — |
| 165 | scientific-clinical-pharmacology | P | — |
| 166 | scientific-clinical-standards | H | loinc, icd |

---

## Phase 15 候補 (v0.23.0 ロードマップ)

### Priority A: 既存スキル TU key 追加候補

| スキル | 候補 TU key | ツール数 |
|---|---|---|
| toxicology-env | ctdbase (CTD 未確認) | 要調査 |
| infectious-disease | card (CARD AMR 未確認) | 要調査 |
| noncoding-rna | rfam (Rfam RNA families) | 要調査 |
| ontology-enrichment | umls (UMLS メタシソーラス) | 要調査 |
| phylogenetics | timetree (TimeTree 分岐年代) | 要調査 |
| protein-design | rosetta (Rosetta 構造予測) | 要調査 |

### Priority B: 新規スキル候補

| 候補名 | カテゴリ | 概要 |
|---|---|---|
| scientific-proteomics-atlas | F | ProteomicsDB/PRIDE プロテオミクスアトラス |
| scientific-food-science | F | FooDB/FoodData Central 食品化学・栄養素 |
| scientific-exosome-ev | F | ExoCarta/Vesiclepedia 細胞外小胞カーゴ |
| scientific-single-cell-multimodal | T | CITE-seq/TEA-seq マルチモーダルシングルセル |
| scientific-long-read-assembly | V | HiFi/ONT ロングリードアセンブリ・Hifiasm |
| scientific-metaproteomics | F | メタプロテオミクス・コミュニティレベルタンパク質 |

### Priority C: パイプライン強化

| パイプライン | 説明 |
|---|---|
| metabolomics → lipidomics → glycomics → multi-omics | 代謝物・脂質・糖鎖のオミクス統合 |
| crispr-design → perturbation-analysis → functional-genomics | CRISPR 設計→摂動→機能ゲノミクス |
| clinical-standards → clinical-decision-support | 臨床標準用語→意思決定支援 |
| metagenome-assembled-genomes → environmental-ecology | MAG → 環境生態系統合 |
| clinical-pharmacology → pharmacogenomics → precision-oncology | 臨床薬理→ゲノム薬理→精密医療 |

---

## カテゴリ分布 (v0.22.0)

| カテゴリ | スキル数 |
|---|:---:|
| A. 基盤・ワークフロー | 17 |
| B. 統計・探索的解析 | 4 |
| C. 機械学習・モデリング | 3 |
| D. 実験計画・プロセス最適化 | 2 |
| E. 信号・スペクトル・時系列 | 4 |
| F. 生命科学・オミクス | 28 |
| G. 化学・材料・イメージング | 9 |
| H. 臨床・疫学・メタ科学 | 7 |
| I. Deep Research・文献検索 | 4 |
| J. 創薬・ファーマコロジー | 9 |
| K. 構造生物学・タンパク質工学 | 7 |
| L. 精密医療・臨床意思決定 | 6 |
| M. 実験室自動化・データ管理 | 3 |
| N. 科学プレゼンテーション・図式 | 2 |
| O. 研究計画・グラント・規制 | 3 |
| P. ファーマコビジランス・薬理ゲノミクス | 4 |
| Q. 腫瘍学・疾患研究 | 10 |
| R. 量子・先端計算 | 7 |
| S. 医用イメージング | 1 |
| T. シングルセル・空間・エピゲノミクス | 13 |
| U. 免疫・感染症 | 2 |
| V. マイクロバイオーム・環境 | 9 |
| W. システム生物学 | 4 |
| X. 疫学・公衆衛生 | 3 |
| Y. 集団遺伝学 | 2 |
| Z. 科学テキストマイニング | 3 |
| **合計** | **166** |

## TU 連携推移

| バージョン | スキル | TU 連携 | カバレッジ |
|---|:---:|:---:|:---:|
| v0.10.0 | 72 | 22 | 30.6% |
| v0.11.0 | 80 | 30 | 37.5% |
| v0.12.0 | 88 | 37 | 42.0% |
| v0.13.0 | 95 | 44 | 46.3% |
| v0.14.0 | 103 | 50 | 48.5% |
| v0.15.0 | 110 | 57 | 51.8% |
| v0.16.0 | 118 | 68 | 57.6% |
| v0.17.0 | 126 | 82 | 65.1% |
| v0.18.0 | 134 | 92 | 68.7% |
| v0.19.0 | 142 | 100 | 70.4% |
| v0.20.0 | 154 | 114 | 74.0% |
| v0.21.0 | 160 | 124 | 77.5% |
| **v0.22.0** | **166** | **131** | **78.9%** |
