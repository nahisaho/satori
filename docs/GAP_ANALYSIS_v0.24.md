# GAP ANALYSIS — SATORI v0.24.0

## 現状サマリー

| 指標 | 値 |
|---|---|
| 総スキル数 | 182 |
| TU 連携スキル | 131 |
| TU カバレッジ | 72.0% (131/182) |
| カテゴリ数 | 26 (A-Z) |
| Phase | 16 完了 |

## v0.24.0 追加分 (Phase 16)

| # | スキル | カテゴリ | TU |
|---|---|---|---|
| 175 | scientific-anomaly-detection | C | — |
| 176 | scientific-causal-ml | C | — |
| 177 | scientific-model-monitoring | C | — |
| 178 | scientific-time-series-forecasting | E | — |
| 179 | scientific-data-profiling | B | — |
| 180 | scientific-geospatial-analysis | B | — |
| 181 | scientific-network-visualization | B | — |
| 182 | scientific-reproducible-reporting | N | — |

## カテゴリ別分布

| カテゴリ | スキル数 | TU連携 | TU率 |
|---|:---:|:---:|:---:|
| A. 基盤・ワークフロー | 17 | 11 | 64.7% |
| B. 統計・探索的解析 | 9 | 0 | 0.0% |
| C. 機械学習・モデリング | 9 | 0 | 0.0% |
| D. 実験計画・プロセス最適化 | 2 | 2 | 100% |
| E. 信号・スペクトル・時系列 | 5 | 3 | 60.0% |
| F. 生命科学・オミクス | 28 | 27 | 96.4% |
| G. 化学・材料・イメージング | 9 | 8 | 88.9% |
| H. 臨床・疫学・メタ科学 | 7 | 6 | 85.7% |
| I. Deep Research・文献検索 | 4 | 4 | 100% |
| J. 創薬・ファーマコロジー | 9 | 9 | 100% |
| K. 構造生物学・タンパク質工学 | 7 | 7 | 100% |
| L. 精密医療・臨床意思決定 | 6 | 6 | 100% |
| M. 実験室自動化・データ管理 | 3 | 2 | 66.7% |
| N. 科学プレゼンテーション・図式 | 4 | 0 | 0.0% |
| O. 研究計画・グラント・規制 | 3 | 3 | 100% |
| P. ファーマコビジランス・薬理ゲノミクス | 4 | 4 | 100% |
| Q. 腫瘍学・疾患研究 | 10 | 10 | 100% |
| R. 量子・先端計算 | 9 | 0 | 0.0% |
| S. 医用イメージング | 1 | 1 | 100% |
| T. シングルセル・空間・エピゲノミクス | 13 | 13 | 100% |
| U. 免疫・感染症 | 2 | 2 | 100% |
| V. マイクロバイオーム・環境 | 9 | 9 | 100% |
| W. システム生物学 | 4 | 4 | 100% |
| X. 疫学・公衆衛生 | 3 | 3 | 100% |
| Y. 集団遺伝学 | 2 | 2 | 100% |
| Z. 科学テキストマイニング | 3 | 3 | 100% |

## Phase 17 候補 (v0.25.0)

### HIGH: TU カバレッジ回復 (72.0% → 目標 75%+)

TU 連携 0% のカテゴリに TU key を追加:
- **B (0/9 TU)**: eda-correlation/statistical-testing/pca-tsne に TU key 追加候補 — biotools, omictools
- **C (0/9 TU)**: ml-regression/ml-classification に TU key 追加候補 — openml, mlflow
- **N (0/4 TU)**: presentation-design に TU key 追加候補 — bio.tools 可視化ツール
- **R (0/9 TU)**: deep-learning/bayesian-statistics に TU key 追加候補 — DagsHub, Weights & Biases

### MEDIUM: 薄いカテゴリ拡充

| カテゴリ | 現在 | 候補スキル |
|---|:---:|---|
| S. 医用イメージング | 1 | radiology-ai, pathology-ai, ophthalmology-ai |
| D. 実験計画 | 2 | adaptive-experiments, ab-testing-design |
| U. 免疫・感染症 | 2 | tcr-repertoire, vaccine-informatics |
| Y. 集団遺伝学 | 2 | ancient-dna, ld-analysis |

### LOW: 新ドメイン追加

- 農業・食品科学 (AgriTech)
- 宇宙科学・天文学 (AstroInformatics)
- ロボティクス・自動制御 (RoboSci)
- 教育データマイニング (EDM)

## TU カバレッジ推移

| Version | Skills | TU | Coverage |
|---|:---:|:---:|:---:|
| v0.20.0 | 139 | 108 | 77.7% |
| v0.21.0 | 160 | 124 | 77.5% |
| v0.22.0 | 166 | 131 | 78.9% |
| v0.23.0 | 174 | 131 | 75.3% |
| **v0.24.0** | **182** | **131** | **72.0%** |

> v0.25.0 では TU key 追加によるカバレッジ回復を最優先とする。
