# GAP ANALYSIS — SATORI v0.25.0

## 現状サマリー

| 指標 | 値 |
|---|---|
| 総スキル数 | 190 |
| TU 連携スキル | 131 |
| TU カバレッジ | 68.9% (131/190) |
| カテゴリ数 | 26 (A-Z) |
| Phase | 17 完了 |

## v0.25.0 追加分 (Phase 17)

| # | スキル | カテゴリ | TU |
|---|---|---|---|
| 183 | scientific-federated-learning | R | — |
| 184 | scientific-neural-architecture-search | R | — |
| 185 | scientific-semi-supervised-learning | C | — |
| 186 | scientific-multi-task-learning | C | — |
| 187 | scientific-statistical-simulation | B | — |
| 188 | scientific-streaming-analytics | B | — |
| 189 | scientific-radiology-ai | S | — |
| 190 | scientific-adaptive-experiments | D | — |

## カテゴリ別分布

| カテゴリ | スキル数 | TU連携 | TU率 |
|---|:---:|:---:|:---:|
| A. 基盤・ワークフロー | 17 | 11 | 64.7% |
| B. 統計・探索的解析 | 11 | 0 | 0.0% |
| C. 機械学習・モデリング | 11 | 0 | 0.0% |
| D. 実験計画・プロセス最適化 | 3 | 2 | 66.7% |
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
| R. 量子・先端計算 | 11 | 0 | 0.0% |
| S. 医用イメージング | 2 | 1 | 50.0% |
| T. シングルセル・空間・エピゲノミクス | 13 | 13 | 100% |
| U. 免疫・感染症 | 2 | 2 | 100% |
| V. マイクロバイオーム・環境 | 9 | 9 | 100% |
| W. システム生物学 | 4 | 4 | 100% |
| X. 疫学・公衆衛生 | 3 | 3 | 100% |
| Y. 集団遺伝学 | 2 | 2 | 100% |
| Z. 科学テキストマイニング | 3 | 3 | 100% |

## Phase 18 候補 (v0.26.0)

### CRITICAL: TU カバレッジ回復 (68.9% → 目標 75%+)

v0.23.0〜v0.25.0 で 24 スキルを追加したが全て TU 非連携 → カバレッジが 78.9% → 68.9% に低下。
v0.26.0 では **TU key 追加のみ** に集中し、カバレッジ回復を最優先とする。

#### Track A: 既存スキルへの TU key 追加候補 (12+ 件)

| スキル | カテゴリ | TU key 候補 | 備考 |
|---|---|---|---|
| scientific-eda-correlation | B | biotools, omictools | EDA ツール連携 |
| scientific-statistical-testing | B | biotools | 統計検定 BI ツール |
| scientific-pca-tsne | B | biotools | 次元削減 |
| scientific-ml-regression | C | openml | ML ベンチマーク |
| scientific-ml-classification | C | openml | ML ベンチマーク |
| scientific-feature-importance | C | openml | 特徴量解析 |
| scientific-deep-learning | R | papers_with_code | 深層学習 |
| scientific-bayesian-statistics | R | biotools | ベイズ統計 |
| scientific-explainable-ai | R | biotools | XAI |
| scientific-presentation-design | N | biotools | 可視化 |
| scientific-medical-imaging | S → 既に TU | — | — |
| scientific-radiology-ai | S | tcia | The Cancer Imaging Archive |

仮に 12 件の TU key を追加すれば: 131+12=143, 143/190 = **75.3%** → 目標達成

### MEDIUM: 薄いカテゴリ拡充

| カテゴリ | 現在 | 候補スキル |
|---|:---:|---|
| U. 免疫・感染症 | 2 | tcr-repertoire, vaccine-informatics |
| Y. 集団遺伝学 | 2 | ancient-dna, ld-analysis |

### LOW: 新ドメイン追加

- 農業・食品科学 (AgriTech)
- 宇宙科学・天文学 (AstroInformatics)
- 教育データマイニング (EDM)

## TU カバレッジ推移

| Version | Skills | TU | Coverage |
|---|:---:|:---:|:---:|
| v0.20.0 | 139 | 108 | 77.7% |
| v0.21.0 | 160 | 124 | 77.5% |
| v0.22.0 | 166 | 131 | 78.9% |
| v0.23.0 | 174 | 131 | 75.3% |
| v0.24.0 | 182 | 131 | 72.0% |
| **v0.25.0** | **190** | **131** | **68.9%** |

> **v0.26.0 では新スキル追加を凍結し、TU key 追加によるカバレッジ回復に専念する。
> 目標: 131 → 143+ TU 連携 (75%+ 回復)**
