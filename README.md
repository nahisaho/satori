# SATORI（悟り）— Agent Skills for Science

**SATORI** は、科学データ解析のための **GitHub Copilot Agent Skills** コレクションです。

## Overview

このディレクトリには、Exp-01〜13 で蓄積した科学データ解析技法を Agent Skills として体系化した **27 個**のスキルを格納しています。Copilot がプロンプトの文脈に応じて適切なスキルを自動ロードし、各実験で確立した解析パターンを再利用します。

スキルは **8 つの中区分**に分類されています。

| 中区分 | スキル数 | 概要 |
|---|:---:|---|
| A. 基盤・ワークフロー | 5 | パイプライン構築・前処理・データ生成・論文図表・学術論文執筆 |
| B. 統計・探索的解析 | 3 | EDA・仮説検定・次元削減 |
| C. 機械学習・モデリング | 3 | 回帰・分類・特徴量重要度 |
| D. 実験計画・プロセス最適化 | 2 | DOE・応答曲面法・ベイズ最適化 |
| E. 信号・スペクトル・時系列 | 3 | スペクトル解析・生体信号・時系列分解 |
| F. 生命科学・オミクス | 5 | バイオインフォ・メタボロ・ゲノム配列・マルチオミクス・ネットワーク |
| G. 化学・材料・イメージング | 3 | ケモインフォ・材料特性評価・画像解析 |
| H. 臨床・疫学・メタ科学 | 3 | 生存解析・因果推論・メタアナリシス |

---

## Skills 一覧

### A. 基盤・ワークフロー（5 種）

全 Exp に共通する横断的な基盤スキル。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 1 | [scientific-pipeline-scaffold](scientific-pipeline-scaffold/SKILL.md) | パイプライン雛形・ディレクトリ構造・StepLogger・JSON 要約 | 全 Exp |
| 2 | [scientific-data-preprocessing](scientific-data-preprocessing/SKILL.md) | 欠損値補完・エンコーディング・スケーリング・外れ値処理 | 全 Exp |
| 3 | [scientific-data-simulation](scientific-data-simulation/SKILL.md) | 物理/化学/生物ベースの合成データ生成 | 06-09, 12, 13 |
| 4 | [scientific-publication-figures](scientific-publication-figures/SKILL.md) | 論文品質図表・rcParams・カラーパレット・マルチパネル | 10, 11-13 |
| 5 | [scientific-academic-writing](scientific-academic-writing/SKILL.md) | 学術論文執筆・ジャーナル別テンプレート・Cover Letter・査読対応 | 汎用 |

### B. 統計・探索的解析（3 種）

データの理解・検定・次元削減を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 6 | [scientific-eda-correlation](scientific-eda-correlation/SKILL.md) | 探索的データ解析・相関ヒートマップ・分布可視化 | 02, 12, 13 |
| 7 | [scientific-statistical-testing](scientific-statistical-testing/SKILL.md) | 仮説検定・多重比較・エンリッチメント・ベイズ推論 | 03, 04, 06, 07 |
| 8 | [scientific-pca-tsne](scientific-pca-tsne/SKILL.md) | PCA / t-SNE / UMAP 次元削減・クラスタリング | 02, 03, 07, 11, 13 |

### C. 機械学習・モデリング（3 種）

教師あり学習と特徴量解釈を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 9 | [scientific-ml-regression](scientific-ml-regression/SKILL.md) | マルチターゲット回帰・モデル比較・レーダーチャート | 05, 12, 13 |
| 10 | [scientific-ml-classification](scientific-ml-classification/SKILL.md) | 分類 ML・ROC・PR 曲線・混同行列・PDP・Volcano | 03, 05 |
| 11 | [scientific-feature-importance](scientific-feature-importance/SKILL.md) | Tree-based & Permutation 特徴量重要度・PDP | 05, 12, 13 |

### D. 実験計画・プロセス最適化（2 種）

実験設計と応答曲面最適化を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 12 | [scientific-doe](scientific-doe/SKILL.md) | 田口直交表・CCD/Box-Behnken・ANOVA 因子効果・ベイズ最適化 | 汎用 |
| 13 | [scientific-process-optimization](scientific-process-optimization/SKILL.md) | 応答曲面法 (ML-RSM)・パレート最適化・プロセスウィンドウ | 12, 13 |

### E. 信号・スペクトル・時系列（3 種）

波形・周波数領域の解析を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 14 | [scientific-spectral-signal](scientific-spectral-signal/SKILL.md) | スペクトル前処理・フィルタリング・ピーク検出 | 11 |
| 15 | [scientific-biosignal-processing](scientific-biosignal-processing/SKILL.md) | ECG R波/HRV・EEG バンドパワー/ERP・EMG バースト・Poincaré | 08 |
| 16 | [scientific-time-series](scientific-time-series/SKILL.md) | STL 分解・SARIMA 予測・変化点検出・FFT 周期解析・Granger 因果 | 汎用 |

### F. 生命科学・オミクス（5 種）

バイオ・オミクス・ネットワーク解析を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 17 | [scientific-bioinformatics](scientific-bioinformatics/SKILL.md) | scRNA-seq・PPI ネットワーク・バルク RNA-seq | 01, 04 |
| 18 | [scientific-metabolomics](scientific-metabolomics/SKILL.md) | PLS-DA/VIP スコア・Pareto スケーリング・パスウェイ濃縮 | 07 |
| 19 | [scientific-sequence-analysis](scientific-sequence-analysis/SKILL.md) | RSCU/CAI コドン解析・アラインメント・系統樹・ORF/CpG 島 | 09 |
| 20 | [scientific-multi-omics](scientific-multi-omics/SKILL.md) | CCA 正準相関・SNF ネットワーク融合・パスウェイ統合・マルチオミクスクラスタ | 汎用 |
| 21 | [scientific-network-analysis](scientific-network-analysis/SKILL.md) | ネットワーク構築・中心性・コミュニティ・PSP パス図 | 04, 07, 13 |

### G. 化学・材料・イメージング（3 種）

化学構造・材料特性評価・画像形態解析を担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 22 | [scientific-cheminformatics](scientific-cheminformatics/SKILL.md) | RDKit 分子記述子・Tanimoto・構造アラート・Lipinski | 02, 05 |
| 23 | [scientific-materials-characterization](scientific-materials-characterization/SKILL.md) | Thornton-Anders SZM・XRD Scherrer・Tauc プロット | 11, 12, 13 |
| 24 | [scientific-image-analysis](scientific-image-analysis/SKILL.md) | Otsu/Watershed セグメンテーション・粒径分布・GLCM テクスチャ・蛍光合成 | 汎用 |

### H. 臨床・疫学・メタ科学（3 種）

臨床試験・因果推論・メタアナリシスを担うスキル群。

| # | Skill | 説明 | 参照 Exp |
|---|---|---|---|
| 25 | [scientific-survival-clinical](scientific-survival-clinical/SKILL.md) | Kaplan-Meier・Cox PH・検出力分析・安全性解析 | 03, 06 |
| 26 | [scientific-causal-inference](scientific-causal-inference/SKILL.md) | PSM 傾向スコア・IPW・DID・RDD・DAG 共変量選択・Rosenbaum 感度分析 | 汎用 |
| 27 | [scientific-meta-analysis](scientific-meta-analysis/SKILL.md) | 固定/ランダム効果モデル・Forest/Funnel プロット・Egger 検定・サブグループ | 汎用 |

---

## 使い方

### GitHub Copilot Agent Mode / Copilot CLI での利用

Skills は `.github/skills/` に配置されているため、Copilot が自動的に検出します。
プロンプトの文脈に応じて、関連する Skill の `SKILL.md` がエージェントに注入されます。

```
# 例：相関解析を依頼すると scientific-eda-correlation が自動ロードされる
> ZnO 薄膜のプロセスパラメータと膜物性の相関ヒートマップを作成して

# 例：分類モデルを依頼すると scientific-ml-classification が自動ロードされる
> がん遺伝子発現データで Random Forest と SVM の ROC 比較をして

# 例：実験計画法を依頼すると scientific-doe が自動ロードされる
> 3因子の田口L9直交表を作成して主効果プロットを描画して

# 例：メタアナリシスを依頼すると scientific-meta-analysis が自動ロードされる
> 5本のRCT論文からランダム効果モデルでForestプロットを作成して
```

### ディレクトリ構造

```
.github/skills/
├── README.md
│
│── [A] 基盤・ワークフロー
│   ├── scientific-pipeline-scaffold/
│   ├── scientific-data-preprocessing/
│   ├── scientific-data-simulation/
│   ├── scientific-publication-figures/
│   └── scientific-academic-writing/
│       └── assets/   ← ジャーナル別テンプレート 7 種
│
│── [B] 統計・探索的解析
│   ├── scientific-eda-correlation/
│   ├── scientific-statistical-testing/
│   └── scientific-pca-tsne/
│
│── [C] 機械学習・モデリング
│   ├── scientific-ml-regression/
│   ├── scientific-ml-classification/
│   └── scientific-feature-importance/
│
│── [D] 実験計画・プロセス最適化
│   ├── scientific-doe/
│   └── scientific-process-optimization/
│
│── [E] 信号・スペクトル・時系列
│   ├── scientific-spectral-signal/
│   ├── scientific-biosignal-processing/
│   └── scientific-time-series/
│
│── [F] 生命科学・オミクス
│   ├── scientific-bioinformatics/
│   ├── scientific-metabolomics/
│   ├── scientific-sequence-analysis/
│   ├── scientific-multi-omics/
│   └── scientific-network-analysis/
│
│── [G] 化学・材料・イメージング
│   ├── scientific-cheminformatics/
│   ├── scientific-materials-characterization/
│   └── scientific-image-analysis/
│
└── [H] 臨床・疫学・メタ科学
    ├── scientific-survival-clinical/
    ├── scientific-causal-inference/
    └── scientific-meta-analysis/
```

> 注: 実際のファイルシステム上ではすべてのスキルディレクトリは `.github/skills/` 直下にフラットに配置されています。上記の中区分グルーピングは論理的な分類です。

---

## 参考

- [SATORI 使い方ガイド](../../docs/qiita-satori-guide.md)
- [GitHub Copilot Agent Skills ドキュメント](https://docs.github.com/en/copilot/concepts/agents/about-agent-skills)
- [Agent Skills オープン標準](https://github.com/agentskills/agentskills)
