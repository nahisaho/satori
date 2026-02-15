---
name: scientific-eda-correlation
description: |
  探索的データ解析（EDA）と相関分析のスキル。データの分布可視化、相関ヒートマップ、
  散布図行列の作成を行う際に使用。Scientific Skills Exp-02, 12, 13 で確立したパターン。
tu_tools:
  - key: biotools
    name: bio.tools
    description: 統計解析ツール検索
---

# Scientific EDA & Correlation Analysis

探索的データ解析（Exploratory Data Analysis）のパイプラインスキル。
データ理解の初期段階で使用し、分布・外れ値・変数間相関を把握する。

## When to Use

- 新しいデータセットを受け取ったとき
- 変数間の関係性を把握したいとき
- 相関ヒートマップを作成したいとき
- 材料別・群別のボックスプロット比較が必要なとき

## Quick Start

## 標準パイプライン

### 1. 記述統計量の算出

```python
import pandas as pd
import numpy as np

def descriptive_statistics(df, numeric_cols, group_col=None):
    """記述統計量を算出して CSV に保存する。"""
    if group_col:
        stats = df.groupby(group_col)[numeric_cols].describe()
    else:
        stats = df[numeric_cols].describe()
    stats.to_csv("results/descriptive_statistics.csv")
    return stats
```

### 2. 分布可視化（ボックスプロット + バイオリンプロット）

```python
import matplotlib.pyplot as plt
import seaborn as sns

def plot_distributions(df, variables, group_col, figsize=(20, 16), ncols=3):
    """群別のボックスプロットを変数ごとに描画する。"""
    nrows = (len(variables) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = axes.flatten()

    for i, var in enumerate(variables):
        sns.boxplot(data=df, x=group_col, y=var, ax=axes[i],
                    palette="Set2", showfliers=True)
        axes[i].set_title(var, fontsize=12, fontweight="bold")
        axes[i].tick_params(axis="x", rotation=45)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    plt.savefig("figures/distribution_boxplots.png", dpi=300, bbox_inches="tight")
    plt.close()
```

### 3. 相関ヒートマップ（Exp-02 / Exp-13 パターン）

```python
def plot_correlation_heatmap(df, numeric_cols, block_boundaries=None,
                              figsize=(14, 12), method="pearson"):
    """
    相関ヒートマップを描画する。
    block_boundaries: PSP などの階層境界を示す線の位置リスト（オプション）。
    """
    corr = df[numeric_cols].corr(method=method)

    fig, ax = plt.subplots(figsize=figsize)
    mask = np.triu(np.ones_like(corr, dtype=bool), k=1)
    sns.heatmap(corr, mask=mask, annot=True, fmt=".2f",
                cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                square=True, linewidths=0.5, ax=ax,
                annot_kws={"size": 8})

    # 階層境界線（PSP ブロック分離）
    if block_boundaries:
        for b in block_boundaries:
            ax.axhline(y=b, color="black", linewidth=2)
            ax.axvline(x=b, color="black", linewidth=2)

    ax.set_title("Correlation Heatmap", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig("figures/correlation_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close()
    return corr
```

### 4. 散布図行列

```python
def plot_scatter_matrix(df, variables, hue_col, figsize=(16, 14)):
    """主要変数の散布図行列を描画する。"""
    g = sns.pairplot(df[variables + [hue_col]], hue=hue_col,
                     diag_kind="kde", palette="Set2",
                     plot_kws={"alpha": 0.6, "s": 30})
    g.fig.suptitle("Scatter Matrix", y=1.02, fontsize=14, fontweight="bold")
    plt.savefig("figures/scatter_matrix.png", dpi=300, bbox_inches="tight")
    plt.close()
```

### 5. PSP ブロック相関分析（Exp-13 独自）

```python
def psp_block_correlation(df, process_cols, structure_cols, property_cols):
    """Process→Structure→Property の 3 ブロック相関を個別に算出する。"""
    ps_corr = df[process_cols + structure_cols].corr().loc[process_cols, structure_cols]
    sp_corr = df[structure_cols + property_cols].corr().loc[structure_cols, property_cols]
    pp_corr = df[process_cols + property_cols].corr().loc[process_cols, property_cols]

    ps_corr.to_csv("results/PSP_process_structure_corr.csv")
    sp_corr.to_csv("results/PSP_structure_property_corr.csv")
    pp_corr.to_csv("results/PSP_process_property_corr.csv")

    return ps_corr, sp_corr, pp_corr
```

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `biotools` | bio.tools | 統計解析ツール検索 |

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/descriptive_statistics.csv` | CSV |
| `figures/distribution_boxplots.png` | PNG (300 DPI) |
| `figures/correlation_heatmap.png` | PNG (300 DPI) |
| `figures/scatter_matrix.png` | PNG (300 DPI) |

#### 参照実験

- **Exp-02**: `sns.heatmap` 相関ヒートマップの基本パターン
- **Exp-12**: 8 プロセスパラメータの EDA
- **Exp-13**: PSP 3 ブロック相関行列
