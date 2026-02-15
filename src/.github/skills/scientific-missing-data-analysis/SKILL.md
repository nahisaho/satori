---
name: scientific-missing-data-analysis
description: |
  欠損データ解析スキル。欠損パターン診断 (MCAR/MAR/MNAR) ・
  Little's MCAR テスト・多重代入法 (MICE) ・KNN 補完・
  MissForest・VAE/GAIN 補完・欠損パターン可視化・Rubin's Rules。
tu_tools:
  - key: biotools
    name: bio.tools
    description: 欠損データ処理ツール検索
---

# Scientific Missing Data Analysis

欠損データの診断・補完・感度分析パイプラインを提供し、
バイアスのない統計推論を実現する。

## When to Use

- データセットの欠損パターンを診断するとき
- MCAR / MAR / MNAR のメカニズムを判定するとき
- 多重代入法 (MICE) で欠損値を補完するとき
- KNN / MissForest / 深層学習ベースの補完をするとき
- 複数の補完結果を Rubin's Rules で統合するとき
- 欠損パターンを可視化するとき

---

## Quick Start

## 1. 欠損パターン診断

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def diagnose_missing_patterns(df, output_prefix="missing"):
    """
    欠損パターン診断 — MCAR/MAR/MNAR 判定支援。

    Parameters:
        df: pd.DataFrame — 入力データ
        output_prefix: str — 出力ファイル接頭辞
    """
    n_rows, n_cols = df.shape
    missing_counts = df.isnull().sum()
    missing_pct = (missing_counts / n_rows * 100).round(2)

    summary = pd.DataFrame({
        "column": df.columns,
        "n_missing": missing_counts.values,
        "pct_missing": missing_pct.values,
        "dtype": df.dtypes.values
    }).sort_values("pct_missing", ascending=False)

    # 欠損パターン行列 (msno 風)
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # (1) 欠損マトリックス
    ax = axes[0, 0]
    missing_matrix = df.isnull().astype(int)
    ax.imshow(missing_matrix.values[:200], aspect="auto", cmap="Greys",
              interpolation="none")
    ax.set_xlabel("Features")
    ax.set_ylabel("Samples")
    ax.set_title("Missing Pattern Matrix (first 200 rows)")

    # (2) 欠損率バー
    ax = axes[0, 1]
    cols_with_missing = summary[summary["pct_missing"] > 0]
    ax.barh(cols_with_missing["column"], cols_with_missing["pct_missing"])
    ax.set_xlabel("Missing %")
    ax.set_title("Missing Rate per Column")

    # (3) 欠損相関ヒートマップ
    ax = axes[1, 0]
    miss_corr = df.isnull().corr()
    sns.heatmap(miss_corr, ax=ax, cmap="RdBu_r", center=0,
                square=True, cbar_kws={"shrink": 0.8})
    ax.set_title("Missing Correlation")

    # (4) 欠損パターン上位
    ax = axes[1, 1]
    patterns = df.isnull().apply(lambda x: tuple(x), axis=1)
    pattern_counts = patterns.value_counts().head(10)
    ax.barh(range(len(pattern_counts)),
            pattern_counts.values)
    ax.set_yticks(range(len(pattern_counts)))
    ax.set_yticklabels([str(p)[:40] for p in pattern_counts.index],
                       fontsize=7)
    ax.set_xlabel("Count")
    ax.set_title("Top 10 Missing Patterns")

    plt.tight_layout()
    path = f"{output_prefix}_diagnosis.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Missing Diagnosis: {n_cols} cols, "
          f"{missing_counts.sum()} total missing ({(missing_counts.sum()/(n_rows*n_cols)*100):.1f}%)")
    return {"summary": summary, "fig": path}
```

## 2. Little's MCAR テスト

```python
def littles_mcar_test(df):
    """
    Little's MCAR テスト — 完全ランダム欠損の検定。

    Parameters:
        df: pd.DataFrame — 数値データのみ
    Returns:
        dict — chi2 統計量, p値, 判定
    """
    from scipy import stats

    numeric_df = df.select_dtypes(include=[np.number])
    n_rows, n_cols = numeric_df.shape

    # 欠損パターンごとにグルーピング
    patterns = numeric_df.isnull().apply(tuple, axis=1)
    unique_patterns = patterns.unique()

    # 全体平均と全体共分散
    global_mean = numeric_df.mean()
    global_cov = numeric_df.cov()

    chi2_stat = 0.0
    df_stat = 0

    for pattern in unique_patterns:
        mask = patterns == pattern
        sub_df = numeric_df[mask]
        n_j = len(sub_df)
        if n_j < 2:
            continue

        # このパターンで観測されているカラム
        obs_cols = [i for i, m in enumerate(pattern) if not m]
        if len(obs_cols) == 0:
            continue

        obs_mean = sub_df.iloc[:, obs_cols].mean().values
        exp_mean = global_mean.iloc[obs_cols].values
        diff = obs_mean - exp_mean

        obs_cov = global_cov.iloc[obs_cols, obs_cols].values
        try:
            cov_inv = np.linalg.pinv(obs_cov / n_j)
        except np.linalg.LinAlgError:
            continue

        chi2_stat += diff @ cov_inv @ diff
        df_stat += len(obs_cols)

    df_stat -= n_cols  # 自由度補正

    if df_stat <= 0:
        return {"chi2": np.nan, "p_value": np.nan,
                "conclusion": "判定不能 (自由度不足)"}

    p_value = 1 - stats.chi2.cdf(chi2_stat, df_stat)
    conclusion = "MCAR (p > 0.05)" if p_value > 0.05 else "Not MCAR (p ≤ 0.05)"

    print(f"Little's MCAR test: χ²={chi2_stat:.2f}, df={df_stat}, "
          f"p={p_value:.4f} → {conclusion}")
    return {"chi2": chi2_stat, "df": df_stat,
            "p_value": p_value, "conclusion": conclusion}
```

## 3. 多重代入法 (MICE)

```python
def mice_imputation(df, n_imputations=5, max_iter=10, random_state=42):
    """
    MICE (Multiple Imputation by Chained Equations)。

    Parameters:
        df: pd.DataFrame — 欠損を含むデータ
        n_imputations: int — 代入データセット数
        max_iter: int — 反復回数
        random_state: int — 乱数シード
    """
    from sklearn.experimental import enable_iterative_imputer  # noqa
    from sklearn.impute import IterativeImputer

    numeric_cols = df.select_dtypes(include=[np.number]).columns
    cat_cols = df.select_dtypes(exclude=[np.number]).columns

    imputed_datasets = []

    for i in range(n_imputations):
        imputer = IterativeImputer(
            max_iter=max_iter,
            random_state=random_state + i,
            sample_posterior=True)

        imputed_numeric = pd.DataFrame(
            imputer.fit_transform(df[numeric_cols]),
            columns=numeric_cols, index=df.index)

        imputed_df = imputed_numeric.copy()
        for col in cat_cols:
            imputed_df[col] = df[col].fillna(df[col].mode().iloc[0]
                                              if not df[col].mode().empty else "UNKNOWN")

        imputed_datasets.append(imputed_df)

    print(f"MICE: {n_imputations} datasets × {max_iter} iterations, "
          f"{len(numeric_cols)} numeric cols")
    return imputed_datasets


def rubins_rules(estimates, variances):
    """
    Rubin's Rules — 多重代入結果の統合。

    Parameters:
        estimates: list[float] — 各代入データセットからの推定値
        variances: list[float] — 各代入データセットからの分散
    """
    m = len(estimates)
    Q_bar = np.mean(estimates)
    U_bar = np.mean(variances)  # Within-imputation variance
    B = np.var(estimates, ddof=1)  # Between-imputation variance
    T = U_bar + (1 + 1 / m) * B  # Total variance

    # 自由度 (Barnard-Rubin)
    r = (1 + 1 / m) * B / U_bar if U_bar > 0 else np.inf
    df_old = (m - 1) * (1 + 1 / r) ** 2 if r > 0 else np.inf

    print(f"Rubin's Rules: Q̄={Q_bar:.4f}, T={T:.4f}, "
          f"within={U_bar:.4f}, between={B:.4f}")
    return {"pooled_estimate": Q_bar, "total_variance": T,
            "within_variance": U_bar, "between_variance": B,
            "df": df_old}
```

## 4. KNN / MissForest 補完

```python
def knn_imputation(df, n_neighbors=5):
    """
    KNN 欠損値補完。

    Parameters:
        df: pd.DataFrame — 欠損を含むデータ
        n_neighbors: int — 近傍数
    """
    from sklearn.impute import KNNImputer

    numeric_cols = df.select_dtypes(include=[np.number]).columns
    imputer = KNNImputer(n_neighbors=n_neighbors)
    imputed = pd.DataFrame(
        imputer.fit_transform(df[numeric_cols]),
        columns=numeric_cols, index=df.index)

    n_imputed = df[numeric_cols].isnull().sum().sum()
    print(f"KNN Imputation (k={n_neighbors}): {n_imputed} values imputed")
    return imputed


def missforest_imputation(df, n_estimators=100, max_iter=10):
    """
    MissForest (Random Forest ベースの反復補完)。

    Parameters:
        df: pd.DataFrame — 欠損を含むデータ
        n_estimators: int — Random Forest の木の数
        max_iter: int — 反復回数
    """
    from sklearn.experimental import enable_iterative_imputer  # noqa
    from sklearn.impute import IterativeImputer
    from sklearn.ensemble import RandomForestRegressor

    numeric_cols = df.select_dtypes(include=[np.number]).columns

    imputer = IterativeImputer(
        estimator=RandomForestRegressor(n_estimators=n_estimators,
                                        random_state=42, n_jobs=-1),
        max_iter=max_iter, random_state=42)

    imputed = pd.DataFrame(
        imputer.fit_transform(df[numeric_cols]),
        columns=numeric_cols, index=df.index)

    n_imputed = df[numeric_cols].isnull().sum().sum()
    print(f"MissForest (n_trees={n_estimators}, iter={max_iter}): "
          f"{n_imputed} values imputed")
    return imputed
```

---

## パイプライン統合

```
eda-correlation → missing-data-analysis → ml-classification
  (探索的解析)      (欠損診断・補完)       (モデリング)
       │                  │                     ↓
 statistical-testing ────┘          advanced-visualization
   (統計検定)                         (結果可視化)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `missing_diagnosis.png` | 欠損パターン可視化 | → reporting |
| `mcar_test_result.json` | Little's MCAR テスト | → 補完戦略選択 |
| `imputed_datasets/` | MICE 多重代入データ | → ml-classification |
| `imputation_comparison.csv` | 補完手法比較 | → 最終選択 |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `biotools` | bio.tools | 欠損データ処理ツール検索 |
