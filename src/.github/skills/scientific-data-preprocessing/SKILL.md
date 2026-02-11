---
name: scientific-data-preprocessing
description: |
  科学データの前処理パイプラインスキル。欠損値補完（KNNImputer/SimpleImputer）、
  エンコーディング（LabelEncoder/OneHot/ダミー変数）、スケーリング（Standard/MinMax/Robust/Pareto）、
  対数変換、外れ値処理のテンプレートを提供。全 Exp-01〜13 に横断的に適用される基盤スキル。
---

# Scientific Data Preprocessing

科学データ解析における前処理の標準パイプライン。データクリーニング、変換、正規化の
一連の手順をトレーサブルかつ再現可能な形で実装するためのテンプレート集。
全 Exp-01〜13 に共通処理として適用される基盤スキル。

## When to Use

- CSV / DataFrame のデータを解析前にクリーニング・変換したいとき
- 欠損値の補完戦略を選択する必要があるとき
- カテゴリカル変数のエンコーディングが必要なとき
- 数値データのスケーリング・正規化が必要なとき
- 外れ値の検出・除去が必要なとき

## Quick Start

## 1. 前処理パイプライン概要

```
Raw Data
  ├─ Step 1: データ品質チェック
  ├─ Step 2: 欠損値処理
  ├─ Step 3: 外れ値処理
  ├─ Step 4: エンコーディング
  ├─ Step 5: スケーリング / 変換
  └─ Step 6: 品質確認レポート
```

## 2. データ品質チェック

```python
import pandas as pd
import numpy as np

def data_quality_report(df, name="dataset"):
    """
    データ品質レポートを出力する。
    前処理方針の決定材料として最初に実行すること。
    """
    report = pd.DataFrame({
        "dtype": df.dtypes,
        "nunique": df.nunique(),
        "missing_count": df.isnull().sum(),
        "missing_pct": (df.isnull().sum() / len(df) * 100).round(2),
        "zeros_count": (df == 0).sum(),
    })

    # 数値カラムの統計
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    for col in numeric_cols:
        report.loc[col, "mean"] = df[col].mean()
        report.loc[col, "std"] = df[col].std()
        report.loc[col, "skewness"] = df[col].skew()
        report.loc[col, "kurtosis"] = df[col].kurtosis()

    print(f"=== Data Quality Report: {name} ===")
    print(f"Shape: {df.shape}")
    print(f"Duplicated rows: {df.duplicated().sum()}")
    print(f"Columns with >50% missing: "
          f"{report[report['missing_pct'] > 50].index.tolist()}")
    print(report.to_string())

    return report
```

## 3. 欠損値処理 — 選択ガイド

| 手法 | ユースケース | 注意点 |
|---|---|---|
| `dropna()` | 完全ケース解析、欠損 <5% | サンプルサイズ減少 |
| `fillna(median)` | 数値データ、外れ値あり | 分布を歪めうる |
| `fillna(mean)` | 数値データ、正規分布 | 外れ値に敏感 |
| `SimpleImputer(strategy='most_frequent')` | カテゴリカル変数 | |
| `KNNImputer(n_neighbors=5)` | 変数間に相関がある場合 | 計算コスト高い |
| `IterativeImputer` | 多変量欠損パターン (MICE) | scikit-learn experimental |

```python
from sklearn.impute import SimpleImputer, KNNImputer

def handle_missing_values(df, strategy="auto", numeric_cols=None, categorical_cols=None):
    """
    欠損値処理パイプライン。

    strategy:
        "auto"    — 欠損率に応じて自動選択
        "drop"    — 欠損行を除去
        "median"  — 中央値で補完
        "knn"     — KNN で補完 (n_neighbors=5)
    """
    df = df.copy()

    if numeric_cols is None:
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if categorical_cols is None:
        categorical_cols = df.select_dtypes(include=["object", "category"]).columns.tolist()

    if strategy == "auto":
        missing_pct = df[numeric_cols].isnull().mean()
        # カラムごとに最適戦略を適用
        low_missing = missing_pct[missing_pct < 0.05].index.tolist()
        mid_missing = missing_pct[(missing_pct >= 0.05) & (missing_pct < 0.30)].index.tolist()
        high_missing = missing_pct[missing_pct >= 0.30].index.tolist()

        # Low: 中央値補完
        if low_missing:
            imp = SimpleImputer(strategy="median")
            df[low_missing] = imp.fit_transform(df[low_missing])

        # Mid: KNN 補完
        if mid_missing:
            imp = KNNImputer(n_neighbors=5)
            df[mid_missing] = imp.fit_transform(df[mid_missing])

        # High: 警告表示のみ (ドロップ推奨)
        if high_missing:
            print(f"WARNING: High missing rate columns ({'>30%'}): {high_missing}")
            print("Consider dropping these columns or using domain-specific imputation.")

    elif strategy == "drop":
        df = df.dropna(subset=numeric_cols)

    elif strategy == "median":
        imp = SimpleImputer(strategy="median")
        df[numeric_cols] = imp.fit_transform(df[numeric_cols])

    elif strategy == "knn":
        imp = KNNImputer(n_neighbors=5)
        df[numeric_cols] = imp.fit_transform(df[numeric_cols])

    # カテゴリカル変数 → 最頻値で補完
    if categorical_cols:
        imp_cat = SimpleImputer(strategy="most_frequent")
        df[categorical_cols] = imp_cat.fit_transform(df[categorical_cols])

    return df
```

## 4. 外れ値処理

```python
def detect_outliers(df, columns, method="iqr", threshold=1.5):
    """
    外れ値を検出する。

    method:
        "iqr"     — IQR法 (Q1 - threshold*IQR, Q3 + threshold*IQR)
        "zscore"  — Z-score 法 (|z| > threshold; default threshold=3)
        "mad"     — MAD (Median Absolute Deviation) 法
    """
    outlier_mask = pd.DataFrame(False, index=df.index, columns=columns)

    for col in columns:
        vals = df[col].dropna()
        if method == "iqr":
            Q1, Q3 = vals.quantile(0.25), vals.quantile(0.75)
            IQR = Q3 - Q1
            lower, upper = Q1 - threshold * IQR, Q3 + threshold * IQR
            outlier_mask[col] = (df[col] < lower) | (df[col] > upper)

        elif method == "zscore":
            z = np.abs((df[col] - vals.mean()) / vals.std())
            outlier_mask[col] = z > (threshold if threshold != 1.5 else 3)

        elif method == "mad":
            median = vals.median()
            mad = np.median(np.abs(vals - median))
            modified_z = 0.6745 * (df[col] - median) / (mad + 1e-10)
            outlier_mask[col] = np.abs(modified_z) > (threshold if threshold != 1.5 else 3.5)

    summary = pd.DataFrame({
        "outlier_count": outlier_mask.sum(),
        "outlier_pct": (outlier_mask.sum() / len(df) * 100).round(2),
    })

    return outlier_mask, summary


def handle_outliers(df, columns, method="iqr", action="clip"):
    """
    外れ値を処理する。

    action:
        "clip"    — 境界値にクリッピング
        "remove"  — 外れ値行を除去
        "nan"     — NaN に置換（後段で補完）
    """
    outlier_mask, summary = detect_outliers(df, columns, method=method)
    df = df.copy()

    if action == "remove":
        any_outlier = outlier_mask.any(axis=1)
        df = df[~any_outlier]

    elif action == "clip":
        for col in columns:
            vals = df[col].dropna()
            Q1, Q3 = vals.quantile(0.25), vals.quantile(0.75)
            IQR = Q3 - Q1
            lower, upper = Q1 - 1.5 * IQR, Q3 + 1.5 * IQR
            df[col] = df[col].clip(lower, upper)

    elif action == "nan":
        for col in columns:
            df.loc[outlier_mask[col], col] = np.nan

    return df, summary
```

## 5. エンコーディング — 選択ガイド

| 手法 | ユースケース | 適用例 |
|---|---|---|
| `LabelEncoder` | 順序カテゴリ、ツリー系モデル | 結晶構造タイプ → 0,1,2 |
| `pd.get_dummies()` | 名義カテゴリ、線形モデル | 材料名 → ダミー変数 |
| `OrdinalEncoder` | 順序のあるカテゴリ（明示的） | 低/中/高 → 0,1,2 |
| `TargetEncoder` | 高カーディナリティ | 化合物 ID |

```python
from sklearn.preprocessing import LabelEncoder, OrdinalEncoder

def encode_categorical(df, columns, method="auto"):
    """
    カテゴリカル変数のエンコーディング。

    method:
        "auto"      — nunique に応じて自動選択
        "label"     — LabelEncoder
        "onehot"    — pd.get_dummies()
        "ordinal"   — OrdinalEncoder
    """
    df = df.copy()
    encoders = {}

    for col in columns:
        n_unique = df[col].nunique()

        if method == "auto":
            chosen = "onehot" if n_unique <= 10 else "label"
        else:
            chosen = method

        if chosen == "label":
            le = LabelEncoder()
            df[col] = le.fit_transform(df[col].astype(str))
            encoders[col] = le

        elif chosen == "onehot":
            dummies = pd.get_dummies(df[col], prefix=col, drop_first=True)
            df = pd.concat([df.drop(columns=[col]), dummies], axis=1)

        elif chosen == "ordinal":
            oe = OrdinalEncoder()
            df[col] = oe.fit_transform(df[[col]])
            encoders[col] = oe

    return df, encoders
```

## 6. スケーリング・変換 — 選択ガイド

| 手法 | 数式 | ユースケース |
|---|---|---|
| `StandardScaler` | $(x - \mu) / \sigma$ | PCA 前、正規分布仮定 |
| `MinMaxScaler` | $(x - x_{min}) / (x_{max} - x_{min})$ | NN 入力、0〜1 に正規化 |
| `RobustScaler` | $(x - Q_2) / (Q_3 - Q_1)$ | 外れ値に頑健 |
| Pareto Scaling | $(x - \bar{x}) / \sqrt{s}$ | メタボロミクス (Exp-07) |
| Log2 Transform | $\log_2(x + 1)$ | カウント・発現量データ |
| Box-Cox | $(x^\lambda - 1)/\lambda$ | 正規性改善（正値のみ） |

```python
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler

def scale_features(df, columns, method="standard"):
    """
    特徴量スケーリング。

    method:
        "standard"  — StandardScaler (z-score)
        "minmax"    — MinMaxScaler (0-1)
        "robust"    — RobustScaler (IQR-based)
        "pareto"    — Pareto scaling (metabolomics 向き)
        "log2"      — log2(x + 1) 変換
    """
    df = df.copy()
    scaler = None

    if method == "standard":
        scaler = StandardScaler()
        df[columns] = scaler.fit_transform(df[columns])

    elif method == "minmax":
        scaler = MinMaxScaler()
        df[columns] = scaler.fit_transform(df[columns])

    elif method == "robust":
        scaler = RobustScaler()
        df[columns] = scaler.fit_transform(df[columns])

    elif method == "pareto":
        # Pareto scaling: (x - mean) / sqrt(std)
        for col in columns:
            mean = df[col].mean()
            std = df[col].std()
            df[col] = (df[col] - mean) / np.sqrt(std + 1e-10)

    elif method == "log2":
        for col in columns:
            df[col] = np.log2(df[col] + 1)

    return df, scaler
```

## 7. 前処理パイプライン統合テンプレート

```python
def preprocessing_pipeline(df, target_col=None, config=None):
    """
    全体前処理パイプライン。設定辞書で制御可能。

    config 例:
        {
            "drop_duplicates": True,
            "missing_strategy": "auto",
            "outlier_method": "iqr",
            "outlier_action": "clip",
            "encoding_method": "auto",
            "scaling_method": "standard",
        }
    """
    if config is None:
        config = {
            "drop_duplicates": True,
            "missing_strategy": "auto",
            "outlier_method": "iqr",
            "outlier_action": "clip",
            "encoding_method": "auto",
            "scaling_method": "standard",
        }

    n_original = len(df)

    # Step 0: 重複除去
    if config.get("drop_duplicates", True):
        df = df.drop_duplicates()
        print(f"  Removed {n_original - len(df)} duplicate rows")

    # Step 1: カラム分類
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    categorical_cols = df.select_dtypes(include=["object", "category"]).columns.tolist()
    if target_col and target_col in numeric_cols:
        numeric_cols.remove(target_col)
    if target_col and target_col in categorical_cols:
        categorical_cols.remove(target_col)

    # Step 2: 欠損値処理
    df = handle_missing_values(
        df,
        strategy=config.get("missing_strategy", "auto"),
        numeric_cols=numeric_cols,
        categorical_cols=categorical_cols,
    )

    # Step 3: 外れ値処理
    if numeric_cols:
        df, outlier_summary = handle_outliers(
            df, numeric_cols,
            method=config.get("outlier_method", "iqr"),
            action=config.get("outlier_action", "clip"),
        )

    # Step 4: エンコーディング
    if categorical_cols:
        df, encoders = encode_categorical(
            df, categorical_cols,
            method=config.get("encoding_method", "auto"),
        )
    else:
        encoders = {}

    # Step 5: スケーリング
    final_numeric = df.select_dtypes(include=[np.number]).columns.tolist()
    if target_col and target_col in final_numeric:
        final_numeric.remove(target_col)
    if final_numeric:
        df, scaler = scale_features(
            df, final_numeric,
            method=config.get("scaling_method", "standard"),
        )
    else:
        scaler = None

    print(f"  Preprocessing complete: {df.shape}")

    return df, {"encoders": encoders, "scaler": scaler}
```

## References

- **Exp-01**: `sc.pp.normalize_total()`, `sc.pp.log1p()` — scRNA-seq 固有の正規化
- **Exp-02**: `LabelEncoder`, `pd.get_dummies()` — 化合物記述子エンコーディング
- **Exp-03**: `StandardScaler` + `log2(x+1)` — がんデータ前処理
- **Exp-05**: `StandardScaler`, `KNNImputer` — 毒性予測データ
- **Exp-07**: Pareto scaling — メタボロミクスデータ
- **Exp-12**: `MinMaxScaler`, `RobustScaler` — プロセスデータ
- **Exp-13**: `LabelEncoder`, ダミー変数 — 材料タイプエンコーディング
