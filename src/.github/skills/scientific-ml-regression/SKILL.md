---
name: scientific-ml-regression
description: |
  マルチターゲット回帰モデルの学習・評価・比較スキル。複数の回帰モデル（Ridge, Lasso,
  Random Forest, Gradient Boosting, Extra Trees）を KFold 交差検証で比較する際に使用。
  Scientific Skills Exp-12, 13 で確立したパターン。
---

# Scientific ML Regression Pipeline

複数の回帰モデルを統一フレームワークで学習・評価・比較するパイプライン。
マルチターゲット回帰（複数の出力変数を同時に予測）に対応。

## When to Use

- 連続値の予測タスク
- 複数モデルの性能比較が必要なとき
- マルチターゲット回帰（複数出力変数の同時予測）
- 交差検証による汎化性能の評価

## Quick Start

## 標準パイプライン

### 1. モデル定義

```python
from sklearn.linear_model import Ridge, Lasso, ElasticNet
from sklearn.ensemble import (RandomForestRegressor,
                              GradientBoostingRegressor,
                              ExtraTreesRegressor)
from sklearn.model_selection import KFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import numpy as np
import pandas as pd

MODEL_DEFS = {
    "Ridge": Ridge(alpha=1.0),
    "Lasso": Lasso(alpha=0.1, max_iter=10000),
    "Random Forest": RandomForestRegressor(n_estimators=200, max_depth=15,
                                           random_state=42),
    "Gradient Boosting": GradientBoostingRegressor(n_estimators=200,
                                                    random_state=42),
    "Extra Trees": ExtraTreesRegressor(n_estimators=200, max_depth=15,
                                       random_state=42),
}
```

### 2. マルチターゲット学習 & 評価

```python
def train_evaluate_models(X_train, X_test, y_train, y_test,
                          target_names, model_defs=None, n_splits=5):
    """
    全モデル × 全ターゲットの学習・評価を一括実行する。
    Returns: results_df (モデル×ターゲットの評価指標), best_models dict
    """
    if model_defs is None:
        model_defs = MODEL_DEFS

    scaler = StandardScaler()
    X_train_sc = scaler.fit_transform(X_train)
    X_test_sc = scaler.transform(X_test)

    results = []
    best_models = {}
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    for target in target_names:
        y_tr = y_train[target].values
        y_te = y_test[target].values
        best_r2 = -np.inf

        for name, model in model_defs.items():
            import copy
            m = copy.deepcopy(model)
            m.fit(X_train_sc, y_tr)
            y_pred = m.predict(X_test_sc)

            r2 = r2_score(y_te, y_pred)
            rmse = np.sqrt(mean_squared_error(y_te, y_pred))
            mae = mean_absolute_error(y_te, y_pred)

            # 交差検証
            cv_scores = cross_val_score(
                copy.deepcopy(model), X_train_sc, y_tr,
                cv=kf, scoring="neg_root_mean_squared_error"
            )
            cv_rmse_mean = -cv_scores.mean()
            cv_rmse_std = cv_scores.std()

            results.append({
                "Target": target,
                "Model": name,
                "R2": r2,
                "RMSE": rmse,
                "MAE": mae,
                "CV_RMSE_mean": cv_rmse_mean,
                "CV_RMSE_std": cv_rmse_std,
            })

            if r2 > best_r2:
                best_r2 = r2
                best_models[target] = {"model": m, "name": name, "r2": r2}

    results_df = pd.DataFrame(results)
    results_df.to_csv("results/model_metrics.csv", index=False)
    return results_df, best_models
```

### 3. モデル比較可視化（R² バーチャート）

```python
import matplotlib.pyplot as plt
import seaborn as sns

def plot_model_comparison_r2(results_df, figsize=(14, 6)):
    """R² の grouped bar chart を描画する。"""
    pivot = results_df.pivot(index="Target", columns="Model", values="R2")
    ax = pivot.plot(kind="bar", figsize=figsize, colormap="Set2", edgecolor="black")
    ax.set_ylabel("R² Score")
    ax.set_title("Model Comparison: R² by Target", fontweight="bold")
    ax.set_ylim(0, 1.05)
    ax.legend(title="Model", bbox_to_anchor=(1.05, 1))
    ax.axhline(y=0.8, color="gray", linestyle="--", alpha=0.5, label="R²=0.8")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig("figures/model_comparison_r2.png", dpi=300, bbox_inches="tight")
    plt.close()
```

### 4. Actual vs Predicted プロット

```python
def plot_actual_vs_predicted(X_test_sc, y_test, best_models, target_names,
                             ncols=3, figsize=(18, 12)):
    """最良モデルの実測値 vs 予測値プロットを描画する。"""
    nrows = (len(target_names) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = axes.flatten()

    for i, target in enumerate(target_names):
        m = best_models[target]["model"]
        y_pred = m.predict(X_test_sc)
        y_true = y_test[target].values
        r2 = best_models[target]["r2"]

        axes[i].scatter(y_true, y_pred, alpha=0.5, s=20, edgecolors="k", linewidth=0.3)
        lims = [min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())]
        axes[i].plot(lims, lims, "r--", linewidth=1.5)
        axes[i].set_xlabel("Actual")
        axes[i].set_ylabel("Predicted")
        axes[i].set_title(f"{target}\n({best_models[target]['name']}, R²={r2:.3f})",
                         fontsize=10, fontweight="bold")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    plt.savefig("figures/actual_vs_predicted.png", dpi=300, bbox_inches="tight")
    plt.close()
```

### 5. レーダーチャート比較（Exp-05 パターン）

```python
def plot_radar_comparison(results_df, target_name, figsize=(8, 8)):
    """指定ターゲットのモデル比較レーダーチャートを描画する。"""
    subset = results_df[results_df["Target"] == target_name]
    metrics = ["R2", "1-RMSE_norm", "1-MAE_norm", "Stability"]

    # 正規化（全指標を 0-1, 高い方が良い に統一）
    subset = subset.copy()
    for col in ["RMSE", "MAE"]:
        max_val = subset[col].max()
        subset[f"1-{col}_norm"] = 1 - subset[col] / max_val if max_val > 0 else 1
    subset["Stability"] = 1 - subset["CV_RMSE_std"] / subset["CV_RMSE_mean"]

    N = len(metrics)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))

    for _, row in subset.iterrows():
        values = [row[m] for m in metrics] + [row[metrics[0]]]
        ax.plot(angles, values, "o-", linewidth=2, label=row["Model"])
        ax.fill(angles, values, alpha=0.1)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metrics)
    ax.set_title(f"Radar: {target_name}", fontweight="bold")
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1))
    plt.tight_layout()
    plt.savefig(f"figures/radar_{target_name}.png", dpi=300, bbox_inches="tight")
    plt.close()
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/model_metrics.csv` | CSV |
| `figures/model_comparison_r2.png` | PNG |
| `figures/actual_vs_predicted.png` | PNG |
| `figures/radar_*.png` | PNG |

#### 参照実験

- **Exp-12**: 6 モデル × 4 ターゲット（エッチング品質特性）
- **Exp-13**: 5 モデル × 6 ターゲット（PSP 構造・物性）
- **Exp-05**: レーダーチャートによるモデル比較
