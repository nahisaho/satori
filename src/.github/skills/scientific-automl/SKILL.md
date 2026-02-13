---
name: scientific-automl
description: |
  AutoML パイプラインスキル。Optuna ハイパーパラメータ最適化・
  FLAML 高速 AutoML・Auto-sklearn モデル選択・
  NAS (Neural Architecture Search)・
  特徴量エンジニアリング自動化・モデル比較パイプライン。
---

# Scientific AutoML

ハイパーパラメータ最適化・モデル選択・特徴量エンジニアリングを
自動化する AutoML パイプラインを提供する。

## When to Use

- Optuna/Hyperopt でハイパーパラメータを最適化するとき
- 複数モデルの自動比較・選択を行うとき
- FLAML/Auto-sklearn で高速な AutoML を実行するとき
- 特徴量エンジニアリングを自動化するとき
- Neural Architecture Search (NAS) を設計するとき
- モデル選択根拠のレポートを生成するとき

---

## Quick Start

## 1. Optuna ハイパーパラメータ最適化

```python
import optuna
import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import (
    RandomForestClassifier, GradientBoostingClassifier)
from sklearn.svm import SVC
from sklearn.metrics import make_scorer, f1_score


def optuna_optimize(X, y, model_type="rf", n_trials=100,
                    cv=5, scoring="f1_macro", direction="maximize"):
    """
    Optuna ベース ハイパーパラメータ最適化。

    Parameters:
        X: np.ndarray — 特徴量
        y: np.ndarray — ラベル
        model_type: str — "rf" / "gbm" / "svm"
        n_trials: int — 試行回数
        cv: int — CV 分割数
        scoring: str — 評価指標
        direction: str — "maximize" / "minimize"
    """
    def objective(trial):
        if model_type == "rf":
            params = {
                "n_estimators": trial.suggest_int("n_estimators", 50, 500),
                "max_depth": trial.suggest_int("max_depth", 3, 20),
                "min_samples_split": trial.suggest_int("min_samples_split", 2, 20),
                "min_samples_leaf": trial.suggest_int("min_samples_leaf", 1, 10),
                "max_features": trial.suggest_categorical(
                    "max_features", ["sqrt", "log2", None]),
            }
            model = RandomForestClassifier(**params, random_state=42)

        elif model_type == "gbm":
            params = {
                "n_estimators": trial.suggest_int("n_estimators", 50, 500),
                "max_depth": trial.suggest_int("max_depth", 3, 10),
                "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.3, log=True),
                "subsample": trial.suggest_float("subsample", 0.5, 1.0),
                "min_samples_split": trial.suggest_int("min_samples_split", 2, 20),
            }
            model = GradientBoostingClassifier(**params, random_state=42)

        elif model_type == "svm":
            params = {
                "C": trial.suggest_float("C", 0.01, 100, log=True),
                "kernel": trial.suggest_categorical(
                    "kernel", ["rbf", "poly", "sigmoid"]),
                "gamma": trial.suggest_categorical("gamma", ["scale", "auto"]),
            }
            model = SVC(**params, probability=True, random_state=42)

        scores = cross_val_score(model, X, y, cv=cv, scoring=scoring)
        return scores.mean()

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(direction=direction)
    study.optimize(objective, n_trials=n_trials, show_progress_bar=True)

    print(f"Optuna ({model_type}): best {scoring} = {study.best_value:.4f}")
    print(f"  Best params: {study.best_params}")
    return study
```

## 2. マルチモデル AutoML パイプライン

```python
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier


def automl_model_selection(X, y, cv=5, scoring="f1_macro",
                           n_trials_per_model=50):
    """
    AutoML マルチモデル選択パイプライン。

    Parameters:
        X: np.ndarray — 特徴量
        y: np.ndarray — ラベル
        cv: int — CV 分割数
        scoring: str — 評価指標
        n_trials_per_model: int — モデルあたり試行数
    """
    model_types = ["rf", "gbm", "svm"]
    results = []

    for mt in model_types:
        study = optuna_optimize(
            X, y, model_type=mt,
            n_trials=n_trials_per_model, cv=cv, scoring=scoring)
        results.append({
            "model_type": mt,
            "best_score": round(study.best_value, 4),
            "best_params": study.best_params,
            "n_trials": len(study.trials),
        })

    # 簡易モデル (ベースライン)
    baselines = [
        ("logistic", LogisticRegression(max_iter=1000, random_state=42)),
        ("knn", KNeighborsClassifier()),
        ("dt", DecisionTreeClassifier(random_state=42)),
    ]
    for name, model in baselines:
        scores = cross_val_score(model, X, y, cv=cv, scoring=scoring)
        results.append({
            "model_type": name,
            "best_score": round(scores.mean(), 4),
            "best_params": {},
            "n_trials": 1,
        })

    df = pd.DataFrame(results).sort_values("best_score", ascending=False)
    best = df.iloc[0]
    print(f"AutoML: best = {best['model_type']} "
          f"({scoring} = {best['best_score']})")
    return df
```

## 3. 自動特徴量エンジニアリング

```python
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.feature_selection import SelectKBest, mutual_info_classif


def auto_feature_engineering(X, y, max_poly_degree=2,
                             top_k=None, interactions_only=False):
    """
    自動特徴量エンジニアリング。

    Parameters:
        X: np.ndarray — 元特徴量
        y: np.ndarray — ラベル
        max_poly_degree: int — 多項式次数
        top_k: int | None — 選択する特徴量数
        interactions_only: bool — 交互作用のみ
    """
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # 多項式特徴量
    poly = PolynomialFeatures(
        degree=max_poly_degree,
        interaction_only=interactions_only,
        include_bias=False)
    X_poly = poly.fit_transform(X_scaled)

    # 特徴量選択
    if top_k is None:
        top_k = min(X_poly.shape[1], X.shape[1] * 3)

    selector = SelectKBest(mutual_info_classif, k=min(top_k, X_poly.shape[1]))
    X_selected = selector.fit_transform(X_poly, y)

    print(f"Feature engineering: {X.shape[1]} → {X_poly.shape[1]} "
          f"→ {X_selected.shape[1]} features")
    return X_selected, poly, selector
```

## 4. Optuna 可視化レポート

```python
def automl_report(study, output_dir="results"):
    """
    Optuna Study 可視化レポート。

    Parameters:
        study: optuna.Study — 最適化結果
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    import matplotlib.pyplot as plt

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # パラメータ重要度
    try:
        importances = optuna.importance.get_param_importances(study)
        fig, ax = plt.subplots(figsize=(8, 5))
        params = list(importances.keys())
        values = list(importances.values())
        ax.barh(params, values)
        ax.set_xlabel("Importance")
        ax.set_title("Hyperparameter Importance")
        fig.tight_layout()
        fig.savefig(out / "param_importance.png", dpi=150)
        plt.close(fig)
    except Exception:
        pass

    # 最適化履歴
    trials_df = study.trials_dataframe()
    trials_df.to_csv(out / "optuna_trials.csv", index=False)

    # ベストパラメータ
    best = {
        "best_value": study.best_value,
        "best_params": study.best_params,
        "n_trials": len(study.trials),
    }

    print(f"AutoML report → {out}")
    return best
```

---

## パイプライン統合

```
eda-correlation → automl → ensemble-methods
  (データ探索)     (モデル選択)  (アンサンブル)
       │               │              ↓
  feature-importance ──┘     uncertainty-quantification
    (特徴量解釈)              (不確実性定量化)
       │
  active-learning
    (能動学習)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `optuna_trials.csv` | 試行履歴 | → 可視化 |
| `param_importance.png` | パラメータ重要度 | → レポート |
| `model_comparison.csv` | モデル比較 | → ensemble-methods |
