---
name: scientific-ensemble-methods
description: |
  アンサンブル学習スキル。Stacking/Blending 多段積層・
  Boosting (XGBoost/LightGBM/CatBoost) 勾配ブースティング・
  Bagging/Random Subspace・Voting 分類器/回帰器・
  アンサンブル多様性評価・モデル統合パイプライン。
tu_tools:
  - key: openml
    name: OpenML
    description: アンサンブル手法ベンチマーク参照
---

# Scientific Ensemble Methods

複数モデルの組み合わせによる予測精度向上・安定化を実現する
アンサンブル学習手法の設計・評価パイプラインを提供する。

## When to Use

- XGBoost/LightGBM/CatBoost で勾配ブースティングを実行するとき
- Stacking/Blending で多段アンサンブルを構築するとき
- 複数モデルの Voting/Averaging で安定予測を得るとき
- アンサンブルの多様性を評価するとき
- Out-of-Fold 予測でリーク防止を行うとき
- モデルの寄与度を分析するとき

---

## Quick Start

## 1. 勾配ブースティング比較

```python
import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score


def compare_boosting(X, y, cv=5, scoring="f1_macro",
                     task="classification"):
    """
    XGBoost / LightGBM / CatBoost 比較。

    Parameters:
        X: np.ndarray — 特徴量
        y: np.ndarray — ラベル
        cv: int — CV 分割数
        scoring: str — 評価指標
        task: str — "classification" / "regression"
    """
    results = []

    try:
        from xgboost import XGBClassifier, XGBRegressor
        model = (XGBClassifier(n_estimators=200, max_depth=6,
                               learning_rate=0.1, random_state=42,
                               use_label_encoder=False, eval_metric="logloss")
                 if task == "classification"
                 else XGBRegressor(n_estimators=200, max_depth=6,
                                   learning_rate=0.1, random_state=42))
        scores = cross_val_score(model, X, y, cv=cv, scoring=scoring)
        results.append({"model": "XGBoost", "mean": scores.mean(),
                        "std": scores.std()})
    except ImportError:
        pass

    try:
        from lightgbm import LGBMClassifier, LGBMRegressor
        model = (LGBMClassifier(n_estimators=200, max_depth=6,
                                learning_rate=0.1, random_state=42, verbose=-1)
                 if task == "classification"
                 else LGBMRegressor(n_estimators=200, max_depth=6,
                                    learning_rate=0.1, random_state=42, verbose=-1))
        scores = cross_val_score(model, X, y, cv=cv, scoring=scoring)
        results.append({"model": "LightGBM", "mean": scores.mean(),
                        "std": scores.std()})
    except ImportError:
        pass

    try:
        from catboost import CatBoostClassifier, CatBoostRegressor
        model = (CatBoostClassifier(iterations=200, depth=6,
                                    learning_rate=0.1, random_seed=42, verbose=0)
                 if task == "classification"
                 else CatBoostRegressor(iterations=200, depth=6,
                                        learning_rate=0.1, random_seed=42, verbose=0))
        scores = cross_val_score(model, X, y, cv=cv, scoring=scoring)
        results.append({"model": "CatBoost", "mean": scores.mean(),
                        "std": scores.std()})
    except ImportError:
        pass

    df = pd.DataFrame(results).sort_values("mean", ascending=False)
    if not df.empty:
        print(f"Boosting: best = {df.iloc[0]['model']} "
              f"({scoring} = {df.iloc[0]['mean']:.4f})")
    return df
```

## 2. Stacking アンサンブル

```python
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC


def stacking_ensemble(X_train, y_train, X_test,
                      base_models=None, meta_model=None,
                      n_folds=5):
    """
    Stacking アンサンブル (Out-of-Fold 予測)。

    Parameters:
        X_train: np.ndarray — 学習データ
        y_train: np.ndarray — 学習ラベル
        X_test: np.ndarray — テストデータ
        base_models: list | None — ベースモデル
        meta_model: classifier | None — メタモデル
        n_folds: int — CV 分割数
    """
    if base_models is None:
        base_models = [
            ("rf", RandomForestClassifier(n_estimators=200, random_state=42)),
            ("gbm", GradientBoostingClassifier(n_estimators=200, random_state=42)),
            ("svm", SVC(probability=True, random_state=42)),
        ]
    if meta_model is None:
        meta_model = LogisticRegression(max_iter=1000, random_state=42)

    kf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    n_classes = len(np.unique(y_train))

    # Out-of-Fold predictions
    oof_preds = np.zeros((len(y_train), len(base_models) * n_classes))
    test_preds = np.zeros((len(X_test), len(base_models) * n_classes))

    for i, (name, model) in enumerate(base_models):
        col_start = i * n_classes
        col_end = (i + 1) * n_classes
        test_fold_preds = np.zeros((len(X_test), n_classes, n_folds))

        for fold, (train_idx, val_idx) in enumerate(kf.split(X_train, y_train)):
            m = model.__class__(**model.get_params()).fit(
                X_train[train_idx], y_train[train_idx])
            oof_preds[val_idx, col_start:col_end] = m.predict_proba(
                X_train[val_idx])
            test_fold_preds[:, :, fold] = m.predict_proba(X_test)

        test_preds[:, col_start:col_end] = test_fold_preds.mean(axis=2)
        print(f"  Stacking base: {name} done")

    # Meta-model
    meta_model.fit(oof_preds, y_train)
    final_pred = meta_model.predict(test_preds)
    final_proba = meta_model.predict_proba(test_preds)

    print(f"Stacking: {len(base_models)} base models → meta-model")
    return final_pred, final_proba, meta_model
```

## 3. Voting アンサンブル

```python
from sklearn.ensemble import VotingClassifier, VotingRegressor


def voting_ensemble(X, y, models=None, voting="soft",
                    cv=5, scoring="f1_macro"):
    """
    Voting アンサンブル。

    Parameters:
        X: np.ndarray — 特徴量
        y: np.ndarray — ラベル
        models: list | None — (name, model) ペア
        voting: str — "soft" / "hard"
        cv: int — CV 分割数
        scoring: str — 評価指標
    """
    if models is None:
        models = [
            ("rf", RandomForestClassifier(n_estimators=200, random_state=42)),
            ("gbm", GradientBoostingClassifier(n_estimators=200, random_state=42)),
            ("lr", LogisticRegression(max_iter=1000, random_state=42)),
        ]

    # 個別モデル評価
    results = []
    for name, model in models:
        scores = cross_val_score(model, X, y, cv=cv, scoring=scoring)
        results.append({"model": name, "mean": scores.mean(), "std": scores.std()})

    # Voting
    vc = VotingClassifier(estimators=models, voting=voting)
    scores = cross_val_score(vc, X, y, cv=cv, scoring=scoring)
    results.append({"model": f"Voting({voting})",
                    "mean": scores.mean(), "std": scores.std()})

    df = pd.DataFrame(results).sort_values("mean", ascending=False)
    print(f"Voting ensemble: {scoring} = {scores.mean():.4f} ± {scores.std():.4f}")
    return df
```

## 4. アンサンブル多様性評価

```python
def ensemble_diversity(models, X, y):
    """
    アンサンブル多様性 (Q-statistic / Disagreement)。

    Parameters:
        models: list — 学習済みモデルリスト
        X: np.ndarray — 評価データ
        y: np.ndarray — 真ラベル
    """
    predictions = np.array([m.predict(X) for m in models])
    n_models = len(models)
    correct = (predictions == y).astype(int)

    # 全ペアの Q-statistic
    q_stats = []
    disagree_rates = []
    for i in range(n_models):
        for j in range(i + 1, n_models):
            n11 = np.sum((correct[i] == 1) & (correct[j] == 1))
            n00 = np.sum((correct[i] == 0) & (correct[j] == 0))
            n10 = np.sum((correct[i] == 1) & (correct[j] == 0))
            n01 = np.sum((correct[i] == 0) & (correct[j] == 1))

            denom = n11 * n00 - n10 * n01
            numer = n11 * n00 + n10 * n01
            q = denom / numer if numer != 0 else 0
            q_stats.append(q)
            disagree_rates.append((n10 + n01) / len(y))

    result = {
        "mean_q_statistic": round(np.mean(q_stats), 4),
        "mean_disagreement": round(np.mean(disagree_rates), 4),
        "n_models": n_models,
    }
    print(f"Diversity: Q={result['mean_q_statistic']:.3f}, "
          f"Disagree={result['mean_disagreement']:.3f}")
    return result
```

---

## パイプライン統合

```
automl → ensemble-methods → uncertainty-quantification
  (モデル選択)  (アンサンブル)     (不確実性定量化)
       │            │                    ↓
  feature-importance ┘         explainable-ai
    (特徴量重要度)               (説明可能 AI)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `stacking_meta.pkl` | Stacking メタモデル | → 予測 |
| `boosting_comparison.csv` | ブースティング比較 | → レポート |
| `ensemble_diversity.json` | 多様性指標 | → モデル改善 |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `openml` | OpenML | アンサンブル手法ベンチマーク参照 |
