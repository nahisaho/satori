---
name: scientific-semi-supervised-learning
description: |
  半教師あり学習スキル。Self-Training・Label Propagation・
  MixMatch/FixMatch・Pseudo-Labeling・ラベル効率評価。
---

# Scientific Semi-Supervised Learning

少量のラベル付きデータと大量の未ラベルデータを活用する
半教師あり学習パイプラインを提供する。

## When to Use

- ラベル付きデータが少量しかないとき
- アノテーションコストが高く全量ラベリングが困難なとき
- Self-Training で反復的にラベルを拡張するとき
- グラフベースの Label Propagation を適用するとき
- Pseudo-Labeling の信頼度閾値を設計するとき

---

## Quick Start

## 1. Self-Training パイプライン

```python
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import accuracy_score, classification_report


def self_training_pipeline(X_labeled, y_labeled, X_unlabeled,
                           base_estimator=None, threshold=0.95,
                           max_iterations=10, batch_size=None,
                           X_test=None, y_test=None):
    """
    Self-Training 半教師あり学習。

    Parameters:
        X_labeled: np.ndarray — ラベル付き特徴量
        y_labeled: np.ndarray — ラベル
        X_unlabeled: np.ndarray — 未ラベル特徴量
        base_estimator: sklearn estimator | None — 基底分類器
        threshold: float — Pseudo-Label 採用閾値
        max_iterations: int — 最大反復回数
        batch_size: int | None — 各反復で追加するサンプル数上限
        X_test: np.ndarray | None — テスト特徴量
        y_test: np.ndarray | None — テストラベル
    """
    from sklearn.ensemble import GradientBoostingClassifier

    if base_estimator is None:
        base_estimator = GradientBoostingClassifier(
            n_estimators=100, random_state=42)

    X_train = X_labeled.copy()
    y_train = y_labeled.copy()
    X_pool = X_unlabeled.copy()
    history = []

    for iteration in range(max_iterations):
        if len(X_pool) == 0:
            print(f"Iteration {iteration}: Pool exhausted")
            break

        model = clone(base_estimator)
        model.fit(X_train, y_train)
        proba = model.predict_proba(X_pool)
        max_proba = proba.max(axis=1)
        pseudo_labels = proba.argmax(axis=1)

        confident_mask = max_proba >= threshold
        n_confident = confident_mask.sum()

        if batch_size and n_confident > batch_size:
            top_idx = np.argsort(max_proba)[-batch_size:]
            confident_mask = np.zeros(len(X_pool), dtype=bool)
            confident_mask[top_idx] = True
            n_confident = batch_size

        if n_confident == 0:
            print(f"Iteration {iteration}: No confident samples")
            break

        X_train = np.vstack([X_train, X_pool[confident_mask]])
        y_train = np.concatenate([
            y_train, pseudo_labels[confident_mask]])
        X_pool = X_pool[~confident_mask]

        record = {"iteration": iteration,
                  "n_labeled": len(X_train),
                  "n_pool": len(X_pool),
                  "n_added": int(n_confident),
                  "mean_confidence": float(max_proba[confident_mask].mean())}

        if X_test is not None and y_test is not None:
            test_acc = accuracy_score(y_test, model.predict(X_test))
            record["test_accuracy"] = test_acc

        history.append(record)
        print(f"Iter {iteration}: +{n_confident} samples, "
              f"total={len(X_train)}, pool={len(X_pool)}")

    final_model = clone(base_estimator)
    final_model.fit(X_train, y_train)
    return final_model, pd.DataFrame(history)
```

## 2. Label Propagation

```python
def label_propagation_ssl(X_all, y_partial, kernel="rbf",
                          gamma=20, n_neighbors=7,
                          max_iter=1000):
    """
    グラフベース Label Propagation。

    Parameters:
        X_all: np.ndarray — 全サンプル特徴量 (ラベル付き+未ラベル)
        y_partial: np.ndarray — ラベル (-1 = 未ラベル)
        kernel: str — "rbf" / "knn"
        gamma: float — RBF カーネルの γ
        n_neighbors: int — KNN カーネルの k
        max_iter: int — 最大反復回数
    """
    from sklearn.semi_supervised import (
        LabelPropagation, LabelSpreading)

    models = {
        "propagation": LabelPropagation(
            kernel=kernel, gamma=gamma,
            n_neighbors=n_neighbors, max_iter=max_iter),
        "spreading": LabelSpreading(
            kernel=kernel, gamma=gamma,
            n_neighbors=n_neighbors, max_iter=max_iter, alpha=0.2),
    }

    results = {}
    for name, model in models.items():
        model.fit(X_all, y_partial)
        y_pred = model.transduction_
        n_propagated = (y_partial == -1).sum()
        results[name] = {
            "model": model,
            "predictions": y_pred,
            "n_propagated": int(n_propagated),
            "label_distributions": model.label_distributions_,
        }
        print(f"{name}: propagated {n_propagated} labels")

    return results
```

## 3. Pseudo-Labeling 品質評価

```python
def evaluate_pseudo_labels(y_true_unlabeled, pseudo_labels,
                           confidences, thresholds=None):
    """
    Pseudo-Label の品質を評価。

    Parameters:
        y_true_unlabeled: np.ndarray — 真のラベル (評価用)
        pseudo_labels: np.ndarray — 予測した疑似ラベル
        confidences: np.ndarray — 各予測の信頼度
        thresholds: list[float] | None — 閾値リスト
    """
    if thresholds is None:
        thresholds = [0.5, 0.7, 0.8, 0.9, 0.95, 0.99]

    records = []
    for t in thresholds:
        mask = confidences >= t
        if mask.sum() == 0:
            continue
        acc = accuracy_score(y_true_unlabeled[mask],
                             pseudo_labels[mask])
        records.append({
            "threshold": t,
            "n_selected": int(mask.sum()),
            "coverage": float(mask.mean()),
            "pseudo_accuracy": acc,
        })
        print(f"τ={t:.2f}: {mask.sum()} samples, "
              f"coverage={mask.mean():.1%}, acc={acc:.3f}")

    return pd.DataFrame(records)
```

---

## パイプライン統合

```
[少量ラベル] → semi-supervised-learning → ml-classification
                   (ラベル拡張)              (本分類)
                        │
                 active-learning ← data-profiling
                  (能動学習)       (データ品質)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `self_training_history.csv` | 反復学習履歴 | → 収束分析 |
| `pseudo_label_quality.csv` | 疑似ラベル品質 | → 閾値選択 |
| `propagated_labels.npy` | 伝播ラベル | → ml-classification |
