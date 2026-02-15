---
name: scientific-active-learning
description: |
  アクティブラーニング (能動学習) スキル。不確実性サンプリング・
  Query-by-Committee・期待モデル変化・プール型/ストリーム型・
  バッチアクティブラーニング・停止基準判定・
  モデル改善パイプライン。
tu_tools:
  - key: openml
    name: OpenML
    description: 能動学習データセット・評価指標
---

# Scientific Active Learning

ラベル付けコストを最小化しながらモデル精度を最大化する
アクティブラーニング戦略の設計・実行・評価パイプラインを提供する。

## When to Use

- 少量のラベル付きデータで効率的にモデルを改善するとき
- ラベル付けコストが高い実験データを扱うとき
- 不確実性サンプリングでモデルの弱点を特定するとき
- Query-by-Committee で意見が分かれるサンプルを選択するとき
- バッチアクティブラーニングで複数サンプルを同時取得するとき
- 停止基準 (パフォーマンス収束) を判定するとき

---

## Quick Start

## 1. 不確実性サンプリング

```python
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score


def uncertainty_sampling(model, X_pool, strategy="entropy"):
    """
    不確実性ベース サンプリング。

    Parameters:
        model: fitted sklearn classifier
        X_pool: np.ndarray — ラベルなしプール
        strategy: str — "entropy" / "margin" / "least_confident"
    Returns:
        indices: np.ndarray — 不確実性降順のインデックス
    """
    proba = model.predict_proba(X_pool)

    if strategy == "entropy":
        scores = -np.sum(proba * np.log(proba + 1e-10), axis=1)
    elif strategy == "margin":
        sorted_p = np.sort(proba, axis=1)
        scores = 1.0 - (sorted_p[:, -1] - sorted_p[:, -2])
    elif strategy == "least_confident":
        scores = 1.0 - np.max(proba, axis=1)
    else:
        raise ValueError(f"Unknown strategy: {strategy}")

    indices = np.argsort(scores)[::-1]
    print(f"Uncertainty ({strategy}): top score = {scores[indices[0]]:.4f}")
    return indices


def query_by_committee(models, X_pool, n_members=5):
    """
    Query-by-Committee サンプリング。

    Parameters:
        models: list — 学習済み分類器リスト
        X_pool: np.ndarray — ラベルなしプール
        n_members: int — 委員会メンバ数
    """
    predictions = np.array([m.predict(X_pool) for m in models[:n_members]])
    # Vote entropy
    n_samples = X_pool.shape[0]
    n_classes = len(np.unique(predictions))
    scores = np.zeros(n_samples)

    for i in range(n_samples):
        votes = predictions[:, i]
        _, counts = np.unique(votes, return_counts=True)
        proba = counts / len(votes)
        scores[i] = -np.sum(proba * np.log(proba + 1e-10))

    indices = np.argsort(scores)[::-1]
    print(f"QBC: top disagreement = {scores[indices[0]]:.4f}")
    return indices
```

## 2. バッチアクティブラーニング

```python
def batch_active_learning(model, X_pool, batch_size=10,
                          strategy="entropy", diversity_weight=0.5):
    """
    多様性を考慮したバッチアクティブラーニング。

    Parameters:
        model: fitted classifier
        X_pool: np.ndarray — ラベルなしプール
        batch_size: int — バッチサイズ
        strategy: str — 不確実性戦略
        diversity_weight: float — 多様性重み (0-1)
    """
    from sklearn.metrics.pairwise import euclidean_distances

    # 不確実性スコア
    indices = uncertainty_sampling(model, X_pool, strategy)

    # 候補プール (上位 batch_size * 3)
    candidate_size = min(batch_size * 3, len(indices))
    candidates = indices[:candidate_size]

    # 多様性ベース選択 (k-center greedy)
    selected = [candidates[0]]
    for _ in range(batch_size - 1):
        remaining = [c for c in candidates if c not in selected]
        if not remaining:
            break

        dists = euclidean_distances(
            X_pool[remaining], X_pool[selected])
        min_dists = dists.min(axis=1)

        # 不確実性ランクを正規化
        uncertainty_ranks = np.array([
            np.where(indices == r)[0][0] for r in remaining])
        uncertainty_scores = 1.0 - uncertainty_ranks / len(indices)

        # 複合スコア
        combined = (diversity_weight * min_dists / (min_dists.max() + 1e-10)
                    + (1 - diversity_weight) * uncertainty_scores)

        best_idx = remaining[np.argmax(combined)]
        selected.append(best_idx)

    print(f"Batch AL: selected {len(selected)} diverse-uncertain samples")
    return np.array(selected)
```

## 3. アクティブラーニングループ

```python
def active_learning_loop(X_labeled, y_labeled, X_pool, y_pool_true,
                         X_test, y_test,
                         model=None, n_rounds=20, batch_size=10,
                         strategy="entropy"):
    """
    アクティブラーニング実験ループ。

    Parameters:
        X_labeled: np.ndarray — 初期ラベル付きデータ
        y_labeled: np.ndarray — 初期ラベル
        X_pool: np.ndarray — ラベルなしプール
        y_pool_true: np.ndarray — プールの真ラベル (Oracle)
        X_test: np.ndarray — テストデータ
        y_test: np.ndarray — テストラベル
        model: sklearn classifier (default: RF)
        n_rounds: int — ラウンド数
        batch_size: int — バッチサイズ
        strategy: str — サンプリング戦略
    """
    if model is None:
        model = RandomForestClassifier(n_estimators=100, random_state=42)

    X_l = X_labeled.copy()
    y_l = y_labeled.copy()
    X_p = X_pool.copy()
    y_p = y_pool_true.copy()

    history = []
    for rnd in range(n_rounds):
        m = clone(model).fit(X_l, y_l)
        acc = accuracy_score(y_test, m.predict(X_test))
        history.append({
            "round": rnd,
            "n_labeled": len(y_l),
            "accuracy": round(acc, 4),
            "pool_size": len(y_p),
        })

        if len(X_p) == 0:
            break

        # バッチ選択
        selected = batch_active_learning(
            m, X_p, batch_size, strategy)

        # Oracle に問い合わせ
        X_l = np.vstack([X_l, X_p[selected]])
        y_l = np.concatenate([y_l, y_p[selected]])

        mask = np.ones(len(X_p), dtype=bool)
        mask[selected] = False
        X_p = X_p[mask]
        y_p = y_p[mask]

    df = pd.DataFrame(history)
    improvement = df["accuracy"].iloc[-1] - df["accuracy"].iloc[0]
    print(f"AL loop: {n_rounds} rounds, "
          f"acc {df['accuracy'].iloc[0]:.3f} → {df['accuracy'].iloc[-1]:.3f} "
          f"(+{improvement:.3f})")
    return df


def compare_strategies(X_labeled, y_labeled, X_pool, y_pool_true,
                       X_test, y_test, n_rounds=20, batch_size=10):
    """
    複数アクティブラーニング戦略の比較。

    Parameters: (同上)
    """
    strategies = ["entropy", "margin", "least_confident"]
    results = {}

    for strat in strategies:
        history = active_learning_loop(
            X_labeled, y_labeled, X_pool, y_pool_true,
            X_test, y_test,
            n_rounds=n_rounds, batch_size=batch_size,
            strategy=strat)
        results[strat] = history

    # ランダムベースライン
    np.random.seed(42)
    random_history = active_learning_loop(
        X_labeled, y_labeled,
        X_pool[np.random.permutation(len(X_pool))],
        y_pool_true[np.random.permutation(len(y_pool_true))],
        X_test, y_test,
        n_rounds=n_rounds, batch_size=batch_size,
        strategy="least_confident")
    results["random"] = random_history

    print(f"Strategy comparison: {len(strategies) + 1} methods evaluated")
    return results
```

## 4. 停止基準判定

```python
def stopping_criterion(history_df, patience=5, min_improvement=0.001):
    """
    アクティブラーニング停止基準判定。

    Parameters:
        history_df: pd.DataFrame — AL 履歴
        patience: int — 改善なしラウンド数
        min_improvement: float — 最小改善幅
    """
    accs = history_df["accuracy"].values
    if len(accs) < patience + 1:
        return False, "insufficient rounds"

    recent = accs[-patience:]
    best_before = accs[:-patience].max()
    improvement = recent.max() - best_before

    if improvement < min_improvement:
        return True, (f"converged: improvement {improvement:.5f} "
                      f"< threshold {min_improvement}")

    return False, f"continuing: improvement {improvement:.5f}"
```

---

## パイプライン統合

```
eda-correlation → active-learning → ml-classification
  (データ探索)     (サンプル選択)     (モデル構築)
       │                │                ↓
  missing-data ─────────┘     ensemble-methods
    (欠損値処理)               (アンサンブル)
                                    ↓
                          uncertainty-quantification
                            (不確実性定量化)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `al_history.csv` | AL ラウンド履歴 | → 停止判定 |
| `selected_samples.csv` | 選択サンプル | → ラベル付け |
| `strategy_comparison.csv` | 戦略比較 | → advanced-visualization |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `openml` | OpenML | 能動学習データセット・評価指標 |
