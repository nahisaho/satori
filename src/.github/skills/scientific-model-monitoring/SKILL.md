---
name: scientific-model-monitoring
description: |
  MLOps モデル監視スキル。データドリフト検出 (Evidently/NannyML)・
  モデル性能劣化検出・特徴量ドリフト・コンセプトドリフト・
  A/B テスト統計・モデルレジストリ・再学習トリガー。
tu_tools:
  - key: openml
    name: OpenML
    description: モデルモニタリング指標・ベンチマーク
---

# Scientific Model Monitoring

本番環境の ML モデル監視パイプラインを提供し、
データドリフト・性能劣化を検出して再学習トリガーを実現する。

## When to Use

- デプロイ済みモデルの予測品質を継続監視するとき
- データドリフト (共変量シフト) を検出するとき
- コンセプトドリフト (P(Y|X) の変化) を検出するとき
- A/B テストで新旧モデルを比較するとき
- 特徴量分布の変化を追跡するとき
- 再学習トリガーの自動化ルールを設定するとき

---

## Quick Start

## 1. データドリフト検出

```python
import numpy as np
import pandas as pd
from scipy import stats


def detect_data_drift(reference_df, current_df,
                      method="ks", threshold=0.05):
    """
    データドリフト検出 — 参照データ vs 現在データ。

    Parameters:
        reference_df: pd.DataFrame — 学習時データ (参照)
        current_df: pd.DataFrame — 推論時データ (現在)
        method: str — "ks" (KS 検定) / "psi" (PSI) / "wasserstein"
        threshold: float — 有意水準 or PSI 閾値
    """
    numeric_cols = reference_df.select_dtypes(include=[np.number]).columns
    common_cols = [c for c in numeric_cols if c in current_df.columns]

    drift_results = []

    for col in common_cols:
        ref_vals = reference_df[col].dropna().values
        cur_vals = current_df[col].dropna().values

        if method == "ks":
            stat, p_value = stats.ks_2samp(ref_vals, cur_vals)
            is_drift = p_value < threshold
            drift_results.append({
                "feature": col, "statistic": stat,
                "p_value": p_value, "is_drift": is_drift})

        elif method == "psi":
            # Population Stability Index
            psi_val = _compute_psi(ref_vals, cur_vals)
            is_drift = psi_val > 0.2  # >0.2 = significant shift
            drift_results.append({
                "feature": col, "psi": psi_val,
                "is_drift": is_drift,
                "severity": "high" if psi_val > 0.25 else
                            "medium" if psi_val > 0.1 else "low"})

        elif method == "wasserstein":
            w_dist = stats.wasserstein_distance(ref_vals, cur_vals)
            ref_std = np.std(ref_vals)
            normalized = w_dist / ref_std if ref_std > 0 else w_dist
            is_drift = normalized > 0.1
            drift_results.append({
                "feature": col, "wasserstein": w_dist,
                "normalized": normalized, "is_drift": is_drift})

    result_df = pd.DataFrame(drift_results)
    n_drift = result_df["is_drift"].sum()
    print(f"Data Drift ({method}): {n_drift}/{len(common_cols)} features drifted")
    return result_df


def _compute_psi(expected, actual, n_bins=10):
    """PSI (Population Stability Index) 計算。"""
    breakpoints = np.quantile(expected, np.linspace(0, 1, n_bins + 1))
    breakpoints[0] = -np.inf
    breakpoints[-1] = np.inf

    expected_pct = np.histogram(expected, bins=breakpoints)[0] / len(expected)
    actual_pct = np.histogram(actual, bins=breakpoints)[0] / len(actual)

    expected_pct = np.clip(expected_pct, 1e-4, None)
    actual_pct = np.clip(actual_pct, 1e-4, None)

    psi = np.sum((actual_pct - expected_pct) * np.log(actual_pct / expected_pct))
    return psi
```

## 2. モデル性能劣化検出

```python
def detect_performance_degradation(y_true_batches, y_pred_batches,
                                    metric="accuracy",
                                    window_size=10, alert_threshold=0.05):
    """
    モデル性能劣化のスライディングウィンドウ検出。

    Parameters:
        y_true_batches: list[np.ndarray] — バッチごとの真値
        y_pred_batches: list[np.ndarray] — バッチごとの予測値
        metric: str — "accuracy" / "f1" / "rmse" / "auc"
        window_size: int — 移動平均ウィンドウ
        alert_threshold: float — 性能低下アラート閾値
    """
    from sklearn.metrics import accuracy_score, f1_score, mean_squared_error
    from sklearn.metrics import roc_auc_score
    import matplotlib.pyplot as plt

    metric_funcs = {
        "accuracy": accuracy_score,
        "f1": lambda y, p: f1_score(y, p, average="macro"),
        "rmse": lambda y, p: -np.sqrt(mean_squared_error(y, p)),
        "auc": lambda y, p: roc_auc_score(y, p)
    }

    func = metric_funcs[metric]
    scores = [func(yt, yp) for yt, yp in zip(y_true_batches, y_pred_batches)]

    # 移動平均
    scores_arr = np.array(scores)
    if len(scores_arr) >= window_size:
        ma = np.convolve(scores_arr, np.ones(window_size)/window_size, mode="valid")
    else:
        ma = scores_arr

    # ベースライン (最初の window_size バッチ)
    baseline = np.mean(scores_arr[:window_size])
    current = np.mean(scores_arr[-window_size:])
    degradation = baseline - current

    is_degraded = degradation > alert_threshold

    # 可視化
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(scores, "b-o", markersize=3, alpha=0.5, label="Batch score")
    if len(ma) > 0:
        ax.plot(range(window_size - 1, window_size - 1 + len(ma)),
                ma, "r-", linewidth=2, label=f"MA({window_size})")
    ax.axhline(baseline, color="g", linestyle="--",
               label=f"Baseline={baseline:.4f}")
    ax.axhline(baseline - alert_threshold, color="orange", linestyle="--",
               label=f"Alert={baseline - alert_threshold:.4f}")
    ax.set_xlabel("Batch")
    ax.set_ylabel(metric)
    ax.set_title(f"Model Performance Monitoring ({metric})")
    ax.legend()

    path = "performance_monitoring.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

    status = "DEGRADED ⚠️" if is_degraded else "OK ✓"
    print(f"Performance ({metric}): baseline={baseline:.4f}, "
          f"current={current:.4f}, Δ={degradation:.4f} → {status}")
    return {"baseline": baseline, "current": current,
            "degradation": degradation, "is_degraded": is_degraded,
            "scores": scores, "fig": path}
```

## 3. A/B テスト統計

```python
def ab_test_models(y_true, preds_a, preds_b, metric="accuracy",
                   n_bootstrap=10000, alpha=0.05):
    """
    A/B テスト — 2 モデルの統計的比較。

    Parameters:
        y_true: np.ndarray — 真値
        preds_a: np.ndarray — モデル A 予測
        preds_b: np.ndarray — モデル B 予測
        metric: str — 評価指標
        n_bootstrap: int — ブートストラップ回数
        alpha: float — 有意水準
    """
    from sklearn.metrics import accuracy_score, f1_score, mean_squared_error

    metric_funcs = {
        "accuracy": accuracy_score,
        "f1": lambda y, p: f1_score(y, p, average="macro"),
        "rmse": lambda y, p: np.sqrt(mean_squared_error(y, p))
    }

    func = metric_funcs[metric]
    score_a = func(y_true, preds_a)
    score_b = func(y_true, preds_b)

    # Bootstrap confidence interval for difference
    diffs = []
    n = len(y_true)
    rng = np.random.RandomState(42)

    for _ in range(n_bootstrap):
        idx = rng.choice(n, n, replace=True)
        sa = func(y_true[idx], preds_a[idx])
        sb = func(y_true[idx], preds_b[idx])
        diffs.append(sb - sa)

    diffs = np.array(diffs)
    ci_lower = np.percentile(diffs, 100 * alpha / 2)
    ci_upper = np.percentile(diffs, 100 * (1 - alpha / 2))
    p_value = np.mean(diffs <= 0)  # P(B ≤ A)

    winner = "B" if ci_lower > 0 else ("A" if ci_upper < 0 else "Tie")

    print(f"A/B Test ({metric}): A={score_a:.4f}, B={score_b:.4f}")
    print(f"  Δ(B-A)={score_b - score_a:.4f}, "
          f"95% CI=[{ci_lower:.4f}, {ci_upper:.4f}], "
          f"p={p_value:.4f} → Winner: {winner}")
    return {"score_a": score_a, "score_b": score_b,
            "diff": score_b - score_a, "ci": (ci_lower, ci_upper),
            "p_value": p_value, "winner": winner}
```

---

## パイプライン統合

```
ensemble-methods → model-monitoring → anomaly-detection
  (モデル構築)      (監視)             (異常検知)
       │                │                  ↓
  automl ──────────────┘          active-learning
   (AutoML)                        (再学習)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `drift_report.csv` | ドリフト検出結果 | → 再学習判断 |
| `performance_monitoring.png` | 性能推移 | → reporting |
| `ab_test_result.json` | A/B テスト結果 | → デプロイ判断 |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `openml` | OpenML | モデルモニタリング指標・ベンチマーク |
