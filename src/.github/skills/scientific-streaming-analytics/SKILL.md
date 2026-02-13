---
name: scientific-streaming-analytics
description: |
  ストリーミング解析スキル。River オンライン学習・
  リアルタイム異常検知・ストリーミング統計・
  増分データ可視化・概念ドリフト検出。
---

# Scientific Streaming Analytics

データストリームに対するリアルタイム学習・異常検知・
統計モニタリングパイプラインを提供する。

## When to Use

- データが逐次的に到着するストリーミング環境のとき
- オンライン学習で逐次モデル更新が必要なとき
- リアルタイム異常検知を実装するとき
- 概念ドリフトを検出・対応するとき
- メモリ制約下で増分的に統計量を計算するとき

---

## Quick Start

## 1. River オンライン学習

```python
from river import (
    compose, preprocessing, linear_model, metrics,
    tree, ensemble, drift)
import pandas as pd
import numpy as np


def online_learning_pipeline(stream_data, model_type="ht",
                             target_col="y",
                             feature_cols=None):
    """
    River オンライン学習パイプライン。

    Parameters:
        stream_data: iterable — (X_dict, y) のストリーム or DataFrame
        model_type: str — "ht" (Hoeffding Tree) / "lr" / "arf"
        target_col: str — 目標変数名
        feature_cols: list[str] | None — 特徴量カラム
    """
    models = {
        "ht": tree.HoeffdingTreeClassifier(),
        "lr": compose.Pipeline(
            preprocessing.StandardScaler(),
            linear_model.LogisticRegression()),
        "arf": ensemble.AdaptiveRandomForestClassifier(
            n_models=10, seed=42),
    }
    model = models.get(model_type, models["ht"])
    metric = metrics.Accuracy()
    history = []

    if isinstance(stream_data, pd.DataFrame):
        if feature_cols is None:
            feature_cols = [c for c in stream_data.columns
                           if c != target_col]
        iterator = (
            (row[feature_cols].to_dict(), row[target_col])
            for _, row in stream_data.iterrows())
    else:
        iterator = stream_data

    for i, (x, y) in enumerate(iterator):
        y_pred = model.predict_one(x)
        if y_pred is not None:
            metric.update(y, y_pred)
        model.learn_one(x, y)

        if (i + 1) % 100 == 0:
            history.append({
                "step": i + 1,
                "accuracy": metric.get(),
            })

    print(f"Online {model_type}: {metric}")
    return model, pd.DataFrame(history)
```

## 2. ストリーミング異常検知

```python
def streaming_anomaly_detection(stream_data, window_size=100,
                                threshold_sigma=3.0,
                                method="zscore"):
    """
    ストリーミング異常検知。

    Parameters:
        stream_data: iterable — 数値ストリーム
        window_size: int — スライディングウィンドウサイズ
        threshold_sigma: float — 異常判定の σ 閾値
        method: str — "zscore" / "iqr" / "ewma"
    """
    from collections import deque

    window = deque(maxlen=window_size)
    results = []
    ewma_mean = None
    ewma_var = None
    alpha = 2.0 / (window_size + 1)

    for i, value in enumerate(stream_data):
        is_anomaly = False

        if method == "zscore" and len(window) >= 10:
            mean = np.mean(window)
            std = np.std(window) + 1e-10
            z = abs(value - mean) / std
            is_anomaly = z > threshold_sigma

        elif method == "iqr" and len(window) >= 10:
            q1, q3 = np.percentile(window, [25, 75])
            iqr = q3 - q1
            lower = q1 - 1.5 * iqr
            upper = q3 + 1.5 * iqr
            is_anomaly = value < lower or value > upper

        elif method == "ewma":
            if ewma_mean is None:
                ewma_mean = value
                ewma_var = 0
            else:
                ewma_mean = alpha * value + (1 - alpha) * ewma_mean
                ewma_var = alpha * (value - ewma_mean) ** 2 + \
                    (1 - alpha) * ewma_var
                ewma_std = np.sqrt(ewma_var) + 1e-10
                is_anomaly = abs(value - ewma_mean) / ewma_std > threshold_sigma

        window.append(value)
        results.append({
            "step": i, "value": value,
            "is_anomaly": is_anomaly,
        })

    df = pd.DataFrame(results)
    n_anomalies = df["is_anomaly"].sum()
    print(f"Streaming anomaly ({method}): "
          f"{n_anomalies}/{len(df)} anomalies detected "
          f"({n_anomalies/len(df):.1%})")
    return df
```

## 3. 概念ドリフト検出

```python
def concept_drift_detection(stream_data, target_col="y",
                            feature_cols=None,
                            detector_type="adwin"):
    """
    概念ドリフト検出。

    Parameters:
        stream_data: pd.DataFrame — ストリームデータ
        target_col: str — 目標変数名
        feature_cols: list[str] | None — 特徴量カラム
        detector_type: str — "adwin" / "ddm" / "eddm"
    """
    detectors = {
        "adwin": drift.ADWIN(delta=0.002),
        "ddm": drift.DDM(min_num_instances=30),
        "eddm": drift.EDDM(),
    }
    detector = detectors.get(detector_type, detectors["adwin"])

    model = tree.HoeffdingTreeClassifier()
    metric = metrics.Accuracy()
    drift_points = []

    if feature_cols is None:
        feature_cols = [c for c in stream_data.columns
                       if c != target_col]

    for i, (_, row) in enumerate(stream_data.iterrows()):
        x = row[feature_cols].to_dict()
        y = row[target_col]
        y_pred = model.predict_one(x)

        if y_pred is not None:
            is_correct = int(y_pred == y)
            metric.update(y, y_pred)
            detector.update(is_correct)

            if detector.drift_detected:
                drift_points.append({
                    "step": i,
                    "accuracy_at_drift": metric.get(),
                })
                print(f"⚠ Drift at step {i}, acc={metric.get():.3f}")

        model.learn_one(x, y)

    print(f"Total drifts: {len(drift_points)}")
    return pd.DataFrame(drift_points)
```

---

## パイプライン統合

```
[データストリーム] → streaming-analytics → model-monitoring
                       (オンライン学習)      (性能監視)
                            │
                    anomaly-detection ← data-profiling
                      (バッチ異常検知)    (データ品質)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `online_model.pkl` | オンラインモデル | → 推論 |
| `stream_anomalies.csv` | 異常検知結果 | → alerting |
| `drift_report.csv` | ドリフト検出点 | → model-monitoring |
