---
name: scientific-anomaly-detection
description: |
  異常検知・外れ値検出スキル。Isolation Forest・LOF・
  One-Class SVM・Autoencoder 異常検知・統計的工程管理 (SPC)・
  多変量異常検知・異常スコアリング・閾値最適化。
---

# Scientific Anomaly Detection

科学データにおける異常値・外れ値・異常パターンの検出と
統計的工程管理 (SPC) パイプラインを提供する。

## When to Use

- 実験データの外れ値を統計的に検出するとき
- 製造プロセスの異常監視 (SPC) をするとき
- 多変量データで異常パターンを発見するとき
- Autoencoder で複雑な異常を検出するとき
- 異常スコアの閾値を最適化するとき
- 複数手法のアンサンブル異常検知をするとき

---

## Quick Start

## 1. 統計的異常検知アンサンブル

```python
import numpy as np
import pandas as pd
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.svm import OneClassSVM
from sklearn.preprocessing import StandardScaler


def anomaly_detection_ensemble(X, contamination=0.05,
                               methods=None, threshold_vote=2):
    """
    複数手法アンサンブル異常検知。

    Parameters:
        X: np.ndarray | pd.DataFrame — 入力データ
        contamination: float — 想定異常率
        methods: list[str] | None — 使用手法 ("iforest", "lof", "ocsvm")
        threshold_vote: int — 最低投票数 (多数決)
    """
    if methods is None:
        methods = ["iforest", "lof", "ocsvm"]

    if isinstance(X, pd.DataFrame):
        feature_names = X.columns.tolist()
        X_arr = X.values
    else:
        feature_names = [f"f{i}" for i in range(X.shape[1])]
        X_arr = X

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_arr)

    results = {}
    predictions = {}

    for method in methods:
        if method == "iforest":
            model = IsolationForest(
                contamination=contamination, random_state=42, n_jobs=-1)
            preds = model.fit_predict(X_scaled)
            scores = -model.score_samples(X_scaled)
        elif method == "lof":
            model = LocalOutlierFactor(
                n_neighbors=20, contamination=contamination)
            preds = model.fit_predict(X_scaled)
            scores = -model.negative_outlier_factor_
        elif method == "ocsvm":
            model = OneClassSVM(kernel="rbf", nu=contamination)
            preds = model.fit_predict(X_scaled)
            scores = -model.decision_function(X_scaled)
        else:
            continue

        is_anomaly = (preds == -1).astype(int)
        predictions[method] = is_anomaly
        results[method] = {
            "n_anomalies": int(is_anomaly.sum()),
            "scores": scores
        }

    # アンサンブル多数決
    vote_matrix = np.column_stack(list(predictions.values()))
    ensemble_votes = vote_matrix.sum(axis=1)
    ensemble_anomaly = (ensemble_votes >= threshold_vote).astype(int)

    result_df = pd.DataFrame(X_arr, columns=feature_names)
    for method, preds in predictions.items():
        result_df[f"anomaly_{method}"] = preds
    result_df["ensemble_votes"] = ensemble_votes
    result_df["is_anomaly"] = ensemble_anomaly

    n_ens = ensemble_anomaly.sum()
    print(f"Anomaly Ensemble ({len(methods)} methods, vote≥{threshold_vote}): "
          f"{n_ens}/{len(X_arr)} anomalies ({n_ens/len(X_arr)*100:.1f}%)")

    for m, r in results.items():
        print(f"  {m}: {r['n_anomalies']} anomalies")

    return result_df, results
```

## 2. Autoencoder 異常検知

```python
def autoencoder_anomaly(X, encoding_dim=8, epochs=100,
                        threshold_percentile=95):
    """
    Autoencoder ベース異常検知。

    Parameters:
        X: np.ndarray — 入力データ (正常データで学習)
        encoding_dim: int — 潜在次元数
        epochs: int — 学習エポック数
        threshold_percentile: float — 再構成誤差の閾値パーセンタイル
    """
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader, TensorDataset

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    n_features = X_scaled.shape[1]

    # Autoencoder 定義
    class AE(nn.Module):
        def __init__(self):
            super().__init__()
            self.encoder = nn.Sequential(
                nn.Linear(n_features, 64), nn.ReLU(),
                nn.Linear(64, 32), nn.ReLU(),
                nn.Linear(32, encoding_dim))
            self.decoder = nn.Sequential(
                nn.Linear(encoding_dim, 32), nn.ReLU(),
                nn.Linear(32, 64), nn.ReLU(),
                nn.Linear(64, n_features))

        def forward(self, x):
            z = self.encoder(x)
            return self.decoder(z)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = AE().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    criterion = nn.MSELoss()

    X_tensor = torch.FloatTensor(X_scaled).to(device)
    dataset = TensorDataset(X_tensor, X_tensor)
    loader = DataLoader(dataset, batch_size=64, shuffle=True)

    model.train()
    for epoch in range(epochs):
        total_loss = 0
        for batch_x, _ in loader:
            optimizer.zero_grad()
            recon = model(batch_x)
            loss = criterion(recon, batch_x)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

    # 再構成誤差
    model.eval()
    with torch.no_grad():
        recon = model(X_tensor).cpu().numpy()

    recon_errors = np.mean((X_scaled - recon) ** 2, axis=1)
    threshold = np.percentile(recon_errors, threshold_percentile)
    is_anomaly = (recon_errors > threshold).astype(int)

    print(f"Autoencoder Anomaly: threshold={threshold:.4f} (P{threshold_percentile}), "
          f"{is_anomaly.sum()} anomalies")
    return {"reconstruction_error": recon_errors, "threshold": threshold,
            "is_anomaly": is_anomaly, "model": model}
```

## 3. 統計的工程管理 (SPC)

```python
def spc_control_chart(data, column, subgroup_size=1,
                      chart_type="individuals"):
    """
    SPC 管理図 (X-bar, R, Individuals-MR)。

    Parameters:
        data: pd.DataFrame | pd.Series — 時系列データ
        column: str — 対象カラム名
        subgroup_size: int — サブグループサイズ
        chart_type: str — "individuals" / "xbar_r" / "cusum"
    """
    import matplotlib.pyplot as plt

    if isinstance(data, pd.DataFrame):
        values = data[column].values
    else:
        values = data.values

    if chart_type == "individuals":
        x_bar = np.mean(values)
        mr = np.abs(np.diff(values))
        mr_bar = np.mean(mr)
        d2 = 1.128  # d2 for n=2

        ucl = x_bar + 3 * (mr_bar / d2)
        lcl = x_bar - 3 * (mr_bar / d2)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

        # Individuals chart
        ax1.plot(values, "b-o", markersize=3)
        ax1.axhline(x_bar, color="g", linestyle="-", label=f"CL={x_bar:.3f}")
        ax1.axhline(ucl, color="r", linestyle="--", label=f"UCL={ucl:.3f}")
        ax1.axhline(lcl, color="r", linestyle="--", label=f"LCL={lcl:.3f}")

        # OOC points
        ooc = np.where((values > ucl) | (values < lcl))[0]
        if len(ooc) > 0:
            ax1.scatter(ooc, values[ooc], c="red", s=50, zorder=5,
                        label=f"OOC ({len(ooc)})")
        ax1.set_title("Individuals Chart")
        ax1.legend(fontsize=8)

        # Moving Range chart
        mr_ucl = 3.267 * mr_bar
        ax2.plot(mr, "b-o", markersize=3)
        ax2.axhline(mr_bar, color="g", linestyle="-")
        ax2.axhline(mr_ucl, color="r", linestyle="--")
        ax2.set_title("Moving Range Chart")

        plt.tight_layout()
        path = "spc_control_chart.png"
        plt.savefig(path, dpi=150, bbox_inches="tight")
        plt.close()

        print(f"SPC Individuals: CL={x_bar:.3f}, UCL={ucl:.3f}, "
              f"LCL={lcl:.3f}, OOC={len(ooc)}")
        return {"cl": x_bar, "ucl": ucl, "lcl": lcl,
                "ooc_indices": ooc, "fig": path}

    elif chart_type == "cusum":
        target = np.mean(values)
        se = np.std(values)
        k = 0.5 * se
        h = 5 * se

        cusum_pos = np.zeros(len(values))
        cusum_neg = np.zeros(len(values))

        for i in range(1, len(values)):
            cusum_pos[i] = max(0, cusum_pos[i-1] + (values[i] - target) - k)
            cusum_neg[i] = min(0, cusum_neg[i-1] + (values[i] - target) + k)

        fig, ax = plt.subplots(figsize=(12, 5))
        ax.plot(cusum_pos, "b-", label="CUSUM+")
        ax.plot(cusum_neg, "r-", label="CUSUM-")
        ax.axhline(h, color="b", linestyle="--", alpha=0.5)
        ax.axhline(-h, color="r", linestyle="--", alpha=0.5)
        ax.set_title("CUSUM Control Chart")
        ax.legend()

        path = "cusum_chart.png"
        plt.savefig(path, dpi=150, bbox_inches="tight")
        plt.close()

        print(f"CUSUM: target={target:.3f}, k={k:.3f}, h={h:.3f}")
        return {"target": target, "k": k, "h": h,
                "cusum_pos": cusum_pos, "cusum_neg": cusum_neg, "fig": path}
```

---

## パイプライン統合

```
eda-correlation → anomaly-detection → ml-classification
  (探索的解析)      (外れ値検出)        (モデリング)
       │                 │                  ↓
 data-profiling ────────┘         model-monitoring
  (データ品質)                     (モデル監視)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `anomaly_ensemble.csv` | アンサンブル異常検知結果 | → EDA |
| `autoencoder_anomaly.json` | AE 異常スコア | → reporting |
| `spc_control_chart.png` | SPC 管理図 | → process-optimization |
