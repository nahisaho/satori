---
name: scientific-uncertainty-quantification
description: |
  不確実性定量化スキル。Conformal Prediction・MC Dropout・
  深層アンサンブル・アレアトリック / エピステミック分離・
  Calibration Curve・予測区間推定・Expected Calibration Error。
tu_tools:
  - key: papers_with_code
    name: Papers with Code
    description: 不確実性定量化手法・ベンチマーク
---

# Scientific Uncertainty Quantification

機械学習予測の不確実性を定量化し、信頼できる予測区間と
校正された確率推定を提供する。

## When to Use

- モデル予測の信頼度を定量化したいとき
- 予測区間・信頼区間を推定するとき
- 校正曲線でモデルの確率推定品質を評価するとき
- Conformal Prediction で分布フリーの保証区間を得るとき
- MC Dropout でベイズ近似推論するとき
- アレアトリック / エピステミック不確実性を分離するとき

---

## Quick Start

## 1. Conformal Prediction

```python
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator


def conformal_prediction(model, X_calib, y_calib, X_test,
                         alpha=0.1, method="quantile"):
    """
    Conformal Prediction — 分布フリーの予測区間推定。

    Parameters:
        model: BaseEstimator — 学習済みモデル
        X_calib: np.ndarray — 校正データ特徴量
        y_calib: np.ndarray — 校正データ目的変数
        X_test: np.ndarray — テストデータ特徴量
        alpha: float — 有意水準 (1-α がカバレッジ保証)
        method: str — "quantile" / "cqr" (Conformalized Quantile Regression)
    """
    # 分類 or 回帰を自動判定
    if hasattr(model, "predict_proba"):
        return _conformal_classification(model, X_calib, y_calib, X_test, alpha)
    else:
        return _conformal_regression(model, X_calib, y_calib, X_test, alpha, method)


def _conformal_regression(model, X_calib, y_calib, X_test, alpha, method):
    preds_calib = model.predict(X_calib)
    residuals = np.abs(y_calib - preds_calib)

    # (1-α)(1 + 1/n) 分位点
    n = len(residuals)
    q = np.quantile(residuals, np.ceil((1 - alpha) * (n + 1)) / n)

    preds_test = model.predict(X_test)
    lower = preds_test - q
    upper = preds_test + q

    coverage = None
    width = np.mean(upper - lower)

    result_df = pd.DataFrame({
        "prediction": preds_test,
        "lower": lower,
        "upper": upper,
        "interval_width": upper - lower
    })

    print(f"Conformal Prediction (α={alpha}): "
          f"mean width = {width:.4f}")
    return result_df


def _conformal_classification(model, X_calib, y_calib, X_test, alpha):
    proba_calib = model.predict_proba(X_calib)
    n_classes = proba_calib.shape[1]

    # Non-conformity score = 1 - 正解クラスの確率
    scores = 1 - proba_calib[np.arange(len(y_calib)), y_calib]
    n = len(scores)
    q = np.quantile(scores, np.ceil((1 - alpha) * (n + 1)) / n)

    proba_test = model.predict_proba(X_test)
    prediction_sets = []
    for i in range(len(X_test)):
        pred_set = np.where(proba_test[i] >= 1 - q)[0].tolist()
        prediction_sets.append(pred_set)

    avg_size = np.mean([len(s) for s in prediction_sets])
    print(f"Conformal Classification (α={alpha}): "
          f"avg set size = {avg_size:.2f}")
    return prediction_sets
```

## 2. MC Dropout

```python
def mc_dropout_predict(model, X, n_forward=100, dropout_rate=0.1):
    """
    MC Dropout — ベイズ近似推論。

    Parameters:
        model: nn.Module — ドロップアウト付きモデル
        X: torch.Tensor — 入力テンソル
        n_forward: int — フォワードパス回数
        dropout_rate: float — ドロップアウト率
    """
    import torch

    model.train()  # Dropout を有効にする

    predictions = []
    with torch.no_grad():
        for _ in range(n_forward):
            pred = model(X)
            if pred.dim() > 1 and pred.shape[1] > 1:
                pred = torch.softmax(pred, dim=1)
            predictions.append(pred.cpu().numpy())

    predictions = np.array(predictions)  # (n_forward, n_samples, ...)
    mean_pred = predictions.mean(axis=0)
    std_pred = predictions.std(axis=0)

    # エピステミック不確実性: 予測のばらつき
    epistemic = std_pred
    # 分類の場合: 予測エントロピー
    if mean_pred.ndim == 2 and mean_pred.shape[1] > 1:
        entropy = -np.sum(mean_pred * np.log(mean_pred + 1e-10), axis=1)
    else:
        entropy = None

    print(f"MC Dropout: {n_forward} passes, "
          f"mean epistemic unc = {epistemic.mean():.4f}")
    return {"mean": mean_pred, "std": std_pred,
            "epistemic": epistemic, "entropy": entropy}
```

## 3. 深層アンサンブル不確実性

```python
def deep_ensemble_uncertainty(models, X, task="classification"):
    """
    深層アンサンブル不確実性推定 (Lakshminarayanan+ 2017)。

    Parameters:
        models: list[nn.Module] — アンサンブルメンバー
        X: torch.Tensor — 入力テンソル
        task: str — "classification" / "regression"
    """
    import torch

    all_preds = []
    for m in models:
        m.eval()
        with torch.no_grad():
            pred = m(X)
            if task == "classification":
                pred = torch.softmax(pred, dim=1)
            all_preds.append(pred.cpu().numpy())

    all_preds = np.array(all_preds)  # (M, N, ...)

    # 平均予測
    mean_pred = all_preds.mean(axis=0)

    if task == "classification":
        # アレアトリック: 平均エントロピー
        aleatoric = np.mean([
            -np.sum(p * np.log(p + 1e-10), axis=1) for p in all_preds
        ], axis=0)
        # 全体エントロピー
        total_entropy = -np.sum(mean_pred * np.log(mean_pred + 1e-10), axis=1)
        # エピステミック = 全体 - アレアトリック (相互情報量)
        epistemic = total_entropy - aleatoric
    else:
        # 回帰: 分散分解
        aleatoric = np.zeros(len(X))  # 各モデルの分散 (ここでは簡略)
        epistemic = all_preds.var(axis=0).squeeze()
        total_entropy = None

    result = {
        "mean_prediction": mean_pred,
        "aleatoric": aleatoric,
        "epistemic": epistemic,
        "total_uncertainty": (aleatoric + epistemic) if task == "regression" else total_entropy
    }

    print(f"Deep Ensemble ({len(models)} models): "
          f"epistemic={epistemic.mean():.4f}, "
          f"aleatoric={aleatoric.mean():.4f}")
    return result
```

## 4. Calibration 評価

```python
def calibration_analysis(y_true, y_proba, n_bins=10):
    """
    校正曲線とECE (Expected Calibration Error) を計算。

    Parameters:
        y_true: np.ndarray — 真のラベル (0/1)
        y_proba: np.ndarray — 予測確率
        n_bins: int — ビン数
    """
    import matplotlib.pyplot as plt
    from sklearn.calibration import calibration_curve

    # 校正曲線
    frac_pos, mean_pred = calibration_curve(
        y_true, y_proba, n_bins=n_bins, strategy="uniform")

    # ECE
    bin_edges = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    bin_stats = []
    for i in range(n_bins):
        mask = (y_proba >= bin_edges[i]) & (y_proba < bin_edges[i + 1])
        if mask.sum() == 0:
            continue
        bin_acc = y_true[mask].mean()
        bin_conf = y_proba[mask].mean()
        bin_count = mask.sum()
        ece += (bin_count / len(y_true)) * abs(bin_acc - bin_conf)
        bin_stats.append({
            "bin": f"[{bin_edges[i]:.1f}, {bin_edges[i+1]:.1f})",
            "count": int(bin_count),
            "accuracy": bin_acc,
            "confidence": bin_conf,
            "gap": abs(bin_acc - bin_conf)
        })

    # 可視化
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # 校正曲線
    ax1.plot([0, 1], [0, 1], "k--", label="Perfect calibration")
    ax1.plot(mean_pred, frac_pos, "o-", label=f"Model (ECE={ece:.4f})")
    ax1.set_xlabel("Mean predicted probability")
    ax1.set_ylabel("Fraction of positives")
    ax1.set_title("Calibration Curve")
    ax1.legend()

    # 信頼度ヒストグラム
    ax2.hist(y_proba, bins=n_bins, edgecolor="black", alpha=0.7)
    ax2.set_xlabel("Predicted probability")
    ax2.set_ylabel("Count")
    ax2.set_title("Confidence Histogram")

    plt.tight_layout()
    plt.savefig("calibration_analysis.png", dpi=150, bbox_inches="tight")
    plt.close()

    stats_df = pd.DataFrame(bin_stats)
    print(f"Calibration: ECE = {ece:.4f}, {n_bins} bins")
    return {"ece": ece, "bin_stats": stats_df, "fig": "calibration_analysis.png"}
```

---

## パイプライン統合

```
ensemble-methods → uncertainty-quantification → explainable-ai
  (アンサンブル)    (不確実性定量化)            (説明可能 AI)
       │                    │                      ↓
 active-learning ──────────┘              bayesian-statistics
  (能動学習)                                (ベイズ統計)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `conformal_intervals.csv` | Conformal 予測区間 | → reporting |
| `mc_dropout_uncertainty.csv` | MC Dropout 不確実性 | → active-learning |
| `calibration_analysis.png` | 校正曲線 | → presentation |
| `uncertainty_decomposition.json` | 分離結果 | → explainable-ai |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `papers_with_code` | Papers with Code | 不確実性定量化手法・ベンチマーク |
