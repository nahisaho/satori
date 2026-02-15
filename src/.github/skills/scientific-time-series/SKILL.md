---
name: scientific-time-series
description: |
  時系列解析・予測スキル。ARIMA/SARIMA/Prophet モデリング、変化点検出（PELT/Bayesian）、
  周期解析（FFT/ウェーブレット）、季節分解（STL）、異常検出、Granger 因果性検定の
  テンプレートを提供。実験データのトレンド解析・予測モデリングに適用。
tu_tools:
  - key: biotools
    name: bio.tools
    description: 時系列解析ツール検索
---

# Scientific Time Series Analysis

時系列データの分解・モデリング・予測・異常検出パイプライン。
プロセスモニタリング、環境測定、臨床バイタルサイン、実験時系列データなどに適用する。

## When to Use

- 時系列トレンドの分解・可視化が必要なとき
- ARIMA / Prophet による将来予測モデルを構築するとき
- 変化点や異常値をデータから検出するとき
- 周期性（季節性）を解析するとき
- Granger 因果性など時系列間の関係を調べるとき

---

## Quick Start

## 1. STL 季節分解

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.seasonal import STL

def stl_decomposition(series, period, robust=True, figsize=(12, 10)):
    """
    STL (Seasonal and Trend decomposition using Loess) による分解。

    Components:
        - Trend: 長期トレンド
        - Seasonal: 周期成分
        - Residual: 残差

    Parameters:
        series: pd.Series with DatetimeIndex
        period: 季節周期（データポイント数）
    """
    stl = STL(series, period=period, robust=robust)
    result = stl.fit()

    fig, axes = plt.subplots(4, 1, figsize=figsize, sharex=True)
    components = [("Observed", series), ("Trend", result.trend),
                  ("Seasonal", result.seasonal), ("Residual", result.resid)]

    for ax, (name, data) in zip(axes, components):
        ax.plot(data, linewidth=1)
        ax.set_ylabel(name)
        ax.grid(alpha=0.3)

    plt.suptitle("STL Decomposition", fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig("figures/stl_decomposition.png", dpi=300, bbox_inches="tight")
    plt.close()

    return result
```

## 2. ARIMA / SARIMA モデリング

```python
from statsmodels.tsa.statespace.sarimax import SARIMAX
from statsmodels.tsa.stattools import adfuller

def adf_stationarity_test(series, significance=0.05):
    """ADF 検定による定常性チェック。"""
    result = adfuller(series.dropna())
    return {
        "adf_statistic": result[0],
        "p_value": result[1],
        "used_lag": result[2],
        "n_obs": result[3],
        "is_stationary": result[1] < significance,
        "critical_values": result[4],
    }


def fit_sarima(series, order=(1, 1, 1), seasonal_order=(1, 1, 1, 12),
               forecast_steps=24):
    """
    SARIMA モデルのフィッティングと予測。

    Parameters:
        order: (p, d, q) — AR次数, 差分次数, MA次数
        seasonal_order: (P, D, Q, s) — 季節 AR/差分/MA 次数 + 周期
        forecast_steps: 予測期間数
    """
    model = SARIMAX(series, order=order, seasonal_order=seasonal_order,
                     enforce_stationarity=False, enforce_invertibility=False)
    fitted = model.fit(disp=False)

    # 予測
    forecast = fitted.get_forecast(steps=forecast_steps)
    pred_mean = forecast.predicted_mean
    conf_int = forecast.conf_int()

    # 診断
    diagnostics = {
        "aic": fitted.aic,
        "bic": fitted.bic,
        "ljung_box_p": fitted.test_serial_correlation("ljungbox")[0][0, 1],
    }

    return fitted, pred_mean, conf_int, diagnostics


def plot_forecast(series, pred_mean, conf_int, title="SARIMA Forecast",
                   figsize=(12, 5)):
    """時系列の実測値と予測を描画する。"""
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(series.index, series.values, "b-", label="Observed", linewidth=1)
    ax.plot(pred_mean.index, pred_mean.values, "r--", label="Forecast",
            linewidth=2)
    ax.fill_between(conf_int.index, conf_int.iloc[:, 0], conf_int.iloc[:, 1],
                    color="red", alpha=0.1, label="95% CI")
    ax.set_title(title, fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("figures/time_series_forecast.png", dpi=300, bbox_inches="tight")
    plt.close()
```

## 3. 変化点検出

```python
def detect_changepoints(series, method="pelt", penalty="bic", min_size=10):
    """
    変化点検出。

    method:
        "pelt"    — PELT アルゴリズム (ruptures)
        "cusum"   — CUSUM ベース
        "bayesian" — ベイズオンライン変化点検出

    Returns:
        list of changepoint indices
    """
    import ruptures as rpt

    signal = series.values

    if method == "pelt":
        algo = rpt.Pelt(model="rbf", min_size=min_size).fit(signal)
        cps = algo.predict(pen=10)
    elif method == "cusum":
        algo = rpt.Binseg(model="l2", min_size=min_size).fit(signal)
        cps = algo.predict(n_bkps=5)
    elif method == "bayesian":
        algo = rpt.BottomUp(model="l2", min_size=min_size).fit(signal)
        cps = algo.predict(n_bkps=5)

    # 最後の要素 (n) を除去
    cps = [cp for cp in cps if cp < len(signal)]

    return cps


def plot_changepoints(series, changepoints, figsize=(12, 4)):
    """変化点を時系列上に可視化する。"""
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(series.index, series.values, "b-", linewidth=1)
    for cp in changepoints:
        ax.axvline(series.index[cp], color="red", linestyle="--",
                   alpha=0.7, linewidth=1.5)
    ax.set_title(f"Changepoint Detection ({len(changepoints)} points found)",
                fontweight="bold")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("figures/changepoints.png", dpi=300, bbox_inches="tight")
    plt.close()
```

## 4. 周期解析（FFT / ウェーブレット）

```python
def fft_periodicity(series, fs=1.0, top_n=5, figsize=(10, 5)):
    """
    FFT によるパワースペクトルと支配的周期の抽出。
    """
    signal = series.values - np.mean(series.values)
    N = len(signal)
    fft_vals = np.fft.rfft(signal)
    fft_power = np.abs(fft_vals)**2
    freqs = np.fft.rfftfreq(N, d=1/fs)

    # DC成分除去
    fft_power[0] = 0

    # Top-N 周期
    top_idx = np.argsort(fft_power)[::-1][:top_n]
    dominant_periods = []
    for idx in top_idx:
        if freqs[idx] > 0:
            dominant_periods.append({
                "frequency": freqs[idx],
                "period": 1 / freqs[idx],
                "power": fft_power[idx],
            })

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(freqs[1:], fft_power[1:], "b-", linewidth=1)
    for dp in dominant_periods:
        ax.axvline(dp["frequency"], color="red", linestyle="--", alpha=0.5)
    ax.set_xlabel("Frequency")
    ax.set_ylabel("Power")
    ax.set_title("FFT Power Spectrum", fontweight="bold")
    plt.tight_layout()
    plt.savefig("figures/fft_spectrum.png", dpi=300, bbox_inches="tight")
    plt.close()

    return dominant_periods
```

## 5. Granger 因果性検定

```python
from statsmodels.tsa.stattools import grangercausalitytests

def granger_causality(df, cause_col, effect_col, max_lag=10,
                       significance=0.05):
    """
    Granger 因果性検定: cause → effect の因果関係を検定する。

    Returns:
        dict with optimal_lag, p_values, is_causal
    """
    data = df[[effect_col, cause_col]].dropna()
    results = grangercausalitytests(data, maxlag=max_lag, verbose=False)

    p_values = {}
    for lag in range(1, max_lag + 1):
        p_val = results[lag][0]["ssr_ftest"][1]
        p_values[lag] = p_val

    optimal_lag = min(p_values, key=p_values.get)

    return {
        "cause": cause_col,
        "effect": effect_col,
        "optimal_lag": optimal_lag,
        "p_value_at_optimal": p_values[optimal_lag],
        "is_causal": p_values[optimal_lag] < significance,
        "all_p_values": p_values,
    }
```

## 6. 異常検出

```python
def anomaly_detection_zscore(series, window=30, threshold=3.0):
    """
    移動平均ベースの Z-score 異常検出。
    """
    rolling_mean = series.rolling(window=window, center=True).mean()
    rolling_std = series.rolling(window=window, center=True).std()
    z_scores = (series - rolling_mean) / (rolling_std + 1e-10)
    anomalies = series[np.abs(z_scores) > threshold]
    return anomalies, z_scores
```

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `biotools` | bio.tools | 時系列解析ツール検索 |

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/stationarity_test.csv` | CSV |
| `results/forecast_results.csv` | CSV |
| `results/changepoints.csv` | CSV |
| `results/dominant_periods.csv` | CSV |
| `results/granger_causality.csv` | CSV |
| `figures/stl_decomposition.png` | PNG |
| `figures/time_series_forecast.png` | PNG |
| `figures/changepoints.png` | PNG |
| `figures/fft_spectrum.png` | PNG |

#### 依存パッケージ

```
statsmodels>=0.14
ruptures>=1.1
```
