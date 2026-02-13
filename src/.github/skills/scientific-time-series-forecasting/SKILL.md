---
name: scientific-time-series-forecasting
description: |
  ML 時系列予測スキル。Prophet/NeuralProphet・N-BEATS・
  Temporal Fusion Transformer (TFT)・時系列特徴量エンジニアリング・
  バックテスト・多段階予測・アンサンブル予測。
---

# Scientific Time Series Forecasting

深層学習・ML ベースの時系列予測パイプラインを提供し、
Prophet から Transformer まで最新手法を網羅する。

## When to Use

- Prophet/NeuralProphet で季節性時系列を予測するとき
- 深層学習 (N-BEATS/TFT) で高精度予測するとき
- 時系列特徴量エンジニアリングでラグ・ローリング特徴を生成するとき
- バックテストで予測性能を厳密に評価するとき
- 複数モデルのアンサンブル予測をするとき
- 多変量・多段階予測をするとき

> **Note**: 古典時系列 (ARIMA/STL/FFT) は `scientific-time-series` を参照。

---

## Quick Start

## 1. Prophet / NeuralProphet

```python
import numpy as np
import pandas as pd


def prophet_forecast(df, date_col, value_col, periods=30,
                     freq="D", yearly=True, weekly=True,
                     changepoint_prior=0.05):
    """
    Prophet 時系列予測。

    Parameters:
        df: pd.DataFrame — 時系列データ
        date_col: str — 日付カラム
        value_col: str — 値カラム
        periods: int — 予測期間
        freq: str — 頻度 ("D" / "H" / "M")
        yearly: bool — 年次季節性
        weekly: bool — 週次季節性
        changepoint_prior: float — 変化点感度
    """
    from prophet import Prophet

    prophet_df = df[[date_col, value_col]].rename(
        columns={date_col: "ds", value_col: "y"})

    model = Prophet(
        yearly_seasonality=yearly,
        weekly_seasonality=weekly,
        changepoint_prior_scale=changepoint_prior)
    model.fit(prophet_df)

    future = model.make_future_dataframe(periods=periods, freq=freq)
    forecast = model.predict(future)

    # 評価
    merged = forecast.merge(prophet_df, on="ds", how="left")
    valid = merged.dropna(subset=["y"])
    mae = np.mean(np.abs(valid["y"] - valid["yhat"]))
    mape = np.mean(np.abs((valid["y"] - valid["yhat"]) / valid["y"])) * 100

    fig1 = model.plot(forecast)
    fig1.savefig("prophet_forecast.png", dpi=150, bbox_inches="tight")

    fig2 = model.plot_components(forecast)
    fig2.savefig("prophet_components.png", dpi=150, bbox_inches="tight")

    print(f"Prophet: {periods} periods, MAE={mae:.4f}, MAPE={mape:.1f}%")
    return {"forecast": forecast, "model": model,
            "mae": mae, "mape": mape}


def neuralprophet_forecast(df, date_col, value_col, periods=30,
                           n_lags=60, n_forecasts=30):
    """
    NeuralProphet 時系列予測 (AR-Net)。

    Parameters:
        df: pd.DataFrame — 時系列データ
        date_col: str — 日付カラム
        value_col: str — 値カラム
        periods: int — 予測期間
        n_lags: int — 自己回帰ラグ数
        n_forecasts: int — 多段階予測ステップ
    """
    from neuralprophet import NeuralProphet

    np_df = df[[date_col, value_col]].rename(
        columns={date_col: "ds", value_col: "y"})

    model = NeuralProphet(
        n_lags=n_lags, n_forecasts=n_forecasts,
        yearly_seasonality=True, weekly_seasonality=True,
        learning_rate=0.01, epochs=100)

    metrics = model.fit(np_df, freq="D")

    future = model.make_future_dataframe(np_df, periods=periods, n_historic_predictions=True)
    forecast = model.predict(future)

    fig = model.plot(forecast)
    fig.savefig("neuralprophet_forecast.png", dpi=150, bbox_inches="tight")

    print(f"NeuralProphet: lags={n_lags}, forecasts={n_forecasts}")
    return {"forecast": forecast, "model": model, "metrics": metrics}
```

## 2. 時系列特徴量エンジニアリング

```python
def create_ts_features(df, date_col, value_col,
                       lags=None, rolling_windows=None):
    """
    時系列特徴量エンジニアリング。

    Parameters:
        df: pd.DataFrame — 時系列データ
        date_col: str — 日付カラム
        value_col: str — 値カラム
        lags: list[int] | None — ラグ特徴量 (e.g., [1,7,14,28])
        rolling_windows: list[int] | None — ローリング窓 (e.g., [7,14,30])
    """
    if lags is None:
        lags = [1, 3, 7, 14, 28]
    if rolling_windows is None:
        rolling_windows = [7, 14, 30]

    result = df.copy()
    result[date_col] = pd.to_datetime(result[date_col])
    result = result.sort_values(date_col)

    # カレンダー特徴量
    result["dayofweek"] = result[date_col].dt.dayofweek
    result["dayofyear"] = result[date_col].dt.dayofyear
    result["month"] = result[date_col].dt.month
    result["quarter"] = result[date_col].dt.quarter
    result["is_weekend"] = (result[date_col].dt.dayofweek >= 5).astype(int)

    # 周期エンコーディング
    result["sin_day"] = np.sin(2 * np.pi * result["dayofyear"] / 365.25)
    result["cos_day"] = np.cos(2 * np.pi * result["dayofyear"] / 365.25)
    result["sin_week"] = np.sin(2 * np.pi * result["dayofweek"] / 7)
    result["cos_week"] = np.cos(2 * np.pi * result["dayofweek"] / 7)

    # ラグ特徴量
    for lag in lags:
        result[f"lag_{lag}"] = result[value_col].shift(lag)

    # ローリング統計量
    for window in rolling_windows:
        result[f"rolling_mean_{window}"] = result[value_col].rolling(window).mean()
        result[f"rolling_std_{window}"] = result[value_col].rolling(window).std()
        result[f"rolling_min_{window}"] = result[value_col].rolling(window).min()
        result[f"rolling_max_{window}"] = result[value_col].rolling(window).max()

    # 差分特徴量
    result["diff_1"] = result[value_col].diff(1)
    result["diff_7"] = result[value_col].diff(7)

    n_features = len(result.columns) - len(df.columns)
    print(f"TS Features: {n_features} features created "
          f"(lags={lags}, windows={rolling_windows})")
    return result


def ts_backtest(df, date_col, value_col, model_fn,
                n_splits=5, horizon=30, gap=0):
    """
    時系列バックテスト (Walk-forward validation)。

    Parameters:
        df: pd.DataFrame — 時系列データ
        date_col: str — 日付カラム
        value_col: str — 値カラム
        model_fn: callable — モデル学習・予測関数 (train_df → forecast_df)
        n_splits: int — 分割数
        horizon: int — 予測ホライズン
        gap: int — 学習-テスト間ギャップ
    """
    from sklearn.metrics import mean_absolute_error, mean_squared_error

    sorted_df = df.sort_values(date_col).reset_index(drop=True)
    n = len(sorted_df)
    fold_size = (n - horizon) // n_splits

    results = []

    for i in range(n_splits):
        train_end = fold_size * (i + 1)
        test_start = train_end + gap
        test_end = min(test_start + horizon, n)

        if test_end > n:
            break

        train_df = sorted_df.iloc[:train_end]
        test_df = sorted_df.iloc[test_start:test_end]

        forecast = model_fn(train_df)
        y_true = test_df[value_col].values[:len(forecast)]
        y_pred = forecast[:len(y_true)]

        mae = mean_absolute_error(y_true, y_pred)
        rmse = np.sqrt(mean_squared_error(y_true, y_pred))
        mape = np.mean(np.abs((y_true - y_pred) / (y_true + 1e-10))) * 100

        results.append({
            "fold": i, "train_size": train_end,
            "test_size": test_end - test_start,
            "mae": mae, "rmse": rmse, "mape": mape})

    results_df = pd.DataFrame(results)
    print(f"Backtest ({n_splits} folds, h={horizon}): "
          f"MAE={results_df['mae'].mean():.4f} ± {results_df['mae'].std():.4f}")
    return results_df
```

---

## パイプライン統合

```
time-series → time-series-forecasting → model-monitoring
  (古典解析)    (ML 予測)                 (監視)
       │              │                       ↓
 spectral-signal ────┘             anomaly-detection
  (周波数解析)                       (異常検知)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `prophet_forecast.png` | Prophet 予測結果 | → presentation |
| `ts_features.csv` | 時系列特徴量 | → ml-regression |
| `backtest_results.csv` | バックテスト結果 | → model selection |
