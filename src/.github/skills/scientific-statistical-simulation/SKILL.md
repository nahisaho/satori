---
name: scientific-statistical-simulation
description: |
  統計シミュレーションスキル。Monte Carlo 法・Bootstrap 推論・
  Permutation Test・統計的検出力分析・確率的リスク評価。
---

# Scientific Statistical Simulation

コンピュータベースの統計シミュレーションにより、
推論の不確実性定量化・検出力設計・リスク評価を行う。

## When to Use

- Monte Carlo シミュレーションで確率分布を推定するとき
- Bootstrap で信頼区間を算出するとき
- Permutation Test でノンパラメトリック検定を行うとき
- 実験前に必要なサンプルサイズを検出力分析で決めるとき
- リスクシナリオの確率的評価を行うとき

---

## Quick Start

## 1. Monte Carlo シミュレーション

```python
import numpy as np
import pandas as pd
from typing import Callable, Dict, Any


def monte_carlo_simulation(func, param_distributions,
                           n_simulations=10000, seed=42,
                           summary_quantiles=None):
    """
    Monte Carlo シミュレーション。

    Parameters:
        func: Callable — シミュレーション対象関数 (dict → float)
        param_distributions: dict — {param_name: scipy.stats distribution}
        n_simulations: int — シミュレーション回数
        seed: int — 乱数シード
        summary_quantiles: list[float] | None — サマリー分位点
    """
    rng = np.random.default_rng(seed)
    if summary_quantiles is None:
        summary_quantiles = [0.025, 0.25, 0.5, 0.75, 0.975]

    results = []
    param_samples = {}

    # パラメータサンプリング
    for name, dist in param_distributions.items():
        param_samples[name] = dist.rvs(
            size=n_simulations, random_state=rng)

    # シミュレーション実行
    for i in range(n_simulations):
        params = {name: samples[i]
                  for name, samples in param_samples.items()}
        results.append(func(params))

    results = np.array(results)

    # サマリー統計
    summary = {
        "mean": np.mean(results),
        "std": np.std(results),
        "min": np.min(results),
        "max": np.max(results),
    }
    for q in summary_quantiles:
        summary[f"q{q:.3f}"] = np.quantile(results, q)

    print(f"Monte Carlo ({n_simulations} runs):")
    print(f"  Mean={summary['mean']:.4f} ± {summary['std']:.4f}")
    print(f"  95% CI: [{summary.get('q0.025', 'N/A'):.4f}, "
          f"{summary.get('q0.975', 'N/A'):.4f}]")

    return {"results": results, "summary": summary,
            "param_samples": param_samples}
```

## 2. Bootstrap 推論

```python
def bootstrap_inference(data, statistic_fn, n_bootstrap=10000,
                        confidence_level=0.95, method="bca",
                        seed=42):
    """
    Bootstrap 信頼区間推定。

    Parameters:
        data: np.ndarray — データ
        statistic_fn: Callable — 統計量計算関数 (data → float)
        n_bootstrap: int — Bootstrap 回数
        confidence_level: float — 信頼水準
        method: str — "percentile" / "bca" / "basic"
        seed: int — 乱数シード
    """
    rng = np.random.default_rng(seed)
    n = len(data)
    observed = statistic_fn(data)

    # Bootstrap 標本
    boot_stats = np.array([
        statistic_fn(data[rng.integers(0, n, size=n)])
        for _ in range(n_bootstrap)])

    alpha = 1 - confidence_level
    if method == "percentile":
        ci_low = np.quantile(boot_stats, alpha / 2)
        ci_high = np.quantile(boot_stats, 1 - alpha / 2)
    elif method == "bca":
        from scipy import stats as sp_stats
        # Bias correction
        z0 = sp_stats.norm.ppf(np.mean(boot_stats < observed))
        # Acceleration (jackknife)
        jack_stats = np.array([
            statistic_fn(np.delete(data, i)) for i in range(n)])
        jack_mean = jack_stats.mean()
        a = (np.sum((jack_mean - jack_stats) ** 3) /
             (6 * np.sum((jack_mean - jack_stats) ** 2) ** 1.5 + 1e-10))
        z_alpha = sp_stats.norm.ppf([alpha / 2, 1 - alpha / 2])
        adjusted = sp_stats.norm.cdf(
            z0 + (z0 + z_alpha) / (1 - a * (z0 + z_alpha)))
        ci_low = np.quantile(boot_stats, adjusted[0])
        ci_high = np.quantile(boot_stats, adjusted[1])
    else:  # basic
        ci_low = 2 * observed - np.quantile(boot_stats, 1 - alpha / 2)
        ci_high = 2 * observed - np.quantile(boot_stats, alpha / 2)

    result = {
        "observed": observed,
        "boot_mean": np.mean(boot_stats),
        "boot_std": np.std(boot_stats),
        "ci_low": ci_low, "ci_high": ci_high,
        "method": method,
        "confidence_level": confidence_level,
    }
    print(f"Bootstrap ({method}): {observed:.4f} "
          f"[{ci_low:.4f}, {ci_high:.4f}] ({confidence_level:.0%} CI)")
    return result, boot_stats
```

## 3. 統計的検出力分析

```python
def power_analysis(effect_size_range=None, n_range=None,
                   alpha=0.05, test_type="two-sample-t",
                   n_simulations=5000, seed=42):
    """
    シミュレーションベース統計的検出力分析。

    Parameters:
        effect_size_range: list[float] | None — 効果量の範囲
        n_range: list[int] | None — サンプルサイズの範囲
        alpha: float — 有意水準
        test_type: str — "two-sample-t" / "paired-t" / "chi-square"
        n_simulations: int — 各条件のシミュレーション回数
        seed: int — 乱数シード
    """
    from scipy import stats as sp_stats

    rng = np.random.default_rng(seed)
    if effect_size_range is None:
        effect_size_range = [0.2, 0.5, 0.8]
    if n_range is None:
        n_range = [10, 20, 30, 50, 100, 200]

    records = []
    for es in effect_size_range:
        for n in n_range:
            rejections = 0
            for _ in range(n_simulations):
                if test_type == "two-sample-t":
                    x = rng.normal(0, 1, n)
                    y = rng.normal(es, 1, n)
                    _, p = sp_stats.ttest_ind(x, y)
                elif test_type == "paired-t":
                    diff = rng.normal(es, 1, n)
                    _, p = sp_stats.ttest_1samp(diff, 0)
                else:
                    raise ValueError(f"Unknown test: {test_type}")

                if p < alpha:
                    rejections += 1

            power = rejections / n_simulations
            records.append({
                "effect_size": es, "n": n,
                "power": power, "alpha": alpha,
            })

    df = pd.DataFrame(records)
    for es in effect_size_range:
        sub = df[df["effect_size"] == es]
        adequate = sub[sub["power"] >= 0.8]
        if len(adequate) > 0:
            min_n = adequate["n"].min()
            print(f"d={es}: min N={min_n} for power≥0.80")
        else:
            print(f"d={es}: power<0.80 for all tested N")

    return df
```

---

## パイプライン統合

```
[仮説設計] → statistical-simulation → statistical-testing
               (シミュレーション)        (本解析)
                      │
               doe ← adaptive-experiments
           (実験計画)   (適応実験)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `mc_results.npz` | Monte Carlo 結果 | → リスク評価 |
| `bootstrap_ci.csv` | 信頼区間 | → 統計レポート |
| `power_analysis.csv` | 検出力カーブ | → DOE サンプルサイズ |
