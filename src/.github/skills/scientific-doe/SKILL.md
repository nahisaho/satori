---
name: scientific-doe
description: |
  実験計画法（DOE）スキル。直交配列表（L9/L16/L27）、中心複合計画（CCD）、
  Box-Behnken 設計、D-最適計画、応答曲面法（RSM）、交互作用解析、
  ベイズ最適化（Gaussian Process）、効果プロット（主効果/交互作用/pareto）の
  テンプレートを提供。
---

# Scientific Design of Experiments (DOE)

体系的な実験計画と最適化のためのスキル。直交表による因子スクリーニングから
RSM による最適条件探索、ベイズ最適化による逐次最適化まで、実験の各段階に
対応するテンプレートを提供する。

## When to Use

- 多因子実験の計画（因子・水準の設計）が必要なとき
- 直交表やCCD で実験回数を最小化したいとき
- 主効果・交互作用の寄与率を定量化するとき
- 応答曲面で最適条件を探索するとき
- ベイズ最適化で逐次実験を行いたいとき

---

## Quick Start

## 1. 因子設計テンプレート

```python
import numpy as np
import pandas as pd

def define_factors(factor_dict):
    """
    因子定義テンプレート。

    factor_dict 例:
    {
        "Temperature": {"levels": [200, 250, 300], "unit": "°C", "type": "continuous"},
        "Pressure": {"levels": [1, 5, 10], "unit": "mTorr", "type": "continuous"},
        "Gas_Ratio": {"levels": [0.2, 0.5, 0.8], "unit": "-", "type": "continuous"},
        "Material": {"levels": ["ZnO", "ITO", "TiO2"], "unit": "-", "type": "categorical"},
    }
    """
    summary = pd.DataFrame([
        {"Factor": k, "Levels": len(v["levels"]), "Values": str(v["levels"]),
         "Unit": v["unit"], "Type": v["type"]}
        for k, v in factor_dict.items()
    ])
    print("=== Factor Design ===")
    print(summary.to_string(index=False))
    return summary
```

## 2. 直交配列表

```python
# 田口 L9 直交表 (3 因子 × 3 水準)
L9 = np.array([
    [0, 0, 0],
    [0, 1, 1],
    [0, 2, 2],
    [1, 0, 1],
    [1, 1, 2],
    [1, 2, 0],
    [2, 0, 2],
    [2, 1, 0],
    [2, 2, 1],
])

def generate_taguchi_design(factor_dict, array="L9"):
    """
    田口直交表から実験計画を生成する。

    Available arrays: L4(2^3), L9(3^3-4), L16(2^15), L27(3^13)
    """
    arrays = {
        "L9": L9,
        "L4": np.array([[0,0,0],[0,1,1],[1,0,1],[1,1,0]]),
    }
    oa = arrays.get(array, L9)
    factors = list(factor_dict.keys())

    runs = []
    for row in oa:
        run = {}
        for i, factor in enumerate(factors[:oa.shape[1]]):
            levels = factor_dict[factor]["levels"]
            run[factor] = levels[row[i] % len(levels)]
        runs.append(run)

    design_df = pd.DataFrame(runs)
    design_df.index.name = "Run"
    design_df.index += 1
    return design_df
```

## 3. 中心複合計画 (CCD)

```python
from itertools import product

def central_composite_design(factor_dict, alpha="rotatable", center_points=3):
    """
    中心複合計画 (Central Composite Design) を生成する。

    Components:
        - 2^k 完全実施要因計画 (cube points)
        - 2k 軸点 (axial/star points)
        - n_c 中心点

    alpha:
        "rotatable"  — α = (2^k)^(1/4) (回転可能)
        "face"       — α = 1 (面心)
        float        — 任意の値
    """
    continuous_factors = {k: v for k, v in factor_dict.items()
                         if v["type"] == "continuous"}
    factor_names = list(continuous_factors.keys())
    k = len(factor_names)

    if alpha == "rotatable":
        alpha_val = (2 ** k) ** 0.25
    elif alpha == "face":
        alpha_val = 1.0
    else:
        alpha_val = float(alpha)

    # コード化: -1, 0, +1
    midpoints = {}
    half_ranges = {}
    for name, info in continuous_factors.items():
        levels = info["levels"]
        mid = (max(levels) + min(levels)) / 2
        half = (max(levels) - min(levels)) / 2
        midpoints[name] = mid
        half_ranges[name] = half

    runs = []

    # Cube points (2^k)
    for combo in product([-1, 1], repeat=k):
        run = {factor_names[i]: midpoints[factor_names[i]] + combo[i] * half_ranges[factor_names[i]]
               for i in range(k)}
        run["_type"] = "cube"
        runs.append(run)

    # Axial points (2k)
    for i in range(k):
        for direction in [-1, 1]:
            run = {name: midpoints[name] for name in factor_names}
            run[factor_names[i]] = midpoints[factor_names[i]] + direction * alpha_val * half_ranges[factor_names[i]]
            run["_type"] = "axial"
            runs.append(run)

    # Center points
    for _ in range(center_points):
        run = {name: midpoints[name] for name in factor_names}
        run["_type"] = "center"
        runs.append(run)

    design_df = pd.DataFrame(runs)
    design_df.index.name = "Run"
    design_df.index += 1
    return design_df
```

## 4. 分散分析 (ANOVA) — 因子効果解析

```python
from scipy.stats import f_oneway

def anova_factor_effects(design_df, response_col, factor_cols):
    """
    各因子の主効果を ANOVA で評価する。

    Returns:
        DataFrame with Factor, SS, DF, MS, F_value, p_value, contribution_pct
    """
    ss_total = np.sum((design_df[response_col] - design_df[response_col].mean())**2)
    results = []

    for factor in factor_cols:
        groups = [group[response_col].values
                  for _, group in design_df.groupby(factor)]
        if len(groups) < 2:
            continue
        f_val, p_val = f_oneway(*groups)

        # SS_factor
        grand_mean = design_df[response_col].mean()
        ss_factor = sum(len(g) * (np.mean(g) - grand_mean)**2 for g in groups)
        df_factor = len(groups) - 1
        ms_factor = ss_factor / df_factor

        results.append({
            "Factor": factor,
            "SS": ss_factor,
            "DF": df_factor,
            "MS": ms_factor,
            "F_value": f_val,
            "p_value": p_val,
            "Contribution_pct": ss_factor / ss_total * 100 if ss_total > 0 else 0,
        })

    return pd.DataFrame(results).sort_values("Contribution_pct", ascending=False)
```

## 5. 主効果プロット

```python
import matplotlib.pyplot as plt

def main_effects_plot(design_df, response_col, factor_cols, figsize=None):
    """全因子の主効果プロットを描画する。"""
    n = len(factor_cols)
    if figsize is None:
        figsize = (4 * n, 4)

    fig, axes = plt.subplots(1, n, figsize=figsize, sharey=True)
    if n == 1:
        axes = [axes]

    grand_mean = design_df[response_col].mean()

    for ax, factor in zip(axes, factor_cols):
        means = design_df.groupby(factor)[response_col].mean()
        ax.plot(range(len(means)), means.values, "bo-", linewidth=2, markersize=8)
        ax.axhline(grand_mean, color="gray", linestyle="--", alpha=0.5)
        ax.set_xticks(range(len(means)))
        ax.set_xticklabels(means.index, rotation=45)
        ax.set_xlabel(factor)
        ax.grid(alpha=0.3)

    axes[0].set_ylabel(response_col)
    plt.suptitle("Main Effects Plot", fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig("figures/main_effects_plot.png", dpi=300, bbox_inches="tight")
    plt.close()
```

## 6. ベイズ最適化（Gaussian Process）

```python
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern

def bayesian_optimization(objective_func, bounds, n_initial=5,
                           n_iterations=20, kappa=2.576):
    """
    ベイズ最適化（Gaussian Process + Expected Improvement）。

    Parameters:
        objective_func: callable f(x) → y (最大化)
        bounds: dict {"param": (low, high)}
        n_initial: 初期ランダムサンプリング数
        n_iterations: 最適化ステップ数
        kappa: 探索-活用トレードオフ (UCB の κ)
    """
    from scipy.optimize import minimize as scipy_minimize
    from scipy.stats import norm

    param_names = list(bounds.keys())
    lows = np.array([bounds[p][0] for p in param_names])
    highs = np.array([bounds[p][1] for p in param_names])

    # 初期サンプリング
    X_init = np.random.uniform(lows, highs, size=(n_initial, len(param_names)))
    y_init = np.array([objective_func(dict(zip(param_names, x))) for x in X_init])

    X_observed = X_init.tolist()
    y_observed = y_init.tolist()

    gp = GaussianProcessRegressor(kernel=Matern(nu=2.5), n_restarts_optimizer=5,
                                   random_state=42)

    for i in range(n_iterations):
        X_arr = np.array(X_observed)
        y_arr = np.array(y_observed)
        gp.fit(X_arr, y_arr)

        # UCB acquisition function
        def neg_ucb(x):
            mu, sigma = gp.predict(x.reshape(1, -1), return_std=True)
            return -(mu + kappa * sigma)

        # 複数の開始点から最適化
        best_x = None
        best_val = float("inf")
        for _ in range(10):
            x0 = np.random.uniform(lows, highs)
            res = scipy_minimize(neg_ucb, x0, bounds=list(zip(lows, highs)),
                                method="L-BFGS-B")
            if res.fun < best_val:
                best_val = res.fun
                best_x = res.x

        # 新しい点を評価
        y_new = objective_func(dict(zip(param_names, best_x)))
        X_observed.append(best_x.tolist())
        y_observed.append(y_new)

    # 最適解
    best_idx = np.argmax(y_observed)
    best_params = dict(zip(param_names, X_observed[best_idx]))
    best_y = y_observed[best_idx]

    return {
        "best_params": best_params,
        "best_value": best_y,
        "X_history": np.array(X_observed),
        "y_history": np.array(y_observed),
        "gp_model": gp,
    }
```

## 7. 交互作用プロット

```python
def interaction_plot(design_df, response_col, factor1, factor2, figsize=(8, 6)):
    """2 因子間の交互作用プロットを描画する。"""
    fig, ax = plt.subplots(figsize=figsize)

    for level2, group in design_df.groupby(factor2):
        means = group.groupby(factor1)[response_col].mean()
        ax.plot(range(len(means)), means.values, "o-", linewidth=2,
                markersize=8, label=f"{factor2}={level2}")

    ax.set_xticks(range(len(means)))
    ax.set_xticklabels(means.index)
    ax.set_xlabel(factor1)
    ax.set_ylabel(response_col)
    ax.set_title(f"Interaction Plot: {factor1} × {factor2}", fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("figures/interaction_plot.png", dpi=300, bbox_inches="tight")
    plt.close()
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/experimental_design.csv` | CSV |
| `results/anova_factor_effects.csv` | CSV |
| `results/bayesian_optimization_history.csv` | CSV |
| `figures/main_effects_plot.png` | PNG |
| `figures/interaction_plot.png` | PNG |
| `figures/bayesian_convergence.png` | PNG |

#### 依存パッケージ

```
scipy>=1.10
scikit-learn>=1.3
```
