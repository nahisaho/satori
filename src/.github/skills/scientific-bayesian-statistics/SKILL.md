---
name: scientific-bayesian-statistics
description: |
  ベイズ統計スキル。PyMC・Stan・ArviZ を活用し、ベイズ回帰・階層モデル・
  MCMC サンプリング・ベイズ最適化・事後予測チェック・モデル比較を支援。
  「ベイズ回帰して」「MCMC で推定して」「事後分布を求めて」で発火。
tu_tools:
  - key: biotools
    name: bio.tools
    description: ベイズ統計ツールレジストリ検索
---

# Scientific Bayesian Statistics

ベイズ統計推論のための解析スキル。確率的プログラミング
フレームワーク（PyMC, Stan, NumPyro）を活用した
パラメータ推定・モデル比較・意思決定支援を行う。

## When to Use

- パラメータの事後分布推定
- 階層ベイズモデル（マルチレベル・混合効果）
- ベイズ回帰・ベイズ分類
- MCMC 診断（Rhat, ESS, トレースプロット）
- モデル比較（WAIC, LOO-CV, ベイズファクター）
- ベイズ最適化（Gaussian Process）
- 事後予測チェック

## Quick Start

### ベイズ解析パイプライン

```
Phase 1: Model Specification
  - 尤度関数の選択 (Normal, Poisson, Binomial, etc.)
  - 事前分布の設定 (弱情報事前分布推奨)
  - パラメータ空間の定義
    ↓
Phase 2: Prior Predictive Check
  - 事前予測シミュレーション
  - 事前分布の妥当性検証
  - ドメイン知識との整合性確認
    ↓
Phase 3: MCMC Sampling
  - NUTS (No-U-Turn Sampler) 実行
  - チェーン数・サンプル数設定 (chains≥4, draws≥2000)
  - 収束診断 (Rhat < 1.01, ESS > 400)
    ↓
Phase 4: Posterior Analysis
  - 事後分布可視化 (forest, trace, pair plot)
  - 95% HDI (Highest Density Interval)
  - 事後予測チェック (PPC)
    ↓
Phase 5: Model Comparison
  - WAIC / LOO-CV 算出
  - Bayes Factor 計算
  - モデル重み付け平均 (Stacking)
    ↓
Phase 6: Reporting
  - 推定値サマリーテーブル
  - 事後分布プロット
  - 結果の解釈・レポート生成
```

## Workflow

### 1. PyMC: ベイズ線形回帰

```python
import pymc as pm
import arviz as az
import numpy as np
import matplotlib.pyplot as plt

# === データ生成 (例) ===
np.random.seed(42)
N = 100
X = np.random.randn(N)
true_alpha, true_beta, true_sigma = 1.0, 2.5, 0.5
y = true_alpha + true_beta * X + np.random.normal(0, true_sigma, N)

# === ベイズ線形回帰 ===
with pm.Model() as linear_model:
    # 事前分布
    alpha = pm.Normal("alpha", mu=0, sigma=10)
    beta = pm.Normal("beta", mu=0, sigma=10)
    sigma = pm.HalfNormal("sigma", sigma=5)

    # 尤度
    mu = alpha + beta * X
    y_obs = pm.Normal("y_obs", mu=mu, sigma=sigma, observed=y)

    # MCMC サンプリング (NUTS)
    trace = pm.sample(
        draws=2000,
        tune=1000,
        chains=4,
        cores=4,
        target_accept=0.95,
        return_inferencedata=True,
    )

# === 収束診断 ===
print(az.summary(trace, var_names=["alpha", "beta", "sigma"]))
# Rhat < 1.01 && ESS > 400 を確認

# === トレースプロット ===
az.plot_trace(trace, var_names=["alpha", "beta", "sigma"])
plt.savefig("figures/bayesian_trace.png", dpi=300, bbox_inches="tight")
plt.show()
```

### 2. 階層ベイズモデル（マルチレベル）

```python
import pandas as pd

# === 階層ベイズ回帰 (ランダム効果) ===
# 例: 複数学校の生徒の成績予測
# school_id: グループ変数
# x: 予測変数, y: 目的変数

def hierarchical_model(df):
    school_ids = df["school_id"].unique()
    school_idx = df["school_id"].map({s: i for i, s in enumerate(school_ids)}).values
    n_schools = len(school_ids)

    with pm.Model() as hier_model:
        # Hyper-priors (全体平均)
        mu_alpha = pm.Normal("mu_alpha", mu=0, sigma=10)
        sigma_alpha = pm.HalfNormal("sigma_alpha", sigma=5)
        mu_beta = pm.Normal("mu_beta", mu=0, sigma=10)
        sigma_beta = pm.HalfNormal("sigma_beta", sigma=5)

        # 学校レベルのパラメータ (部分プーリング)
        alpha = pm.Normal("alpha", mu=mu_alpha, sigma=sigma_alpha, shape=n_schools)
        beta = pm.Normal("beta", mu=mu_beta, sigma=sigma_beta, shape=n_schools)

        # 残差
        sigma = pm.HalfNormal("sigma", sigma=5)

        # 尤度
        mu = alpha[school_idx] + beta[school_idx] * df["x"].values
        y_obs = pm.Normal("y_obs", mu=mu, sigma=sigma, observed=df["y"].values)

        # サンプリング
        trace = pm.sample(2000, tune=1000, chains=4, target_accept=0.95,
                          return_inferencedata=True)

    return trace

# === Forest Plot（学校間パラメータ比較）===
# az.plot_forest(trace, var_names=["alpha", "beta"], combined=True)
```

### 3. 事後予測チェック (PPC)

```python
def posterior_predictive_check(model, trace, observed_y):
    """事後予測チェック: モデルの適合度検証"""
    with model:
        ppc = pm.sample_posterior_predictive(trace, random_seed=42)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # 1. 密度オーバーレイ
    az.plot_ppc(az.from_pymc3(posterior_predictive=ppc, model=model),
                ax=axes[0], kind="kde")
    axes[0].set_title("Posterior Predictive: KDE")

    # 2. 統計量比較 (mean)
    ppc_mean = ppc.posterior_predictive["y_obs"].mean(dim=("chain", "draw")).values
    axes[1].hist(ppc_mean, bins=30, alpha=0.7, label="PPC mean")
    axes[1].axvline(observed_y.mean(), color="red", linestyle="--", label="Observed mean")
    axes[1].set_title("PPC: Mean Comparison")
    axes[1].legend()

    # 3. 残差チェック
    ppc_median = np.median(ppc.posterior_predictive["y_obs"].values, axis=(0, 1))
    residuals = observed_y - ppc_median
    axes[2].scatter(ppc_median, residuals, alpha=0.5, s=10)
    axes[2].axhline(0, color="red", linestyle="--")
    axes[2].set_title("PPC Residuals")

    plt.tight_layout()
    plt.savefig("figures/bayesian_ppc.png", dpi=300, bbox_inches="tight")
    plt.show()
```

### 4. モデル比較 (WAIC / LOO-CV)

```python
def bayesian_model_comparison(models_dict):
    """
    WAIC / LOO-CV でベイズモデルを比較。
    models_dict: {"model_name": (model, trace)} の辞書
    """
    comparison_data = {}

    for name, (model, trace) in models_dict.items():
        # WAIC
        waic = az.waic(trace)
        # LOO-CV (PSIS-LOO)
        loo = az.loo(trace)

        comparison_data[name] = {
            "waic": round(waic.waic, 2),
            "waic_se": round(waic.waic_se, 2),
            "p_waic": round(waic.p_waic, 2),
            "loo": round(loo.loo, 2),
            "loo_se": round(loo.loo_se, 2),
            "p_loo": round(loo.p_loo, 2),
        }

    # ArviZ モデル比較
    comp_df = az.compare(
        {name: trace for name, (_, trace) in models_dict.items()},
        ic="loo",
    )
    print("Model Comparison (LOO-CV):")
    print(comp_df)

    # 比較プロット
    az.plot_compare(comp_df)
    plt.savefig("figures/bayesian_model_comparison.png", dpi=300, bbox_inches="tight")
    plt.show()

    return comparison_data
```

### 5. ベイズ最適化 (Gaussian Process)

```python
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern
from scipy.optimize import minimize
from scipy.stats import norm

def bayesian_optimization(objective_func, bounds, n_init=5, n_iter=25):
    """
    Gaussian Process ベースのベイズ最適化。
    Expected Improvement (EI) を獲得関数として使用。
    """
    dim = len(bounds)

    # 初期サンプリング (Latin Hypercube)
    X_samples = np.random.uniform(
        [b[0] for b in bounds],
        [b[1] for b in bounds],
        size=(n_init, dim)
    )
    y_samples = np.array([objective_func(x) for x in X_samples])

    # GP モデル
    kernel = Matern(nu=2.5)
    gp = GaussianProcessRegressor(kernel=kernel, alpha=1e-6, normalize_y=True)

    history = {"X": list(X_samples), "y": list(y_samples)}

    for i in range(n_iter):
        gp.fit(np.array(history["X"]), np.array(history["y"]))

        # Expected Improvement
        def neg_ei(x):
            mu, sigma = gp.predict(x.reshape(1, -1), return_std=True)
            best_y = np.min(history["y"])
            z = (best_y - mu) / (sigma + 1e-8)
            ei = (best_y - mu) * norm.cdf(z) + sigma * norm.pdf(z)
            return -ei[0]

        # 次の探索点
        best_x = None
        best_ei = np.inf
        for _ in range(20):
            x0 = np.random.uniform([b[0] for b in bounds], [b[1] for b in bounds])
            result = minimize(neg_ei, x0, bounds=bounds, method="L-BFGS-B")
            if result.fun < best_ei:
                best_ei = result.fun
                best_x = result.x

        new_y = objective_func(best_x)
        history["X"].append(best_x)
        history["y"].append(new_y)

        print(f"Iter {i+1}: x = {best_x}, y = {new_y:.4f}, best = {min(history['y']):.4f}")

    return history
```

### 6. Stan によるモデリング

```python
# === CmdStanPy インターフェース ===
import cmdstanpy

def fit_stan_model(stan_code, data_dict, **kwargs):
    """Stan モデルのコンパイルとフィット"""
    model = cmdstanpy.CmdStanModel(stan_file=None, stan_file_path=None,
                                     model_code=stan_code)
    fit = model.sample(data=data_dict, chains=4, iter_sampling=2000,
                       iter_warmup=1000, **kwargs)

    # ArviZ 変換
    idata = az.from_cmdstanpy(fit)
    return idata

# Stan モデルコード例: 階層正規モデル
STAN_HIERARCHICAL = """
data {
  int<lower=0> N;
  int<lower=0> J;
  array[N] int<lower=1,upper=J> group;
  vector[N] x;
  vector[N] y;
}
parameters {
  real mu_alpha;
  real mu_beta;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma;
  vector[J] alpha;
  vector[J] beta;
}
model {
  mu_alpha ~ normal(0, 10);
  mu_beta ~ normal(0, 10);
  sigma_alpha ~ half_normal(5);
  sigma_beta ~ half_normal(5);
  sigma ~ half_normal(5);
  alpha ~ normal(mu_alpha, sigma_alpha);
  beta ~ normal(mu_beta, sigma_beta);
  y ~ normal(alpha[group] + beta[group] .* x, sigma);
}
"""
```

---

## Best Practices

1. **弱情報事前分布**: データに語らせる。過度に情報的な事前分布は避ける
2. **収束診断必須**: Rhat < 1.01, ESS > 400/chain, 発散なし
3. **事後予測チェック**: モデルの適合度を必ず PPC で検証
4. **LOO-CV 優先**: モデル比較は WAIC より LOO-CV (PSIS-LOO) が安定
5. **中心化パラメトリゼーション**: 階層モデルでは non-centered parameterization を検討
6. **再現性**: random_seed を固定
7. **感度分析**: 事前分布を変えた結果の安定性を確認

## Completeness Checklist

- [ ] モデル仕様 (尤度・事前分布) の定義
- [ ] Prior predictive check
- [ ] MCMC サンプリング完了
- [ ] 収束診断 (Rhat, ESS, divergences)
- [ ] 事後分布の要約 (mean, HDI)
- [ ] 事後予測チェック (PPC)
- [ ] モデル比較 (LOO-CV / WAIC)
- [ ] 結果レポート・プロット生成

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `biotools` | bio.tools | ベイズ統計ツールレジストリ検索 |

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/bayesian_summary.json` | 事後分布サマリー（JSON） | サンプリング完了時 |
| `figures/bayesian_trace.png` | トレースプロット | MCMC 完了時 |
| `figures/bayesian_ppc.png` | 事後予測チェック図 | PPC 完了時 |
| `figures/bayesian_model_comparison.png` | モデル比較プロット | 比較完了時 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-statistical-testing` | ← 頻度論検定との比較基盤 |
| `scientific-ml-regression` | ← 回帰モデルのベイズ拡張 |
| `scientific-doe` | ← 実験計画のベイズ最適化 |
| `scientific-causal-inference` | ← ベイズ因果モデル |
| `scientific-quantum-computing` | → VQE パラメータのベイズ最適化 |
| `scientific-meta-analysis` | ← ベイズメタ分析 |
