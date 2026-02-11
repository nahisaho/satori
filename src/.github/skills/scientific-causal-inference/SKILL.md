---
name: scientific-causal-inference
description: |
  因果推論スキル。傾向スコアマッチング（PSM）、逆確率重み付け（IPW / IPTW）、
  操作変数法（2SLS）、差分の差分法（DID）、回帰不連続デザイン（RDD）、
  DAG ベースの共変量選択（backdoor criterion）、感度分析テンプレートを提供。
---

# Scientific Causal Inference

観察データから因果効果を推定するための統計的手法パイプライン。
RCT が実施できない状況で交絡因子を調整し、因果的解釈を可能にする。

## When to Use

- 観察データから因果効果（ATE / ATT）を推定したいとき
- 交絡因子の調整が必要なとき
- 傾向スコアによるマッチングや重み付けが必要なとき
- 自然実験データ（DID / RDD）を分析するとき
- DAG を描いて因果構造を明示化するとき

---

## Quick Start

## 1. DAG（有向非巡回グラフ）の定義

```python
import networkx as nx
import matplotlib.pyplot as plt

def define_causal_dag(edges, treatment, outcome, figsize=(10, 6)):
    """
    因果 DAG を定義し可視化する。

    Parameters:
        edges: list of (cause, effect) tuples
        treatment: 処置変数名
        outcome: アウトカム変数名

    Example:
        edges = [("Age", "Treatment"), ("Age", "Outcome"),
                 ("Treatment", "Outcome"), ("Gender", "Treatment")]
    """
    G = nx.DiGraph()
    G.add_edges_from(edges)

    # 色分け
    color_map = []
    for node in G.nodes():
        if node == treatment:
            color_map.append("#FF6B6B")
        elif node == outcome:
            color_map.append("#4ECDC4")
        else:
            color_map.append("#95E1D3")

    fig, ax = plt.subplots(figsize=figsize)
    pos = nx.spring_layout(G, k=2, seed=42)
    nx.draw(G, pos, ax=ax, with_labels=True, node_color=color_map,
            node_size=2000, font_size=11, font_weight="bold",
            edge_color="gray", arrows=True, arrowsize=20, width=2)
    ax.set_title("Causal DAG", fontweight="bold", fontsize=14)
    plt.tight_layout()
    plt.savefig("figures/causal_dag.png", dpi=300, bbox_inches="tight")
    plt.close()

    return G


def identify_confounders(dag, treatment, outcome):
    """
    Backdoor criterion に基づく共変量の同定。
    処置→アウトカムのバックドアパスを遮断するために調整すべき変数を返す。
    """
    # 簡易版: treatment の親ノードのうち、outcome にもパスがある変数
    parents_of_treatment = set(dag.predecessors(treatment))
    confounders = set()

    for parent in parents_of_treatment:
        if nx.has_path(dag, parent, outcome):
            confounders.add(parent)

    return confounders
```

## 2. 傾向スコアマッチング (PSM)

```python
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from scipy.spatial.distance import cdist

def propensity_score_matching(df, treatment_col, covariates, outcome_col,
                                caliper=0.2, n_matches=1):
    """
    傾向スコアマッチング。

    Steps:
        1. ロジスティック回帰で傾向スコア P(T=1|X) を推定
        2. 最近傍マッチング（キャリパー制約付き）
        3. マッチング後の共変量バランスチェック
        4. ATE / ATT の推定

    Parameters:
        caliper: マッチング許容距離（傾向スコアの標準偏差の倍率）
    """
    # Step 1: 傾向スコア推定
    X = df[covariates].values
    T = df[treatment_col].values

    lr = LogisticRegression(max_iter=1000, random_state=42)
    lr.fit(X, T)
    ps = lr.predict_proba(X)[:, 1]
    df = df.copy()
    df["propensity_score"] = ps

    # Step 2: マッチング
    treated_idx = df[df[treatment_col] == 1].index
    control_idx = df[df[treatment_col] == 0].index

    ps_treated = ps[treated_idx]
    ps_control = ps[control_idx]

    caliper_val = caliper * np.std(ps)
    matched_pairs = []

    for i, t_idx in enumerate(treated_idx):
        distances = np.abs(ps_treated[i] - ps_control)
        within_caliper = np.where(distances <= caliper_val)[0]
        if len(within_caliper) > 0:
            best = within_caliper[np.argmin(distances[within_caliper])]
            matched_pairs.append((t_idx, control_idx[best]))

    print(f"  Matched {len(matched_pairs)} / {len(treated_idx)} treated units")

    # Step 3: バランスチェック (SMD)
    matched_treated = df.loc[[p[0] for p in matched_pairs]]
    matched_control = df.loc[[p[1] for p in matched_pairs]]

    balance = []
    for cov in covariates:
        smd_before = _standardized_mean_diff(
            df[df[treatment_col]==1][cov], df[df[treatment_col]==0][cov])
        smd_after = _standardized_mean_diff(
            matched_treated[cov], matched_control[cov])
        balance.append({
            "covariate": cov,
            "SMD_before": smd_before,
            "SMD_after": smd_after,
            "balanced": abs(smd_after) < 0.1,
        })

    balance_df = pd.DataFrame(balance)

    # Step 4: ATT 推定
    att = matched_treated[outcome_col].mean() - matched_control[outcome_col].mean()

    return {
        "ATT": att,
        "n_matched": len(matched_pairs),
        "balance": balance_df,
        "propensity_scores": ps,
        "matched_pairs": matched_pairs,
    }


def _standardized_mean_diff(x1, x2):
    """Standardized Mean Difference (SMD) = |μ1 - μ2| / sqrt((s1² + s2²)/2)"""
    return abs(x1.mean() - x2.mean()) / np.sqrt((x1.var() + x2.var()) / 2 + 1e-10)
```

## 3. 逆確率重み付け (IPW / IPTW)

```python
def inverse_probability_weighting(df, treatment_col, covariates, outcome_col):
    """
    逆確率重み付け推定量 (IPTW - Inverse Probability of Treatment Weighting)。

    ATE = E[Y(1)] - E[Y(0)]
        = Σ (T·Y/PS) / Σ (T/PS) - Σ ((1-T)·Y/(1-PS)) / Σ ((1-T)/(1-PS))
    """
    X = df[covariates].values
    T = df[treatment_col].values
    Y = df[outcome_col].values

    # 傾向スコア推定
    lr = LogisticRegression(max_iter=1000, random_state=42)
    lr.fit(X, T)
    ps = lr.predict_proba(X)[:, 1]

    # 重み計算（安定化重み）
    w_treated = T / (ps + 1e-10)
    w_control = (1 - T) / (1 - ps + 1e-10)

    # ATE
    E_Y1 = np.sum(w_treated * Y) / np.sum(w_treated)
    E_Y0 = np.sum(w_control * Y) / np.sum(w_control)
    ate = E_Y1 - E_Y0

    # ATT
    att = np.mean(T * Y) - np.sum((1-T) * ps / (1-ps+1e-10) * Y) / np.sum((1-T) * ps / (1-ps+1e-10))

    return {
        "ATE": ate,
        "ATT": att,
        "E_Y1": E_Y1,
        "E_Y0": E_Y0,
        "propensity_scores": ps,
        "weights_treated": w_treated,
        "weights_control": w_control,
    }
```

## 4. 差分の差分法 (DID)

```python
import statsmodels.api as sm

def difference_in_differences(df, time_col, treatment_col, outcome_col,
                                covariates=None):
    """
    差分の差分法 (Difference-in-Differences)。

    Y = β0 + β1·Post + β2·Treat + β3·(Post × Treat) + ε
    β3 が因果効果（DID 推定量）
    """
    df = df.copy()
    df["interaction"] = df[time_col] * df[treatment_col]

    X_cols = [time_col, treatment_col, "interaction"]
    if covariates:
        X_cols += covariates

    X = sm.add_constant(df[X_cols])
    model = sm.OLS(df[outcome_col], X).fit()

    return {
        "DID_estimate": model.params["interaction"],
        "DID_se": model.bse["interaction"],
        "DID_pvalue": model.pvalues["interaction"],
        "DID_ci_95": model.conf_int().loc["interaction"].tolist(),
        "model_summary": model.summary2().tables[1],
    }
```

## 5. 回帰不連続デザイン (RDD)

```python
def regression_discontinuity(df, running_var, outcome_col, cutoff,
                               bandwidth=None, kernel="triangular"):
    """
    回帰不連続デザイン (Sharp RDD)。

    カットオフ前後での局所回帰による処置効果の推定。

    Parameters:
        running_var: 強制変数（running variable）列名
        cutoff: カットオフ値
        bandwidth: バンド幅（None = IK 最適バンド幅）
    """
    df = df.copy()
    df["centered"] = df[running_var] - cutoff
    df["treated"] = (df[running_var] >= cutoff).astype(int)

    if bandwidth is None:
        bandwidth = 1.5 * df["centered"].std()

    # バンド幅内のデータ
    in_band = df[df["centered"].abs() <= bandwidth]

    # カーネル重み
    if kernel == "triangular":
        weights = 1 - np.abs(in_band["centered"]) / bandwidth
    else:
        weights = np.ones(len(in_band))

    # 局所線形回帰
    X = sm.add_constant(in_band[["centered", "treated"]])
    X["interaction"] = in_band["centered"] * in_band["treated"]
    model = sm.WLS(in_band[outcome_col], X, weights=weights).fit()

    return {
        "RDD_estimate": model.params["treated"],
        "RDD_se": model.bse["treated"],
        "RDD_pvalue": model.pvalues["treated"],
        "RDD_ci_95": model.conf_int().loc["treated"].tolist(),
        "bandwidth": bandwidth,
        "n_in_bandwidth": len(in_band),
    }
```

## 6. 感度分析（Rosenbaum Bounds）

```python
def rosenbaum_sensitivity(matched_outcomes_treated, matched_outcomes_control,
                           gamma_range=None):
    """
    Rosenbaum 感度分析。
    隠れた交絡の影響度 Γ に対する因果推定の頑健性を評価。

    Γ = 1: 交絡なし仮定
    Γ > 1: 隠れた交絡がΓ倍のオッズ比を持つ場合
    """
    from scipy.stats import norm

    if gamma_range is None:
        gamma_range = np.arange(1.0, 3.1, 0.1)

    diffs = matched_outcomes_treated - matched_outcomes_control
    n = len(diffs)
    signs = (diffs > 0).astype(int)
    T_obs = np.sum(signs)

    results = []
    for gamma in gamma_range:
        p_upper = gamma / (1 + gamma)
        E_T = n * p_upper
        Var_T = n * p_upper * (1 - p_upper)
        z = (T_obs - E_T) / np.sqrt(Var_T + 1e-10)
        p_value = 1 - norm.cdf(z)
        results.append({"gamma": gamma, "z_statistic": z, "p_value": p_value})

    return pd.DataFrame(results)
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/causal_estimates.csv` | CSV |
| `results/covariate_balance.csv` | CSV |
| `results/sensitivity_analysis.csv` | CSV |
| `figures/causal_dag.png` | PNG |
| `figures/propensity_distribution.png` | PNG |
| `figures/rdd_plot.png` | PNG |

#### 依存パッケージ

```
statsmodels>=0.14
scikit-learn>=1.3
networkx>=3.0
```
