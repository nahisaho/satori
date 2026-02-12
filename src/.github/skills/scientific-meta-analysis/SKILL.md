---
name: scientific-meta-analysis
description: |
  メタ解析スキル。固定効果・ランダム効果モデル（DerSimonian-Laird）、Forest プロット、
  異質性評価（I²/Q 検定/τ²）、出版バイアス検出（Funnel プロット/Egger/Begg 検定）、
  サブグループ解析、メタ回帰、累積メタ解析のテンプレートを提供。
---

# Scientific Meta-Analysis

複数の独立した研究結果を統合し、全体的なエビデンスを定量化するためのスキル。
効果量（SMD / OR / RR / MD）の統合、異質性評価、出版バイアス検出のパイプラインを
提供する。

## When to Use

- 複数の研究・実験結果を統合的に評価するとき
- 効果量（Hedges' g / Cohen's d / SMD）を算出するとき
- Forest プロット / Funnel プロットを描画するとき
- 研究間の異質性を定量化するとき（I² / Q / τ²）
- 出版バイアスの有無を検定するとき

---

## Quick Start

## 1. 効果量の算出

```python
import numpy as np
import pandas as pd
from scipy.stats import norm

def compute_effect_sizes(studies_df, effect_type="SMD"):
    """
    各研究の効果量と分散を算出する。

    effect_type:
        "SMD"  — Standardized Mean Difference (Hedges' g)
        "MD"   — Mean Difference (同スケール)
        "OR"   — Odds Ratio (log 変換)
        "RR"   — Risk Ratio (log 変換)

    Input columns (SMD/MD):
        mean1, sd1, n1, mean2, sd2, n2

    Input columns (OR/RR):
        events1, total1, events2, total2
    """
    df = studies_df.copy()

    if effect_type == "SMD":
        # Cohen's d → Hedges' g (小標本補正)
        pooled_sd = np.sqrt(
            ((df["n1"]-1)*df["sd1"]**2 + (df["n2"]-1)*df["sd2"]**2) /
            (df["n1"] + df["n2"] - 2)
        )
        d = (df["mean1"] - df["mean2"]) / pooled_sd
        # Hedges' correction factor
        J = 1 - 3 / (4*(df["n1"]+df["n2"]-2) - 1)
        df["effect_size"] = d * J
        df["variance"] = (df["n1"]+df["n2"])/(df["n1"]*df["n2"]) + \
                          df["effect_size"]**2 / (2*(df["n1"]+df["n2"]))

    elif effect_type == "MD":
        df["effect_size"] = df["mean1"] - df["mean2"]
        df["variance"] = df["sd1"]**2/df["n1"] + df["sd2"]**2/df["n2"]

    elif effect_type == "OR":
        a = df["events1"]; b = df["total1"] - df["events1"]
        c = df["events2"]; d_val = df["total2"] - df["events2"]
        df["effect_size"] = np.log((a * d_val) / (b * c + 1e-10) + 1e-10)
        df["variance"] = 1/a + 1/b + 1/c + 1/d_val

    elif effect_type == "RR":
        p1 = df["events1"] / df["total1"]
        p2 = df["events2"] / df["total2"]
        df["effect_size"] = np.log(p1 / (p2 + 1e-10) + 1e-10)
        df["variance"] = (1-p1)/(df["events1"]+1e-10) + \
                          (1-p2)/(df["events2"]+1e-10)

    df["se"] = np.sqrt(df["variance"])
    df["ci_lower"] = df["effect_size"] - 1.96 * df["se"]
    df["ci_upper"] = df["effect_size"] + 1.96 * df["se"]
    df["weight"] = 1 / df["variance"]

    return df
```

## 2. 固定効果 / ランダム効果モデル

```python
def meta_analysis(studies_df, model="random"):
    """
    メタ解析統合。

    model:
        "fixed"   — 固定効果モデル (Inverse-Variance weighted)
        "random"  — ランダム効果モデル (DerSimonian-Laird)

    Input: DataFrame with columns: study, effect_size, variance
    """
    es = studies_df["effect_size"].values
    var = studies_df["variance"].values
    w = 1 / var
    k = len(es)

    # 固定効果
    theta_fixed = np.sum(w * es) / np.sum(w)
    se_fixed = 1 / np.sqrt(np.sum(w))

    # 異質性
    Q = np.sum(w * (es - theta_fixed)**2)
    df_Q = k - 1
    p_Q = 1 - __import__("scipy").stats.chi2.cdf(Q, df_Q)
    I2 = max(0, (Q - df_Q) / Q * 100) if Q > 0 else 0

    if model == "random":
        # DerSimonian-Laird τ² 推定
        C = np.sum(w) - np.sum(w**2) / np.sum(w)
        tau2 = max(0, (Q - df_Q) / C)

        # ランダム効果重み
        w_re = 1 / (var + tau2)
        theta_random = np.sum(w_re * es) / np.sum(w_re)
        se_random = 1 / np.sqrt(np.sum(w_re))

        summary_effect = theta_random
        summary_se = se_random
    else:
        tau2 = 0
        summary_effect = theta_fixed
        summary_se = se_fixed

    z = summary_effect / summary_se
    p_val = 2 * (1 - norm.cdf(abs(z)))

    return {
        "model": model,
        "summary_effect": summary_effect,
        "se": summary_se,
        "ci_lower": summary_effect - 1.96 * summary_se,
        "ci_upper": summary_effect + 1.96 * summary_se,
        "z_value": z,
        "p_value": p_val,
        "Q_statistic": Q,
        "Q_p_value": p_Q,
        "I_squared": I2,
        "tau_squared": tau2,
        "k_studies": k,
    }
```

## 3. Forest プロット

```python
import matplotlib.pyplot as plt

def forest_plot(studies_df, meta_result, effect_label="SMD",
                figsize=(10, None)):
    """
    Forest プロットを描画する。

    studies_df: study, effect_size, ci_lower, ci_upper, weight
    meta_result: meta_analysis() の出力
    """
    k = len(studies_df)
    if figsize[1] is None:
        figsize = (figsize[0], max(4, k * 0.4 + 2))

    fig, ax = plt.subplots(figsize=figsize)

    y_positions = range(k, 0, -1)

    # 個別研究
    for i, (_, row) in enumerate(studies_df.iterrows()):
        y = list(y_positions)[i]
        ax.plot([row["ci_lower"], row["ci_upper"]], [y, y],
               "b-", linewidth=1.5)
        size = row.get("weight", 1) / studies_df["weight"].max() * 200 + 20
        ax.plot(row["effect_size"], y, "bs", markersize=np.sqrt(size),
               markerfacecolor="steelblue")
        ax.text(-0.05, y, row.get("study", f"Study {i+1}"),
               ha="right", va="center", fontsize=9,
               transform=ax.get_yaxis_transform())

    # サマリーダイヤモンド
    y_summary = 0
    diamond_x = [meta_result["ci_lower"], meta_result["summary_effect"],
                 meta_result["ci_upper"], meta_result["summary_effect"]]
    diamond_y = [y_summary, y_summary + 0.3, y_summary, y_summary - 0.3]
    ax.fill(diamond_x, diamond_y, color="red", alpha=0.7)
    ax.text(-0.05, y_summary, "Summary",
           ha="right", va="center", fontsize=9, fontweight="bold",
           transform=ax.get_yaxis_transform())

    # 参照線
    ax.axvline(0, color="black", linestyle="-", linewidth=0.5)

    ax.set_xlabel(effect_label)
    ax.set_yticks([])
    ax.set_title(f"Forest Plot (I²={meta_result['I_squared']:.1f}%, "
                f"p={meta_result['p_value']:.4f})", fontweight="bold")

    # 右側に数値
    for i, (_, row) in enumerate(studies_df.iterrows()):
        y = list(y_positions)[i]
        ax.text(1.02, y, f"{row['effect_size']:.2f} [{row['ci_lower']:.2f}, "
               f"{row['ci_upper']:.2f}]",
               ha="left", va="center", fontsize=8,
               transform=ax.get_yaxis_transform())

    plt.tight_layout()
    plt.savefig("figures/forest_plot.png", dpi=300, bbox_inches="tight")
    plt.close()
```

## 4. 出版バイアス検出

### 4.1 Funnel プロット

```python
def funnel_plot(studies_df, meta_result, figsize=(8, 6)):
    """
    Funnel プロットを描画する。
    非対称なら出版バイアスの存在を示唆。
    """
    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(studies_df["effect_size"], studies_df["se"],
              c="steelblue", s=50, edgecolors="black", zorder=5)

    # 参照線
    ax.axvline(meta_result["summary_effect"], color="red", linestyle="--")

    # 95% 擬似信頼区間
    se_range = np.linspace(0, studies_df["se"].max() * 1.1, 100)
    ax.plot(meta_result["summary_effect"] - 1.96 * se_range, se_range,
           "gray", linestyle="--", alpha=0.5)
    ax.plot(meta_result["summary_effect"] + 1.96 * se_range, se_range,
           "gray", linestyle="--", alpha=0.5)

    ax.set_xlabel("Effect Size")
    ax.set_ylabel("Standard Error")
    ax.set_title("Funnel Plot", fontweight="bold")
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig("figures/funnel_plot.png", dpi=300, bbox_inches="tight")
    plt.close()
```

### 4.2 Egger 検定

```python
import statsmodels.api as sm

def egger_test(studies_df):
    """
    Egger 回帰検定 — Funnel プロットの非対称性を統計的に検定。

    y = effect_size / se
    x = 1 / se
    切片 ≠ 0 → 出版バイアスあり
    """
    precision = 1 / studies_df["se"]
    z_score = studies_df["effect_size"] / studies_df["se"]

    X = sm.add_constant(precision)
    model = sm.OLS(z_score, X).fit()

    return {
        "intercept": model.params["const"],
        "intercept_se": model.bse["const"],
        "intercept_p": model.pvalues["const"],
        "publication_bias": model.pvalues["const"] < 0.05,
    }
```

## 5. サブグループ解析

```python
def subgroup_analysis(studies_df, subgroup_col, model="random"):
    """サブグループごとにメタ解析を行い、グループ間差を検定する。"""
    subgroups = studies_df[subgroup_col].unique()
    results = []

    for sg in subgroups:
        subset = studies_df[studies_df[subgroup_col] == sg]
        if len(subset) >= 2:
            ma = meta_analysis(subset, model=model)
            ma["subgroup"] = sg
            ma["k"] = len(subset)
            results.append(ma)

    results_df = pd.DataFrame(results)

    # グループ間検定 (Q_between)
    overall = meta_analysis(studies_df, model=model)
    Q_within = sum(r["Q_statistic"] for r in results)
    Q_between = overall["Q_statistic"] - Q_within
    df_between = len(results) - 1
    p_between = 1 - __import__("scipy").stats.chi2.cdf(Q_between, df_between)

    return {
        "subgroup_results": results_df,
        "Q_between": Q_between,
        "df_between": df_between,
        "p_between": p_between,
    }
```

## 6. 累積メタ解析

```python
def cumulative_meta_analysis(studies_df, sort_by="year", model="random"):
    """
    研究を順に追加しながらメタ解析を実行する。
    エビデンスの蓄積過程を可視化。
    """
    sorted_df = studies_df.sort_values(sort_by)
    cumulative_results = []

    for i in range(2, len(sorted_df) + 1):
        subset = sorted_df.iloc[:i]
        ma = meta_analysis(subset, model=model)
        ma["n_studies"] = i
        ma["last_added"] = sorted_df.iloc[i-1].get("study", f"Study {i}")
        cumulative_results.append(ma)

    return pd.DataFrame(cumulative_results)
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/meta_analysis_summary.csv` | CSV |
| `results/effect_sizes.csv` | CSV |
| `results/publication_bias_tests.csv` | CSV |
| `results/subgroup_analysis.csv` | CSV |
| `figures/forest_plot.png` | PNG |
| `figures/funnel_plot.png` | PNG |
| `figures/cumulative_meta.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| PubMed | `PubMed_search_articles` | メタアナリシス対象文献検索 |
| PubMed | `PubMed_get_article` | 論文メタデータ取得 |
| EuropePMC | `EuropePMC_search_articles` | ヨーロッパ文献検索 |
| Crossref | `Crossref_search_works` | 出版情報検索 |

#### 依存パッケージ

```
scipy>=1.10
statsmodels>=0.14
```
