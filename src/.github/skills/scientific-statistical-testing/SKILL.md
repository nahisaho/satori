---
name: scientific-statistical-testing
description: |
  統計検定・多重比較・エンリッチメント解析のスキル。t検定、カイ二乗検定、ANOVA、
  Bonferroni/BH 補正、Fisher 正確検定、ベイズ推論を行う際に使用。
  Scientific Skills Exp-03, 04, 06, 07 で確立したパターン。
---

# Scientific Statistical Testing & Enrichment Analysis

仮説検定、多重比較補正、エンリッチメント解析のための統計パイプラインスキル。
頻度論的検定とベイズ推論の両方のアプローチを提供する。

## When to Use

- 2 群間の有意差を検定したいとき（t 検定、Mann-Whitney U）
- 多群間の比較（ANOVA、Kruskal-Wallis）
- 多重比較の補正（Bonferroni、Benjamini-Hochberg）
- パスウェイエンリッチメント解析（Fisher 正確検定）
- ベイズ推論（Beta-Binomial 共役モデル）

## Quick Start

## 標準パイプライン

### 1. 2 群間検定

```python
from scipy import stats
import numpy as np
import pandas as pd

def two_group_test(group1, group2, test="auto", alternative="two-sided"):
    """
    2 群間の検定を実行する。
    test='auto' の場合、正規性検定に基づいて t 検定 or Mann-Whitney U を選択。
    """
    # 正規性検定（Shapiro-Wilk）
    if test == "auto":
        _, p1 = stats.shapiro(group1) if len(group1) <= 5000 else (0, 0.05)
        _, p2 = stats.shapiro(group2) if len(group2) <= 5000 else (0, 0.05)
        test = "ttest" if (p1 > 0.05 and p2 > 0.05) else "mannwhitney"

    if test == "ttest":
        stat, pval = stats.ttest_ind(group1, group2, alternative=alternative)
        test_name = "Welch's t-test"
    elif test == "mannwhitney":
        stat, pval = stats.mannwhitneyu(group1, group2, alternative=alternative)
        test_name = "Mann-Whitney U"
    else:
        raise ValueError(f"Unknown test: {test}")

    # 効果量（Cohen's d）
    pooled_std = np.sqrt((np.var(group1, ddof=1) + np.var(group2, ddof=1)) / 2)
    cohens_d = (np.mean(group1) - np.mean(group2)) / pooled_std if pooled_std > 0 else 0

    return {
        "test": test_name,
        "statistic": stat,
        "p_value": pval,
        "cohens_d": cohens_d,
        "effect_size": ("large" if abs(cohens_d) > 0.8 else
                       "medium" if abs(cohens_d) > 0.5 else "small"),
    }
```

### 2. 多重比較補正

```python
from statsmodels.stats.multitest import multipletests

def multiple_testing_correction(p_values, method="fdr_bh", alpha=0.05):
    """
    多重比較補正を適用する。
    method: 'bonferroni', 'fdr_bh' (Benjamini-Hochberg), 'holm'
    """
    reject, p_corrected, _, _ = multipletests(p_values, alpha=alpha, method=method)
    return reject, p_corrected
```

### 3. ANOVA / Kruskal-Wallis

```python
def multi_group_test(groups, test="auto"):
    """
    多群間の比較を実行する。
    groups: [array1, array2, ...] のリスト
    """
    # 正規性チェック
    normal = all(stats.shapiro(g)[1] > 0.05 for g in groups if len(g) <= 5000)

    if test == "auto":
        test = "anova" if normal else "kruskal"

    if test == "anova":
        stat, pval = stats.f_oneway(*groups)
        test_name = "One-way ANOVA"
    elif test == "kruskal":
        stat, pval = stats.kruskal(*groups)
        test_name = "Kruskal-Wallis"

    return {"test": test_name, "statistic": stat, "p_value": pval}
```

### 4. Fisher 正確検定パスウェイエンリッチメント（Exp-04, 07）

```python
def pathway_enrichment(deg_list, pathway_dict, background_size,
                       method="fisher", correction="fdr_bh"):
    """
    Fisher 正確検定によるパスウェイエンリッチメント解析。
    deg_list: 差次発現遺伝子リスト
    pathway_dict: {pathway_name: [gene1, gene2, ...]} の辞書
    background_size: バックグラウンド遺伝子数
    """
    results = []
    deg_set = set(deg_list)

    for pathway, genes in pathway_dict.items():
        gene_set = set(genes)
        overlap = deg_set & gene_set

        # 2×2 分割表
        a = len(overlap)                           # DEG ∩ Pathway
        b = len(deg_set) - a                       # DEG ∩ ~Pathway
        c = len(gene_set) - a                      # ~DEG ∩ Pathway
        d = background_size - a - b - c            # ~DEG ∩ ~Pathway

        _, pval = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        fold_enrichment = (a / len(deg_set)) / (len(gene_set) / background_size) \
                         if len(gene_set) > 0 and len(deg_set) > 0 else 0

        results.append({
            "Pathway": pathway,
            "Overlap": a,
            "Pathway_Size": len(gene_set),
            "Fold_Enrichment": fold_enrichment,
            "p_value": pval,
            "Genes": ", ".join(sorted(overlap)),
        })

    results_df = pd.DataFrame(results)

    # 多重検定補正
    if len(results_df) > 0:
        reject, p_adj = multiple_testing_correction(
            results_df["p_value"].values, method=correction
        )
        results_df["p_adjusted"] = p_adj
        results_df["Significant"] = reject

    results_df = results_df.sort_values("p_value")
    results_df.to_csv("results/pathway_enrichment.csv", index=False)
    return results_df
```

### 5. ベイズ推論（Beta-Binomial, Exp-06 パターン）

```python
def bayesian_beta_binomial(successes, trials, prior_alpha=1, prior_beta=1):
    """
    Beta-Binomial 共役モデルによるベイズ推論。
    事前分布: Beta(alpha, beta), デフォルトは一様事前分布。
    """
    post_alpha = prior_alpha + successes
    post_beta = prior_beta + (trials - successes)

    from scipy.stats import beta
    posterior = beta(post_alpha, post_beta)

    return {
        "posterior_mean": posterior.mean(),
        "posterior_std": posterior.std(),
        "95%_CI": (posterior.ppf(0.025), posterior.ppf(0.975)),
        "MAP": (post_alpha - 1) / (post_alpha + post_beta - 2)
               if post_alpha > 1 and post_beta > 1 else posterior.mean(),
        "posterior_alpha": post_alpha,
        "posterior_beta": post_beta,
    }
```

### 6. 生存解析（Kaplan-Meier + Cox PH, Exp-03/06）

```python
def survival_analysis(df, time_col, event_col, group_col):
    """
    Kaplan-Meier 生存曲線と Log-rank 検定を実行する。
    lifelines ライブラリが必要。
    """
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test
    import matplotlib.pyplot as plt

    groups = df[group_col].unique()
    fig, ax = plt.subplots(figsize=(8, 6))
    kmf = KaplanMeierFitter()

    for group in sorted(groups):
        mask = df[group_col] == group
        kmf.fit(df.loc[mask, time_col], event_observed=df.loc[mask, event_col],
                label=str(group))
        kmf.plot_survival_function(ax=ax)

    # Log-rank 検定（2 群の場合）
    if len(groups) == 2:
        g1 = df[df[group_col] == groups[0]]
        g2 = df[df[group_col] == groups[1]]
        lr = logrank_test(g1[time_col], g2[time_col],
                         event_observed_A=g1[event_col],
                         event_observed_B=g2[event_col])
        ax.text(0.7, 0.9, f"Log-rank p={lr.p_value:.4f}",
               transform=ax.transAxes, fontsize=10)

    ax.set_xlabel("Time")
    ax.set_ylabel("Survival Probability")
    ax.set_title("Kaplan-Meier Survival Curves", fontweight="bold")
    plt.tight_layout()
    plt.savefig("figures/kaplan_meier.png", dpi=300, bbox_inches="tight")
    plt.close()

    return lr.p_value if len(groups) == 2 else None
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/pathway_enrichment.csv` | CSV |
| `results/statistical_tests.csv` | CSV |
| `figures/kaplan_meier.png` | PNG |
| `figures/enrichment_dotplot.png` | PNG |

#### 参照実験

- **Exp-03**: Mann-Whitney U + Volcano Plot + 生存解析
- **Exp-04**: Fisher パスウェイエンリッチメント + Louvain コミュニティ
- **Exp-06**: 頻度論 + ベイズ推論 + 検出力分析
- **Exp-07**: Welch t 検定 + BH 補正 + PLS-DA VIP
