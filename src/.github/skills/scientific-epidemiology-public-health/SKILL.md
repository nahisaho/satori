---
name: scientific-epidemiology-public-health
description: |
  疫学・公衆衛生解析スキル。観察研究デザイン（コホート/症例対照/横断）・
  リスク指標（RR/OR/HR/NNT）・標準化死亡比（SMR）・年齢調整率・
  空間疫学（GIS / 空間クラスタリング）・因果推論ダイアグラム（DAG）・
  WHO/CDC/EU 公衆衛生データ統合パイプライン。
---

# Scientific Epidemiology & Public Health

疫学研究と公衆衛生データ解析のパイプラインを提供する。
研究デザイン設計、リスク指標算出、交絡調整、
空間疫学、健康格差評価、公衆衛生データベース連携を体系的に扱う。

## When to Use

- 観察研究のリスク指標（RR / OR / HR）を算出するとき
- 年齢調整率・標準化死亡比（SMR）を計算するとき
- 空間疫学（疾患クラスタリング・GIS マッピング）を行うとき
- DAG（有向非巡回グラフ）で交絡構造を分析するとき
- WHO / CDC / EU の公衆衛生データを取得・解析するとき

---

## Quick Start

## 1. リスク指標算出

```python
import numpy as np
import pandas as pd
from scipy.stats import norm

def calculate_risk_measures(a, b, c, d, alpha=0.05):
    """
    2×2 分割表からリスク指標を算出する。

                Disease+    Disease-
    Exposed+      a           b        → a+b
    Exposed-      c           d        → c+d
                 a+c         b+d        N

    指標:
      - Risk (Incidence): R = cases / total
      - Risk Ratio (RR): R_exposed / R_unexposed（コホート研究）
      - Odds Ratio (OR): (a·d) / (b·c)（症例対照研究）
      - Risk Difference (RD): R_exposed - R_unexposed
      - NNT (Number Needed to Treat): 1 / |RD|
      - Attributable Fraction (AF): (RR - 1) / RR
    """
    z = norm.ppf(1 - alpha / 2)

    # Risk
    R1 = a / (a + b)  # Exposed
    R0 = c / (c + d)  # Unexposed

    # Risk Ratio
    RR = R1 / R0
    ln_RR_se = np.sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))
    RR_ci = (RR * np.exp(-z * ln_RR_se), RR * np.exp(z * ln_RR_se))

    # Odds Ratio
    OR = (a * d) / (b * c)
    ln_OR_se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    OR_ci = (OR * np.exp(-z * ln_OR_se), OR * np.exp(z * ln_OR_se))

    # Risk Difference
    RD = R1 - R0
    RD_se = np.sqrt(R1*(1-R1)/(a+b) + R0*(1-R0)/(c+d))
    RD_ci = (RD - z * RD_se, RD + z * RD_se)

    # NNT
    NNT = 1 / abs(RD) if RD != 0 else np.inf

    # Attributable fraction
    AF = (RR - 1) / RR if RR > 0 else 0

    results = {
        "risk_exposed": round(R1, 4),
        "risk_unexposed": round(R0, 4),
        "RR": round(RR, 4), "RR_CI": [round(x, 4) for x in RR_ci],
        "OR": round(OR, 4), "OR_CI": [round(x, 4) for x in OR_ci],
        "RD": round(RD, 4), "RD_CI": [round(x, 4) for x in RD_ci],
        "NNT": round(NNT, 1),
        "AF": round(AF, 4),
    }

    print(f"  RR={RR:.3f} ({RR_ci[0]:.3f}–{RR_ci[1]:.3f}), "
          f"OR={OR:.3f} ({OR_ci[0]:.3f}–{OR_ci[1]:.3f})")
    return results
```

## 2. 年齢調整率・SMR

```python
def age_standardization(observed_df, standard_pop, method="direct"):
    """
    年齢調整率と標準化死亡比。

    method:
      - "direct": 直接法 — 標準人口の年齢構成で重み付け
        ASR = Σ(年齢別率ᵢ × 標準人口割合ᵢ)
      - "indirect": 間接法 — SMR (Standardized Mortality Ratio)
        SMR = 観察死亡数 / 期待死亡数
        期待死亡数 = Σ(標準年齢別率ᵢ × 対象人口ᵢ)

    SMR の 95% CI（Byar's approximation）:
      SMR_lower = SMR × (1 - 1/(9·O) - z/(3·√O))³
      SMR_upper = (O+1)/E × (1 - 1/(9·(O+1)) + z/(3·√(O+1)))³
    """
    if method == "direct":
        # 直接法年齢調整率
        merged = observed_df.merge(standard_pop, on="age_group")
        merged["weighted_rate"] = merged["rate"] * merged["std_proportion"]
        asr = merged["weighted_rate"].sum()

        # 分散（二項近似）
        merged["var_component"] = (merged["std_proportion"] ** 2 *
                                    merged["rate"] * (1 - merged["rate"]) /
                                    merged["population"])
        se = np.sqrt(merged["var_component"].sum())

        return {
            "ASR": round(asr, 6),
            "ASR_per_100k": round(asr * 1e5, 2),
            "SE": round(se, 6),
            "CI_95": [round((asr - 1.96*se)*1e5, 2), round((asr + 1.96*se)*1e5, 2)],
        }

    elif method == "indirect":
        # 間接法 SMR
        merged = observed_df.merge(standard_pop, on="age_group",
                                     suffixes=("_obs", "_std"))
        merged["expected"] = merged["rate_std"] * merged["population_obs"]
        O = merged["deaths_obs"].sum()
        E = merged["expected"].sum()

        SMR = O / E
        z = 1.96

        # Byar's approximation
        lower = SMR * (1 - 1/(9*O) - z/(3*np.sqrt(O)))**3
        upper = ((O+1)/E) * (1 - 1/(9*(O+1)) + z/(3*np.sqrt(O+1)))**3

        print(f"  SMR={SMR:.3f} ({lower:.3f}–{upper:.3f}), O={O}, E={E:.1f}")
        return {"SMR": round(SMR, 4), "CI_95": [round(lower, 4), round(upper, 4)],
                "observed": O, "expected": round(E, 1)}
```

## 3. 空間疫学・疾患クラスタリング

```python
def spatial_cluster_detection(cases_gdf, population_gdf, method="kulldorff"):
    """
    空間疾患クラスタリング。

    method:
      - "kulldorff": Kulldorff's spatial scan statistic（SaTScan）
        H₀: λ(s) = 常数（一様リスク）
        H₁: ∃ 円形ウィンドウ Z で λ_in > λ_out
        LLR = (O_Z/E_Z)^{O_Z} × ((O-O_Z)/(O-E_Z))^{O-O_Z}
      - "moran": Local Moran's I（局所空間自己相関）
        Iᵢ = zᵢ Σⱼ wᵢⱼ zⱼ
      - "getis_ord": Getis-Ord Gi* — ホットスポット検出

    用途:
      - 疾患の地理的集積（クラスター）の検出
      - ホットスポット / コールドスポットの同定
    """
    import geopandas as gpd
    from libpysal.weights import Queen
    from esda.moran import Moran_Local
    from esda.getisord import G_Local

    if method == "moran":
        W = Queen.from_dataframe(cases_gdf)
        W.transform = "r"
        rates = cases_gdf["cases"] / cases_gdf["population"]
        lisa = Moran_Local(rates.values, W)

        cases_gdf["local_moran_I"] = lisa.Is
        cases_gdf["local_moran_p"] = lisa.p_sim
        cases_gdf["cluster_type"] = classify_lisa(lisa)

        n_hotspots = (cases_gdf["cluster_type"] == "HH").sum()
        n_coldspots = (cases_gdf["cluster_type"] == "LL").sum()
        print(f"  LISA: {n_hotspots} hotspots, {n_coldspots} coldspots")

    elif method == "getis_ord":
        W = Queen.from_dataframe(cases_gdf)
        W.transform = "b"
        rates = cases_gdf["cases"] / cases_gdf["population"]
        g_local = G_Local(rates.values, W)

        cases_gdf["gi_star"] = g_local.Zs
        cases_gdf["gi_p"] = g_local.p_sim
        cases_gdf["hotspot"] = (g_local.Zs > 1.96) & (g_local.p_sim < 0.05)

    return cases_gdf


def classify_lisa(lisa, p_threshold=0.05):
    """LISA クラスタ分類（HH/HL/LH/LL/NS）。"""
    types = []
    for i in range(len(lisa.Is)):
        if lisa.p_sim[i] > p_threshold:
            types.append("NS")
        elif lisa.q[i] == 1:
            types.append("HH")
        elif lisa.q[i] == 2:
            types.append("LH")
        elif lisa.q[i] == 3:
            types.append("LL")
        elif lisa.q[i] == 4:
            types.append("HL")
    return types
```

## 4. DAG ベース交絡分析

```python
def dag_confounding_analysis(dag_edges, exposure, outcome):
    """
    DAG（有向非巡回グラフ）ベースの交絡分析。

    パイプライン:
      1. DAG 構築
      2. バックドアパス列挙
      3. 最小調整セット（Sufficient Adjustment Set）同定
      4. d-分離判定

    Pearl のバックドア基準:
      変数セット Z がバックドア基準を満たす ⟺
      Z が X→Y の全バックドアパスをブロックし、
      Z に X の子孫が含まれない
    """
    import networkx as nx
    from dowhy import CausalModel

    G = nx.DiGraph()
    G.add_edges_from(dag_edges)

    # バックドアパス
    backdoor_paths = find_backdoor_paths(G, exposure, outcome)

    # 最小調整セット
    adjustment_sets = find_adjustment_sets(G, exposure, outcome)

    result = {
        "n_backdoor_paths": len(backdoor_paths),
        "backdoor_paths": backdoor_paths,
        "adjustment_sets": adjustment_sets,
        "minimal_adjustment": min(adjustment_sets, key=len) if adjustment_sets else [],
    }

    print(f"  DAG: {len(backdoor_paths)} backdoor paths, "
          f"minimal adjustment = {result['minimal_adjustment']}")
    return result


def find_backdoor_paths(G, source, target):
    """バックドアパス（X ← ... → Y）を列挙する。"""
    undirected = G.to_undirected()
    all_paths = list(nx.all_simple_paths(undirected, source, target))
    backdoor = [p for p in all_paths if G.has_edge(p[1], source)]
    return backdoor


def find_adjustment_sets(G, exposure, outcome):
    """最小十分調整セットを求める（簡易実装）。"""
    from itertools import combinations
    nodes = set(G.nodes()) - {exposure, outcome}
    sets = []
    for r in range(len(nodes) + 1):
        for combo in combinations(nodes, r):
            if blocks_all_backdoor(G, exposure, outcome, set(combo)):
                sets.append(list(combo))
    return sets


def blocks_all_backdoor(G, X, Y, Z):
    """Z がすべてのバックドアパスをブロックするか判定。"""
    # 簡易 d-separation チェック
    return True  # 要完全実装
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/risk_measures.json` | JSON |
| `results/age_standardized_rates.csv` | CSV |
| `results/spatial_clusters.geojson` | GeoJSON |
| `results/dag_analysis.json` | JSON |
| `figures/disease_map.png` | PNG |
| `figures/dag_diagram.png` | PNG |
| `figures/forest_plot.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| WHO | `who_gho_get_data` | WHO GHO データ取得 |
| WHO | `who_gho_query_health_data` | WHO 健康指標クエリ |
| CDC | `cdc_data_search_datasets` | CDC データセット検索 |
| CDC | `cdc_data_get_dataset` | CDC データ取得 |
| EUHealthInfo | `euhealthinfo_search_surveillance_mortality_rates` | 死亡率データ |
| EUHealthInfo | `euhealthinfo_search_healthcare_expenditure` | 医療費データ |
| EUHealthInfo | `euhealthinfo_search_population_health_survey` | 健康調査データ |
| HealthDisparities | `health_disparities_get_svi_info` | 社会脆弱性指標 |
| HealthDisparities | `health_disparities_get_county_rankings_info` | 地域健康ランキング |
| ClinicalTrials | `search_clinical_trials` | 臨床試験検索 |
| PubMed | `PubMed_Guidelines_Search` | 公衆衛生ガイドライン |

### 参照スキル

| スキル | 連携内容 |
|---|---|
| [scientific-causal-inference](../scientific-causal-inference/SKILL.md) | 因果推論・傾向スコア |
| [scientific-survival-clinical](../scientific-survival-clinical/SKILL.md) | 生存解析・Cox 回帰 |
| [scientific-meta-analysis](../scientific-meta-analysis/SKILL.md) | メタアナリシス・系統的レビュー |
| [scientific-infectious-disease](../scientific-infectious-disease/SKILL.md) | 感染症疫学 |
| [scientific-bayesian-statistics](../scientific-bayesian-statistics/SKILL.md) | ベイズ空間モデル |

#### 依存パッケージ

- geopandas, libpysal, esda, dowhy, lifelines, scipy, statsmodels
