---
name: scientific-environmental-ecology
description: |
  環境科学・生態学解析スキル。種分布モデリング（SDM / MaxEnt）・
  生物多様性指標（α/β/γ 多様性）・群集構造解析（NMDS/CCA/RDA）・
  生態学的ニッチモデリング・保全優先順位評価・OBIS/GBIF データ統合パイプライン。
  ToolUniverse 連携: gbif。
tu_tools:
  - key: gbif
    name: GBIF
    description: 地球規模生物多様性情報ファシリティ
---

# Scientific Environmental Ecology

環境科学・生態学に特化した解析パイプラインを提供する。
種分布モデリング、生物多様性評価、群集構造解析、
保全優先順位付け、海洋/陸域の生態系データ統合を扱う。

## When to Use

- 種分布モデル（SDM）を構築して生息適地を推定するとき
- 群集の生物多様性指標を算出・比較するとき
- 群集構造の環境要因への応答を解析するとき（CCA / RDA）
- GBIF / OBIS から出現データを取得して空間解析を行うとき
- 保全優先区域の評価・ランキングを行うとき

---

## Quick Start

## 1. 種分布モデリング（SDM）

```python
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_auc_score

def species_distribution_model(occurrences, background, env_layers,
                                 method="maxent", n_folds=5):
    """
    種分布モデリング（SDM）パイプライン。

    method:
      - "maxent": MaxEnt — 最大エントロピーモデル（在データのみ可）
      - "rf": Random Forest — 在/不在データ
      - "gbm": Gradient Boosting — アンサンブル学習
      - "ensemble": 複数モデルの加重平均

    MaxEnt 原理:
      P(x) を環境変数 x の関数として推定。
      情報エントロピーを最大化する分布を選択:
        H(P) = -Σ P(x) log P(x) → maximize
      制約: E_P[fⱼ] = E_data[fⱼ] (特徴量の期待値一致)

    入力:
      - occurrences: 種の出現座標 (lon, lat)
      - background: 疑似不在点 (lon, lat)
      - env_layers: 環境変数ラスタ（Bio1-Bio19 等）
    """
    # 環境変数を出現/不在点で抽出
    X_pres = extract_env_values(occurrences, env_layers)
    X_bg = extract_env_values(background, env_layers)
    X = np.vstack([X_pres, X_bg])
    y = np.concatenate([np.ones(len(X_pres)), np.zeros(len(X_bg))])

    if method == "maxent":
        from elapid import MaxentModel
        model = MaxentModel()
        model.fit(X_pres, X_bg)
        pred = model.predict(env_layers)

    elif method == "rf":
        model = RandomForestClassifier(n_estimators=500, random_state=42)
        model.fit(X, y)
        auc_scores = cross_val_score(model, X, y, cv=n_folds, scoring="roc_auc")
        print(f"  RF AUC: {np.mean(auc_scores):.3f} ± {np.std(auc_scores):.3f}")
        pred = model.predict_proba(env_layers.reshape(-1, env_layers.shape[-1]))[:, 1]

    elif method == "gbm":
        model = GradientBoostingClassifier(n_estimators=300, max_depth=5,
                                            random_state=42)
        model.fit(X, y)
        auc_scores = cross_val_score(model, X, y, cv=n_folds, scoring="roc_auc")
        print(f"  GBM AUC: {np.mean(auc_scores):.3f} ± {np.std(auc_scores):.3f}")

    return model, pred


def extract_env_values(coords, env_layers):
    """座標から環境変数値を抽出する。"""
    import rasterio
    values = []
    for lon, lat in coords:
        row, col = env_layers.index(lon, lat)
        values.append(env_layers.read()[:, row, col])
    return np.array(values)
```

## 2. 生物多様性指標

```python
from scipy.stats import entropy

def biodiversity_indices(community_matrix, metadata=None):
    """
    群集ベースの生物多様性指標算出。

    α 多様性（サイト内）:
      - Species richness: S = 種数
      - Shannon: H' = -Σ pᵢ ln(pᵢ)
      - Simpson: D = 1 - Σ pᵢ²
      - Pielou's Evenness: J = H' / ln(S)
      - Chao1: S_est = S_obs + f₁²/(2·f₂)

    β 多様性（サイト間）:
      - Bray-Curtis dissimilarity
      - Jaccard distance
      - Sørensen index
      - Whittaker's β: γ/ᾱ - 1

    γ 多様性（景観全体）:
      - 全サイトの合計種数
    """
    results = []
    for idx, row in community_matrix.iterrows():
        counts = row[row > 0].values
        freqs = counts / counts.sum()
        S = len(counts)

        H = entropy(freqs)
        D_simpson = 1 - np.sum(freqs ** 2)
        J = H / np.log(S) if S > 1 else 0

        f1 = np.sum(counts == 1)
        f2 = max(np.sum(counts == 2), 1)
        chao1 = S + (f1 ** 2) / (2 * f2)

        results.append({
            "site": idx,
            "richness": S,
            "shannon": round(H, 4),
            "simpson": round(D_simpson, 4),
            "evenness": round(J, 4),
            "chao1": round(chao1, 1),
            "total_abundance": int(counts.sum()),
        })

    alpha_df = pd.DataFrame(results).set_index("site")

    # γ 多様性
    gamma = (community_matrix > 0).any(axis=0).sum()
    mean_alpha = alpha_df["richness"].mean()
    beta_whittaker = gamma / mean_alpha - 1

    summary = {
        "gamma_diversity": gamma,
        "mean_alpha": round(mean_alpha, 2),
        "beta_whittaker": round(beta_whittaker, 3),
    }

    print(f"  Biodiversity: γ={gamma}, ᾱ={mean_alpha:.1f}, β_w={beta_whittaker:.3f}")
    return alpha_df, summary
```

## 3. 群集構造解析（NMDS / CCA / RDA）

```python
def community_ordination(community_matrix, env_df=None, method="nmds",
                           n_dims=2, distance="bray"):
    """
    群集構造の序列化（Ordination）。

    method:
      - "nmds": Non-metric Multidimensional Scaling — ランクベース
      - "cca": Canonical Correspondence Analysis — 制約付き（単峰型応答）
      - "rda": Redundancy Analysis — 制約付き（線形応答）
      - "dca": Detrended Correspondence Analysis — 勾配長評価

    NMDS stress 基準:
      - < 0.05: Excellent
      - < 0.10: Good
      - < 0.20: Acceptable
      - > 0.20: Poor（次元数増加を検討）
    """
    from skbio.stats.ordination import pcoa
    from skbio.diversity import beta_diversity
    from scipy.spatial.distance import squareform

    if method == "nmds":
        from sklearn.manifold import MDS
        dm = beta_diversity(distance, community_matrix.values,
                             community_matrix.index)
        mds = MDS(n_components=n_dims, dissimilarity="precomputed",
                   metric=False, random_state=42, max_iter=500)
        coords = mds.fit_transform(squareform(dm.data))
        stress = mds.stress_
        print(f"  NMDS: stress={stress:.4f} ({n_dims}D)")
        return coords, stress

    elif method == "pcoa":
        dm = beta_diversity(distance, community_matrix.values,
                             community_matrix.index)
        result = pcoa(dm)
        return result.samples.values[:, :n_dims], result.proportion_explained[:n_dims]
```

## 4. 種の保全優先順位評価

```python
def conservation_priority(species_data, criteria_weights=None):
    """
    保全優先順位の多基準評価。

    IUCN レッドリスト基準:
      - CR: Critically Endangered
      - EN: Endangered
      - VU: Vulnerable
      - NT: Near Threatened

    評価基準:
      1. 絶滅リスク（IUCN カテゴリ）
      2. 系統的独自性（Evolutionary Distinctiveness）
      3. 生息地面積減少率
      4. Endemic 性（固有種かどうか）
      5. 生態系サービス寄与
    """
    if criteria_weights is None:
        criteria_weights = {
            "iucn_score": 0.30,
            "evolutionary_distinctiveness": 0.20,
            "habitat_loss_rate": 0.20,
            "endemism": 0.15,
            "ecosystem_service": 0.15,
        }

    iucn_mapping = {"CR": 5, "EN": 4, "VU": 3, "NT": 2, "LC": 1, "DD": 0}
    species_data["iucn_score"] = species_data["iucn_category"].map(iucn_mapping)

    # 正規化
    for col in criteria_weights:
        if col in species_data.columns:
            min_v = species_data[col].min()
            max_v = species_data[col].max()
            if max_v > min_v:
                species_data[f"{col}_norm"] = (species_data[col] - min_v) / (max_v - min_v)

    # Composite score
    species_data["priority_score"] = sum(
        w * species_data.get(f"{col}_norm", 0) for col, w in criteria_weights.items()
    )
    species_data = species_data.sort_values("priority_score", ascending=False)

    print(f"  Conservation: {len(species_data)} species ranked")
    return species_data
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/sdm_predictions.tif` | GeoTIFF |
| `results/biodiversity_indices.csv` | CSV |
| `results/ordination_scores.csv` | CSV |
| `results/conservation_priority.csv` | CSV |
| `figures/sdm_map.png` | PNG |
| `figures/nmds_plot.png` | PNG |
| `figures/diversity_comparison.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| OBIS | `OBIS_search_taxa` | 海洋生物分類検索 |
| OBIS | `OBIS_search_occurrences` | 海洋生物出現データ |
| GBIF | `GBIF_search_species` | 種名検索 |
| GBIF | `GBIF_search_occurrences` | 出現記録検索 |
| Paleobiology | `Paleobiology_get_fossils` | 化石記録データ |
| OLS | `ols_search_terms` | 生態学オントロジー検索 |
| PubMed | `PubMed_search_articles` | 生態学文献検索 |

### 参照スキル

| スキル | 連携内容 |
|---|---|
| [scientific-statistical-testing](../scientific-statistical-testing/SKILL.md) | 多様性有意差検定 |
| [scientific-pca-tsne](../scientific-pca-tsne/SKILL.md) | 次元削減・序列化 |
| [scientific-ml-classification](../scientific-ml-classification/SKILL.md) | SDM モデル（RF/GBM） |
| [scientific-image-analysis](../scientific-image-analysis/SKILL.md) | リモートセンシング画像解析 |
| [scientific-time-series](../scientific-time-series/SKILL.md) | 生態系時系列トレンド |

#### 依存パッケージ

- scikit-bio, rasterio, geopandas, elapid, shapely, pygbif
