---
name: scientific-squidpy-advanced
description: |
  高度 Squidpy 空間解析スキル。空間自己相関・共起解析・空間
  近傍・リガンド受容体空間マッピング・ニッチ同定。
  K-Dense 連携: squidpy-advanced。
  ToolUniverse 連携: cellxgene。
tu_tools: []
kdense_ref: squidpy-advanced
tu_tools:
  - key: cellxgene
    name: CellxGene
    description: 空間トランスクリプトミクスデータ検索
---

# Scientific Squidpy Advanced

Squidpy の高度な空間統計解析機能を活用した空間トランスクリプ
トミクス解析パイプラインを提供する。

## When to Use

- 空間自己相関 (Moran's I, Geary's C) を計算するとき
- 細胞タイプの空間共起パターンを分析するとき
- 近傍エンリッチメント解析を実施するとき
- リガンド-受容体のペア空間マッピングを行うとき
- 空間ニッチ (組織微小環境) を同定するとき
- 空間グラフ中心性指標を計算するとき

---

## Quick Start

## 1. 空間自己相関解析

```python
import squidpy as sq
import scanpy as sc
import pandas as pd
import numpy as np


def squidpy_spatial_autocorrelation(adata, genes=None,
                                      mode="moran", n_perms=100,
                                      n_jobs=4):
    """
    Squidpy — 空間自己相関解析 (Moran's I / Geary's C)。

    Parameters:
        adata: AnnData — 空間トランスクリプトミクスデータ
        genes: list[str] — 対象遺伝子 (None=HVG 使用)
        mode: str — "moran" or "geary"
        n_perms: int — 並べ替え検定回数
        n_jobs: int — 並列ジョブ数
    """
    # 空間グラフ構築 (未構築の場合)
    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="generic",
                                  n_neighs=6)

    if genes is None:
        sc.pp.highly_variable_genes(adata, n_top_genes=200)
        genes = adata.var_names[adata.var["highly_variable"]].tolist()

    sq.gr.spatial_autocorr(adata, mode=mode, genes=genes,
                            n_perms=n_perms, n_jobs=n_jobs)

    result_key = f"{mode}I" if mode == "moran" else f"{mode}C"
    result = adata.uns[result_key].copy()
    sig = result[result["pval_norm"] < 0.05]

    print(f"Spatial autocorrelation ({mode}): "
          f"{len(sig)}/{len(result)} significant genes")
    return result.sort_values("pval_norm")
```

## 2. 空間共起解析

```python
def squidpy_co_occurrence(adata, cluster_key="cell_type",
                            interval=50, n_splits=None):
    """
    Squidpy — 細胞タイプ空間共起解析。

    Parameters:
        adata: AnnData — 空間トランスクリプトミクスデータ
        cluster_key: str — 細胞タイプアノテーションカラム
        interval: int — 距離インターバル数
        n_splits: int — 分割数 (大規模データ用, None=自動)
    """
    sq.gr.co_occurrence(adata, cluster_key=cluster_key,
                         interval=interval, n_splits=n_splits)

    co_occ = adata.uns[f"{cluster_key}_co_occurrence"]
    occ_scores = co_occ["occ"]
    interval_vals = co_occ["interval"]
    categories = adata.obs[cluster_key].cat.categories.tolist()

    # 最近距離での共起スコアマトリクス
    near_idx = 0
    scores = occ_scores[near_idx]
    co_df = pd.DataFrame(scores, index=categories, columns=categories)

    # 有意な共起ペア
    pairs = []
    for i, cat_i in enumerate(categories):
        for j, cat_j in enumerate(categories):
            if i < j:
                pairs.append({
                    "cell_type_a": cat_i,
                    "cell_type_b": cat_j,
                    "co_occurrence_score": scores[i, j],
                })
    pairs_df = pd.DataFrame(pairs).sort_values(
        "co_occurrence_score", ascending=False)

    print(f"Co-occurrence: {len(categories)} cell types, "
          f"top pair: {pairs_df.iloc[0]['cell_type_a']} — "
          f"{pairs_df.iloc[0]['cell_type_b']}")
    return co_df, pairs_df
```

## 3. 近傍エンリッチメント & ニッチ同定

```python
def squidpy_niche_identification(adata, cluster_key="cell_type",
                                    n_neighs=15, n_niches=5,
                                    seed=42):
    """
    Squidpy — 空間ニッチ同定。

    Parameters:
        adata: AnnData — 空間トランスクリプトミクスデータ
        cluster_key: str — 細胞タイプカラム
        n_neighs: int — 近傍数
        n_niches: int — ニッチクラスタ数
        seed: int — 乱数シード
    """
    # 近傍エンリッチメント
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
    nhood = adata.uns[f"{cluster_key}_nhood_enrichment"]
    zscore_matrix = nhood["zscore"]
    categories = adata.obs[cluster_key].cat.categories.tolist()
    nhood_df = pd.DataFrame(zscore_matrix, index=categories,
                              columns=categories)

    # 空間近傍グラフ
    sq.gr.spatial_neighbors(adata, coord_type="generic",
                              n_neighs=n_neighs)

    # 近傍組成プロファイル
    from sklearn.cluster import KMeans
    cell_types = adata.obs[cluster_key]
    type_dummies = pd.get_dummies(cell_types)

    # 各細胞の近傍組成
    conn = adata.obsp["spatial_connectivities"]
    nhood_composition = conn.dot(type_dummies.values)
    nhood_comp_df = pd.DataFrame(
        nhood_composition, columns=type_dummies.columns,
        index=adata.obs_names)

    # ニッチクラスタリング
    km = KMeans(n_clusters=n_niches, random_state=seed, n_init=10)
    adata.obs["spatial_niche"] = km.fit_predict(
        nhood_comp_df.values).astype(str)

    niche_summary = adata.obs.groupby("spatial_niche")[cluster_key]\
        .value_counts(normalize=True).unstack(fill_value=0)

    print(f"Niche identification: {n_niches} niches, "
          f"{len(categories)} cell types")
    return nhood_df, niche_summary
```

## 4. Squidpy 高度空間解析統合パイプライン

```python
def squidpy_advanced_pipeline(adata, cluster_key="cell_type",
                                output_dir="results"):
    """
    Squidpy 高度空間解析統合パイプライン。

    Parameters:
        adata: AnnData — 空間トランスクリプトミクスデータ
        cluster_key: str — 細胞タイプカラム
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 空間自己相関
    autocorr = squidpy_spatial_autocorrelation(adata)
    autocorr.to_csv(output_dir / "spatial_autocorrelation.csv")

    # 2) 空間共起
    co_df, pairs_df = squidpy_co_occurrence(adata,
        cluster_key=cluster_key)
    co_df.to_csv(output_dir / "co_occurrence_matrix.csv")
    pairs_df.to_csv(output_dir / "co_occurrence_pairs.csv",
                      index=False)

    # 3) ニッチ同定
    nhood_df, niche_summary = squidpy_niche_identification(
        adata, cluster_key=cluster_key)
    nhood_df.to_csv(output_dir / "nhood_enrichment.csv")
    niche_summary.to_csv(output_dir / "niche_composition.csv")

    # 4) 中心性スコア
    sq.gr.centrality_scores(adata, cluster_key=cluster_key)
    centrality = adata.uns[f"{cluster_key}_centrality_scores"]
    centrality_df = pd.DataFrame(centrality)
    centrality_df.to_csv(output_dir / "centrality_scores.csv")

    adata.write(output_dir / "adata_spatial.h5ad")

    print(f"Squidpy advanced pipeline: {output_dir}")
    return {
        "autocorrelation": autocorr,
        "co_occurrence": co_df,
        "niches": niche_summary,
        "centrality": centrality_df,
    }
```

---

## K-Dense 連携

| K-Dense Key | 参照内容 |
|-------------|---------|
| `squidpy-advanced` | 高度空間統計・ニッチ解析手法 |

## パイプライン統合

```
spatial-transcriptomics → squidpy-advanced → single-cell-genomics
  (Visium/MERFISH)       (空間統計)          (scRNA-seq)
       │                       │                   ↓
  image-analysis ─────────────┘            cell-communication
  (H&E/IF 画像)          │                (リガンド-受容体)
                          ↓
                    multi-omics
                    (空間マルチモーダル)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/spatial_autocorrelation.csv` | 空間自己相関 | → gene-expression |
| `results/co_occurrence_matrix.csv` | 共起マトリクス | → cell-communication |
| `results/niche_composition.csv` | ニッチ組成 | → single-cell-genomics |
| `results/centrality_scores.csv` | 中心性スコア | → spatial-transcriptomics |
| `results/adata_spatial.h5ad` | AnnData (全結果) | → multi-omics |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `cellxgene` | CellxGene | 空間トランスクリプトミクスデータ検索 |
