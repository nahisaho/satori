---
name: scientific-single-cell-genomics
description: |
  シングルセルゲノミクス解析スキル。scRNA-seq データの品質管理・正規化・
  次元削減（PCA/UMAP）・クラスタリング（Leiden）・差次発現遺伝子（DEG）同定・
  セルタイプアノテーション・RNA velocity・細胞間コミュニケーション推定パイプライン。
  Scanpy / AnnData フレームワークに準拠。
---

# Scientific Single-Cell Genomics

scRNA-seq / snRNA-seq データを対象に、QC → 正規化 → 高変動遺伝子選択 →
次元削減 → クラスタリング → DEG → セルタイプアノテーションの標準パイプラインを提供する。
CELLxGENE Census・HCA データポータルとの連携も組み込む。

## When to Use

- scRNA-seq / snRNA-seq データの解析パイプラインが必要なとき
- セルタイプのクラスタリングとアノテーションを行うとき
- クラスタ間の差次発現遺伝子を同定するとき
- RNA velocity による細胞分化軌跡を解析するとき
- CellChat / CellPhoneDB による細胞間コミュニケーション推定を行うとき

---

## Quick Start

## 1. QC・前処理パイプライン

```python
import scanpy as sc
import numpy as np
import pandas as pd

def sc_qc_preprocessing(adata, min_genes=200, max_genes=5000,
                         max_pct_mito=20, min_cells=3):
    """
    scRNA-seq 標準 QC パイプライン。

    QC メトリクス:
      - n_genes_by_counts: 細胞あたり検出遺伝子数
      - total_counts: 細胞あたり総 UMI カウント
      - pct_counts_mt: ミトコンドリア遺伝子比率（%）

    フィルタリング基準:
      - min_genes ≤ n_genes ≤ max_genes
      - pct_mito ≤ max_pct_mito
      - 遺伝子は min_cells 以上の細胞で発現
    """
    # ミトコンドリア遺伝子のアノテーション
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    n_before = adata.n_obs
    # 細胞フィルタ
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["n_genes_by_counts"] <= max_genes].copy()
    adata = adata[adata.obs["pct_counts_mt"] <= max_pct_mito].copy()
    # 遺伝子フィルタ
    sc.pp.filter_genes(adata, min_cells=min_cells)

    n_after = adata.n_obs
    print(f"  QC: {n_before} → {n_after} cells ({n_before - n_after} removed)")

    return adata


def sc_normalize(adata, target_sum=1e4, n_top_genes=2000):
    """
    正規化・HVG 選択パイプライン。

    手順:
      1. Library size normalization (target_sum)
      2. Log1p 変換
      3. Highly Variable Gene (HVG) 選択 (Seurat v3 法)
    """
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers["log_normalized"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat_v3",
                                 layer="counts")
    print(f"  HVG: {adata.var['highly_variable'].sum()} genes selected")
    return adata
```

## 2. 次元削減・クラスタリング

```python
def sc_clustering(adata, n_pcs=50, n_neighbors=15, resolution=1.0,
                   random_state=42):
    """
    PCA → 近傍グラフ → UMAP → Leiden クラスタリング。

    Parameters:
        n_pcs: PCA 主成分数
        n_neighbors: k-NN グラフの近傍数
        resolution: Leiden クラスタリングの解像度
    """
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, random_state=random_state)

    # Elbow plot 用の分散説明率
    variance_ratio = adata.uns["pca"]["variance_ratio"]
    cumvar = np.cumsum(variance_ratio)
    n_pcs_use = int(np.argmax(cumvar >= 0.9)) + 1
    n_pcs_use = max(n_pcs_use, 15)
    print(f"  PCA: using {n_pcs_use} PCs (cumulative variance ≥ 90%)")

    sc.pp.neighbors(adata, n_pcs=n_pcs_use, n_neighbors=n_neighbors,
                     random_state=random_state)
    sc.tl.umap(adata, random_state=random_state)
    sc.tl.leiden(adata, resolution=resolution, random_state=random_state)

    n_clusters = adata.obs["leiden"].nunique()
    print(f"  Leiden: {n_clusters} clusters (resolution={resolution})")
    return adata
```

## 3. 差次発現遺伝子（DEG）同定

```python
def sc_deg_analysis(adata, groupby="leiden", method="wilcoxon",
                     n_genes=200, min_logfc=0.25, max_pval=0.05):
    """
    クラスタ間の差次発現遺伝子を同定する。

    method:
      - "wilcoxon": Wilcoxon rank-sum test（推奨）
      - "t-test_overestim_var": Welch's t-test（高速）
      - "logreg": Logistic regression

    Returns:
        DataFrame with top DEGs per cluster
    """
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=n_genes)

    deg_results = []
    for cluster in adata.obs[groupby].unique():
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df["cluster"] = cluster
        df_sig = df[(df["logfoldchanges"].abs() >= min_logfc) &
                     (df["pvals_adj"] < max_pval)]
        deg_results.append(df_sig)

    deg_df = pd.concat(deg_results, ignore_index=True)
    print(f"  DEG: {len(deg_df)} significant genes across {adata.obs[groupby].nunique()} clusters")
    return deg_df
```

## 4. セルタイプアノテーション

```python
def annotate_celltypes_marker(adata, marker_dict, groupby="leiden",
                               threshold=0.5):
    """
    マーカー遺伝子ベースのセルタイプアノテーション。

    marker_dict 例:
        {
            "T cells": ["CD3D", "CD3E", "CD8A", "CD4"],
            "B cells": ["CD19", "MS4A1", "CD79A"],
            "Monocytes": ["CD14", "LYZ", "FCGR3A"],
            "NK cells": ["NKG7", "GNLY", "KLRD1"],
            "Dendritic": ["FCER1A", "CST3", "CLEC10A"],
        }

    各クラスタのマーカー遺伝子発現スコアを計算し、最も高いスコアに対応する
    セルタイプをアサインする。
    """
    sc.tl.score_genes(adata, gene_list=[], score_name="_dummy")

    scores = {}
    for cell_type, markers in marker_dict.items():
        valid_markers = [m for m in markers if m in adata.var_names]
        if valid_markers:
            score_name = f"score_{cell_type.replace(' ', '_')}"
            sc.tl.score_genes(adata, gene_list=valid_markers, score_name=score_name)
            scores[cell_type] = score_name

    # クラスタごとに最高スコアのセルタイプを割り当て
    cluster_annotations = {}
    for cluster in adata.obs[groupby].unique():
        mask = adata.obs[groupby] == cluster
        best_type = "Unknown"
        best_score = -np.inf
        for cell_type, score_col in scores.items():
            mean_score = adata.obs.loc[mask, score_col].mean()
            if mean_score > best_score and mean_score > threshold:
                best_score = mean_score
                best_type = cell_type
        cluster_annotations[cluster] = best_type

    adata.obs["cell_type"] = adata.obs[groupby].map(cluster_annotations)
    print(f"  Annotation: {len(set(cluster_annotations.values()))} cell types assigned")
    return adata, cluster_annotations
```

## 5. RNA Velocity

```python
def rna_velocity_analysis(adata_loom_path, adata, basis="umap"):
    """
    scVelo による RNA velocity 解析。

    RNA velocity は、スプライシング状態（unspliced/spliced）の比率から
    遺伝子発現の時間的変化方向を推定する。

    Modes:
      - stochastic: 確率的モデル（推奨、高速）
      - dynamical: 動的モデル（精度高、低速）
    """
    import scvelo as scv

    # 前処理
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    # Velocity 推定
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata)

    # 可視化
    scv.pl.velocity_embedding_stream(adata, basis=basis,
                                      save="figures/velocity_stream.png")

    # Latent time（擬似時間）
    scv.tl.latent_time(adata)
    print(f"  Velocity: latent time range [{adata.obs['latent_time'].min():.3f}, "
          f"{adata.obs['latent_time'].max():.3f}]")

    return adata
```

## 6. 細胞間コミュニケーション推定

```python
def cell_communication_analysis(adata, groupby="cell_type",
                                  database="CellChatDB.human"):
    """
    CellChat によるリガンド-レセプター相互作用解析。

    パイプライン:
      1. L-R ペアデータベースのロード
      2. 発現ベースのコミュニケーション確率計算
      3. シグナリングネットワーク推定
      4. 経路レベル集約
    """
    import cellchat

    cc = cellchat.CellChat(adata, groupby=groupby)
    cc.preprocess()
    cc.identify_overexpressed_genes()
    cc.identify_overexpressed_interactions()

    cc.compute_communication_prob(database=database)
    cc.filter_communication(min_cells=10)
    cc.compute_communication_prob_pathway()

    # ネットワーク可視化
    cc.aggregate_net()
    cc.net_analysis()

    # 結果取得
    lr_pairs = cc.get_significant_interactions()
    print(f"  CellChat: {len(lr_pairs)} significant L-R interactions")

    return cc, lr_pairs
```

## 7. 可視化

```python
def sc_visualization_panel(adata, deg_df=None, save_dir="figures"):
    """
    シングルセル解析結果の可視化パネル。

    生成図:
      1. QC violin plot (n_genes, total_counts, pct_mito)
      2. UMAP — leiden clusters
      3. UMAP — cell types
      4. DEG dot plot (top markers per cluster)
      5. Marker heatmap
    """
    import matplotlib.pyplot as plt
    import os
    os.makedirs(save_dir, exist_ok=True)

    # 1. QC violin
    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                  jitter=0.4, multi_panel=True, save="_qc.png")

    # 2. UMAP — clusters
    sc.pl.umap(adata, color="leiden", legend_loc="on data",
                title="Leiden Clusters", save="_leiden.png")

    # 3. UMAP — cell types
    if "cell_type" in adata.obs.columns:
        sc.pl.umap(adata, color="cell_type",
                    title="Cell Type Annotation", save="_celltypes.png")

    # 4. DEG dot plot
    if deg_df is not None:
        top_markers = deg_df.groupby("cluster").head(5)["names"].unique().tolist()
        sc.pl.dotplot(adata, var_names=top_markers[:30],
                       groupby="leiden", save="_markers.png")

    # 5. Stacked violin
    if deg_df is not None:
        top5 = deg_df.groupby("cluster").head(3)["names"].unique().tolist()[:20]
        sc.pl.stacked_violin(adata, var_names=top5,
                              groupby="leiden", save="_stacked.png")

    print(f"  Figures saved to {save_dir}/")
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/sc_qc_summary.json` | JSON |
| `results/sc_deg_results.csv` | CSV |
| `results/sc_celltype_annotations.json` | JSON |
| `results/sc_velocity_summary.json` | JSON |
| `results/sc_cellchat_interactions.csv` | CSV |
| `figures/umap_leiden.png` | PNG |
| `figures/umap_celltypes.png` | PNG |
| `figures/velocity_stream.png` | PNG |
| `figures/deg_dotplot.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| CELLxGENE | `CELLxGENE_get_expression_data` | シングルセル発現データ取得 |
| CELLxGENE | `CELLxGENE_get_cell_metadata` | 細胞メタデータ取得 |
| CELLxGENE | `CELLxGENE_get_gene_metadata` | 遺伝子メタデータ取得 |
| CELLxGENE | `CELLxGENE_get_presence_matrix` | 遺伝子存在マトリクス |
| CELLxGENE | `CELLxGENE_get_embeddings` | 埋め込みベクトル取得 |
| CELLxGENE | `CELLxGENE_download_h5ad` | H5AD ファイルダウンロード |
| HCA | `hca_search_projects` | Human Cell Atlas プロジェクト検索 |
| HCA | `hca_get_file_manifest` | HCA ファイルマニフェスト取得 |
| HPA | `HPA_get_rna_expression_by_source` | 組織別 RNA 発現データ |

### 参照スキル

| スキル | 連携内容 |
|---|---|
| [scientific-bioinformatics](../scientific-bioinformatics/SKILL.md) | 遺伝子アノテーション・パスウェイ解析 |
| [scientific-multi-omics](../scientific-multi-omics/SKILL.md) | マルチオミクス統合 |
| [scientific-network-analysis](../scientific-network-analysis/SKILL.md) | 遺伝子制御ネットワーク |
| [scientific-deep-learning](../scientific-deep-learning/SKILL.md) | scVI / scGPT 等の深層学習モデル |
| [scientific-pca-tsne](../scientific-pca-tsne/SKILL.md) | 次元削減手法 |

#### 依存パッケージ

- scanpy, anndata, scvelo, cellchat, leidenalg, umap-learn
