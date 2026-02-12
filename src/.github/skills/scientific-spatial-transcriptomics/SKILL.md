---
name: scientific-spatial-transcriptomics
description: |
  空間トランスクリプトミクス解析スキル。10x Visium / MERFISH / Slide-seq データの
  前処理・空間的遺伝子発現パターン検出（Moran's I / SpatialDE）・
  空間ドメイン同定（BayesSpace / STAGATE）・細胞-細胞近接解析・
  Deconvolution（cell2location）パイプライン。Squidpy フレームワーク準拠。
---

# Scientific Spatial Transcriptomics

空間トランスクリプトミクス技術（10x Visium, MERFISH, Slide-seq, STARmap 等）で
取得したデータの空間的遺伝子発現解析パイプラインを提供する。
空間統計・空間ドメイン検出・deconvolution・リガンド-レセプター空間共発現を扱う。

## When to Use

- Visium / MERFISH 等の空間トランスクリプトミクスデータの解析が必要なとき
- 空間的に変動する遺伝子（SVG: Spatially Variable Genes）を同定するとき
- 組織内の空間ドメインを自動検出するとき
- スポットの細胞型 deconvolution を行うとき
- 空間的なリガンド-レセプター相互作用を評価するとき

---

## Quick Start

## 1. 空間データ前処理

```python
import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd

def spatial_preprocessing(adata, min_counts=500, min_genes=200,
                           max_pct_mito=25, n_top_genes=3000):
    """
    空間トランスクリプトミクス前処理パイプライン。

    手順:
      1. QC メトリクス計算
      2. 低品質スポット/セルフィルタリング
      3. 正規化 + log1p
      4. HVG 選択
      5. 空間近傍グラフ構築
    """
    # ミトコンドリア遺伝子
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["pct_counts_mt"] <= max_pct_mito].copy()
    sc.pp.filter_genes(adata, min_cells=10)
    print(f"  QC: {n_before} → {adata.n_obs} spots")

    # 正規化
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat_v3",
                                 layer="counts")

    # 空間近傍グラフ
    sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)
    print(f"  Spatial neighbors: {adata.obsp['spatial_connectivities'].nnz} edges")

    return adata
```

## 2. 空間的変動遺伝子（SVG）検出

```python
def spatially_variable_genes(adata, method="moran", n_perms=100,
                               fdr_threshold=0.05):
    """
    空間的に変動する遺伝子を同定する。

    method:
      - "moran": Moran's I 空間自己相関統計量
        I = (N / W) * Σᵢ Σⱼ wᵢⱼ (xᵢ - x̄)(xⱼ - x̄) / Σᵢ (xᵢ - x̄)²
      - "sepal": SEPAL（拡散ベース）

    Moran's I:
      - I ≈ 1: 強い正の空間自己相関（同値クラスタリング）
      - I ≈ 0: ランダム分布
      - I ≈ -1: 負の空間自己相関（チェッカーボード）
    """
    if method == "moran":
        sq.gr.spatial_autocorr(adata, mode="moran", n_perms=n_perms, n_jobs=-1)
        svg_df = adata.uns["moranI"].copy()
        svg_df["significant"] = svg_df["pval_norm_fdr_bh"] < fdr_threshold
        svg_df = svg_df.sort_values("I", ascending=False)
    elif method == "sepal":
        sq.gr.spatial_autocorr(adata, mode="geary", n_perms=n_perms, n_jobs=-1)
        svg_df = adata.uns["gearyC"].copy()
        svg_df["significant"] = svg_df["pval_norm_fdr_bh"] < fdr_threshold

    n_sig = svg_df["significant"].sum()
    print(f"  SVG ({method}): {n_sig} significant spatially variable genes")
    return svg_df
```

## 3. 空間ドメイン検出

```python
def spatial_domain_detection(adata, n_domains=7, method="leiden", resolution=0.8):
    """
    空間ドメイン（組織領域）を自動検出する。

    method:
      - "leiden": 空間グラフベース Leiden クラスタリング
      - "bayesspace": BayesSpace（ベイズモデル、Visium 特化）

    空間的制約付きクラスタリングにより、遺伝子発現が類似し空間的にも
    隣接するスポットを同一ドメインにグループ化する。
    """
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30)

    if method == "leiden":
        # 空間近傍グラフを利用
        sc.tl.leiden(adata, resolution=resolution, adjacency=adata.obsp["spatial_connectivities"])
        adata.obs["spatial_domain"] = adata.obs["leiden"]
    elif method == "bayesspace":
        from bayesspace import BayesSpace
        bs = BayesSpace(adata, n_clusters=n_domains)
        bs.fit()
        adata.obs["spatial_domain"] = bs.labels_.astype(str)

    n_domains_found = adata.obs["spatial_domain"].nunique()
    print(f"  Domains: {n_domains_found} spatial domains detected ({method})")
    return adata
```

## 4. Cell-type Deconvolution

```python
def spatial_deconvolution(adata_spatial, adata_sc, cell_type_col="cell_type",
                            method="cell2location"):
    """
    空間スポットの細胞型組成を scRNA-seq 参照データから推定する。

    method:
      - "cell2location": Bayesian 推定（推奨、高精度）
      - "rctd": RCTD（ロバスト deconvolution）

    cell2location モデル:
      y_s ~ NB(μ_s, α)
      μ_s = Σ_f w_sf * g_f
      w_sf: スポット s における細胞型 f の存在量
      g_f: 細胞型 f の遺伝子発現シグネチャ
    """
    import cell2location

    # リファレンスシグネチャの学習
    cell2location.models.RegressionModel.setup_anndata(adata_sc, labels_key=cell_type_col)
    ref_model = cell2location.models.RegressionModel(adata_sc)
    ref_model.train(max_epochs=250)
    inf_aver = ref_model.export_posterior(adata_sc)

    # Spatial mapping
    cell2location.models.Cell2location.setup_anndata(adata_spatial)
    mod = cell2location.models.Cell2location(adata_spatial,
                                              cell_state_df=inf_aver)
    mod.train(max_epochs=30000)

    adata_spatial = mod.export_posterior(adata_spatial)
    print(f"  Deconvolution: {len(inf_aver.columns)} cell types mapped to "
          f"{adata_spatial.n_obs} spots")
    return adata_spatial
```

## 5. 空間リガンド-レセプター解析

```python
def spatial_ligand_receptor(adata, cluster_key="spatial_domain",
                              n_perms=1000, fdr_threshold=0.01):
    """
    空間的に隣接するクラスタ間のリガンド-レセプター相互作用を検定する。

    Squidpy の permutation test:
      各 L-R ペアについて、空間的に隣接するスポット間の共発現が
      ランダム置換と比較して有意に高いかを検定する。
    """
    sq.gr.ligrec(adata, n_perms=n_perms, cluster_key=cluster_key,
                  copy=False, use_raw=False, transmitter_params={"categories": "ligand"},
                  receiver_params={"categories": "receptor"})

    lr_results = adata.uns["ligrec"]
    pvals = lr_results["pvalues"]
    sig_count = (pvals < fdr_threshold).sum().sum()
    print(f"  L-R analysis: {sig_count} significant spatial interactions")
    return lr_results
```

## 6. 空間可視化

```python
def spatial_visualization_panel(adata, svg_df=None, save_dir="figures"):
    """
    空間トランスクリプトミクス可視化パネル。

    生成図:
      1. 空間遺伝子発現マップ（top SVG）
      2. 空間ドメインマップ
      3. Deconvolution 結果マップ
      4. Spatial autocorrelation ヒートマップ
    """
    import matplotlib.pyplot as plt
    import os
    os.makedirs(save_dir, exist_ok=True)

    # 1. 空間上の遺伝子発現
    if svg_df is not None:
        top_svgs = svg_df[svg_df["significant"]].head(6).index.tolist()
        sq.pl.spatial_scatter(adata, color=top_svgs, ncols=3,
                               save=f"{save_dir}/spatial_svgs.png")

    # 2. 空間ドメイン
    if "spatial_domain" in adata.obs.columns:
        sq.pl.spatial_scatter(adata, color="spatial_domain",
                               save=f"{save_dir}/spatial_domains.png")

    # 3. Moran's I 分布
    if svg_df is not None:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(svg_df["I"], bins=50, color="steelblue", alpha=0.7)
        ax.axvline(x=0, color="red", linestyle="--", label="Random (I=0)")
        ax.set_xlabel("Moran's I")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of Spatial Autocorrelation")
        ax.legend()
        plt.tight_layout()
        plt.savefig(f"{save_dir}/morans_i_dist.png", dpi=300, bbox_inches="tight")
        plt.close()

    print(f"  Figures saved to {save_dir}/")
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/svg_results.csv` | CSV |
| `results/spatial_domains.json` | JSON |
| `results/deconvolution_proportions.csv` | CSV |
| `results/ligrec_results.json` | JSON |
| `figures/spatial_svgs.png` | PNG |
| `figures/spatial_domains.png` | PNG |
| `figures/morans_i_dist.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| CELLxGENE | `CELLxGENE_get_expression_data` | 参照シングルセルデータ取得 |
| CELLxGENE | `CELLxGENE_get_cell_metadata` | 細胞メタデータ参照 |
| CELLxGENE | `CELLxGENE_download_h5ad` | H5AD データダウンロード |
| HCA | `hca_search_projects` | HCA 空間データ検索 |
| HPA | `HPA_get_rna_expression_by_source` | 組織発現参照データ |

### 参照スキル

| スキル | 連携内容 |
|---|---|
| [scientific-single-cell-genomics](../scientific-single-cell-genomics/SKILL.md) | scRNA-seq 参照データ・deconvolution |
| [scientific-bioinformatics](../scientific-bioinformatics/SKILL.md) | 遺伝子アノテーション・パスウェイ濃縮 |
| [scientific-image-analysis](../scientific-image-analysis/SKILL.md) | 組織画像処理・セグメンテーション |
| [scientific-network-analysis](../scientific-network-analysis/SKILL.md) | 空間近傍ネットワーク解析 |
| [scientific-bayesian-statistics](../scientific-bayesian-statistics/SKILL.md) | BayesSpace ドメイン検出 |

#### 依存パッケージ

- squidpy, scanpy, anndata, cell2location, bayesspace, SpatialDE
