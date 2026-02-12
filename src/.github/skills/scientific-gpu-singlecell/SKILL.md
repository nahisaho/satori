---
name: scientific-gpu-singlecell
description: |
  GPU アクセラレーション シングルセル解析スキル。
  rapids-singlecell / cuML / cuGraph による GPU 並列処理。
  大規模 (>1M cells) データの高速前処理・クラスタリング・
  次元削減。K-Dense: rapids-singlecell。
---

# Scientific GPU Single-Cell

rapids-singlecell / cuML / cuGraph を活用した GPU アクセラレー
ション対応シングルセル解析パイプラインを提供する。100万細胞超
の大規模データセットの高速処理。

## When to Use

- 大規模シングルセルデータ (>100k cells) の高速前処理が必要なとき
- GPU クラスタリング (Leiden/Louvain) を実行するとき
- GPU UMAP/t-SNE で次元削減を高速化するとき
- CPU 版 scanpy では処理時間が実用的でないとき
- 複数サンプル統合に GPU を活用するとき
- ベンチマーク (CPU vs GPU) で性能比較を行うとき

---

## Quick Start

## 1. rapids-singlecell 前処理

```python
import rapids_singlecell as rsc
import scanpy as sc
import anndata as ad
import cupy as cp
import numpy as np
import pandas as pd
import time
from pathlib import Path


def gpu_preprocessing(adata, min_genes=200, min_cells=3,
                       n_top_genes=2000, target_sum=10000):
    """
    rapids-singlecell — GPU 前処理パイプライン。

    Parameters:
        adata: AnnData — 入力データ
        min_genes: int — 最小遺伝子数
        min_cells: int — 最小細胞数
        n_top_genes: int — HVG 数
        target_sum: float — 正規化ターゲット
    """
    t0 = time.time()

    # GPU メモリにデータ転送
    rsc.get.anndata_to_GPU(adata)

    # QC
    rsc.pp.calculate_qc_metrics(adata)
    rsc.pp.filter_cells(adata, min_genes=min_genes)
    rsc.pp.filter_genes(adata, min_cells=min_cells)

    # 正規化
    rsc.pp.normalize_total(adata, target_sum=target_sum)
    rsc.pp.log1p(adata)

    # HVG 選択
    rsc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor="seurat_v3",
    )

    # スケーリング
    rsc.pp.scale(adata, max_value=10)

    elapsed = time.time() - t0
    print(f"GPU preprocessing: {adata.n_obs} cells × {adata.n_vars} genes "
          f"({elapsed:.1f}s)")
    return adata
```

## 2. GPU PCA & 近傍グラフ

```python
def gpu_pca_neighbors(adata, n_comps=50, n_neighbors=15):
    """
    GPU PCA + 近傍グラフ構築。

    Parameters:
        adata: AnnData — 前処理済みデータ (GPU)
        n_comps: int — PCA 成分数
        n_neighbors: int — kNN 近傍数
    """
    t0 = time.time()

    # GPU PCA
    rsc.pp.pca(adata, n_comps=n_comps)

    # GPU kNN
    rsc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps)

    elapsed = time.time() - t0
    print(f"GPU PCA + kNN: {n_comps} PCs, k={n_neighbors} ({elapsed:.1f}s)")
    return adata
```

## 3. GPU クラスタリング (Leiden/Louvain)

```python
def gpu_clustering(adata, method="leiden", resolution=1.0):
    """
    cuGraph — GPU Leiden/Louvain クラスタリング。

    Parameters:
        adata: AnnData — 近傍グラフ付きデータ
        method: str — "leiden" or "louvain"
        resolution: float — クラスタリング解像度
    """
    t0 = time.time()

    if method == "leiden":
        rsc.tl.leiden(adata, resolution=resolution)
    else:
        rsc.tl.louvain(adata, resolution=resolution)

    n_clusters = adata.obs[method].nunique()
    elapsed = time.time() - t0

    print(f"GPU {method}: {n_clusters} clusters, "
          f"resolution={resolution} ({elapsed:.1f}s)")
    return adata
```

## 4. GPU UMAP / t-SNE

```python
def gpu_embedding(adata, method="umap", n_components=2, **kwargs):
    """
    GPU UMAP / t-SNE 次元削減。

    Parameters:
        adata: AnnData — 近傍グラフ付きデータ
        method: str — "umap" or "tsne"
        n_components: int — 出力次元数
    """
    t0 = time.time()

    if method == "umap":
        rsc.tl.umap(adata, n_components=n_components, **kwargs)
    else:
        rsc.tl.tsne(adata, n_pcs=n_components, **kwargs)

    elapsed = time.time() - t0
    print(f"GPU {method.upper()}: {n_components}D ({elapsed:.1f}s)")
    return adata
```

## 5. CPU vs GPU ベンチマーク

```python
def benchmark_cpu_vs_gpu(adata_path, n_top_genes=2000, n_comps=50):
    """
    CPU (scanpy) vs GPU (rapids-singlecell) ベンチマーク。

    Parameters:
        adata_path: str — h5ad ファイルパス
        n_top_genes: int — HVG 数
        n_comps: int — PCA 成分数
    """
    results = {}

    # === CPU (scanpy) ===
    adata_cpu = sc.read_h5ad(adata_path)
    t0 = time.time()
    sc.pp.normalize_total(adata_cpu, target_sum=1e4)
    sc.pp.log1p(adata_cpu)
    sc.pp.highly_variable_genes(adata_cpu, n_top_genes=n_top_genes)
    adata_cpu = adata_cpu[:, adata_cpu.var["highly_variable"]].copy()
    sc.pp.scale(adata_cpu, max_value=10)
    sc.pp.pca(adata_cpu, n_comps=n_comps)
    sc.pp.neighbors(adata_cpu)
    sc.tl.leiden(adata_cpu)
    sc.tl.umap(adata_cpu)
    cpu_time = time.time() - t0
    results["cpu_seconds"] = cpu_time
    results["cpu_clusters"] = adata_cpu.obs["leiden"].nunique()

    # === GPU (rapids-singlecell) ===
    adata_gpu = sc.read_h5ad(adata_path)
    t0 = time.time()
    rsc.get.anndata_to_GPU(adata_gpu)
    rsc.pp.normalize_total(adata_gpu, target_sum=1e4)
    rsc.pp.log1p(adata_gpu)
    rsc.pp.highly_variable_genes(adata_gpu, n_top_genes=n_top_genes,
                                  flavor="seurat_v3")
    adata_gpu = adata_gpu[:, adata_gpu.var["highly_variable"]].copy()
    rsc.pp.scale(adata_gpu, max_value=10)
    rsc.pp.pca(adata_gpu, n_comps=n_comps)
    rsc.pp.neighbors(adata_gpu)
    rsc.tl.leiden(adata_gpu)
    rsc.tl.umap(adata_gpu)
    gpu_time = time.time() - t0
    results["gpu_seconds"] = gpu_time
    results["gpu_clusters"] = adata_gpu.obs["leiden"].nunique()

    results["speedup"] = cpu_time / gpu_time
    results["n_cells"] = adata_cpu.n_obs

    print(f"Benchmark ({results['n_cells']} cells):")
    print(f"  CPU: {cpu_time:.1f}s ({results['cpu_clusters']} clusters)")
    print(f"  GPU: {gpu_time:.1f}s ({results['gpu_clusters']} clusters)")
    print(f"  Speedup: {results['speedup']:.1f}x")
    return results
```

## 6. GPU シングルセル統合パイプライン

```python
def gpu_singlecell_pipeline(input_files, output_dir="results",
                              n_top_genes=3000, resolution=1.0):
    """
    大規模 GPU シングルセル統合パイプライン。

    Parameters:
        input_files: list[str] — h5ad ファイルリスト
        output_dir: str — 出力ディレクトリ
        n_top_genes: int — HVG 数
        resolution: float — Leiden 解像度
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    t_total = time.time()

    # 1) データ読み込み・結合
    adatas = []
    for i, f in enumerate(input_files):
        a = sc.read_h5ad(f)
        a.obs["sample"] = f"sample_{i}"
        adatas.append(a)
    adata = ad.concat(adatas, join="inner")
    print(f"Combined: {adata.n_obs} cells from {len(input_files)} samples")

    # 2) GPU 前処理
    adata = gpu_preprocessing(adata, n_top_genes=n_top_genes)

    # 3) GPU PCA + kNN
    adata = gpu_pca_neighbors(adata)

    # 4) GPU クラスタリング
    adata = gpu_clustering(adata, resolution=resolution)

    # 5) GPU UMAP
    adata = gpu_embedding(adata)

    # 6) CPU に戻して marker 検出
    rsc.get.anndata_to_CPU(adata)
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    # 保存
    adata.write(output_dir / "gpu_singlecell.h5ad")

    # マーカー遺伝子エクスポート
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(output_dir / "markers.csv", index=False)

    total_time = time.time() - t_total
    print(f"GPU pipeline: {adata.n_obs} cells, "
          f"{adata.obs['leiden'].nunique()} clusters ({total_time:.1f}s)")
    return adata
```

---

## パイプライン統合

```
single-cell-genomics → gpu-singlecell → scvi-integration
  (scanpy 標準)        (GPU 高速化)     (深層学習統合)
       │                     │               ↓
  batch-correction ─────────┘         cell-type-annotation
  (Harmony/scVI)       │              (自動アノテーション)
                       ↓
                  atlas-construction
                  (大規模アトラス)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/gpu_singlecell.h5ad` | GPU 処理済み AnnData | → scvi-integration |
| `results/markers.csv` | マーカー遺伝子 | → cell-type-annotation |
| `results/benchmark.json` | CPU/GPU 比較結果 | → atlas-construction |
