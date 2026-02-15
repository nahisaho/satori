---
name: scientific-scatac-signac
description: |
  scATAC-seq 解析スキル (Signac/SnapATAC2/episcanpy)。
  ピークコーリング・モチーフ解析・Gene Activity スコア・
  RNA+ATAC マルチモーダル統合 (WNN)。K-Dense: signac。
tu_tools:
  - key: encode
    name: ENCODE
    description: scATAC-seq 参照エピゲノムデータ
---

# Scientific scATAC-seq / Signac

Signac / SnapATAC2 / episcanpy を活用した scATAC-seq (single-cell
ATAC-seq) 解析パイプラインを提供する。クロマチンアクセシビリティ
解析、モチーフエンリッチメント、マルチモーダル統合。

## When to Use

- scATAC-seq データの前処理とピークコーリングを行うとき
- 転写因子モチーフのエンリッチメント解析を実行するとき
- Gene Activity スコアで scATAC を発現レベルで解釈するとき
- scRNA-seq + scATAC-seq のマルチモーダル統合を行うとき
- クロマチンアクセシビリティの細胞型特異的パターンを同定するとき
- エピゲノム (ヒストン修飾/DNAメチル化) とクロマチンを統合するとき

---

## Quick Start

## 1. scATAC-seq 前処理 (SnapATAC2)

```python
import snapatac2 as snap
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path


def scatac_preprocessing(fragment_file, genome="hg38",
                          min_tsse=5, min_fragments=1000):
    """
    SnapATAC2 — scATAC-seq 前処理パイプライン。

    Parameters:
        fragment_file: str — fragments.tsv.gz パス
        genome: str — リファレンスゲノム ("hg38", "mm10")
        min_tsse: float — 最小 TSS enrichment スコア
        min_fragments: int — 最小フラグメント数
    """
    # フラグメントファイル読み込み
    adata = snap.pp.import_data(
        fragment_file,
        chrom_sizes=snap.genome(genome),
        sorted_by_barcode=False,
    )

    # QC メトリクス
    snap.metrics.tsse(adata, snap.genome(genome))
    snap.metrics.frag_size_distr(adata)

    # フィルタリング
    snap.pp.filter_cells(
        adata,
        min_counts=min_fragments,
        min_tsse=min_tsse,
    )

    print(f"scATAC preprocessing: {adata.n_obs} cells, "
          f"TSS enrichment ≥ {min_tsse}")
    return adata
```

## 2. ピークコーリング & Tile Matrix

```python
def scatac_peak_calling(adata, genome="hg38", peak_method="macs2",
                         n_features=50000):
    """
    scATAC-seq ピークコーリング & アクセシビリティマトリクス作成。

    Parameters:
        adata: AnnData — 前処理済み scATAC データ
        genome: str — リファレンスゲノム
        peak_method: str — "macs2" or "tile"
        n_features: int — feature selection 上位数
    """
    if peak_method == "tile":
        # Tile matrix (500bp bins)
        snap.pp.add_tile_matrix(adata, bin_size=500)
    else:
        # MACS2 peak calling (クラスタ別)
        snap.pp.make_peak_matrix(adata)

    # Feature selection
    snap.pp.select_features(adata, n_features=n_features)

    # 次元削減
    snap.tl.spectral(adata)
    snap.tl.umap(adata)

    # クラスタリング
    snap.pp.knn(adata)
    snap.tl.leiden(adata)

    n_clusters = adata.obs["leiden"].nunique()
    print(f"Peak calling ({peak_method}): {n_clusters} clusters, "
          f"{n_features} features")
    return adata
```

## 3. モチーフエンリッチメント (chromVAR)

```python
def motif_enrichment(adata, genome="hg38", motif_db="JASPAR2022"):
    """
    chromVAR モチーフエンリッチメント解析。

    Parameters:
        adata: AnnData — ピークマトリクス付き scATAC データ
        genome: str — リファレンスゲノム
        motif_db: str — モチーフデータベース
    """
    # モチーフスキャン
    snap.tl.motif_enrichment(
        adata,
        motifs=motif_db,
        genome=genome,
    )

    # クラスタ別差分モチーフ
    cluster_motifs = {}
    for cluster in adata.obs["leiden"].unique():
        mask = adata.obs["leiden"] == cluster
        # chromVAR deviation を取得
        if "chromvar" in adata.obsm:
            deviations = adata.obsm["chromvar"][mask].mean(axis=0)
            top_idx = np.argsort(deviations)[-10:]
            top_motifs = [adata.uns["motif_names"][i] for i in top_idx]
            cluster_motifs[cluster] = top_motifs

    results = []
    for cluster, motifs in cluster_motifs.items():
        for rank, motif in enumerate(motifs, 1):
            results.append({
                "cluster": cluster,
                "rank": rank,
                "motif": motif,
            })

    df = pd.DataFrame(results)
    print(f"Motif enrichment: {len(cluster_motifs)} clusters, "
          f"{len(df)} motif-cluster pairs")
    return df
```

## 4. Gene Activity スコア

```python
def gene_activity_score(adata, genome="hg38", upstream=2000, body=True):
    """
    Gene Activity スコア計算 — ATAC → 擬似発現量。

    Parameters:
        adata: AnnData — scATAC データ
        genome: str — リファレンスゲノム
        upstream: int — プロモーター上流距離 (bp)
        body: bool — 遺伝子本体を含むか
    """
    snap.pp.make_gene_matrix(
        adata,
        gene_anno=snap.genome(genome),
        upstream=upstream,
        include_body=body,
    )

    # Gene Activity を .layers に保存
    if hasattr(adata, "uns") and "gene_activity" in adata.uns:
        gene_act = adata.uns["gene_activity"]
    else:
        gene_act = adata.X  # Gene matrix mode

    # 正規化
    gene_act_norm = gene_act / gene_act.sum(axis=1, keepdims=True) * 10000
    gene_act_log = np.log1p(gene_act_norm)

    print(f"Gene activity: {gene_act_log.shape[1]} genes, "
          f"{gene_act_log.shape[0]} cells")
    return gene_act_log
```

## 5. RNA + ATAC マルチモーダル統合 (WNN)

```python
import scanpy as sc


def multimodal_wnn(adata_atac, adata_rna, n_neighbors=20):
    """
    Weighted Nearest Neighbor (WNN) — RNA + ATAC 統合。

    Parameters:
        adata_atac: AnnData — scATAC データ (LSI 済み)
        adata_rna: AnnData — scRNA データ (PCA 済み)
        n_neighbors: int — 近傍数
    """
    # 共通バーコード
    common_bc = list(
        set(adata_atac.obs_names) & set(adata_rna.obs_names)
    )
    atac_sub = adata_atac[common_bc].copy()
    rna_sub = adata_rna[common_bc].copy()

    print(f"Common barcodes: {len(common_bc)}")

    # 各モダリティの kNN
    sc.pp.neighbors(rna_sub, use_rep="X_pca", key_added="rna")
    sc.pp.neighbors(atac_sub, use_rep="X_spectral", key_added="atac")

    # WNN 統合 (muon)
    try:
        import muon as mu
        mdata = mu.MuData({"rna": rna_sub, "atac": atac_sub})
        mu.pp.neighbors(mdata, key="wnn")
        mu.tl.umap(mdata, neighbors_key="wnn")
        mu.tl.leiden(mdata, neighbors_key="wnn")

        n_clusters = mdata.obs["leiden"].nunique()
        print(f"WNN integration: {n_clusters} clusters")
        return mdata
    except ImportError:
        print("muon not installed — falling back to concatenation")
        # Fallback: 単純連結
        combined = ad.concat([rna_sub, atac_sub], axis=1, merge="same")
        sc.pp.neighbors(combined)
        sc.tl.leiden(combined)
        return combined
```

## 6. scATAC-seq 統合パイプライン

```python
def scatac_pipeline(fragment_file, rna_h5ad=None, genome="hg38",
                     output_dir="results"):
    """
    scATAC-seq 統合解析パイプライン。

    Parameters:
        fragment_file: str — fragments.tsv.gz
        rna_h5ad: str — scRNA-seq h5ad (マルチモーダル用, optional)
        genome: str — リファレンスゲノム
        output_dir: str — 出力ディレクトリ
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 前処理
    adata = scatac_preprocessing(fragment_file, genome=genome)

    # 2) ピーク & クラスタリング
    adata = scatac_peak_calling(adata, genome=genome)
    adata.write(output_dir / "scatac_clustered.h5ad")

    # 3) モチーフ
    motifs = motif_enrichment(adata, genome=genome)
    motifs.to_csv(output_dir / "motif_enrichment.csv", index=False)

    # 4) Gene Activity
    gene_act = gene_activity_score(adata, genome=genome)

    # 5) マルチモーダル統合
    if rna_h5ad:
        adata_rna = sc.read_h5ad(rna_h5ad)
        mdata = multimodal_wnn(adata, adata_rna)
        mdata.write(output_dir / "multimodal_wnn.h5mu")

    print(f"scATAC pipeline: {output_dir}")
    return adata
```

---

## パイプライン統合

```
epigenomics-chromatin → scatac-signac → single-cell-genomics
  (ChIP/ATAC bulk)     (scATAC-seq)    (scRNA 統合)
       │                     │               ↓
  peak-annotation ──────────┘         spatial-transcriptomics
  (ENCODE/ChIPAtlas)   │              (Visium/MERFISH)
                       ↓
                  gene-regulatory-network
                  (GRN 推定)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/scatac_clustered.h5ad` | クラスタリング済み scATAC | → single-cell-genomics |
| `results/motif_enrichment.csv` | モチーフエンリッチメント | → gene-regulatory-network |
| `results/multimodal_wnn.h5mu` | RNA+ATAC 統合 | → spatial-transcriptomics |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `encode` | ENCODE | scATAC-seq 参照エピゲノムデータ |
