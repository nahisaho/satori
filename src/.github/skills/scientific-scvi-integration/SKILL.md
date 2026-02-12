---
name: scientific-scvi-integration
description: |
  scvi-tools シングルセル統合スキル。scVI 変分オートエンコーダ統合・
  scANVI 半教師有りアノテーション・totalVI CITE-seq
  RNA+タンパク質結合解析・SOLO ダブレット検出・潜在空間解析。
---

# Scientific scVI Integration

scvi-tools を活用したシングルセル確率的モデルベース統合
パイプラインを提供する。scVI/scANVI/totalVI/SOLO による
バッチ統合、半教師有りアノテーション、マルチモーダル解析。

## When to Use

- 複数バッチの scRNA-seq データを統合するとき (scVI)
- 半教師有りで細胞型アノテーションを転移するとき (scANVI)
- CITE-seq (RNA + ADT) データを結合解析するとき (totalVI)
- ダブレット (doublet) を検出・除去するとき (SOLO)
- 差次的発現を確率的にテストするとき
- 潜在空間を用いたクラスタリングを行うとき

---

## Quick Start

## 1. scVI モデルセットアップ & 訓練

```python
import scvi
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd


def setup_scvi(adata, batch_key="batch", layer=None, n_latent=30,
               n_hidden=128, n_layers=2):
    """
    scVI 変分オートエンコーダのセットアップ & 訓練。

    Parameters:
        adata: AnnData — 入力データ (raw counts)
        batch_key: str — バッチキー
        layer: str — カウント格納レイヤー
        n_latent: int — 潜在次元数
        n_hidden: int — 隠れユニット数
        n_layers: int — ニューラルネットレイヤー数

    K-Dense: scvi-tools
    """
    # scVI データ登録
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key=batch_key,
        layer=layer,
    )

    # モデル構築
    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_hidden=n_hidden,
        n_layers=n_layers,
        gene_likelihood="zinb",  # zero-inflated negative binomial
    )

    # 訓練
    model.train(max_epochs=400, early_stopping=True)

    # 潜在空間取得
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent

    print(f"scVI trained: {len(adata)} cells → {n_latent}D latent space")
    print(f"  Batches: {adata.obs[batch_key].nunique()}")
    return model
```

## 2. scVI バッチ統合 & UMAP

```python
def scvi_integration(adata, model, resolution=1.0):
    """
    scVI 潜在空間でバッチ統合 & クラスタリング。

    Parameters:
        adata: AnnData — scVI 登録済みデータ
        model: scvi.model.SCVI — 訓練済みモデル
        resolution: float — Leiden 解像度
    """
    # 潜在空間から近傍グラフ
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)

    n_clusters = adata.obs["leiden"].nunique()
    print(f"scVI integration: {n_clusters} clusters (resolution={resolution})")
    return adata
```

## 3. scANVI 半教師有りアノテーション

```python
def scanvi_annotation(adata, scvi_model, labels_key="cell_type",
                      unlabeled_category="Unknown", n_epochs=20):
    """
    scANVI で半教師有り細胞型アノテーション転移。

    Parameters:
        adata: AnnData — scVI 登録済みデータ
        scvi_model: scvi.model.SCVI — 訓練済み scVI
        labels_key: str — 既知ラベルカラム
        unlabeled_category: str — 未知ラベル値
        n_epochs: int — 追加訓練エポック
    """
    # scANVI = scVI + 半教師有り
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category=unlabeled_category,
        labels_key=labels_key,
    )

    scanvi_model.train(max_epochs=n_epochs)

    # 予測ラベル
    predictions = scanvi_model.predict()
    adata.obs["scANVI_prediction"] = predictions

    # 予測確率
    soft_predictions = scanvi_model.predict(soft=True)
    adata.obs["scANVI_confidence"] = soft_predictions.max(axis=1)

    n_labeled = (adata.obs[labels_key] != unlabeled_category).sum()
    n_unlabeled = (adata.obs[labels_key] == unlabeled_category).sum()
    mean_conf = adata.obs["scANVI_confidence"].mean()

    print(f"scANVI annotation:")
    print(f"  Labeled: {n_labeled}, Unlabeled: {n_unlabeled}")
    print(f"  Mean confidence: {mean_conf:.3f}")
    return scanvi_model
```

## 4. totalVI CITE-seq 統合

```python
def totalvi_citeseq(adata, protein_expression_obsm="protein_expression",
                    batch_key="batch", n_latent=20, n_epochs=400):
    """
    totalVI で RNA + ADT (CITE-seq) 結合解析。

    Parameters:
        adata: AnnData — RNA counts + protein expression
        protein_expression_obsm: str — ADT 格納 obsm キー
        batch_key: str — バッチキー
        n_latent: int — 潜在次元数
    """
    scvi.model.TOTALVI.setup_anndata(
        adata,
        batch_key=batch_key,
        protein_expression_obsm_key=protein_expression_obsm,
    )

    model = scvi.model.TOTALVI(
        adata,
        n_latent=n_latent,
        latent_distribution="normal",
    )

    model.train(max_epochs=n_epochs, early_stopping=True)

    # 潜在空間
    latent = model.get_latent_representation()
    adata.obsm["X_totalVI"] = latent

    # タンパク質前景確率 (denoised)
    _, protein_fore = model.get_normalized_expression(
        n_samples=25,
        return_mean=True,
    )

    n_proteins = protein_fore.shape[1] if hasattr(protein_fore, 'shape') else 0
    print(f"totalVI: {len(adata)} cells, {n_proteins} proteins")
    print(f"  Latent dim: {n_latent}")
    return model
```

## 5. SOLO ダブレット検出

```python
def solo_doublet_detection(adata, scvi_model, threshold=0.5):
    """
    SOLO でダブレット検出。

    Parameters:
        adata: AnnData — scVI 登録済みデータ
        scvi_model: scvi.model.SCVI — 訓練済み scVI
        threshold: float — ダブレット判定閾値
    """
    solo_model = scvi.external.SOLO.from_scvi_model(scvi_model)
    solo_model.train()

    # ダブレット予測
    predictions = solo_model.predict()
    predictions["label"] = predictions.apply(
        lambda x: "doublet" if x["doublet"] > threshold else "singlet",
        axis=1,
    )

    adata.obs["solo_doublet"] = predictions["label"].values
    adata.obs["solo_score"] = predictions["doublet"].values

    n_doublets = (adata.obs["solo_doublet"] == "doublet").sum()
    doublet_rate = n_doublets / len(adata) * 100

    print(f"SOLO doublet detection:")
    print(f"  Doublets: {n_doublets} ({doublet_rate:.1f}%)")
    print(f"  Singlets: {len(adata) - n_doublets}")
    return predictions
```

## 6. scVI 差次的発現

```python
def scvi_differential_expression(model, adata, groupby="leiden",
                                  group1="0", group2="1",
                                  delta=0.25, batch_size=256):
    """
    scVI による確率的差次的発現。

    Parameters:
        model: scvi.model.SCVI — 訓練済みモデル
        adata: AnnData — 入力データ
        groupby: str — グループカラム
        group1: str — 比較グループ1
        group2: str — 比較グループ2
        delta: float — LFC 閾値
    """
    de_results = model.differential_expression(
        groupby=groupby,
        group1=group1,
        group2=group2,
        delta=delta,
        batch_size=batch_size,
    )

    # 有意な遺伝子フィルタ
    sig_genes = de_results[
        (de_results["is_de_fdr_0.05"])
        & (de_results["lfc_mean"].abs() > delta)
    ]

    sig_genes = sig_genes.sort_values("lfc_mean", ascending=False)

    print(f"scVI DE: {group1} vs {group2}")
    print(f"  Significant genes: {len(sig_genes)}")
    print(f"  Up in {group1}: {(sig_genes['lfc_mean'] > 0).sum()}")
    print(f"  Up in {group2}: {(sig_genes['lfc_mean'] < 0).sum()}")
    return sig_genes
```

## 7. 統合パイプライン

```python
def scvi_pipeline(adata_paths, batch_labels=None, labels_key="cell_type",
                  output_dir="results"):
    """
    scVI → scANVI → SOLO 統合パイプライン。

    Parameters:
        adata_paths: list — AnnData ファイルパスリスト
        batch_labels: list — バッチラベル
        labels_key: str — 細胞型ラベルカラム
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) データ結合
    adatas = []
    for i, path in enumerate(adata_paths):
        a = sc.read_h5ad(path)
        a.obs["batch"] = batch_labels[i] if batch_labels else f"batch_{i}"
        adatas.append(a)

    adata = ad.concat(adatas, join="inner")
    adata.obs_names_make_unique()

    # 2) 前処理
    sc.pp.highly_variable_genes(
        adata, n_top_genes=2000, batch_key="batch", flavor="seurat_v3"
    )
    adata = adata[:, adata.var["highly_variable"]].copy()

    # 3) scVI 訓練
    scvi_model = setup_scvi(adata, batch_key="batch")

    # 4) SOLO ダブレット除去
    solo_results = solo_doublet_detection(adata, scvi_model)
    adata = adata[adata.obs["solo_doublet"] == "singlet"].copy()

    # 5) 再訓練 (ダブレット除去後)
    scvi_model = setup_scvi(adata, batch_key="batch")

    # 6) scANVI アノテーション
    if labels_key in adata.obs.columns:
        scanvi_model = scanvi_annotation(adata, scvi_model, labels_key=labels_key)

    # 7) UMAP & クラスタリング
    scvi_integration(adata, scvi_model)

    # 保存
    adata.write(output_dir / "integrated.h5ad")
    scvi_model.save(str(output_dir / "scvi_model"))

    print(f"Pipeline complete: {len(adata)} cells integrated")
    return adata, scvi_model
```

---

## パイプライン統合

```
single-cell-genomics → scvi-integration → spatial-transcriptomics
  (scRNA-seq QC)       (scVI/scANVI/totalVI)   (Visium/MERFISH)
       │                       │                      ↓
perturbation-analysis ────────┘               gene-expression
  (Perturb-seq)          │                    (DEG/マーカー)
                         ↓
                  expression-comparison
                  (Expression Atlas 比較)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/integrated.h5ad` | 統合 AnnData | → spatial-transcriptomics |
| `results/scvi_model/` | scVI 訓練済みモデル | → perturbation-analysis |
| `results/de_results.csv` | 差次的発現結果 | → gene-expression |
| `results/annotations.csv` | scANVI アノテーション | → expression-comparison |
