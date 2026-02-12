---
name: scientific-perturbation-analysis
description: |
  シングルセル摂動解析スキル。pertpy による CRISPR スクリーン解析・
  薬剤応答分析・scGen 摂動予測・Augur 摂動応答性スコアリング・
  scIB 統合ベンチマーク・差次的摂動応答パイプライン。
---

# Scientific Perturbation Analysis

pertpy / Augur / scIB を活用したシングルセルレベルの摂動解析
パイプラインを提供する。CRISPR スクリーン、薬剤処理、
遺伝子ノックダウンなどの摂動データの統合解析。

## When to Use

- CRISPR スクリーンデータ (Perturb-seq) を解析するとき
- 薬剤処理前後のシングルセル発現変動を評価するとき
- 摂動応答の細胞型特異性を定量するとき
- 複数のバッチ統合手法をベンチマークするとき (scIB)
- 摂動の効果を in silico で予測するとき (scGen)
- 差次的優先度 (Augur) で摂動応答性の高い細胞型を特定するとき

---

## Quick Start

## 1. pertpy セットアップ & データ読込み

```python
import pertpy as pt
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np


def load_perturbation_data(adata_path, perturbation_key="perturbation",
                           control_label="control"):
    """
    摂動実験 AnnData 読込み & 前処理。

    Parameters:
        adata_path: str — AnnData ファイルパス
        perturbation_key: str — 摂動ラベルカラム
        control_label: str — コントロールラベル

    K-Dense: pertpy
    """
    adata = sc.read_h5ad(adata_path)

    # 基本前処理
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    n_perturbations = adata.obs[perturbation_key].nunique()
    n_control = (adata.obs[perturbation_key] == control_label).sum()
    n_perturbed = len(adata) - n_control

    print(f"Loaded: {len(adata)} cells, {n_perturbations} perturbations")
    print(f"Control: {n_control}, Perturbed: {n_perturbed}")
    return adata
```

## 2. 差次的遺伝子発現 (摂動 vs コントロール)

```python
def differential_perturbation(adata, perturbation_key="perturbation",
                               control="control", target=None):
    """
    摂動-コントロール間差次的発現解析。

    Parameters:
        adata: AnnData — 摂動データ
        perturbation_key: str — 摂動ラベル
        control: str — コントロールラベル
        target: str — 比較対象摂動 (None で全摂動)
    """
    if target:
        mask = adata.obs[perturbation_key].isin([control, target])
        adata_sub = adata[mask].copy()
    else:
        adata_sub = adata.copy()

    sc.tl.rank_genes_groups(
        adata_sub,
        groupby=perturbation_key,
        reference=control,
        method="wilcoxon",
    )

    results = {}
    for group in adata_sub.obs[perturbation_key].unique():
        if group == control:
            continue
        try:
            degs = sc.get.rank_genes_groups_df(adata_sub, group=group)
            degs_sig = degs[degs["pvals_adj"] < 0.05]
            results[group] = {
                "n_degs": len(degs_sig),
                "n_up": (degs_sig["logfoldchanges"] > 0).sum(),
                "n_down": (degs_sig["logfoldchanges"] < 0).sum(),
                "top_genes": degs_sig.head(10)["names"].tolist(),
            }
        except Exception:
            continue

    print(f"DE results: {len(results)} perturbations analyzed")
    return results
```

## 3. Augur 摂動応答性スコアリング

```python
def augur_prioritization(adata, perturbation_key="perturbation",
                         cell_type_key="cell_type", control="control"):
    """
    Augur で細胞型ごとの摂動応答性をスコアリング。

    Parameters:
        adata: AnnData — 摂動データ
        perturbation_key: str — 摂動ラベル
        cell_type_key: str — 細胞型ラベル
        control: str — コントロールラベル

    K-Dense: augur (via pertpy)
    """
    ag = pt.tl.Augur(estimator="random_forest_classifier")

    # 摂動 vs コントロールで各細胞型のAUC計算
    adata_augur, results = ag.predict(
        adata,
        condition_key=perturbation_key,
        cell_type_key=cell_type_key,
        control_label=control,
    )

    # 結果をDataFrameに
    auc_df = results["summary_metrics"]
    auc_df = auc_df.sort_values("auc", ascending=False)

    print(f"Augur prioritization:")
    for _, row in auc_df.head(5).iterrows():
        print(f"  {row['cell_type']}: AUC={row['auc']:.3f}")

    return auc_df
```

## 4. scGen 摂動予測

```python
def scgen_perturbation_prediction(adata, perturbation_key="perturbation",
                                   cell_type_key="cell_type",
                                   control="control", target_perturbation=None,
                                   target_cell_type=None):
    """
    scGen による摂動効果の in silico 予測。

    Parameters:
        adata: AnnData — 訓練データ
        target_perturbation: str — 予測対象の摂動
        target_cell_type: str — 予測対象の細胞型
    """
    import scgen

    # モデル訓練
    scg = scgen.SCGEN(adata)
    scg.train(max_epochs=100, batch_size=32)

    # 予測
    pred, delta = scg.predict(
        ctrl_key=control,
        stim_key=target_perturbation,
        celltype_to_predict=target_cell_type,
    )

    print(f"scGen prediction: {target_cell_type} under {target_perturbation}")
    print(f"  Predicted cells: {pred.shape[0]}")
    return pred, delta
```

## 5. scIB 統合ベンチマーク

```python
def benchmark_integration(adata, batch_key="batch", label_key="cell_type",
                           methods=None):
    """
    scIB でバッチ統合手法をベンチマーク。

    Parameters:
        adata: AnnData — バッチ混在データ
        batch_key: str — バッチラベル
        label_key: str — 細胞型ラベル
        methods: list — 評価するメトリクス

    K-Dense: scib
    """
    import scib

    if methods is None:
        methods = ["scib"]

    # 基本メトリクス
    metrics = {}

    # batch correction metrics
    metrics["batch_kbet"] = scib.me.kBET(
        adata, batch_key=batch_key, label_key=label_key
    )
    metrics["batch_silhouette"] = scib.me.silhouette_batch(
        adata, batch_key=batch_key, label_key=label_key, embed="X_pca"
    )

    # bio conservation metrics
    metrics["bio_nmi"] = scib.me.nmi(adata, label_key, "leiden")
    metrics["bio_ari"] = scib.me.ari(adata, label_key, "leiden")
    metrics["bio_silhouette"] = scib.me.silhouette(
        adata, label_key=label_key, embed="X_pca"
    )

    # 総合スコア
    metrics["overall"] = 0.6 * np.mean([
        metrics["bio_nmi"], metrics["bio_ari"], metrics["bio_silhouette"]
    ]) + 0.4 * np.mean([
        metrics["batch_kbet"], metrics["batch_silhouette"]
    ])

    print(f"scIB benchmark:")
    for k, v in metrics.items():
        print(f"  {k}: {v:.4f}")
    return metrics
```

## 6. 摂動シグネチャ解析

```python
def perturbation_signature(adata, perturbation_key="perturbation",
                            control="control", n_top_genes=50):
    """
    摂動特異的遺伝子シグネチャ抽出。

    Parameters:
        adata: AnnData — 摂動データ
        perturbation_key: str — 摂動ラベル
        control: str — コントロールラベル
        n_top_genes: int — トップ遺伝子数
    """
    perturbations = [p for p in adata.obs[perturbation_key].unique()
                     if p != control]

    signatures = {}
    ctrl_mean = adata[adata.obs[perturbation_key] == control].X.mean(axis=0)
    ctrl_mean = np.asarray(ctrl_mean).flatten()

    for pert in perturbations:
        pert_mask = adata.obs[perturbation_key] == pert
        pert_mean = adata[pert_mask].X.mean(axis=0)
        pert_mean = np.asarray(pert_mean).flatten()

        delta = pert_mean - ctrl_mean
        gene_indices = np.argsort(np.abs(delta))[::-1][:n_top_genes]

        signatures[pert] = {
            "top_genes": adata.var_names[gene_indices].tolist(),
            "deltas": delta[gene_indices].tolist(),
            "n_cells": int(pert_mask.sum()),
        }

    print(f"Signatures extracted: {len(signatures)} perturbations, "
          f"{n_top_genes} genes each")
    return signatures
```

---

## パイプライン統合

```
single-cell-genomics → perturbation-analysis → pathway-enrichment
  (scRNA-seq QC)        (摂動 DE/Augur/scGen)   (KEGG/Reactome)
        │                        │                     ↓
spatial-transcriptomics ──┘      │            disease-research
  (Visium/MERFISH)               ↓              (GWAS/DisGeNET)
                       drug-target-profiling
                       (標的候補評価)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/perturbation_de.json` | 差次的発現結果 | → pathway-enrichment |
| `results/augur_scores.csv` | Augur 応答性スコア | → single-cell-genomics |
| `results/perturbation_signatures.json` | 摂動シグネチャ | → drug-target-profiling |
| `results/scib_benchmark.json` | 統合ベンチマーク | → spatial-transcriptomics |
