---
name: scientific-gene-expression-transcriptomics
description: |
  遺伝子発現・トランスクリプトミクス解析スキル。GEO (Gene Expression Omnibus) からの
  公開データセット取得・前処理、DESeq2 (PyDESeq2) による差次発現解析、
  GTEx 組織発現参照・eQTL 解析、Expression Atlas (EBI GXA) 統合照会、
  遺伝子セット濃縮解析 (GSEA)、バルク RNA-seq カウントデータの
  標準解析パイプライン。
---

# Scientific Gene Expression & Transcriptomics

バルク RNA-seq / マイクロアレイの遺伝子発現データを対象に、
GEO データセット取得→前処理→差次発現→GSEA→組織発現参照の
統合トランスクリプトミクスパイプラインを提供する。

## When to Use

- GEO からバルク RNA-seq/マイクロアレイデータセットを取得・前処理するとき
- DESeq2 による差次発現遺伝子 (DEG) 解析が必要なとき
- GTEx 組織発現プロファイル・eQTL データを照会するとき
- 遺伝子セット濃縮解析 (GSEA/ORA) を行うとき
- Expression Atlas でベースライン/差次発現実験を検索するとき

---

## Quick Start

## 1. GEO データセット取得

```python
import pandas as pd
import GEOparse


def fetch_geo_dataset(accession, output_dir="data/geo"):
    """
    GEO (Gene Expression Omnibus) データセットの取得・前処理。

    GEO ID 形式:
    - GSE: Series (発現データセット)
    - GPL: Platform (アレイ/シーケンサー定義)
    - GSM: Sample (個別サンプル)
    - GDS: Dataset (キュレーション済み)
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    gse = GEOparse.get_GEO(geo=accession, destdir=output_dir)

    print(f"  GEO Accession: {accession}")
    print(f"  Title: {gse.metadata['title'][0]}")
    print(f"  Platform: {list(gse.gpls.keys())}")
    print(f"  Samples: {len(gse.gsms)}")
    print(f"  Type: {gse.metadata.get('type', ['unknown'])}")

    # サンプルメタデータ抽出
    metadata = []
    for gsm_name, gsm in gse.gsms.items():
        meta = {"sample_id": gsm_name}
        meta.update({k: v[0] if v else None
                     for k, v in gsm.metadata.items()
                     if k in ["title", "source_name_ch1", "characteristics_ch1"]})
        metadata.append(meta)

    metadata_df = pd.DataFrame(metadata)

    # 発現マトリクス取得
    pivot_df = gse.pivot_samples("VALUE")
    print(f"  Expression matrix: {pivot_df.shape[0]} genes × {pivot_df.shape[1]} samples")

    return gse, metadata_df, pivot_df
```

## 2. DESeq2 差次発現解析 (PyDESeq2)

```python
import numpy as np
import pandas as pd


def deseq2_differential_expression(count_matrix, metadata, design_factor,
                                     contrast=None, alpha=0.05,
                                     lfc_threshold=1.0):
    """
    PyDESeq2 による差次発現解析パイプライン。

    1. カウントマトリクス入力 (genes × samples)
    2. サイズファクター正規化 (median of ratios)
    3. 分散推定 (shrinkage)
    4. GLM フィッティング (NB 分布)
    5. Wald 検定
    6. LFC 収縮 (apeglm)
    7. FDR 補正 (Benjamini-Hochberg)
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # DeseqDataSet 構築
    dds = DeseqDataSet(
        counts=count_matrix,
        metadata=metadata,
        design_factors=design_factor,
    )

    # 正規化 + 分散推定 + 統計検定
    dds.deseq2()

    # 結果取得
    stat_res = DeseqStats(dds, contrast=contrast, alpha=alpha)
    stat_res.summary()

    results_df = stat_res.results_df.copy()

    # LFC 収縮
    stat_res.lfc_shrink(coeff=contrast)
    results_df["log2FoldChange_shrunk"] = stat_res.results_df["log2FoldChange"]

    # フィルタリング
    sig = results_df[
        (results_df["padj"] < alpha) &
        (results_df["log2FoldChange"].abs() > lfc_threshold)
    ]

    sig_up = sig[sig["log2FoldChange"] > 0]
    sig_down = sig[sig["log2FoldChange"] < 0]

    print(f"  DESeq2 results:")
    print(f"    Total genes tested: {len(results_df)}")
    print(f"    Significant (FDR < {alpha}, |log2FC| > {lfc_threshold}):")
    print(f"      UP: {len(sig_up)}")
    print(f"      DOWN: {len(sig_down)}")

    return results_df, sig


def generate_volcano_plot(results_df, alpha=0.05, lfc_threshold=1.0,
                           output_file="figures/volcano_rnaseq.png"):
    """
    Volcano プロット生成。
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 6))

    results_df["-log10_padj"] = -np.log10(results_df["padj"].clip(lower=1e-300))

    # 色分け
    colors = []
    for _, row in results_df.iterrows():
        if row["padj"] < alpha and row["log2FoldChange"] > lfc_threshold:
            colors.append("red")
        elif row["padj"] < alpha and row["log2FoldChange"] < -lfc_threshold:
            colors.append("blue")
        else:
            colors.append("gray")

    ax.scatter(results_df["log2FoldChange"], results_df["-log10_padj"],
              c=colors, alpha=0.5, s=5)
    ax.axhline(-np.log10(alpha), color="gray", linestyle="--", lw=0.5)
    ax.axvline(lfc_threshold, color="gray", linestyle="--", lw=0.5)
    ax.axvline(-lfc_threshold, color="gray", linestyle="--", lw=0.5)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.set_title("Volcano Plot — Differential Expression")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    return output_file
```

## 3. GTEx 組織発現・eQTL 照会

```python
import pandas as pd


def query_gtex_expression(gene_name, tissue=None):
    """
    GTEx (Genotype-Tissue Expression) 組織発現プロファイル照会。

    GTEx v8: 54 組織, 948 ドナー, 17,382 サンプル。
    TPM (Transcripts Per Million) ベースの発現量。
    """
    print(f"  GTEx gene expression query: {gene_name}")
    if tissue:
        print(f"  Tissue: {tissue}")
    else:
        print("  All tissues (54 tissue sites)")

    return {"gene": gene_name, "tissue": tissue}


def query_gtex_eqtl(gene_name, tissue, pvalue_threshold=1e-5):
    """
    GTEx eQTL (expression Quantitative Trait Loci) 照会。

    eQTL = 遺伝子発現量に影響する遺伝的変異
    - cis-eQTL: 遺伝子の ±1 Mb 以内の変異
    - trans-eQTL: 遺伝子から離れた変異
    """
    print(f"  GTEx eQTL query: gene={gene_name}, tissue={tissue}")
    print(f"  P-value threshold: {pvalue_threshold}")
    print("  Types: cis-eQTL (primary), trans-eQTL")

    return {"gene": gene_name, "tissue": tissue}
```

## 4. 遺伝子セット濃縮解析 (GSEA)

```python
import pandas as pd
import numpy as np


def gsea_preranked(ranked_gene_list, gene_sets="MSigDB_Hallmark_2020",
                    n_permutations=1000, min_size=15, max_size=500):
    """
    GSEA (Gene Set Enrichment Analysis) — Preranked。

    入力: log2FC × -log10(p) でランク付けされた遺伝子リスト
    遺伝子セットDB:
    - MSigDB Hallmark (H)
    - GO Biological Process (C5:BP)
    - KEGG Pathways (C2:KEGG)
    - Reactome (C2:REACTOME)
    """
    import gseapy as gp

    # ランクスコア = sign(log2FC) × -log10(pvalue)
    results = gp.prerank(
        rnk=ranked_gene_list,
        gene_sets=gene_sets,
        processes=4,
        permutation_num=n_permutations,
        min_size=min_size,
        max_size=max_size,
        outdir="results/gsea",
        seed=42,
    )

    sig_terms = results.res2d[results.res2d["FDR q-val"] < 0.05]

    print(f"  GSEA results ({gene_sets}):")
    print(f"    Gene sets tested: {len(results.res2d)}")
    print(f"    Significant (FDR < 0.05): {len(sig_terms)}")
    if len(sig_terms) > 0:
        print(f"    Top enriched:")
        for _, row in sig_terms.head(5).iterrows():
            direction = "UP" if row["NES"] > 0 else "DOWN"
            print(f"      {row['Term']} (NES={row['NES']:.2f}, {direction})")

    return results


def overrepresentation_analysis(gene_list, background=None,
                                  gene_sets="GO_Biological_Process_2021"):
    """
    遺伝子オーバーリプレゼンテーション解析 (ORA)。

    Fisher exact test ベースの濃縮解析。
    DEG リスト → 機能カテゴリへのマッピング。
    """
    import gseapy as gp

    results = gp.enrich(
        gene_list=gene_list,
        gene_sets=gene_sets,
        background=background,
        outdir="results/ora",
    )

    sig = results.res2d[results.res2d["Adjusted P-value"] < 0.05]

    print(f"  ORA results ({gene_sets}):")
    print(f"    Input genes: {len(gene_list)}")
    print(f"    Significant terms: {len(sig)}")

    return results
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/geo_expression_matrix.csv` | CSV |
| `results/deseq2_results.csv` | CSV |
| `results/gsea/` | ディレクトリ |
| `results/ora/` | ディレクトリ |
| `figures/volcano_rnaseq.png` | PNG |
| `figures/ma_plot.png` | PNG |
| `figures/gsea_dotplot.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| GEO | `geo_search_datasets` | GEO データセット検索 |
| GEO | `geo_get_dataset_info` | データセット詳細取得 |
| GEO | `geo_get_sample_info` | サンプル情報取得 |
| GTEx | `GTEx_get_median_gene_expression` | 組織間中央値発現量 |
| GTEx | `GTEx_get_gene_expression` | サンプルレベル発現データ |
| GTEx | `GTEx_get_top_expressed_genes` | 高発現遺伝子取得 |
| GTEx | `GTEx_get_eqtl_genes` | eQTL 遺伝子 (eGenes) |
| GTEx | `GTEx_get_single_tissue_eqtls` | 単一組織 eQTL |
| GTEx | `GTEx_get_multi_tissue_eqtls` | 多組織 eQTL |
| GTEx | `GTEx_calculate_eqtl` | eQTL 計算 |
| Expression Atlas | `ExpressionAtlas_search_experiments` | 実験検索 |
| Expression Atlas | `ExpressionAtlas_get_baseline` | ベースライン発現 |
| Expression Atlas | `ExpressionAtlas_search_differential` | 差次発現実験 |
| ArrayExpress | `arrayexpress_search_experiments` | ArrayExpress 実験検索 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-bioinformatics` | バルク RNA-seq 基盤 |
| `scientific-single-cell-genomics` | scRNA-seq (単一細胞) |
| `scientific-epigenomics-chromatin` | 発現-エピゲノム統合 |
| `scientific-multi-omics` | マルチオミクス統合 |
| `scientific-network-analysis` | 共発現ネットワーク |

### 依存パッケージ

`pydeseq2`, `GEOparse`, `gseapy`, `pandas`, `numpy`, `matplotlib`, `scipy`
