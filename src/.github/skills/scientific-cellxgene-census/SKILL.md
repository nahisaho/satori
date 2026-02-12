---
name: scientific-cellxgene-census
description: |
  CELLxGENE Census 大規模シングルセルアトラススキル。
  CZ CELLxGENE Census API によるヒト/マウス全アトラスの
  メタデータ検索・遺伝子発現クエリ・セルタイプ分布解析・
  データセット横断統合パイプライン。
  ToolUniverse 連携: cellxgene_census。
tu_tools:
  - key: cellxgene_census
    name: CELLxGENE Census
    description: 大規模シングルセルアトラスデータアクセス API
---

# Scientific CELLxGENE Census

CZ CELLxGENE Census API を活用した大規模シングルセルアトラスの
メタデータ検索・遺伝子発現クエリ・セルタイプ分布解析・
データセット横断統合パイプラインを提供する。

## When to Use

- 数千万細胞規模のシングルセルアトラスから特定組織/疾患のデータを抽出するとき
- 組織横断的なセルタイプ分布を比較するとき
- 特定遺伝子の全アトラスにわたる発現パターンを検索するとき
- Census データセットをメタデータベースでフィルタリングするとき
- AnnData/Sparse 行列として大規模データを効率的に取得するとき

---

## Quick Start

## 1. Census メタデータ検索

```python
import cellxgene_census
import pandas as pd


def census_datasets(organism="Homo sapiens",
                      tissue=None, disease=None):
    """
    CELLxGENE Census — データセットメタデータ検索。

    Parameters:
        organism: str — 生物種
        tissue: str — 組織フィルタ
        disease: str — 疾患フィルタ
    """
    with cellxgene_census.open_soma() as census:
        datasets = census["census_info"]["datasets"].read(
        ).concat().to_pandas()

    if tissue:
        datasets = datasets[
            datasets["dataset_title"].str.contains(
                tissue, case=False, na=False)]
    if disease:
        datasets = datasets[
            datasets["dataset_title"].str.contains(
                disease, case=False, na=False)]

    print(f"Census datasets: {len(datasets)} matched")
    return datasets


def census_cell_metadata(organism="Homo sapiens",
                           tissue=None,
                           cell_type=None,
                           max_cells=10000):
    """
    CELLxGENE Census — セルメタデータ取得。

    Parameters:
        organism: str — 生物種
        tissue: str — 組織フィルタ
        cell_type: str — セルタイプフィルタ
        max_cells: int — 最大セル数
    """
    obs_filters = []
    if tissue:
        obs_filters.append(
            f"tissue_general == '{tissue}'")
    if cell_type:
        obs_filters.append(
            f"cell_type == '{cell_type}'")
    value_filter = " and ".join(obs_filters) \
        if obs_filters else None

    with cellxgene_census.open_soma() as census:
        obs_df = cellxgene_census.get_obs(
            census, organism,
            value_filter=value_filter,
            column_names=[
                "cell_type", "tissue_general",
                "disease", "sex", "development_stage",
                "dataset_id", "assay"],
        )

    df = obs_df.head(max_cells)
    print(f"Census cells: {len(df)} retrieved "
          f"(filter: {value_filter})")
    return df
```

## 2. 遺伝子発現クエリ

```python
def census_gene_expression(organism="Homo sapiens",
                             gene_symbols=None,
                             tissue=None,
                             max_cells=5000):
    """
    CELLxGENE Census — 遺伝子発現データ取得。

    Parameters:
        organism: str — 生物種
        gene_symbols: list[str] — 遺伝子シンボルリスト
        tissue: str — 組織フィルタ
        max_cells: int — 最大セル数
    """
    obs_filter = (f"tissue_general == '{tissue}'"
                  if tissue else None)
    var_filter = None
    if gene_symbols:
        genes_str = "', '".join(gene_symbols)
        var_filter = f"feature_name in ['{genes_str}']"

    with cellxgene_census.open_soma() as census:
        adata = cellxgene_census.get_anndata(
            census, organism,
            obs_value_filter=obs_filter,
            var_value_filter=var_filter,
            obs_column_names=[
                "cell_type", "tissue_general",
                "disease"],
        )

    if max_cells and len(adata) > max_cells:
        import numpy as np
        idx = np.random.choice(
            len(adata), max_cells, replace=False)
        adata = adata[idx]

    print(f"Census expression: {adata.shape[0]} cells × "
          f"{adata.shape[1]} genes")
    return adata
```

## 3. セルタイプ分布解析

```python
def celltype_distribution(organism="Homo sapiens",
                            tissue=None):
    """
    組織別セルタイプ分布解析。

    Parameters:
        organism: str — 生物種
        tissue: str — 組織フィルタ
    """
    obs_filter = (f"tissue_general == '{tissue}'"
                  if tissue else None)

    with cellxgene_census.open_soma() as census:
        obs_df = cellxgene_census.get_obs(
            census, organism,
            value_filter=obs_filter,
            column_names=["cell_type",
                          "tissue_general"],
        )

    # セルタイプカウント
    ct_counts = obs_df["cell_type"].value_counts()
    ct_df = pd.DataFrame({
        "cell_type": ct_counts.index,
        "count": ct_counts.values,
        "fraction": ct_counts.values / ct_counts.sum(),
    })

    print(f"Cell types: {len(ct_df)} types, "
          f"total {ct_counts.sum()} cells")
    return ct_df
```

## 4. Census 統合パイプライン

```python
def census_pipeline(organism="Homo sapiens",
                      tissue=None,
                      gene_symbols=None,
                      output_dir="results"):
    """
    CELLxGENE Census 統合パイプライン。

    Parameters:
        organism: str — 生物種
        tissue: str — 組織
        gene_symbols: list[str] — 対象遺伝子
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) データセットメタデータ
    datasets = census_datasets(organism, tissue)
    datasets.to_csv(output_dir / "census_datasets.csv",
                    index=False)

    # 2) セルタイプ分布
    ct_dist = celltype_distribution(organism, tissue)
    ct_dist.to_csv(
        output_dir / "celltype_distribution.csv",
        index=False)

    # 3) 遺伝子発現 (指定時)
    adata = None
    if gene_symbols:
        adata = census_gene_expression(
            organism, gene_symbols, tissue)
        adata.write_h5ad(
            output_dir / "census_expression.h5ad")

    print(f"Census pipeline → {output_dir}")
    return {
        "datasets": datasets,
        "celltype_dist": ct_dist,
        "adata": adata,
    }
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `cellxgene_census` | CELLxGENE Census | 大規模シングルセルアトラス API |

## パイプライン統合

```
human-cell-atlas → cellxgene-census → single-cell-genomics
  (HCA Portal)     (Census 大規模)     (Scanpy 解析)
       │                  │                  ↓
  spatial-multiomics ────┘      scvi-integration
    (空間統合)                   (scVI/scANVI)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/census_datasets.csv` | データセット一覧 | → human-cell-atlas |
| `results/celltype_distribution.csv` | セルタイプ分布 | → single-cell-genomics |
| `results/census_expression.h5ad` | 発現行列 (AnnData) | → scvi-integration |
