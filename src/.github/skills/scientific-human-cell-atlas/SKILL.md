---
name: scientific-human-cell-atlas
description: |
  Human Cell Atlas (HCA) データポータルスキル。HCA Data Portal API
  プロジェクト検索・ファイルダウンロード・CELLxGENE Census 統合・
  細胞型アノテーション・アトラス構築。ToolUniverse 連携: hca_tools。
tu_tools:
  - key: hca_tools
    name: Human Cell Atlas Tools
    description: HCA データポータル プロジェクト・バンドル・ファイル検索
  - key: cellxgene_census
    name: CELLxGENE Census
    description: 大規模シングルセルアトラスデータアクセス API
---

# Scientific Human Cell Atlas

HCA Data Portal / CELLxGENE Census を活用した大規模シングルセル
アトラスデータアクセス・解析パイプラインを提供する。

## When to Use

- HCA Data Portal からプロジェクト・実験データを検索するとき
- CELLxGENE Census で大規模シングルセルアトラスを照会するとき
- 特定組織/疾患の細胞型構成を調べるとき
- 複数 HCA プロジェクト間で細胞型を比較するとき
- シングルセルアトラスのリファレンスマッピングを行うとき
- 希少細胞型の発見・アノテーションを実施するとき

---

## Quick Start

## 1. HCA Data Portal プロジェクト検索

```python
import requests
import pandas as pd
import json

HCA_BASE = "https://service.azul.data.humancellatlas.org"
HCA_CATALOG = "dcp44"


def hca_search_projects(keyword=None, organ=None, disease=None,
                         species="Homo sapiens", limit=25):
    """
    HCA Data Portal — プロジェクト検索。

    Parameters:
        keyword: str — キーワード検索
        organ: str — 臓器 (例: "lung", "heart")
        disease: str — 疾患 (例: "COVID-19")
        species: str — 生物種
        limit: int — 最大結果数
    """
    url = f"{HCA_BASE}/index/projects"
    params = {"catalog": HCA_CATALOG, "size": limit}

    filters = {}
    if organ:
        filters["organ"] = {"is": [organ]}
    if disease:
        filters["disease"] = {"is": [disease]}
    if species:
        filters["genusSpecies"] = {"is": [species]}
    if keyword:
        params["q"] = keyword
    if filters:
        params["filters"] = json.dumps(filters)

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    projects = []
    for hit in data.get("hits", []):
        proj = hit.get("projects", [{}])[0]
        samples = hit.get("samples", [{}])[0]
        protocols = hit.get("protocols", [{}])[0]
        projects.append({
            "project_id": proj.get("projectId", ""),
            "title": proj.get("projectTitle", ""),
            "organ": ", ".join(samples.get("organ", [])),
            "disease": ", ".join(samples.get("disease", [])),
            "species": ", ".join(samples.get("genusSpecies", [])),
            "cell_count": hit.get("cellSuspensions", [{}])[0].get(
                "totalCells", 0),
            "library_method": ", ".join(protocols.get(
                "libraryConstructionApproach", [])),
            "donor_count": samples.get("donorCount", 0),
        })

    df = pd.DataFrame(projects)
    print(f"HCA: {len(df)} projects found")
    return df
```

## 2. HCA ファイル取得

```python
def hca_get_project_files(project_id, file_format=None):
    """
    HCA — プロジェクトのファイル一覧取得。

    Parameters:
        project_id: str — プロジェクト UUID
        file_format: str — ファイル形式 (例: "h5ad", "loom", "csv")
    """
    url = f"{HCA_BASE}/index/files"
    filters = {"projectId": {"is": [project_id]}}
    if file_format:
        filters["fileFormat"] = {"is": [file_format]}

    params = {
        "catalog": HCA_CATALOG,
        "filters": json.dumps(filters),
        "size": 100,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    files = []
    for hit in data.get("hits", []):
        for f in hit.get("files", []):
            files.append({
                "file_id": f.get("uuid", ""),
                "name": f.get("name", ""),
                "format": f.get("format", ""),
                "size_bytes": f.get("size", 0),
                "url": f.get("url", ""),
            })

    df = pd.DataFrame(files)
    print(f"HCA files ({project_id[:8]}): {len(df)} files")
    return df
```

## 3. CELLxGENE Census アトラスクエリ

```python
def cellxgene_census_query(organism="homo_sapiens", tissue=None,
                            disease=None, cell_type=None,
                            gene_list=None, max_cells=50000):
    """
    CELLxGENE Census — 大規模シングルセルアトラスクエリ。

    Parameters:
        organism: str — 生物種
        tissue: str — 組織
        disease: str — 疾患
        cell_type: str — 細胞型
        gene_list: list[str] — 取得遺伝子リスト
        max_cells: int — 最大細胞数
    """
    import cellxgene_census

    census = cellxgene_census.open_soma()

    # 観察フィルタ構築
    obs_filters = []
    if tissue:
        obs_filters.append(f"tissue == '{tissue}'")
    if disease:
        obs_filters.append(f"disease == '{disease}'")
    if cell_type:
        obs_filters.append(f"cell_type == '{cell_type}'")

    obs_filter = " and ".join(obs_filters) if obs_filters else None

    # 遺伝子フィルタ
    var_filter = None
    if gene_list:
        genes_str = "', '".join(gene_list)
        var_filter = f"feature_name in ['{genes_str}']"

    adata = cellxgene_census.get_anndata(
        census,
        organism=organism,
        obs_value_filter=obs_filter,
        var_value_filter=var_filter,
        obs_column_names=[
            "cell_type", "tissue", "disease",
            "donor_id", "dataset_id", "assay",
        ],
    )

    if adata.n_obs > max_cells:
        import numpy as np
        idx = np.random.choice(adata.n_obs, max_cells, replace=False)
        adata = adata[idx].copy()

    census.close()
    print(f"CELLxGENE Census: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata
```

## 4. 細胞型構成解析

```python
import scanpy as sc


def cell_type_composition(adata, groupby="tissue", cell_type_col="cell_type"):
    """
    細胞型構成の定量比較。

    Parameters:
        adata: AnnData — シングルセルデータ
        groupby: str — グループ変数
        cell_type_col: str — 細胞型カラム名
    """
    # 構成比計算
    composition = (
        adata.obs.groupby([groupby, cell_type_col])
        .size()
        .unstack(fill_value=0)
    )
    composition_pct = composition.div(composition.sum(axis=1), axis=0) * 100

    print(f"Cell type composition: {composition.shape[0]} groups × "
          f"{composition.shape[1]} cell types")
    return composition_pct
```

## 5. HCA 統合パイプライン

```python
def hca_atlas_pipeline(organ, disease=None, output_dir="results"):
    """
    HCA + CELLxGENE 統合アトラスパイプライン。

    Parameters:
        organ: str — 対象臓器
        disease: str — 対象疾患
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) HCA プロジェクト検索
    projects = hca_search_projects(organ=organ, disease=disease)
    projects.to_csv(output_dir / "hca_projects.csv", index=False)

    # 2) CELLxGENE Census クエリ
    adata = cellxgene_census_query(tissue=organ, disease=disease)

    # 3) 前処理
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # 4) 細胞型構成
    composition = cell_type_composition(adata)
    composition.to_csv(output_dir / "cell_type_composition.csv")

    # 5) 保存
    adata.write(output_dir / "hca_atlas.h5ad")

    print(f"HCA atlas pipeline: {output_dir}")
    return {"projects": projects, "adata": adata, "composition": composition}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `hca_tools` | HCA Tools | プロジェクト検索・ファイルダウンロード |

## パイプライン統合

```
single-cell-genomics → human-cell-atlas → scvi-integration
  (scanpy 標準)         (HCA/CELLxGENE)   (scVI 統合)
       │                      │                ↓
  spatial-transcriptomics ───┘         cell-type-annotation
  (Visium/MERFISH)       │              (リファレンスマッピング)
                         ↓
                   gpu-singlecell
                   (大規模処理)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/hca_projects.csv` | HCA プロジェクト一覧 | → single-cell-genomics |
| `results/hca_atlas.h5ad` | アトラス AnnData | → scvi-integration |
| `results/cell_type_composition.csv` | 細胞型構成比 | → spatial-transcriptomics |
