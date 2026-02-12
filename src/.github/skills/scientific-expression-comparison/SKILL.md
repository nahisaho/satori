---
name: scientific-expression-comparison
description: |
  Expression Atlas / GTEx / HPA 統合発現比較スキル。EBI Expression Atlas
  ベースライン/差次的発現検索、実験アクセション取得、組織間・条件間
  発現比較、マルチソース統合発現プロファイリングパイプライン。
---

# Scientific Expression Comparison

EBI Expression Atlas API を中核として GTEx/HPA データと統合した
遺伝子発現比較パイプラインを提供する。ベースライン発現・差次的発現・
実験メタデータを横断的に検索・比較。

## When to Use

- EBI Expression Atlas で遺伝子のベースライン発現パターンを調べるとき
- 疾患 vs 正常の差次的発現データを検索するとき
- 複数の組織/細胞型にわたる発現比較を行うとき
- Expression Atlas の実験メタデータを取得するとき
- GTEx/HPA の発現データと統合して比較分析するとき

---

## Quick Start

## 1. ベースライン発現検索

```python
import requests
import pandas as pd

ATLAS_API = "https://www.ebi.ac.uk/gxa/json"
ATLAS_REST = "https://www.ebi.ac.uk/gxa"


def get_baseline_expression(gene, species="homo sapiens"):
    """
    Expression Atlas ベースライン発現プロファイル取得。

    Parameters:
        gene: str — 遺伝子シンボルまたは Ensembl ID
        species: str — 生物種

    ToolUniverse:
        ExpressionAtlas_get_baseline(gene=gene, species=species)
    """
    url = f"{ATLAS_REST}/json/baseline_expression"
    params = {"gene": gene, "species": species}
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    profiles = data.get("profiles", {}).get("rows", [])
    rows = []
    for profile in profiles:
        gene_name = profile.get("name", gene)
        for exp in profile.get("expressions", []):
            rows.append({
                "gene": gene_name,
                "factor_value": exp.get("factorValue", ""),
                "expression_level": exp.get("value"),
                "unit": "TPM",
            })

    df = pd.DataFrame(rows)
    print(f"Baseline expression for {gene}: {len(df)} tissue/cell profiles")
    return df
```

## 2. 差次的発現検索

```python
def search_differential_expression(gene, condition=None, species="homo sapiens"):
    """
    差次的発現実験の検索。

    Parameters:
        gene: str — 遺伝子シンボルまたは Ensembl ID
        condition: str — 条件 (例: "cancer", "inflammation")
        species: str — 生物種

    ToolUniverse:
        ExpressionAtlas_search_differential(
            gene=gene, condition=condition, species=species
        )
    """
    url = f"{ATLAS_REST}/json/search"
    params = {"geneQuery": gene, "species": species}
    if condition:
        params["conditionQuery"] = condition

    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    results = data.get("results", [])
    rows = []
    for r in results:
        rows.append({
            "experiment_accession": r.get("experimentAccession"),
            "experiment_description": r.get("experimentDescription", "")[:200],
            "experiment_type": r.get("experimentType"),
            "species": r.get("species"),
            "contrast": r.get("contrastId", ""),
            "log2_fold_change": r.get("foldChange"),
            "p_value": r.get("pValue"),
        })

    df = pd.DataFrame(rows)
    print(f"Differential expression for {gene}: {len(df)} contrasts found")
    return df
```

## 3. 実験メタデータ取得

```python
def get_experiment_details(accession):
    """
    Expression Atlas 実験詳細取得。

    Parameters:
        accession: str — 実験 ID (例: "E-MTAB-5214")

    ToolUniverse:
        ExpressionAtlas_get_experiment(accession=accession)
    """
    url = f"{ATLAS_REST}/json/experiments/{accession}"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    experiment = data.get("experiment", {})
    info = {
        "accession": experiment.get("accession"),
        "description": experiment.get("description"),
        "type": experiment.get("type"),
        "species": experiment.get("species", []),
        "pubmed_ids": experiment.get("pubmedIds", []),
        "n_assays": experiment.get("numberOfAssays"),
        "n_contrasts": experiment.get("numberOfContrasts"),
        "last_updated": experiment.get("lastUpdate"),
    }

    print(f"Experiment: {info['accession']} — {info['description'][:100]}")
    return info
```

## 4. 実験検索

```python
def search_experiments(gene=None, condition=None, species="homo sapiens"):
    """
    Expression Atlas 実験検索。

    Parameters:
        gene: str — 遺伝子クエリ
        condition: str — 条件クエリ
        species: str — 生物種

    ToolUniverse:
        ExpressionAtlas_search_experiments(
            gene=gene, condition=condition, species=species
        )
    """
    url = f"{ATLAS_REST}/json/search"
    params = {"species": species}
    if gene:
        params["geneQuery"] = gene
    if condition:
        params["conditionQuery"] = condition

    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    experiments = data.get("matchingExperiments", [])
    rows = []
    for e in experiments:
        rows.append({
            "accession": e.get("experimentAccession"),
            "description": e.get("experimentDescription", "")[:200],
            "type": e.get("experimentType"),
            "species": e.get("species"),
            "n_assays": e.get("numberOfAssays"),
        })

    df = pd.DataFrame(rows)
    print(f"Experiments found: {len(df)}")
    return df
```

## 5. 組織横断発現比較

```python
def cross_tissue_comparison(genes, species="homo sapiens"):
    """
    複数遺伝子の組織横断発現比較。

    Parameters:
        genes: list — 遺伝子リスト
        species: str — 生物種

    Returns:
        DataFrame — genes × tissues 発現マトリクス
    """
    all_data = []
    for gene in genes:
        df = get_baseline_expression(gene, species)
        if not df.empty:
            df["gene_query"] = gene
            all_data.append(df)

    if not all_data:
        print("No expression data found")
        return pd.DataFrame()

    combined = pd.concat(all_data, ignore_index=True)

    # ピボットテーブル (genes × tissues)
    matrix = combined.pivot_table(
        index="gene",
        columns="factor_value",
        values="expression_level",
        aggfunc="mean",
    )

    print(f"Expression matrix: {matrix.shape[0]} genes × "
          f"{matrix.shape[1]} tissues/conditions")
    return matrix
```

## 6. 発現ヒートマップ

```python
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def plot_expression_heatmap(matrix, title="Gene Expression Comparison",
                            figsize=(14, 8), save_path=None):
    """
    発現マトリクスのヒートマップ描画。

    Parameters:
        matrix: DataFrame — cross_tissue_comparison の出力
        title: str — 図タイトル
        figsize: tuple — 図サイズ
        save_path: str — 保存パス
    """
    log_matrix = np.log2(matrix.fillna(0) + 1)

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        log_matrix,
        cmap="viridis",
        xticklabels=True,
        yticklabels=True,
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel("Tissue / Condition")
    ax.set_ylabel("Gene")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.show()
```

---

## パイプライン統合

```
gene-expression-transcriptomics → expression-comparison → multi-omics
  (GEO/GTEx/DESeq2)                (Atlas 発現比較)        (統合解析)
        │                                │                     ↓
human-protein-atlas ────────────┘        │            pathway-enrichment
  (HPA 組織/がん発現)                     ↓              (KEGG/GO)
                              ontology-enrichment
                              (EFO 形質マッピング)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/baseline_expression.csv` | ベースライン発現 | → gene-expression |
| `results/differential_expression.csv` | 差次的発現 | → pathway-enrichment |
| `results/expression_matrix.csv` | 発現マトリクス | → multi-omics |
| `figures/expression_heatmap.png` | ヒートマップ | → publication-figures |

## 利用可能ツール (ToolUniverse SMCP)

| ツール名 | 用途 |
|---------|------|
| `ExpressionAtlas_get_baseline` | ベースライン発現 |
| `ExpressionAtlas_search_differential` | 差次的発現検索 |
| `ExpressionAtlas_search_experiments` | 実験検索 |
| `ExpressionAtlas_get_experiment` | 実験詳細 |
