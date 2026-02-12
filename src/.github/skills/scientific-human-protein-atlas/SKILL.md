---
name: scientific-human-protein-atlas
description: |
  Human Protein Atlas (HPA) 統合スキル。組織/細胞タンパク質発現、
  がん予後バイオマーカー、RNA 発現プロファイル、細胞内局在、
  タンパク質相互作用の包括的検索・解析パイプライン。
tu_tools:
  - key: hpa
    name: Human Protein Atlas
    description: 組織/細胞タンパク質発現・RNA 発現・がん予後
---

# Scientific Human Protein Atlas

HPA REST API を活用した組織・細胞レベルの
タンパク質発現プロファイリングパイプラインを提供する。

## When to Use

- 遺伝子/タンパク質の組織発現パターンを調べるとき
- がん予後バイオマーカー候補を評価するとき
- 細胞内局在 (subcellular localization) を確認するとき
- 細胞株間の発現比較を行うとき
- RNA 発現データ (HPA/GTEx/FANTOM5) を統合するとき

---

## Quick Start

## 1. HPA 遺伝子基本情報取得

```python
import requests
import pandas as pd

HPA_API = "https://www.proteinatlas.org/api"


def get_hpa_gene_info(ensembl_id):
    """
    HPA 遺伝子基本情報取得。

    Parameters:
        ensembl_id: str — Ensembl gene ID (e.g., "ENSG00000141510")

    ToolUniverse:
        HPA_get_gene_basic_info_by_ensembl_id(ensembl_id=ensembl_id)
        HPA_get_comprehensive_gene_details_by_ensembl_id(ensembl_id=ensembl_id)
    """
    url = f"https://www.proteinatlas.org/{ensembl_id}.json"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    info = {
        "ensembl_id": ensembl_id,
        "gene_name": data.get("Gene", ""),
        "gene_description": data.get("Gene description", ""),
        "uniprot_id": data.get("Uniprot", []),
        "chromosome": data.get("Chromosome", ""),
        "protein_class": data.get("Protein class", []),
        "evidence": data.get("Evidence", ""),
    }

    print(f"HPA gene: {info['gene_name']} ({ensembl_id})")
    return info, data
```

## 2. 組織 RNA 発現プロファイル

```python
def get_tissue_rna_expression(gene_name):
    """
    HPA 組織別 RNA 発現データ取得。

    ToolUniverse:
        HPA_get_rna_expression_by_source(gene=gene_name, source="HPA")
        HPA_get_rna_expression_in_specific_tissues(gene=gene_name, tissues=tissues)
    """
    url = f"https://www.proteinatlas.org/{gene_name}.json"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    rna_data = data.get("RNA tissue specific nTPM", [])
    results = []
    for entry in rna_data:
        results.append({
            "tissue": entry.get("Tissue", ""),
            "cell_type": entry.get("Cell type", ""),
            "ntpm": float(entry.get("nTPM", 0)),
            "detection": entry.get("Detection", ""),
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("ntpm", ascending=False)

    print(f"HPA RNA expression '{gene_name}': {len(df)} tissue entries")
    return df
```

## 3. がん予後バイオマーカー解析

```python
def get_cancer_prognostics(gene_name):
    """
    HPA がん予後データ取得。

    ToolUniverse:
        HPA_get_cancer_prognostics_by_gene(gene=gene_name)
    """
    url = f"https://www.proteinatlas.org/{gene_name}.json"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    prognostics = data.get("Pathology prognostics", [])
    results = []
    for entry in prognostics:
        results.append({
            "cancer_type": entry.get("Cancer type", ""),
            "prognostic_type": entry.get("Prognostic type", ""),
            "is_prognostic": entry.get("Is prognostic", False),
            "p_value": float(entry.get("p-value", 1.0)),
            "high_expression_favorable": entry.get(
                "High expression is favorable", None
            ),
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("p_value")
        significant = df[df["p_value"] < 0.05]
        print(f"HPA cancer prognostics '{gene_name}': "
              f"{len(significant)}/{len(df)} significant")
    else:
        print(f"HPA cancer prognostics '{gene_name}': no data")
    return df
```

## 4. 細胞内局在

```python
def get_subcellular_location(gene_name):
    """
    HPA 細胞内局在データ取得。

    ToolUniverse:
        HPA_get_subcellular_location(gene=gene_name)
    """
    url = f"https://www.proteinatlas.org/{gene_name}.json"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    sc = data.get("Subcellular location", [])
    results = []
    for entry in sc:
        results.append({
            "location": entry.get("Location", ""),
            "reliability": entry.get("Reliability", ""),
            "enhanced": entry.get("Enhanced", False),
            "supported": entry.get("Supported", False),
            "cell_lines": entry.get("Cell lines", []),
        })

    df = pd.DataFrame(results)
    print(f"HPA subcellular '{gene_name}': {len(df)} locations")
    return df
```

## 5. タンパク質相互作用ネットワーク (HPA)

```python
def get_hpa_protein_interactions(gene_name):
    """
    HPA タンパク質相互作用データ取得。

    ToolUniverse:
        HPA_get_protein_interactions_by_gene(gene=gene_name)
        HPA_get_biological_processes_by_gene(gene=gene_name)
        HPA_get_contextual_biological_process_analysis(gene=gene_name)
    """
    url = f"https://www.proteinatlas.org/{gene_name}.json"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    interactions = data.get("Protein interaction partners", [])
    results = []
    for partner in interactions:
        results.append({
            "partner_gene": partner.get("Gene", ""),
            "partner_ensembl": partner.get("Ensembl", ""),
            "confidence": partner.get("Confidence", ""),
            "source": partner.get("Source", ""),
        })

    df = pd.DataFrame(results)
    print(f"HPA interactions '{gene_name}': {len(df)} partners")
    return df
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/hpa_gene_info.json` | JSON |
| `results/hpa_tissue_expression.csv` | CSV |
| `results/hpa_cancer_prognostics.csv` | CSV |
| `results/hpa_subcellular.csv` | CSV |
| `results/hpa_interactions.csv` | CSV |

### 利用可能ツール

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| HPA | `HPA_generic_search` | 汎用検索 |
| HPA | `HPA_get_gene_basic_info_by_ensembl_id` | 遺伝子基本情報 |
| HPA | `HPA_get_comprehensive_gene_details_by_ensembl_id` | 包括的詳細 |
| HPA | `HPA_get_rna_expression_by_source` | RNA 発現 |
| HPA | `HPA_get_rna_expression_in_specific_tissues` | 組織別発現 |
| HPA | `HPA_get_cancer_prognostics_by_gene` | がん予後 |
| HPA | `HPA_get_subcellular_location` | 細胞内局在 |
| HPA | `HPA_get_protein_interactions_by_gene` | PPI |
| HPA | `HPA_get_biological_processes_by_gene` | 生物学的プロセス |
| HPA | `HPA_get_contextual_biological_process_analysis` | プロセス解析 |
| HPA | `HPA_get_disease_expression_by_gene_tissue_disease` | 疾患発現 |
| HPA | `HPA_get_comparative_expression_by_gene_and_cellline` | 細胞株比較 |
| HPA | `HPA_get_gene_tsv_data_by_ensembl_id` | TSV データ |
| HPA | `HPA_search_genes_by_query` | 遺伝子検索 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-gene-expression-transcriptomics` | GEO/GTEx 発現解析 |
| `scientific-proteomics-mass-spectrometry` | プロテオミクス |
| `scientific-cancer-genomics` | がんゲノミクス |
| `scientific-protein-interaction-network` | PPI ネットワーク |
| `scientific-pathway-enrichment` | パスウェイ濃縮 |

### 依存パッケージ

`requests`, `pandas`
