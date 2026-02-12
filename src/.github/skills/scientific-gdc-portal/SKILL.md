---
name: scientific-gdc-portal
description: |
  NCI Genomic Data Commons ポータルスキル。GDC REST API
  を用いたがんゲノムプロジェクト横断検索・ケースメタデータ・
  体細胞変異 (SSM)・遺伝子発現・ファイル取得。
  ToolUniverse 連携: gdc。
tu_tools:
  - key: gdc
    name: GDC
    description: NCI Genomic Data Commons REST API
---

# Scientific GDC Portal

NCI Genomic Data Commons (GDC) REST API を活用した
がんゲノムプロジェクト横断検索・ケースメタデータ取得・
体細胞変異 (SSM)・遺伝子発現パイプラインを提供する。

## When to Use

- TCGA/TARGET 等のがんゲノムデータを横断検索するとき
- がん種別のケースメタデータを取得するとき
- 特定遺伝子の体細胞変異 (SSM) 頻度を調べるとき
- がんプロジェクトの統計サマリーを取得するとき
- GDC ファイルメタデータを検索してダウンロード URL を取得するとき

---

## Quick Start

## 1. プロジェクト検索・統計

```python
import requests
import pandas as pd

GDC_API = "https://api.gdc.cancer.gov"


def gdc_projects(disease_type=None, limit=50):
    """
    GDC — プロジェクト検索。

    Parameters:
        disease_type: str — 疾患タイプフィルタ
            (例: "Breast Invasive Carcinoma")
        limit: int — 最大結果数
    """
    url = f"{GDC_API}/projects"
    params = {
        "size": limit,
        "fields": ("project_id,name,primary_site,"
                    "disease_type,summary.case_count,"
                    "summary.file_count"),
    }

    if disease_type:
        params["filters"] = (
            '{"op":"=","content":{"field":'
            f'"disease_type","value":"{disease_type}"}}}}'
        )

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for hit in data.get("data", {}).get("hits", []):
        summary = hit.get("summary", {})
        rows.append({
            "project_id": hit.get("project_id", ""),
            "name": hit.get("name", ""),
            "primary_site": "; ".join(
                hit.get("primary_site", [])),
            "disease_type": "; ".join(
                hit.get("disease_type", [])),
            "case_count": summary.get(
                "case_count", 0),
            "file_count": summary.get(
                "file_count", 0),
        })

    df = pd.DataFrame(rows)
    print(f"GDC projects: {len(df)}")
    return df
```

## 2. ケースメタデータ

```python
def gdc_cases(project_id, limit=100):
    """
    GDC — ケースメタデータ取得。

    Parameters:
        project_id: str — プロジェクト ID
            (例: "TCGA-BRCA")
        limit: int — 最大結果数
    """
    url = f"{GDC_API}/cases"
    filters = {
        "op": "=",
        "content": {
            "field": "project.project_id",
            "value": project_id,
        },
    }
    params = {
        "filters": str(filters).replace("'", '"'),
        "fields": ("case_id,submitter_id,"
                    "demographic.gender,"
                    "demographic.race,"
                    "demographic.vital_status,"
                    "diagnoses.primary_diagnosis,"
                    "diagnoses.tumor_stage,"
                    "diagnoses.age_at_diagnosis"),
        "size": limit,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for hit in data.get("data", {}).get("hits", []):
        demo = hit.get("demographic", {}) or {}
        diag = (hit.get("diagnoses", [{}]) or [{}])[0]
        rows.append({
            "case_id": hit.get("case_id", ""),
            "submitter_id": hit.get("submitter_id", ""),
            "gender": demo.get("gender", ""),
            "race": demo.get("race", ""),
            "vital_status": demo.get(
                "vital_status", ""),
            "diagnosis": diag.get(
                "primary_diagnosis", ""),
            "stage": diag.get("tumor_stage", ""),
            "age_at_diagnosis": diag.get(
                "age_at_diagnosis", ""),
        })

    df = pd.DataFrame(rows)
    print(f"GDC cases: {project_id} → {len(df)}")
    return df
```

## 3. 体細胞変異 (SSM) 検索

```python
def gdc_ssm_by_gene(gene_symbol, project_id=None,
                       limit=100):
    """
    GDC — 遺伝子別体細胞変異検索。

    Parameters:
        gene_symbol: str — 遺伝子シンボル (例: "TP53")
        project_id: str — プロジェクト ID フィルタ
        limit: int — 最大結果数
    """
    url = f"{GDC_API}/ssms"
    filters = {
        "op": "and",
        "content": [
            {
                "op": "=",
                "content": {
                    "field":
                        "consequence.transcript."
                        "gene.symbol",
                    "value": gene_symbol,
                },
            }
        ],
    }

    if project_id:
        filters["content"].append({
            "op": "=",
            "content": {
                "field": "cases.project.project_id",
                "value": project_id,
            },
        })

    params = {
        "filters": str(filters).replace("'", '"'),
        "fields": ("ssm_id,genomic_dna_change,"
                    "consequence.transcript.aa_change,"
                    "consequence.transcript."
                    "consequence_type,"
                    "consequence.transcript."
                    "gene.symbol"),
        "size": limit,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for hit in data.get("data", {}).get("hits", []):
        for csq in hit.get("consequence", []):
            tx = csq.get("transcript", {})
            rows.append({
                "ssm_id": hit.get("ssm_id", ""),
                "genomic_change": hit.get(
                    "genomic_dna_change", ""),
                "gene": tx.get("gene", {}).get(
                    "symbol", ""),
                "aa_change": tx.get("aa_change", ""),
                "consequence_type": tx.get(
                    "consequence_type", ""),
            })

    df = pd.DataFrame(rows)
    print(f"GDC SSM: {gene_symbol} → {len(df)} variants")
    return df
```

## 4. GDC 統合パイプライン

```python
def gdc_pipeline(project_id, gene_symbol=None,
                    output_dir="results"):
    """
    GDC 統合パイプライン。

    Parameters:
        project_id: str — プロジェクト ID
        gene_symbol: str — 遺伝子フィルタ
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) プロジェクト情報
    projects = gdc_projects()
    projects.to_csv(output_dir / "gdc_projects.csv",
                    index=False)

    # 2) ケースメタデータ
    cases = gdc_cases(project_id)
    cases.to_csv(output_dir / "gdc_cases.csv",
                 index=False)

    # 3) 体細胞変異
    if gene_symbol:
        ssm = gdc_ssm_by_gene(gene_symbol, project_id)
        ssm.to_csv(output_dir / "gdc_ssm.csv",
                   index=False)

    print(f"GDC pipeline: {project_id} → {output_dir}")
    return {"cases": cases}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `gdc` | GDC | NCI Genomic Data Commons REST API |

## パイプライン統合

```
cancer-genomics → gdc-portal → precision-oncology
  (COSMIC/DepMap)   (GDC API)    (MTB レポート)
        │                │              ↓
icgc-cancer-data ───────┘    variant-interpretation
  (ICGC DCC)                 (ClinVar/ACMG)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/gdc_projects.csv` | プロジェクト一覧 | → cancer-genomics |
| `results/gdc_cases.csv` | ケースメタデータ | → precision-oncology |
| `results/gdc_ssm.csv` | 体細胞変異 | → variant-interpretation |
