---
name: scientific-icgc-cancer-data
description: |
  ICGC がんゲノムデータスキル。ICGC ARGO DCC API および
  レガシー API による国際がんゲノムデータ検索・ドナー/
  検体/変異解析。直接 API (ToolUniverse 非連携)。
tu_tools: []
---

# Scientific ICGC Cancer Data

ICGC (International Cancer Genome Consortium) ARGO DCC API を
活用した国際がんゲノムデータ検索・変異統計・がん種横断解析
パイプラインを提供する。

## When to Use

- 国際がんゲノムプロジェクトのデータを検索するとき
- がん種ごとの体細胞変異プロファイルを調べるとき
- ドナー・検体・変異の統計情報を取得するとき
- がんゲノムの変異シグネチャを分析するとき
- PCAWG (Pan-Cancer Analysis of Whole Genomes) データを活用するとき
- がん遺伝子変異の国際比較データが必要なとき

---

## Quick Start

## 1. ICGC プロジェクト・ドナー検索

```python
import requests
import pandas as pd

ICGC_BASE = "https://dcc.icgc.org/api/v1"


def icgc_search_projects(query=None, limit=50):
    """
    ICGC — がんゲノムプロジェクト検索。

    Parameters:
        query: str — 検索キーワード (例: "lung", "BRCA")
        limit: int — 最大結果数
    """
    url = f"{ICGC_BASE}/projects"
    params = {"size": limit, "from": 1}
    if query:
        params["filters"] = (
            f'{{"project":{{"primarySite":'
            f'{{"is":["{query}"]}}}}}}'
        )

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for hit in data.get("hits", []):
        results.append({
            "project_id": hit.get("id", ""),
            "project_name": hit.get("name", ""),
            "primary_site": hit.get("primarySite", ""),
            "tumour_type": hit.get("tumourType", ""),
            "tumour_subtype": hit.get("tumourSubtype", ""),
            "primary_country": "; ".join(
                hit.get("primaryCountries", [])),
            "total_donors": hit.get("totalDonorCount", 0),
            "ssm_count": hit.get("ssmCount", 0),
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("total_donors", ascending=False)

    total = data.get("pagination", {}).get("total", 0)
    print(f"ICGC projects: {len(df)}/{total} "
          f"(query='{query}')")
    return df


def icgc_search_donors(project_id, limit=100):
    """
    ICGC — プロジェクト内ドナー検索。

    Parameters:
        project_id: str — プロジェクト ID (例: "BRCA-US")
        limit: int — 最大結果数
    """
    url = f"{ICGC_BASE}/donors"
    params = {
        "size": limit,
        "filters": (
            f'{{"donor":{{"projectId":'
            f'{{"is":["{project_id}"]}}}}}}'
        ),
    }

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for hit in data.get("hits", []):
        results.append({
            "donor_id": hit.get("id", ""),
            "project_id": project_id,
            "primary_site": hit.get("primarySite", ""),
            "gender": hit.get("gender", ""),
            "vital_status": hit.get("vitalStatus", ""),
            "age_at_diagnosis": hit.get("ageAtDiagnosis"),
            "disease_status": hit.get(
                "diseaseStatusLastFollowup", ""),
            "ssm_count": hit.get("ssmCount", 0),
        })

    df = pd.DataFrame(results)
    total = data.get("pagination", {}).get("total", 0)
    print(f"ICGC donors: {len(df)}/{total} "
          f"(project={project_id})")
    return df
```

## 2. 体細胞変異 (SSM) 検索

```python
def icgc_search_mutations(gene_symbol=None,
                             project_id=None,
                             consequence_type=None,
                             limit=100):
    """
    ICGC — 体細胞変異 (Simple Somatic Mutation) 検索。

    Parameters:
        gene_symbol: str — 遺伝子シンボル (例: "TP53")
        project_id: str — プロジェクト ID
        consequence_type: str — 変異タイプ
            (例: "missense_variant")
        limit: int — 最大結果数
    """
    url = f"{ICGC_BASE}/mutations"
    filters = {}

    if gene_symbol:
        filters["gene"] = {"symbol": {"is": [gene_symbol]}}
    if project_id:
        filters["donor"] = {"projectId": {"is": [project_id]}}
    if consequence_type:
        filters["mutation"] = {
            "consequenceType": {"is": [consequence_type]}
        }

    import json
    params = {
        "size": limit,
        "filters": json.dumps(filters) if filters else "{}",
    }

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for hit in data.get("hits", []):
        # 主要な consequence 取得
        consequences = hit.get("consequences", [])
        top_cons = consequences[0] if consequences else {}

        results.append({
            "mutation_id": hit.get("id", ""),
            "chromosome": hit.get("chromosome", ""),
            "start": hit.get("start"),
            "end": hit.get("end"),
            "mutation": hit.get("mutation", ""),
            "type": hit.get("type", ""),
            "gene_symbol": top_cons.get("geneSymbol", ""),
            "consequence_type": top_cons.get("type", ""),
            "aa_mutation": top_cons.get("aaMutation", ""),
            "affected_donors": hit.get(
                "affectedDonorCountTotal", 0),
            "affected_projects": hit.get(
                "affectedProjectCount", 0),
            "functional_impact": hit.get(
                "functionalImpact", ""),
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("affected_donors",
                            ascending=False)

    total = data.get("pagination", {}).get("total", 0)
    print(f"ICGC mutations: {len(df)}/{total} "
          f"(gene={gene_symbol}, project={project_id})")
    return df
```

## 3. がん種統計・変異サマリー

```python
def icgc_cancer_stats(project_id=None):
    """
    ICGC — がん種統計サマリー。

    Parameters:
        project_id: str — プロジェクト ID (None で全体統計)
    """
    if project_id:
        url = f"{ICGC_BASE}/projects/{project_id}"
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        stats = {
            "project_id": project_id,
            "project_name": data.get("name", ""),
            "primary_site": data.get("primarySite", ""),
            "total_donors": data.get("totalDonorCount", 0),
            "total_specimens": data.get(
                "totalSpecimenCount", 0),
            "ssm_count": data.get("ssmCount", 0),
            "repository": "; ".join(
                data.get("repository", [])),
        }
        print(f"ICGC stats: {project_id} — "
              f"{stats['total_donors']} donors, "
              f"{stats['ssm_count']} mutations")
        return stats
    else:
        # 全プロジェクト概要
        projects = icgc_search_projects(limit=200)
        summary = {
            "total_projects": len(projects),
            "total_donors": projects[
                "total_donors"].sum(),
            "total_ssm": projects["ssm_count"].sum(),
            "top_sites": projects.groupby(
                "primary_site")["total_donors"].sum(
                ).sort_values(ascending=False).head(10
                ).to_dict(),
        }
        print(f"ICGC summary: {summary['total_projects']} "
              f"projects, {summary['total_donors']} donors")
        return summary


def icgc_gene_mutation_frequency(gene_symbol, top_n=20):
    """
    ICGC — 遺伝子別がん種変異頻度。

    Parameters:
        gene_symbol: str — 遺伝子シンボル
        top_n: int — 上位がん種数
    """
    mutations = icgc_search_mutations(
        gene_symbol=gene_symbol, limit=500)

    if mutations.empty:
        return pd.DataFrame()

    # プロジェクト別集計
    freq = mutations.groupby("gene_symbol").agg(
        total_mutations=("mutation_id", "count"),
        total_affected_donors=("affected_donors", "sum"),
        mutation_types=("consequence_type",
                        lambda x: "; ".join(x.unique()[:5])),
    ).reset_index()

    print(f"ICGC gene frequency: {gene_symbol} — "
          f"{len(freq)} entries")
    return freq
```

## 4. ICGC 統合パイプライン

```python
def icgc_pipeline(gene_symbols, cancer_site=None,
                     output_dir="results"):
    """
    ICGC 統合パイプライン。

    Parameters:
        gene_symbols: list[str] — 遺伝子リスト
        cancer_site: str — がん部位フィルタ
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) プロジェクト検索
    projects = icgc_search_projects(query=cancer_site)
    projects.to_csv(output_dir / "projects.csv", index=False)

    # 2) 遺伝子別変異検索
    all_mutations = []
    for gene in gene_symbols:
        try:
            muts = icgc_search_mutations(
                gene_symbol=gene, limit=200)
            muts["query_gene"] = gene
            all_mutations.append(muts)
        except Exception as e:
            print(f"  Warning: {gene} — {e}")
            continue

    if all_mutations:
        combined = pd.concat(all_mutations,
                              ignore_index=True)
        combined.to_csv(output_dir / "mutations.csv",
                        index=False)

    # 3) がん種統計
    if not projects.empty:
        top_project = projects.iloc[0]["project_id"]
        stats = icgc_cancer_stats(project_id=top_project)
        pd.DataFrame([stats]).to_csv(
            output_dir / "cancer_stats.csv", index=False)

    print(f"ICGC pipeline: {output_dir}")
    return {"projects": projects}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| (direct) | ICGC DCC API | 直接 REST API — TU 非連携 |

## パイプライン統合

```
cancer-genomics → icgc-cancer-data → precision-oncology
  (がんゲノム全般)  (ICGC DCC API)   (精密腫瘍学)
       │                 │                 ↓
  tcga-data ────────────┘        clinical-decision-support
  (TCGA データ)      │            (臨床意思決定)
                     ↓
          variant-interpretation
          (変異臨床解釈)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/projects.csv` | プロジェクト一覧 | → cancer-genomics |
| `results/mutations.csv` | 体細胞変異 | → variant-interpretation |
| `results/cancer_stats.csv` | がん種統計 | → precision-oncology |
