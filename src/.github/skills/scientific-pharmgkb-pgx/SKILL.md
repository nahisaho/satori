---
name: scientific-pharmgkb-pgx
description: |
  PharmGKB 薬理ゲノミクススキル。PharmGKB REST API による
  臨床アノテーション・薬物遺伝子関連・投与量ガイドライン・
  スターアレル解析。ToolUniverse 連携: pharmgkb。
tu_tools:
  - key: pharmgkb
    name: PharmGKB
    description: 臨床アノテーション・薬物遺伝子関連・PGx ガイドライン
---

# Scientific PharmGKB PGx

PharmGKB (Pharmacogenomics Knowledgebase) REST API を活用した
薬理ゲノミクス臨床アノテーション・薬物遺伝子相互作用・投与量
ガイドライン検索パイプラインを提供する。

## When to Use

- 薬物と遺伝子変異の関連を調べるとき
- 臨床アノテーション (エビデンスレベル付き) を検索するとき
- 投与量調整ガイドライン (CPIC/DPWG) を取得するとき
- スターアレルと表現型の対応を確認するとき
- 特定薬物の薬理ゲノミクス情報を包括的に取得するとき
- 精密医療の薬物選択を支援するとき

---

## Quick Start

## 1. 薬物・遺伝子検索

```python
import requests
import pandas as pd

PGKB_BASE = "https://api.pharmgkb.org/v1/data"


def pharmgkb_search_drugs(query, limit=50):
    """
    PharmGKB — 薬物検索。

    Parameters:
        query: str — 薬物名 (例: "warfarin", "clopidogrel")
        limit: int — 最大結果数
    """
    url = f"{PGKB_BASE}/chemical"
    params = {"name": query, "view": "max"}

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for item in data.get("data", []):
        results.append({
            "pharmgkb_id": item.get("id", ""),
            "name": item.get("name", ""),
            "generic_names": "; ".join(
                item.get("genericNames", [])),
            "trade_names": "; ".join(
                item.get("tradeNames", [])[:5]),
            "type": item.get("type", ""),
            "cross_references": len(
                item.get("crossReferences", [])),
        })

    df = pd.DataFrame(results)
    print(f"PharmGKB drugs: {len(df)} results "
          f"(query='{query}')")
    return df


def pharmgkb_search_genes(query, limit=50):
    """
    PharmGKB — 遺伝子検索。

    Parameters:
        query: str — 遺伝子シンボル (例: "CYP2D6")
        limit: int — 最大結果数
    """
    url = f"{PGKB_BASE}/gene"
    params = {"symbol": query, "view": "max"}

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for item in data.get("data", []):
        results.append({
            "pharmgkb_id": item.get("id", ""),
            "symbol": item.get("symbol", ""),
            "name": item.get("name", ""),
            "chromosome": item.get("chromosomeFormatted", ""),
            "cpic_gene": item.get("cpicGene", False),
            "has_prescribing_info": item.get(
                "hasPrescribingInfo", False),
        })

    df = pd.DataFrame(results)
    print(f"PharmGKB genes: {len(df)} results "
          f"(query='{query}')")
    return df
```

## 2. 臨床アノテーション取得

```python
def pharmgkb_clinical_annotations(gene_or_drug,
                                     search_type="gene"):
    """
    PharmGKB — 臨床アノテーション検索。

    Parameters:
        gene_or_drug: str — 遺伝子シンボル or 薬物名
        search_type: str — "gene" or "drug"
    """
    url = f"{PGKB_BASE}/clinicalAnnotation"
    params = {"view": "max"}

    if search_type == "gene":
        # 遺伝子で検索
        gene_url = f"{PGKB_BASE}/gene"
        g_resp = requests.get(gene_url,
                              params={"symbol": gene_or_drug},
                              timeout=30)
        g_resp.raise_for_status()
        genes = g_resp.json().get("data", [])
        if genes:
            params["relatedGenes.id"] = genes[0].get("id", "")
    else:
        # 薬物で検索
        drug_url = f"{PGKB_BASE}/chemical"
        d_resp = requests.get(drug_url,
                              params={"name": gene_or_drug},
                              timeout=30)
        d_resp.raise_for_status()
        drugs = d_resp.json().get("data", [])
        if drugs:
            params["relatedChemicals.id"] = drugs[0].get("id", "")

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for item in data.get("data", []):
        genes = [g.get("symbol", "")
                 for g in item.get("relatedGenes", [])]
        drugs = [c.get("name", "")
                 for c in item.get("relatedChemicals", [])]
        results.append({
            "annotation_id": item.get("id", ""),
            "level": item.get("level", ""),
            "score": item.get("score", ""),
            "genes": "; ".join(genes),
            "drugs": "; ".join(drugs),
            "phenotype_category": item.get(
                "phenotypeCategory", ""),
            "sentences": (item.get("textHtml") or "")[:300],
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("level")

    print(f"PharmGKB annotations: {len(df)} "
          f"({search_type}={gene_or_drug})")
    return df
```

## 3. 投与量ガイドライン取得

```python
def pharmgkb_dosing_guidelines(drug_name=None, gene=None):
    """
    PharmGKB — 投与量ガイドライン (CPIC/DPWG) 検索。

    Parameters:
        drug_name: str — 薬物名
        gene: str — 遺伝子シンボル
    """
    url = f"{PGKB_BASE}/guideline"
    params = {"view": "max"}

    if drug_name:
        params["relatedChemicals.name"] = drug_name
    if gene:
        params["relatedGenes.symbol"] = gene

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for item in data.get("data", []):
        genes = [g.get("symbol", "")
                 for g in item.get("relatedGenes", [])]
        drugs = [c.get("name", "")
                 for c in item.get("relatedChemicals", [])]
        results.append({
            "guideline_id": item.get("id", ""),
            "name": item.get("name", ""),
            "source": item.get("source", ""),
            "genes": "; ".join(genes),
            "drugs": "; ".join(drugs),
            "recommendation": (item.get("textHtml") or "")[:500],
        })

    df = pd.DataFrame(results)
    print(f"PharmGKB guidelines: {len(df)} "
          f"(drug={drug_name}, gene={gene})")
    return df
```

## 4. PharmGKB 統合パイプライン

```python
def pharmgkb_pipeline(drug_name, genes=None,
                         output_dir="results"):
    """
    PharmGKB 統合パイプライン。

    Parameters:
        drug_name: str — 薬物名
        genes: list[str] — 関連遺伝子リスト
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 薬物検索
    drugs = pharmgkb_search_drugs(drug_name)
    drugs.to_csv(output_dir / "drugs.csv", index=False)

    # 2) 薬物の臨床アノテーション
    annotations = pharmgkb_clinical_annotations(
        drug_name, search_type="drug")
    annotations.to_csv(output_dir / "annotations.csv",
                       index=False)

    # 3) 投与量ガイドライン
    guidelines = pharmgkb_dosing_guidelines(
        drug_name=drug_name)
    guidelines.to_csv(output_dir / "guidelines.csv",
                      index=False)

    # 4) 関連遺伝子解析
    if genes:
        gene_results = []
        for g in genes:
            try:
                g_ann = pharmgkb_clinical_annotations(
                    g, search_type="gene")
                g_ann["query_gene"] = g
                gene_results.append(g_ann)
            except Exception:
                continue
        if gene_results:
            gene_df = pd.concat(gene_results,
                                ignore_index=True)
            gene_df.to_csv(
                output_dir / "gene_annotations.csv",
                index=False)

    print(f"PharmGKB pipeline: {output_dir}")
    return {
        "drugs": drugs,
        "annotations": annotations,
        "guidelines": guidelines,
    }
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `pharmgkb` | PharmGKB | 臨床アノテーション・薬物遺伝子・PGx ガイドライン |

## パイプライン統合

```
pharmacogenomics → pharmgkb-pgx → clinical-decision-support
  (PGx 解析全般)  (PharmGKB API)  (臨床意思決定)
       │                │                ↓
  drug-discovery ──────┘         precision-oncology
  (薬物開発)         │            (精密腫瘍学)
                     ↓
           variant-interpretation
           (変異臨床解釈)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/drugs.csv` | 薬物情報 | → drug-discovery |
| `results/annotations.csv` | 臨床アノテーション | → variant-interpretation |
| `results/guidelines.csv` | 投与量ガイドライン | → clinical-decision-support |
| `results/gene_annotations.csv` | 遺伝子別アノテーション | → pharmacogenomics |
