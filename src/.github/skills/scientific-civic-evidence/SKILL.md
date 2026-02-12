---
name: scientific-civic-evidence
description: |
  CIViC 臨床エビデンススキル。CIViC (Clinical Interpretation
  of Variants in Cancer) REST API を用いたバリアント臨床解釈・
  エビデンスアイテム・分子プロファイル・アサーション検索。
  ToolUniverse 連携: civic。
tu_tools:
  - key: civic
    name: CIViC
    description: がんバリアント臨床解釈データベース
---

# Scientific CIViC Evidence

CIViC (Clinical Interpretation of Variants in Cancer) REST API
を活用したバリアント臨床解釈・エビデンスアイテム取得・
分子プロファイル・アサーションパイプラインを提供する。

## When to Use

- がんバリアントの臨床的解釈を検索するとき
- エビデンスアイテム (薬剤応答・予後・診断) を取得するとき
- 遺伝子ごとのバリアントサマリーを確認するとき
- 分子プロファイル (Molecular Profile) を検索するとき
- アサーション (ガイドライン推奨) を取得するとき

---

## Quick Start

## 1. バリアント検索・臨床解釈

```python
import requests
import pandas as pd

CIVIC_API = "https://civicdb.org/api"


def civic_variant_search(gene_name, variant_name=None,
                            limit=50):
    """
    CIViC — バリアント検索。

    Parameters:
        gene_name: str — 遺伝子名 (例: "BRAF")
        variant_name: str — バリアント名
            (例: "V600E")
        limit: int — 最大結果数
    """
    url = f"{CIVIC_API}/variants"
    params = {"count": limit}

    # 遺伝子名で検索
    gene_url = f"{CIVIC_API}/genes/{gene_name}"
    try:
        resp = requests.get(gene_url, timeout=30)
        if resp.status_code == 200:
            gene_data = resp.json()
        else:
            # 検索 API フォールバック
            search_url = f"{CIVIC_API}/genes"
            params_g = {"name": gene_name, "count": 5}
            resp = requests.get(search_url,
                                params=params_g,
                                timeout=30)
            resp.raise_for_status()
            records = resp.json().get("records", [])
            gene_data = records[0] if records else {}
    except Exception as e:
        print(f"  CIViC gene lookup: {e}")
        gene_data = {}

    if not gene_data:
        return pd.DataFrame()

    variants = gene_data.get("variants", [])
    rows = []
    for v in variants[:limit]:
        name = v.get("name", "")
        if variant_name and variant_name.lower() \
                not in name.lower():
            continue
        rows.append({
            "variant_id": v.get("id", ""),
            "gene": gene_name,
            "variant_name": name,
            "description": (v.get("description", "")
                            [:200]),
            "evidence_count": len(
                v.get("evidence_items", [])),
        })

    df = pd.DataFrame(rows)
    print(f"CIViC variants: {gene_name} → {len(df)}")
    return df


def civic_gene_summary(gene_name):
    """
    CIViC — 遺伝子サマリー取得。

    Parameters:
        gene_name: str — 遺伝子名 (例: "EGFR")
    """
    url = f"{CIVIC_API}/genes/{gene_name}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    result = {
        "gene_id": data.get("id", ""),
        "name": data.get("name", ""),
        "description": data.get("description", ""),
        "n_variants": len(data.get("variants", [])),
        "aliases": "; ".join(
            data.get("aliases", [])),
    }
    return result
```

## 2. エビデンスアイテム取得

```python
def civic_evidence_items(variant_id, limit=50):
    """
    CIViC — エビデンスアイテム取得。

    Parameters:
        variant_id: int — バリアント ID
        limit: int — 最大結果数
    """
    url = f"{CIVIC_API}/variants/{variant_id}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for ev in data.get("evidence_items", [])[:limit]:
        drugs = [d.get("name", "")
                 for d in ev.get("drugs", [])]
        rows.append({
            "evidence_id": ev.get("id", ""),
            "variant_id": variant_id,
            "evidence_type": ev.get(
                "evidence_type", ""),
            "evidence_level": ev.get(
                "evidence_level", ""),
            "evidence_direction": ev.get(
                "evidence_direction", ""),
            "clinical_significance": ev.get(
                "clinical_significance", ""),
            "disease": ev.get("disease", {}).get(
                "name", ""),
            "drugs": "; ".join(drugs),
            "rating": ev.get("rating", ""),
            "status": ev.get("status", ""),
            "source_citation": ev.get(
                "source", {}).get("citation", ""),
        })

    df = pd.DataFrame(rows)
    print(f"CIViC evidence: variant {variant_id} "
          f"→ {len(df)} items")
    return df
```

## 3. アサーション取得

```python
def civic_assertions(gene_name=None, limit=50):
    """
    CIViC — アサーション (ガイドライン推奨) 取得。

    Parameters:
        gene_name: str — 遺伝子名フィルタ
        limit: int — 最大結果数
    """
    url = f"{CIVIC_API}/assertions"
    params = {"count": limit}

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for a in data.get("records", []):
        genes = [g.get("name", "")
                 for g in a.get("genes", [])]
        if gene_name and gene_name not in genes:
            continue
        drugs = [d.get("name", "")
                 for d in a.get("drugs", [])]
        rows.append({
            "assertion_id": a.get("id", ""),
            "genes": "; ".join(genes),
            "variant": a.get("variant", {}).get(
                "name", ""),
            "disease": a.get("disease", {}).get(
                "name", ""),
            "drugs": "; ".join(drugs),
            "assertion_type": a.get(
                "assertion_type", ""),
            "assertion_direction": a.get(
                "assertion_direction", ""),
            "clinical_significance": a.get(
                "clinical_significance", ""),
            "amp_level": a.get("amp_level", ""),
            "status": a.get("status", ""),
        })

    df = pd.DataFrame(rows)
    print(f"CIViC assertions: {len(df)}")
    return df
```

## 4. CIViC 統合パイプライン

```python
def civic_pipeline(gene_name, variant_name=None,
                      output_dir="results"):
    """
    CIViC 統合パイプライン。

    Parameters:
        gene_name: str — 遺伝子名 (例: "BRAF")
        variant_name: str — バリアント名 (例: "V600E")
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 遺伝子サマリー
    summary = civic_gene_summary(gene_name)
    pd.DataFrame([summary]).to_csv(
        output_dir / "civic_gene.csv", index=False)

    # 2) バリアント検索
    variants = civic_variant_search(gene_name,
                                    variant_name)
    variants.to_csv(output_dir / "civic_variants.csv",
                    index=False)

    # 3) トップバリアントのエビデンス
    if not variants.empty:
        top_vid = variants.iloc[0]["variant_id"]
        evidence = civic_evidence_items(top_vid)
        evidence.to_csv(
            output_dir / "civic_evidence.csv",
            index=False)

    # 4) アサーション
    assertions = civic_assertions(gene_name)
    assertions.to_csv(
        output_dir / "civic_assertions.csv",
        index=False)

    print(f"CIViC pipeline: {gene_name} → {output_dir}")
    return {"variants": variants}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `civic` | CIViC | がんバリアント臨床解釈 (~12 tools) |

## パイプライン統合

```
variant-interpretation → civic-evidence → precision-oncology
  (ClinVar バリアント)    (CIViC REST)    (精密腫瘍学)
          │                    │                ↓
  gnomad-variants ────────────┘     drug-target-profiling
  (集団頻度)           │           (標的プロファイリング)
                       ↓
             opentargets-genetics
             (OT 標的-疾患)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/civic_gene.csv` | 遺伝子サマリー | → cancer-genomics |
| `results/civic_variants.csv` | バリアント一覧 | → variant-interpretation |
| `results/civic_evidence.csv` | エビデンス | → precision-oncology |
| `results/civic_assertions.csv` | アサーション | → pharmacogenomics |
