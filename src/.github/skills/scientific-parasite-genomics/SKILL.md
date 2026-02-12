---
name: scientific-parasite-genomics
description: |
  寄生虫ゲノミクススキル。PlasmoDB/VectorBase/ToxoDB REST API
  による寄生虫ゲノム検索・遺伝子情報・薬剤標的同定・比較
  ゲノミクス。直接 REST API 連携 (TU 外)。
tu_tools: []
---

# Scientific Parasite Genomics

VEuPathDB ファミリー (PlasmoDB, VectorBase, ToxoDB, TriTrypDB)
の REST API を活用した寄生虫ゲノミクス解析パイプラインを提供
する。

## When to Use

- マラリア原虫ゲノム (PlasmoDB) を検索するとき
- 蚊・ダニ等の媒介生物ゲノム (VectorBase) を検索するとき
- トキソプラズマゲノム (ToxoDB) を検索するとき
- トリパノソーマ/リーシュマニアゲノム (TriTrypDB) を検索するとき
- 寄生虫の薬剤標的候補を同定するとき
- 寄生虫間の比較ゲノミクスを実施するとき

---

## Quick Start

## 1. VEuPathDB 遺伝子検索

```python
import requests
import pandas as pd
import numpy as np

VEUPATHDB_SITES = {
    "plasmo": "https://plasmodb.org/plasmo/service",
    "vector": "https://vectorbase.org/vectorbase/service",
    "toxo": "https://toxodb.org/toxo/service",
    "tritryp": "https://tritrypdb.org/tritrypdb/service",
}


def veupathdb_search_genes(organism, query, db="plasmo",
                             limit=100):
    """
    VEuPathDB — 遺伝子検索。

    Parameters:
        organism: str — 生物種名 (例: "Plasmodium falciparum 3D7")
        query: str — 検索キーワード (例: "kinase", "transporter")
        db: str — データベース ("plasmo", "vector", "toxo", "tritryp")
        limit: int — 最大結果数
    """
    base = VEUPATHDB_SITES.get(db, VEUPATHDB_SITES["plasmo"])
    url = f"{base}/record-types/gene/searches/GenesByTextSearch"

    payload = {
        "searchConfig": {
            "parameters": {
                "text_expression": query,
                "text_fields": "Gene ID,Gene Name or Symbol,"
                                "Gene product",
                "organism": [organism],
            }
        },
        "reportConfig": {
            "attributes": ["primary_key", "gene_name",
                           "gene_product", "gene_type",
                           "chromosome", "start_min",
                           "end_max", "strand"],
            "pagination": {"offset": 0, "numRecords": limit},
        },
    }
    headers = {"Content-Type": "application/json"}
    resp = requests.post(url, json=payload, headers=headers,
                          timeout=60)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for rec in data.get("records", []):
        attrs = rec.get("attributes", {})
        results.append({
            "gene_id": attrs.get("primary_key", ""),
            "gene_name": attrs.get("gene_name", ""),
            "product": attrs.get("gene_product", ""),
            "gene_type": attrs.get("gene_type", ""),
            "chromosome": attrs.get("chromosome", ""),
            "start": attrs.get("start_min", None),
            "end": attrs.get("end_max", None),
            "strand": attrs.get("strand", ""),
        })

    df = pd.DataFrame(results)
    print(f"VEuPathDB ({db}) genes: {len(df)} results "
          f"(organism={organism}, query={query})")
    return df
```

## 2. 遺伝子機能アノテーション

```python
def veupathdb_gene_annotation(gene_id, db="plasmo"):
    """
    VEuPathDB — 遺伝子機能アノテーション取得。

    Parameters:
        gene_id: str — 遺伝子 ID (例: "PF3D7_1133400")
        db: str — データベース
    """
    base = VEUPATHDB_SITES.get(db, VEUPATHDB_SITES["plasmo"])
    url = f"{base}/record-types/gene/records/{gene_id}"

    params = {
        "attributes": "all",
        "tables": "GoTerms,InterPro,MetabolicPathways,"
                  "PubMed,EcNumber",
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    attrs = data.get("attributes", {})
    tables = data.get("tables", {})

    annotation = {
        "gene_id": gene_id,
        "gene_name": attrs.get("gene_name", ""),
        "product": attrs.get("gene_product", ""),
        "molecular_weight": attrs.get("molecular_weight", ""),
        "isoelectric_point": attrs.get("isoelectric_point", ""),
        "signal_peptide": attrs.get("signal_peptide", ""),
        "transmembrane_domains": attrs.get("transmembrane_domains", ""),
    }

    # GO Term 取得
    go_terms = []
    for go_rec in tables.get("GoTerms", []):
        go_terms.append({
            "go_id": go_rec.get("go_id", ""),
            "go_term": go_rec.get("go_term_name", ""),
            "ontology": go_rec.get("ontology", ""),
            "evidence": go_rec.get("evidence_code", ""),
        })
    annotation["go_terms"] = go_terms

    # InterPro ドメイン
    domains = []
    for d in tables.get("InterPro", []):
        domains.append({
            "interpro_id": d.get("interpro_primary_id", ""),
            "name": d.get("interpro_name", ""),
            "description": d.get("interpro_description", ""),
        })
    annotation["domains"] = domains

    print(f"VEuPathDB annotation: {gene_id}, "
          f"{len(go_terms)} GO terms, {len(domains)} domains")
    return annotation
```

## 3. 薬剤標的候補スクリーニング

```python
def parasite_drug_target_screen(organism, db="plasmo",
                                  essentiality_threshold=0.5):
    """
    寄生虫ゲノム — 薬剤標的候補スクリーニング。

    Parameters:
        organism: str — 生物種
        db: str — データベース
        essentiality_threshold: float — 必須性スコア閾値
    """
    # キナーゼ検索
    kinases = veupathdb_search_genes(organism, "kinase", db=db)
    # プロテアーゼ検索
    proteases = veupathdb_search_genes(organism, "protease", db=db)
    # トランスポーター検索
    transporters = veupathdb_search_genes(
        organism, "transporter", db=db)

    all_targets = pd.concat([kinases, proteases, transporters],
                              ignore_index=True)
    all_targets = all_targets.drop_duplicates(subset=["gene_id"])

    # 薬剤標的性スコア (ヒューリスティック)
    all_targets["target_class"] = "unknown"
    all_targets.loc[
        all_targets["gene_id"].isin(kinases["gene_id"]),
        "target_class"] = "kinase"
    all_targets.loc[
        all_targets["gene_id"].isin(proteases["gene_id"]),
        "target_class"] = "protease"
    all_targets.loc[
        all_targets["gene_id"].isin(transporters["gene_id"]),
        "target_class"] = "transporter"

    print(f"Drug target screen: {len(all_targets)} candidates "
          f"(kinases={len(kinases)}, proteases={len(proteases)}, "
          f"transporters={len(transporters)})")
    return all_targets
```

## 4. 寄生虫ゲノミクス統合パイプライン

```python
def parasite_genomics_pipeline(organism, query,
                                  db="plasmo",
                                  output_dir="results"):
    """
    寄生虫ゲノミクス統合パイプライン。

    Parameters:
        organism: str — 生物種 (例: "Plasmodium falciparum 3D7")
        query: str — 検索クエリ
        db: str — データベース
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 遺伝子検索
    genes = veupathdb_search_genes(organism, query, db=db)
    genes.to_csv(output_dir / "genes.csv", index=False)

    # 2) トップ遺伝子のアノテーション
    annotations = []
    for gene_id in genes["gene_id"].head(10):
        try:
            ann = veupathdb_gene_annotation(gene_id, db=db)
            annotations.append(ann)
        except Exception:
            continue
    ann_df = pd.DataFrame([{
        k: v for k, v in a.items()
        if not isinstance(v, list)
    } for a in annotations])
    ann_df.to_csv(output_dir / "annotations.csv", index=False)

    # 3) 薬剤標的スクリーニング
    targets = parasite_drug_target_screen(organism, db=db)
    targets.to_csv(output_dir / "drug_targets.csv", index=False)

    print(f"Parasite genomics pipeline: {output_dir}")
    return {
        "genes": genes,
        "annotations": annotations,
        "drug_targets": targets,
    }
```

---

## ToolUniverse 連携

直接 REST API 使用 (VEuPathDB は ToolUniverse 外)。

## パイプライン統合

```
infectious-disease → parasite-genomics → phylogenetics
  (病原体情報)      (寄生虫ゲノム)       (系統解析)
       │                   │                  ↓
  drug-discovery ─────────┘           comparative-genomics
  (薬剤探索)           │              (比較ゲノミクス)
                        ↓
                  pathway-enrichment
                  (パスウェイ解析)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/genes.csv` | 遺伝子一覧 | → phylogenetics |
| `results/annotations.csv` | 機能アノテーション | → pathway-enrichment |
| `results/drug_targets.csv` | 薬剤標的候補 | → drug-discovery |
