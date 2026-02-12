---
name: scientific-uniprot-proteome
description: |
  UniProt プロテオームスキル。UniProt REST API による
  タンパク質検索・ID マッピング・配列取得・機能アノテーション・
  UniRef/UniParc 横断検索。ToolUniverse 連携: uniprot。
tu_tools:
  - key: uniprot
    name: UniProt
    description: タンパク質配列・機能アノテーション・ID マッピング
---

# Scientific UniProt Proteome

UniProt REST API を活用したタンパク質検索・ID マッピング・
機能アノテーション取得・プロテオーム解析パイプラインを提供する。

## When to Use

- タンパク質のアミノ酸配列・機能情報を検索するとき
- UniProt ID と他データベース ID を相互変換するとき
- タンパク質のドメイン・GO アノテーション・疾患関連を調べるとき
- プロテオーム規模でのタンパク質機能解析を行うとき
- UniRef クラスター・UniParc アーカイブを横断検索するとき
- 生物種ごとのリファレンスプロテオームを取得するとき

---

## Quick Start

## 1. タンパク質検索・エントリ取得

```python
import requests
import pandas as pd

UNIPROT_BASE = "https://rest.uniprot.org"


def uniprot_search(query, organism=None, reviewed=True,
                      limit=50, fields=None):
    """
    UniProt — タンパク質検索。

    Parameters:
        query: str — 検索クエリ (例: "BRCA1", "kinase")
        organism: str — 生物種 (例: "9606" for Human)
        reviewed: bool — Swiss-Prot のみ (True) / TrEMBL 含む (False)
        limit: int — 最大結果数
        fields: list[str] — 返却フィールド
    """
    url = f"{UNIPROT_BASE}/uniprotkb/search"
    default_fields = [
        "accession", "id", "protein_name", "gene_names",
        "organism_name", "length", "reviewed",
        "go_p", "go_f", "go_c",
    ]
    params = {
        "query": query,
        "size": min(limit, 500),
        "fields": ",".join(fields or default_fields),
        "format": "json",
    }
    if organism:
        params["query"] += f" AND organism_id:{organism}"
    if reviewed:
        params["query"] += " AND reviewed:true"

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for entry in data.get("results", []):
        protein_name = ""
        if entry.get("proteinDescription"):
            rec = entry["proteinDescription"].get(
                "recommendedName")
            if rec:
                protein_name = rec.get("fullName", {}).get(
                    "value", "")

        genes = [g.get("geneName", {}).get("value", "")
                 for g in entry.get("genes", [])]

        results.append({
            "accession": entry.get("primaryAccession", ""),
            "entry_name": entry.get("uniProtkbId", ""),
            "protein_name": protein_name,
            "gene_names": "; ".join(genes),
            "organism": entry.get("organism", {}).get(
                "scientificName", ""),
            "length": entry.get("sequence", {}).get(
                "length", 0),
            "reviewed": entry.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)",
        })

    df = pd.DataFrame(results)
    print(f"UniProt search: {len(df)} entries "
          f"(query='{query}')")
    return df


def uniprot_get_entry(accession, format="json"):
    """
    UniProt — エントリ詳細取得。

    Parameters:
        accession: str — UniProt ID (例: "P38398")
        format: str — "json" or "fasta"
    """
    url = f"{UNIPROT_BASE}/uniprotkb/{accession}"
    resp = requests.get(url, params={"format": format},
                        timeout=30)
    resp.raise_for_status()
    if format == "fasta":
        return resp.text
    return resp.json()
```

## 2. ID マッピング

```python
import time


def uniprot_id_mapping(from_db, to_db, ids):
    """
    UniProt — ID マッピング (非同期ジョブ)。

    Parameters:
        from_db: str — 変換元 DB (例: "UniProtKB_AC-ID")
        to_db: str — 変換先 DB (例: "Gene_Name", "PDB",
            "Ensembl", "RefSeq_Protein")
        ids: list[str] — ID リスト
    """
    # ジョブ投入
    url = f"{UNIPROT_BASE}/idmapping/run"
    resp = requests.post(url, data={
        "from": from_db, "to": to_db,
        "ids": ",".join(ids),
    }, timeout=30)
    resp.raise_for_status()
    job_id = resp.json()["jobId"]

    # ポーリング
    status_url = f"{UNIPROT_BASE}/idmapping/status/{job_id}"
    for _ in range(30):
        s = requests.get(status_url, timeout=30).json()
        if "results" in s or "redirectURL" in s:
            break
        time.sleep(2)

    # 結果取得
    result_url = (f"{UNIPROT_BASE}/idmapping/results/"
                  f"{job_id}")
    r = requests.get(result_url, timeout=30)
    r.raise_for_status()
    data = r.json()

    results = []
    for item in data.get("results", []):
        results.append({
            "from_id": item.get("from", ""),
            "to_id": item.get("to", ""),
        })

    df = pd.DataFrame(results)
    print(f"UniProt ID mapping: {len(df)} mappings "
          f"({from_db} → {to_db})")
    return df
```

## 3. 機能アノテーション取得

```python
def uniprot_get_features(accession):
    """
    UniProt — タンパク質機能フィーチャー取得。

    Parameters:
        accession: str — UniProt ID
    """
    entry = uniprot_get_entry(accession)

    features = []
    for f in entry.get("features", []):
        loc = f.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        features.append({
            "type": f.get("type", ""),
            "description": f.get("description", ""),
            "start": start,
            "end": end,
            "evidence": len(f.get("evidences", [])),
        })

    df = pd.DataFrame(features)
    print(f"UniProt features: {accession} → "
          f"{len(df)} features")
    return df
```

## 4. UniProt 統合パイプライン

```python
def uniprot_pipeline(gene_names, organism="9606",
                        output_dir="results"):
    """
    UniProt 統合パイプライン。

    Parameters:
        gene_names: list[str] — 遺伝子シンボルリスト
        organism: str — 生物種 Taxonomy ID
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_entries = []
    for gene in gene_names:
        try:
            df = uniprot_search(gene, organism=organism)
            all_entries.append(df)
        except Exception as e:
            print(f"  Warning: {gene} — {e}")

    if all_entries:
        combined = pd.concat(all_entries, ignore_index=True)
        combined.to_csv(output_dir / "uniprot_entries.csv",
                        index=False)

        # トップエントリの機能フィーチャー
        if not combined.empty:
            top_acc = combined.iloc[0]["accession"]
            features = uniprot_get_features(top_acc)
            features.to_csv(
                output_dir / "protein_features.csv",
                index=False)

    print(f"UniProt pipeline: {output_dir}")
    return {"entries": combined if all_entries else pd.DataFrame()}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `uniprot` | UniProt | タンパク質検索・ID マッピング・配列・機能アノテーション |

## パイプライン統合

```
protein-structure-analysis → uniprot-proteome → protein-design
  (PDB/AlphaFold 構造)     (UniProt REST API)  (de novo 設計)
          │                       │                    ↓
  alphafold-structures ──────────┘           drug-target-profiling
  (AlphaFold DB)            │               (標的プロファイリング)
                            ↓
                  protein-domain-family
                  (InterPro/Pfam)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/uniprot_entries.csv` | タンパク質検索結果 | → protein-structure-analysis |
| `results/protein_features.csv` | 機能フィーチャー | → protein-domain-family |
