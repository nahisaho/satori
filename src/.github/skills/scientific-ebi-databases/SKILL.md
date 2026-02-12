---
name: scientific-ebi-databases
description: |
  EBI データベース群統合アクセススキル。EBI Search 横断検索、ENA Browser
  ヌクレオチドアーカイブ、BioStudies 研究データ、dbfetch エントリ取得、
  MetaboLights メタボロミクスリポジトリの統合パイプライン。
---

# Scientific EBI Databases

EBI Search / ENA Browser / BioStudies / dbfetch / MetaboLights を統合した
EBI データベース群アクセスパイプラインを提供する。

## When to Use

- EBI Search で複数データベースを横断検索するとき
- ENA (European Nucleotide Archive) で配列データを検索するとき
- BioStudies で研究プロジェクトデータを探すとき
- dbfetch でエントリを一括取得するとき
- MetaboLights でメタボロミクス実験データにアクセスするとき

---

## Quick Start

## 1. EBI Search 横断検索

```python
import requests
import pandas as pd

EBI_SEARCH_API = "https://www.ebi.ac.uk/ebisearch/ws/rest"


def search_ebi(query, domain="allebi", size=25, fields=None):
    """
    EBI Search 横断検索 — 複数 EBI データベースを一括検索。

    Parameters:
        query: str — 検索クエリ
        domain: str — 検索ドメイン ("allebi", "uniprot", "pdb", "ena", etc.)
        size: int — 最大取得数
        fields: list — 返却フィールド

    ToolUniverse:
        EBI_Search_query(query=query, domain=domain)
        EBI_Search_get_entry(domain=domain, entry_id=entry_id)
    """
    params = {
        "query": query,
        "size": size,
        "format": "json",
    }
    if fields:
        params["fields"] = ",".join(fields)

    resp = requests.get(f"{EBI_SEARCH_API}/{domain}", params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for entry in data.get("entries", []):
        row = {"id": entry.get("id", ""), "source": entry.get("source", "")}
        for field in entry.get("fields", {}):
            row[field] = entry["fields"][field][0] if entry["fields"][field] else ""
        results.append(row)

    df = pd.DataFrame(results)
    total = data.get("hitCount", 0)
    print(f"EBI Search [{domain}] '{query}': {total} total hits, {len(df)} returned")
    return df
```

## 2. ENA (European Nucleotide Archive) 配列検索

```python
ENA_API = "https://www.ebi.ac.uk/ena/browser/api"


def search_ena(query, result_type="sequence", limit=100):
    """
    ENA ヌクレオチドアーカイブ検索。

    Parameters:
        query: str — 検索クエリ or Taxon ID
        result_type: str — "sequence", "read_run", "analysis", "study"
        limit: int — 最大取得数

    ToolUniverse:
        ENA_search(query=query, result=result_type)
        ENA_get_entry(accession=accession)
    """
    params = {
        "query": query,
        "result": result_type,
        "limit": limit,
        "format": "json",
    }
    resp = requests.get(f"{ENA_API}/search", params=params)
    resp.raise_for_status()
    data = resp.json()

    df = pd.DataFrame(data) if isinstance(data, list) else pd.DataFrame()
    print(f"ENA search '{query}' [{result_type}]: {len(df)} entries")
    return df


def get_ena_entry(accession, display="json"):
    """
    ENA アクセッション番号によるエントリ取得。

    Parameters:
        accession: str — ENA accession (e.g., "ERS000001", "ERR000001")
    """
    resp = requests.get(
        f"{ENA_API}/entry/{accession}",
        params={"display": display}
    )
    resp.raise_for_status()
    print(f"ENA entry {accession}: retrieved")
    return resp.json() if display == "json" else resp.text
```

## 3. BioStudies 研究データ検索

```python
BIOSTUDIES_API = "https://www.ebi.ac.uk/biostudies/api/v1"


def search_biostudies(query, page_size=25):
    """
    BioStudies 研究プロジェクトデータ検索。

    Parameters:
        query: str — 検索クエリ
        page_size: int — ページサイズ

    ToolUniverse:
        BioStudies_search(query=query)
        BioStudies_get_study(accession=accession)
    """
    params = {"query": query, "pageSize": page_size}
    resp = requests.get(f"{BIOSTUDIES_API}/search", params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for hit in data.get("hits", []):
        results.append({
            "accession": hit.get("accno", ""),
            "title": hit.get("title", ""),
            "author": hit.get("author", ""),
            "release_date": hit.get("rtime", ""),
            "type": hit.get("type", ""),
            "files": hit.get("files", 0),
            "links": hit.get("links", 0),
        })

    df = pd.DataFrame(results)
    total = data.get("totalHits", 0)
    print(f"BioStudies search '{query}': {total} total, {len(df)} returned")
    return df
```

## 4. dbfetch エントリ一括取得

```python
DBFETCH_API = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch"


def dbfetch(db, ids, format_type="json", style="raw"):
    """
    dbfetch — EBI データベースエントリ一括取得。

    Parameters:
        db: str — データベース名 (e.g., "uniprotkb", "embl", "pdb")
        ids: list — ID リスト
        format_type: str — 出力形式 ("json", "fasta", "xml")
        style: str — スタイル ("raw", "html")

    ToolUniverse:
        dbfetch_get_entries(db=db, ids=ids, format=format_type)
    """
    ids_str = ",".join(ids) if isinstance(ids, list) else ids
    params = {
        "db": db,
        "id": ids_str,
        "format": format_type,
        "style": style,
    }
    resp = requests.get(DBFETCH_API, params=params)
    resp.raise_for_status()

    print(f"dbfetch [{db}]: {len(ids) if isinstance(ids, list) else 1} entries, "
          f"format={format_type}")
    if format_type == "json":
        return resp.json()
    return resp.text
```

## 5. MetaboLights メタボロミクスリポジトリ

```python
METABOLIGHTS_API = "https://www.ebi.ac.uk/metabolights/ws"


def search_metabolights(query):
    """
    MetaboLights メタボロミクス実験データ検索。

    Parameters:
        query: str — 検索クエリ (化合物名、疾患名、生物種)

    ToolUniverse:
        MetaboLights_search_studies(query=query)
        MetaboLights_get_study(study_id=study_id)
    """
    resp = requests.get(
        f"{METABOLIGHTS_API}/studies/search",
        params={"query": query}
    )
    resp.raise_for_status()
    data = resp.json()

    results = []
    for study in data.get("content", []):
        results.append({
            "study_id": study.get("studyIdentifier", ""),
            "title": study.get("title", ""),
            "organism": study.get("organism", ""),
            "description": (study.get("description") or "")[:200],
            "submission_date": study.get("submissionDate", ""),
            "status": study.get("studyStatus", ""),
        })

    df = pd.DataFrame(results)
    print(f"MetaboLights search '{query}': {len(df)} studies")
    return df


def get_metabolights_study(study_id):
    """MetaboLights 個別研究取得。"""
    resp = requests.get(f"{METABOLIGHTS_API}/studies/{study_id}")
    resp.raise_for_status()
    data = resp.json()
    print(f"MetaboLights {study_id}: {data.get('title', '')[:80]}")
    return data
```

---

## 利用可能ツール

| ToolUniverse カテゴリ | 主なツール |
|---|---|
| `ebi_search` | `EBI_Search_query`, `EBI_Search_get_entry` |
| `ena_browser` | `ENA_search`, `ENA_get_entry` |
| `biostudies` | `BioStudies_search`, `BioStudies_get_study` |
| `dbfetch` | `dbfetch_get_entries` |
| `metabolights` | `MetaboLights_search_studies`, `MetaboLights_get_study` |

## パイプライン出力

| 出力ファイル | 説明 | 連携先スキル |
|---|---|---|
| `results/ebi_search.csv` | EBI 横断検索結果 | → bioinformatics, literature-search |
| `results/ena_sequences.fasta` | ENA 配列データ | → genome-sequence-tools, sequence-analysis |
| `results/biostudies_metadata.json` | 研究プロジェクト情報 | → multi-omics, systematic-review |
| `results/metabolights_study.json` | メタボロミクスデータ | → metabolomics, metabolomics-databases |

## パイプライン統合

```
genome-sequence-tools ──→ ebi-databases ──→ metabolomics-databases
  (NCBI/BLAST)           (ENA/EBI Search)   (HMDB/MetaCyc)
                              │
                              ├──→ bioinformatics (配列データ)
                              ├──→ sequence-analysis (FASTA)
                              └──→ structural-proteomics (PDBe cross-ref)
```
