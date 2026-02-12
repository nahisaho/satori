---
name: scientific-rcsb-pdb-search
description: |
  RCSB PDB 構造検索スキル。RCSB PDB Search API および
  Data API によるタンパク質立体構造検索・メタデータ取得・
  リガンド情報・解像度フィルタリング。ToolUniverse 連携:
  rcsb_pdb, rcsb_search。
tu_tools:
  - key: rcsb_pdb
    name: RCSB PDB Data
    description: PDB エントリデータ取得・構造メタデータ
  - key: rcsb_search
    name: RCSB PDB Search
    description: PDB 構造検索・テキスト/配列/構造類似検索
---

# Scientific RCSB PDB Search

RCSB PDB Search API および Data API を活用したタンパク質立体構造
検索・メタデータ取得・リガンド情報パイプラインを提供する。

## When to Use

- PDB のタンパク質立体構造をテキスト検索するとき
- 解像度・実験手法でフィルタリングするとき
- リガンド結合構造を検索するとき
- 構造のメタデータ (著者・引用・解像度) を取得するとき
- 配列類似性で構造を検索するとき
- PDB エントリからリガンド・結合サイト情報を取得するとき

---

## Quick Start

## 1. テキスト検索・構造メタデータ

```python
import requests
import pandas as pd

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA = "https://data.rcsb.org/rest/v1/core"


def rcsb_text_search(query, method=None,
                        resolution_max=None, limit=50):
    """
    RCSB PDB — テキスト検索。

    Parameters:
        query: str — 検索クエリ (例: "BRCA1", "kinase")
        method: str — 実験手法フィルタ
            (例: "X-RAY DIFFRACTION", "ELECTRON MICROSCOPY")
        resolution_max: float — 最大解像度 (Å)
        limit: int — 最大結果数
    """
    search_query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {"value": query},
                }
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }

    if method:
        search_query["query"]["nodes"].append({
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "exptl.method",
                "operator": "exact_match",
                "value": method,
            },
        })

    if resolution_max:
        search_query["query"]["nodes"].append({
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info."
                             "resolution_combined",
                "operator": "less_or_equal",
                "value": resolution_max,
            },
        })

    resp = requests.post(RCSB_SEARCH, json=search_query,
                         timeout=30)
    resp.raise_for_status()
    data = resp.json()

    pdb_ids = [r["identifier"]
               for r in data.get("result_set", [])]
    total = data.get("total_count", 0)
    print(f"RCSB PDB search: {len(pdb_ids)}/{total} "
          f"(query='{query}')")
    return pdb_ids


def rcsb_get_entry(pdb_id):
    """
    RCSB PDB — エントリメタデータ取得。

    Parameters:
        pdb_id: str — PDB ID (例: "1BRS", "7S4O")
    """
    url = f"{RCSB_DATA}/entry/{pdb_id}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    info = data.get("rcsb_entry_info", {})
    citation = (data.get("rcsb_primary_citation") or {})

    result = {
        "pdb_id": pdb_id,
        "title": data.get("struct", {}).get("title", ""),
        "method": info.get("experimental_method", ""),
        "resolution": info.get("resolution_combined", [None])[0],
        "deposition_date": info.get("deposition_date", ""),
        "polymer_count": info.get(
            "deposited_polymer_entity_count", 0),
        "nonpolymer_count": info.get(
            "deposited_nonpolymer_entity_count", 0),
        "citation_title": citation.get("title", ""),
        "citation_doi": citation.get(
            "pdbx_database_id_doi", ""),
    }
    return result
```

## 2. 構造バッチ取得・比較

```python
def rcsb_batch_metadata(pdb_ids):
    """
    RCSB PDB — バッチメタデータ取得。

    Parameters:
        pdb_ids: list[str] — PDB ID リスト
    """
    results = []
    for pid in pdb_ids:
        try:
            meta = rcsb_get_entry(pid)
            results.append(meta)
        except Exception as e:
            print(f"  Warning: {pid} — {e}")
            continue

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("resolution")
    print(f"RCSB batch: {len(df)} entries")
    return df
```

## 3. リガンド・結合サイト情報

```python
def rcsb_get_ligands(pdb_id):
    """
    RCSB PDB — リガンド情報取得。

    Parameters:
        pdb_id: str — PDB ID
    """
    url = (f"https://data.rcsb.org/rest/v1/core/"
           f"nonpolymer_entity/{pdb_id}")
    # まずエントリのnonpolymerエンティティを取得
    entry = rcsb_get_entry(pdb_id)
    n_ligands = entry.get("nonpolymer_count", 0)

    ligands = []
    for i in range(1, n_ligands + 1):
        try:
            lig_url = (f"https://data.rcsb.org/rest/v1/core/"
                       f"nonpolymer_entity/{pdb_id}/{i}")
            r = requests.get(lig_url, timeout=15)
            if r.status_code == 200:
                ld = r.json()
                comp_id = ld.get(
                    "pdbx_entity_nonpoly", {}).get(
                    "comp_id", "")
                ligands.append({
                    "pdb_id": pdb_id,
                    "entity_id": i,
                    "comp_id": comp_id,
                    "name": ld.get(
                        "rcsb_nonpolymer_entity", {}).get(
                        "pdbx_description", ""),
                    "formula": ld.get(
                        "rcsb_nonpolymer_entity", {}).get(
                        "formula_weight", ""),
                })
        except Exception:
            continue

    df = pd.DataFrame(ligands)
    print(f"RCSB ligands: {pdb_id} → {len(df)} ligands")
    return df
```

## 4. RCSB PDB 統合パイプライン

```python
def rcsb_pipeline(query, resolution_max=3.0,
                     output_dir="results"):
    """
    RCSB PDB 統合パイプライン。

    Parameters:
        query: str — 検索クエリ
        resolution_max: float — 最大解像度
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) テキスト検索
    pdb_ids = rcsb_text_search(
        query, resolution_max=resolution_max)

    # 2) バッチメタデータ
    metadata = rcsb_batch_metadata(pdb_ids[:20])
    metadata.to_csv(output_dir / "pdb_entries.csv",
                    index=False)

    # 3) トップ構造のリガンド
    if not metadata.empty:
        top = metadata.iloc[0]["pdb_id"]
        ligands = rcsb_get_ligands(top)
        ligands.to_csv(output_dir / "ligands.csv",
                       index=False)

    print(f"RCSB pipeline: {output_dir}")
    return {"entries": metadata}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `rcsb_pdb` | RCSB PDB Data | エントリデータ・構造メタデータ |
| `rcsb_search` | RCSB PDB Search | テキスト/配列/構造類似検索 |

## パイプライン統合

```
protein-structure-analysis → rcsb-pdb-search → molecular-docking
  (PDB/AlphaFold 構造)     (RCSB Search API)   (Vina/DiffDock)
          │                       │                    ↓
  uniprot-proteome ──────────────┘           drug-target-profiling
  (UniProt配列)             │               (標的プロファイリング)
                            ↓
                  alphafold-structures
                  (AlphaFold DB)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/pdb_entries.csv` | 構造メタデータ | → protein-structure-analysis |
| `results/ligands.csv` | リガンド情報 | → molecular-docking |
