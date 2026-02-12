---
name: scientific-reactome-pathways
description: |
  Reactome パスウェイスキル。Reactome Content Service
  REST API によるパスウェイ検索・階層取得・UniProt マッピング・
  パスウェイ図データ取得。ToolUniverse 連携: reactome。
tu_tools:
  - key: reactome
    name: Reactome
    description: パスウェイデータベース REST API
---

# Scientific Reactome Pathways

Reactome Content Service REST API を活用したパスウェイ検索・
階層構造取得・UniProt アクセッション→パスウェイマッピング・
パスウェイ図データ取得パイプラインを提供する。

## When to Use

- 生物学パスウェイを名前やキーワードで検索するとき
- パスウェイ階層ツリーを取得するとき
- UniProt アクセッションからパスウェイをマッピングするとき
- パスウェイの参加者 (タンパク質/化合物) を列挙するとき
- パスウェイ図のレイアウトデータを取得するとき

---

## Quick Start

## 1. パスウェイ検索・詳細取得

```python
import requests
import pandas as pd

REACTOME = "https://reactome.org/ContentService"


def reactome_search(query, species="Homo sapiens",
                       limit=25):
    """
    Reactome — パスウェイ検索。

    Parameters:
        query: str — 検索クエリ (例: "apoptosis", "MAPK")
        species: str — 種名 (例: "Homo sapiens")
        limit: int — 最大結果数
    """
    url = f"{REACTOME}/search/query"
    params = {
        "query": query,
        "species": species,
        "types": "Pathway",
        "cluster": "true",
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for group in data.get("results", []):
        for entry in group.get("entries", [])[:limit]:
            rows.append({
                "stId": entry.get("stId", ""),
                "name": entry.get("name", ""),
                "species": entry.get("species", ""),
                "exact_type": entry.get(
                    "exactType", ""),
                "compartments": "; ".join(
                    entry.get("compartmentNames", [])),
            })

    df = pd.DataFrame(rows[:limit])
    print(f"Reactome search: '{query}' → {len(df)} "
          f"pathways")
    return df


def reactome_pathway_detail(pathway_id):
    """
    Reactome — パスウェイ詳細取得。

    Parameters:
        pathway_id: str — Reactome Stable ID
            (例: "R-HSA-109581")
    """
    url = f"{REACTOME}/data/pathway/{pathway_id}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    result = {
        "stId": data.get("stId", ""),
        "name": data.get("displayName", ""),
        "species": data.get("speciesName", ""),
        "is_inferred": data.get("isInferred", False),
        "has_diagram": data.get("hasDiagram", False),
        "n_sub_events": len(
            data.get("hasEvent", [])),
        "n_compartments": len(
            data.get("compartment", [])),
        "release_date": data.get("releaseDate", ""),
    }
    return result
```

## 2. UniProt→パスウェイ マッピング

```python
def reactome_uniprot_pathways(uniprot_id,
                                  species="Homo sapiens"):
    """
    Reactome — UniProt → パスウェイマッピング。

    Parameters:
        uniprot_id: str — UniProt アクセッション
            (例: "P38398" = BRCA1)
        species: str — 種名
    """
    url = f"{REACTOME}/data/pathways/low/entity/{uniprot_id}"
    params = {"species": species}
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for pw in data:
        rows.append({
            "pathway_id": pw.get("stId", ""),
            "pathway_name": pw.get("displayName", ""),
            "species": pw.get("speciesName", ""),
            "has_diagram": pw.get("hasDiagram", False),
        })

    df = pd.DataFrame(rows)
    print(f"Reactome UniProt→pathway: {uniprot_id} "
          f"→ {len(df)} pathways")
    return df
```

## 3. パスウェイ参加者取得

```python
def reactome_participants(pathway_id):
    """
    Reactome — パスウェイ参加者一覧。

    Parameters:
        pathway_id: str — Reactome Stable ID
    """
    url = (f"{REACTOME}/data/participants/"
           f"{pathway_id}")
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for item in data:
        pe_name = item.get("displayName", "")
        for ref in item.get("refEntities", []):
            rows.append({
                "pathway_id": pathway_id,
                "participant": pe_name,
                "db_name": ref.get("databaseName", ""),
                "identifier": ref.get("identifier", ""),
                "name": ref.get("displayName", ""),
            })

    df = pd.DataFrame(rows)
    print(f"Reactome participants: {pathway_id} "
          f"→ {len(df)} entities")
    return df
```

## 4. Reactome 統合パイプライン

```python
def reactome_pipeline(query_or_uniprot,
                         output_dir="results"):
    """
    Reactome 統合パイプライン。

    Parameters:
        query_or_uniprot: str — 検索クエリ or UniProt ID
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    is_uniprot = (len(query_or_uniprot) == 6
                  and query_or_uniprot[0].isalpha())

    if is_uniprot:
        # UniProt → パスウェイ
        pathways = reactome_uniprot_pathways(
            query_or_uniprot)
    else:
        # テキスト検索
        pathways = reactome_search(query_or_uniprot)

    pathways.to_csv(output_dir / "reactome_pathways.csv",
                    index=False)

    # トップパスウェイの参加者
    if not pathways.empty:
        top_id = (pathways.iloc[0].get("pathway_id")
                  or pathways.iloc[0].get("stId"))
        parts = reactome_participants(top_id)
        parts.to_csv(
            output_dir / "reactome_participants.csv",
            index=False)

    print(f"Reactome pipeline: {output_dir}")
    return {"pathways": pathways}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `reactome` | Reactome | パスウェイデータベース REST API |

## パイプライン統合

```
pathway-enrichment → reactome-pathways → systems-biology
  (GO/パスウェイ)     (Reactome API)     (ネットワーク解析)
        │                   │                 ↓
uniprot-proteome ──────────┘       metabolomics-databases
  (UniProt ID)                     (MetaCyc 代謝)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/reactome_pathways.csv` | パスウェイ一覧 | → pathway-enrichment |
| `results/reactome_participants.csv` | 参加者 | → protein-interaction-network |
