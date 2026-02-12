---
name: scientific-compound-screening
description: |
  化合物スクリーニングスキル。ZINC データベースを活用した購入可能化合物検索、
  SMILES/名前ベースの類似性検索、カタログフィルタリング、
  バーチャルスクリーニング前処理パイプライン。
---

# Scientific Compound Screening

ZINC データベースを活用した化合物ライブラリ検索・
バーチャルスクリーニング前処理パイプラインを提供する。

## When to Use

- 購入可能な化合物ライブラリを検索するとき
- SMILES 構造式から類似化合物を探すとき
- 化合物名からデータベースレコードを取得するとき
- ベンダーカタログの絞り込みを行うとき
- バーチャルスクリーニング用の化合物セットを準備するとき

---

## Quick Start

## 1. ZINC 化合物名検索

```python
import requests
import pandas as pd

ZINC_API = "https://zinc15.docking.org"


def zinc_search_by_name(name, max_results=20):
    """
    ZINC データベースで化合物名による検索。

    Parameters:
        name: str — compound name (e.g., "aspirin")
        max_results: int — maximum results

    ToolUniverse:
        ZINC_search_by_name(name=name)
    """
    url = f"{ZINC_API}/substances/search"
    params = {"q": name, "count": max_results}
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for item in data:
        results.append({
            "zinc_id": item.get("zinc_id", ""),
            "name": item.get("name", ""),
            "smiles": item.get("smiles", ""),
            "mwt": item.get("mwt", ""),
            "logp": item.get("logp", ""),
            "purchasable": item.get("purchasability", ""),
        })

    df = pd.DataFrame(results)
    print(f"ZINC search '{name}': {len(df)} compounds")
    return df
```

## 2. ZINC SMILES 類似性検索

```python
def zinc_search_by_smiles(smiles, similarity=0.7, max_results=20):
    """
    ZINC で SMILES 構造式による類似性検索。

    Parameters:
        smiles: str — SMILES string
        similarity: float — Tanimoto similarity threshold (0-1)

    ToolUniverse:
        ZINC_search_by_smiles(smiles=smiles)
    """
    url = f"{ZINC_API}/substances/search"
    params = {
        "smiles": smiles,
        "similarity": similarity,
        "count": max_results,
    }
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for item in data:
        results.append({
            "zinc_id": item.get("zinc_id", ""),
            "smiles": item.get("smiles", ""),
            "similarity": item.get("similarity", ""),
            "mwt": item.get("mwt", ""),
            "logp": item.get("logp", ""),
            "purchasable": item.get("purchasability", ""),
        })

    df = pd.DataFrame(results)
    print(f"ZINC SMILES search: {len(df)} similar compounds "
          f"(threshold={similarity})")
    return df
```

## 3. ZINC 化合物詳細取得

```python
def zinc_get_substance(zinc_id):
    """
    ZINC ID から化合物の完全情報を取得。

    Parameters:
        zinc_id: str — ZINC ID (e.g., "ZINC000000000001")

    ToolUniverse:
        ZINC_get_substance(zinc_id=zinc_id)
    """
    url = f"{ZINC_API}/substances/{zinc_id}.json"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    info = {
        "zinc_id": data.get("zinc_id", ""),
        "name": data.get("name", ""),
        "smiles": data.get("smiles", ""),
        "inchikey": data.get("inchikey", ""),
        "mwt": data.get("mwt", ""),
        "logp": data.get("logp", ""),
        "num_rotatable_bonds": data.get("num_rotatable_bonds", ""),
        "num_hba": data.get("num_hba", ""),
        "num_hbd": data.get("num_hbd", ""),
        "tpsa": data.get("tpsa", ""),
        "purchasable": data.get("purchasability", ""),
    }

    print(f"ZINC {zinc_id}: {info['name']} (MW={info['mwt']})")
    return info, data
```

## 4. ZINC カタログ一覧

```python
def zinc_get_catalogs():
    """
    ZINC の利用可能カタログ (ベンダー) 一覧を取得。

    ToolUniverse:
        ZINC_get_catalogs()
    """
    url = f"{ZINC_API}/catalogs.json"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for cat in data:
        results.append({
            "catalog_name": cat.get("name", ""),
            "short_name": cat.get("short_name", ""),
            "num_substances": cat.get("num_substances", 0),
            "url": cat.get("url", ""),
        })

    df = pd.DataFrame(results)
    print(f"ZINC catalogs: {len(df)} vendors")
    return df
```

## 5. バーチャルスクリーニング前処理パイプライン

```python
def virtual_screening_prep(query_smiles, lipinski=True, max_compounds=100):
    """
    バーチャルスクリーニング用の化合物セット準備。
    Lipinski's Rule of Five フィルタリング含む。

    ToolUniverse (横断):
        ZINC_search_by_smiles(smiles=query_smiles) → ZINC_get_substance(zinc_id)
    """
    # Step 1: Similar compound search
    df = zinc_search_by_smiles(query_smiles, similarity=0.6,
                               max_results=max_compounds)

    if df.empty:
        print("No similar compounds found")
        return df

    # Step 2: Lipinski filter
    if lipinski:
        df["mwt"] = pd.to_numeric(df["mwt"], errors="coerce")
        df["logp"] = pd.to_numeric(df["logp"], errors="coerce")
        before = len(df)
        df = df[
            (df["mwt"] <= 500)
            & (df["logp"] <= 5)
        ]
        print(f"Lipinski filter: {before} → {len(df)} compounds")

    # Step 3: Sort by similarity
    df["similarity"] = pd.to_numeric(df["similarity"], errors="coerce")
    df = df.sort_values("similarity", ascending=False)

    print(f"VS prep: {len(df)} compounds ready for screening")
    return df
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/zinc_search.csv` | CSV |
| `results/zinc_similar.csv` | CSV |
| `results/zinc_substance.json` | JSON |
| `results/zinc_catalogs.csv` | CSV |
| `results/vs_library.csv` | CSV |

### 利用可能ツール

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| ZINC | `ZINC_search_by_name` | 化合物名検索 |
| ZINC | `ZINC_search_by_smiles` | SMILES 類似性検索 |
| ZINC | `ZINC_get_substance` | 化合物詳細 |
| ZINC | `ZINC_get_catalogs` | カタログ一覧 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-compound-similarity` | 化合物類似性 |
| `scientific-pharmacology-targets` | 薬理学ターゲット |
| `scientific-molecular-docking` | 分子ドッキング |
| `scientific-drug-target-interaction` | DTI 解析 |
| `scientific-admet-toxicity` | ADMET 毒性 |

### 依存パッケージ

`requests`, `pandas`
