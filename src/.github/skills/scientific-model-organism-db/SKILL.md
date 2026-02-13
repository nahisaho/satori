---
name: scientific-model-organism-db
description: |
  モデル生物データベース統合スキル。FlyBase (ショウジョウバエ)、
  WormBase (線虫)、ZFIN (ゼブラフィッシュ)、RGD (ラット)、
  MGI (マウス) の REST API を統合した
  モデル生物遺伝子・表現型・疾患モデル横断検索パイプライン。
  ToolUniverse 連携: impc, mpd。
tu_tools:
  - key: impc
    name: IMPC
    description: 国際マウス表現型解析コンソーシアム
  - key: mpd
    name: MPD
    description: Mouse Phenome Database マウス表現型
---

# Scientific Model Organism Database

主要 5 モデル生物データベース (FlyBase / WormBase / ZFIN / RGD / MGI) を
統合した遺伝子・表現型・疾患モデル横断検索パイプラインを提供する。

## When to Use

- ヒト遺伝子のモデル生物オルソログを検索するとき
- モデル生物の表現型データを疾患研究に活用するとき
- 遺伝子改変動物モデルの表現型情報を取得するとき
- 複数のモデル生物間で機能保存性を比較するとき
- IMPC (既存スキル) を補完してラット/魚/ハエ/線虫データが必要なとき

---

## Quick Start

## 1. MGI (Mouse Genome Informatics) 遺伝子検索

```python
import requests
import pandas as pd

MGI_API = "http://www.informatics.jax.org/api"


def search_mgi_gene(query, limit=20):
    """
    MGI マウス遺伝子検索。

    Parameters:
        query: str — 遺伝子名またはシンボル
        limit: int — 最大取得数
    """
    url = f"{MGI_API}/gene/search"
    params = {"query": query, "limit": limit}
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for gene in data.get("results", []):
        rows.append({
            "mgi_id": gene.get("mgiId"),
            "symbol": gene.get("symbol"),
            "name": gene.get("name"),
            "chromosome": gene.get("chromosome"),
            "feature_type": gene.get("featureType"),
            "organism": "Mus musculus",
        })

    df = pd.DataFrame(rows[:limit])
    print(f"MGI search '{query}': {len(df)} genes")
    return df
```

## 2. RGD (Rat Genome Database) 遺伝子検索

```python
RGD_API = "https://rest.rgd.mcw.edu/rgdws"


def search_rgd_gene(query, species="rat"):
    """
    RGD ラット遺伝子検索。

    Parameters:
        query: str — 遺伝子シンボル
        species: str — "rat", "mouse", "human"
    """
    species_map = {"rat": 3, "mouse": 2, "human": 1}
    species_key = species_map.get(species, 3)

    url = f"{RGD_API}/genes/{query}/{species_key}"
    resp = requests.get(url, headers={"Accept": "application/json"})
    resp.raise_for_status()
    data = resp.json()

    if isinstance(data, dict):
        data = [data]

    rows = []
    for gene in data:
        rows.append({
            "rgd_id": gene.get("rgdId"),
            "symbol": gene.get("symbol"),
            "name": gene.get("name"),
            "chromosome": gene.get("chromosome"),
            "type": gene.get("type"),
            "organism": species,
        })

    df = pd.DataFrame(rows)
    print(f"RGD search '{query}': {len(df)} genes ({species})")
    return df
```

## 3. ZFIN (Zebrafish Information Network)

```python
ZFIN_API = "https://zfin.org/action/api"


def search_zfin_gene(query, limit=20):
    """
    ZFIN ゼブラフィッシュ遺伝子検索。

    Parameters:
        query: str — 遺伝子名またはシンボル
        limit: int — 最大取得数
    """
    url = f"{ZFIN_API}/marker/search"
    params = {"name": query, "limit": limit, "type": "gene"}
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for gene in data.get("results", []):
        rows.append({
            "zfin_id": gene.get("zdbID"),
            "symbol": gene.get("abbreviation"),
            "name": gene.get("name"),
            "type": gene.get("markerType"),
            "organism": "Danio rerio",
        })

    df = pd.DataFrame(rows[:limit])
    print(f"ZFIN search '{query}': {len(df)} genes")
    return df
```

## 4. FlyBase (Drosophila)

```python
FLYBASE_API = "https://api.flybase.org/api/v1.0"


def search_flybase_gene(query, limit=20):
    """
    FlyBase ショウジョウバエ遺伝子検索。

    Parameters:
        query: str — 遺伝子名またはシンボル
        limit: int — 最大取得数
    """
    url = f"{FLYBASE_API}/gene/search/{query}"
    params = {"limit": limit}
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    genes = data.get("resultset", {}).get("result", [])
    rows = []
    for gene in genes:
        rows.append({
            "flybase_id": gene.get("id"),
            "symbol": gene.get("symbol"),
            "name": gene.get("name"),
            "chromosome": gene.get("location", {}).get("chromosome"),
            "organism": "Drosophila melanogaster",
        })

    df = pd.DataFrame(rows[:limit])
    print(f"FlyBase search '{query}': {len(df)} genes")
    return df
```

## 5. WormBase (C. elegans)

```python
WORMBASE_API = "https://wormbase.org/rest"


def search_wormbase_gene(query, limit=20):
    """
    WormBase 線虫遺伝子検索。

    Parameters:
        query: str — 遺伝子名またはシンボル
        limit: int — 最大取得数
    """
    url = f"{WORMBASE_API}/field/gene/{query}/overview"
    resp = requests.get(url, headers={"Accept": "application/json"})
    resp.raise_for_status()
    data = resp.json()

    overview = data.get("overview", {})
    info = {
        "wormbase_id": overview.get("name", {}).get("data", {}).get("id"),
        "symbol": overview.get("name", {}).get("data", {}).get("label"),
        "concise_description": overview.get("concise_description", {}).get(
            "data", ""
        ),
        "organism": "Caenorhabditis elegans",
    }

    print(f"WormBase: {info['symbol']} ({info['wormbase_id']})")
    return info
```

## 6. モデル生物横断オルソログ検索

```python
def cross_species_ortholog_search(human_gene):
    """
    ヒト遺伝子の 5 モデル生物オルソログ横断検索。

    Parameters:
        human_gene: str — ヒト遺伝子シンボル (例: "TP53")

    Pipeline:
        MGI → RGD → ZFIN → FlyBase → WormBase
    """
    results = []

    # Mouse (MGI)
    try:
        mgi = search_mgi_gene(human_gene, limit=3)
        if not mgi.empty:
            results.append({"organism": "Mouse", "db": "MGI",
                            "symbol": mgi.iloc[0]["symbol"],
                            "id": mgi.iloc[0]["mgi_id"]})
    except Exception as e:
        print(f"MGI error: {e}")

    # Rat (RGD)
    try:
        rgd = search_rgd_gene(human_gene, "rat")
        if not rgd.empty:
            results.append({"organism": "Rat", "db": "RGD",
                            "symbol": rgd.iloc[0]["symbol"],
                            "id": str(rgd.iloc[0]["rgd_id"])})
    except Exception as e:
        print(f"RGD error: {e}")

    # Zebrafish (ZFIN)
    try:
        zfin = search_zfin_gene(human_gene.lower(), limit=3)
        if not zfin.empty:
            results.append({"organism": "Zebrafish", "db": "ZFIN",
                            "symbol": zfin.iloc[0]["symbol"],
                            "id": zfin.iloc[0]["zfin_id"]})
    except Exception as e:
        print(f"ZFIN error: {e}")

    # Drosophila (FlyBase)
    try:
        fly = search_flybase_gene(human_gene, limit=3)
        if not fly.empty:
            results.append({"organism": "Drosophila", "db": "FlyBase",
                            "symbol": fly.iloc[0]["symbol"],
                            "id": fly.iloc[0]["flybase_id"]})
    except Exception as e:
        print(f"FlyBase error: {e}")

    # C. elegans (WormBase)
    try:
        worm = search_wormbase_gene(human_gene.lower())
        if worm.get("wormbase_id"):
            results.append({"organism": "C. elegans", "db": "WormBase",
                            "symbol": worm["symbol"],
                            "id": worm["wormbase_id"]})
    except Exception as e:
        print(f"WormBase error: {e}")

    df = pd.DataFrame(results)
    print(f"\nOrthologs of {human_gene}: {len(df)} model organisms")
    return df
```

## 7. 表現型データ統合

```python
def get_mgi_phenotypes(mgi_id):
    """
    MGI 遺伝子の表現型アノテーション取得。

    Parameters:
        mgi_id: str — MGI ID (例: "MGI:98834")
    """
    url = f"{MGI_API}/gene/{mgi_id}/phenotypes"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for pheno in data.get("phenotypes", []):
        rows.append({
            "mp_id": pheno.get("mpId"),
            "mp_term": pheno.get("mpTerm"),
            "allele_symbol": pheno.get("alleleSymbol"),
            "genetic_background": pheno.get("geneticBackground"),
        })

    df = pd.DataFrame(rows)
    print(f"Phenotypes for {mgi_id}: {len(df)} MP annotations")
    return df
```

---

## パイプライン統合

```
ensembl-genomics ──→ model-organism-db ──→ disease-research
  (オルソログ ID)      (表現型データ)        (疾患モデル)
        │                     │                   ↓
biothings-idmapping ──┘      ↓          rare-disease-genetics
  (ID マッピング)     phylogenetics        (OMIM/Orphanet)
                    (種間系統解析)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/model_orthologs.csv` | モデル生物オルソログ | → ensembl-genomics |
| `results/mgi_phenotypes.csv` | マウス表現型 | → disease-research |
| `results/cross_species.json` | 横断比較結果 | → phylogenetics |
