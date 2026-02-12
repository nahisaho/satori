---
name: scientific-marine-ecology
description: |
  海洋生態学統合スキル。OBIS 海洋生物分布・WoRMS 海洋分類体系・
  GBIF 生物多様性レコード・FishBase 魚類データ。ToolUniverse
  連携: obis, worms, gbif。
tu_tools:
  - key: obis
    name: OBIS (Ocean Biodiversity Information System)
    description: 海洋生物の出現・分布データを提供する国際プラットフォーム
  - key: worms
    name: WoRMS (World Register of Marine Species)
    description: 海洋生物種の権威ある分類学的参照データベース
  - key: gbif
    name: GBIF (Global Biodiversity Information Facility)
    description: 生物多様性データの国際的なオープンアクセスインフラ
---

# Scientific Marine Ecology

OBIS / WoRMS / GBIF / FishBase を活用した海洋生物多様性・
分布解析パイプラインを提供する。

## When to Use

- 海洋生物の地理的分布データを取得するとき
- 海洋生物の分類学的情報 (WoRMS) を検証するとき
- 生物多様性ホットスポットを解析するとき
- 魚類の生態・形態データ (FishBase) を取得するとき
- 海洋保全区域の生物多様性評価を行うとき
- 海洋環境変動と種分布の関係を分析するとき

---

## Quick Start

## 1. OBIS 海洋生物分布

```python
import requests
import pandas as pd
import numpy as np

OBIS_BASE = "https://api.obis.org/v3"


def obis_occurrence_search(taxon_name=None, taxon_id=None,
                            geometry=None, year_range=None, limit=1000):
    """
    OBIS — 海洋生物出現記録検索。

    Parameters:
        taxon_name: str — 学名 (例: "Delphinidae")
        taxon_id: int — AphiaID
        geometry: str — WKT ジオメトリ (例: "POLYGON((...))")
        year_range: tuple — (start_year, end_year)
        limit: int — 最大取得件数
    """
    url = f"{OBIS_BASE}/occurrence"
    params = {"size": min(limit, 5000)}

    if taxon_name:
        params["scientificname"] = taxon_name
    if taxon_id:
        params["taxonid"] = taxon_id
    if geometry:
        params["geometry"] = geometry
    if year_range:
        params["startdate"] = f"{year_range[0]}-01-01"
        params["enddate"] = f"{year_range[1]}-12-31"

    resp = requests.get(url, params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    records = []
    for rec in data.get("results", []):
        records.append({
            "scientific_name": rec.get("scientificName", ""),
            "aphia_id": rec.get("aphiaID", ""),
            "latitude": rec.get("decimalLatitude", None),
            "longitude": rec.get("decimalLongitude", None),
            "depth": rec.get("depth", None),
            "date": rec.get("date_mid", ""),
            "dataset_id": rec.get("dataset_id", ""),
            "basis_of_record": rec.get("basisOfRecord", ""),
        })

    df = pd.DataFrame(records)
    print(f"OBIS: '{taxon_name or taxon_id}' → {len(df)} occurrences")
    return df


def obis_checklist(geometry=None, area_id=None):
    """
    OBIS — 地域別種チェックリスト。

    Parameters:
        geometry: str — WKT ジオメトリ
        area_id: int — OBIS エリア ID
    """
    url = f"{OBIS_BASE}/checklist"
    params = {"size": 5000}
    if geometry:
        params["geometry"] = geometry
    if area_id:
        params["areaid"] = area_id

    resp = requests.get(url, params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    species = []
    for sp in data.get("results", []):
        species.append({
            "scientific_name": sp.get("scientificName", ""),
            "aphia_id": sp.get("taxonID", ""),
            "records": sp.get("records", 0),
            "kingdom": sp.get("kingdom", ""),
            "phylum": sp.get("phylum", ""),
            "class": sp.get("class", ""),
            "order": sp.get("order", ""),
            "family": sp.get("family", ""),
        })

    df = pd.DataFrame(species)
    print(f"OBIS checklist: {len(df)} species")
    return df
```

## 2. WoRMS 分類学検索

```python
WORMS_BASE = "https://www.marinespecies.org/rest"


def worms_taxon_search(name, fuzzy=True, marine_only=True):
    """
    WoRMS — 海洋生物分類学的検索。

    Parameters:
        name: str — 種名/属名
        fuzzy: bool — ファジー検索
        marine_only: bool — 海洋種のみ
    """
    url = f"{WORMS_BASE}/AphiaRecordsByName/{name}"
    params = {
        "like": str(fuzzy).lower(),
        "marine_only": str(marine_only).lower(),
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for taxon in data if isinstance(data, list) else [data]:
        results.append({
            "aphia_id": taxon.get("AphiaID", ""),
            "scientific_name": taxon.get("scientificname", ""),
            "authority": taxon.get("authority", ""),
            "status": taxon.get("status", ""),
            "rank": taxon.get("rank", ""),
            "valid_name": taxon.get("valid_name", ""),
            "kingdom": taxon.get("kingdom", ""),
            "phylum": taxon.get("phylum", ""),
            "class": taxon.get("class", ""),
            "order": taxon.get("order", ""),
            "family": taxon.get("family", ""),
            "genus": taxon.get("genus", ""),
            "is_marine": taxon.get("isMarine", 0),
            "is_brackish": taxon.get("isBrackish", 0),
            "is_freshwater": taxon.get("isFreshwater", 0),
        })

    df = pd.DataFrame(results)
    print(f"WoRMS: '{name}' → {len(df)} taxa")
    return df


def worms_classification(aphia_id):
    """
    WoRMS — 完全分類階層取得。

    Parameters:
        aphia_id: int — AphiaID
    """
    url = f"{WORMS_BASE}/AphiaClassificationByAphiaID/{aphia_id}"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    hierarchy = []
    node = data
    while node:
        hierarchy.append({
            "aphia_id": node.get("AphiaID", ""),
            "name": node.get("scientificname", ""),
            "rank": node.get("rank", ""),
        })
        node = node.get("child")

    df = pd.DataFrame(hierarchy)
    print(f"WoRMS classification: {len(df)} levels")
    return df
```

## 3. GBIF 生物多様性レコード

```python
GBIF_BASE = "https://api.gbif.org/v1"


def gbif_species_search(name, limit=20):
    """
    GBIF — 種名検索・分類マッチング。

    Parameters:
        name: str — 種名
        limit: int — 結果上限
    """
    url = f"{GBIF_BASE}/species/search"
    params = {"q": name, "limit": limit}
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for sp in data.get("results", []):
        results.append({
            "taxon_key": sp.get("key", ""),
            "scientific_name": sp.get("scientificName", ""),
            "canonical_name": sp.get("canonicalName", ""),
            "status": sp.get("taxonomicStatus", ""),
            "rank": sp.get("rank", ""),
            "kingdom": sp.get("kingdom", ""),
            "phylum": sp.get("phylum", ""),
            "class": sp.get("class", ""),
            "order": sp.get("order", ""),
            "family": sp.get("family", ""),
            "num_occurrences": sp.get("numOccurrences", 0),
        })

    df = pd.DataFrame(results)
    print(f"GBIF species: '{name}' → {len(df)} taxa")
    return df


def gbif_occurrence_search(taxon_key=None, country=None,
                            year_range=None, limit=300):
    """
    GBIF — 出現記録検索。

    Parameters:
        taxon_key: int — GBIF taxon key
        country: str — ISO 国コード (例: "JP")
        year_range: tuple — (start, end)
        limit: int — 最大件数
    """
    url = f"{GBIF_BASE}/occurrence/search"
    params = {"limit": min(limit, 300)}

    if taxon_key:
        params["taxonKey"] = taxon_key
    if country:
        params["country"] = country
    if year_range:
        params["year"] = f"{year_range[0]},{year_range[1]}"

    resp = requests.get(url, params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()

    records = []
    for rec in data.get("results", []):
        records.append({
            "gbif_id": rec.get("gbifID", ""),
            "scientific_name": rec.get("scientificName", ""),
            "latitude": rec.get("decimalLatitude", None),
            "longitude": rec.get("decimalLongitude", None),
            "country": rec.get("country", ""),
            "year": rec.get("year", ""),
            "basis_of_record": rec.get("basisOfRecord", ""),
            "institution": rec.get("institutionCode", ""),
        })

    df = pd.DataFrame(records)
    print(f"GBIF occurrences: {len(df)} records (total: {data.get('count', 0)})")
    return df
```

## 4. FishBase 魚類データ

```python
FISHBASE_BASE = "https://fishbase.ropensci.org"


def fishbase_species(genus=None, species=None, family=None):
    """
    FishBase — 魚類種データ取得。

    Parameters:
        genus: str — 属名
        species: str — 種小名
        family: str — 科名
    """
    url = f"{FISHBASE_BASE}/species"
    params = {"limit": 100}
    if genus:
        params["Genus"] = genus
    if species:
        params["Species"] = species
    if family:
        params["Family"] = family

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    records = []
    for fish in data.get("data", []):
        records.append({
            "spec_code": fish.get("SpecCode", ""),
            "genus": fish.get("Genus", ""),
            "species": fish.get("Species", ""),
            "family": fish.get("Family", ""),
            "body_shape": fish.get("BodyShapeI", ""),
            "max_length": fish.get("Length", None),
            "vulnerability": fish.get("Vulnerability", None),
            "importance": fish.get("Importance", ""),
            "habitat": fish.get("DemersPelag", ""),
            "depth_range": f"{fish.get('DepthRangeShallow', '')}-{fish.get('DepthRangeDeep', '')}",
        })

    df = pd.DataFrame(records)
    print(f"FishBase: {len(df)} species")
    return df
```

## 5. 海洋生態学統合パイプライン

```python
def marine_ecology_pipeline(taxon_name, region_wkt=None,
                             output_dir="results"):
    """
    OBIS + WoRMS + GBIF + FishBase 統合パイプライン。

    Parameters:
        taxon_name: str — 分類群名
        region_wkt: str — 調査海域 WKT
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) WoRMS 分類検証
    taxonomy = worms_taxon_search(taxon_name)
    taxonomy.to_csv(output_dir / "worms_taxonomy.csv", index=False)

    aphia_id = None
    if len(taxonomy) > 0:
        aphia_id = taxonomy.iloc[0]["aphia_id"]
        classification = worms_classification(aphia_id)
        classification.to_csv(output_dir / "classification.csv", index=False)

    # 2) OBIS 海洋分布
    obis_data = obis_occurrence_search(
        taxon_name=taxon_name,
        taxon_id=aphia_id,
        geometry=region_wkt,
    )
    obis_data.to_csv(output_dir / "obis_occurrences.csv", index=False)

    # 3) GBIF 補完
    gbif_sp = gbif_species_search(taxon_name)
    if len(gbif_sp) > 0:
        taxon_key = gbif_sp.iloc[0]["taxon_key"]
        gbif_occ = gbif_occurrence_search(taxon_key=taxon_key)
        gbif_occ.to_csv(output_dir / "gbif_occurrences.csv", index=False)
    else:
        gbif_occ = pd.DataFrame()

    # 4) 多様性指標
    if len(obis_data) > 0:
        n_species = obis_data["scientific_name"].nunique()
        lat_range = (obis_data["latitude"].min(), obis_data["latitude"].max())
        depth_range = (obis_data["depth"].min(), obis_data["depth"].max())
        print(f"Diversity: {n_species} species, "
              f"lat {lat_range}, depth {depth_range}")

    print(f"Marine ecology pipeline: {output_dir}")
    return {
        "taxonomy": taxonomy,
        "obis_occurrences": obis_data,
        "gbif_occurrences": gbif_occ,
    }
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `obis` | OBIS | 海洋生物出現・分布レコード検索 |
| `worms` | WoRMS | 海洋種分類学的検証・階層取得 |
| `gbif` | GBIF | 生物多様性出現レコード検索 |

## パイプライン統合

```
environmental-ecology → marine-ecology → phylogenetics
  (陸域生態学)          (OBIS/WoRMS/GBIF)  (系統解析)
       │                      │               ↓
  biodiversity-db ───────────┘         species-distribution
  (ENA/BOLD)           │               (空間分布モデル)
                       ↓
                  fisheries-management
                  (水産資源管理)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/worms_taxonomy.csv` | WoRMS 分類情報 | → phylogenetics |
| `results/obis_occurrences.csv` | OBIS 出現記録 | → species-distribution |
| `results/gbif_occurrences.csv` | GBIF 出現記録 | → biodiversity-db |
| `results/classification.csv` | 分類階層 | → environmental-ecology |
