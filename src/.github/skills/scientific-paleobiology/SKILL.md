---
name: scientific-paleobiology
description: |
  古生物学データベーススキル。Paleobiology Database (PBDB) REST
  API による化石産出記録・分類群・コレクション検索、地質年代
  多様性曲線・古地理解析。ToolUniverse 連携: paleobiology。
tu_tools:
  - key: paleobiology
    name: Paleobiology Database
    description: PBDB 化石産出記録・分類群・コレクション検索
---

# Scientific Paleobiology

Paleobiology Database (PBDB) REST API を活用した古生物学的
多様性解析パイプラインを提供する。

## When to Use

- 化石産出記録 (occurrence) を検索するとき
- 分類群 (taxa) の地質年代分布を調べるとき
- 化石コレクション/産地情報を検索するとき
- 地質年代を通じた多様性曲線を作成するとき
- 大量絶滅イベントのパターンを分析するとき
- 古地理的分布を解析するとき

---

## Quick Start

## 1. PBDB 化石産出記録検索

```python
import requests
import pandas as pd
import numpy as np

PBDB_BASE = "https://paleobiodb.org/data1.2"


def pbdb_search_occurrences(taxon=None, interval=None,
                              lngmin=None, lngmax=None,
                              latmin=None, latmax=None, limit=1000):
    """
    PBDB — 化石産出記録検索。

    Parameters:
        taxon: str — 分類群名 (例: "Dinosauria", "Trilobita")
        interval: str — 地質年代区間 (例: "Cretaceous", "Permian")
        lngmin: float — 経度最小値
        lngmax: float — 経度最大値
        latmin: float — 緯度最小値
        latmax: float — 緯度最大値
        limit: int — 最大結果数
    """
    url = f"{PBDB_BASE}/occs/list.json"
    params = {
        "show": "coords,phylo,time",
        "limit": limit,
    }
    if taxon:
        params["base_name"] = taxon
    if interval:
        params["interval"] = interval
    if all(v is not None for v in [lngmin, lngmax, latmin, latmax]):
        params.update({
            "lngmin": lngmin, "lngmax": lngmax,
            "latmin": latmin, "latmax": latmax,
        })

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    records = resp.json().get("records", [])

    results = []
    for r in records:
        results.append({
            "occurrence_no": r.get("oid", ""),
            "taxon_name": r.get("tna", ""),
            "taxon_rank": r.get("rnk", ""),
            "phylum": r.get("phl", ""),
            "class": r.get("cll", ""),
            "order": r.get("odl", ""),
            "family": r.get("fml", ""),
            "early_interval": r.get("oei", ""),
            "late_interval": r.get("oli", ""),
            "max_ma": r.get("eag", None),
            "min_ma": r.get("lag", None),
            "lng": r.get("lng", None),
            "lat": r.get("lat", None),
            "collection_no": r.get("cid", ""),
            "reference_no": r.get("rid", ""),
        })

    df = pd.DataFrame(results)
    print(f"PBDB occurrences: {len(df)} records "
          f"(taxon={taxon}, interval={interval})")
    return df
```

## 2. PBDB 分類群情報検索

```python
def pbdb_search_taxa(name=None, rank=None, interval=None, limit=500):
    """
    PBDB — 分類群検索。

    Parameters:
        name: str — 分類群名 (例: "Dinosauria")
        rank: str — ランク (例: "genus", "family", "order")
        interval: str — 地質年代区間
        limit: int — 最大結果数
    """
    url = f"{PBDB_BASE}/taxa/list.json"
    params = {
        "show": "attr,app,size",
        "limit": limit,
    }
    if name:
        params["base_name"] = name
    if rank:
        params["rank"] = rank
    if interval:
        params["interval"] = interval

    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    records = resp.json().get("records", [])

    results = []
    for r in records:
        results.append({
            "taxon_no": r.get("oid", ""),
            "taxon_name": r.get("nam", ""),
            "rank": r.get("rnk", ""),
            "parent_name": r.get("prl", ""),
            "n_occs": r.get("noc", 0),
            "first_appearance": r.get("fea", ""),
            "last_appearance": r.get("lla", ""),
            "extant": r.get("ext", ""),
        })

    df = pd.DataFrame(results)
    print(f"PBDB taxa: {len(df)} records (name={name})")
    return df
```

## 3. 地質年代多様性曲線

```python
def pbdb_diversity_curve(taxon, time_resolution="stage",
                           rank="genus"):
    """
    PBDB — 地質年代多様性曲線生成。

    Parameters:
        taxon: str — 分類群名
        time_resolution: str — "stage" or "epoch" or "period"
        rank: str — カウントするランク ("genus", "family")
    """
    url = f"{PBDB_BASE}/occs/diversity.json"
    params = {
        "base_name": taxon,
        "count": rank,
        "time_reso": time_resolution,
    }
    resp = requests.get(url, params=params, timeout=60)
    resp.raise_for_status()
    records = resp.json().get("records", [])

    results = []
    for r in records:
        results.append({
            "interval_name": r.get("idn", ""),
            "max_ma": r.get("eag", None),
            "min_ma": r.get("lag", None),
            "mid_ma": (float(r.get("eag", 0)) +
                       float(r.get("lag", 0))) / 2,
            "sampled_in_bin": r.get("dsb", 0),
            "n_originations": r.get("dor", 0),
            "n_extinctions": r.get("dex", 0),
            "range_through": r.get("drt", 0),
        })

    df = pd.DataFrame(results)
    print(f"PBDB diversity: {len(df)} intervals, "
          f"max diversity={df['sampled_in_bin'].max()} {rank}")
    return df
```

## 4. 古生物学統合パイプライン

```python
def paleobiology_pipeline(taxon, interval=None,
                            output_dir="results"):
    """
    古生物学統合パイプライン。

    Parameters:
        taxon: str — 分類群名 (例: "Dinosauria")
        interval: str — 地質年代区間 (オプション)
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 産出記録
    occ = pbdb_search_occurrences(taxon=taxon, interval=interval)
    occ.to_csv(output_dir / "occurrences.csv", index=False)

    # 2) 分類群情報
    taxa = pbdb_search_taxa(name=taxon)
    taxa.to_csv(output_dir / "taxa.csv", index=False)

    # 3) 多様性曲線
    diversity = pbdb_diversity_curve(taxon)
    diversity.to_csv(output_dir / "diversity.csv", index=False)

    # 4) 地理的サマリ
    if "lat" in occ.columns and "lng" in occ.columns:
        geo_summary = occ.groupby("early_interval").agg(
            n_records=("occurrence_no", "count"),
            mean_lat=("lat", "mean"),
            mean_lng=("lng", "mean"),
        ).reset_index()
        geo_summary.to_csv(output_dir / "geo_summary.csv", index=False)

    print(f"Paleobiology pipeline: {output_dir}")
    return {
        "occurrences": occ,
        "taxa": taxa,
        "diversity": diversity,
    }
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `paleobiology` | Paleobiology Database | 化石産出・分類群・コレクション検索 |

## パイプライン統合

```
phylogenetics → paleobiology → environmental-ecology
  (系統解析)   (化石記録)      (GBIF/生態)
      │              │                ↓
  taxonomy ─────────┘         environmental-geodata
  (分類体系)    │              (環境モデリング)
                ↓
         macroevolution
         (大進化パターン)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/occurrences.csv` | 化石産出記録 | → environmental-ecology |
| `results/taxa.csv` | 分類群情報 | → phylogenetics |
| `results/diversity.csv` | 多様性曲線 | → macroevolution |
| `results/geo_summary.csv` | 古地理サマリ | → environmental-geodata |
