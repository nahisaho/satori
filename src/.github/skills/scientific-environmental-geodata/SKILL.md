---
name: scientific-environmental-geodata
description: |
  環境地理空間データスキル。SoilGrids REST API による土壌特性
  取得、WorldClim/CHELSA 気候データ、生物多様性-環境モデリング
  統合。直接 REST API 連携 (TU 外)。
tu_tools: []
---

# Scientific Environmental Geodata

SoilGrids・WorldClim 等の地球観測/環境データ API を活用した
生態学的環境モデリングパイプラインを提供する。

## When to Use

- グローバル土壌特性 (pH, SOC, 粘土含量) を取得するとき
- バイオクリマティック変数 (BIO1-BIO19) を取得するとき
- 種分布モデル (SDM) の環境変数を準備するとき
- 気候変動シナリオの生息地適性を評価するとき
- 環境ニッチモデリングを実施するとき
- 土壌-植生-気候の相互作用を解析するとき

---

## Quick Start

## 1. SoilGrids 土壌特性取得

```python
import requests
import pandas as pd
import numpy as np

SOILGRIDS_BASE = "https://rest.isric.org/soilgrids/v2.0"


def soilgrids_get_properties(lat, lon, properties=None,
                              depths=None, values=None):
    """
    SoilGrids — 地点の土壌特性取得。

    Parameters:
        lat: float — 緯度
        lon: float — 経度
        properties: list[str] — 土壌特性 (例: ["phh2o", "soc", "clay"])
        depths: list[str] — 深度 (例: ["0-5cm", "5-15cm"])
        values: list[str] — 値の種類 (例: ["mean", "Q0.05", "Q0.95"])
    """
    if properties is None:
        properties = ["phh2o", "soc", "clay", "sand", "nitrogen",
                       "bdod", "cec", "ocd"]
    if depths is None:
        depths = ["0-5cm", "5-15cm", "15-30cm", "30-60cm"]
    if values is None:
        values = ["mean", "Q0.05", "Q0.95"]

    url = f"{SOILGRIDS_BASE}/properties/query"
    params = {
        "lat": lat,
        "lon": lon,
        "property": properties,
        "depth": depths,
        "value": values,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for layer in data.get("properties", {}).get("layers", []):
        prop_name = layer.get("name", "")
        unit = layer.get("unit_measure", {})
        conversion = unit.get("mapped_units", "")
        for depth_info in layer.get("depths", []):
            row = {
                "property": prop_name,
                "depth": depth_info.get("label", ""),
                "unit": conversion,
            }
            for val_key, val_val in depth_info.get("values", {}).items():
                row[val_key] = val_val
            results.append(row)

    df = pd.DataFrame(results)
    print(f"SoilGrids ({lat}, {lon}): {len(df)} records, "
          f"{len(properties)} properties")
    return df
```

## 2. WorldClim バイオクリマティック変数

```python
import rasterio
from rasterio.sample import sample_gen


def worldclim_get_bioclim(lat, lon, resolution="2.5m",
                           data_dir="worldclim"):
    """
    WorldClim — バイオクリマティック変数取得。

    Parameters:
        lat: float — 緯度
        lon: float — 経度
        resolution: str — 空間解像度 ("30s", "2.5m", "5m", "10m")
        data_dir: str — WorldClim データディレクトリ
    """
    from pathlib import Path
    bio_dir = Path(data_dir) / f"wc2.1_{resolution}_bio"

    bioclim_names = {
        1: "Annual Mean Temperature",
        2: "Mean Diurnal Range",
        3: "Isothermality",
        4: "Temperature Seasonality",
        5: "Max Temperature Warmest Month",
        6: "Min Temperature Coldest Month",
        7: "Temperature Annual Range",
        8: "Mean Temperature Wettest Quarter",
        9: "Mean Temperature Driest Quarter",
        10: "Mean Temperature Warmest Quarter",
        11: "Mean Temperature Coldest Quarter",
        12: "Annual Precipitation",
        13: "Precipitation Wettest Month",
        14: "Precipitation Driest Month",
        15: "Precipitation Seasonality",
        16: "Precipitation Wettest Quarter",
        17: "Precipitation Driest Quarter",
        18: "Precipitation Warmest Quarter",
        19: "Precipitation Coldest Quarter",
    }

    results = []
    for bio_num, bio_name in bioclim_names.items():
        tif_path = bio_dir / f"wc2.1_{resolution}_bio_{bio_num}.tif"
        if not tif_path.exists():
            continue
        with rasterio.open(tif_path) as src:
            vals = list(sample_gen(src, [(lon, lat)]))
            value = vals[0][0] if vals else None
        results.append({
            "variable": f"BIO{bio_num}",
            "name": bio_name,
            "value": value,
        })

    df = pd.DataFrame(results)
    print(f"WorldClim ({lat}, {lon}): {len(df)} bioclim variables")
    return df
```

## 3. 種分布モデル環境変数統合

```python
def sdm_environmental_stack(occurrences_df, lat_col="latitude",
                              lon_col="longitude", buffer_deg=0.5):
    """
    SDM — 種の出現記録に対する環境変数スタック生成。

    Parameters:
        occurrences_df: pd.DataFrame — 種出現記録
        lat_col: str — 緯度カラム名
        lon_col: str — 経度カラム名
        buffer_deg: float — バッファ距離 (度)
    """
    results = []
    for _, row in occurrences_df.iterrows():
        lat, lon = row[lat_col], row[lon_col]

        # SoilGrids
        soil = soilgrids_get_properties(lat, lon,
            properties=["phh2o", "soc", "clay"])
        soil_mean = {}
        for _, s in soil.iterrows():
            if s.get("depth") == "0-5cm":
                soil_mean[f"soil_{s['property']}"] = s.get("mean", None)

        # WorldClim (if available)
        bioclim = {}
        try:
            bio_df = worldclim_get_bioclim(lat, lon)
            bioclim = {r["variable"]: r["value"]
                        for _, r in bio_df.iterrows()}
        except Exception:
            pass

        combined = {
            lat_col: lat,
            lon_col: lon,
            **soil_mean,
            **bioclim,
        }
        results.append(combined)

    df = pd.DataFrame(results)
    print(f"SDM env stack: {len(df)} points, {len(df.columns)} variables")
    return df
```

## 4. 環境地理空間統合パイプライン

```python
def environmental_geodata_pipeline(occurrences_csv, output_dir="results"):
    """
    環境地理空間統合パイプライン。

    Parameters:
        occurrences_csv: str — 種出現記録 CSV パス
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    occ = pd.read_csv(occurrences_csv)
    print(f"Occurrences: {len(occ)} records")

    # 環境変数スタック
    env_df = sdm_environmental_stack(occ)
    env_df.to_csv(output_dir / "env_stack.csv", index=False)

    # 環境空間要約
    summary = env_df.describe().T
    summary.to_csv(output_dir / "env_summary.csv")

    print(f"Environmental pipeline: {output_dir}")
    return {"occurrences": occ, "env_stack": env_df, "summary": summary}
```

---

## ToolUniverse 連携

直接 REST API 使用 (SoilGrids, WorldClim は ToolUniverse 外)。

## パイプライン統合

```
environmental-ecology → environmental-geodata → marine-ecology
  (GBIF/iNaturalist)   (SoilGrids/WorldClim)   (OBIS/WoRMS)
       │                        │                     ↓
  phylogenetics ───────────────┘              biodiversity-indices
  (系統情報)             │                    (多様性指標)
                         ↓
                  species-distribution-model
                  (SDM 統合)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/env_stack.csv` | 環境変数スタック | → species-distribution-model |
| `results/env_summary.csv` | 環境空間要約 | → environmental-ecology |
