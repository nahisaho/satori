---
name: scientific-geospatial-analysis
description: |
  地理空間データ解析スキル。GeoPandas ベクターデータ処理・
  Rasterio ラスター解析・Folium/Kepler.gl インタラクティブ地図・
  空間自己相関 (Moran's I)・クリギング補間・CRS 変換。
tu_tools:
  - key: gbif
    name: GBIF
    description: 地理空間生物多様性データ取得
---

# Scientific Geospatial Analysis

地理空間データの前処理・空間統計・インタラクティブ地図可視化
パイプラインを提供する。

## When to Use

- GeoPandas でベクターデータ (Shapefile/GeoJSON) を処理するとき
- ラスターデータ (GeoTIFF) を読み込み解析するとき
- 空間自己相関 (Moran's I / LISA) を検定するとき
- クリギング (Kriging) で空間補間するとき
- Folium/Kepler.gl でインタラクティブ地図を作成するとき
- CRS (座標参照系) 変換・空間結合をするとき

> **Note**: 環境特化 GIS (SoilGrids/WorldClim) は `scientific-environmental-geodata` を参照。

---

## Quick Start

## 1. GeoPandas ベクターデータ処理

```python
import numpy as np
import pandas as pd


def load_and_process_geodata(filepath, target_crs="EPSG:4326"):
    """
    GeoPandas ベクター/ポイントデータ読み込み・CRS 変換。

    Parameters:
        filepath: str — Shapefile / GeoJSON / GPKG パス
        target_crs: str — 変換先座標系
    """
    import geopandas as gpd

    gdf = gpd.read_file(filepath)
    original_crs = gdf.crs

    if gdf.crs != target_crs:
        gdf = gdf.to_crs(target_crs)

    # 基本統計
    bounds = gdf.total_bounds  # [minx, miny, maxx, maxy]
    geom_types = gdf.geometry.geom_type.value_counts().to_dict()

    print(f"GeoData: {len(gdf)} features, CRS: {original_crs} → {target_crs}")
    print(f"  Bounds: [{bounds[0]:.4f}, {bounds[1]:.4f}] "
          f"to [{bounds[2]:.4f}, {bounds[3]:.4f}]")
    print(f"  Geometry types: {geom_types}")
    return gdf


def spatial_join(gdf_left, gdf_right, how="inner", predicate="intersects"):
    """
    空間結合 (Spatial Join)。

    Parameters:
        gdf_left: GeoDataFrame — 左テーブル
        gdf_right: GeoDataFrame — 右テーブル
        how: str — "inner" / "left" / "right"
        predicate: str — "intersects" / "within" / "contains"
    """
    import geopandas as gpd

    if gdf_left.crs != gdf_right.crs:
        gdf_right = gdf_right.to_crs(gdf_left.crs)

    joined = gpd.sjoin(gdf_left, gdf_right, how=how, predicate=predicate)

    print(f"Spatial Join ({predicate}, {how}): "
          f"{len(gdf_left)} × {len(gdf_right)} → {len(joined)}")
    return joined
```

## 2. 空間自己相関 (Moran's I / LISA)

```python
def spatial_autocorrelation(gdf, value_col, weight_type="queen"):
    """
    空間自己相関検定 — Global Moran's I + LISA。

    Parameters:
        gdf: GeoDataFrame — ジオメトリ + 属性データ
        value_col: str — 解析対象カラム
        weight_type: str — "queen" / "rook" / "knn"
    """
    from libpysal.weights import Queen, Rook, KNN
    from esda.moran import Moran, Moran_Local
    import matplotlib.pyplot as plt

    # 空間重み行列
    if weight_type == "queen":
        w = Queen.from_dataframe(gdf)
    elif weight_type == "rook":
        w = Rook.from_dataframe(gdf)
    elif weight_type == "knn":
        w = KNN.from_dataframe(gdf, k=5)

    w.transform = "r"
    y = gdf[value_col].values

    # Global Moran's I
    moran_global = Moran(y, w)

    # LISA (Local Indicators of Spatial Association)
    moran_local = Moran_Local(y, w)

    gdf = gdf.copy()
    gdf["lisa_cluster"] = moran_local.q  # 1=HH, 2=LH, 3=LL, 4=HL
    gdf["lisa_significant"] = moran_local.p_sim < 0.05
    gdf["local_moran_i"] = moran_local.Is

    # 可視化
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    gdf.plot(column=value_col, ax=ax1, legend=True,
             cmap="RdYlBu_r", edgecolor="gray", linewidth=0.3)
    ax1.set_title(f"{value_col} (Moran's I={moran_global.I:.4f}, "
                  f"p={moran_global.p_sim:.4f})")

    cluster_labels = {1: "High-High", 2: "Low-High",
                      3: "Low-Low", 4: "High-Low", 0: "Not Significant"}
    sig_gdf = gdf[gdf["lisa_significant"]]
    if len(sig_gdf) > 0:
        sig_gdf.plot(column="lisa_cluster", ax=ax2,
                     categorical=True, legend=True,
                     edgecolor="gray", linewidth=0.3)
    ax2.set_title("LISA Clusters (p < 0.05)")

    plt.tight_layout()
    path = "spatial_autocorrelation.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Moran's I = {moran_global.I:.4f}, p = {moran_global.p_sim:.4f}")
    print(f"LISA: {gdf['lisa_significant'].sum()} significant clusters")
    return {"moran_i": moran_global.I, "p_value": moran_global.p_sim,
            "gdf": gdf, "fig": path}
```

## 3. クリギング空間補間

```python
def kriging_interpolation(points_df, x_col, y_col, value_col,
                          grid_resolution=100,
                          variogram_model="spherical"):
    """
    Ordinary Kriging 空間補間。

    Parameters:
        points_df: pd.DataFrame — 観測点データ
        x_col, y_col: str — 座標カラム
        value_col: str — 補間対象カラム
        grid_resolution: int — グリッド解像度
        variogram_model: str — "spherical" / "exponential" / "gaussian"
    """
    from pykrige.ok import OrdinaryKriging
    import matplotlib.pyplot as plt

    x = points_df[x_col].values
    y = points_df[y_col].values
    z = points_df[value_col].values

    ok = OrdinaryKriging(
        x, y, z,
        variogram_model=variogram_model,
        verbose=False, enable_plotting=False)

    grid_x = np.linspace(x.min(), x.max(), grid_resolution)
    grid_y = np.linspace(y.min(), y.max(), grid_resolution)
    z_pred, ss_pred = ok.execute("grid", grid_x, grid_y)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    im1 = ax1.imshow(z_pred, origin="lower",
                     extent=[x.min(), x.max(), y.min(), y.max()],
                     cmap="viridis")
    ax1.scatter(x, y, c="red", s=10, edgecolors="black", linewidths=0.5)
    ax1.set_title(f"Kriging Prediction ({variogram_model})")
    plt.colorbar(im1, ax=ax1)

    im2 = ax2.imshow(ss_pred, origin="lower",
                     extent=[x.min(), x.max(), y.min(), y.max()],
                     cmap="Reds")
    ax2.set_title("Kriging Variance (Uncertainty)")
    plt.colorbar(im2, ax=ax2)

    plt.tight_layout()
    path = "kriging_result.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Kriging ({variogram_model}): {grid_resolution}×{grid_resolution} grid, "
          f"{len(x)} observation points")
    return {"z_pred": z_pred, "variance": ss_pred,
            "grid_x": grid_x, "grid_y": grid_y, "fig": path}
```

## 4. Folium インタラクティブ地図

```python
def interactive_map(gdf, value_col=None, popup_cols=None,
                    tiles="CartoDB positron",
                    output="interactive_map.html"):
    """
    Folium インタラクティブ地図。

    Parameters:
        gdf: GeoDataFrame — 地理空間データ
        value_col: str | None — Choropleth カラム
        popup_cols: list[str] | None — ポップアップ表示カラム
        tiles: str — タイル名
        output: str — 出力 HTML
    """
    import folium

    center = [gdf.geometry.centroid.y.mean(),
              gdf.geometry.centroid.x.mean()]
    m = folium.Map(location=center, zoom_start=8, tiles=tiles)

    if value_col and gdf.geometry.geom_type.iloc[0] in ["Polygon", "MultiPolygon"]:
        folium.Choropleth(
            geo_data=gdf.__geo_interface__,
            data=gdf, columns=[gdf.index.name or "index", value_col],
            key_on="feature.id",
            fill_color="YlOrRd",
            legend_name=value_col
        ).add_to(m)
    else:
        for _, row in gdf.iterrows():
            popup_text = ""
            if popup_cols:
                popup_text = "<br>".join(
                    [f"<b>{c}</b>: {row[c]}" for c in popup_cols])
            folium.CircleMarker(
                location=[row.geometry.centroid.y, row.geometry.centroid.x],
                radius=5, popup=popup_text,
                color="blue", fill=True
            ).add_to(m)

    m.save(output)
    print(f"Interactive map → {output} ({len(gdf)} features)")
    return output
```

---

## パイプライン統合

```
environmental-geodata → geospatial-analysis → advanced-visualization
  (環境データ取得)        (空間解析)             (高度可視化)
         │                     │                      ↓
  epidemiology ───────────────┘          interactive-dashboard
    (空間疫学)                              (ダッシュボード)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `spatial_autocorrelation.png` | Moran's I + LISA | → reporting |
| `kriging_result.png` | クリギング補間 | → visualization |
| `interactive_map.html` | Folium 地図 | → dashboard |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `gbif` | GBIF | 地理空間生物多様性データ取得 |
