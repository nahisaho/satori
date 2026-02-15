---
name: scientific-advanced-imaging
description: |
  高度バイオイメージング解析スキル。CellProfiler によるモフォロジカル
  プロファイリング・Cell Painting 解析、Cellpose による深層学習
  セルセグメンテーション、napari によるインタラクティブ 3D 可視化。
tu_tools:
  - key: biotools
    name: bio.tools
    description: 高度イメージングツール検索
---

# Scientific Advanced Imaging

CellProfiler / Cellpose / napari を活用した高度バイオイメージング
解析パイプラインを提供する。セルセグメンテーション、形態学的
特徴抽出、Cell Painting アッセイ解析、3D 画像可視化。

## When to Use

- Cell Painting アッセイ画像を解析するとき
- 深層学習によるセルセグメンテーションが必要なとき (Cellpose)
- 形態学的特徴量 (面積、真円度、テクスチャ等) を定量するとき
- 蛍光顕微鏡スタック画像を 3D 可視化するとき
- カスタム Cellpose モデルを微調整するとき
- 高コンテンツスクリーニング (HCS) データを処理するとき

---

## Quick Start

## 1. Cellpose セルセグメンテーション

```python
from cellpose import models, io
import numpy as np
from pathlib import Path


def cellpose_segmentation(image_path, model_type="cyto2", diameter=None,
                          channels=None, gpu=True):
    """
    Cellpose で細胞セグメンテーション。

    Parameters:
        image_path: str — 画像ファイルパス
        model_type: str — モデル ("cyto", "cyto2", "nuclei")
        diameter: float — 細胞直径 (None で自動)
        channels: list — チャネル [cytoplasm, nuclei]
        gpu: bool — GPU 使用

    K-Dense: cellpose
    """
    model = models.Cellpose(model_type=model_type, gpu=gpu)

    if channels is None:
        channels = [0, 0]  # grayscale

    img = io.imread(str(image_path))
    masks, flows, styles, diams = model.eval(
        img, diameter=diameter, channels=channels
    )

    n_cells = len(np.unique(masks)) - 1  # background除外
    print(f"Cellpose: {n_cells} cells detected (model={model_type})")
    print(f"  Auto diameter: {diams:.1f}")
    return masks, flows, styles
```

## 2. バッチセグメンテーション

```python
def batch_segmentation(image_dir, output_dir, model_type="cyto2",
                       image_pattern="*.tif", diameter=None):
    """
    ディレクトリ内画像のバッチセグメンテーション。

    Parameters:
        image_dir: str — 入力画像ディレクトリ
        output_dir: str — 出力ディレクトリ
        model_type: str — Cellpose モデル
    """
    model = models.Cellpose(model_type=model_type, gpu=True)
    image_dir = Path(image_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    image_files = sorted(image_dir.glob(image_pattern))
    results = []

    for img_path in image_files:
        img = io.imread(str(img_path))
        masks, flows, styles, diams = model.eval(
            img, diameter=diameter, channels=[0, 0]
        )
        n_cells = len(np.unique(masks)) - 1

        # マスク保存
        out_path = output_dir / f"{img_path.stem}_masks.npy"
        np.save(str(out_path), masks)

        results.append({
            "image": img_path.name,
            "n_cells": n_cells,
            "diameter": float(diams),
        })

    print(f"Batch: {len(results)} images, "
          f"{sum(r['n_cells'] for r in results)} total cells")
    return results
```

## 3. CellProfiler 形態学的プロファイリング

```python
import skimage.measure as measure
import pandas as pd


def morphological_profiling(image, masks, channel_names=None):
    """
    CellProfiler スタイルの形態学的特徴量抽出。

    Parameters:
        image: np.ndarray — 画像 (H, W) or (H, W, C)
        masks: np.ndarray — セグメンテーションマスク
        channel_names: list — チャネル名

    K-Dense: cellprofiler
    """
    if image.ndim == 2:
        image = image[..., np.newaxis]
    if channel_names is None:
        channel_names = [f"ch{i}" for i in range(image.shape[-1])]

    props = measure.regionprops(masks, intensity_image=image[..., 0])

    features = []
    for prop in props:
        feat = {
            "cell_id": prop.label,
            "area": prop.area,
            "perimeter": prop.perimeter,
            "eccentricity": prop.eccentricity,
            "solidity": prop.solidity,
            "major_axis": prop.major_axis_length,
            "minor_axis": prop.minor_axis_length,
        }
        # 真円度
        if prop.perimeter > 0:
            feat["circularity"] = (4 * np.pi * prop.area) / (prop.perimeter ** 2)
        else:
            feat["circularity"] = 0.0

        # チャネルごと蛍光強度
        for ch_idx, ch_name in enumerate(channel_names):
            ch_props = measure.regionprops(
                masks, intensity_image=image[..., ch_idx]
            )
            for cp in ch_props:
                if cp.label == prop.label:
                    feat[f"{ch_name}_mean"] = cp.mean_intensity
                    feat[f"{ch_name}_max"] = cp.max_intensity
                    feat[f"{ch_name}_min"] = cp.min_intensity
                    break

        features.append(feat)

    df = pd.DataFrame(features)
    print(f"Morphological profiling: {len(df)} cells, {len(df.columns)} features")
    return df
```

## 4. Cell Painting 解析

```python
def cell_painting_analysis(image_5ch, masks, normalize=True):
    """
    Cell Painting 5 チャネル解析。

    Channels:
      ch0: DNA (Hoechst), ch1: ER (ConA), ch2: RNA (SYTO14),
      ch3: AGP (Phalloidin), ch4: Mito (MitoTracker)

    Parameters:
        image_5ch: np.ndarray — (H, W, 5) 画像
        masks: np.ndarray — セグメンテーションマスク
        normalize: bool — Z-score 正規化
    """
    channel_names = ["DNA", "ER", "RNA", "AGP", "Mito"]

    # 形態 + 蛍光特徴量抽出
    features_df = morphological_profiling(image_5ch, masks, channel_names)

    # テクスチャ特徴量 (Haralick-like)
    from skimage.feature import graycomatrix, graycoprops

    texture_features = []
    for prop in measure.regionprops(masks):
        cell_mask = masks == prop.label
        bbox = prop.bbox
        cell_region = image_5ch[bbox[0]:bbox[2], bbox[1]:bbox[3], 0]
        cell_mask_crop = cell_mask[bbox[0]:bbox[2], bbox[1]:bbox[3]]
        cell_region = (cell_region * cell_mask_crop).astype(np.uint8)

        if cell_region.max() > 0:
            glcm = graycomatrix(cell_region, [1], [0], levels=256, symmetric=True)
            texture_features.append({
                "cell_id": prop.label,
                "contrast": graycoprops(glcm, "contrast")[0, 0],
                "homogeneity": graycoprops(glcm, "homogeneity")[0, 0],
                "energy": graycoprops(glcm, "energy")[0, 0],
                "correlation": graycoprops(glcm, "correlation")[0, 0],
            })

    texture_df = pd.DataFrame(texture_features)
    combined = features_df.merge(texture_df, on="cell_id", how="left")

    if normalize:
        numeric_cols = combined.select_dtypes(include=[np.number]).columns
        numeric_cols = [c for c in numeric_cols if c != "cell_id"]
        combined[numeric_cols] = (
            (combined[numeric_cols] - combined[numeric_cols].mean())
            / combined[numeric_cols].std()
        )

    print(f"Cell Painting: {len(combined)} cells, {len(combined.columns)} features")
    return combined
```

## 5. napari 3D 可視化

```python
def napari_visualization(image, masks=None, points=None, labels=None):
    """
    napari でインタラクティブ 3D 可視化。

    Parameters:
        image: np.ndarray — 画像 (Z, Y, X) or (Z, Y, X, C)
        masks: np.ndarray — セグメンテーションマスク
        points: np.ndarray — 座標点 (N, 3)
        labels: list — ポイントラベル

    K-Dense: napari
    """
    import napari

    viewer = napari.Viewer()

    # 画像レイヤー
    if image.ndim == 4:
        for ch in range(image.shape[-1]):
            viewer.add_image(
                image[..., ch],
                name=f"Channel-{ch}",
                blending="additive",
            )
    else:
        viewer.add_image(image, name="Image")

    # マスクレイヤー
    if masks is not None:
        viewer.add_labels(masks, name="Segmentation")

    # ポイントレイヤー
    if points is not None:
        properties = {"label": labels} if labels else None
        viewer.add_points(
            points,
            properties=properties,
            name="Points",
            size=5,
        )

    print(f"napari viewer: {image.shape}, "
          f"masks={'Yes' if masks is not None else 'No'}")
    return viewer
```

## 6. Cellpose カスタムモデル微調整

```python
def finetune_cellpose(train_dir, model_type="cyto2", n_epochs=100,
                      learning_rate=0.1):
    """
    Cellpose モデルの微調整。

    Parameters:
        train_dir: str — 訓練データ (画像 + _masks)
        model_type: str — ベースモデル
        n_epochs: int — エポック数
        learning_rate: float — 学習率
    """
    from cellpose import train

    train_dir = Path(train_dir)
    image_files = sorted(train_dir.glob("*.tif"))
    mask_files = sorted(train_dir.glob("*_masks.tif"))

    images = [io.imread(str(f)) for f in image_files]
    labels = [io.imread(str(f)) for f in mask_files]

    model = models.CellposeModel(model_type=model_type, gpu=True)
    model_path = model.train(
        images, labels,
        channels=[0, 0],
        n_epochs=n_epochs,
        learning_rate=learning_rate,
        save_path=str(train_dir / "models"),
    )

    print(f"Fine-tuned model saved: {model_path}")
    print(f"  Base: {model_type}, Epochs: {n_epochs}")
    return model_path
```

## 7. 統合イメージングパイプライン

```python
def imaging_pipeline(image_dir, output_dir, model_type="cyto2",
                     channels_5=True):
    """
    セグメンテーション + 特徴量抽出 統合パイプライン。

    Parameters:
        image_dir: str — 入力画像ディレクトリ
        output_dir: str — 出力ディレクトリ
        model_type: str — Cellpose モデル
        channels_5: bool — 5ch Cell Painting
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) バッチセグメンテーション
    seg_results = batch_segmentation(
        image_dir, output_dir / "masks", model_type=model_type
    )

    # 2) 特徴量抽出
    all_features = []
    for res in seg_results:
        img_path = Path(image_dir) / res["image"]
        mask_path = output_dir / "masks" / f"{img_path.stem}_masks.npy"

        img = io.imread(str(img_path))
        masks = np.load(str(mask_path))

        if channels_5 and img.ndim == 3 and img.shape[-1] == 5:
            feats = cell_painting_analysis(img, masks)
        else:
            feats = morphological_profiling(img, masks)

        feats["image"] = res["image"]
        all_features.append(feats)

    combined = pd.concat(all_features, ignore_index=True)
    combined.to_csv(output_dir / "features.csv", index=False)

    print(f"Pipeline complete: {len(seg_results)} images, "
          f"{len(combined)} cells profiled")
    return combined
```

---

## パイプライン統合

```
image-analysis → advanced-imaging → medical-imaging
  (OpenCV/基本)   (Cellpose/CellProfiler)   (DICOM/MONAI)
       │                   │                      ↓
fluorescence-microscopy ──┘               drug-target-profiling
  (蛍光画像取得)          │                 (Cell Painting hit)
                          ↓
                   cheminformatics
                   (Cell Painting → 化合物活性)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/masks/` | セグメンテーションマスク群 | → medical-imaging |
| `results/features.csv` | 形態学的特徴量マトリクス | → cheminformatics |
| `results/cell_painting.csv` | Cell Painting プロファイル | → drug-target-profiling |
| `results/model/` | 微調整 Cellpose モデル | — |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `biotools` | bio.tools | 高度イメージングツール検索 |
