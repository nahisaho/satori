---
name: scientific-medical-imaging
description: |
  医用イメージングスキル。DICOM/NIfTI 処理・WSI (Whole Slide Image) 解析・
  PathML・MONAI・3D Slicer 連携・放射線画像解析・病理組織画像解析を支援。
  「DICOM を解析して」「WSI を処理して」「医用画像をセグメンテーションして」で発火。
---

# Scientific Medical Imaging

医用イメージングのための専門解析スキル。
DICOM / NIfTI の放射線画像と WSI (Whole Slide Image) の
病理組織画像を対象とした解析パイプラインを提供する。

## When to Use

- DICOM / NIfTI 医用画像の読み込み・前処理
- CT / MRI / PET 画像のセグメンテーション
- WSI (Whole Slide Image) 病理組織画像解析
- 放射線画像の特徴量抽出 (Radiomics)
- 医用画像の深層学習モデル (MONAI)
- 3D ボリュームレンダリング・可視化

## Quick Start

### 医用イメージングパイプライン

```
Phase 1: Data Ingestion
  - DICOM / NIfTI ファイル読み込み
  - メタデータ抽出 (患者情報・撮影条件)
  - 匿名化処理
    ↓
Phase 2: Preprocessing
  - リサンプリング (Isotropic spacing)
  - 正規化 (Windowing / Z-score)
  - Registration (位置合わせ)
    ↓
Phase 3: Segmentation
  - U-Net / nnU-Net / Swin UNETR
  - 臓器/病変セグメンテーション
  - 後処理 (形態学的操作)
    ↓
Phase 4: Feature Extraction
  - Radiomics (PyRadiomics)
  - 形状・テクスチャ・強度特徴量
  - WSI パッチレベル特徴量
    ↓
Phase 5: Analysis
  - 腫瘍体積・成長率計算
  - 放射線治療計画支援
  - 病理組織分類
    ↓
Phase 6: Reporting
  - DICOM-SR (Structured Report)
  - 解析レポート生成
  - 3D 可視化
```

## Workflow

### 1. DICOM 画像処理

```python
import pydicom
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# === DICOM 読み込み ===
def load_dicom_series(dicom_dir):
    """DICOM シリーズ → 3D ボリューム"""
    dicom_dir = Path(dicom_dir)
    slices = []

    for dcm_file in sorted(dicom_dir.glob("*.dcm")):
        ds = pydicom.dcmread(dcm_file)
        slices.append(ds)

    # スライス位置でソート
    slices.sort(key=lambda s: float(s.ImagePositionPatient[2]))

    # メタデータ
    metadata = {
        "patient_id": slices[0].PatientID,
        "modality": slices[0].Modality,
        "rows": slices[0].Rows,
        "cols": slices[0].Columns,
        "n_slices": len(slices),
        "pixel_spacing": list(slices[0].PixelSpacing),
        "slice_thickness": float(slices[0].SliceThickness),
    }

    # 3D volume 構築
    volume = np.stack([s.pixel_array for s in slices])

    # HU (Hounsfield Unit) 変換
    slope = float(slices[0].RescaleSlope)
    intercept = float(slices[0].RescaleIntercept)
    volume_hu = volume * slope + intercept

    print(f"Volume shape: {volume_hu.shape}")
    print(f"HU range: [{volume_hu.min()}, {volume_hu.max()}]")

    return volume_hu, metadata


# === Windowing (ウィンドウ処理) ===
def apply_window(volume, window_center, window_width):
    """CT ウィンドウ設定"""
    # 標準ウィンドウ設定:
    # Lung: WC=-600, WW=1500
    # Bone: WC=300, WW=1500
    # Brain: WC=40, WW=80
    # Liver: WC=50, WW=350
    lower = window_center - window_width / 2
    upper = window_center + window_width / 2
    volume_windowed = np.clip(volume, lower, upper)
    volume_windowed = (volume_windowed - lower) / (upper - lower)
    return volume_windowed
```

### 2. NIfTI 処理 (SimpleITK)

```python
import SimpleITK as sitk

def load_nifti(filepath):
    """NIfTI 読み込み"""
    img = sitk.ReadImage(filepath)
    array = sitk.GetArrayFromImage(img)
    spacing = img.GetSpacing()
    origin = img.GetOrigin()
    direction = img.GetDirection()

    return {
        "array": array,
        "spacing": spacing,
        "origin": origin,
        "direction": direction,
        "shape": array.shape,
    }


def resample_isotropic(image, new_spacing=(1.0, 1.0, 1.0)):
    """Isotropic リサンプリング"""
    original_spacing = image.GetSpacing()
    original_size = image.GetSize()

    new_size = [
        int(round(osz * osp / nsp))
        for osz, osp, nsp in zip(original_size, original_spacing, new_spacing)
    ]

    resampler = sitk.ResampleImageFilter()
    resampler.SetOutputSpacing(new_spacing)
    resampler.SetSize(new_size)
    resampler.SetOutputDirection(image.GetDirection())
    resampler.SetOutputOrigin(image.GetOrigin())
    resampler.SetTransform(sitk.Transform())
    resampler.SetDefaultPixelValue(image.GetPixelIDValue())
    resampler.SetInterpolator(sitk.sitkBSpline)

    return resampler.Execute(image)


def registration(fixed_image, moving_image):
    """画像レジストレーション (位置合わせ)"""
    registration_method = sitk.ImageRegistrationMethod()

    # Similarity metric
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.01)

    # Optimizer
    registration_method.SetOptimizerAsGradientDescent(
        learningRate=1.0, numberOfIterations=100,
        convergenceMinimumValue=1e-6, convergenceWindowSize=10)

    # Transform
    initial_transform = sitk.CenteredTransformInitializer(
        fixed_image, moving_image, sitk.Euler3DTransform(),
        sitk.CenteredTransformInitializerFilter.GEOMETRY)
    registration_method.SetInitialTransform(initial_transform, inPlace=False)

    final_transform = registration_method.Execute(fixed_image, moving_image)
    return final_transform
```

### 3. MONAI: 医用画像深層学習

```python
import monai
from monai.networks.nets import UNet, SwinUNETR
from monai.losses import DiceLoss, DiceCELoss
from monai.metrics import DiceMetric
from monai.transforms import (
    Compose, LoadImaged, EnsureChannelFirstd, Spacingd,
    ScaleIntensityRanged, CropForegroundd, RandCropByPosNegLabeld,
    RandFlipd, RandRotate90d,
)

# === MONAI Transform パイプライン ===
train_transforms = Compose([
    LoadImaged(keys=["image", "label"]),
    EnsureChannelFirstd(keys=["image", "label"]),
    Spacingd(keys=["image", "label"], pixdim=(1.5, 1.5, 2.0), mode=("bilinear", "nearest")),
    ScaleIntensityRanged(keys=["image"], a_min=-175, a_max=250, b_min=0.0, b_max=1.0, clip=True),
    CropForegroundd(keys=["image", "label"], source_key="image"),
    RandCropByPosNegLabeld(
        keys=["image", "label"], label_key="label",
        spatial_size=(96, 96, 96), pos=1, neg=1, num_samples=4,
    ),
    RandFlipd(keys=["image", "label"], prob=0.5, spatial_axis=0),
    RandRotate90d(keys=["image", "label"], prob=0.5),
])

# === U-Net for 3D Segmentation ===
model = UNet(
    spatial_dims=3,
    in_channels=1,
    out_channels=14,  # 臓器数
    channels=(16, 32, 64, 128, 256),
    strides=(2, 2, 2, 2),
    num_res_units=2,
    dropout=0.2,
)

# === Swin UNETR (Transformer ベース) ===
model_swin = SwinUNETR(
    img_size=(96, 96, 96),
    in_channels=1,
    out_channels=14,
    feature_size=48,
    use_checkpoint=True,
)

# Loss + Optimizer
loss_fn = DiceCELoss(to_onehot_y=True, softmax=True)
dice_metric = DiceMetric(include_background=False, reduction="mean")
```

### 4. WSI (Whole Slide Image) 解析

```python
# === openslide で WSI 読み込み ===
import openslide

def load_wsi(wsi_path):
    """WSI (svs/ndpi/tiff) 読み込み"""
    slide = openslide.OpenSlide(wsi_path)

    metadata = {
        "dimensions": slide.dimensions,  # (width, height) at level 0
        "level_count": slide.level_count,
        "level_dimensions": slide.level_dimensions,
        "level_downsamples": slide.level_downsamples,
        "properties": dict(slide.properties),
    }

    print(f"WSI size: {slide.dimensions[0]:,} x {slide.dimensions[1]:,} pixels")
    print(f"Levels: {slide.level_count}")

    return slide, metadata


def extract_tissue_patches(slide, patch_size=256, level=0, threshold=0.7):
    """
    組織領域からパッチを抽出。
    Otsu 閾値で背景を除去し、組織含有率が threshold 以上のパッチのみ取得。
    """
    from skimage.filters import threshold_otsu
    from skimage.color import rgb2gray

    # サムネイル取得 (組織検出用)
    thumb_level = min(slide.level_count - 1, 4)
    thumbnail = slide.get_thumbnail(
        (slide.level_dimensions[thumb_level][0], slide.level_dimensions[thumb_level][1])
    )
    thumb_gray = rgb2gray(np.array(thumbnail))
    thresh = threshold_otsu(thumb_gray)
    tissue_mask = thumb_gray < thresh

    # パッチ座標計算
    downsample = slide.level_downsamples[thumb_level]
    patches = []
    w, h = slide.level_dimensions[level]

    for y in range(0, h, patch_size):
        for x in range(0, w, patch_size):
            # サムネイル上の対応位置
            tx = int(x / downsample)
            ty = int(y / downsample)
            tw = int(patch_size / downsample)
            th = int(patch_size / downsample)

            if tx + tw <= tissue_mask.shape[1] and ty + th <= tissue_mask.shape[0]:
                tissue_ratio = tissue_mask[ty:ty+th, tx:tx+tw].mean()
                if tissue_ratio >= threshold:
                    patch = slide.read_region((x, y), level, (patch_size, patch_size))
                    patches.append({
                        "image": np.array(patch.convert("RGB")),
                        "x": x, "y": y,
                        "tissue_ratio": round(tissue_ratio, 3),
                    })

    print(f"Extracted {len(patches)} tissue patches from WSI")
    return patches
```

### 5. Radiomics 特徴量抽出

```python
from radiomics import featureextractor

def extract_radiomics(image_path, mask_path, params=None):
    """
    PyRadiomics で放射線画像特徴量を抽出。
    Shape, First-Order, GLCM, GLRLM, GLSZM, GLDM, NGTDM
    """
    if params is None:
        params = {
            "setting": {
                "binWidth": 25,
                "resampledPixelSpacing": [1, 1, 1],
                "interpolator": "sitkBSpline",
                "normalize": True,
            },
            "featureClass": {
                "shape": None,
                "firstorder": None,
                "glcm": None,
                "glrlm": None,
                "glszm": None,
            },
        }

    extractor = featureextractor.RadiomicsFeatureExtractor()
    for key, val in params.get("setting", {}).items():
        extractor.settings[key] = val

    features = extractor.execute(image_path, mask_path)

    # diagnostics を除外
    radiomic_features = {
        k: float(v) for k, v in features.items()
        if not k.startswith("diagnostics_")
    }

    print(f"Extracted {len(radiomic_features)} radiomic features")
    return radiomic_features
```

### 6. レポート生成

```python
import json

def generate_imaging_report(patient_id, modality, findings,
                            segmentation_results=None, radiomics=None,
                            output_dir="results"):
    """医用画像解析レポート"""
    report = {
        "patient_id": patient_id,
        "modality": modality,
        "analysis_date": pd.Timestamp.now().isoformat() if "pd" in dir() else "",
        "findings": findings,
    }

    if segmentation_results:
        report["segmentation"] = {
            "model": segmentation_results.get("model", ""),
            "dice_scores": segmentation_results.get("dice_scores", {}),
            "volumes_ml": segmentation_results.get("volumes", {}),
        }

    if radiomics:
        report["radiomics"] = {
            "n_features": len(radiomics),
            "top_features": dict(sorted(
                radiomics.items(), key=lambda x: abs(x[1]), reverse=True
            )[:20]),
        }

    with open(f"{output_dir}/imaging_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    md = f"# Medical Imaging Report\n\n"
    md += f"**Patient**: {patient_id} | **Modality**: {modality}\n\n"
    md += f"## Findings\n\n"
    for finding in findings:
        md += f"- {finding}\n"

    if segmentation_results:
        md += f"\n## Segmentation Results\n\n"
        md += "| Structure | Dice Score | Volume (ml) |\n|---|---|---|\n"
        for struct, dice in segmentation_results.get("dice_scores", {}).items():
            vol = segmentation_results.get("volumes", {}).get(struct, "")
            md += f"| {struct} | {dice:.3f} | {vol} |\n"

    with open(f"{output_dir}/imaging_report.md", "w") as f:
        f.write(md)

    return report
```

---

## Best Practices

1. **Isotropic リサンプリング**: 異方性ボクセルは必ず等方性に統一してから解析
2. **Windowing**: CT は組織に適したウィンドウ設定を適用 (肺・骨・軟部組織)
3. **nnU-Net**: セグメンテーションタスクでは nnU-Net が自動設定で SOTA レベル
4. **WSI 解析はパッチベース**: ギガピクセル画像は直接処理できないためタイル化
5. **Dice + Hausdorff**: セグメンテーション評価は Dice と Hausdorff Distance を併用
6. **匿名化**: DICOM メタデータの患者情報を必ず匿名化
7. **3D context**: 2D スライスではなく 3D ボリュームベースの解析を推奨

## Completeness Checklist

- [ ] DICOM/NIfTI 読み込み・メタデータ抽出
- [ ] 前処理（リサンプリング・正規化・ウィンドウ処理）
- [ ] セグメンテーション実行・Dice スコア評価
- [ ] Radiomics 特徴量抽出（該当時）
- [ ] WSI パッチ抽出・組織検出（病理画像の場合）
- [ ] 3D 可視化
- [ ] レポート（JSON + Markdown）生成

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/imaging_report.json` | 画像解析レポート（JSON） | 解析完了時 |
| `results/imaging_report.md` | 画像解析レポート（Markdown） | レポート生成時 |
| `results/radiomics_features.json` | Radiomics 特徴量（JSON） | 特徴量抽出時 |
| `figures/segmentation_overlay.png` | セグメンテーション結果図 | セグメンテーション時 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-image-analysis` | ← 一般的な画像処理・セグメンテーション基盤 |
| `scientific-deep-learning` | ← U-Net/SwinUNETR のトレーニング基盤 |
| `scientific-ml-classification` | ← Radiomics 特徴量を用いた分類 |
| `scientific-precision-oncology` | → 腫瘍画像の精密医療連携 |
| `scientific-clinical-decision-support` | → 画像所見の臨床意思決定反映 |
| `scientific-explainable-ai` | → 医用 AI の説明可能性 (GradCAM) |
