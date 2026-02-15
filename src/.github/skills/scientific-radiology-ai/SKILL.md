---
name: scientific-radiology-ai
description: |
  放射線診断支援 AI スキル。CADe/CADx パイプライン・
  CT/MRI 分類・セグメンテーション・Grad-CAM 説明可能性・
  構造化レポート・AI-RADS グレーディング。
  ※ scientific-medical-imaging (DICOM/WSI/Radiomics) の
  放射線診断 AI 特化拡張。
tu_tools:
  - key: tcia
    name: TCIA
    description: 放射線画像データセット検索
---

# Scientific Radiology AI

放射線画像（CT/MRI/X 線）に対する AI 診断支援
パイプラインを提供する。MONAI ベースの学習・推論・
説明可能性・構造化レポート生成を含む。

## When to Use

- CT/MRI/X 線画像の AI 分類・セグメンテーションを行うとき
- CADe (検出) / CADx (診断) パイプラインを構築するとき
- Grad-CAM で AI 判断の説明可能性を付与するとき
- 構造化放射線レポートを自動生成するとき
- AI-RADS スコアリングを実装するとき

---

## Quick Start

## 1. MONAI 放射線 AI 分類パイプライン

```python
import numpy as np
import torch
import torch.nn as nn


def build_radiology_classifier(in_channels=1, num_classes=2,
                                spatial_dims=3,
                                architecture="densenet121"):
    """
    MONAI ベース放射線画像分類モデル。

    Parameters:
        in_channels: int — 入力チャネル数 (CT=1, MRI multimodal=4)
        num_classes: int — クラス数
        spatial_dims: int — 2 (2D スライス) or 3 (3D ボリューム)
        architecture: str — "densenet121" / "resnet50" / "efficientnet"
    """
    import monai.networks.nets as nets

    models = {
        "densenet121": nets.DenseNet121(
            spatial_dims=spatial_dims,
            in_channels=in_channels,
            out_channels=num_classes),
        "resnet50": nets.ResNet(
            block="bottleneck", layers=[3, 4, 6, 3],
            block_inplanes=[64, 128, 256, 512],
            spatial_dims=spatial_dims,
            n_input_channels=in_channels,
            num_classes=num_classes),
        "efficientnet": nets.EfficientNetBN(
            "efficientnet-b0",
            spatial_dims=spatial_dims,
            in_channels=in_channels,
            num_classes=num_classes),
    }
    model = models.get(architecture, models["densenet121"])
    total_params = sum(p.numel() for p in model.parameters())
    print(f"Radiology classifier: {architecture} | "
          f"{total_params:,} params | {spatial_dims}D")
    return model


def train_radiology_model(model, train_loader, val_loader,
                          epochs=50, lr=1e-4, device="cuda"):
    """
    放射線 AI モデル学習。

    Parameters:
        model: nn.Module — 分類モデル
        train_loader: DataLoader — 訓練データ
        val_loader: DataLoader — 検証データ
        epochs: int — 学習エポック数
        lr: float — 学習率
        device: str — デバイス
    """
    import pandas as pd
    from monai.utils import set_determinism
    set_determinism(seed=42)

    model.to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr,
                                  weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=epochs)
    criterion = nn.CrossEntropyLoss()
    history = []

    best_val_acc = 0
    for epoch in range(epochs):
        model.train()
        train_loss, correct, total = 0, 0, 0
        for batch in train_loader:
            images = batch["image"].to(device)
            labels = batch["label"].to(device)
            optimizer.zero_grad()
            outputs = model(images)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
            correct += (outputs.argmax(1) == labels).sum().item()
            total += len(labels)

        scheduler.step()

        # Validation
        model.eval()
        val_loss, val_correct, val_total = 0, 0, 0
        with torch.no_grad():
            for batch in val_loader:
                images = batch["image"].to(device)
                labels = batch["label"].to(device)
                outputs = model(images)
                val_loss += criterion(outputs, labels).item()
                val_correct += (outputs.argmax(1) == labels).sum().item()
                val_total += len(labels)

        val_acc = val_correct / val_total
        if val_acc > best_val_acc:
            best_val_acc = val_acc
            torch.save(model.state_dict(), "best_radiology_model.pt")

        history.append({
            "epoch": epoch + 1,
            "train_loss": train_loss / len(train_loader),
            "train_acc": correct / total,
            "val_loss": val_loss / len(val_loader),
            "val_acc": val_acc,
        })

        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}: train_acc={correct/total:.3f}, "
                  f"val_acc={val_acc:.3f}")

    print(f"Best val_acc: {best_val_acc:.4f}")
    return pd.DataFrame(history)
```

## 2. Grad-CAM 説明可能性

```python
def radiology_gradcam(model, image_tensor, target_layer=None,
                      target_class=None, device="cuda"):
    """
    放射線画像に対する Grad-CAM 可視化。

    Parameters:
        model: nn.Module — 学習済み分類モデル
        image_tensor: torch.Tensor — 入力画像 [1, C, H, W] or [1, C, D, H, W]
        target_layer: nn.Module | None — CAM 対象層
        target_class: int | None — 対象クラス (None=予測クラス)
        device: str — デバイス
    """
    import matplotlib.pyplot as plt
    from monai.visualize import GradCAM

    model.to(device).eval()
    image_tensor = image_tensor.to(device)

    if target_layer is None:
        # DenseNet の最終 features 層を使用
        for name, module in model.named_modules():
            if "features" in name or "layer4" in name:
                target_layer = name
        if target_layer is None:
            target_layer = list(model.named_modules())[-2][0]

    cam = GradCAM(nn_module=model, target_layers=target_layer)

    if target_class is None:
        with torch.no_grad():
            target_class = model(image_tensor).argmax(1).item()

    result = cam(x=image_tensor, class_idx=target_class)
    cam_map = result.squeeze().cpu().numpy()

    # 2D スライス可視化
    if cam_map.ndim == 3:
        mid_slice = cam_map.shape[0] // 2
        cam_map_2d = cam_map[mid_slice]
        img_2d = image_tensor.squeeze().cpu().numpy()[mid_slice]
    else:
        cam_map_2d = cam_map
        img_2d = image_tensor.squeeze().cpu().numpy()

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].imshow(img_2d, cmap="gray")
    axes[0].set_title("Original")
    axes[1].imshow(cam_map_2d, cmap="jet")
    axes[1].set_title(f"Grad-CAM (class={target_class})")
    axes[2].imshow(img_2d, cmap="gray")
    axes[2].imshow(cam_map_2d, cmap="jet", alpha=0.4)
    axes[2].set_title("Overlay")
    for ax in axes:
        ax.axis("off")
    plt.tight_layout()
    plt.savefig("gradcam_radiology.png", dpi=150, bbox_inches="tight")
    print(f"Grad-CAM saved → gradcam_radiology.png (class={target_class})")
    return cam_map
```

## 3. 構造化放射線レポート

```python
def generate_structured_report(predictions, patient_info=None,
                               modality="CT", body_part="Chest"):
    """
    AI 支援構造化放射線レポート生成。

    Parameters:
        predictions: dict — {"finding": str, "probability": float, ...}
        patient_info: dict | None — 患者情報
        modality: str — "CT" / "MRI" / "XR"
        body_part: str — 検査部位
    """
    if patient_info is None:
        patient_info = {"id": "ANON", "age": "N/A", "sex": "N/A"}

    findings = []
    for finding, prob in predictions.items():
        if prob >= 0.5:
            confidence = "High" if prob >= 0.8 else "Moderate"
            findings.append(f"- {finding}: {prob:.1%} ({confidence} confidence)")

    report = f"""## Structured Radiology Report (AI-Assisted)

**Patient**: {patient_info.get('id', 'N/A')} | \
Age: {patient_info.get('age', 'N/A')} | Sex: {patient_info.get('sex', 'N/A')}
**Modality**: {modality} | **Body Part**: {body_part}

### AI Findings

{chr(10).join(findings) if findings else '- No significant findings detected'}

### AI Confidence Summary

| Finding | Probability | AI-RADS |
|---------|:-----------:|:-------:|
"""
    for finding, prob in sorted(predictions.items(),
                                key=lambda x: x[1], reverse=True):
        rads = 5 if prob >= 0.9 else 4 if prob >= 0.7 else \
            3 if prob >= 0.5 else 2 if prob >= 0.3 else 1
        report += f"| {finding} | {prob:.1%} | {rads} |\n"

    report += """
### Disclaimer
> This report was generated with AI assistance and requires
> review by a qualified radiologist before clinical use.
"""
    print(report)
    return report
```

---

## パイプライン統合

```
[DICOM 取得] → medical-imaging → radiology-ai → clinical-report
               (前処理/Radiomics)  (AI 診断)     (臨床レポート)
                                       │
                              explainable-ai ← deep-learning
                               (説明可能性)     (基盤学習)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `best_radiology_model.pt` | 学習済み分類モデル | → 推論 |
| `gradcam_radiology.png` | Grad-CAM 可視化 | → レポート |
| `structured_report.md` | 構造化レポート | → clinical-report |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `tcia` | TCIA | 放射線画像データセット検索 |
