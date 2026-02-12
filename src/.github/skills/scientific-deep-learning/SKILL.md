---
name: scientific-deep-learning
description: |
  深層学習スキル。PyTorch Lightning・Hugging Face Transformers・timm を活用し、
  NN アーキテクチャ設計・転移学習・分散トレーニング・ハイパーパラメータ最適化・
  モデルデプロイを支援。
  「ニューラルネットで学習して」「Transformer を Fine-tune して」「深層学習モデルを構築して」で発火。
---

# Scientific Deep Learning

深層学習のための統合スキル。科学データに対する
NN アーキテクチャ設計・トレーニング・評価・デプロイメントを支援する。

## When to Use

- CNN / RNN / Transformer アーキテクチャ設計
- 転移学習・Fine-tuning（科学ドメイン特化モデル）
- 分散トレーニング（マルチ GPU / マルチノード）
- ハイパーパラメータ最適化 (Optuna / Ray Tune)
- 科学データ向けモデル（タンパク質言語モデル・化学 Transformer 等）
- モデルエクスポート・推論パイプライン

## Quick Start

### 深層学習パイプライン

```
Phase 1: Data Preparation
  - データセット構築 (Dataset / DataLoader)
  - データ拡張 (Augmentation)
  - 正規化・前処理
    ↓
Phase 2: Architecture Design
  - バックボーン選定 (ResNet/ViT/BERT/etc.)
  - ヘッド設計 (分類/回帰/生成)
  - 損失関数選択
    ↓
Phase 3: Training
  - Lightning Trainer 設定
  - オプティマイザ・スケジューラ
  - Mixed Precision (AMP)
  - 早期停止・チェックポイント
    ↓
Phase 4: Hyperparameter Optimization
  - Optuna / Ray Tune
  - 探索空間定義
  - Pruning (Median / Hyperband)
    ↓
Phase 5: Evaluation
  - テストセット評価
  - 学習曲線・混同行列
  - 解釈性分析 (→ XAI スキル)
    ↓
Phase 6: Deployment
  - ONNX / TorchScript エクスポート
  - バッチ推論パイプライン
  - モデルカード作成
```

## Workflow

### 1. PyTorch Lightning: 構造化パイプライン

```python
import torch
import torch.nn as nn
import torch.nn.functional as F
import pytorch_lightning as pl
from torch.utils.data import DataLoader, Dataset
from torchmetrics import Accuracy, AUROC, MeanSquaredError

# === Lightning Module ===
class ScientificModel(pl.LightningModule):
    def __init__(self, model, lr=1e-3, weight_decay=1e-5, task="classification"):
        super().__init__()
        self.save_hyperparameters(ignore=["model"])
        self.model = model
        self.lr = lr
        self.weight_decay = weight_decay
        self.task = task

        if task == "classification":
            self.criterion = nn.CrossEntropyLoss()
            self.train_acc = Accuracy(task="binary")
            self.val_acc = Accuracy(task="binary")
            self.val_auroc = AUROC(task="binary")
        else:
            self.criterion = nn.MSELoss()
            self.val_rmse = MeanSquaredError(squared=False)

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        x, y = batch
        logits = self(x)
        loss = self.criterion(logits, y)
        self.log("train/loss", loss, prog_bar=True)
        if self.task == "classification":
            preds = torch.argmax(logits, dim=1)
            self.train_acc(preds, y)
            self.log("train/acc", self.train_acc, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        x, y = batch
        logits = self(x)
        loss = self.criterion(logits, y)
        self.log("val/loss", loss, prog_bar=True)
        if self.task == "classification":
            preds = torch.argmax(logits, dim=1)
            self.val_acc(preds, y)
            self.val_auroc(logits[:, 1], y)
            self.log("val/acc", self.val_acc, prog_bar=True)
            self.log("val/auroc", self.val_auroc, prog_bar=True)
        else:
            self.val_rmse(logits.squeeze(), y)
            self.log("val/rmse", self.val_rmse, prog_bar=True)

    def configure_optimizers(self):
        optimizer = torch.optim.AdamW(self.parameters(), lr=self.lr,
                                       weight_decay=self.weight_decay)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer, T_max=self.trainer.max_epochs)
        return {"optimizer": optimizer, "lr_scheduler": scheduler}


# === トレーニング ===
from pytorch_lightning.callbacks import EarlyStopping, ModelCheckpoint, LearningRateMonitor

trainer = pl.Trainer(
    max_epochs=100,
    accelerator="auto",
    devices="auto",
    precision="16-mixed",  # AMP
    callbacks=[
        EarlyStopping(monitor="val/loss", patience=10, mode="min"),
        ModelCheckpoint(monitor="val/loss", save_top_k=3, mode="min"),
        LearningRateMonitor(logging_interval="epoch"),
    ],
    gradient_clip_val=1.0,
    deterministic=True,
)
# trainer.fit(model, train_loader, val_loader)
```

### 2. CNN アーキテクチャ (画像分類)

```python
import timm

# === timm: 事前学習済み CNN ===
def create_cnn_model(model_name="resnet50", num_classes=10, pretrained=True):
    """timm ライブラリから事前学習済みモデルをロード"""
    model = timm.create_model(
        model_name,
        pretrained=pretrained,
        num_classes=num_classes,
    )
    return model

# 利用可能モデル例:
# - "resnet50", "resnet101"
# - "efficientnet_b0" ~ "efficientnet_b7"
# - "vit_base_patch16_224" (Vision Transformer)
# - "swin_base_patch4_window7_224" (Swin Transformer)
# - "convnext_base"

# === カスタム CNN ===
class ConvBlock(nn.Module):
    def __init__(self, in_ch, out_ch, kernel_size=3, stride=1):
        super().__init__()
        self.conv = nn.Conv2d(in_ch, out_ch, kernel_size, stride, padding=kernel_size // 2)
        self.bn = nn.BatchNorm2d(out_ch)
        self.act = nn.GELU()

    def forward(self, x):
        return self.act(self.bn(self.conv(x)))


class ScientificCNN(nn.Module):
    def __init__(self, in_channels=1, num_classes=10):
        super().__init__()
        self.features = nn.Sequential(
            ConvBlock(in_channels, 32),
            ConvBlock(32, 64, stride=2),
            ConvBlock(64, 128, stride=2),
            ConvBlock(128, 256, stride=2),
            nn.AdaptiveAvgPool2d(1),
        )
        self.classifier = nn.Sequential(
            nn.Flatten(),
            nn.Dropout(0.3),
            nn.Linear(256, num_classes),
        )

    def forward(self, x):
        return self.classifier(self.features(x))
```

### 3. Transformer / Fine-tuning

```python
from transformers import AutoModel, AutoTokenizer, AutoModelForSequenceClassification
from transformers import Trainer, TrainingArguments

# === Hugging Face Transformer Fine-tuning ===
def finetune_transformer(model_name, train_dataset, val_dataset,
                          num_labels=2, num_epochs=5, lr=2e-5):
    """科学テキスト分類の Transformer Fine-tuning"""
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForSequenceClassification.from_pretrained(
        model_name, num_labels=num_labels
    )

    training_args = TrainingArguments(
        output_dir="./results/transformer_ft",
        num_train_epochs=num_epochs,
        per_device_train_batch_size=16,
        per_device_eval_batch_size=32,
        learning_rate=lr,
        weight_decay=0.01,
        warmup_ratio=0.1,
        evaluation_strategy="epoch",
        save_strategy="epoch",
        load_best_model_at_end=True,
        metric_for_best_model="eval_loss",
        fp16=True,
        dataloader_num_workers=4,
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=val_dataset,
        tokenizer=tokenizer,
    )

    trainer.train()
    return trainer

# 科学ドメイン特化モデル:
# - "allenai/scibert_scivocab_uncased" (SciBERT)
# - "microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract" (BiomedBERT)
# - "ChemBERTa" (化学 SMILES)
# - "facebook/esm2_t33_650M_UR50D" (タンパク質言語モデル)
```

### 4. Optuna ハイパーパラメータ最適化

```python
import optuna
from optuna.integration import PyTorchLightningPruningCallback

def optuna_optimization(train_loader, val_loader, n_trials=50):
    """Optuna によるハイパーパラメータ最適化"""

    def objective(trial):
        # 探索空間
        lr = trial.suggest_float("lr", 1e-5, 1e-2, log=True)
        weight_decay = trial.suggest_float("weight_decay", 1e-6, 1e-3, log=True)
        hidden_dim = trial.suggest_categorical("hidden_dim", [128, 256, 512])
        n_layers = trial.suggest_int("n_layers", 2, 6)
        dropout = trial.suggest_float("dropout", 0.1, 0.5)

        # モデル構築
        model = build_model(hidden_dim, n_layers, dropout)
        lit_model = ScientificModel(model, lr=lr, weight_decay=weight_decay)

        # トレーナー
        trainer = pl.Trainer(
            max_epochs=30,
            callbacks=[
                PyTorchLightningPruningCallback(trial, monitor="val/loss"),
                EarlyStopping(monitor="val/loss", patience=5),
            ],
            enable_progress_bar=False,
        )

        trainer.fit(lit_model, train_loader, val_loader)
        return trainer.callback_metrics["val/loss"].item()

    study = optuna.create_study(
        direction="minimize",
        pruner=optuna.pruners.HyperbandPruner(),
    )
    study.optimize(objective, n_trials=n_trials)

    print(f"Best trial: {study.best_trial.value:.4f}")
    print(f"Best params: {study.best_trial.params}")

    # 可視化
    optuna.visualization.plot_optimization_history(study).write_image(
        "figures/optuna_history.png")
    optuna.visualization.plot_param_importances(study).write_image(
        "figures/optuna_importance.png")

    return study
```

### 5. モデルエクスポート

```python
def export_model(model, sample_input, export_dir="models"):
    """モデルエクスポート (ONNX + TorchScript)"""
    import os
    os.makedirs(export_dir, exist_ok=True)

    model.eval()

    # TorchScript
    scripted = torch.jit.trace(model, sample_input)
    scripted.save(f"{export_dir}/model.pt")

    # ONNX
    torch.onnx.export(
        model, sample_input,
        f"{export_dir}/model.onnx",
        input_names=["input"],
        output_names=["output"],
        dynamic_axes={"input": {0: "batch_size"}, "output": {0: "batch_size"}},
        opset_version=17,
    )

    print(f"Exported: {export_dir}/model.pt, {export_dir}/model.onnx")
```

---

## Best Practices

1. **PyTorch Lightning で構造化**: 再現性と可読性のため Lightning Module を使用
2. **Mixed Precision**: `precision="16-mixed"` で GPU メモリ効率とスピードを改善
3. **Gradient Clipping**: `gradient_clip_val=1.0` で勾配爆発を防止
4. **Cosine Annealing**: 学習率スケジューラとして安定的
5. **Optuna + Pruning**: 不良トライアルを早期に枝刈りして効率化
6. **Scaffold Split** (分子データ): ランダム分割ではなく scaffold split を使用
7. **Model Card**: 学習データ・性能・制約を documentation
8. **Seed 固定**: `pl.seed_everything(42)` で再現性確保

## Completeness Checklist

- [ ] Dataset / DataLoader 構築
- [ ] アーキテクチャ設計・選定
- [ ] トレーニング収束確認
- [ ] 学習曲線プロット
- [ ] テストセット評価メトリクス
- [ ] ハイパーパラメータ最適化
- [ ] モデルエクスポート (ONNX)
- [ ] モデルカード作成

## References

### Output Files

| ファイル | 形式 | 生成タイミング |
|---|---|---|
| `results/dl_training_log.json` | トレーニングログ（JSON） | トレーニング完了時 |
| `models/model.onnx` | ONNX モデル | エクスポート時 |
| `figures/dl_learning_curve.png` | 学習曲線 | トレーニング完了時 |
| `figures/optuna_history.png` | Optuna 最適化履歴 | HPO 完了時 |

### 参照スキル

| スキル | 連携 |
|---|---|
| `scientific-ml-classification` | ← 従来 ML との比較ベースライン |
| `scientific-ml-regression` | ← 回帰タスクの DL 拡張 |
| `scientific-explainable-ai` | → DeepSHAP/Captum による説明 |
| `scientific-graph-neural-networks` | → グラフデータの DL |
| `scientific-image-analysis` | ← 画像データの CNN 応用 |
| `scientific-medical-imaging` | → 医用画像の DL モデル |
| `scientific-quantum-computing` | ← 量子-古典ハイブリッド ML |
