---
name: scientific-transfer-learning
description: |
  転移学習・ドメイン適応スキル。事前学習モデルファインチューニング・
  Few-shot / Zero-shot 学習・ドメイン適応 (DA)・
  知識蒸留・マルチタスク学習・科学ドメイン特化モデル転移。
---

# Scientific Transfer Learning

事前学習モデルの科学データへの転移・ドメイン適応・
Few-shot 学習パイプラインを提供する。

## When to Use

- 事前学習済みモデル (ImageNet/BERT) をファインチューニングするとき
- 小規模科学データセットで高精度を実現したいとき
- ドメイン適応で異なるデータ分布間のギャップを埋めるとき
- Few-shot 学習で数例から分類するとき
- 知識蒸留で大規模モデルを軽量化するとき
- マルチタスク学習で複数タスクを共同学習するとき

---

## Quick Start

## 1. Vision モデルファインチューニング

```python
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import numpy as np


def finetune_vision_model(train_loader, val_loader,
                          model_name="resnet50",
                          num_classes=10, epochs=20,
                          lr=1e-4, freeze_backbone=True):
    """
    Vision モデルファインチューニング。

    Parameters:
        train_loader: DataLoader — 学習データ
        val_loader: DataLoader — 検証データ
        model_name: str — "resnet50" / "vit_b_16" / "efficientnet_b0"
        num_classes: int — クラス数
        epochs: int — エポック数
        lr: float — 学習率
        freeze_backbone: bool — バックボーン凍結
    """
    import torchvision.models as models

    # モデルロード
    model_fn = getattr(models, model_name)
    weights_name = model_name.replace("_", "").title() + "_Weights"
    try:
        weights = getattr(models, weights_name).DEFAULT
    except AttributeError:
        weights = "DEFAULT"
    model = model_fn(weights=weights)

    # 最終層置換
    if hasattr(model, "fc"):
        in_features = model.fc.in_features
        model.fc = nn.Linear(in_features, num_classes)
    elif hasattr(model, "classifier"):
        if isinstance(model.classifier, nn.Sequential):
            in_features = model.classifier[-1].in_features
            model.classifier[-1] = nn.Linear(in_features, num_classes)
        else:
            in_features = model.classifier.in_features
            model.classifier = nn.Linear(in_features, num_classes)
    elif hasattr(model, "heads"):
        in_features = model.heads.head.in_features
        model.heads.head = nn.Linear(in_features, num_classes)

    # バックボーン凍結
    if freeze_backbone:
        for name, param in model.named_parameters():
            if "fc" not in name and "classifier" not in name and "heads" not in name:
                param.requires_grad = False

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)

    optimizer = torch.optim.AdamW(
        filter(lambda p: p.requires_grad, model.parameters()), lr=lr)
    criterion = nn.CrossEntropyLoss()
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, epochs)

    best_acc = 0.0
    history = []

    for epoch in range(epochs):
        model.train()
        train_loss = 0.0
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            optimizer.zero_grad()
            outputs = model(X_batch)
            loss = criterion(outputs, y_batch)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
        scheduler.step()

        # Validation
        model.eval()
        correct = total = 0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                outputs = model(X_batch)
                _, predicted = outputs.max(1)
                total += y_batch.size(0)
                correct += predicted.eq(y_batch).sum().item()

        val_acc = correct / total
        history.append({"epoch": epoch, "train_loss": train_loss / len(train_loader),
                        "val_acc": val_acc})
        if val_acc > best_acc:
            best_acc = val_acc

    print(f"Finetune {model_name}: best val acc = {best_acc:.4f}")
    return model, history
```

## 2. NLP モデルファインチューニング

```python
def finetune_text_classifier(train_texts, train_labels,
                             val_texts, val_labels,
                             model_name="dmis-lab/biobert-base-cased-v1.2",
                             num_labels=2, epochs=5, lr=2e-5):
    """
    BERT/BioBERT テキスト分類ファインチューニング。

    Parameters:
        train_texts: list[str] — 学習テキスト
        train_labels: list[int] — 学習ラベル
        val_texts: list[str] — 検証テキスト
        val_labels: list[int] — 検証ラベル
        model_name: str — HuggingFace モデル名
        num_labels: int — ラベル数
        epochs: int — エポック数
        lr: float — 学習率
    """
    from transformers import (
        AutoTokenizer, AutoModelForSequenceClassification,
        TrainingArguments, Trainer)
    from datasets import Dataset

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForSequenceClassification.from_pretrained(
        model_name, num_labels=num_labels)

    def tokenize(examples):
        return tokenizer(examples["text"], truncation=True,
                         padding="max_length", max_length=512)

    train_ds = Dataset.from_dict({"text": train_texts, "label": train_labels})
    val_ds = Dataset.from_dict({"text": val_texts, "label": val_labels})
    train_ds = train_ds.map(tokenize, batched=True)
    val_ds = val_ds.map(tokenize, batched=True)

    args = TrainingArguments(
        output_dir="./ft_output", num_train_epochs=epochs,
        per_device_train_batch_size=16, learning_rate=lr,
        evaluation_strategy="epoch", save_strategy="epoch",
        load_best_model_at_end=True, metric_for_best_model="accuracy")

    def compute_metrics(eval_pred):
        preds = np.argmax(eval_pred.predictions, axis=-1)
        acc = (preds == eval_pred.label_ids).mean()
        return {"accuracy": acc}

    trainer = Trainer(model=model, args=args, train_dataset=train_ds,
                      eval_dataset=val_ds, compute_metrics=compute_metrics)
    trainer.train()

    metrics = trainer.evaluate()
    print(f"Finetune {model_name}: val acc = {metrics['eval_accuracy']:.4f}")
    return model, tokenizer, metrics
```

## 3. Few-shot 学習

```python
def prototypical_network(support_X, support_y, query_X,
                         feature_extractor=None):
    """
    Prototypical Network — Few-shot 分類。

    Parameters:
        support_X: np.ndarray — サポートセット特徴量
        support_y: np.ndarray — サポートラベル
        query_X: np.ndarray — クエリセット特徴量
        feature_extractor: callable | None — 特徴量抽出器
    """
    if feature_extractor is not None:
        support_emb = feature_extractor(support_X)
        query_emb = feature_extractor(query_X)
    else:
        support_emb = support_X
        query_emb = query_X

    classes = np.unique(support_y)
    prototypes = np.array([
        support_emb[support_y == c].mean(axis=0) for c in classes])

    # ユークリッド距離
    dists = np.array([
        np.linalg.norm(query_emb - p, axis=1) for p in prototypes]).T

    predictions = classes[np.argmin(dists, axis=1)]
    confidences = np.exp(-dists.min(axis=1))

    print(f"Few-shot: {len(classes)} classes, "
          f"{len(support_y)} support → {len(query_X)} query")
    return predictions, confidences
```

## 4. 知識蒸留

```python
def knowledge_distillation(teacher, student, train_loader,
                           epochs=20, temperature=4.0, alpha=0.7,
                           lr=1e-3):
    """
    知識蒸留 (Teacher → Student)。

    Parameters:
        teacher: nn.Module — 教師モデル (frozen)
        student: nn.Module — 生徒モデル
        train_loader: DataLoader — 学習データ
        epochs: int — エポック数
        temperature: float — 蒸留温度
        alpha: float — soft loss の重み
        lr: float — 学習率
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    teacher = teacher.to(device).eval()
    student = student.to(device)

    optimizer = torch.optim.AdamW(student.parameters(), lr=lr)
    ce_loss = nn.CrossEntropyLoss()
    kl_loss = nn.KLDivLoss(reduction="batchmean")

    for epoch in range(epochs):
        student.train()
        total_loss = 0.0
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)

            with torch.no_grad():
                teacher_logits = teacher(X_batch)

            student_logits = student(X_batch)

            soft_loss = kl_loss(
                nn.functional.log_softmax(student_logits / temperature, dim=1),
                nn.functional.softmax(teacher_logits / temperature, dim=1)
            ) * (temperature ** 2)

            hard_loss = ce_loss(student_logits, y_batch)
            loss = alpha * soft_loss + (1 - alpha) * hard_loss

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        print(f"  Epoch {epoch}: loss = {total_loss / len(train_loader):.4f}")

    print(f"Distillation: T={temperature}, α={alpha}, {epochs} epochs")
    return student
```

---

## パイプライン統合

```
deep-learning → transfer-learning → active-learning
  (モデル設計)    (転移・適応)       (効率的ラベル付け)
       │               │                 ↓
  healthcare-ai ───────┘       ensemble-methods
    (臨床 AI)                    (アンサンブル)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `ft_model.pt` | ファインチューニング済みモデル | → 推論 |
| `ft_history.csv` | 学習履歴 | → visualization |
| `few_shot_predictions.csv` | Few-shot 予測 | → 評価 |
