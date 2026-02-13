---
name: scientific-neural-architecture-search
description: |
  ニューラルアーキテクチャ探索 (NAS) スキル。DARTS 微分可能 NAS・
  Optuna NAS 統合・効率的ネットワーク設計・探索空間定義・
  Pareto 最適化 (精度 vs FLOPS)。
---

# Scientific Neural Architecture Search

ニューラルネットワーク構造の自動探索・最適化パイプラインを提供する。

## When to Use

- NN アーキテクチャを自動的に最適化したいとき
- 精度と計算コストの Pareto 最適解を探索するとき
- 探索空間を定義して効率的な構造探索を行うとき
- DARTS 系の微分可能 NAS を実装するとき
- モデル圧縮・効率化を自動化するとき

---

## Quick Start

## 1. Optuna NAS — ネットワーク構造探索

```python
import optuna
import torch
import torch.nn as nn
from typing import Dict, Any


def optuna_nas(train_loader, val_loader, search_space=None,
               n_trials=50, n_epochs=10, device="cpu",
               direction="maximize", metric="accuracy"):
    """
    Optuna によるニューラルアーキテクチャ探索。

    Parameters:
        train_loader: DataLoader — 訓練データ
        val_loader: DataLoader — 検証データ
        search_space: dict | None — カスタム探索空間
        n_trials: int — 試行回数
        n_epochs: int — 各試行の学習エポック数
        device: str — "cpu" / "cuda"
        direction: str — "maximize" / "minimize"
        metric: str — 最適化指標名
    """
    input_dim = next(iter(train_loader))[0].shape[1]
    n_classes = len(torch.unique(
        torch.cat([y for _, y in train_loader])))

    def build_model(trial):
        n_layers = trial.suggest_int("n_layers", 1, 5)
        layers = []
        in_features = input_dim

        for i in range(n_layers):
            out_features = trial.suggest_int(
                f"n_units_l{i}", 16, 512, log=True)
            layers.append(nn.Linear(in_features, out_features))

            activation = trial.suggest_categorical(
                f"activation_l{i}", ["relu", "gelu", "silu"])
            act_map = {"relu": nn.ReLU(), "gelu": nn.GELU(),
                       "silu": nn.SiLU()}
            layers.append(act_map[activation])

            dropout = trial.suggest_float(
                f"dropout_l{i}", 0.0, 0.5)
            if dropout > 0:
                layers.append(nn.Dropout(dropout))

            use_bn = trial.suggest_categorical(
                f"batchnorm_l{i}", [True, False])
            if use_bn:
                layers.append(nn.BatchNorm1d(out_features))

            in_features = out_features

        layers.append(nn.Linear(in_features, n_classes))
        return nn.Sequential(*layers)

    def objective(trial):
        model = build_model(trial).to(device)
        lr = trial.suggest_float("lr", 1e-5, 1e-1, log=True)
        optimizer_name = trial.suggest_categorical(
            "optimizer", ["Adam", "AdamW", "SGD"])

        opt_cls = getattr(torch.optim, optimizer_name)
        optimizer = opt_cls(model.parameters(), lr=lr)
        criterion = nn.CrossEntropyLoss()

        for epoch in range(n_epochs):
            model.train()
            for X, y in train_loader:
                X, y = X.to(device), y.to(device)
                optimizer.zero_grad()
                loss = criterion(model(X), y)
                loss.backward()
                optimizer.step()

            # 枝刈り
            model.eval()
            correct, total = 0, 0
            with torch.no_grad():
                for X, y in val_loader:
                    X, y = X.to(device), y.to(device)
                    correct += (model(X).argmax(1) == y).sum().item()
                    total += len(y)

            trial.report(correct / total, epoch)
            if trial.should_prune():
                raise optuna.TrialPruned()

        return correct / total

    study = optuna.create_study(
        direction=direction,
        pruner=optuna.pruners.MedianPruner(n_warmup_steps=3))
    study.optimize(objective, n_trials=n_trials)

    print(f"Best {metric}: {study.best_value:.4f}")
    print(f"Best params: {study.best_params}")
    return study


def nas_pareto_search(train_loader, val_loader,
                      n_trials=100, device="cpu"):
    """
    多目的 NAS — 精度 vs モデルサイズの Pareto 最適化。

    Parameters:
        train_loader: DataLoader — 訓練データ
        val_loader: DataLoader — 検証データ
        n_trials: int — 試行回数
        device: str — デバイス
    """

    def objective(trial):
        n_layers = trial.suggest_int("n_layers", 1, 4)
        total_params = 0
        in_f = next(iter(train_loader))[0].shape[1]
        n_classes = len(torch.unique(
            torch.cat([y for _, y in val_loader])))

        layers = []
        for i in range(n_layers):
            out_f = trial.suggest_int(f"units_{i}", 8, 256, log=True)
            layers.append(nn.Linear(in_f, out_f))
            layers.append(nn.ReLU())
            total_params += in_f * out_f + out_f
            in_f = out_f
        layers.append(nn.Linear(in_f, n_classes))
        total_params += in_f * n_classes + n_classes

        model = nn.Sequential(*layers).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
        criterion = nn.CrossEntropyLoss()

        for _ in range(5):
            model.train()
            for X, y in train_loader:
                X, y = X.to(device), y.to(device)
                optimizer.zero_grad()
                criterion(model(X), y).backward()
                optimizer.step()

        model.eval()
        correct, total = 0, 0
        with torch.no_grad():
            for X, y in val_loader:
                X, y = X.to(device), y.to(device)
                correct += (model(X).argmax(1) == y).sum().item()
                total += len(y)

        return correct / total, total_params

    study = optuna.create_study(
        directions=["maximize", "minimize"])
    study.optimize(objective, n_trials=n_trials)

    print(f"Pareto front: {len(study.best_trials)} solutions")
    return study
```

---

## パイプライン統合

```
[タスク定義] → neural-architecture-search → deep-learning
                     (構造探索)                 (本格学習)
                          │
                   automl ← ensemble-methods
                  (HPO)      (アンサンブル)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `nas_study.pkl` | Optuna Study | → 最良構造抽出 |
| `pareto_front.csv` | Pareto 最適解群 | → モデル選択 |
| `best_architecture.json` | 最良アーキテクチャ | → deep-learning |
