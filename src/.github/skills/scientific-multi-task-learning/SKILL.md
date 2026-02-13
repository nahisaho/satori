---
name: scientific-multi-task-learning
description: |
  マルチタスク学習スキル。Hard/Soft Parameter Sharing・
  GradNorm 勾配正規化・PCGrad 勾配投影・
  タスクバランシング・補助タスク設計。
---

# Scientific Multi-Task Learning

複数の関連タスクを同時に学習し、共有表現を活用して
各タスクの汎化性能を向上させるパイプラインを提供する。

## When to Use

- 複数の関連予測タスクを同時に実行するとき
- 共有表現を学習してデータ効率を高めたいとき
- 主タスク + 補助タスクの構成で学習するとき
- タスク間の勾配干渉を解消するとき
- マルチ出力回帰・分類を設計するとき

---

## Quick Start

## 1. Hard Parameter Sharing MTL

```python
import torch
import torch.nn as nn
from typing import Dict, List, Tuple


class HardSharingMTL(nn.Module):
    """
    Hard Parameter Sharing マルチタスクモデル。

    共有エンコーダ + タスク別ヘッドの構成。
    """

    def __init__(self, input_dim, shared_dims, task_configs):
        """
        Parameters:
            input_dim: int — 入力次元
            shared_dims: list[int] — 共有層のユニット数
            task_configs: dict — {task_name: {"output_dim": int, "head_dims": [int]}}
        """
        super().__init__()

        # 共有エンコーダ
        layers = []
        in_d = input_dim
        for d in shared_dims:
            layers.extend([nn.Linear(in_d, d), nn.ReLU(),
                           nn.BatchNorm1d(d), nn.Dropout(0.2)])
            in_d = d
        self.shared_encoder = nn.Sequential(*layers)

        # タスク別ヘッド
        self.task_heads = nn.ModuleDict()
        for name, config in task_configs.items():
            head_layers = []
            h_in = in_d
            for h_d in config.get("head_dims", [64]):
                head_layers.extend([nn.Linear(h_in, h_d), nn.ReLU()])
                h_in = h_d
            head_layers.append(nn.Linear(h_in, config["output_dim"]))
            self.task_heads[name] = nn.Sequential(*head_layers)

    def forward(self, x):
        shared = self.shared_encoder(x)
        return {name: head(shared) for name, head in self.task_heads.items()}


def train_mtl_model(model, train_loader, task_losses,
                    task_weights=None, epochs=50,
                    lr=1e-3, device="cpu"):
    """
    MTL モデルの学習。

    Parameters:
        model: HardSharingMTL — MTL モデル
        train_loader: DataLoader — {task_name: (X, y)} バッチ
        task_losses: dict — {task_name: loss_fn}
        task_weights: dict | None — {task_name: float} タスク重み
        epochs: int — 学習エポック数
        lr: float — 学習率
        device: str — デバイス
    """
    import pandas as pd

    if task_weights is None:
        task_weights = {name: 1.0 for name in task_losses}

    model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    history = []

    for epoch in range(epochs):
        model.train()
        epoch_losses = {name: 0.0 for name in task_losses}

        for batch in train_loader:
            X = batch["X"].to(device)
            outputs = model(X)
            optimizer.zero_grad()

            total_loss = 0
            for name, loss_fn in task_losses.items():
                y = batch[name].to(device)
                task_loss = loss_fn(outputs[name], y)
                total_loss += task_weights[name] * task_loss
                epoch_losses[name] += task_loss.item()

            total_loss.backward()
            optimizer.step()

        record = {"epoch": epoch + 1}
        for name in task_losses:
            record[f"loss_{name}"] = epoch_losses[name] / len(train_loader)
        history.append(record)

        if (epoch + 1) % 10 == 0:
            losses_str = " | ".join(
                f"{n}={epoch_losses[n]/len(train_loader):.4f}"
                for n in task_losses)
            print(f"Epoch {epoch+1}: {losses_str}")

    return pd.DataFrame(history)
```

## 2. GradNorm — 動的タスクバランシング

```python
def gradnorm_balance(model, task_losses, train_loader,
                     alpha=1.5, epochs=50, lr=1e-3, device="cpu"):
    """
    GradNorm による動的タスク重みバランシング。

    Parameters:
        model: HardSharingMTL — MTL モデル
        task_losses: dict — {task_name: loss_fn}
        train_loader: DataLoader
        alpha: float — GradNorm 非対称度パラメータ
        epochs: int — 学習エポック
        lr: float — 学習率
        device: str — デバイス
    """
    import pandas as pd

    task_names = list(task_losses.keys())
    n_tasks = len(task_names)
    log_weights = torch.zeros(n_tasks, requires_grad=True, device=device)

    model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    weight_optimizer = torch.optim.Adam([log_weights], lr=0.025)

    initial_losses = None
    history = []

    for epoch in range(epochs):
        model.train()
        epoch_losses = {n: 0.0 for n in task_names}

        for batch in train_loader:
            X = batch["X"].to(device)
            outputs = model(X)

            weights = torch.softmax(log_weights, dim=0) * n_tasks
            losses = []
            for i, name in enumerate(task_names):
                y = batch[name].to(device)
                task_loss = task_losses[name](outputs[name], y)
                losses.append(task_loss)
                epoch_losses[name] += task_loss.item()

            if initial_losses is None:
                initial_losses = [l.item() for l in losses]

            total_loss = sum(w * l for w, l in zip(weights, losses))

            optimizer.zero_grad()
            weight_optimizer.zero_grad()
            total_loss.backward(retain_graph=True)

            # GradNorm 更新
            shared_params = list(model.shared_encoder.parameters())
            norms = []
            for l in losses:
                g = torch.autograd.grad(l, shared_params[-1],
                                        retain_graph=True)[0]
                norms.append(torch.norm(g))

            avg_norm = torch.stack(norms).mean()
            loss_ratios = torch.tensor(
                [l.item() / il for l, il in
                 zip(losses, initial_losses)], device=device)
            relative_inv = loss_ratios / loss_ratios.mean()
            target_norms = avg_norm * (relative_inv ** alpha)

            gradnorm_loss = sum(
                torch.abs(n - t) for n, t in
                zip(norms, target_norms))
            gradnorm_loss.backward()

            optimizer.step()
            weight_optimizer.step()

        record = {"epoch": epoch + 1}
        for i, name in enumerate(task_names):
            record[f"loss_{name}"] = epoch_losses[name] / len(train_loader)
            record[f"weight_{name}"] = (
                torch.softmax(log_weights, 0) * n_tasks)[i].item()
        history.append(record)

    return pd.DataFrame(history)
```

---

## パイプライン統合

```
[複数タスク定義] → multi-task-learning → feature-importance
                     (共有表現学習)        (特徴量解釈)
                          │
                   deep-learning ← transfer-learning
                    (基盤 NN)        (転移学習)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `mtl_model.pt` | MTL モデル | → 推論 |
| `mtl_history.csv` | タスク別学習履歴 | → 可視化 |
| `gradnorm_weights.csv` | 動的タスク重み推移 | → バランシング分析 |
