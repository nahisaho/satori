---
name: scientific-federated-learning
description: |
  連合学習スキル。Flower フレームワークによる FL パイプライン・
  FedAvg/FedProx/FedOpt 集約戦略・差分プライバシー (DP-SGD)・
  非 IID データ分割・通信効率化。
---

# Scientific Federated Learning

プライバシー保護型分散機械学習を実現する連合学習パイプラインを提供する。

## When to Use

- 複数施設・組織のデータを集約せずにモデル学習するとき
- 医療データ・個人情報を含むデータで ML を行うとき
- 差分プライバシーを適用した学習が必要なとき
- 非 IID データ分割下での連合学習を設計するとき
- 通信効率を考慮した分散学習を構築するとき

---

## Quick Start

## 1. Flower 連合学習パイプライン

```python
import flwr as fl
import numpy as np
from typing import Dict, List, Tuple, Optional


def create_fl_client(model, train_loader, val_loader,
                     device="cpu"):
    """
    Flower クライアント生成。

    Parameters:
        model: nn.Module — PyTorch モデル
        train_loader: DataLoader — 訓練データ
        val_loader: DataLoader — 検証データ
        device: str — "cpu" / "cuda"
    """
    import torch

    class SatoriFlClient(fl.client.NumPyClient):
        def get_parameters(self, config):
            return [val.cpu().numpy()
                    for val in model.parameters()]

        def set_parameters(self, parameters):
            for param, new_val in zip(model.parameters(), parameters):
                param.data = torch.tensor(new_val).to(device)

        def fit(self, parameters, config):
            self.set_parameters(parameters)
            model.train()
            optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
            criterion = torch.nn.CrossEntropyLoss()

            epochs = config.get("local_epochs", 1)
            for _ in range(epochs):
                for X, y in train_loader:
                    X, y = X.to(device), y.to(device)
                    optimizer.zero_grad()
                    loss = criterion(model(X), y)
                    loss.backward()
                    optimizer.step()

            return self.get_parameters(config), len(train_loader.dataset), {}

        def evaluate(self, parameters, config):
            self.set_parameters(parameters)
            model.eval()
            criterion = torch.nn.CrossEntropyLoss()
            total_loss, correct, total = 0.0, 0, 0

            with torch.no_grad():
                for X, y in val_loader:
                    X, y = X.to(device), y.to(device)
                    preds = model(X)
                    total_loss += criterion(preds, y).item() * len(y)
                    correct += (preds.argmax(1) == y).sum().item()
                    total += len(y)

            return total_loss / total, total, {"accuracy": correct / total}

    return SatoriFlClient()


def create_fl_strategy(algorithm="fedavg", min_clients=2,
                       fraction_fit=1.0, fraction_evaluate=1.0,
                       proximal_mu=0.1):
    """
    連合学習集約戦略の選択。

    Parameters:
        algorithm: str — "fedavg" / "fedprox" / "fedopt" / "fedadam"
        min_clients: int — 最小クライアント数
        fraction_fit: float — 学習参加率
        fraction_evaluate: float — 評価参加率
        proximal_mu: float — FedProx 近接項の強度
    """
    common = dict(
        min_fit_clients=min_clients,
        min_evaluate_clients=min_clients,
        min_available_clients=min_clients,
        fraction_fit=fraction_fit,
        fraction_evaluate=fraction_evaluate,
    )

    strategies = {
        "fedavg": fl.server.strategy.FedAvg(**common),
        "fedprox": fl.server.strategy.FedProx(
            proximal_mu=proximal_mu, **common),
        "fedadam": fl.server.strategy.FedAdam(
            eta=1e-1, eta_l=1e-1, tau=1e-9, **common),
    }

    strategy = strategies.get(algorithm, strategies["fedavg"])
    print(f"FL Strategy: {algorithm} | min_clients={min_clients}")
    return strategy
```

## 2. 差分プライバシー (DP-SGD)

```python
def apply_differential_privacy(model, train_loader,
                               target_epsilon=1.0,
                               target_delta=1e-5,
                               max_grad_norm=1.0,
                               noise_multiplier=1.1,
                               epochs=10, lr=1e-3):
    """
    Opacus DP-SGD による差分プライバシー学習。

    Parameters:
        model: nn.Module — PyTorch モデル
        train_loader: DataLoader — 訓練データ
        target_epsilon: float — プライバシーバジェット ε
        target_delta: float — プライバシーパラメータ δ
        max_grad_norm: float — 勾配クリッピングノルム
        noise_multiplier: float — ノイズ乗数 σ
        epochs: int — 学習エポック数
        lr: float — 学習率
    """
    import torch
    from opacus import PrivacyEngine

    optimizer = torch.optim.SGD(model.parameters(), lr=lr)
    privacy_engine = PrivacyEngine()

    model, optimizer, train_loader = privacy_engine.make_private_with_epsilon(
        module=model,
        optimizer=optimizer,
        data_loader=train_loader,
        epochs=epochs,
        target_epsilon=target_epsilon,
        target_delta=target_delta,
        max_grad_norm=max_grad_norm,
    )

    criterion = torch.nn.CrossEntropyLoss()
    history = []

    for epoch in range(epochs):
        model.train()
        total_loss = 0
        for X, y in train_loader:
            optimizer.zero_grad()
            loss = criterion(model(X), y)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        epsilon = privacy_engine.get_epsilon(delta=target_delta)
        history.append({"epoch": epoch + 1,
                        "loss": total_loss / len(train_loader),
                        "epsilon": epsilon})
        print(f"Epoch {epoch+1}: loss={total_loss/len(train_loader):.4f}, "
              f"ε={epsilon:.2f}")

    import pandas as pd
    return pd.DataFrame(history)
```

## 3. 非 IID データ分割

```python
def create_non_iid_splits(dataset_labels, n_clients=5,
                          alpha=0.5, seed=42):
    """
    Dirichlet 分布ベースの非 IID データ分割。

    Parameters:
        dataset_labels: np.ndarray — 全データのラベル配列
        n_clients: int — クライアント数
        alpha: float — Dirichlet α (小さいほど偏りが大きい)
        seed: int — 乱数シード
    """
    rng = np.random.default_rng(seed)
    n_classes = len(np.unique(dataset_labels))
    client_indices = [[] for _ in range(n_clients)]

    for c in range(n_classes):
        class_idx = np.where(dataset_labels == c)[0]
        proportions = rng.dirichlet(np.repeat(alpha, n_clients))
        split_points = (np.cumsum(proportions) * len(class_idx)).astype(int)
        splits = np.split(class_idx, split_points[:-1])
        for i, split in enumerate(splits):
            client_indices[i].extend(split.tolist())

    # 分布サマリー
    for i, indices in enumerate(client_indices):
        labels = dataset_labels[indices]
        unique, counts = np.unique(labels, return_counts=True)
        dist = dict(zip(unique.tolist(), counts.tolist()))
        print(f"Client {i}: {len(indices)} samples, dist={dist}")

    return client_indices
```

---

## パイプライン統合

```
[プライバシー要件] → federated-learning → model-monitoring
                      (連合学習)              (モデル監視)
                           │
                    deep-learning ← transfer-learning
                      (基盤 NN)       (転移学習)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `fl_strategy_config.json` | FL 集約設定 | → サーバー起動 |
| `dp_training_history.csv` | DP 学習履歴 | → model-monitoring |
| `client_splits.json` | 非 IID 分割情報 | → FL クライアント |
