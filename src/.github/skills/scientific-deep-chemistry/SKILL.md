---
name: scientific-deep-chemistry
description: |
  深層学習分子特性予測スキル。DeepChem による GCN/MPNN/AttentiveFP
  分子特性予測・MoleculeNet ベンチマーク・ChemBERTa/GROVER
  事前学習モデル・分子フィンガープリントフィーチャライザ。
tu_tools:
  - key: chembl
    name: ChEMBL
    description: 化学的活性・化合物データ検索
---

# Scientific Deep Chemistry

DeepChem を活用した深層学習ベース分子特性予測パイプラインを提供する。
グラフニューラルネットワーク (GCN/MPNN/AttentiveFP)、MoleculeNet
ベンチマーク、事前学習モデル (ChemBERTa/GROVER)。

## When to Use

- 分子の ADMET/物性を深層学習で予測するとき
- MoleculeNet ベンチマークデータセットを使うとき
- GCN / MPNN / AttentiveFP モデルを訓練するとき
- ChemBERTa で分子表現学習を行うとき
- 毒性予測 (Tox21, ToxCast) を行うとき
- 薬理活性予測の分子特徴量を生成するとき

---

## Quick Start

## 1. MoleculeNet データセット読込み

```python
import deepchem as dc
import numpy as np
import pandas as pd


def load_moleculenet(dataset_name="delaney", featurizer="GraphConv",
                     split="scaffold"):
    """
    MoleculeNet ベンチマークデータセット読込み。

    Parameters:
        dataset_name: str — データセット名
          ("delaney", "tox21", "bbbp", "hiv", "muv", "pcba",
           "sider", "clintox", "freesolv", "lipo")
        featurizer: str — 特徴量化手法
          ("GraphConv", "ECFP", "Weave", "MolGraphConv")
        split: str — 分割方法 ("scaffold", "random", "stratified")

    K-Dense: deepchem
    """
    loader_map = {
        "delaney": dc.molnet.load_delaney,
        "tox21": dc.molnet.load_tox21,
        "bbbp": dc.molnet.load_bbbp,
        "hiv": dc.molnet.load_hiv,
        "muv": dc.molnet.load_muv,
        "pcba": dc.molnet.load_pcba,
        "sider": dc.molnet.load_sider,
        "clintox": dc.molnet.load_clintox,
        "freesolv": dc.molnet.load_freesolv,
        "lipo": dc.molnet.load_lipo,
    }

    if dataset_name not in loader_map:
        raise ValueError(f"Unknown dataset: {dataset_name}")

    tasks, datasets, transformers = loader_map[dataset_name](
        featurizer=featurizer, splitter=split
    )
    train, valid, test = datasets

    print(f"MoleculeNet '{dataset_name}':")
    print(f"  Tasks: {len(tasks)}")
    print(f"  Train: {len(train)}, Valid: {len(valid)}, Test: {len(test)}")
    print(f"  Featurizer: {featurizer}, Split: {split}")
    return tasks, (train, valid, test), transformers
```

## 2. GCN モデル訓練

```python
def train_gcn(train_data, valid_data, tasks, n_epochs=50,
              learning_rate=0.001, batch_size=64):
    """
    Graph Convolutional Network (GCN) モデル訓練。

    Parameters:
        train_data: dc.data.Dataset — 訓練データ
        valid_data: dc.data.Dataset — 検証データ
        tasks: list — タスク名リスト
        n_epochs: int — エポック数
    """
    model = dc.models.GraphConvModel(
        n_tasks=len(tasks),
        mode="classification" if len(tasks) > 1 else "regression",
        batch_size=batch_size,
        learning_rate=learning_rate,
    )

    for epoch in range(n_epochs):
        loss = model.fit(train_data, nb_epoch=1)
        if (epoch + 1) % 10 == 0:
            metric = dc.metrics.Metric(
                dc.metrics.roc_auc_score if len(tasks) > 1
                else dc.metrics.pearson_r2_score
            )
            train_score = model.evaluate(train_data, [metric])
            valid_score = model.evaluate(valid_data, [metric])
            print(f"  Epoch {epoch+1}: "
                  f"train={list(train_score.values())[0]:.4f}, "
                  f"valid={list(valid_score.values())[0]:.4f}")

    return model
```

## 3. MPNN モデル訓練

```python
def train_mpnn(train_data, valid_data, tasks, n_epochs=50,
               learning_rate=0.001):
    """
    Message Passing Neural Network (MPNN) 訓練。

    Parameters:
        train_data: dc.data.Dataset — GraphConv 特徴量訓練データ
        valid_data: dc.data.Dataset — 検証データ
        tasks: list — タスク名リスト
    """
    model = dc.models.MPNNModel(
        n_tasks=len(tasks),
        mode="classification" if len(tasks) > 1 else "regression",
        learning_rate=learning_rate,
        node_out_feats=64,
        edge_hidden_feats=128,
        num_step_message_passing=3,
    )

    model.fit(train_data, nb_epoch=n_epochs)

    metric = dc.metrics.Metric(
        dc.metrics.roc_auc_score if len(tasks) > 1
        else dc.metrics.pearson_r2_score
    )
    valid_score = model.evaluate(valid_data, [metric])
    print(f"MPNN: valid score = {list(valid_score.values())[0]:.4f}")
    return model
```

## 4. AttentiveFP モデル訓練

```python
def train_attentivefp(train_data, valid_data, tasks, n_epochs=50,
                      learning_rate=0.001, num_layers=2):
    """
    AttentiveFP (Attention-based Fingerprint) 訓練。

    Parameters:
        train_data: dc.data.Dataset — 訓練データ
        valid_data: dc.data.Dataset — 検証データ
        tasks: list — タスク名
        num_layers: int — GATレイヤー数
    """
    model = dc.models.AttentiveFPModel(
        n_tasks=len(tasks),
        mode="classification" if len(tasks) > 1 else "regression",
        learning_rate=learning_rate,
        num_layers=num_layers,
        graph_feat_size=200,
        num_timesteps=2,
    )

    model.fit(train_data, nb_epoch=n_epochs)

    metric = dc.metrics.Metric(
        dc.metrics.roc_auc_score if len(tasks) > 1
        else dc.metrics.pearson_r2_score
    )
    valid_score = model.evaluate(valid_data, [metric])
    print(f"AttentiveFP: valid score = {list(valid_score.values())[0]:.4f}")
    return model
```

## 5. ChemBERTa 分子表現学習

```python
def chemberta_embeddings(smiles_list, model_name="seyonec/ChemBERTa-zinc-base-v1"):
    """
    ChemBERTa で SMILES → 分子埋込みベクトル。

    Parameters:
        smiles_list: list — SMILES 文字列リスト
        model_name: str — HuggingFace モデル名
    """
    from transformers import AutoTokenizer, AutoModel
    import torch

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name)
    model.eval()

    embeddings = []
    batch_size = 32

    for i in range(0, len(smiles_list), batch_size):
        batch = smiles_list[i:i+batch_size]
        inputs = tokenizer(batch, padding=True, truncation=True,
                          max_length=512, return_tensors="pt")

        with torch.no_grad():
            outputs = model(**inputs)
            # CLS トークン埋込み
            cls_emb = outputs.last_hidden_state[:, 0, :].numpy()
            embeddings.append(cls_emb)

    embeddings = np.vstack(embeddings)
    print(f"ChemBERTa: {len(smiles_list)} molecules → "
          f"{embeddings.shape[1]}D embeddings")
    return embeddings
```

## 6. モデル比較ベンチマーク

```python
def benchmark_models(dataset_name="tox21", models_to_test=None,
                     n_epochs=30):
    """
    複数モデルのベンチマーク比較。

    Parameters:
        dataset_name: str — MoleculeNet データセット
        models_to_test: list — テストモデル名
        n_epochs: int — エポック数
    """
    if models_to_test is None:
        models_to_test = ["GCN", "MPNN", "AttentiveFP"]

    results = {}

    for model_name in models_to_test:
        featurizer = "GraphConv" if model_name != "ECFP_RF" else "ECFP"
        tasks, (train, valid, test), transformers = load_moleculenet(
            dataset_name, featurizer=featurizer
        )

        is_classification = len(tasks) > 1 or dataset_name in [
            "tox21", "bbbp", "hiv", "sider", "clintox"
        ]

        if model_name == "GCN":
            model = train_gcn(train, valid, tasks, n_epochs=n_epochs)
        elif model_name == "MPNN":
            model = train_mpnn(train, valid, tasks, n_epochs=n_epochs)
        elif model_name == "AttentiveFP":
            model = train_attentivefp(train, valid, tasks, n_epochs=n_epochs)
        else:
            continue

        metric = dc.metrics.Metric(
            dc.metrics.roc_auc_score if is_classification
            else dc.metrics.pearson_r2_score
        )
        test_score = model.evaluate(test, [metric])
        results[model_name] = list(test_score.values())[0]

    print(f"\nBenchmark on '{dataset_name}':")
    for name, score in sorted(results.items(), key=lambda x: -x[1]):
        print(f"  {name}: {score:.4f}")
    return results
```

## 7. 分子特性予測パイプライン

```python
def molecular_prediction_pipeline(smiles_list, property_name="solubility",
                                   model_type="AttentiveFP"):
    """
    SMILES → 分子特性予測 統合パイプライン。

    Parameters:
        smiles_list: list — SMILES リスト
        property_name: str — 予測対象物性
        model_type: str — 使用モデル
    """
    # データセットマッピング
    property_dataset = {
        "solubility": "delaney",
        "toxicity": "tox21",
        "bbb_penetration": "bbbp",
        "hiv_activity": "hiv",
        "lipophilicity": "lipo",
        "solvation_energy": "freesolv",
    }

    dataset_name = property_dataset.get(property_name, "delaney")

    # 1) ベンチマークデータで訓練
    tasks, (train, valid, test), transformers = load_moleculenet(
        dataset_name, featurizer="GraphConv"
    )

    if model_type == "GCN":
        model = train_gcn(train, valid, tasks)
    elif model_type == "AttentiveFP":
        model = train_attentivefp(train, valid, tasks)
    else:
        model = train_mpnn(train, valid, tasks)

    # 2) 新規分子を予測
    featurizer = dc.feat.MolGraphConvFeaturizer()
    features = featurizer.featurize(smiles_list)
    pred_dataset = dc.data.NumpyDataset(X=features)
    predictions = model.predict(pred_dataset)

    results = []
    for smi, pred in zip(smiles_list, predictions):
        results.append({
            "smiles": smi,
            "prediction": float(pred[0]) if pred.ndim > 1 else float(pred),
            "property": property_name,
            "model": model_type,
        })

    df = pd.DataFrame(results)
    print(f"Predictions: {len(df)} molecules, property='{property_name}'")
    return df
```

---

## パイプライン統合

```
cheminformatics → deep-chemistry → drug-target-profiling
  (RDKit/SMILES)   (GCN/MPNN/FP)   (ChEMBL/標的)
       │                  │                ↓
molecular-docking ───────┘         admet-pharmacokinetics
  (AutoDock/Vina)      │            (ADMET予測)
                       ↓
                  md-simulation
                  (分子動力学検証)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/predictions.csv` | 分子特性予測値 | → drug-target-profiling |
| `results/benchmark.json` | モデルベンチマーク結果 | — |
| `results/embeddings.npy` | ChemBERTa 埋込み | → cheminformatics |
| `results/model/` | 訓練済みモデル | → admet-pharmacokinetics |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `chembl` | ChEMBL | 化学的活性・化合物データ検索 |
