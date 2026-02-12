---
name: scientific-healthcare-ai
description: |
  ヘルスケア AI スキル。PyHealth 臨床 ML パイプライン、
  フローサイトメトリー (FlowIO) 解析、電子健康記録 (EHR) 処理、
  臨床予測モデル構築のガイダンス。
---

# Scientific Healthcare AI

臨床データ解析・ヘルスケア機械学習パイプラインを提供する。
PyHealth フレームワーク、フローサイトメトリー解析ツールを活用。

## When to Use

- 臨床予測モデル (再入院予測, 死亡率予測等) を構築するとき
- EHR (電子健康記録) データの前処理・特徴量エンジニアリングを行うとき
- フローサイトメトリー (FACS) データを読み込み・解析するとき
- 臨床タスク向けの ML パイプラインを設計するとき
- 医療コード (ICD-10, SNOMED, ATC) のマッピングを行うとき

---

## Quick Start

## 1. PyHealth 臨床予測パイプライン

```python
"""
PyHealth による臨床予測モデル構築。
pip install pyhealth

K-Dense-AI 参照: pyhealth — 臨床 ML フレームワーク
"""
from pyhealth.datasets import MIMIC3Dataset
from pyhealth.tasks import readmission_prediction_mimic3_fn


def build_clinical_pipeline(
    mimic3_root,
    tables=("DIAGNOSES_ICD", "PROCEDURES_ICD", "PRESCRIPTIONS"),
    code_mapping=None,
):
    """
    MIMIC-III データセットから臨床予測パイプラインを構築。

    Parameters:
        mimic3_root: str — MIMIC-III CSV ディレクトリパス
        tables: tuple — 使用するテーブル
        code_mapping: dict | None — コードマッピング設定
    """
    # Step 1: Dataset loading
    if code_mapping is None:
        code_mapping = {
            "NDC": ("ATC", {"target_kwargs": {"level": 3}}),
            "ICD9CM": "CCSCM",
            "ICD9PROC": "CCSPROC",
        }

    dataset = MIMIC3Dataset(
        root=mimic3_root,
        tables=tables,
        code_mapping=code_mapping,
    )

    print(f"MIMIC-III dataset: {len(dataset.patients)} patients")
    return dataset


def apply_clinical_task(dataset, task_fn=None):
    """
    臨床タスク関数を適用しサンプルを生成。
    """
    from pyhealth.datasets import split_by_patient

    if task_fn is None:
        task_fn = readmission_prediction_mimic3_fn

    samples = dataset.set_task(task_fn)
    train, val, test = split_by_patient(samples, [0.8, 0.1, 0.1])

    print(f"Clinical task samples: "
          f"train={len(train)}, val={len(val)}, test={len(test)}")
    return train, val, test
```

## 2. PyHealth モデル学習

```python
def train_clinical_model(
    train_dataset,
    val_dataset,
    model_type="Transformer",
    epochs=20,
    batch_size=64,
):
    """
    PyHealth モデルの学習。

    Parameters:
        train_dataset: SampleDataset
        val_dataset: SampleDataset
        model_type: str — "Transformer", "RETAIN", "GRU", "CNN"
        epochs: int — 学習エポック数
    """
    from pyhealth.models import Transformer
    from pyhealth.trainer import Trainer

    model_classes = {
        "Transformer": Transformer,
    }
    ModelClass = model_classes.get(model_type, Transformer)

    model = ModelClass(
        dataset=train_dataset,
        feature_keys=["conditions", "procedures", "drugs"],
        label_key="label",
        mode="binary",
    )

    trainer = Trainer(model=model)
    trainer.train(
        train_dataloader=train_dataset,
        val_dataloader=val_dataset,
        epochs=epochs,
        monitor="pr_auc",
    )

    print(f"Clinical model ({model_type}): trained for {epochs} epochs")
    return model, trainer
```

## 3. フローサイトメトリー解析

```python
def read_fcs_file(fcs_path):
    """
    FCS ファイルの読み込みと前処理。
    pip install flowio

    K-Dense-AI 参照: flowio — FCS file I/O

    Parameters:
        fcs_path: str — FCS ファイルパス
    """
    import flowio
    import numpy as np
    import pandas as pd

    fcs_data = flowio.FlowData(fcs_path)

    # Extract channel names
    channels = []
    for i in range(1, fcs_data.channel_count + 1):
        name = fcs_data.channels.get(f"P{i}N", f"Channel_{i}")
        channels.append(name)

    # Convert to DataFrame
    events = np.array(fcs_data.events).reshape(-1, fcs_data.channel_count)
    df = pd.DataFrame(events, columns=channels)

    print(f"FCS '{fcs_path}': {len(df)} events x {len(channels)} channels")
    return df, fcs_data


def gate_fcs_data(df, channel, low=None, high=None):
    """
    単純な矩形ゲーティング。

    Parameters:
        df: pd.DataFrame — FCS データ
        channel: str — チャネル名
        low: float | None — 下限
        high: float | None — 上限
    """
    mask = pd.Series([True] * len(df))
    if low is not None:
        mask &= df[channel] >= low
    if high is not None:
        mask &= df[channel] <= high

    gated = df[mask]
    pct = len(gated) / len(df) * 100
    print(f"Gate '{channel}' [{low},{high}]: "
          f"{len(gated)}/{len(df)} events ({pct:.1f}%)")
    return gated
```

## 4. 医療コードマッピング

```python
def map_medical_codes(codes, source_system, target_system):
    """
    医療コード間のマッピング。

    Parameters:
        codes: list[str] — ソースコードのリスト
        source_system: str — "ICD9CM", "ICD10CM", "NDC", "ATC", "SNOMED"
        target_system: str — 変換先コード体系
    """
    try:
        from pyhealth.medcode import CrossMap

        mapper = CrossMap(source_system, target_system)
        results = {}
        for code in codes:
            mapped = mapper.map(code)
            results[code] = mapped

        mapped_count = sum(1 for v in results.values() if v)
        print(f"Code mapping {source_system}→{target_system}: "
              f"{mapped_count}/{len(codes)} mapped")
        return results

    except ImportError:
        print("pyhealth.medcode not available; install pyhealth")
        return {}
```

## 5. 臨床モデル評価

```python
def evaluate_clinical_model(trainer, test_dataset):
    """
    臨床予測モデルの評価。

    Parameters:
        trainer: Trainer — 学習済み Trainer
        test_dataset: SampleDataset — テストデータ
    """
    metrics = trainer.evaluate(test_dataset)

    print("Clinical model evaluation:")
    for metric_name, value in metrics.items():
        print(f"  {metric_name}: {value:.4f}")

    return metrics
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/clinical_predictions.csv` | CSV |
| `results/clinical_metrics.json` | JSON |
| `results/fcs_processed.csv` | CSV |
| `results/code_mapping.json` | JSON |

### 利用可能ツール

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| (K-Dense) | `pyhealth` | 臨床 ML フレームワーク |
| (K-Dense) | `flowio` | FCS ファイル I/O |

> **注**: 本スキルは ToolUniverse ツールを持たず、
> K-Dense-AI Scientific Skills からの参照のみ。

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-clinical-nlp` | 臨床 NLP |
| `scientific-biostatistics-survival` | 生存時間解析 |
| `scientific-single-cell-rnaseq` | 単一細胞解析 |
| `scientific-machine-learning-omics` | ML x オミクス |
| `scientific-biothings-idmapping` | ID マッピング |

### 依存パッケージ

`pyhealth`, `flowio`, `numpy`, `pandas`, `scikit-learn`
