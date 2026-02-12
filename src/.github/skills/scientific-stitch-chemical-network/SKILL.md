---
name: scientific-stitch-chemical-network
description: |
  STITCH 化学-タンパク質相互作用ネットワークスキル。STITCH
  REST API を用いた化学物質-タンパク質インタラクション検索・
  信頼度スコアリング・ネットワーク薬理学・ポリファーマコロジー解析。
  ToolUniverse 連携: stitch。
tu_tools:
  - key: stitch
    name: STITCH
    description: 化学物質-タンパク質相互作用ネットワーク (EMBL)
---

# Scientific STITCH Chemical Network

STITCH (Search Tool for Interactions of Chemicals) REST API
を活用した化学物質-タンパク質相互作用検索・信頼度スコアリング・
ネットワーク薬理学・ポリファーマコロジー解析パイプラインを提供する。

## When to Use

- 化学物質とタンパク質の相互作用エビデンスを検索するとき
- 薬物の標的タンパク質ネットワークを構築するとき
- ポリファーマコロジー (多標的薬理作用) を解析するとき
- 化学物質間の類似ネットワークを構築するとき
- ネットワーク薬理学 (Network Pharmacology) を実施するとき

---

## Quick Start

## 1. 化学物質-タンパク質相互作用検索

```python
import requests
import pandas as pd

STITCH_API = "http://stitch.embl.de/api"


def stitch_interactions(chemical, species=9606,
                          required_score=400, limit=50):
    """
    STITCH — 化学物質-タンパク質相互作用検索。

    Parameters:
        chemical: str — 化学物質名または CID
            (例: "aspirin", "CIDm00002244")
        species: int — NCBI Taxonomy ID
            (9606=ヒト)
        required_score: int — 最低信頼度スコア
            (0-1000, 400=medium)
        limit: int — 最大結果数
    """
    url = f"{STITCH_API}/tsv/interactionsList"
    params = {
        "identifiers": chemical,
        "species": species,
        "required_score": required_score,
        "limit": limit,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    lines = resp.text.strip().split("\n")
    if len(lines) < 2:
        return pd.DataFrame()

    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    df = pd.DataFrame(rows, columns=header)

    if "score" in df.columns:
        df["score"] = pd.to_numeric(
            df["score"], errors="coerce")

    print(f"STITCH: {chemical} → {len(df)} interactions")
    return df


def stitch_resolve(identifiers, species=9606):
    """
    STITCH — 化学物質/タンパク質 ID 解決。

    Parameters:
        identifiers: list[str] — 化学物質/タンパク質名リスト
        species: int — NCBI Taxonomy ID
    """
    url = f"{STITCH_API}/tsv/resolveList"
    params = {
        "identifiers": "\r".join(identifiers),
        "species": species,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    lines = resp.text.strip().split("\n")
    if len(lines) < 2:
        return pd.DataFrame()

    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    df = pd.DataFrame(rows, columns=header)

    print(f"STITCH resolve: {len(identifiers)} queries → "
          f"{len(df)} results")
    return df
```

## 2. ネットワーク薬理学

```python
def stitch_network(chemicals, species=9606,
                     required_score=400):
    """
    STITCH — 多化学物質ネットワーク構築。

    Parameters:
        chemicals: list[str] — 化学物質名リスト
        species: int — NCBI Taxonomy ID
        required_score: int — 最低信頼度スコア
    """
    url = f"{STITCH_API}/tsv/network"
    params = {
        "identifiers": "\r".join(chemicals),
        "species": species,
        "required_score": required_score,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    lines = resp.text.strip().split("\n")
    if len(lines) < 2:
        return pd.DataFrame()

    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    df = pd.DataFrame(rows, columns=header)

    nodes = set()
    if "stringId_A" in df.columns:
        nodes.update(df["stringId_A"].unique())
    if "stringId_B" in df.columns:
        nodes.update(df["stringId_B"].unique())

    print(f"STITCH network: {len(nodes)} nodes, "
          f"{len(df)} edges")
    return df


def polypharmacology_analysis(drug_list, species=9606,
                                required_score=700):
    """
    ポリファーマコロジー解析。

    Parameters:
        drug_list: list[str] — 薬物名リスト
        species: int — NCBI Taxonomy ID
        required_score: int — 高信頼度スコア閾値
    """
    all_targets = {}
    for drug in drug_list:
        interactions = stitch_interactions(
            drug, species, required_score)
        if interactions.empty:
            continue

        targets = set()
        for col in ["stringId_A", "stringId_B"]:
            if col in interactions.columns:
                targets.update(
                    interactions[col].unique())
        # 化学物質自身を除外
        targets = {t for t in targets
                   if not t.startswith("CID")}
        all_targets[drug] = targets

    # 共通標的計算
    if len(all_targets) < 2:
        return pd.DataFrame()

    pairs = []
    drugs = list(all_targets.keys())
    for i in range(len(drugs)):
        for j in range(i + 1, len(drugs)):
            shared = (all_targets[drugs[i]]
                      & all_targets[drugs[j]])
            union = (all_targets[drugs[i]]
                     | all_targets[drugs[j]])
            jaccard = (len(shared) / len(union)
                       if union else 0)
            pairs.append({
                "drug_a": drugs[i],
                "drug_b": drugs[j],
                "shared_targets": len(shared),
                "jaccard_index": jaccard,
                "shared_list": "; ".join(
                    sorted(shared)),
            })

    df = pd.DataFrame(pairs)
    df.sort_values("jaccard_index", ascending=False,
                   inplace=True)
    print(f"Polypharmacology: {len(drugs)} drugs, "
          f"{len(df)} pairs")
    return df
```

## 3. エンリッチメント解析

```python
def stitch_enrichment(identifiers, species=9606):
    """
    STITCH — 機能エンリッチメント解析。

    Parameters:
        identifiers: list[str] — タンパク質/化学物質リスト
        species: int — NCBI Taxonomy ID
    """
    url = f"{STITCH_API}/tsv/enrichment"
    params = {
        "identifiers": "\r".join(identifiers),
        "species": species,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()

    lines = resp.text.strip().split("\n")
    if len(lines) < 2:
        return pd.DataFrame()

    header = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    df = pd.DataFrame(rows, columns=header)

    if "p_value" in df.columns:
        df["p_value"] = pd.to_numeric(
            df["p_value"], errors="coerce")
        df.sort_values("p_value", inplace=True)

    print(f"STITCH enrichment: {len(df)} terms")
    return df
```

## 4. STITCH 統合パイプライン

```python
def stitch_pipeline(chemicals, species=9606,
                      output_dir="results"):
    """
    STITCH 統合パイプライン。

    Parameters:
        chemicals: list[str] — 化学物質名リスト
        species: int — NCBI Taxonomy ID
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 個別相互作用
    all_interactions = []
    for chem in chemicals:
        ixns = stitch_interactions(chem, species)
        ixns["query_chemical"] = chem
        all_interactions.append(ixns)
    ixn_df = pd.concat(all_interactions, ignore_index=True)
    ixn_df.to_csv(
        output_dir / "stitch_interactions.csv",
        index=False)

    # 2) ネットワーク
    network = stitch_network(chemicals, species)
    network.to_csv(
        output_dir / "stitch_network.csv",
        index=False)

    # 3) ポリファーマコロジー
    polypharm = polypharmacology_analysis(
        chemicals, species)
    polypharm.to_csv(
        output_dir / "polypharmacology.csv",
        index=False)

    print(f"STITCH pipeline → {output_dir}")
    return {
        "interactions": ixn_df,
        "network": network,
        "polypharmacology": polypharm,
    }
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `stitch` | STITCH | 化学物質-タンパク質相互作用 (EMBL) |

## パイプライン統合

```
cheminformatics → stitch-chemical-network → drug-target-profiling
  (化合物記述子)     (STITCH 相互作用)     (DGIdb 標的)
        │                   │                    ↓
string-network-api ────────┘          pharmacology-targets
  (STRING PPI)                        (BindingDB/GtoPdb)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/stitch_interactions.csv` | 化学物質-標的 | → drug-target-profiling |
| `results/stitch_network.csv` | ネットワーク | → string-network-api |
| `results/polypharmacology.csv` | 多標的解析 | → pharmacology-targets |
