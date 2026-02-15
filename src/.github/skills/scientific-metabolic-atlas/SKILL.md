---
name: scientific-metabolic-atlas
description: |
  代謝アトラススキル。Metabolic Atlas / Human-GEM REST API による
  代謝反応・代謝産物・コンパートメント検索、フラックス解析統合、
  代謝ネットワーク可視化。K-Dense 連携: metabolic-atlas。
  ToolUniverse 連携: metabolic_atlas。
tu_tools: []
kdense_ref: metabolic-atlas
tu_tools:
  - key: metabolic_atlas
    name: Metabolic Atlas
    description: ヒトゲノム規模代謝モデル検索
---

# Scientific Metabolic Atlas

Metabolic Atlas REST API を活用したゲノムスケール代謝モデル
(GEM) 解析パイプラインを提供する。

## When to Use

- ヒト代謝反応・代謝産物を検索するとき
- Human-GEM のコンパートメント情報を取得するとき
- 代謝経路のネットワーク構造を解析するとき
- フラックスバランス解析 (FBA) の入力を準備するとき
- 代謝産物コネクティビティを可視化するとき
- 組織特異的代謝モデルを構築するとき

---

## Quick Start

## 1. 代謝反応検索

```python
import requests
import pandas as pd
import numpy as np

MA_BASE = "https://metabolicatlas.org/api/v2"


def metabolic_atlas_search_reactions(query, model="Human-GEM",
                                       compartment=None, limit=50):
    """
    Metabolic Atlas — 代謝反応検索。

    Parameters:
        query: str — 検索クエリ (例: "glycolysis", "citrate")
        model: str — GEM モデル名
        compartment: str — コンパートメント (例: "cytosol", "mitochondria")
        limit: int — 最大結果数
    """
    url = f"{MA_BASE}/search"
    params = {
        "query": query,
        "model": model,
        "type": "reaction",
        "limit": limit,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for r in data.get("results", data) if isinstance(data, dict) else data:
        rxn = r if isinstance(r, dict) else {}
        row = {
            "reaction_id": rxn.get("id", ""),
            "name": rxn.get("name", ""),
            "equation": rxn.get("equation", ""),
            "subsystem": rxn.get("subsystem", ""),
            "compartment": rxn.get("compartment", ""),
            "gene_rule": rxn.get("geneRule", ""),
            "lower_bound": rxn.get("lowerBound", None),
            "upper_bound": rxn.get("upperBound", None),
        }
        if compartment and compartment.lower() not in str(
                row.get("compartment", "")).lower():
            continue
        results.append(row)

    df = pd.DataFrame(results[:limit])
    print(f"Metabolic Atlas reactions: {len(df)} results "
          f"(query={query})")
    return df
```

## 2. 代謝産物検索

```python
def metabolic_atlas_search_metabolites(query, model="Human-GEM",
                                          limit=50):
    """
    Metabolic Atlas — 代謝産物検索。

    Parameters:
        query: str — 検索クエリ (例: "glucose", "ATP")
        model: str — GEM モデル名
        limit: int — 最大結果数
    """
    url = f"{MA_BASE}/search"
    params = {
        "query": query,
        "model": model,
        "type": "metabolite",
        "limit": limit,
    }
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for m in data.get("results", data) if isinstance(data, dict) else data:
        met = m if isinstance(m, dict) else {}
        results.append({
            "metabolite_id": met.get("id", ""),
            "name": met.get("name", ""),
            "formula": met.get("formula", ""),
            "charge": met.get("charge", None),
            "compartment": met.get("compartment", ""),
            "chebi_id": met.get("chebiId", ""),
            "kegg_id": met.get("keggId", ""),
        })

    df = pd.DataFrame(results[:limit])
    print(f"Metabolic Atlas metabolites: {len(df)} results "
          f"(query={query})")
    return df
```

## 3. 代謝ネットワーク解析

```python
import networkx as nx


def metabolic_atlas_network(subsystem, model="Human-GEM"):
    """
    Metabolic Atlas — サブシステム代謝ネットワーク構築。

    Parameters:
        subsystem: str — サブシステム名 (例: "Glycolysis")
        model: str — GEM モデル名
    """
    reactions = metabolic_atlas_search_reactions(
        subsystem, model=model, limit=200)

    G = nx.DiGraph()

    for _, rxn in reactions.iterrows():
        rxn_id = rxn["reaction_id"]
        equation = str(rxn.get("equation", ""))

        # 簡易パーサ: "A + B => C + D"
        if "=>" in equation:
            substrates_str, products_str = equation.split("=>", 1)
        elif "=" in equation:
            substrates_str, products_str = equation.split("=", 1)
        else:
            continue

        substrates = [s.strip() for s in substrates_str.split("+")
                       if s.strip()]
        products = [p.strip() for p in products_str.split("+")
                     if p.strip()]

        G.add_node(rxn_id, type="reaction",
                    name=rxn.get("name", ""))

        for s in substrates:
            G.add_node(s, type="metabolite")
            G.add_edge(s, rxn_id)

        for p in products:
            G.add_node(p, type="metabolite")
            G.add_edge(rxn_id, p)

    # ネットワーク統計
    n_reactions = sum(1 for _, d in G.nodes(data=True)
                      if d.get("type") == "reaction")
    n_metabolites = sum(1 for _, d in G.nodes(data=True)
                         if d.get("type") == "metabolite")

    print(f"Metabolic network: {n_reactions} reactions, "
          f"{n_metabolites} metabolites, {G.number_of_edges()} edges")
    return G
```

## 4. 代謝アトラス統合パイプライン

```python
def metabolic_atlas_pipeline(query, model="Human-GEM",
                               output_dir="results"):
    """
    代謝アトラス統合パイプライン。

    Parameters:
        query: str — 代謝経路/サブシステム名
        model: str — GEM モデル名
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 反応検索
    reactions = metabolic_atlas_search_reactions(query, model=model)
    reactions.to_csv(output_dir / "reactions.csv", index=False)

    # 2) 代謝産物検索
    metabolites = metabolic_atlas_search_metabolites(query, model=model)
    metabolites.to_csv(output_dir / "metabolites.csv", index=False)

    # 3) ネットワーク構築
    G = metabolic_atlas_network(query, model=model)
    nx.write_graphml(G, str(output_dir / "metabolic_network.graphml"))

    # 4) ハブ代謝産物
    met_nodes = [n for n, d in G.nodes(data=True)
                  if d.get("type") == "metabolite"]
    hub_scores = {n: G.degree(n) for n in met_nodes}
    hub_df = pd.DataFrame([
        {"metabolite": k, "degree": v}
        for k, v in sorted(hub_scores.items(),
                           key=lambda x: -x[1])[:20]
    ])
    hub_df.to_csv(output_dir / "hub_metabolites.csv", index=False)

    print(f"Metabolic Atlas pipeline: {output_dir}")
    return {
        "reactions": reactions,
        "metabolites": metabolites,
        "network": G,
        "hubs": hub_df,
    }
```

---

## K-Dense 連携

| K-Dense Key | 参照内容 |
|-------------|---------|
| `metabolic-atlas` | 代謝モデル構造・反応データベース |

## パイプライン統合

```
metabolic-modeling → metabolic-atlas → systems-biology
  (COBRA/FBA)        (Human-GEM)       (統合モデリング)
       │                   │                  ↓
  pathway-enrichment ─────┘          gene-expression
  (KEGG/Reactome)     │              (発現データ)
                       ↓
                 multi-omics
                 (メタボロミクス統合)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/reactions.csv` | 代謝反応一覧 | → metabolic-modeling |
| `results/metabolites.csv` | 代謝産物一覧 | → pathway-enrichment |
| `results/metabolic_network.graphml` | 代謝ネットワーク | → systems-biology |
| `results/hub_metabolites.csv` | ハブ代謝産物 | → multi-omics |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `metabolic_atlas` | Metabolic Atlas | ヒトゲノム規模代謝モデル検索 |
