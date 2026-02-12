---
name: scientific-protein-interaction-network
description: |
  タンパク質-タンパク質相互作用 (PPI) ネットワーク解析スキル。STRING、IntAct、
  BioGRID、STITCH (化学-タンパク質) 相互作用データベースを統合した
  ネットワーク構築・解析パイプライン。GO/KEGG 富化、相互作用パートナー発見、
  組織特異的ネットワーク (HumanBase)、化合物-標的ネットワーク対応。
  14 の ToolUniverse SMCP ツールと連携。
tu_tools:
  - key: intact
    name: IntAct
    description: 分子相互作用データベース (EBI)
---

# Scientific Protein Interaction Network

STRING / IntAct / BioGRID / STITCH の 4 大 PPI データベースを統合した
タンパク質相互作用ネットワーク解析パイプラインを提供する。

## When to Use

- DEG や変異遺伝子の PPI ネットワークを構築するとき
- ハブタンパク質やボトルネックの特定が必要なとき
- 化合物と標的タンパク質の相互作用を調べるとき
- 組織特異的な相互作用ネットワークを評価するとき
- PPI データに基づく GO/KEGG 富化やモジュール解析を行うとき

---

## Quick Start

## 1. STRING PPI ネットワーク取得

```python
import requests
import pandas as pd
import networkx as nx


def string_get_interactions(proteins, species=9606,
                             score_threshold=400,
                             network_type="functional"):
    """
    STRING API v12 による PPI ネットワーク取得。

    Parameters:
        proteins: list — タンパク質名/UniProt ID リスト
        species: int — NCBI Taxonomy ID (9606=Homo sapiens)
        score_threshold: int — 信頼スコア閾値 (0-1000)
        network_type: "functional" or "physical"
    """
    base = "https://string-db.org/api/json"

    # タンパク質 ID 解決
    resolve_url = f"{base}/get_string_ids"
    resolved = []
    for batch in [proteins[i:i+10] for i in range(0, len(proteins), 10)]:
        params = {
            "identifiers": "\r".join(batch),
            "species": species,
            "limit": 1,
        }
        resp = requests.get(resolve_url, params=params)
        for r in resp.json():
            resolved.append(r["stringId"])

    if not resolved:
        print("No proteins resolved")
        return pd.DataFrame(), nx.Graph()

    # 相互作用取得
    interaction_url = f"{base}/network"
    params = {
        "identifiers": "\r".join(resolved),
        "species": species,
        "required_score": score_threshold,
        "network_type": network_type,
    }
    resp = requests.get(interaction_url, params=params)
    interactions = resp.json()

    edges = []
    for i in interactions:
        edges.append({
            "protein_a": i["preferredName_A"],
            "protein_b": i["preferredName_B"],
            "score": i["score"],
            "nscore": i.get("nscore", 0),
            "fscore": i.get("fscore", 0),
            "pscore": i.get("pscore", 0),
            "ascore": i.get("ascore", 0),
            "escore": i.get("escore", 0),
            "dscore": i.get("dscore", 0),
            "tscore": i.get("tscore", 0),
        })

    df = pd.DataFrame(edges)

    # NetworkX グラフ構築
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(row["protein_a"], row["protein_b"],
                   weight=row["score"] / 1000.0)

    print(f"STRING network: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges (score ≥ {score_threshold})")
    return df, G
```

## 2. IntAct 分子相互作用検索

```python
def intact_search_interactions(query, species="human",
                                interaction_type=None,
                                max_results=200):
    """
    IntAct REST API による分子相互作用検索。

    Parameters:
        query: str — タンパク質名/UniProt ID
        species: str or int — "human" or taxonomy ID
        interaction_type: str — MI term (e.g., "MI:0407" physical association)
    """
    url = "https://www.ebi.ac.uk/intact/ws/interaction/findInteractions"
    params = {
        "query": query,
        "maxResults": max_results,
    }
    if species:
        params["species"] = species

    resp = requests.get(url, params=params)
    if resp.status_code != 200:
        print(f"IntAct error: {resp.status_code}")
        return pd.DataFrame()

    data = resp.json()
    interactions = data.get("content", [])

    results = []
    for ix in interactions:
        interactor_a = ix.get("interactorA", {})
        interactor_b = ix.get("interactorB", {})
        results.append({
            "interactor_a": interactor_a.get("preferredIdentifier", ""),
            "interactor_a_name": interactor_a.get("shortLabel", ""),
            "interactor_b": interactor_b.get("preferredIdentifier", ""),
            "interactor_b_name": interactor_b.get("shortLabel", ""),
            "interaction_type": ix.get("interactionType", ""),
            "detection_method": ix.get("detectionMethod", ""),
            "confidence": ix.get("confidenceValue", 0),
            "publication": ix.get("pubmedId", ""),
        })

    df = pd.DataFrame(results)
    print(f"IntAct: {len(df)} interactions for '{query}'")
    return df
```

## 3. STITCH 化合物-タンパク質相互作用

```python
def stitch_chemical_protein(chemicals, species=9606,
                             score_threshold=400):
    """
    STITCH API による化合物-タンパク質相互作用検索。

    Parameters:
        chemicals: list — 化合物名/CID リスト
        species: int — NCBI Taxonomy ID
        score_threshold: int — 信頼スコア閾値
    """
    url = "http://stitch.embl.de/api/json/interactionsList"
    params = {
        "identifiers": "\r".join(chemicals),
        "species": species,
        "required_score": score_threshold,
    }

    resp = requests.get(url, params=params)
    interactions = resp.json()

    results = []
    for i in interactions:
        results.append({
            "chemical": i.get("preferredName_A", ""),
            "protein": i.get("preferredName_B", ""),
            "score": i.get("score", 0),
            "type_a": "chemical" if i.get("ncbiTaxonId_A") == -1 else "protein",
        })

    df = pd.DataFrame(results)
    print(f"STITCH: {len(df)} chemical-protein interactions")
    return df
```

## 4. PPI ネットワーク解析 (中心性・モジュール)

```python
def analyze_ppi_network(G, community_method="louvain"):
    """
    PPI ネットワークのトポロジー解析。

    Parameters:
        G: nx.Graph — PPI ネットワーク
        community_method: "louvain" or "label_propagation"
    """
    if G.number_of_nodes() == 0:
        return {}

    # 中心性指標
    degree_cent = nx.degree_centrality(G)
    betweenness = nx.betweenness_centrality(G)
    closeness = nx.closeness_centrality(G)

    # ハブタンパク質 (degree top 10)
    hubs = sorted(degree_cent.items(), key=lambda x: -x[1])[:10]

    # ボトルネック (betweenness top 10)
    bottlenecks = sorted(betweenness.items(), key=lambda x: -x[1])[:10]

    # コミュニティ検出
    if community_method == "louvain":
        from community import community_louvain
        partition = community_louvain.best_partition(G)
    else:
        communities = nx.community.label_propagation_communities(G)
        partition = {}
        for i, comm in enumerate(communities):
            for node in comm:
                partition[node] = i

    n_communities = len(set(partition.values()))

    stats = {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "density": round(nx.density(G), 4),
        "avg_clustering": round(nx.average_clustering(G), 4),
        "connected_components": nx.number_connected_components(G),
        "communities": n_communities,
        "hub_proteins": [h[0] for h in hubs],
        "bottleneck_proteins": [b[0] for b in bottlenecks],
    }

    centrality_df = pd.DataFrame({
        "protein": list(degree_cent.keys()),
        "degree_centrality": list(degree_cent.values()),
        "betweenness": [betweenness[n] for n in degree_cent.keys()],
        "closeness": [closeness[n] for n in degree_cent.keys()],
        "community": [partition.get(n, -1) for n in degree_cent.keys()],
    }).sort_values("degree_centrality", ascending=False)

    print(f"PPI analysis: {stats['nodes']} nodes, {stats['edges']} edges, "
          f"{n_communities} communities")
    return stats, centrality_df, partition
```

## 5. PPI ネットワーク可視化

```python
def visualize_ppi_network(G, partition=None, hub_proteins=None,
                           output="figures/ppi_network.png",
                           layout="spring"):
    """
    PPI ネットワークの可視化。
    """
    import matplotlib.pyplot as plt
    import os
    os.makedirs(os.path.dirname(output), exist_ok=True)

    fig, ax = plt.subplots(figsize=(14, 14))

    if layout == "spring":
        pos = nx.spring_layout(G, k=1.5, seed=42)
    elif layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(G)

    # ノードサイズ = 次数
    node_sizes = [300 + 100 * G.degree(n) for n in G.nodes()]

    # コミュニティカラー
    if partition:
        import matplotlib.cm as cm
        n_comm = len(set(partition.values()))
        colors = [cm.Set3(partition.get(n, 0) / max(n_comm, 1)) for n in G.nodes()]
    else:
        colors = "steelblue"

    nx.draw_networkx_edges(G, pos, alpha=0.2, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes,
                           node_color=colors, alpha=0.8, ax=ax)

    # ハブタンパク質のラベル
    if hub_proteins:
        labels = {n: n for n in G.nodes() if n in hub_proteins}
    else:
        labels = {n: n for n in G.nodes() if G.degree(n) >= 5}
    nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax)

    ax.set_title(f"PPI Network ({G.number_of_nodes()} proteins, "
                 f"{G.number_of_edges()} interactions)")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {output}")
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/string_interactions.csv` | CSV |
| `results/intact_interactions.csv` | CSV |
| `results/stitch_interactions.csv` | CSV |
| `results/ppi_centrality.csv` | CSV |
| `results/ppi_network.graphml` | GraphML |
| `figures/ppi_network.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| IntAct | `intact_search_interactions` | 分子相互作用検索 |
| IntAct | `intact_get_interactions` | 相互作用データ取得 |
| IntAct | `intact_get_interactor` | 相互作用因子詳細 |
| IntAct | `intact_get_interaction_details` | 相互作用詳細 |
| IntAct | `intact_get_interaction_network` | ネットワーク取得 |
| IntAct | `intact_get_interactions_by_organism` | 生物種別相互作用 |
| IntAct | `intact_get_interactions_by_complex` | 複合体別相互作用 |
| IntAct | `intact_get_complex_details` | 複合体詳細 |
| STRING/BioGRID | `STRING_get_protein_interactions` | STRING PPI 取得 |
| STRING/BioGRID | `BioGRID_get_interactions` | BioGRID 相互作用取得 |
| STITCH | `STITCH_get_chemical_protein_interactions` | 化合物-タンパク質相互作用 |
| STITCH | `STITCH_get_interaction_partners` | 相互作用パートナー |
| STITCH | `STITCH_resolve_identifier` | 化合物 ID 解決 |
| HumanBase | `humanbase_ppi_analysis` | 組織特異的 PPI 解析 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-drug-target-profiling` | 標的タンパク質 → PPI 拡張 |
| `scientific-network-analysis` | 汎用ネットワーク解析手法 |
| `scientific-pathway-enrichment` | PPI モジュール → パスウェイ富化 |
| `scientific-protein-structure-analysis` | 構造情報 → 相互作用界面 |
| `scientific-systems-biology` | GRN ↔ PPI 統合 |

### 依存パッケージ

`networkx`, `requests`, `pandas`, `matplotlib`, `python-louvain` (community)
