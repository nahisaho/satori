---
name: scientific-network-visualization
description: |
  ネットワーク解析・可視化スキル。NetworkX グラフ構築・
  コミュニティ検出 (Louvain/Leiden)・中心性指標・
  PyVis インタラクティブ・ネットワーク統計量・動的ネットワーク。
tu_tools:
  - key: string
    name: STRING
    description: ネットワーク可視化・相互作用データ
---

# Scientific Network Visualization

ネットワーク/グラフデータの解析・コミュニティ検出・
インタラクティブ可視化パイプラインを提供する。

## When to Use

- 関連性・相互作用のネットワーク構造を解析するとき
- コミュニティ (クラスタ) を検出するとき
- 中心性指標でハブ・ブリッジノードを特定するとき
- PyVis でインタラクティブなネットワーク図を作成するとき
- 相関行列からネットワークを構築するとき
- 時間発展するネットワークを解析するとき

> **Note**: タンパク質 PPI ネットワークは `scientific-protein-interaction-network` を参照。

---

## Quick Start

## 1. ネットワーク構築・基本統計

```python
import numpy as np
import pandas as pd
import networkx as nx


def build_network_from_edgelist(edges_df, source_col, target_col,
                                weight_col=None, directed=False):
    """
    エッジリストからネットワーク構築 + 基本統計。

    Parameters:
        edges_df: pd.DataFrame — エッジリスト
        source_col: str — ソースノードカラム
        target_col: str — ターゲットノードカラム
        weight_col: str | None — 重みカラム
        directed: bool — 有向グラフ
    """
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()

    for _, row in edges_df.iterrows():
        kwargs = {}
        if weight_col and pd.notna(row[weight_col]):
            kwargs["weight"] = row[weight_col]
        G.add_edge(row[source_col], row[target_col], **kwargs)

    stats = {
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
        "density": nx.density(G),
        "is_connected": nx.is_connected(G) if not directed else nx.is_weakly_connected(G),
        "n_components": nx.number_connected_components(G) if not directed
                        else nx.number_weakly_connected_components(G),
        "avg_degree": np.mean([d for _, d in G.degree()]),
        "avg_clustering": nx.average_clustering(G) if not directed else None,
    }

    if stats["is_connected"] and not directed:
        stats["avg_path_length"] = nx.average_shortest_path_length(G)
        stats["diameter"] = nx.diameter(G)

    print(f"Network: {stats['n_nodes']} nodes, {stats['n_edges']} edges, "
          f"density={stats['density']:.4f}")
    return G, stats


def build_network_from_correlation(df, threshold=0.5,
                                    method="pearson"):
    """
    相関行列からネットワーク構築。

    Parameters:
        df: pd.DataFrame — 数値データ
        threshold: float — 相関閾値 (|r| ≥ threshold でエッジ)
        method: str — "pearson" / "spearman"
    """
    corr = df.corr(method=method)
    G = nx.Graph()

    for i, col_i in enumerate(corr.columns):
        G.add_node(col_i)
        for j, col_j in enumerate(corr.columns):
            if i < j and abs(corr.iloc[i, j]) >= threshold:
                G.add_edge(col_i, col_j,
                           weight=abs(corr.iloc[i, j]),
                           correlation=corr.iloc[i, j])

    print(f"Correlation Network (|r|≥{threshold}): "
          f"{G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G
```

## 2. コミュニティ検出

```python
def detect_communities(G, method="louvain", resolution=1.0):
    """
    コミュニティ検出。

    Parameters:
        G: nx.Graph — ネットワーク
        method: str — "louvain" / "leiden" / "label_propagation" / "girvan_newman"
        resolution: float — 解像度パラメータ (Louvain/Leiden)
    """
    import matplotlib.pyplot as plt

    if method == "louvain":
        communities = nx.community.louvain_communities(
            G, resolution=resolution, seed=42)
    elif method == "leiden":
        try:
            import leidenalg
            import igraph as ig
            ig_graph = ig.Graph.from_networkx(G)
            partition = leidenalg.find_partition(
                ig_graph, leidenalg.RBConfigurationVertexPartition,
                resolution_parameter=resolution, seed=42)
            communities = [set(ig_graph.vs[c]["_nx_name"] for c in comm)
                           for comm in partition]
        except ImportError:
            communities = nx.community.louvain_communities(
                G, resolution=resolution, seed=42)
    elif method == "label_propagation":
        communities = list(nx.community.label_propagation_communities(G))
    elif method == "girvan_newman":
        comp = nx.community.girvan_newman(G)
        communities = next(comp)  # 最初の分割

    # ノードにコミュニティ ID 割当
    node_community = {}
    for i, comm in enumerate(communities):
        for node in comm:
            node_community[node] = i
    nx.set_node_attributes(G, node_community, "community")

    # モジュラリティ
    modularity = nx.community.modularity(G, communities)

    # 可視化
    fig, ax = plt.subplots(figsize=(12, 10))
    pos = nx.spring_layout(G, k=1/np.sqrt(G.number_of_nodes()), seed=42)
    colors = [node_community.get(n, 0) for n in G.nodes()]

    nx.draw_networkx(G, pos, ax=ax, node_color=colors,
                     cmap=plt.cm.Set3, node_size=100,
                     font_size=6, edge_color="gray", alpha=0.7,
                     with_labels=G.number_of_nodes() < 100)
    ax.set_title(f"Communities ({method}): {len(communities)} clusters, "
                 f"Q={modularity:.4f}")
    plt.tight_layout()

    path = "network_communities.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"Communities ({method}): {len(communities)} clusters, "
          f"modularity={modularity:.4f}")
    return {"communities": communities, "modularity": modularity,
            "node_community": node_community, "fig": path}
```

## 3. 中心性指標

```python
def centrality_analysis(G, top_n=20):
    """
    多面的中心性解析。

    Parameters:
        G: nx.Graph — ネットワーク
        top_n: int — 上位ノード数
    """
    centralities = {
        "degree": nx.degree_centrality(G),
        "betweenness": nx.betweenness_centrality(G),
        "closeness": nx.closeness_centrality(G),
        "eigenvector": nx.eigenvector_centrality(G, max_iter=1000),
        "pagerank": nx.pagerank(G)
    }

    # DataFrame 化
    cent_df = pd.DataFrame(centralities)
    cent_df.index.name = "node"

    # ランキング
    rankings = {}
    for metric in centralities:
        top = cent_df[metric].nlargest(top_n)
        rankings[metric] = top.index.tolist()

    # ハブスコア (複数指標の統合)
    for metric in centralities:
        cent_df[f"{metric}_rank"] = cent_df[metric].rank(ascending=False)
    cent_df["hub_score"] = cent_df[[f"{m}_rank" for m in centralities]].mean(axis=1)
    cent_df = cent_df.sort_values("hub_score")

    print(f"Centrality: {len(G.nodes())} nodes analyzed")
    print(f"  Top hubs: {cent_df.head(5).index.tolist()}")
    return {"centrality_df": cent_df, "rankings": rankings}
```

## 4. PyVis インタラクティブ可視化

```python
def interactive_network(G, output="network_interactive.html",
                        height="700px", width="100%"):
    """
    PyVis インタラクティブネットワーク図。

    Parameters:
        G: nx.Graph — ネットワーク
        output: str — 出力 HTML パス
        height: str — 高さ
        width: str — 幅
    """
    from pyvis.network import Network

    nt = Network(height=height, width=width, notebook=False,
                 bgcolor="#ffffff", font_color="black")

    # コミュニティカラーリング
    community_map = nx.get_node_attributes(G, "community")
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
              "#ff7f00", "#ffff33", "#a65628", "#f781bf"]

    for node in G.nodes():
        comm = community_map.get(node, 0)
        degree = G.degree(node)
        nt.add_node(str(node), label=str(node),
                    color=colors[comm % len(colors)],
                    size=max(5, min(degree * 3, 50)),
                    title=f"{node}\nDegree: {degree}\nCommunity: {comm}")

    for u, v, data in G.edges(data=True):
        weight = data.get("weight", 1)
        nt.add_edge(str(u), str(v), value=weight)

    nt.toggle_physics(True)
    nt.show_buttons(filter_=["physics"])
    nt.save_graph(output)

    print(f"Interactive Network → {output} "
          f"({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")
    return output
```

---

## パイプライン統合

```
eda-correlation → network-visualization → advanced-visualization
  (相関解析)        (ネットワーク解析)       (高度可視化)
       │                   │                      ↓
 graph-neural-networks ───┘          interactive-dashboard
  (GNN)                                (ダッシュボード)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `network_communities.png` | コミュニティ構造 | → presentation |
| `centrality_analysis.csv` | 中心性指標 | → feature-importance |
| `network_interactive.html` | PyVis 図 | → dashboard |

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|--------|
| `string` | STRING | ネットワーク可視化・相互作用データ |
