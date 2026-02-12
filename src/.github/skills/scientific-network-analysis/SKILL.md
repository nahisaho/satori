---
name: scientific-network-analysis
description: |
  ネットワーク解析・相関ネットワーク構築のスキル。NetworkX を用いたグラフ構築、
  中心性解析、コミュニティ検出、ネットワーク可視化を行う際に使用。
  Scientific Skills Exp-04, 07 で確立したパターン。PSP パスダイアグラムにも適用。
---

# Scientific Network Analysis

NetworkX を用いたネットワーク解析パイプラインスキル。PPI ネットワーク、
相関ネットワーク、PSP パスダイアグラムなどの構築・解析・可視化を提供する。

## When to Use

- タンパク質相互作用ネットワークを構築・解析したいとき
- 相関行列からネットワークを構築したいとき
- ノードの重要性（ハブ、ボトルネック）を評価したいとき
- コミュニティ（モジュール）を検出したいとき
- PSP パスダイアグラム（因果連鎖）を可視化したいとき

## Quick Start

## 標準パイプライン

### 1. ネットワーク構築

```python
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def build_network_from_edgelist(edges_df, source_col, target_col,
                                weight_col=None, directed=False):
    """エッジリスト DataFrame からネットワークを構築する。"""
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()

    for _, row in edges_df.iterrows():
        kwargs = {}
        if weight_col and weight_col in row:
            kwargs["weight"] = row[weight_col]
        G.add_edge(row[source_col], row[target_col], **kwargs)

    return G

def build_correlation_network(corr_matrix, threshold=0.5, absolute=True):
    """
    相関行列から閾値以上のエッジを持つネットワークを構築する。
    メタボロミクスの相関ネットワーク等に使用（Exp-07）。
    """
    G = nx.Graph()
    variables = corr_matrix.columns

    for i, var1 in enumerate(variables):
        for j, var2 in enumerate(variables):
            if i < j:
                r = corr_matrix.iloc[i, j]
                if absolute:
                    if abs(r) >= threshold:
                        G.add_edge(var1, var2, weight=abs(r),
                                  correlation=r, sign="+" if r > 0 else "-")
                else:
                    if r >= threshold:
                        G.add_edge(var1, var2, weight=r)

    return G
```

### 2. 中心性解析

```python
def comprehensive_centrality(G):
    """4 種の中心性指標を一括計算する。"""
    centrality = pd.DataFrame(index=list(G.nodes()))
    centrality["Degree"] = pd.Series(dict(G.degree()))
    centrality["Betweenness"] = pd.Series(nx.betweenness_centrality(G))
    centrality["Closeness"] = pd.Series(nx.closeness_centrality(G))

    try:
        centrality["Eigenvector"] = pd.Series(
            nx.eigenvector_centrality(G, max_iter=1000)
        )
    except nx.PowerIterationFailedConvergence:
        centrality["Eigenvector"] = np.nan

    centrality = centrality.sort_values("Degree", ascending=False)
    centrality.to_csv("results/centrality_measures.csv")
    return centrality

def identify_hubs(centrality_df, top_n=10):
    """ハブノード（高次数+高媒介中心性）を同定する。"""
    # Degree と Betweenness の両方で上位のノード
    degree_top = set(centrality_df.nlargest(top_n, "Degree").index)
    between_top = set(centrality_df.nlargest(top_n, "Betweenness").index)
    hubs = degree_top & between_top
    return list(hubs), centrality_df.loc[list(hubs)]
```

### 3. コミュニティ検出

```python
def detect_communities(G, method="louvain"):
    """
    コミュニティ検出を実行する。
    method: 'louvain', 'greedy', 'label_propagation'
    """
    if method == "louvain":
        try:
            import community as community_louvain
            partition = community_louvain.best_partition(G, random_state=42)
        except ImportError:
            # フォールバック
            from networkx.algorithms.community import greedy_modularity_communities
            communities = list(greedy_modularity_communities(G))
            partition = {node: i for i, comm in enumerate(communities)
                        for node in comm}
    elif method == "greedy":
        from networkx.algorithms.community import greedy_modularity_communities
        communities = list(greedy_modularity_communities(G))
        partition = {node: i for i, comm in enumerate(communities)
                    for node in comm}
    elif method == "label_propagation":
        from networkx.algorithms.community import label_propagation_communities
        communities = list(label_propagation_communities(G))
        partition = {node: i for i, comm in enumerate(communities)
                    for node in comm}

    nx.set_node_attributes(G, partition, "community")
    modularity = nx.community.modularity(
        G, [{n for n, c in partition.items() if c == i}
            for i in set(partition.values())]
    )

    return partition, modularity
```

### 4. ネットワーク可視化

```python
def visualize_network(G, partition=None, node_size_attr="Degree",
                       title="Network", figsize=(12, 12)):
    """ネットワークを spring layout で可視化する。"""
    fig, ax = plt.subplots(figsize=figsize)
    pos = nx.spring_layout(G, seed=42, k=2/np.sqrt(len(G.nodes())))

    # ノードサイズ
    if node_size_attr == "Degree":
        sizes = np.array([G.degree(n) for n in G.nodes()])
    else:
        sizes = np.array([G.nodes[n].get(node_size_attr, 1) for n in G.nodes()])
    sizes = 100 + sizes / sizes.max() * 500

    # ノード色（コミュニティ）
    if partition:
        colors = [partition.get(n, 0) for n in G.nodes()]
        cmap = plt.cm.Set2
    else:
        colors = "steelblue"
        cmap = None

    nx.draw_networkx_edges(G, pos, alpha=0.2, ax=ax)
    nodes = nx.draw_networkx_nodes(G, pos, node_size=sizes,
                                    node_color=colors, cmap=cmap,
                                    alpha=0.8, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig("figures/network_visualization.png", dpi=300, bbox_inches="tight")
    plt.close()
```

### 5. PSP パスダイアグラム（Exp-13 独自）

```python
def psp_path_diagram(ps_corr, sp_corr, pp_corr,
                      threshold=0.3, figsize=(14, 10)):
    """
    Process → Structure → Property のパスダイアグラムを描画する。
    矢印の色: 赤=正の相関、青=負の相関、太さ=|r| に比例。
    """
    fig, ax = plt.subplots(figsize=figsize)

    # ノード配置（3 列）
    process_vars = list(ps_corr.index)
    structure_vars = list(ps_corr.columns)
    property_vars = list(sp_corr.columns)

    def place_nodes(names, x, color):
        positions = {}
        n = len(names)
        for i, name in enumerate(names):
            y = (n - 1 - i) / max(n - 1, 1)
            positions[name] = (x, y)
            ax.scatter(x, y, s=300, c=color, zorder=5, edgecolors="black")
            ax.text(x, y, name, ha="center", va="center", fontsize=7,
                   fontweight="bold", zorder=6)
        return positions

    pos_p = place_nodes(process_vars, 0, "#FFB3BA")
    pos_s = place_nodes(structure_vars, 1, "#BAE1FF")
    pos_pr = place_nodes(property_vars, 2, "#BAFFC9")

    # エッジ描画
    def draw_edges(corr_df, pos_from, pos_to):
        for var1 in corr_df.index:
            for var2 in corr_df.columns:
                r = corr_df.loc[var1, var2]
                if abs(r) >= threshold:
                    color = "red" if r > 0 else "blue"
                    width = abs(r) * 3
                    ax.annotate("", xy=pos_to[var2], xytext=pos_from[var1],
                              arrowprops=dict(arrowstyle="->", color=color,
                                            lw=width, alpha=0.6))

    draw_edges(ps_corr, pos_p, pos_s)
    draw_edges(sp_corr, pos_s, pos_pr)

    # ラベル
    ax.text(0, -0.1, "Process", ha="center", fontsize=12,
           fontweight="bold", color="#FF6B6B")
    ax.text(1, -0.1, "Structure", ha="center", fontsize=12,
           fontweight="bold", color="#4ECDC4")
    ax.text(2, -0.1, "Property", ha="center", fontsize=12,
           fontweight="bold", color="#45B7D1")

    ax.set_xlim(-0.3, 2.3)
    ax.set_ylim(-0.2, 1.1)
    ax.axis("off")
    ax.set_title("PSP Linkage Path Diagram", fontweight="bold", fontsize=14)
    plt.tight_layout()
    plt.savefig("figures/psp_path_diagram.png", dpi=300, bbox_inches="tight")
    plt.close()
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/centrality_measures.csv` | CSV |
| `results/edge_list.csv` | CSV |
| `results/node_attributes.csv` | CSV |
| `figures/network_visualization.png` | PNG |
| `figures/psp_path_diagram.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| STRING | `STRING_get_protein_interactions` | タンパク質相互作用ネットワーク |
| IntAct | `intact_search_interactions` | 分子相互作用検索 |
| IntAct | `intact_get_interaction_network` | PPI ネットワーク取得 |
| Reactome | `Reactome_map_uniprot_to_pathways` | パスウェイマッピング |
| KEGG | `kegg_get_pathway_info` | パスウェイ情報 |
| GO | `GO_get_annotations_for_gene` | GO アノテーション |

#### 参照実験

- **Exp-04**: PPI ネットワーク（71 タンパク質、4 種中心性、Louvain コミュニティ）
- **Exp-07**: Spearman 相関ネットワーク（メタボロミクス）
- **Exp-13**: PSP パスダイアグラム（Process→Structure→Property 因果連鎖）
