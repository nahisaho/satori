---
name: scientific-bioinformatics
description: |
  バイオインフォマティクス解析パイプラインのスキル。scRNA-seq（Scanpy）、ゲノム配列解析
  （BioPython）、PPI ネットワーク解析（NetworkX）、メタボロミクスの前処理を行う際に使用。
  Scientific Skills Exp-01, 04, 07, 09 で確立したパターン。
---

# Scientific Bioinformatics Pipelines

生命科学データの解析パイプラインスキル。scRNA-seq、ゲノム配列解析、PPI ネットワーク、
メタボロミクスの 4 種のワークフローを統合的に提供する。

## When to Use

- scRNA-seq データの QC・前処理・クラスタリング・DEG 解析
- ゲノム/タンパク質配列のアラインメント・系統解析・コドン使用解析
- タンパク質相互作用ネットワークの構築・中心性解析・コミュニティ検出
- メタボロミクスの前処理（補完・変換・スケーリング）・PLS-DA

## Quick Start

## 1. scRNA-seq パイプライン（Exp-01）

### Scanpy 標準ワークフロー

```python
import scanpy as sc

def scrnaseq_pipeline(adata, min_genes=200, max_genes=5000,
                       max_pct_mito=20, n_hvg=2000, resolution=0.5):
    """
    Scanpy の標準 scRNA-seq 解析パイプライン。
    QC → 正規化 → HVG → PCA → Neighbors → Leiden → UMAP
    """
    # QC フィルタリング
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"],
                                percent_top=None, inplace=True)
    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mito].copy()
    adata = adata[adata.obs["n_genes_by_counts"] < max_genes].copy()

    # 正規化 & Log 変換
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG 選択
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]].copy()

    # PCA → Neighbors → UMAP → Leiden
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)

    return adata
```

### 細胞タイプアノテーション

```python
CELL_MARKERS = {
    "CD4+ T": ["IL7R", "CD4"],
    "CD8+ T": ["CD8A", "CD8B"],
    "B cell": ["MS4A1", "CD79A"],
    "NK": ["GNLY", "NKG7"],
    "Monocyte": ["CD14", "LYZ"],
    "DC": ["FCER1A", "CST3"],
    "Platelet": ["PPBP"],
}

def annotate_clusters(adata, marker_dict=None):
    """マーカー遺伝子ベースのクラスタアノテーション。"""
    if marker_dict is None:
        marker_dict = CELL_MARKERS
    # マーカーごとの平均発現量でクラスタに細胞型を割り当て
    cluster_annotations = {}
    for cluster in adata.obs["leiden"].unique():
        mask = adata.obs["leiden"] == cluster
        best_score = -1
        best_type = "Unknown"
        for cell_type, markers in marker_dict.items():
            available = [m for m in markers if m in adata.raw.var_names]
            if available:
                score = adata.raw[mask][:, available].X.mean()
                if score > best_score:
                    best_score = score
                    best_type = cell_type
        cluster_annotations[cluster] = best_type
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_annotations)
    return adata
```

## 2. ゲノム配列解析パイプライン（Exp-09）

```python
def compute_sequence_statistics(sequence):
    """DNA 配列の基本統計量を算出する。"""
    from collections import Counter
    seq_str = str(sequence).upper()
    counts = Counter(seq_str)
    length = len(seq_str)
    gc_content = (counts.get("G", 0) + counts.get("C", 0)) / length * 100
    return {
        "Length": length,
        "GC_Content": gc_content,
        "A": counts.get("A", 0), "T": counts.get("T", 0),
        "G": counts.get("G", 0), "C": counts.get("C", 0),
    }

def codon_usage_rscu(coding_sequence):
    """RSCU（Relative Synonymous Codon Usage）を算出する。"""
    from collections import Counter
    codons = [coding_sequence[i:i+3] for i in range(0, len(coding_sequence)-2, 3)]
    codon_counts = Counter(codons)
    # RSCU = observed / expected_if_uniform
    # 同義コドン群ごとに計算
    return codon_counts  # 簡略版、完全版は Exp-09 参照
```

## 3. PPI ネットワーク解析（Exp-04）

```python
import networkx as nx

def build_ppi_network(interactions_df, source_col, target_col,
                       weight_col=None):
    """PPI ネットワークを構築する。"""
    G = nx.Graph()
    for _, row in interactions_df.iterrows():
        if weight_col:
            G.add_edge(row[source_col], row[target_col],
                      weight=row[weight_col])
        else:
            G.add_edge(row[source_col], row[target_col])
    return G

def centrality_analysis(G):
    """4 種の中心性指標を一括計算する。"""
    centralities = pd.DataFrame({
        "Degree": dict(G.degree()),
        "Betweenness": nx.betweenness_centrality(G),
        "Closeness": nx.closeness_centrality(G),
        "Eigenvector": nx.eigenvector_centrality(G, max_iter=1000),
    })
    centralities.to_csv("results/centrality_measures.csv")
    return centralities

def community_detection(G, method="louvain"):
    """コミュニティ検出。"""
    if method == "louvain":
        import community as community_louvain
        partition = community_louvain.best_partition(G)
        nx.set_node_attributes(G, partition, "community")
        return partition
    elif method == "greedy":
        from networkx.algorithms.community import greedy_modularity_communities
        communities = list(greedy_modularity_communities(G))
        partition = {}
        for i, comm in enumerate(communities):
            for node in comm:
                partition[node] = i
        return partition
```

## 4. メタボロミクス前処理（Exp-07）

```python
from sklearn.impute import KNNImputer

def metabolomics_preprocessing(df, metabolite_cols, group_col=None):
    """
    メタボロミクスの標準前処理パイプライン。
    KNN 補完 → Log2 変換 → Pareto スケーリング
    """
    X = df[metabolite_cols].values.copy()

    # KNN 補完
    imputer = KNNImputer(n_neighbors=5)
    X_imputed = imputer.fit_transform(X)

    # Log2 変換（ゼロ値に小さな正の値を加算）
    X_log = np.log2(X_imputed + 1)

    # Pareto スケーリング
    means = X_log.mean(axis=0)
    stds = X_log.std(axis=0)
    X_pareto = (X_log - means) / np.sqrt(stds + 1e-10)

    result_df = pd.DataFrame(X_pareto, columns=metabolite_cols, index=df.index)
    if group_col:
        result_df[group_col] = df[group_col].values

    return result_df
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/centrality_measures.csv` | CSV |
| `results/community_pathway_mapping.csv` | CSV |
| `results/codon_usage.csv` | CSV |
| `figures/umap_clusters.png` | PNG |
| `figures/network_visualization.png` | PNG |

#### 参照実験

- **Exp-01**: scRNA-seq（Scanpy + PyDESeq2 + KEGG）
- **Exp-04**: PPI ネットワーク（NetworkX + Louvain）
- **Exp-07**: メタボロミクス（KNN + PLS-DA + VIP）
- **Exp-09**: ゲノム配列（BioPython + 系統樹 + RSCU）
