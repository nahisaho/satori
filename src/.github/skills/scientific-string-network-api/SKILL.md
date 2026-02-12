---
name: scientific-string-network-api
description: |
  STRING/BioGRID/STITCH ネットワーク解析スキル。STRING タンパク質相互作用
  ネットワーク直接 API、BioGRID 実験的 PPI、STITCH 化学-タンパク質ネットワーク、
  ネットワークトポロジー解析・コミュニティ検出・機能濃縮統合パイプライン。
tu_tools:
  - key: ppi
    name: STRING/BioGRID PPI
    description: タンパク質・化学物質相互作用ネットワーク
---

# Scientific STRING Network API

STRING v12 / BioGRID / STITCH API を活用した PPI・化合物-タンパク質
ネットワーク解析パイプラインを提供する。既存の protein-interaction-network
スキル (IntAct/HumanBase) を補完し、STRING 直接 API ベースの高度な
ネットワーク分析を統合。

## When to Use

- STRING API でタンパク質相互作用ネットワークを直接構築するとき
- BioGRID から実験的エビデンスベースの PPI を取得するとき
- STITCH で化合物-タンパク質間ネットワークを検索するとき
- ネットワークトポロジー指標 (次数分布・媒介中心性) を計算するとき
- PPI ネットワーク上でコミュニティ検出を行うとき
- 機能濃縮解析 (STRING enrichment) をネットワーク上で実行するとき

---

## Quick Start

## 1. STRING PPI ネットワーク取得

```python
import requests
import pandas as pd
import networkx as nx

STRING_API = "https://string-db.org/api"
OUTPUT_FORMAT = "json"


def get_string_network(proteins, species=9606, score_threshold=400,
                       network_type="functional", limit=50):
    """
    STRING PPI ネットワーク取得。

    Parameters:
        proteins: list — タンパク質名リスト (例: ["TP53", "MDM2", "BRCA1"])
        species: int — NCBI Taxonomy ID (9606=human)
        score_threshold: int — 信頼スコア閾値 (0-1000)
        network_type: str — "functional" or "physical"
        limit: int — interactor 最大数

    ToolUniverse:
        STRING_get_protein_interactions(
            protein_ids=proteins, species=species,
            confidence_score=score_threshold/1000,
            network_type=network_type, limit=limit
        )
    """
    url = f"{STRING_API}/{OUTPUT_FORMAT}/network"
    params = {
        "identifiers": "\r".join(proteins),
        "species": species,
        "required_score": score_threshold,
        "network_type": network_type,
        "limit": limit,
    }

    resp = requests.post(url, data=params)
    resp.raise_for_status()
    interactions = resp.json()

    rows = []
    for i in interactions:
        rows.append({
            "protein_a": i.get("preferredName_A"),
            "protein_b": i.get("preferredName_B"),
            "combined_score": i.get("score"),
            "nscore": i.get("nscore"),
            "fscore": i.get("fscore"),
            "pscore": i.get("pscore"),
            "ascore": i.get("ascore"),
            "escore": i.get("escore"),
            "dscore": i.get("dscore"),
            "tscore": i.get("tscore"),
        })

    df = pd.DataFrame(rows)
    print(f"STRING network: {len(df)} interactions "
          f"(score ≥ {score_threshold/1000})")
    return df
```

## 2. BioGRID 実験的 PPI 取得

```python
def get_biogrid_interactions(genes, organism=9606, evidence_type=None,
                             api_key="YOUR_KEY", limit=500):
    """
    BioGRID 実験的 PPI データ取得。

    Parameters:
        genes: list — 遺伝子名リスト
        organism: int — NCBI Taxonomy ID
        evidence_type: str — "physical" or "genetic"
        api_key: str — BioGRID API key (https://webservice.thebiogrid.org)
        limit: int — 最大取得数

    ToolUniverse:
        BioGRID_get_interactions(
            gene_names=genes, organism=organism,
            interaction_type=evidence_type, limit=limit
        )
    """
    url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": api_key,
        "geneList": "|".join(genes),
        "organism": organism,
        "format": "json",
        "max": limit,
        "searchNames": "true",
        "includeInteractors": "true",
    }
    if evidence_type:
        params["interSpeciesExcluded"] = "true"

    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    rows = []
    for _, interaction in data.items():
        rows.append({
            "gene_a": interaction.get("OFFICIAL_SYMBOL_A"),
            "gene_b": interaction.get("OFFICIAL_SYMBOL_B"),
            "experimental_system": interaction.get("EXPERIMENTAL_SYSTEM"),
            "throughput": interaction.get("THROUGHPUT"),
            "pubmed_id": interaction.get("PUBMED_ID"),
            "source_db": "BioGRID",
        })

    df = pd.DataFrame(rows)
    print(f"BioGRID: {len(df)} interactions for {genes}")
    return df
```

## 3. STITCH 化合物-タンパク質ネットワーク

```python
def get_stitch_interactions(identifiers, species=9606, score=400, limit=20):
    """
    STITCH 化合物-タンパク質相互作用取得。

    Parameters:
        identifiers: list — CID (化合物) または遺伝子名リスト
        species: int — NCBI Taxonomy ID
        score: int — 信頼スコア閾値
        limit: int — 最大結果数

    ToolUniverse:
        STITCH_get_chemical_protein_interactions(
            identifiers=identifiers, species=species,
            required_score=score, limit=limit
        )
        STITCH_get_interaction_partners(identifiers=identifiers)
        STITCH_resolve_identifier(identifiers=identifiers)
    """
    url = f"https://stitch.embl.de/api/{OUTPUT_FORMAT}/interactionsList"
    params = {
        "identifiers": "\r".join(identifiers),
        "species": species,
        "required_score": score,
        "limit": limit,
    }

    resp = requests.post(url, data=params)
    resp.raise_for_status()
    interactions = resp.json()

    rows = []
    for i in interactions:
        rows.append({
            "interactor_a": i.get("preferredName_A", i.get("stringId_A")),
            "interactor_b": i.get("preferredName_B", i.get("stringId_B")),
            "combined_score": i.get("score"),
            "is_chemical": "CID" in str(i.get("stringId_A", ""))
                or "CID" in str(i.get("stringId_B", "")),
        })

    df = pd.DataFrame(rows)
    print(f"STITCH: {len(df)} chemical-protein interactions")
    return df
```

## 4. ネットワーク構築 & トポロジー解析

```python
def build_network(interaction_df, source_col="protein_a", target_col="protein_b",
                  weight_col="combined_score"):
    """
    NetworkX グラフ構築 & トポロジー解析。

    Parameters:
        interaction_df: DataFrame — 相互作用データ
        source_col, target_col: str — ノードカラム名
        weight_col: str — エッジ重みカラム名
    """
    G = nx.Graph()
    for _, row in interaction_df.iterrows():
        G.add_edge(
            row[source_col], row[target_col],
            weight=row.get(weight_col, 1.0),
        )

    # トポロジー指標
    degree = dict(G.degree())
    betweenness = nx.betweenness_centrality(G)
    closeness = nx.closeness_centrality(G)
    clustering = nx.clustering(G)

    metrics = pd.DataFrame({
        "node": list(degree.keys()),
        "degree": list(degree.values()),
        "betweenness": [betweenness[n] for n in degree],
        "closeness": [closeness[n] for n in degree],
        "clustering": [clustering[n] for n in degree],
    }).sort_values("betweenness", ascending=False)

    print(f"Network: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges, "
          f"density={nx.density(G):.4f}")
    return G, metrics
```

## 5. コミュニティ検出

```python
from networkx.algorithms.community import greedy_modularity_communities


def detect_communities(G, resolution=1.0):
    """
    ネットワーク上のコミュニティ (モジュール) 検出。

    Parameters:
        G: nx.Graph — ネットワークグラフ
        resolution: float — 解像度パラメータ
    """
    communities = list(greedy_modularity_communities(G, resolution=resolution))
    modularity = nx.algorithms.community.modularity(G, communities)

    comm_data = []
    for i, comm in enumerate(communities):
        for node in comm:
            comm_data.append({"node": node, "community": i})

    df = pd.DataFrame(comm_data)
    print(f"Communities: {len(communities)} detected, "
          f"modularity={modularity:.4f}")
    return df, modularity
```

## 6. STRING 機能濃縮解析

```python
def string_enrichment(proteins, species=9606):
    """
    STRING API 機能濃縮解析 (GO/KEGG/Reactome/InterPro)。

    Parameters:
        proteins: list — タンパク質名リスト
        species: int — NCBI Taxonomy ID
    """
    url = f"{STRING_API}/{OUTPUT_FORMAT}/enrichment"
    params = {
        "identifiers": "\r".join(proteins),
        "species": species,
    }

    resp = requests.post(url, data=params)
    resp.raise_for_status()
    enrichment = resp.json()

    rows = []
    for e in enrichment:
        rows.append({
            "category": e.get("category"),
            "term": e.get("term"),
            "description": e.get("description"),
            "p_value": e.get("p_value"),
            "fdr": e.get("fdr"),
            "number_of_genes": e.get("number_of_genes"),
            "input_genes": e.get("inputGenes", ""),
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("fdr")
    print(f"Enrichment: {len(df)} terms, "
          f"{df[df['fdr'] < 0.05].shape[0]} significant (FDR<0.05)")
    return df
```

## 7. 統合 PPI 解析パイプライン

```python
def integrated_ppi_pipeline(genes, species=9606, score=700):
    """
    STRING + BioGRID + STITCH 統合 PPI パイプライン。

    Pipeline:
        STRING network → BioGRID validation → topology → communities →
        enrichment
    """
    # STRING ネットワーク
    string_df = get_string_network(genes, species, score)

    # ネットワーク構築 & トポロジー
    G, metrics = build_network(string_df)

    # コミュニティ検出
    comm_df, modularity = detect_communities(G)

    # STRING 濃縮解析
    all_nodes = list(G.nodes())
    enrichment = string_enrichment(all_nodes[:500], species)

    result = {
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
        "density": round(nx.density(G), 4),
        "n_communities": comm_df["community"].nunique(),
        "modularity": round(modularity, 4),
        "hub_genes": metrics.head(10)["node"].tolist(),
        "n_enriched_terms": len(enrichment[enrichment["fdr"] < 0.05])
            if not enrichment.empty else 0,
    }

    print(f"\n=== Integrated PPI Pipeline ===")
    print(f"Nodes: {result['n_nodes']}, Edges: {result['n_edges']}")
    print(f"Hub genes: {', '.join(result['hub_genes'][:5])}")
    return result
```

---

## パイプライン統合

```
drug-target-profiling → string-network-api → pathway-enrichment
  (候補ターゲット)       (STRING PPI 構築)     (GO/KEGG 濃縮)
        │                       │                    ↓
protein-interaction ───┘        │           ontology-enrichment
  (IntAct/HumanBase)            ↓             (EFO/Enrichr)
                       network-analysis
                       (既存スキル補完)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/string_network.csv` | STRING PPI ネットワーク | → network-analysis |
| `results/ppi_topology.csv` | トポロジー指標 | → drug-target-profiling |
| `results/ppi_communities.csv` | コミュニティ割当 | → pathway-enrichment |
| `results/string_enrichment.csv` | 機能濃縮結果 | → ontology-enrichment |

## 利用可能ツール (ToolUniverse SMCP)

| ツール名 | 用途 |
|---------|------|
| `STRING_get_protein_interactions` | STRING PPI 取得 |
| `BioGRID_get_interactions` | BioGRID 実験的 PPI |
| `STITCH_get_chemical_protein_interactions` | STITCH 化合物-タンパク質 |
| `STITCH_get_interaction_partners` | STITCH 相互作用パートナー |
| `STITCH_resolve_identifier` | STITCH ID 解決 |
