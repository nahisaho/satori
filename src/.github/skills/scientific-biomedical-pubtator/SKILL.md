---
name: scientific-biomedical-pubtator
description: |
  バイオメディカルテキストマイニングスキル。PubTator3 API による
  遺伝子・疾患・化合物・変異・種のエンティティ認識、関係抽出、
  バイオメディカル文献アノテーション自動化パイプライン。
---

# Scientific Biomedical PubTator

PubTator3 API を活用したバイオメディカル文献エンティティ認識・
関係抽出パイプラインを提供する。

## When to Use

- PubMed 論文から遺伝子・疾患・化合物のエンティティを自動抽出するとき
- バイオメディカル NER (Named Entity Recognition) を実行するとき
- 遺伝子-疾患・薬物-標的の関係を文献から抽出するとき
- 大規模文献コーパスのバイオアノテーションを行うとき
- テキストマイニング結果を知識グラフに統合するとき

---

## Quick Start

## 1. PubTator3 エンティティアノテーション

```python
import requests
import pandas as pd
import json
import time

PUBTATOR_API = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"


def annotate_pmids(pmids, concepts=None):
    """
    PubTator3 — PMID リストのバイオメディカルエンティティアノテーション。

    Parameters:
        pmids: list — PMID リスト (e.g., [12345678, 23456789])
        concepts: list — エンティティタイプ
            "gene", "disease", "chemical", "mutation", "species", "cellline"

    ToolUniverse:
        PubTator_annotate(pmids=pmids, concepts=concepts)
        PubTator_search(query=query)
    """
    if concepts is None:
        concepts = ["gene", "disease", "chemical", "mutation", "species"]

    pmid_str = ",".join(str(p) for p in pmids)
    params = {
        "pmids": pmid_str,
        "concepts": ",".join(concepts),
        "format": "biocjson",
    }

    resp = requests.get(f"{PUBTATOR_API}/publications/export/biocjson", params=params)
    resp.raise_for_status()
    data = resp.json()

    # Parse annotations
    all_annotations = []
    for doc in data.get("PubTator3", []) if isinstance(data, dict) else [data]:
        pmid = doc.get("pmid", doc.get("id", ""))
        for passage in doc.get("passages", []):
            for annotation in passage.get("annotations", []):
                infons = annotation.get("infons", {})
                all_annotations.append({
                    "pmid": pmid,
                    "text": annotation.get("text", ""),
                    "type": infons.get("type", ""),
                    "identifier": infons.get("identifier", ""),
                    "offset": annotation.get("locations", [{}])[0].get("offset", ""),
                    "length": annotation.get("locations", [{}])[0].get("length", ""),
                    "passage_type": passage.get("infons", {}).get("type", ""),
                })

    df = pd.DataFrame(all_annotations)
    type_counts = df["type"].value_counts().to_dict() if not df.empty else {}
    print(f"PubTator annotation: {len(pmids)} PMIDs → "
          f"{len(df)} entities {type_counts}")
    return df
```

## 2. PubTator3 テキスト検索

```python
def search_pubtator(query, max_results=100):
    """
    PubTator3 テキスト検索 — バイオメディカルエンティティ付き論文検索。

    Parameters:
        query: str — 検索クエリ (遺伝子名、疾患名、化合物名)
        max_results: int — 最大取得数
    """
    params = {
        "text": query,
        "sort": "score",
        "page_size": min(max_results, 100),
    }
    resp = requests.get(f"{PUBTATOR_API}/search/", params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for hit in data.get("results", []):
        results.append({
            "pmid": hit.get("pmid", ""),
            "title": hit.get("title", ""),
            "journal": hit.get("journal", ""),
            "year": hit.get("year", ""),
            "score": hit.get("score", 0),
            "genes": hit.get("genes", []),
            "diseases": hit.get("diseases", []),
            "chemicals": hit.get("chemicals", []),
            "mutations": hit.get("mutations", []),
        })

    df = pd.DataFrame(results)
    total = data.get("count", 0)
    print(f"PubTator search '{query}': {total} total, {len(df)} returned")
    return df
```

## 3. エンティティ関係抽出

```python
def extract_entity_relations(pmids, relation_types=None):
    """
    PubTator3 — エンティティ間関係 (gene-disease, drug-target 等) 抽出。

    Parameters:
        pmids: list — PMID リスト
        relation_types: list — 関係タイプフィルタ
            "GDA" (gene-disease), "CDA" (chemical-disease),
            "CGA" (chemical-gene), "PPI" (protein-protein)
    """
    # Get annotations with relations
    df_annotations = annotate_pmids(pmids)

    # Extract co-occurrences within same passage
    relations = []
    for pmid in df_annotations["pmid"].unique():
        pmid_df = df_annotations[df_annotations["pmid"] == pmid]

        # Gene-Disease relations
        genes = pmid_df[pmid_df["type"] == "Gene"]
        diseases = pmid_df[pmid_df["type"] == "Disease"]
        chemicals = pmid_df[pmid_df["type"] == "Chemical"]

        if not relation_types or "GDA" in relation_types:
            for _, gene in genes.iterrows():
                for _, disease in diseases.iterrows():
                    relations.append({
                        "pmid": pmid,
                        "relation_type": "GDA",
                        "entity1_type": "Gene",
                        "entity1_text": gene["text"],
                        "entity1_id": gene["identifier"],
                        "entity2_type": "Disease",
                        "entity2_text": disease["text"],
                        "entity2_id": disease["identifier"],
                    })

        if not relation_types or "CGA" in relation_types:
            for _, chem in chemicals.iterrows():
                for _, gene in genes.iterrows():
                    relations.append({
                        "pmid": pmid,
                        "relation_type": "CGA",
                        "entity1_type": "Chemical",
                        "entity1_text": chem["text"],
                        "entity1_id": chem["identifier"],
                        "entity2_type": "Gene",
                        "entity2_text": gene["text"],
                        "entity2_id": gene["identifier"],
                    })

        if not relation_types or "CDA" in relation_types:
            for _, chem in chemicals.iterrows():
                for _, disease in diseases.iterrows():
                    relations.append({
                        "pmid": pmid,
                        "relation_type": "CDA",
                        "entity1_type": "Chemical",
                        "entity1_text": chem["text"],
                        "entity1_id": chem["identifier"],
                        "entity2_type": "Disease",
                        "entity2_text": disease["text"],
                        "entity2_id": disease["identifier"],
                    })

    rel_df = pd.DataFrame(relations)
    rel_counts = rel_df["relation_type"].value_counts().to_dict() if not rel_df.empty else {}
    print(f"Entity relations: {len(rel_df)} total {rel_counts}")
    return rel_df
```

## 4. バイオアノテーション集計ダッシュボード

```python
def annotation_summary_dashboard(pmids, output_prefix="pubtator"):
    """
    PubTator アノテーション集計・可視化。

    Parameters:
        pmids: list — PMID リスト
        output_prefix: str — 出力ファイルプレフィックス
    """
    import matplotlib.pyplot as plt

    # Get annotations
    df = annotate_pmids(pmids)
    if df.empty:
        print("No annotations found")
        return {}

    # Entity type distribution
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # 1. Entity type counts
    type_counts = df["type"].value_counts()
    type_counts.plot(kind="bar", ax=axes[0], color="#2196F3")
    axes[0].set_title("Entity Type Distribution")
    axes[0].set_ylabel("Count")

    # 2. Top entities per type
    for entity_type in ["Gene", "Disease", "Chemical"]:
        sub = df[df["type"] == entity_type]
        top = sub["text"].value_counts().head(10)
        if not top.empty:
            print(f"\nTop {entity_type}s: {top.to_dict()}")

    # 3. Articles per entity count
    per_article = df.groupby("pmid")["type"].count()
    per_article.hist(ax=axes[1], bins=20, color="#4CAF50")
    axes[1].set_title("Entities per Article")
    axes[1].set_xlabel("Number of entities")

    # Entity type per article
    pivot = df.groupby(["pmid", "type"]).size().unstack(fill_value=0)
    pivot.plot(kind="box", ax=axes[2])
    axes[2].set_title("Entity Types per Article")

    plt.tight_layout()
    fig_path = f"figures/{output_prefix}_dashboard.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()

    # Save results
    df.to_csv(f"results/{output_prefix}_annotations.csv", index=False)

    summary = {
        "total_pmids": df["pmid"].nunique(),
        "total_annotations": len(df),
        "entity_types": type_counts.to_dict(),
        "unique_entities": df.groupby("type")["text"].nunique().to_dict(),
    }
    print(f"\nSummary: {summary}")
    return summary
```

## 5. 知識グラフ構築用エンティティネットワーク

```python
def build_entity_network(pmids, min_cooccurrence=2):
    """
    PubTator エンティティ共起ネットワーク構築。

    Parameters:
        pmids: list — PMID リスト
        min_cooccurrence: int — 最小共起回数
    """
    import networkx as nx
    from collections import Counter

    rel_df = extract_entity_relations(pmids)
    if rel_df.empty:
        return nx.Graph()

    # Count co-occurrences
    edge_counter = Counter()
    for _, row in rel_df.iterrows():
        key = tuple(sorted([
            f"{row['entity1_type']}:{row['entity1_text']}",
            f"{row['entity2_type']}:{row['entity2_text']}",
        ]))
        edge_counter[key] += 1

    # Build network
    G = nx.Graph()
    for (node1, node2), count in edge_counter.items():
        if count >= min_cooccurrence:
            G.add_edge(node1, node2, weight=count)

    print(f"Entity network: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges "
          f"(min cooccurrence = {min_cooccurrence})")
    return G
```

---

## 利用可能ツール

| ToolUniverse カテゴリ | 主なツール |
|---|---|
| `pubtator` | `PubTator_annotate`, `PubTator_search` |

## パイプライン出力

| 出力ファイル | 説明 | 連携先スキル |
|---|---|---|
| `results/pubtator_annotations.csv` | エンティティアノテーション | → text-mining-nlp, knowledge-graph |
| `results/entity_relations.csv` | エンティティ間関係 | → network-analysis, disease-research |
| `results/entity_network.graphml` | エンティティ共起ネットワーク | → graph-neural-networks |
| `figures/pubtator_dashboard.png` | アノテーション集計 | → publication-figures |

## パイプライン統合

```
literature-search ──→ biomedical-pubtator ──→ text-mining-nlp
  (PubMed/OpenAlex)   (PubTator NER)         (KG 構築)
                            │
                            ├──→ disease-research (GDA 関係)
                            ├──→ drug-target-profiling (CGA 関係)
                            └──→ preprint-archive (プレプリント NER)
```
