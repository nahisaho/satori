---
name: scientific-text-mining-nlp
description: |
  科学テキストマイニング・NLP スキル。生物医学 NER（遺伝子/疾患/薬物/化合物）・
  関係抽出（PPI / DDI / GDA）・文献ベースナレッジグラフ構築・
  エビデンス要約・トピックモデリング・引用ネットワーク解析パイプライン。
  PubTator / SemanticScholar / EuropePMC データ統合。
---

# Scientific Text Mining & NLP

科学文献に対する自然言語処理（NLP）パイプラインを提供する。
生物医学エンティティ認識、関係抽出、ナレッジグラフ構築、
トピックモデリング、自動エビデンス要約を体系的に扱う。

## When to Use

- 大量の科学文献から遺伝子・疾患・薬物名を自動抽出するとき
- タンパク質-タンパク質相互作用（PPI）等の関係を文献から抽出するとき
- 文献ベースのナレッジグラフを構築するとき
- 研究トレンドのトピックモデリングを行うとき
- 引用ネットワーク分析で影響力のある論文を同定するとき

---

## Quick Start

## 1. 生物医学 NER（Named Entity Recognition）

```python
import numpy as np
import pandas as pd

def biomedical_ner(texts, model="biobert", entity_types=None):
    """
    生物医学テキストからのエンティティ認識。

    model:
      - "biobert": BioBERT — PubMed 事前学習 BERT
      - "scispacy": SciSpaCy — 科学テキスト特化 spaCy
      - "pubtator": PubTator3 API — NCBI の NER サービス

    entity_types:
      - Gene/Protein: 遺伝子・タンパク質名
      - Disease: 疾患名（MESH / OMIM ID）
      - Chemical/Drug: 化合物・薬物名（MeSH / DrugBank ID）
      - Species: 生物種
      - Mutation: 変異（tmVar 形式）
      - Cell Line / Cell Type
    """
    if entity_types is None:
        entity_types = ["Gene", "Disease", "Chemical", "Species", "Mutation"]

    if model == "scispacy":
        import spacy
        nlp = spacy.load("en_core_sci_lg")
        from scispacy.linking import EntityLinker
        nlp.add_pipe("scispacy_linker", config={
            "resolve_abbreviations": True,
            "linker_name": "umls"
        })

        all_entities = []
        for i, text in enumerate(texts):
            doc = nlp(text)
            for ent in doc.ents:
                all_entities.append({
                    "doc_id": i,
                    "text": ent.text,
                    "label": ent.label_,
                    "start": ent.start_char,
                    "end": ent.end_char,
                    "kb_id": ent._.kb_ents[0][0] if ent._.kb_ents else None,
                    "confidence": ent._.kb_ents[0][1] if ent._.kb_ents else None,
                })

        df = pd.DataFrame(all_entities)
        print(f"  NER: {len(df)} entities from {len(texts)} documents")
        return df

    elif model == "biobert":
        from transformers import pipeline
        ner_pipeline = pipeline("ner", model="dmis-lab/biobert-large-cased-v1.1-ner",
                                 aggregation_strategy="simple")

        all_entities = []
        for i, text in enumerate(texts):
            entities = ner_pipeline(text)
            for ent in entities:
                all_entities.append({
                    "doc_id": i, "text": ent["word"],
                    "label": ent["entity_group"],
                    "score": ent["score"],
                })

        return pd.DataFrame(all_entities)
```

## 2. 関係抽出

```python
def relation_extraction(texts, relation_type="ppi", model="biobert_re"):
    """
    科学文献からの関係抽出。

    relation_type:
      - "ppi": Protein-Protein Interaction
      - "ddi": Drug-Drug Interaction
      - "gda": Gene-Disease Association
      - "chem_disease": Chemical-Disease Relation
      - "chem_gene": Chemical-Gene Interaction

    パイプライン:
      1. NER でエンティティ抽出
      2. 同一文内のエンティティペアを候補として列挙
      3. 各ペアの関係分類（BERT ベース）
      4. 信頼度フィルタリング
    """
    from transformers import pipeline

    if relation_type == "ppi":
        re_model = "dmis-lab/biobert-v1.1"  # Fine-tuned for PPI
    elif relation_type == "ddi":
        re_model = "dmis-lab/biobert-v1.1"

    classifier = pipeline("text-classification", model=re_model)

    relations = []
    for i, text in enumerate(texts):
        # エンティティペア候補をマーキング
        ner_results = biomedical_ner([text], model="scispacy")
        entities = ner_results[ner_results["doc_id"] == 0]

        # 全ペアの関係分類
        for idx_a, ent_a in entities.iterrows():
            for idx_b, ent_b in entities.iterrows():
                if idx_a < idx_b:
                    # コンテキスト付きテキスト
                    marked_text = mark_entities(text, ent_a, ent_b)
                    pred = classifier(marked_text[:512])

                    if pred[0]["score"] > 0.7:
                        relations.append({
                            "doc_id": i,
                            "entity_a": ent_a["text"],
                            "entity_b": ent_b["text"],
                            "relation": pred[0]["label"],
                            "confidence": pred[0]["score"],
                        })

    df = pd.DataFrame(relations)
    print(f"  RE: {len(df)} relations from {len(texts)} documents")
    return df


def mark_entities(text, ent_a, ent_b):
    """エンティティをマーキングしたテキストを生成。"""
    # 簡易実装: @ENTITY_A@ / @ENTITY_B@ でマーク
    marked = text.replace(ent_a["text"], f"@ENTITY_A@ {ent_a['text']} @/ENTITY_A@")
    marked = marked.replace(ent_b["text"], f"@ENTITY_B@ {ent_b['text']} @/ENTITY_B@")
    return marked
```

## 3. ナレッジグラフ構築

```python
def build_knowledge_graph(entities_df, relations_df, min_confidence=0.7):
    """
    文献ベースのナレッジグラフ構築。

    ノード: エンティティ（遺伝子、疾患、薬物、経路 etc.）
    エッジ: 関係（interacts_with, treats, causes, associated_with etc.）

    パイプライン:
      1. エンティティ正規化（UMLS CUI / MeSH 統一）
      2. 重複エンティティマージ
      3. 関係集約（頻度 + 最大信頼度）
      4. グラフ構築 + コミュニティ検出
    """
    import networkx as nx
    from collections import Counter

    # 信頼度フィルタ
    rel_filtered = relations_df[relations_df["confidence"] >= min_confidence]

    # グラフ構築
    G = nx.MultiDiGraph()

    # エンティティノード追加
    for _, ent in entities_df.iterrows():
        G.add_node(ent["text"], type=ent["label"],
                    kb_id=ent.get("kb_id", None))

    # 関係エッジ追加
    edge_counts = Counter()
    for _, rel in rel_filtered.iterrows():
        key = (rel["entity_a"], rel["entity_b"], rel["relation"])
        edge_counts[key] += 1
        G.add_edge(rel["entity_a"], rel["entity_b"],
                    relation=rel["relation"],
                    confidence=rel["confidence"],
                    frequency=edge_counts[key])

    # コミュニティ検出
    G_simple = nx.Graph(G)
    from networkx.algorithms.community import louvain_communities
    communities = louvain_communities(G_simple, resolution=1.0)

    print(f"  KG: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
          f"{len(communities)} communities")
    return G, communities
```

## 4. トピックモデリング

```python
def topic_modeling(abstracts, n_topics=10, method="bertopic"):
    """
    科学文献のトピックモデリング。

    method:
      - "bertopic": BERTopic — BERT 埋め込み + HDBSCAN + c-TF-IDF
      - "lda": LDA (Latent Dirichlet Allocation) — 確率的トピックモデル
      - "nmf": NMF (Non-negative Matrix Factorization)

    BERTopic パイプライン:
      1. BERT / SPECTER で文書埋め込み
      2. UMAP で次元削減
      3. HDBSCAN でクラスタリング
      4. c-TF-IDF でトピックワード抽出
    """
    if method == "bertopic":
        from bertopic import BERTopic
        from sentence_transformers import SentenceTransformer

        embedding_model = SentenceTransformer("allenai-specter")
        topic_model = BERTopic(embedding_model=embedding_model,
                                nr_topics=n_topics,
                                calculate_probabilities=True)

        topics, probs = topic_model.fit_transform(abstracts)

        topic_info = topic_model.get_topic_info()
        print(f"  Topics: {len(topic_info) - 1} topics from {len(abstracts)} documents")
        return topic_model, topics, probs

    elif method == "lda":
        from sklearn.decomposition import LatentDirichletAllocation
        from sklearn.feature_extraction.text import CountVectorizer

        vectorizer = CountVectorizer(max_df=0.95, min_df=2, stop_words="english")
        dtm = vectorizer.fit_transform(abstracts)

        lda = LatentDirichletAllocation(n_components=n_topics, random_state=42)
        lda.fit(dtm)

        feature_names = vectorizer.get_feature_names_out()
        topics = {}
        for i, topic_dist in enumerate(lda.components_):
            top_words = [feature_names[j] for j in topic_dist.argsort()[-10:]]
            topics[f"Topic_{i}"] = top_words

        return lda, topics
```

## 5. 引用ネットワーク分析

```python
def citation_network_analysis(papers_df, citations_df):
    """
    引用ネットワーク分析。

    指標:
      - In-degree: 被引用数 → 影響力
      - PageRank: 引用の質を加味した影響力
      - Hub/Authority (HITS): Hub=多数引用、Authority=多数被引用
      - Citation burst: 急激な被引用増加（新興トピック）
      - Bibliographic coupling: 同じ論文を引用するペア
      - Co-citation: 同時に引用されるペア
    """
    import networkx as nx

    G = nx.DiGraph()
    for _, paper in papers_df.iterrows():
        G.add_node(paper["paper_id"], title=paper["title"],
                    year=paper["year"])

    for _, cite in citations_df.iterrows():
        G.add_edge(cite["citing"], cite["cited"])

    # PageRank
    pagerank = nx.pagerank(G, alpha=0.85)

    # HITS
    hubs, authorities = nx.hits(G, max_iter=100)

    # 結果集約
    metrics_df = pd.DataFrame({
        "paper_id": list(G.nodes()),
        "in_degree": [G.in_degree(n) for n in G.nodes()],
        "out_degree": [G.out_degree(n) for n in G.nodes()],
        "pagerank": [pagerank.get(n, 0) for n in G.nodes()],
        "hub_score": [hubs.get(n, 0) for n in G.nodes()],
        "authority_score": [authorities.get(n, 0) for n in G.nodes()],
    })
    metrics_df = metrics_df.sort_values("pagerank", ascending=False)

    print(f"  Citation network: {G.number_of_nodes()} papers, "
          f"{G.number_of_edges()} citations")
    return G, metrics_df
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/ner_entities.csv` | CSV |
| `results/relations.csv` | CSV |
| `results/knowledge_graph.json` | JSON |
| `results/topic_model_info.csv` | CSV |
| `results/citation_metrics.csv` | CSV |
| `figures/kg_visualization.png` | PNG |
| `figures/topic_distribution.png` | PNG |
| `figures/citation_network.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| PubTator | `PubTator3_LiteratureSearch` | NER 付き文献検索 |
| PubTator | `PubTator3_EntityAutocomplete` | エンティティ補完 |
| PubMed | `PubMed_search_articles` | PubMed 文献検索 |
| PubMed | `PubMed_get_cited_by` | 被引用論文取得 |
| PubMed | `PubMed_get_related` | 関連論文取得 |
| SemanticScholar | `SemanticScholar_search_papers` | 学術論文検索 |
| EuropePMC | `EuropePMC_search_articles` | EuropePMC 検索 |
| EuropePMC | `EuropePMC_get_fulltext` | 全文テキスト取得 |
| EuropePMC | `EuropePMC_get_references` | 引用文献取得 |
| OpenAlex | `openalex_search_works` | OpenAlex 検索 |
| DBLP | `DBLP_search_publications` | CS 文献検索 |

### 参照スキル

| スキル | 連携内容 |
|---|---|
| [scientific-deep-research](../scientific-deep-research/SKILL.md) | 深層文献調査 |
| [scientific-citation-checker](../scientific-citation-checker/SKILL.md) | 引用検証 |
| [scientific-network-analysis](../scientific-network-analysis/SKILL.md) | ネットワーク解析 |
| [scientific-meta-analysis](../scientific-meta-analysis/SKILL.md) | 系統的文献レビュー |
| [scientific-graph-neural-networks](../scientific-graph-neural-networks/SKILL.md) | ナレッジグラフ推論 |

#### 依存パッケージ

- scispacy, spacy, transformers, bertopic, sentence-transformers, networkx
