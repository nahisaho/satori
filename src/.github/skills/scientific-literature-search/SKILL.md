---
name: scientific-literature-search
description: |
  学術文献検索・取得スキル。PubMed E-utilities、Semantic Scholar、
  OpenAlex、EuropePMC、CrossRef の 5 大学術データベース API を統合した
  文献検索パイプライン。MeSH 構造化検索、引用ネットワーク分析、
  著者/機関メトリクス、全文取得、PICO ベース検索戦略対応。
  29 の ToolUniverse SMCP ツールと連携。
---

# Scientific Literature Search

PubMed / Semantic Scholar / OpenAlex / EuropePMC / CrossRef の
5 大学術 DB を統合した文献検索・メタデータ取得パイプラインを提供する。

## When to Use

- PubMed で MeSH 用語を用いた構造化検索が必要なとき
- Semantic Scholar でセマンティック類似論文を発見するとき
- OpenAlex で著者・機関・ジャーナルメトリクスを分析するとき
- 特定論文の引用/被引用ネットワークを構築するとき
- 系統的レビューのためのマルチ DB 横断検索が必要なとき

---

## Quick Start

## 1. PubMed E-utilities 検索

```python
import requests
import xml.etree.ElementTree as ET
import pandas as pd
import time


def pubmed_search(query, max_results=100, sort="relevance",
                  date_from=None, date_to=None, rettype="xml"):
    """
    PubMed E-utilities による構造化検索。

    Parameters:
        query: str — PubMed クエリ (MeSH 構文対応)
        max_results: int — 最大取得件数
        sort: "relevance" or "date"
        date_from/to: "YYYY/MM/DD" フォーマット
    """
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Step 1: esearch — PMID 取得
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "sort": sort,
        "retmode": "json",
    }
    if date_from:
        params["mindate"] = date_from
        params["datetype"] = "pdat"
    if date_to:
        params["maxdate"] = date_to

    resp = requests.get(f"{base}/esearch.fcgi", params=params)
    data = resp.json()
    pmids = data["esearchresult"]["idlist"]
    total = int(data["esearchresult"]["count"])
    print(f"PubMed search: {total} total results, fetching {len(pmids)}")

    if not pmids:
        return pd.DataFrame()

    # Step 2: efetch — 詳細取得
    time.sleep(0.4)
    fetch_params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "rettype": "xml",
        "retmode": "xml",
    }
    resp = requests.get(f"{base}/efetch.fcgi", params=fetch_params)
    root = ET.fromstring(resp.text)

    articles = []
    for article in root.findall(".//PubmedArticle"):
        medline = article.find(".//MedlineCitation")
        pmid = medline.findtext("PMID", "")
        title = medline.findtext(".//ArticleTitle", "")
        abstract = medline.findtext(".//AbstractText", "")
        journal = medline.findtext(".//Journal/Title", "")
        year = medline.findtext(".//PubDate/Year", "")

        authors = []
        for author in medline.findall(".//Author"):
            last = author.findtext("LastName", "")
            fore = author.findtext("ForeName", "")
            if last:
                authors.append(f"{last} {fore}")

        mesh_terms = [m.findtext("DescriptorName", "")
                      for m in medline.findall(".//MeshHeading")]

        articles.append({
            "pmid": pmid,
            "title": title,
            "abstract": abstract[:500],
            "journal": journal,
            "year": year,
            "authors": "; ".join(authors[:5]),
            "mesh_terms": "; ".join(mesh_terms[:10]),
        })

    return pd.DataFrame(articles)
```

## 2. Semantic Scholar API 検索

```python
def semantic_scholar_search(query, max_results=50, year_from=None,
                             fields=None):
    """
    Semantic Scholar Academic Graph API 検索。

    Parameters:
        query: str — 検索クエリ
        fields: list — 取得フィールド
    """
    url = "https://api.semanticscholar.org/graph/v1/paper/search"

    if fields is None:
        fields = ["title", "abstract", "year", "citationCount",
                  "influentialCitationCount", "authors", "url",
                  "openAccessPdf"]

    params = {
        "query": query,
        "limit": min(max_results, 100),
        "fields": ",".join(fields),
    }
    if year_from:
        params["year"] = f"{year_from}-"

    resp = requests.get(url, params=params)
    data = resp.json()
    papers = data.get("data", [])

    results = []
    for p in papers:
        results.append({
            "paper_id": p.get("paperId", ""),
            "title": p.get("title", ""),
            "abstract": (p.get("abstract") or "")[:300],
            "year": p.get("year"),
            "citations": p.get("citationCount", 0),
            "influential_citations": p.get("influentialCitationCount", 0),
            "authors": "; ".join([a.get("name", "")
                                  for a in (p.get("authors") or [])[:5]]),
            "url": p.get("url", ""),
            "open_access": bool(p.get("openAccessPdf")),
        })

    df = pd.DataFrame(results)
    print(f"Semantic Scholar: {data.get('total', 0)} total, "
          f"{len(df)} fetched")
    return df
```

## 3. OpenAlex 著者・機関メトリクス

```python
def openalex_search(query, entity_type="works", max_results=50,
                     filters=None):
    """
    OpenAlex API 検索 (250M+ works, 90M+ authors)。

    Parameters:
        query: str — 検索クエリ
        entity_type: "works", "authors", "institutions", "concepts"
        filters: dict — OpenAlex フィルタ (e.g., {"from_publication_date": "2020-01-01"})
    """
    base = "https://api.openalex.org"
    url = f"{base}/{entity_type}"

    params = {
        "search": query,
        "per-page": min(max_results, 200),
        "mailto": "research@example.com",
    }
    if filters:
        filter_parts = [f"{k}:{v}" for k, v in filters.items()]
        params["filter"] = ",".join(filter_parts)

    resp = requests.get(url, params=params)
    data = resp.json()

    results = data.get("results", [])
    total = data.get("meta", {}).get("count", 0)
    print(f"OpenAlex {entity_type}: {total} total, {len(results)} fetched")

    if entity_type == "works":
        return pd.DataFrame([{
            "id": r.get("id", ""),
            "title": r.get("title", ""),
            "year": r.get("publication_year"),
            "citations": r.get("cited_by_count", 0),
            "doi": r.get("doi", ""),
            "type": r.get("type", ""),
            "open_access": r.get("open_access", {}).get("is_oa", False),
        } for r in results])
    elif entity_type == "authors":
        return pd.DataFrame([{
            "id": r.get("id", ""),
            "name": r.get("display_name", ""),
            "works_count": r.get("works_count", 0),
            "citations": r.get("cited_by_count", 0),
            "h_index": r.get("summary_stats", {}).get("h_index", 0),
            "institution": (r.get("last_known_institutions") or [{}])[0].get(
                "display_name", "") if r.get("last_known_institutions") else "",
        } for r in results])

    return results
```

## 4. EuropePMC 全文検索

```python
def europepmc_search(query, max_results=50, source="MED",
                      open_access_only=False):
    """
    EuropePMC REST API 検索 (39M+ records)。

    Parameters:
        query: str — 検索クエリ
        source: "MED" (PubMed), "PMC" (PubMed Central), "AGR", "CBA"
        open_access_only: True でオープンアクセスのみ
    """
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

    q = query
    if open_access_only:
        q += " AND OPEN_ACCESS:y"

    params = {
        "query": q,
        "resultType": "core",
        "pageSize": min(max_results, 100),
        "format": "json",
        "src": source,
    }

    resp = requests.get(url, params=params)
    data = resp.json()
    results = data.get("resultList", {}).get("result", [])
    total = data.get("hitCount", 0)

    articles = []
    for r in results:
        articles.append({
            "pmid": r.get("pmid", ""),
            "pmcid": r.get("pmcid", ""),
            "title": r.get("title", ""),
            "journal": r.get("journalTitle", ""),
            "year": r.get("pubYear", ""),
            "citations": r.get("citedByCount", 0),
            "has_fulltext": r.get("hasTextMinedTerms", "N") == "Y",
            "is_open_access": r.get("isOpenAccess", "N") == "Y",
        })

    df = pd.DataFrame(articles)
    print(f"EuropePMC: {total} total, {len(df)} fetched")
    return df
```

## 5. CrossRef メタデータ取得

```python
def crossref_search(query, max_results=50, filter_type=None,
                     from_date=None, sort="relevance"):
    """
    CrossRef REST API 検索 (DOI メタデータ)。

    Parameters:
        query: str — 検索クエリ
        filter_type: str — DOI タイプフィルタ ("journal-article", "book-chapter")
        from_date: str — "YYYY-MM-DD"
        sort: "relevance", "published", "is-referenced-by-count"
    """
    url = "https://api.crossref.org/works"
    params = {
        "query": query,
        "rows": min(max_results, 100),
        "sort": sort,
        "mailto": "research@example.com",
    }

    filters = []
    if filter_type:
        filters.append(f"type:{filter_type}")
    if from_date:
        filters.append(f"from-pub-date:{from_date}")
    if filters:
        params["filter"] = ",".join(filters)

    resp = requests.get(url, params=params)
    data = resp.json()
    items = data.get("message", {}).get("items", [])
    total = data.get("message", {}).get("total-results", 0)

    results = []
    for item in items:
        title = item.get("title", [""])[0] if item.get("title") else ""
        results.append({
            "doi": item.get("DOI", ""),
            "title": title,
            "journal": item.get("container-title", [""])[0] if item.get("container-title") else "",
            "year": item.get("published", {}).get("date-parts", [[None]])[0][0],
            "citations": item.get("is-referenced-by-count", 0),
            "type": item.get("type", ""),
            "publisher": item.get("publisher", ""),
        })

    df = pd.DataFrame(results)
    print(f"CrossRef: {total} total, {len(df)} fetched")
    return df
```

## 6. 引用ネットワーク構築

```python
def build_citation_network(seed_pmids, depth=1, max_per_level=20):
    """
    PubMed 論文の引用/被引用ネットワーク構築。

    Parameters:
        seed_pmids: list — 起点となる PMID リスト
        depth: int — 探索深度 (1=直接引用のみ)
        max_per_level: int — 各レベルの最大ノード数
    """
    import networkx as nx

    G = nx.DiGraph()
    visited = set()
    current_level = set(seed_pmids)

    for level in range(depth + 1):
        next_level = set()
        for pmid in list(current_level)[:max_per_level]:
            if pmid in visited:
                continue
            visited.add(pmid)
            G.add_node(pmid, level=level, is_seed=(level == 0))

            # elink で引用関係取得
            time.sleep(0.4)
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
            params = {
                "dbfrom": "pubmed", "db": "pubmed",
                "id": pmid, "linkname": "pubmed_pubmed_citedin",
                "retmode": "json",
            }
            try:
                resp = requests.get(url, params=params)
                data = resp.json()
                links = data.get("linksets", [{}])[0].get("linksetdbs", [])
                if links:
                    cited_by = [l["id"] for l in links[0].get("links", [])[:max_per_level]]
                    for citing in cited_by:
                        G.add_edge(citing, pmid)
                        next_level.add(citing)
            except Exception:
                continue

        current_level = next_level - visited

    print(f"Citation network: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges, depth={depth}")
    return G
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/pubmed_search.csv` | CSV |
| `results/semantic_scholar_results.csv` | CSV |
| `results/openalex_results.csv` | CSV |
| `results/europepmc_results.csv` | CSV |
| `results/crossref_results.csv` | CSV |
| `results/citation_network.graphml` | GraphML |
| `figures/citation_network.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| PubMed | `PubMed_search_articles` | PubMed 論文検索 |
| PubMed | `PubMed_get_article` | 論文詳細取得 |
| PubMed | `PubMed_get_cited_by` | 被引用論文取得 |
| PubMed | `PubMed_get_related` | 関連論文検索 |
| PubMed | `PubMed_get_links` | 関連リンク取得 |
| PubMed | `PubMed_Guidelines_Search` | ガイドライン検索 |
| EuropePMC | `EuropePMC_search_articles` | EuropePMC 検索 |
| EuropePMC | `EuropePMC_get_article` | 記事詳細取得 |
| EuropePMC | `EuropePMC_get_references` | 参考文献取得 |
| EuropePMC | `EuropePMC_get_citations` | 引用取得 |
| EuropePMC | `EuropePMC_get_fulltext` | 全文取得 |
| EuropePMC | `EuropePMC_get_fulltext_snippets` | 全文スニペット |
| EuropePMC | `EuropePMC_Guidelines_Search` | ガイドライン検索 |
| Semantic Scholar | `SemanticScholar_search_papers` | セマンティック論文検索 |
| Semantic Scholar | `SemanticScholar_get_pdf_snippets` | PDF スニペット |
| OpenAlex | `openalex_literature_search` | 文献検索 |
| OpenAlex | `openalex_get_work` | 論文詳細 |
| OpenAlex | `openalex_get_work_by_doi` | DOI → 論文 |
| OpenAlex | `openalex_search_works` | 論文検索 (REST) |
| OpenAlex | `openalex_get_author` | 著者詳細 |
| OpenAlex | `openalex_search_authors` | 著者検索 |
| OpenAlex | `openalex_get_institution` | 機関詳細 |
| OpenAlex | `openalex_search_institutions` | 機関検索 |
| CrossRef | `Crossref_search_works` | DOI メタデータ検索 |
| CrossRef | `Crossref_get_work` | 論文メタデータ取得 |
| CrossRef | `Crossref_get_funder` | 助成機関情報 |
| CrossRef | `Crossref_get_journal` | ジャーナル情報 |
| CrossRef | `Crossref_list_funders` | 助成機関一覧 |
| CrossRef | `Crossref_list_types` | 出版タイプ一覧 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-deep-research` | 文献検索 → 深層リサーチ |
| `scientific-citation-checker` | 引用検証連携 |
| `scientific-systematic-review` | 系統的検索戦略 |
| `scientific-meta-analysis` | メタアナリシス文献収集 |
| `scientific-academic-writing` | 関連文献の自動発見 |
| `scientific-text-mining-nlp` | 抽出テキストの NLP 解析 |

### 依存パッケージ

`requests`, `pandas`, `networkx`, `xml.etree.ElementTree` (stdlib)
