---
name: scientific-semantic-scholar
description: |
  Semantic Scholar 学術グラフスキル。Semantic Scholar Academic
  Graph API による論文検索・著者プロファイル・引用グラフ・
  推薦・TLDR 要約。ToolUniverse 連携: semantic_scholar。
tu_tools:
  - key: semantic_scholar
    name: Semantic Scholar
    description: 学術論文検索・引用解析・著者プロファイル
---

# Scientific Semantic Scholar

Semantic Scholar Academic Graph API を活用した学術論文検索・
引用ネットワーク解析・著者プロファイル・論文推薦パイプライン
を提供する。

## When to Use

- 学術論文を高精度で検索するとき
- 引用・被引用ネットワークを解析するとき
- 著者の h-index・論文数・研究領域を調べるとき
- 関連論文の推薦を受けるとき
- TLDR (自動要約) を取得するとき
- 特定分野の引用傾向を分析するとき
- PubMed/OpenAlex 以外の学術検索エンジンを使うとき

---

## Quick Start

## 1. 論文検索

```python
import requests
import pandas as pd

S2_BASE = "https://api.semanticscholar.org/graph/v1"
S2_HEADERS = {}  # API key: {"x-api-key": "YOUR_KEY"}


def semantic_scholar_search(query, limit=50,
                              year_range=None,
                              fields_of_study=None):
    """
    Semantic Scholar — 論文検索。

    Parameters:
        query: str — 検索クエリ
        limit: int — 最大結果数
        year_range: str — 年範囲 (例: "2020-2024")
        fields_of_study: list[str] — 分野フィルタ
    """
    url = f"{S2_BASE}/paper/search"
    params = {
        "query": query,
        "limit": min(limit, 100),
        "fields": ("paperId,title,year,citationCount,"
                    "influentialCitationCount,authors,"
                    "journal,tldr,openAccessPdf,fieldsOfStudy"),
    }
    if year_range:
        params["year"] = year_range
    if fields_of_study:
        params["fieldsOfStudy"] = ",".join(fields_of_study)

    resp = requests.get(url, params=params,
                        headers=S2_HEADERS, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for p in data.get("data", []):
        authors = [a.get("name", "") for a in p.get("authors", [])]
        tldr_text = ""
        if p.get("tldr"):
            tldr_text = p["tldr"].get("text", "")
        results.append({
            "paper_id": p.get("paperId", ""),
            "title": p.get("title", ""),
            "year": p.get("year"),
            "citation_count": p.get("citationCount", 0),
            "influential_citations": p.get(
                "influentialCitationCount", 0),
            "authors": "; ".join(authors[:5]),
            "journal": (p.get("journal") or {}).get("name", ""),
            "fields": ", ".join(p.get("fieldsOfStudy") or []),
            "tldr": tldr_text[:300],
            "pdf_url": (p.get("openAccessPdf") or {}).get("url", ""),
        })

    df = pd.DataFrame(results)
    print(f"Semantic Scholar: {len(df)} papers "
          f"(query='{query}')")
    return df


def semantic_scholar_get_paper(paper_id):
    """
    Semantic Scholar — 論文詳細取得。

    Parameters:
        paper_id: str — S2 Paper ID / DOI / ArXiv ID
    """
    url = f"{S2_BASE}/paper/{paper_id}"
    params = {
        "fields": ("paperId,title,year,abstract,citationCount,"
                    "influentialCitationCount,authors,references,"
                    "citations,journal,tldr,openAccessPdf,"
                    "fieldsOfStudy,publicationDate,venue"),
    }
    resp = requests.get(url, params=params,
                        headers=S2_HEADERS, timeout=30)
    resp.raise_for_status()
    return resp.json()
```

## 2. 著者プロファイル・引用解析

```python
def semantic_scholar_author(author_id, paper_limit=100):
    """
    Semantic Scholar — 著者プロファイル取得。

    Parameters:
        author_id: str — S2 Author ID
        paper_limit: int — 取得論文数上限
    """
    url = f"{S2_BASE}/author/{author_id}"
    params = {
        "fields": ("authorId,name,affiliations,homepage,"
                    "paperCount,citationCount,hIndex"),
    }
    resp = requests.get(url, params=params,
                        headers=S2_HEADERS, timeout=30)
    resp.raise_for_status()
    profile = resp.json()

    # 論文一覧
    papers_url = f"{S2_BASE}/author/{author_id}/papers"
    p_params = {
        "fields": "paperId,title,year,citationCount,venue",
        "limit": min(paper_limit, 1000),
    }
    p_resp = requests.get(papers_url, params=p_params,
                          headers=S2_HEADERS, timeout=30)
    p_resp.raise_for_status()

    papers = []
    for p in p_resp.json().get("data", []):
        papers.append({
            "paper_id": p.get("paperId", ""),
            "title": p.get("title", ""),
            "year": p.get("year"),
            "citations": p.get("citationCount", 0),
            "venue": p.get("venue", ""),
        })

    papers_df = pd.DataFrame(papers)

    print(f"Author {profile.get('name', '')}: "
          f"h-index={profile.get('hIndex', 0)}, "
          f"{profile.get('paperCount', 0)} papers, "
          f"{profile.get('citationCount', 0)} citations")
    return profile, papers_df
```

## 3. 引用ネットワーク・影響度分析

```python
def semantic_scholar_citation_graph(paper_id,
                                       direction="both",
                                       limit=100):
    """
    Semantic Scholar — 引用グラフ取得。

    Parameters:
        paper_id: str — S2 Paper ID
        direction: str — "citations", "references", "both"
        limit: int — 各方向の上限
    """
    graphs = {}
    fields = "paperId,title,year,citationCount,authors"

    if direction in ("citations", "both"):
        url = f"{S2_BASE}/paper/{paper_id}/citations"
        resp = requests.get(url, params={"fields": fields,
                            "limit": limit},
                            headers=S2_HEADERS, timeout=30)
        resp.raise_for_status()
        cites = []
        for c in resp.json().get("data", []):
            cp = c.get("citingPaper", {})
            cites.append({
                "paper_id": cp.get("paperId", ""),
                "title": cp.get("title", ""),
                "year": cp.get("year"),
                "citations": cp.get("citationCount", 0),
            })
        graphs["citations"] = pd.DataFrame(cites)

    if direction in ("references", "both"):
        url = f"{S2_BASE}/paper/{paper_id}/references"
        resp = requests.get(url, params={"fields": fields,
                            "limit": limit},
                            headers=S2_HEADERS, timeout=30)
        resp.raise_for_status()
        refs = []
        for r in resp.json().get("data", []):
            rp = r.get("citedPaper", {})
            refs.append({
                "paper_id": rp.get("paperId", ""),
                "title": rp.get("title", ""),
                "year": rp.get("year"),
                "citations": rp.get("citationCount", 0),
            })
        graphs["references"] = pd.DataFrame(refs)

    for k, v in graphs.items():
        print(f"  {k}: {len(v)} papers")
    return graphs
```

## 4. 学術文献統合パイプライン

```python
def semantic_scholar_pipeline(query, year_range=None,
                                 output_dir="results"):
    """
    Semantic Scholar 統合パイプライン。

    Parameters:
        query: str — 検索クエリ
        year_range: str — 年範囲
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 論文検索
    papers = semantic_scholar_search(query,
                                     year_range=year_range)
    papers.to_csv(output_dir / "papers.csv", index=False)

    # 2) トップ被引用論文の引用グラフ
    if not papers.empty:
        top = papers.sort_values("citation_count",
                                  ascending=False).iloc[0]
        pid = top["paper_id"]
        graphs = semantic_scholar_citation_graph(pid)
        for k, df in graphs.items():
            df.to_csv(output_dir / f"{k}.csv", index=False)

    # 3) 年次引用傾向
    if not papers.empty and "year" in papers.columns:
        yearly = papers.groupby("year").agg(
            papers_count=("paper_id", "count"),
            total_citations=("citation_count", "sum"),
            avg_citations=("citation_count", "mean"),
        ).reset_index()
        yearly.to_csv(output_dir / "yearly_trend.csv",
                      index=False)

    print(f"Semantic Scholar pipeline: {output_dir}")
    return {"papers": papers}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `semantic_scholar` | Semantic Scholar | 論文検索・引用解析・著者・TLDR |

## パイプライン統合

```
literature-search → semantic-scholar → deep-research
  (PubMed/NCBI)   (Academic Graph API) (knowledge synthesis)
       │                  │                    ↓
  crossref-metadata ─────┘             citation-checker
  (DOI/metadata)     │                 (引用品質検証)
                     ↓
            gene-expression-transcriptomics
            (論文引用データからの解析)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/papers.csv` | 論文検索結果 | → deep-research |
| `results/citations.csv` | 被引用論文 | → citation-checker |
| `results/references.csv` | 引用論文 | → meta-analysis |
| `results/yearly_trend.csv` | 年次引用傾向 | → bibliometrics |
