---
name: scientific-crossref-metadata
description: |
  CrossRef メタデータスキル。CrossRef REST API による
  DOI 解決・論文メタデータ・引用数・ジャーナル情報・
  助成金情報検索。ToolUniverse 連携: crossref。
tu_tools:
  - key: crossref
    name: CrossRef
    description: DOI 解決・論文メタデータ・引用数・ジャーナル情報
---

# Scientific CrossRef Metadata

CrossRef REST API を活用した学術文献 DOI 解決・メタデータ検索・
引用分析・ジャーナル情報・助成金レジストリ検索パイプラインを
提供する。

## When to Use

- DOI から論文メタデータを取得するとき
- 学術文献をタイトル・著者で検索するとき
- ジャーナルの ISSN やインパクト情報を調べるとき
- 論文の引用数・被引用数を確認するとき
- 研究助成金の情報を検索するとき
- 参考文献リストのメタデータを一括取得するとき
- 特定出版社やジャーナルの出版傾向を分析するとき

---

## Quick Start

## 1. DOI 解決・論文メタデータ

```python
import requests
import pandas as pd

CR_BASE = "https://api.crossref.org"
CR_HEADERS = {
    "User-Agent": "SATORI/0.18.0 (mailto:your@email.com)",
}


def crossref_resolve_doi(doi):
    """
    CrossRef — DOI からメタデータ取得。

    Parameters:
        doi: str — DOI (例: "10.1038/s41586-020-2649-2")
    """
    url = f"{CR_BASE}/works/{doi}"
    resp = requests.get(url, headers=CR_HEADERS, timeout=30)
    resp.raise_for_status()
    item = resp.json().get("message", {})

    authors = []
    for a in item.get("author", []):
        name = f"{a.get('given', '')} {a.get('family', '')}"
        authors.append(name.strip())

    result = {
        "doi": item.get("DOI", ""),
        "title": " ".join(item.get("title", [])),
        "authors": "; ".join(authors[:10]),
        "journal": " ".join(
            item.get("container-title", [])),
        "publisher": item.get("publisher", ""),
        "type": item.get("type", ""),
        "published_date": _cr_date(
            item.get("published-print") or
            item.get("published-online")),
        "citation_count": item.get(
            "is-referenced-by-count", 0),
        "reference_count": item.get("reference-count", 0),
        "issn": ", ".join(item.get("ISSN", [])),
        "url": item.get("URL", ""),
        "abstract": (item.get("abstract") or "")[:500],
        "funder": "; ".join(
            f.get("name", "")
            for f in item.get("funder", [])),
        "license": _cr_license(item),
    }

    print(f"CrossRef DOI: {doi}")
    print(f"  {result['title'][:80]}")
    print(f"  Citations: {result['citation_count']}")
    return result


def _cr_date(date_obj):
    if not date_obj:
        return ""
    parts = date_obj.get("date-parts", [[]])[0]
    return "-".join(str(p) for p in parts)


def _cr_license(item):
    licenses = item.get("license", [])
    if licenses:
        return licenses[0].get("content-version", "")
    return ""
```

## 2. 論文検索

```python
def crossref_search_works(query, limit=50,
                             sort="relevance",
                             filter_type=None,
                             from_date=None):
    """
    CrossRef — 論文検索。

    Parameters:
        query: str — 検索クエリ
        limit: int — 最大結果数
        sort: str — ソート ("relevance", "published",
            "is-referenced-by-count")
        filter_type: str — 文献タイプフィルタ
            (例: "journal-article")
        from_date: str — 開始日 (例: "2020-01-01")
    """
    url = f"{CR_BASE}/works"
    params = {
        "query": query,
        "rows": min(limit, 1000),
        "sort": sort,
    }

    filters = []
    if filter_type:
        filters.append(f"type:{filter_type}")
    if from_date:
        filters.append(f"from-pub-date:{from_date}")
    if filters:
        params["filter"] = ",".join(filters)

    resp = requests.get(url, params=params,
                        headers=CR_HEADERS, timeout=30)
    resp.raise_for_status()
    data = resp.json().get("message", {})

    results = []
    for item in data.get("items", []):
        authors = []
        for a in item.get("author", []):
            name = f"{a.get('given', '')} {a.get('family', '')}"
            authors.append(name.strip())
        results.append({
            "doi": item.get("DOI", ""),
            "title": " ".join(item.get("title", [])),
            "authors": "; ".join(authors[:5]),
            "journal": " ".join(
                item.get("container-title", [])),
            "year": _cr_date(
                item.get("published-print") or
                item.get("published-online")),
            "citations": item.get(
                "is-referenced-by-count", 0),
            "type": item.get("type", ""),
        })

    df = pd.DataFrame(results)
    total = data.get("total-results", 0)
    print(f"CrossRef search: {len(df)}/{total} works "
          f"(query='{query}')")
    return df
```

## 3. ジャーナル情報・助成金検索

```python
def crossref_journal_info(issn):
    """
    CrossRef — ジャーナル情報取得。

    Parameters:
        issn: str — ISSN (例: "0028-0836")
    """
    url = f"{CR_BASE}/journals/{issn}"
    resp = requests.get(url, headers=CR_HEADERS, timeout=30)
    resp.raise_for_status()
    data = resp.json().get("message", {})

    counts = data.get("counts", {})
    result = {
        "issn": issn,
        "title": data.get("title", ""),
        "publisher": data.get("publisher", ""),
        "subjects": "; ".join(
            s.get("name", "")
            for s in data.get("subjects", [])),
        "total_dois": counts.get("total-dois", 0),
        "current_dois": counts.get("current-dois", 0),
        "backfile_dois": counts.get("backfile-dois", 0),
    }

    print(f"CrossRef journal: {result['title']} "
          f"({result['total_dois']} DOIs)")
    return result


def crossref_search_funders(query, limit=20):
    """
    CrossRef — 助成金機関検索。

    Parameters:
        query: str — 機関名 (例: "NIH", "JSPS")
        limit: int — 最大結果数
    """
    url = f"{CR_BASE}/funders"
    params = {"query": query, "rows": limit}
    resp = requests.get(url, params=params,
                        headers=CR_HEADERS, timeout=30)
    resp.raise_for_status()
    data = resp.json().get("message", {})

    results = []
    for item in data.get("items", []):
        results.append({
            "funder_id": item.get("id", ""),
            "name": item.get("name", ""),
            "location": item.get("location", ""),
            "alt_names": "; ".join(
                item.get("alt-names", [])[:3]),
            "work_count": item.get("work-count", 0),
        })

    df = pd.DataFrame(results)
    print(f"CrossRef funders: {len(df)} (query='{query}')")
    return df
```

## 4. CrossRef 統合パイプライン

```python
def crossref_pipeline(query, dois=None,
                         output_dir="results"):
    """
    CrossRef 統合パイプライン。

    Parameters:
        query: str — 検索クエリ
        dois: list[str] — DOI リスト (直接解決)
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 論文検索
    works = crossref_search_works(query)
    works.to_csv(output_dir / "works.csv", index=False)

    # 2) DOI 一括解決
    if dois:
        resolved = []
        for doi in dois:
            try:
                meta = crossref_resolve_doi(doi)
                resolved.append(meta)
            except Exception as e:
                print(f"  Warning: {doi} — {e}")
                continue
        resolved_df = pd.DataFrame(resolved)
        resolved_df.to_csv(output_dir / "doi_resolved.csv",
                           index=False)

    # 3) 引用分析
    if not works.empty:
        stats = {
            "total_works": len(works),
            "total_citations": works["citations"].sum(),
            "mean_citations": works["citations"].mean(),
            "median_citations": works["citations"].median(),
            "max_citations": works["citations"].max(),
        }
        pd.DataFrame([stats]).to_csv(
            output_dir / "citation_stats.csv", index=False)

    print(f"CrossRef pipeline: {output_dir}")
    return {"works": works}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `crossref` | CrossRef | DOI 解決・メタデータ・引用・ジャーナル情報 |

## パイプライン統合

```
literature-search → crossref-metadata → citation-checker
  (PubMed/NCBI)   (CrossRef REST API)  (引用品質検証)
       │                  │                  ↓
  semantic-scholar ──────┘          deep-research
  (S2 Academic Graph)   │          (知識統合)
                        ↓
               bibliometrics
               (書誌計量分析)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/works.csv` | 論文検索結果 | → semantic-scholar |
| `results/doi_resolved.csv` | DOI メタデータ | → citation-checker |
| `results/citation_stats.csv` | 引用統計 | → bibliometrics |
