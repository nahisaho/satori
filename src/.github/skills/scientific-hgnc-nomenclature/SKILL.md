---
name: scientific-hgnc-nomenclature
description: |
  HGNC 遺伝子命名法スキル。HUGO Gene Nomenclature Committee
  REST API による公式遺伝子シンボル検索・エイリアス解決・
  遺伝子ファミリー/グループクエリ・ID クロスリファレンス
  パイプライン。
  TU 外スキル (直接 REST API)。
---

# Scientific HGNC Nomenclature

HGNC (HUGO Gene Nomenclature Committee) REST API を活用した
公式遺伝子シンボル検索・エイリアス/旧シンボル解決・
遺伝子ファミリー照会・マルチデータベース ID 相互参照
パイプラインを提供する。

## When to Use

- 遺伝子エイリアスから公式 HGNC シンボルを取得するとき
- 旧遺伝子シンボル (previous symbol) を最新名に変換するとき
- 遺伝子ファミリー/グループのメンバーリストを取得するとき
- HGNC ID ↔ Ensembl / NCBI Gene / UniProt のクロスリファレンスを行うとき
- 遺伝子座タイプ (protein-coding, ncRNA 等) でフィルタするとき

---

## Quick Start

## 1. HGNC シンボル検索

```python
import requests
import pandas as pd

HGNC_BASE = "https://rest.genenames.org"
HEADERS = {"Accept": "application/json"}


def hgnc_search(query):
    """
    HGNC — 遺伝子シンボル/名前検索。

    Parameters:
        query: str — 検索クエリ (シンボル/名前)
    """
    url = f"{HGNC_BASE}/search/{query}"
    resp = requests.get(url, headers=HEADERS,
                        timeout=30)
    resp.raise_for_status()
    data = resp.json().get("response", {})
    docs = data.get("docs", [])

    rows = []
    for doc in docs:
        rows.append({
            "hgnc_id": doc.get("hgnc_id", ""),
            "symbol": doc.get("symbol", ""),
            "name": doc.get("name", ""),
            "locus_type": doc.get("locus_type", ""),
            "status": doc.get("status", ""),
        })

    df = pd.DataFrame(rows)
    print(f"HGNC search '{query}': {len(df)} hits")
    return df


def hgnc_fetch_symbol(symbol):
    """
    HGNC — 公式シンボルで遺伝子詳細取得。

    Parameters:
        symbol: str — 公式遺伝子シンボル (例: "BRCA1")
    """
    url = f"{HGNC_BASE}/fetch/symbol/{symbol}"
    resp = requests.get(url, headers=HEADERS,
                        timeout=30)
    resp.raise_for_status()
    docs = resp.json().get("response", {}).get(
        "docs", [])

    if not docs:
        print(f"HGNC: {symbol} not found")
        return {}

    doc = docs[0]
    info = {
        "hgnc_id": doc.get("hgnc_id", ""),
        "symbol": doc.get("symbol", ""),
        "name": doc.get("name", ""),
        "locus_type": doc.get("locus_type", ""),
        "location": doc.get("location", ""),
        "alias_symbol": doc.get("alias_symbol", []),
        "prev_symbol": doc.get("prev_symbol", []),
        "ensembl_gene_id": doc.get(
            "ensembl_gene_id", ""),
        "entrez_id": doc.get("entrez_id", ""),
        "uniprot_ids": doc.get("uniprot_ids", []),
        "gene_group": doc.get("gene_group", []),
    }

    print(f"HGNC: {symbol} → {info['name']} "
          f"({info['locus_type']})")
    return info
```

## 2. エイリアス/旧シンボル解決

```python
def hgnc_resolve_alias(alias):
    """
    HGNC — エイリアスから公式シンボルへ解決。

    Parameters:
        alias: str — エイリアスまたは旧シンボル
    """
    # 1) alias_symbol で検索
    url = f"{HGNC_BASE}/fetch/alias_symbol/{alias}"
    resp = requests.get(url, headers=HEADERS,
                        timeout=30)
    resp.raise_for_status()
    docs = resp.json().get("response", {}).get(
        "docs", [])

    if docs:
        symbols = [d["symbol"] for d in docs]
        print(f"HGNC alias '{alias}' → "
              f"{', '.join(symbols)}")
        return symbols

    # 2) prev_symbol で検索
    url2 = f"{HGNC_BASE}/fetch/prev_symbol/{alias}"
    resp2 = requests.get(url2, headers=HEADERS,
                         timeout=30)
    resp2.raise_for_status()
    docs2 = resp2.json().get("response", {}).get(
        "docs", [])

    if docs2:
        symbols = [d["symbol"] for d in docs2]
        print(f"HGNC prev '{alias}' → "
              f"{', '.join(symbols)}")
        return symbols

    print(f"HGNC: '{alias}' not resolved")
    return []


def hgnc_resolve_batch(aliases):
    """
    HGNC — バッチエイリアス解決。

    Parameters:
        aliases: list[str] — エイリアス/旧シンボルリスト
    """
    results = []
    for alias in aliases:
        resolved = hgnc_resolve_alias(alias)
        results.append({
            "input": alias,
            "resolved": resolved[0] if resolved
                        else "UNRESOLVED",
            "ambiguous": len(resolved) > 1,
        })

    df = pd.DataFrame(results)
    n_resolved = (df["resolved"] != "UNRESOLVED").sum()
    print(f"HGNC batch: {n_resolved}/{len(df)} "
          f"resolved")
    return df
```

## 3. 遺伝子ファミリー/グループ

```python
def hgnc_gene_group(group_name):
    """
    HGNC — 遺伝子ファミリー/グループメンバー取得。

    Parameters:
        group_name: str — グループ名
                          (例: "Kinases", "Ion channels")
    """
    url = (f"{HGNC_BASE}/search/"
           f"gene_group:%22{group_name}%22")
    resp = requests.get(url, headers=HEADERS,
                        timeout=30)
    resp.raise_for_status()
    docs = resp.json().get("response", {}).get(
        "docs", [])

    rows = []
    for doc in docs:
        rows.append({
            "symbol": doc.get("symbol", ""),
            "name": doc.get("name", ""),
            "locus_type": doc.get("locus_type", ""),
            "location": doc.get("location", ""),
        })

    df = pd.DataFrame(rows)
    print(f"HGNC group '{group_name}': "
          f"{len(df)} members")
    return df
```

## 4. HGNC 統合パイプライン

```python
def hgnc_pipeline(symbols, aliases=None,
                    output_dir="results"):
    """
    HGNC 統合命名法パイプライン。

    Parameters:
        symbols: list[str] — 公式シンボルリスト
        aliases: list[str] | None — 解決するエイリアス
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) シンボル詳細
    details = []
    for sym in symbols:
        info = hgnc_fetch_symbol(sym)
        if info:
            details.append(info)
    detail_df = pd.DataFrame(details)
    detail_df.to_csv(
        output_dir / "hgnc_details.csv",
        index=False)

    # 2) エイリアス解決
    if aliases:
        alias_df = hgnc_resolve_batch(aliases)
        alias_df.to_csv(
            output_dir / "hgnc_alias_resolved.csv",
            index=False)

    # 3) ID クロスリファレンス
    xref_rows = []
    for d in details:
        xref_rows.append({
            "symbol": d.get("symbol", ""),
            "hgnc_id": d.get("hgnc_id", ""),
            "ensembl": d.get("ensembl_gene_id", ""),
            "entrez": d.get("entrez_id", ""),
            "uniprot": (d.get("uniprot_ids", [""])[0]
                        if d.get("uniprot_ids")
                        else ""),
        })
    xref_df = pd.DataFrame(xref_rows)
    xref_df.to_csv(
        output_dir / "hgnc_xref.csv",
        index=False)

    print(f"HGNC pipeline → {output_dir}")
    return {"details": detail_df, "xref": xref_df}
```

---

## パイプライン統合

```
biothings-idmapping → hgnc-nomenclature → genome-sequence-tools
  (MyGene/MyVariant)     (公式シンボル)       (配列解析)
         │                     │                  ↓
  gene-expression ────────────┘       variant-interpretation
    (RNA-seq)                          (バリアント解釈)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/hgnc_details.csv` | 遺伝子詳細 | → gene-expression |
| `results/hgnc_alias_resolved.csv` | エイリアス解決 | → biothings-idmapping |
| `results/hgnc_xref.csv` | ID 相互参照 | → genome-sequence-tools |
