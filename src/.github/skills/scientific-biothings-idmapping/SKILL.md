---
name: scientific-biothings-idmapping
description: |
  BioThings API (MyGene.info, MyVariant.info, MyChem.info) を活用した
  遺伝子・変異・化合物の横断的 ID マッピングおよびアノテーション統合スキル。
---

# Scientific BioThings ID Mapping

BioThings API スイート (MyGene, MyVariant, MyChem) を活用した
多データベース横断の ID 変換・アノテーション取得パイプラインを提供する。

## When to Use

- 遺伝子 ID 間の変換 (Entrez ↔ Ensembl ↔ Symbol ↔ UniProt) を行うとき
- 変異 ID のアノテーション (ClinVar, dbSNP, CADD 等) を取得するとき
- 化合物 ID の変換 (DrugBank ↔ ChEMBL ↔ InChIKey ↔ PubChem) を行うとき
- バッチクエリで多数の ID を一括アノテーションするとき
- 複数データベースのメタ情報を統合するとき

---

## Quick Start

## 1. MyGene.info 遺伝子アノテーション

```python
import requests
import pandas as pd

MYGENE_API = "https://mygene.info/v3"


def mygene_query(query, fields=None, species="human", size=10):
    """
    MyGene.info で遺伝子検索。

    Parameters:
        query: str — gene symbol, Entrez ID, or keyword
        fields: str | None — comma-separated fields
        species: str — "human", "mouse", etc.

    ToolUniverse:
        MyGene_query_genes(q=query, fields=fields, species=species)
    """
    params = {
        "q": query,
        "species": species,
        "size": size,
    }
    if fields:
        params["fields"] = fields

    resp = requests.get(f"{MYGENE_API}/query", params=params)
    resp.raise_for_status()
    data = resp.json()

    hits = data.get("hits", [])
    print(f"MyGene query '{query}': {data.get('total', 0)} total, "
          f"{len(hits)} returned")
    return hits


def mygene_get_gene(gene_id, fields=None):
    """
    MyGene.info 遺伝子詳細アノテーション取得。

    ToolUniverse:
        MyGene_get_gene_annotation(gene_id=gene_id, fields=fields)
    """
    params = {}
    if fields:
        params["fields"] = fields

    resp = requests.get(f"{MYGENE_API}/gene/{gene_id}", params=params)
    resp.raise_for_status()
    data = resp.json()

    print(f"MyGene gene {gene_id}: {data.get('symbol', '?')} "
          f"({data.get('name', '')})")
    return data


def mygene_batch_query(gene_ids, fields=None, species="human"):
    """
    MyGene.info バッチ遺伝子アノテーション。

    ToolUniverse:
        MyGene_batch_query(ids=gene_ids, fields=fields, species=species)
    """
    payload = {
        "ids": ",".join(str(g) for g in gene_ids),
        "species": species,
    }
    if fields:
        payload["fields"] = fields

    resp = requests.post(f"{MYGENE_API}/gene", json=payload)
    resp.raise_for_status()
    data = resp.json()

    print(f"MyGene batch: {len(gene_ids)} queried → {len(data)} results")
    return data
```

## 2. MyVariant.info 変異アノテーション

```python
MYVARIANT_API = "https://myvariant.info/v1"


def myvariant_get(variant_id, fields=None):
    """
    MyVariant.info 変異アノテーション取得。

    Parameters:
        variant_id: str — HGVS notation (e.g., "chr17:g.7674220C>T")

    ToolUniverse:
        MyVariant_get_variant_annotation(variant_id=variant_id, fields=fields)
    """
    params = {}
    if fields:
        params["fields"] = fields

    resp = requests.get(f"{MYVARIANT_API}/variant/{variant_id}", params=params)
    resp.raise_for_status()
    data = resp.json()

    clinvar = data.get("clinvar", {})
    cadd = data.get("cadd", {})
    print(f"MyVariant {variant_id}: "
          f"ClinVar={clinvar.get('clinical_significance', 'N/A')}, "
          f"CADD={cadd.get('phred', 'N/A')}")
    return data


def myvariant_query(query, fields=None, size=10):
    """
    MyVariant.info 変異検索。

    ToolUniverse:
        MyVariant_query_variants(q=query, fields=fields, size=size)
    """
    params = {"q": query, "size": size}
    if fields:
        params["fields"] = fields

    resp = requests.get(f"{MYVARIANT_API}/query", params=params)
    resp.raise_for_status()
    data = resp.json()

    hits = data.get("hits", [])
    print(f"MyVariant query '{query}': {data.get('total', 0)} total")
    return hits
```

## 3. MyChem.info 化合物アノテーション

```python
MYCHEM_API = "https://mychem.info/v1"


def mychem_get(chem_id, fields=None):
    """
    MyChem.info 化合物アノテーション取得。

    Parameters:
        chem_id: str — InChIKey, DrugBank ID, ChEMBL ID, etc.

    ToolUniverse:
        MyChem_get_chemical_annotation(chem_id=chem_id, fields=fields)
    """
    params = {}
    if fields:
        params["fields"] = fields

    resp = requests.get(f"{MYCHEM_API}/chem/{chem_id}", params=params)
    resp.raise_for_status()
    data = resp.json()

    drugbank = data.get("drugbank", {})
    print(f"MyChem {chem_id}: {drugbank.get('name', 'N/A')}")
    return data


def mychem_query(query, fields=None, size=10):
    """
    MyChem.info 化合物検索。

    ToolUniverse:
        MyChem_query_chemicals(q=query, fields=fields, size=size)
    """
    params = {"q": query, "size": size}
    if fields:
        params["fields"] = fields

    resp = requests.get(f"{MYCHEM_API}/query", params=params)
    resp.raise_for_status()
    data = resp.json()

    hits = data.get("hits", [])
    print(f"MyChem query '{query}': {data.get('total', 0)} total")
    return hits
```

## 4. クロスデータベース ID マッピング

```python
def cross_db_id_mapping(gene_symbol):
    """
    遺伝子シンボルから Entrez, Ensembl, UniProt, RefSeq を一括取得。

    ToolUniverse (横断):
        MyGene_query_genes(q=gene_symbol, fields="entrezgene,ensembl.gene,uniprot,refseq")
    """
    fields = "entrezgene,ensembl.gene,uniprot.Swiss-Prot,refseq.rna,symbol,name"
    hits = mygene_query(gene_symbol, fields=fields)

    results = []
    for hit in hits:
        ensembl = hit.get("ensembl", {})
        if isinstance(ensembl, list):
            ensembl = ensembl[0] if ensembl else {}
        uniprot = hit.get("uniprot", {})

        results.append({
            "symbol": hit.get("symbol", ""),
            "name": hit.get("name", ""),
            "entrez_id": hit.get("entrezgene", ""),
            "ensembl_gene": ensembl.get("gene", ""),
            "uniprot_swissprot": uniprot.get("Swiss-Prot", ""),
            "refseq_rna": hit.get("refseq", {}).get("rna", []),
        })

    df = pd.DataFrame(results)
    print(f"ID mapping '{gene_symbol}': {len(df)} entries")
    return df
```

## 5. バッチ統合アノテーション

```python
def batch_integrated_annotation(gene_symbols, include_variants=False):
    """
    複数遺伝子のバッチ統合アノテーション。

    ToolUniverse (横断):
        MyGene_batch_query(ids=entrez_ids, fields=fields)
        MyVariant_query_variants(q=gene_query) [optional]
    """
    # Step 1: Batch gene annotation
    all_hits = []
    for symbol in gene_symbols:
        hits = mygene_query(symbol, fields="entrezgene,symbol,name,summary")
        all_hits.extend(hits[:1])  # top hit per symbol

    df = pd.DataFrame(all_hits)
    print(f"Batch annotation: {len(gene_symbols)} genes → {len(df)} results")
    return df
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/mygene_annotation.json` | JSON |
| `results/myvariant_annotation.json` | JSON |
| `results/mychem_annotation.json` | JSON |
| `results/id_mapping.csv` | CSV |

### 利用可能ツール

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| BioThings | `MyGene_query_genes` | 遺伝子検索 |
| BioThings | `MyGene_get_gene_annotation` | 遺伝子詳細 |
| BioThings | `MyGene_batch_query` | バッチアノテーション |
| BioThings | `MyVariant_get_variant_annotation` | 変異アノテーション |
| BioThings | `MyVariant_query_variants` | 変異検索 |
| BioThings | `MyChem_get_chemical_annotation` | 化合物アノテーション |
| BioThings | `MyChem_query_chemicals` | 化合物検索 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-variant-interpretation` | 変異アノテーション |
| `scientific-gene-expression-transcriptomics` | 遺伝子発現 |
| `scientific-drug-target-interaction` | DTI 解析 |
| `scientific-rare-disease-genetics` | 希少疾患 |
| `scientific-pathway-enrichment` | パスウェイ解析 |

### 依存パッケージ

`requests`, `pandas`
