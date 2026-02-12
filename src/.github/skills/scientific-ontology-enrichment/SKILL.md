---
name: scientific-ontology-enrichment
description: |
  オントロジー・エンリッチメント解析スキル。EFO 実験ファクターオントロジー、
  OLS オントロジー検索サービス、Enrichr 遺伝子セット濃縮解析、
  UMLS メタシソーラス統一医学言語体系の統合パイプライン。
---

# Scientific Ontology Enrichment

EFO / OLS / Enrichr / UMLS を統合した
オントロジー検索・エンリッチメント解析パイプラインを提供する。

## When to Use

- EFO で実験条件 (疾患・細胞型・組織) のオントロジー ID を取得するとき
- OLS で複数オントロジー横断検索 (HP, MONDO, DOID, GO, CHEBI) するとき
- Enrichr で遺伝子リストの濃縮解析を行うとき
- UMLS CUI で異なる用語体系間のマッピングを行うとき
- GWAS Catalog の trait を EFO 用語で標準化するとき

---

## Quick Start

## 1. EFO 実験ファクターオントロジー

```python
import requests
import pandas as pd

OLS_API = "https://www.ebi.ac.uk/ols4/api"


def search_efo(query, exact=False):
    """
    EFO (Experimental Factor Ontology) 検索。

    Parameters:
        query: str — 検索語 (疾患名、細胞型、組織名等)
        exact: bool — 完全一致検索

    ToolUniverse:
        EFO_search(query=query, exact=exact)
    """
    params = {
        "q": query,
        "ontology": "efo",
        "exact": str(exact).lower(),
        "rows": 30,
    }
    resp = requests.get(f"{OLS_API}/search", params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for doc in data.get("response", {}).get("docs", []):
        results.append({
            "efo_id": doc.get("obo_id", ""),
            "label": doc.get("label", ""),
            "description": (doc.get("description") or [""])[0][:200],
            "iri": doc.get("iri", ""),
            "ontology": doc.get("ontology_name", ""),
            "is_defining_ontology": doc.get("is_defining_ontology", False),
            "synonyms": doc.get("synonym", []),
        })

    df = pd.DataFrame(results)
    print(f"EFO search '{query}': {len(df)} terms")
    return df
```

## 2. OLS マルチオントロジー検索

```python
def search_ols(query, ontologies=None, type_filter=None):
    """
    OLS (Ontology Lookup Service) マルチオントロジー横断検索。

    Parameters:
        query: str — 検索語
        ontologies: list — オントロジー ID リスト (e.g., ["hp", "mondo", "go"])
        type_filter: str — "class", "property", "individual"

    ToolUniverse:
        OLS_search(query=query, ontology=ontology)
        OLS_get_term(ontology=ontology, iri=iri)
        OLS_get_ancestors(ontology=ontology, iri=iri)
    """
    params = {"q": query, "rows": 50}
    if ontologies:
        params["ontology"] = ",".join(ontologies)
    if type_filter:
        params["type"] = type_filter

    resp = requests.get(f"{OLS_API}/search", params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for doc in data.get("response", {}).get("docs", []):
        results.append({
            "obo_id": doc.get("obo_id", ""),
            "label": doc.get("label", ""),
            "ontology": doc.get("ontology_name", ""),
            "description": (doc.get("description") or [""])[0][:200],
            "iri": doc.get("iri", ""),
            "synonyms": doc.get("synonym", []),
            "has_children": doc.get("has_children", False),
        })

    df = pd.DataFrame(results)
    print(f"OLS search '{query}' "
          f"[{','.join(ontologies) if ontologies else 'all'}]: "
          f"{len(df)} terms")
    return df


def get_ols_term_hierarchy(ontology, term_id):
    """
    OLS 用語の階層構造 (ancestors/descendants) 取得。

    Parameters:
        ontology: str — オントロジー ID (e.g., "hp", "go")
        term_id: str — OBO ID (e.g., "HP:0001250")
    """
    iri = f"http://purl.obolibrary.org/obo/{term_id.replace(':', '_')}"
    encoded_iri = requests.utils.quote(requests.utils.quote(iri, safe=""), safe="")

    # Ancestors
    anc_resp = requests.get(
        f"{OLS_API}/ontologies/{ontology}/terms/{encoded_iri}/ancestors"
    )

    # Descendants
    desc_resp = requests.get(
        f"{OLS_API}/ontologies/{ontology}/terms/{encoded_iri}/descendants"
    )

    hierarchy = {"ancestors": [], "descendants": []}

    if anc_resp.status_code == 200:
        for t in anc_resp.json().get("_embedded", {}).get("terms", []):
            hierarchy["ancestors"].append({
                "id": t.get("obo_id", ""),
                "label": t.get("label", ""),
            })

    if desc_resp.status_code == 200:
        for t in desc_resp.json().get("_embedded", {}).get("terms", []):
            hierarchy["descendants"].append({
                "id": t.get("obo_id", ""),
                "label": t.get("label", ""),
            })

    print(f"OLS hierarchy {term_id}: "
          f"{len(hierarchy['ancestors'])} ancestors, "
          f"{len(hierarchy['descendants'])} descendants")
    return hierarchy
```

## 3. Enrichr 遺伝子セット濃縮解析

```python
ENRICHR_API = "https://maayanlab.cloud/Enrichr"


def run_enrichr(gene_list, description="", gene_set_libraries=None):
    """
    Enrichr 遺伝子リスト濃縮解析。

    Parameters:
        gene_list: list — 遺伝子シンボルリスト (e.g., ["TP53", "BRCA1", "EGFR"])
        description: str — 解析の説明
        gene_set_libraries: list — 使用する遺伝子セットライブラリ

    ToolUniverse:
        Enrichr_submit_gene_list(genes=gene_list)
        Enrichr_get_enrichment(user_list_id=id, library=library)
    """
    if gene_set_libraries is None:
        gene_set_libraries = [
            "GO_Biological_Process_2023",
            "GO_Molecular_Function_2023",
            "KEGG_2021_Human",
            "Reactome_2022",
            "WikiPathway_2023_Human",
            "DisGeNET",
        ]

    # Submit gene list
    genes_str = "\n".join(gene_list)
    submit_resp = requests.post(
        f"{ENRICHR_API}/addList",
        files={"list": (None, genes_str), "description": (None, description)},
    )
    submit_resp.raise_for_status()
    user_list_id = submit_resp.json().get("userListId")
    print(f"Enrichr: submitted {len(gene_list)} genes (ID={user_list_id})")

    # Get enrichment results per library
    all_results = {}
    for library in gene_set_libraries:
        enrich_resp = requests.get(
            f"{ENRICHR_API}/enrich",
            params={"userListId": user_list_id, "backgroundType": library},
        )
        enrich_resp.raise_for_status()
        data = enrich_resp.json()

        results = []
        for term_data in data.get(library, []):
            results.append({
                "rank": term_data[0],
                "term": term_data[1],
                "p_value": term_data[2],
                "z_score": term_data[3],
                "combined_score": term_data[4],
                "overlap_genes": term_data[5],
                "adjusted_p": term_data[6],
            })

        df = pd.DataFrame(results)
        if not df.empty:
            df = df.sort_values("adjusted_p")
        all_results[library] = df
        sig_count = (df["adjusted_p"] < 0.05).sum() if not df.empty else 0
        print(f"  {library}: {sig_count} significant terms (FDR < 0.05)")

    return all_results
```

## 4. UMLS メタシソーラスマッピング

```python
UMLS_API = "https://uts-ws.nlm.nih.gov/rest"


def search_umls(query, api_key, search_type="words"):
    """
    UMLS メタシソーラス検索。

    Parameters:
        query: str — 検索語 (疾患名、症状、薬剤名)
        api_key: str — UMLS API キー
        search_type: str — "words", "exact", "leftTruncation"

    ToolUniverse:
        UMLS_search(query=query, search_type=search_type)
        UMLS_get_concept(cui=cui)
    """
    params = {
        "string": query,
        "searchType": search_type,
        "apiKey": api_key,
        "pageSize": 25,
    }
    resp = requests.get(f"{UMLS_API}/search/current", params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for item in data.get("result", {}).get("results", []):
        results.append({
            "cui": item.get("ui", ""),
            "name": item.get("name", ""),
            "root_source": item.get("rootSource", ""),
            "uri": item.get("uri", ""),
        })

    df = pd.DataFrame(results)
    print(f"UMLS search '{query}': {len(df)} concepts")
    return df


def get_umls_crosswalk(cui, api_key, target_source=None):
    """
    UMLS CUI からの用語体系間マッピング。

    Parameters:
        cui: str — UMLS CUI (e.g., "C0023264")
        api_key: str — UMLS API キー
        target_source: str — ターゲット用語体系 (e.g., "SNOMEDCT_US", "ICD10CM", "MeSH")
    """
    params = {"apiKey": api_key, "pageSize": 100}
    if target_source:
        params["sabs"] = target_source

    resp = requests.get(f"{UMLS_API}/content/current/CUI/{cui}/atoms", params=params)
    resp.raise_for_status()
    data = resp.json()

    mappings = []
    for atom in data.get("result", []):
        mappings.append({
            "source": atom.get("rootSource", ""),
            "code": atom.get("sourceConcept", ""),
            "name": atom.get("name", ""),
            "term_type": atom.get("termType", ""),
        })

    df = pd.DataFrame(mappings)
    if target_source:
        df = df[df["source"] == target_source]

    print(f"UMLS crosswalk {cui}: {len(df)} mappings "
          f"({target_source or 'all sources'})")
    return df
```

---

## 利用可能ツール

| ToolUniverse カテゴリ | 主なツール |
|---|---|
| `efo` | `EFO_search` |
| `ols` | `OLS_search`, `OLS_get_term`, `OLS_get_ancestors` |
| `enrichr` | `Enrichr_submit_gene_list`, `Enrichr_get_enrichment` |
| `umls` | `UMLS_search`, `UMLS_get_concept` |

## パイプライン出力

| 出力ファイル | 説明 | 連携先スキル |
|---|---|---|
| `results/efo_terms.csv` | EFO 標準化用語 | → disease-research, gene-expression |
| `results/enrichr_results/` | 遺伝子セット濃縮結果 | → pathway-enrichment, multi-omics |
| `results/umls_mapping.json` | UMLS 用語マッピング | → clinical-decision-support, public-health-data |
| `results/ontology_hierarchy.json` | オントロジー階層 | → text-mining-nlp, knowledge-graph |

## パイプライン統合

```
disease-research ──→ ontology-enrichment ──→ pathway-enrichment
  (GWAS/DisGeNET)    (EFO/OLS/UMLS/Enrichr)   (KEGG/Reactome/GO)
                              │
                              ├──→ biothings-idmapping (CUI→Gene→Protein)
                              ├──→ public-health-data (UMLS→RxNorm)
                              └──→ clinical-reporting (SNOMED/ICD マッピング)
```
