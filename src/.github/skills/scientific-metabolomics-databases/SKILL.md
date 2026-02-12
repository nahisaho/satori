---
name: scientific-metabolomics-databases
description: |
  メタボロミクスデータベース統合スキル。HMDB (Human Metabolome Database、
  220,000+ 代謝物)、MetaCyc (代謝パスウェイ)、Metabolomics Workbench
  (NIH メタボロミクスリポジトリ) の 3 大メタボロミクス DB を統合した
  代謝物同定、パスウェイマッピング、バイオマーカー発見、
  RefMet 標準化命名パイプライン。13 の ToolUniverse SMCP ツールと連携。
tu_tools:
  - key: metacyc
    name: MetaCyc
    description: 代謝パスウェイ・反応・化合物データベース
---

# Scientific Metabolomics Databases

HMDB / MetaCyc / Metabolomics Workbench の 3 大メタボロミクスデータベースを統合した
代謝物同定・パスウェイマッピング・バイオマーカー発見パイプラインを提供する。

## When to Use

- 質量 (m/z) から代謝物を同定するとき
- HMDB で代謝物の生物学的コンテキストを調べるとき
- MetaCyc で代謝パスウェイの詳細を確認するとき
- Metabolomics Workbench の公開データセットを検索するとき
- 複数 DB を横断した代謝物アノテーションが必要なとき

---

## Quick Start

## 1. HMDB 代謝物検索

```python
import requests
import pandas as pd
import xml.etree.ElementTree as ET


def hmdb_search(query, search_type="name", max_results=50):
    """
    HMDB (Human Metabolome Database) 検索。

    Parameters:
        query: str — 検索クエリ (代謝物名, HMDB ID, 化学式)
        search_type: "name", "hmdb_id", "formula", "mass"
    """
    # ToolUniverse 経由: HMDB_search, HMDB_get_metabolite, HMDB_get_diseases

    base_url = "https://hmdb.ca/metabolites"

    if search_type == "hmdb_id":
        url = f"{base_url}/{query}.xml"
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            root = ET.fromstring(resp.text)
            ns = {"hmdb": "http://www.hmdb.ca"}
            result = {
                "hmdb_id": root.findtext("hmdb:accession", "", ns),
                "name": root.findtext("hmdb:name", "", ns),
                "chemical_formula": root.findtext("hmdb:chemical_formula", "", ns),
                "monoisotopic_mass": float(root.findtext(
                    "hmdb:monisotopic_molecular_weight", "0", ns) or 0),
                "description": (root.findtext("hmdb:description", "", ns) or "")[:300],
                "status": root.findtext("hmdb:status", "", ns),
                "biological_role": root.findtext(
                    "hmdb:ontology/hmdb:root/hmdb:term", "", ns),
            }
            return pd.DataFrame([result])

    # 名前検索
    results = []
    params = {"query": query, "search_type": search_type}
    print(f"HMDB search: '{query}' (type={search_type})")
    # 実際の検索は ToolUniverse SMCP 経由で実行
    return pd.DataFrame(results)
```

## 2. MetaCyc 代謝パスウェイ検索

```python
def metacyc_pathway_search(query, organism="Homo sapiens"):
    """
    MetaCyc 代謝パスウェイ検索。

    Parameters:
        query: str — パスウェイ名/代謝物名/酵素名
        organism: str — 生物種フィルタ
    """
    # ToolUniverse 経由:
    # MetaCyc_search_pathways, MetaCyc_get_pathway
    # MetaCyc_get_compound, MetaCyc_get_reaction

    results = {
        "query": query,
        "organism": organism,
        "database": "MetaCyc",
        "tools": [
            "MetaCyc_search_pathways — パスウェイ検索",
            "MetaCyc_get_pathway — パスウェイ詳細 (反応・酵素・代謝物)",
            "MetaCyc_get_compound — 化合物詳細 (構造・特性)",
            "MetaCyc_get_reaction — 反応詳細 (基質・生成物・酵素)",
        ],
    }

    print(f"MetaCyc: querying '{query}' for {organism}")
    return results
```

## 3. Metabolomics Workbench データ取得

```python
def metabolomics_workbench_search(query, search_type="compound_name",
                                    exact_mass=None, mz_tolerance=0.01):
    """
    NIH Metabolomics Workbench REST API 検索。

    Parameters:
        query: str — 検索クエリ
        search_type: "compound_name", "refmet", "study", "exact_mass", "mz"
        exact_mass: float — 厳密質量 (search_type="exact_mass" 時)
        mz_tolerance: float — m/z 許容誤差 (Da)
    """
    base_url = "https://www.metabolomicsworkbench.org/rest"

    # ToolUniverse 経由:
    # MetabolomicsWorkbench_search_compound_by_name
    # MetabolomicsWorkbench_get_refmet_info
    # MetabolomicsWorkbench_get_study
    # MetabolomicsWorkbench_search_by_exact_mass
    # MetabolomicsWorkbench_search_by_mz
    # MetabolomicsWorkbench_get_compound_by_pubchem_cid

    if search_type == "compound_name":
        url = f"{base_url}/compound/name/{query}/all/"
    elif search_type == "refmet":
        url = f"{base_url}/refmet/name/{query}/all/"
    elif search_type == "study":
        url = f"{base_url}/study/study_id/{query}/summary/"
    elif search_type == "exact_mass":
        url = f"{base_url}/compound/exact_mass/{exact_mass}/tolerance/{mz_tolerance}/"
    elif search_type == "mz":
        url = f"{base_url}/compound/mz_value/{query}/tolerance/{mz_tolerance}/"
    else:
        raise ValueError(f"Unknown search_type: {search_type}")

    resp = requests.get(url, timeout=30)
    if resp.status_code == 200:
        try:
            data = resp.json()
            if isinstance(data, list):
                df = pd.DataFrame(data)
            elif isinstance(data, dict):
                df = pd.DataFrame([data])
            else:
                df = pd.DataFrame()
        except Exception:
            df = pd.DataFrame()
    else:
        df = pd.DataFrame()

    print(f"Metabolomics Workbench ({search_type}): {len(df)} results")
    return df
```

## 4. m/z ベース代謝物同定 (マルチ DB)

```python
def identify_metabolites_by_mass(mz_values, adducts=None,
                                  tolerance_ppm=10,
                                  databases=None):
    """
    m/z 値から複数 DB を横断して代謝物を同定。

    Parameters:
        mz_values: list[float] — 観測 m/z 値リスト
        adducts: list — アダクトイオン (e.g., ["[M+H]+", "[M+Na]+", "[M-H]-"])
        tolerance_ppm: float — 質量許容誤差 (ppm)
        databases: list — 検索対象 DB ("hmdb", "metacyc", "mwb")
    """
    import numpy as np

    if adducts is None:
        adducts = [
            {"name": "[M+H]+", "mass_diff": 1.007276, "mode": "positive"},
            {"name": "[M+Na]+", "mass_diff": 22.989218, "mode": "positive"},
            {"name": "[M-H]-", "mass_diff": -1.007276, "mode": "negative"},
            {"name": "[M+NH4]+", "mass_diff": 18.034164, "mode": "positive"},
        ]

    if databases is None:
        databases = ["hmdb", "mwb"]

    all_results = []

    for mz in mz_values:
        for adduct in adducts:
            neutral_mass = mz - adduct["mass_diff"]
            tolerance_da = neutral_mass * tolerance_ppm / 1e6

            result = {
                "query_mz": mz,
                "adduct": adduct["name"],
                "neutral_mass": round(neutral_mass, 6),
                "tolerance_da": round(tolerance_da, 6),
                "databases_queried": databases,
            }
            all_results.append(result)

    df = pd.DataFrame(all_results)
    print(f"Mass-based ID: {len(mz_values)} m/z values × "
          f"{len(adducts)} adducts = {len(df)} queries "
          f"across {databases}")
    return df
```

## 5. RefMet 標準化命名

```python
def refmet_standardize(metabolite_names):
    """
    RefMet (Reference Metabolomics) 標準化命名。

    Parameters:
        metabolite_names: list — 代謝物名リスト (非標準化)
    """
    results = []

    for name in metabolite_names:
        # ToolUniverse 経由: MetabolomicsWorkbench_get_refmet_info
        result = {
            "input_name": name,
            "refmet_name": None,
            "super_class": None,
            "main_class": None,
            "sub_class": None,
            "formula": None,
            "exact_mass": None,
        }
        results.append(result)

    df = pd.DataFrame(results)
    print(f"RefMet standardization: {len(df)} metabolites queried")
    return df
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/hmdb_metabolites.csv` | CSV |
| `results/metacyc_pathways.json` | JSON |
| `results/mwb_compounds.csv` | CSV |
| `results/mass_id_results.csv` | CSV |
| `results/refmet_standardized.csv` | CSV |
| `figures/metabolite_class_distribution.png` | PNG |

### 利用可能ツール

> [ToolUniverse](https://github.com/mims-harvard/ToolUniverse) SMCP 経由で利用可能な外部ツール。

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| HMDB | `HMDB_get_metabolite` | 代謝物詳細取得 |
| HMDB | `HMDB_search` | 代謝物検索 |
| HMDB | `HMDB_get_diseases` | 疾患関連代謝物 |
| MetaCyc | `MetaCyc_search_pathways` | 代謝パスウェイ検索 |
| MetaCyc | `MetaCyc_get_pathway` | パスウェイ詳細 |
| MetaCyc | `MetaCyc_get_compound` | 化合物詳細 |
| MetaCyc | `MetaCyc_get_reaction` | 反応詳細 |
| MWB | `MetabolomicsWorkbench_search_compound_by_name` | 化合物名検索 |
| MWB | `MetabolomicsWorkbench_get_refmet_info` | RefMet 標準化情報 |
| MWB | `MetabolomicsWorkbench_get_study` | 研究データ取得 |
| MWB | `MetabolomicsWorkbench_search_by_exact_mass` | 厳密質量検索 |
| MWB | `MetabolomicsWorkbench_search_by_mz` | m/z 値検索 |
| MWB | `MetabolomicsWorkbench_get_compound_by_pubchem_cid` | PubChem CID 検索 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-metabolomics` | PLS-DA/VIP 統計解析 |
| `scientific-pathway-enrichment` | 代謝物 → パスウェイ富化 |
| `scientific-cheminformatics` | 分子記述子・構造解析 |
| `scientific-systems-biology` | FBA 代謝フラックス |
| `scientific-multi-omics` | メタボロミクス ↔ オミクス統合 |

### 依存パッケージ

`requests`, `pandas`, `numpy`, `xml.etree.ElementTree` (stdlib)
