---
name: scientific-metabolic-modeling
description: |
  代謝モデリングスキル。BiGG Models ゲノムスケール代謝モデル、
  BioModels SBML リポジトリを統合した代謝ネットワーク解析・
  モデル検索パイプライン。
---

# Scientific Metabolic Modeling

BiGG Models と BioModels を活用したゲノムスケール代謝モデルの
検索・探索・解析パイプラインを提供する。

## When to Use

- ゲノムスケール代謝モデル (GEM) を検索・取得するとき
- BiGG Models の反応・代謝物データを調べるとき
- BioModels リポジトリから SBML モデルを取得するとき
- 代謝パスウェイのフラックス解析の準備を行うとき
- 複数生物種の代謝モデルを比較するとき

---

## Quick Start

## 1. BiGG Models 検索

```python
import requests
import pandas as pd

BIGG_API = "http://bigg.ucsd.edu/api/v2"


def bigg_search(query, search_type="models"):
    """
    BiGG Models データベースを検索。

    Parameters:
        query: str — search term
        search_type: str — "models", "reactions", "metabolites"

    ToolUniverse:
        BiGG_search(query=query, search_type=search_type)
        BiGG_list_models()
    """
    url = f"{BIGG_API}/search"
    params = {
        "query": query,
        "search_type": search_type,
    }
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    results = data.get("results", [])
    df = pd.DataFrame(results)
    print(f"BiGG search '{query}' ({search_type}): {len(df)} results")
    return df
```

## 2. BiGG モデル詳細取得

```python
def bigg_get_model(model_id):
    """
    BiGG Models からゲノムスケール代謝モデルの詳細を取得。

    Parameters:
        model_id: str — BiGG model ID (e.g., "iJO1366")

    ToolUniverse:
        BiGG_get_model(model_id=model_id)
        BiGG_get_model_reactions(model_id=model_id)
        BiGG_get_database_version()
    """
    url = f"{BIGG_API}/models/{model_id}"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    info = {
        "model_id": data.get("bigg_id", ""),
        "organism": data.get("organism", ""),
        "genome_name": data.get("genome_name", ""),
        "num_reactions": data.get("reaction_count", 0),
        "num_metabolites": data.get("metabolite_count", 0),
        "num_genes": data.get("gene_count", 0),
    }

    print(f"BiGG model {model_id}: {info['organism']}, "
          f"{info['num_reactions']} reactions, "
          f"{info['num_metabolites']} metabolites, "
          f"{info['num_genes']} genes")
    return info, data
```

## 3. BiGG 反応・代謝物データ

```python
def bigg_get_reaction(reaction_id):
    """
    BiGG 反応の詳細 (反応式, 関連モデル) を取得。

    ToolUniverse:
        BiGG_get_reaction(reaction_id=reaction_id)
    """
    url = f"{BIGG_API}/universal/reactions/{reaction_id}"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    info = {
        "reaction_id": data.get("bigg_id", ""),
        "name": data.get("name", ""),
        "reaction_string": data.get("reaction_string", ""),
        "pseudoreaction": data.get("pseudoreaction", False),
        "model_count": len(data.get("models_containing_reaction", [])),
    }

    print(f"BiGG reaction {reaction_id}: {info['name']}")
    return info, data


def bigg_get_metabolite(metabolite_id):
    """
    BiGG 代謝物の詳細を取得。

    ToolUniverse:
        BiGG_get_metabolite(metabolite_id=metabolite_id)
    """
    url = f"{BIGG_API}/universal/metabolites/{metabolite_id}"
    resp = requests.get(url)
    resp.raise_for_status()
    data = resp.json()

    info = {
        "metabolite_id": data.get("bigg_id", ""),
        "name": data.get("name", ""),
        "formulae": data.get("formulae", []),
        "charges": data.get("charges", []),
        "model_count": len(data.get("models_containing_metabolite", [])),
    }

    print(f"BiGG metabolite {metabolite_id}: {info['name']}")
    return info, data
```

## 4. BioModels リポジトリ検索

```python
BIOMODELS_API = "https://www.ebi.ac.uk/biomodels"


def biomodels_search(query, num_results=10):
    """
    BioModels (SBML) リポジトリからモデルを検索。

    Parameters:
        query: str — search term

    ToolUniverse:
        biomodels_search(query=query)
        BioModels_search_parameters(query=query)
    """
    url = f"{BIOMODELS_API}/search"
    params = {"query": query, "numResults": num_results}
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()

    models = data.get("models", [])
    results = []
    for m in models:
        results.append({
            "model_id": m.get("id", ""),
            "name": m.get("name", ""),
            "format": m.get("format", {}).get("name", ""),
            "submission_date": m.get("submissionDate", ""),
            "publication": m.get("publication", {}).get("title", ""),
        })

    df = pd.DataFrame(results)
    print(f"BioModels '{query}': {data.get('matches', 0)} total, "
          f"{len(df)} returned")
    return df


def biomodels_get_model(model_id):
    """
    BioModels モデル詳細取得。

    ToolUniverse:
        BioModels_get_model(model_id=model_id)
        BioModels_list_files(model_id=model_id)
        BioModels_download_model(model_id=model_id)
    """
    url = f"{BIOMODELS_API}/{model_id}"
    resp = requests.get(url, headers={"Accept": "application/json"})
    resp.raise_for_status()
    data = resp.json()

    info = {
        "model_id": data.get("publicationId", model_id),
        "name": data.get("name", ""),
        "description": data.get("description", ""),
        "format": data.get("format", {}).get("name", ""),
    }

    print(f"BioModels {model_id}: {info['name']}")
    return info, data
```

## 5. 統合代謝モデル探索パイプライン

```python
def metabolic_model_exploration(organism_query):
    """
    BiGG + BioModels を横断した代謝モデル探索。

    ToolUniverse (横断):
        BiGG_search(query=organism_query) → BiGG_get_model(model_id)
        biomodels_search(query=organism_query) → BioModels_get_model(model_id)
    """
    pipeline = {"query": organism_query}

    # Step 1: BiGG search
    bigg_df = bigg_search(organism_query, search_type="models")
    pipeline["bigg_models"] = len(bigg_df)

    if not bigg_df.empty:
        top_model = bigg_df.iloc[0]
        model_id = top_model.get("bigg_id", "")
        if model_id:
            info, _ = bigg_get_model(model_id)
            pipeline["bigg_top_model"] = info

    # Step 2: BioModels search
    bm_df = biomodels_search(organism_query)
    pipeline["biomodels_models"] = len(bm_df)

    print(f"Metabolic models '{organism_query}': "
          f"BiGG={pipeline['bigg_models']}, "
          f"BioModels={pipeline['biomodels_models']}")
    return pipeline
```

## References

### Output Files

| ファイル | 形式 |
|---|---|
| `results/bigg_search.csv` | CSV |
| `results/bigg_model.json` | JSON |
| `results/bigg_reaction.json` | JSON |
| `results/biomodels_search.csv` | CSV |
| `results/biomodels_model.json` | JSON |

### 利用可能ツール

| カテゴリ | 主要ツール | 用途 |
|---|---|---|
| BiGG | `BiGG_search` | モデル/反応/代謝物検索 |
| BiGG | `BiGG_list_models` | モデル一覧 |
| BiGG | `BiGG_get_model` | モデル詳細 |
| BiGG | `BiGG_get_model_reactions` | モデル反応一覧 |
| BiGG | `BiGG_get_reaction` | 反応詳細 |
| BiGG | `BiGG_get_metabolite` | 代謝物詳細 |
| BiGG | `BiGG_get_database_version` | DB バージョン |
| BioModels | `biomodels_search` | モデル検索 |
| BioModels | `BioModels_get_model` | モデル詳細 |
| BioModels | `BioModels_list_files` | ファイル一覧 |
| BioModels | `BioModels_download_model` | モデル DL |
| BioModels | `BioModels_search_parameters` | パラメータ検索 |

### 参照スキル

| スキル | 関連 |
|---|---|
| `scientific-pathway-enrichment` | パスウェイ解析 |
| `scientific-systems-biology` | システム生物学 |
| `scientific-gene-expression-transcriptomics` | 発現データ |
| `scientific-biothings-idmapping` | ID マッピング |

### 依存パッケージ

`requests`, `pandas`
