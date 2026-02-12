---
name: scientific-public-health-data
description: |
  公衆衛生データアクセススキル。NHANES 疫学調査データ、MedlinePlus 一般向け
  健康情報、RxNorm 薬剤標準語彙、ODPHP 健康目標・ガイドライン、
  Health Disparities 健康格差データ統合パイプライン。
---

# Scientific Public Health Data

NHANES / MedlinePlus / RxNorm / ODPHP / Health Disparities /
Guidelines を統合した公衆衛生データアクセスパイプラインを提供する。

## When to Use

- NHANES 疫学調査データ (検査値・アンケート) を取得するとき
- MedlinePlus で一般向け健康情報を検索するとき
- RxNorm で薬剤名の標準化・マッピングを行うとき
- ODPHP Healthy People 目標や健康ガイドラインを参照するとき
- 健康格差 (Health Disparities) データを分析するとき
- 臨床ガイドライン (USPSTF/WHO) を検索するとき

---

## Quick Start

## 1. NHANES 疫学調査データ取得

```python
import requests
import pandas as pd
import io

NHANES_BASE = "https://wwwn.cdc.gov/nchs/nhanes"


def get_nhanes_dataset(cycle, dataset_name):
    """
    NHANES データセット (XPT/SAS 形式) 取得。

    Parameters:
        cycle: str — 調査サイクル (e.g., "2017-2018", "2019-2020")
        dataset_name: str — データセット名 (e.g., "DEMO_J", "BIOPRO_J")

    ToolUniverse:
        NHANES_get_dataset(cycle=cycle, dataset=dataset_name)
        NHANES_list_datasets(cycle=cycle)
    """
    cycle_code = cycle.replace("-", "_")
    url = f"{NHANES_BASE}/search/DataPage.aspx"

    # XPT ファイルの直接ダウンロード
    xpt_url = f"https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/{dataset_name}.XPT"
    resp = requests.get(xpt_url)
    resp.raise_for_status()

    df = pd.read_sas(io.BytesIO(resp.content), format="xport")
    print(f"NHANES {cycle} {dataset_name}: {df.shape[0]} rows × {df.shape[1]} columns")
    return df


def search_nhanes_variables(keyword):
    """
    NHANES 変数検索。

    Parameters:
        keyword: str — 変数名/説明の検索語

    ToolUniverse:
        NHANES_search_variables(keyword=keyword)
    """
    url = f"{NHANES_BASE}/search/variablelist.aspx"
    params = {"SearchTarget": keyword}
    resp = requests.get(url, params=params)
    resp.raise_for_status()

    print(f"NHANES variable search '{keyword}': response received")
    return resp.text
```

## 2. MedlinePlus 健康情報検索

```python
MEDLINEPLUS_API = "https://connect.medlineplus.gov/service"
MEDLINEPLUS_WS = "https://wsearch.nlm.nih.gov/ws/query"


def search_medlineplus_health_topics(query, language="English"):
    """
    MedlinePlus 健康トピック検索。

    ToolUniverse:
        MedlinePlus_search_health_topics(query=query)
        MedlinePlus_get_health_topic(topic_id=topic_id)
        MedlinePlus_search_drugs(query=query)
        MedlinePlus_search_labs(query=query)
        MedlinePlus_connect(code=code, code_system=system)
    """
    params = {
        "db": "healthTopics",
        "term": query,
    }
    resp = requests.get(MEDLINEPLUS_WS, params=params)
    resp.raise_for_status()

    # XML response parsing
    import xml.etree.ElementTree as ET
    root = ET.fromstring(resp.text)

    results = []
    for doc in root.findall(".//document"):
        results.append({
            "title": doc.find(".//content[@name='title']").text
            if doc.find(".//content[@name='title']") is not None else "",
            "url": doc.get("url", ""),
            "summary": doc.find(".//content[@name='FullSummary']").text[:300]
            if doc.find(".//content[@name='FullSummary']") is not None else "",
            "rank": doc.get("rank", ""),
        })

    df = pd.DataFrame(results)
    print(f"MedlinePlus search '{query}': {len(df)} health topics")
    return df
```

## 3. RxNorm 薬剤標準語彙

```python
RXNORM_API = "https://rxnav.nlm.nih.gov/REST"


def rxnorm_lookup(drug_name):
    """
    RxNorm 薬剤名正規化・コードマッピング。

    Parameters:
        drug_name: str — 薬剤名 (商品名 or 一般名)

    ToolUniverse:
        RxNorm_get_rxcui(name=drug_name)
    """
    resp = requests.get(
        f"{RXNORM_API}/rxcui.json",
        params={"name": drug_name}
    )
    resp.raise_for_status()
    data = resp.json()

    rxcui = data.get("idGroup", {}).get("rxnormId", [None])[0]
    if not rxcui:
        print(f"RxNorm: '{drug_name}' not found")
        return None

    # Get properties
    props_resp = requests.get(f"{RXNORM_API}/rxcui/{rxcui}/properties.json")
    props_resp.raise_for_status()
    props = props_resp.json().get("properties", {})

    # Get related concepts
    related_resp = requests.get(
        f"{RXNORM_API}/rxcui/{rxcui}/related.json",
        params={"tty": "IN+BN+SBD+SCD"}
    )
    related_resp.raise_for_status()
    related = related_resp.json()

    result = {
        "rxcui": rxcui,
        "name": props.get("name", ""),
        "tty": props.get("tty", ""),
        "synonym": props.get("synonym", ""),
        "related_concepts": [
            {
                "rxcui": c.get("rxcui"),
                "name": c.get("name"),
                "tty": c.get("tty"),
            }
            for group in related.get("relatedGroup", {}).get("conceptGroup", [])
            for c in group.get("conceptProperties", [])
        ],
    }
    print(f"RxNorm '{drug_name}': RXCUI={rxcui}, TTY={result['tty']}")
    return result
```

## 4. Health Disparities データ取得

```python
HD_API = "https://data.cdc.gov/resource"


def get_health_disparities(indicator, dataset_id="pqnx-3xr5"):
    """
    CDC 健康格差データ取得。

    Parameters:
        indicator: str — 健康指標名
        dataset_id: str — CDC Socrata データセット ID

    ToolUniverse:
        HealthDisparities_search(query=indicator)
        HealthDisparities_get_indicators(category=category)
    """
    params = {
        "$where": f"indicator LIKE '%{indicator}%'",
        "$limit": 1000,
    }
    resp = requests.get(f"{HD_API}/{dataset_id}.json", params=params)
    resp.raise_for_status()
    data = resp.json()

    df = pd.DataFrame(data)
    print(f"Health Disparities '{indicator}': {len(df)} records")
    return df
```

## 5. ODPHP 健康ガイドライン

```python
ODPHP_API = "https://health.gov/myhealthfinder/api/v3"


def search_health_guidelines(keyword, category=None):
    """
    ODPHP MyHealthfinder ガイドライン検索。

    ToolUniverse:
        ODPHP_search_topics(keyword=keyword)
        ODPHP_get_topic(topic_id=topic_id)
    """
    params = {"keyword": keyword}
    if category:
        params["categoryId"] = category
    resp = requests.get(f"{ODPHP_API}/topicsearch.json", params=params)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for topic in data.get("Result", {}).get("Resources", {}).get("Resource", []):
        results.append({
            "title": topic.get("Title", ""),
            "categories": topic.get("Categories", ""),
            "url": topic.get("AccessibleVersion", ""),
            "sections": [
                s.get("Title", "") for s in topic.get("Sections", {}).get("section", [])
            ],
        })

    df = pd.DataFrame(results)
    print(f"ODPHP search '{keyword}': {len(df)} guidelines")
    return df
```

## 6. 臨床ガイドライン検索 (USPSTF)

```python
def search_clinical_guidelines(query, source="uspstf"):
    """
    USPSTF/WHO 臨床ガイドライン検索。

    ToolUniverse:
        Guidelines_search(query=query, source=source)
        Guidelines_get_recommendations(topic_id=topic_id)
    """
    sources = {
        "uspstf": "https://www.uspreventiveservicestaskforce.org/uspstf/api",
        "who": "https://app.magicapp.org/api",
    }
    base_url = sources.get(source, sources["uspstf"])

    resp = requests.get(f"{base_url}/search", params={"q": query})
    if resp.status_code == 200:
        data = resp.json()
        results = []
        for item in data.get("results", []):
            results.append({
                "title": item.get("title", ""),
                "grade": item.get("grade", ""),
                "population": item.get("population", ""),
                "date": item.get("date", ""),
                "recommendation": item.get("recommendation", ""),
            })
        df = pd.DataFrame(results)
    else:
        df = pd.DataFrame()

    print(f"Guidelines ({source}) search '{query}': {len(df)} recommendations")
    return df
```

---

## 利用可能ツール

| ToolUniverse カテゴリ | 主なツール |
|---|---|
| `nhanes` | `NHANES_get_dataset`, `NHANES_list_datasets`, `NHANES_search_variables` |
| `health_disparities` | `HealthDisparities_search`, `HealthDisparities_get_indicators` |
| `medlineplus` | `MedlinePlus_search_health_topics`, `MedlinePlus_get_health_topic`, `MedlinePlus_search_drugs`, `MedlinePlus_search_labs`, `MedlinePlus_connect` |
| `odphp` | `ODPHP_search_topics`, `ODPHP_get_topic` |
| `rxnorm` | `RxNorm_get_rxcui` |
| `guidelines_tools` | `Guidelines_search`, `Guidelines_get_recommendations` |

## パイプライン出力

| 出力ファイル | 説明 | 連携先スキル |
|---|---|---|
| `results/nhanes_data.csv` | NHANES 疫学データ | → epidemiology-public-health, survival-clinical |
| `results/drug_mapping.json` | RxNorm 薬剤マッピング | → pharmacovigilance, pharmacogenomics |
| `results/health_guidelines.json` | 臨床ガイドライン | → clinical-decision-support |
| `results/health_disparities.csv` | 健康格差指標 | → epidemiology-public-health, causal-inference |

## パイプライン統合

```
epidemiology-public-health ──→ public-health-data ──→ clinical-decision-support
  (RR/OR/DAG)                  (NHANES/CDC/ODPHP)    (GRADE エビデンス)
                                      │
                                      ├──→ pharmacovigilance (RxNorm + 安全性)
                                      ├──→ pharmacogenomics (RxNorm + PGx)
                                      └──→ survival-clinical (NHANES コホート)
```
