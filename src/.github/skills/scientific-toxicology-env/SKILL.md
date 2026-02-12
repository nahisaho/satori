---
name: scientific-toxicology-env
description: |
  毒性学・環境衛生スキル。CTD (Comparative Toxicogenomics Database)
  化学-遺伝子-疾患関連・ToxCast/Tox21 高スループット毒性スクリーニング・
  IRIS ヒトリスク評価・T3DB 食品/環境毒性物質・PubChem BioAssay 毒性データ。
---

# Scientific Toxicology & Environmental Health

CTD / Tox21 / ToxCast / T3DB / IRIS を活用した毒性学・環境衛生
パイプラインを提供する。化学物質-遺伝子-疾患関連、高スループット
毒性スクリーニング、ヒト健康リスク評価。

## When to Use

- 化学物質が関連する遺伝子・疾患を調べるとき (CTD)
- Tox21/ToxCast アッセイデータで毒性メカニズムを解析するとき
- IRIS リスク評価データ (RfD/RfC/UR) を参照するとき
- 食品・環境毒性物質の詳細データベース (T3DB) を検索するとき
- PubChem BioAssay から毒性関連ハイスループットデータを取得するとき
- ADMET 毒性予測と実験毒性データを併用するとき

---

## Quick Start

## 1. CTD 化学-遺伝子-疾患関連検索

```python
import requests
import pandas as pd
import json

CTD_BASE = "https://ctdbase.org/tools"


def ctd_chemical_gene(chemical_name, limit=100):
    """
    CTD — 化学物質-遺伝子相互作用を検索。

    Parameters:
        chemical_name: str — 化学物質名 (例: "Bisphenol A")
        limit: int — 最大取得数
    """
    url = f"{CTD_BASE}/batchQuery.go"
    params = {
        "inputType": "chem",
        "inputTerms": chemical_name,
        "report": "genes_curated",
        "format": "json",
    }
    resp = requests.get(url, params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()[:limit]

    results = []
    for entry in data:
        results.append({
            "chemical": entry.get("ChemicalName", ""),
            "gene": entry.get("GeneSymbol", ""),
            "organism": entry.get("Organism", ""),
            "interaction": entry.get("Interaction", ""),
            "pubmed_ids": entry.get("PubMedIDs", ""),
        })

    df = pd.DataFrame(results)
    print(f"CTD: {chemical_name} → {len(df)} gene interactions")
    return df
```

## 2. CTD 化学-疾患関連検索

```python
def ctd_chemical_disease(chemical_name, limit=100):
    """
    CTD — 化学物質-疾患関連を検索。

    Parameters:
        chemical_name: str — 化学物質名
        limit: int — 最大取得数
    """
    url = f"{CTD_BASE}/batchQuery.go"
    params = {
        "inputType": "chem",
        "inputTerms": chemical_name,
        "report": "diseases_curated",
        "format": "json",
    }
    resp = requests.get(url, params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()[:limit]

    results = []
    for entry in data:
        results.append({
            "chemical": entry.get("ChemicalName", ""),
            "disease": entry.get("DiseaseName", ""),
            "disease_id": entry.get("DiseaseID", ""),
            "direct_evidence": entry.get("DirectEvidence", ""),
            "inference_score": float(entry.get("InferenceScore", 0)),
        })

    df = pd.DataFrame(results)
    df = df.sort_values("inference_score", ascending=False)
    print(f"CTD: {chemical_name} → {len(df)} disease associations")
    return df
```

## 3. Tox21/ToxCast アッセイデータ取得

```python
COMPTOX_BASE = "https://comptox.epa.gov/dashboard/api"


def tox21_assay_search(chemical_identifier, assay_source="Tox21"):
    """
    CompTox Dashboard — Tox21/ToxCast アッセイ結果を取得。

    Parameters:
        chemical_identifier: str — DTXSID or CAS
        assay_source: str — "Tox21" or "ToxCast"
    """
    # CompTox Dashboard API で化学物質のアッセイデータ取得
    url = f"{COMPTOX_BASE}/chemical/search"
    params = {"query": chemical_identifier}
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    chem_data = resp.json()

    dtxsid = chem_data.get("dtxsid", chemical_identifier)

    # アッセイエンドポイント取得
    url_assay = f"{COMPTOX_BASE}/chemical/{dtxsid}/assays"
    resp_assay = requests.get(url_assay, timeout=30)
    resp_assay.raise_for_status()
    assays = resp_assay.json()

    results = []
    for assay in assays:
        if assay_source.lower() in assay.get("assaySource", "").lower():
            results.append({
                "assay_name": assay.get("assayName", ""),
                "assay_source": assay.get("assaySource", ""),
                "endpoint": assay.get("assayEndpoint", ""),
                "activity": assay.get("activity", ""),
                "ac50_um": assay.get("ac50", None),
                "hit_call": assay.get("hitCall", ""),
            })

    df = pd.DataFrame(results)
    n_active = (df["hit_call"] == "Active").sum() if len(df) > 0 else 0
    print(f"Tox21/ToxCast: {dtxsid} → {len(df)} assays, {n_active} active")
    return df
```

## 4. T3DB 毒性物質検索

```python
T3DB_BASE = "https://t3db.ca/api"


def t3db_search(query, search_type="name"):
    """
    T3DB — 食品/環境毒性物質データベース検索。

    Parameters:
        query: str — 検索語
        search_type: str — "name", "cas", "category"
    """
    url = f"{T3DB_BASE}/toxins/search"
    params = {"query": query, "search_type": search_type}
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    results = []
    for toxin in data.get("toxins", []):
        results.append({
            "name": toxin.get("name", ""),
            "t3db_id": toxin.get("t3db_id", ""),
            "cas_number": toxin.get("cas_number", ""),
            "category": toxin.get("category", ""),
            "toxicity_class": toxin.get("toxicity_class", ""),
            "ld50_oral": toxin.get("ld50_oral", ""),
            "target_organs": toxin.get("target_organs", []),
        })

    df = pd.DataFrame(results)
    print(f"T3DB: '{query}' → {len(df)} toxins")
    return df
```

## 5. EPA IRIS リスク評価

```python
def iris_risk_assessment(chemical_name):
    """
    EPA IRIS — ヒト健康リスク評価データ取得。

    Parameters:
        chemical_name: str — 化学物質名
    """
    url = "https://iris.epa.gov/AtoZ"
    resp = requests.get(url, timeout=30)

    # IRIS は構造化 API なし — スクレイピングまたはローカルデータ
    # 代替: CompTox Dashboard API 経由
    url_comptox = f"{COMPTOX_BASE}/chemical/search"
    params = {"query": chemical_name}
    resp = requests.get(url_comptox, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    risk_data = {
        "chemical": chemical_name,
        "dtxsid": data.get("dtxsid", ""),
        "rfd_oral_mg_kg_day": data.get("rfdOral", None),
        "rfc_inhalation_mg_m3": data.get("rfcInhalation", None),
        "cancer_classification": data.get("cancerClassification", ""),
        "oral_slope_factor": data.get("oralSlopeFactor", None),
        "inhalation_unit_risk": data.get("inhalationUnitRisk", None),
    }

    print(f"IRIS: {chemical_name}")
    for k, v in risk_data.items():
        if v and k != "chemical":
            print(f"  {k}: {v}")
    return risk_data
```

## 6. 毒性パスウェイ解析

```python
def toxicity_pathway_analysis(chemical_name, species="Homo sapiens"):
    """
    CTD + パスウェイ統合毒性解析。

    Parameters:
        chemical_name: str — 化学物質名
        species: str — 生物種
    """
    # 1) CTD 遺伝子取得
    gene_df = ctd_chemical_gene(chemical_name, limit=500)
    if species:
        gene_df = gene_df[gene_df["organism"] == species]

    gene_list = gene_df["gene"].unique().tolist()

    # 2) CTD 疾患取得
    disease_df = ctd_chemical_disease(chemical_name, limit=100)

    # 3) パスウェイ濃縮 (KEGG enrichment via Enrichr)
    enrichr_url = "https://maayanlab.cloud/Enrichr"
    add_resp = requests.post(
        f"{enrichr_url}/addList",
        files={"list": (None, "\n".join(gene_list))}
    )
    user_list_id = add_resp.json()["userListId"]

    enrich_resp = requests.get(
        f"{enrichr_url}/enrich",
        params={"userListId": user_list_id, "backgroundType": "KEGG_2021_Human"}
    )
    pathways = enrich_resp.json().get("KEGG_2021_Human", [])

    pathway_results = []
    for pw in pathways[:20]:
        pathway_results.append({
            "pathway": pw[1],
            "p_value": pw[2],
            "adj_p_value": pw[6],
            "genes": pw[5],
        })

    print(f"Toxicity pathway: {chemical_name}")
    print(f"  Target genes: {len(gene_list)}")
    print(f"  Diseases: {len(disease_df)}")
    print(f"  Pathways: {len(pathway_results)}")
    return {
        "genes": gene_df,
        "diseases": disease_df,
        "pathways": pd.DataFrame(pathway_results),
    }
```

---

## パイプライン統合

```
admet-pharmacokinetics → toxicology-env → pharmacovigilance
  (ADMET 毒性予測)       (CTD/Tox21/IRIS)   (市販後安全性)
        │                       │                  ↓
cheminformatics ───────────────┘           disease-research
  (RDKit 構造アラート)    │                (疾患-遺伝子関連)
                          ↓
                   public-health-data
                   (CDC/WHO 公衆衛生)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/ctd_gene_interactions.csv` | CTD 化学-遺伝子関連 | → pathway-enrichment |
| `results/ctd_disease_associations.csv` | CTD 化学-疾患関連 | → disease-research |
| `results/tox21_assays.csv` | Tox21/ToxCast アッセイ結果 | → admet-pharmacokinetics |
| `results/toxicity_pathways.json` | 毒性パスウェイ解析結果 | → pharmacovigilance |
