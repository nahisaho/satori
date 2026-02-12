---
name: scientific-opentargets-genetics
description: |
  Open Targets Platform 遺伝学スキル。Open Targets Platform
  GraphQL API を用いた標的-疾患アソシエーション・薬剤
  エビデンス・L2G 遺伝的関連・ファーマコゲノミクス検索。
  ToolUniverse 連携: opentarget。
tu_tools:
  - key: opentarget
    name: Open Targets
    description: 標的-疾患アソシエーション GraphQL API
---

# Scientific Open Targets Genetics

Open Targets Platform GraphQL API を活用した標的-疾患
アソシエーションスコア取得・薬剤エビデンス検索・L2G
遺伝的関連パイプラインを提供する。

## When to Use

- 遺伝子 (標的) と疾患のアソシエーションスコアを検索するとき
- 薬剤エビデンスデータを取得するとき
- GWAS バリアントから遺伝子を L2G スコアでマッピングするとき
- 標的の安全性プロファイルを確認するとき
- ファーマコゲノミクスデータを検索するとき

---

## Quick Start

## 1. 標的-疾患アソシエーション

```python
import requests
import pandas as pd

OT_API = ("https://api.platform.opentargets.org"
          "/api/v4/graphql")


def ot_target_disease_assoc(target_id, limit=25):
    """
    Open Targets — 標的-疾患アソシエーション。

    Parameters:
        target_id: str — Ensembl Gene ID
            (例: "ENSG00000012048" = BRCA1)
        limit: int — 最大結果数
    """
    query = """
    query targetDisease($id: String!, $size: Int!) {
      target(ensemblId: $id) {
        id
        approvedSymbol
        associatedDiseases(page: {size: $size, index: 0}) {
          count
          rows {
            disease { id name }
            score
            datatypeScores {
              componentId: id
              score
            }
          }
        }
      }
    }
    """
    variables = {"id": target_id, "size": limit}
    resp = requests.post(OT_API,
                         json={"query": query,
                               "variables": variables},
                         timeout=30)
    resp.raise_for_status()
    data = resp.json()["data"]["target"]

    rows = []
    for r in data["associatedDiseases"]["rows"]:
        row = {
            "target_id": target_id,
            "target_symbol": data["approvedSymbol"],
            "disease_id": r["disease"]["id"],
            "disease_name": r["disease"]["name"],
            "overall_score": r["score"],
        }
        for dt in r["datatypeScores"]:
            row[dt["componentId"]] = dt["score"]
        rows.append(row)

    df = pd.DataFrame(rows)
    total = data["associatedDiseases"]["count"]
    print(f"OT associations: {data['approvedSymbol']} "
          f"→ {len(df)}/{total} diseases")
    return df
```

## 2. 薬剤エビデンス

```python
def ot_drug_evidence(target_id, disease_id, limit=50):
    """
    Open Targets — 薬剤エビデンス。

    Parameters:
        target_id: str — Ensembl Gene ID
        disease_id: str — EFO Disease ID
            (例: "EFO_0000305" = breast carcinoma)
        limit: int — 最大結果数
    """
    query = """
    query drugEvidence($ensemblId: String!,
                       $efoId: String!,
                       $size: Int!) {
      disease(efoId: $efoId) {
        id
        name
        evidences(
          ensemblIds: [$ensemblId]
          datasourceIds: ["chembl"]
          size: $size
        ) {
          count
          rows {
            id
            score
            drug {
              id name drugType
              maximumClinicalTrialPhase
              mechanismsOfAction {
                rows { actionType }
              }
            }
            clinicalPhase
            clinicalStatus
            urls { niceName url }
          }
        }
      }
    }
    """
    variables = {"ensemblId": target_id,
                 "efoId": disease_id,
                 "size": limit}
    resp = requests.post(OT_API,
                         json={"query": query,
                               "variables": variables},
                         timeout=30)
    resp.raise_for_status()
    data = resp.json()["data"]["disease"]

    results = []
    for ev in data["evidences"]["rows"]:
        drug = ev.get("drug", {})
        moas = drug.get("mechanismsOfAction", {})
        moa_list = [m["actionType"]
                    for m in moas.get("rows", [])]
        results.append({
            "disease": data["name"],
            "drug_id": drug.get("id", ""),
            "drug_name": drug.get("name", ""),
            "drug_type": drug.get("drugType", ""),
            "max_phase": drug.get(
                "maximumClinicalTrialPhase", 0),
            "clinical_phase": ev.get("clinicalPhase", ""),
            "clinical_status": ev.get(
                "clinicalStatus", ""),
            "moa": "; ".join(moa_list),
            "score": ev.get("score", 0),
        })

    df = pd.DataFrame(results)
    print(f"OT drug evidence: {len(df)} entries")
    return df
```

## 3. L2G 遺伝的関連 (Locus-to-Gene)

```python
def ot_l2g_variants(study_id, limit=50):
    """
    Open Targets Genetics — L2G バリアント-遺伝子マッピング。

    Parameters:
        study_id: str — GWAS Study ID
            (例: "GCST004988")
        limit: int — 最大結果数
    """
    # OT Genetics API
    OT_GENETICS = ("https://api.genetics.opentargets.org"
                   "/graphql")
    query = """
    query l2g($studyId: String!, $size: Int!) {
      studyLocus2GeneTable(studyId: $studyId,
                           pageSize: $size) {
        rows {
          gene { id symbol }
          variant { id rsId }
          yProbaModel
          yProbaDistance
          yProbaInteraction
          yProbaMolecularQTL
          yProbaPathogenicity
          hasColoc
          distanceToLocus
        }
      }
    }
    """
    variables = {"studyId": study_id, "size": limit}
    resp = requests.post(OT_GENETICS,
                         json={"query": query,
                               "variables": variables},
                         timeout=30)
    resp.raise_for_status()
    data = resp.json()["data"]["studyLocus2GeneTable"]

    rows = []
    for r in data["rows"]:
        rows.append({
            "gene_id": r["gene"]["id"],
            "gene_symbol": r["gene"]["symbol"],
            "variant_id": r["variant"]["id"],
            "rsid": r["variant"]["rsId"],
            "l2g_score": r["yProbaModel"],
            "distance_score": r["yProbaDistance"],
            "interaction_score": r["yProbaInteraction"],
            "qtl_score": r["yProbaMolecularQTL"],
            "pathogenicity": r["yProbaPathogenicity"],
            "has_coloc": r["hasColoc"],
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("l2g_score", ascending=False)
    print(f"OT L2G: {study_id} → {len(df)} gene mappings")
    return df
```

## 4. Open Targets 統合パイプライン

```python
def ot_pipeline(gene_symbol, ensembl_id,
                   output_dir="results"):
    """
    Open Targets 統合パイプライン。

    Parameters:
        gene_symbol: str — 遺伝子シンボル (例: "BRCA1")
        ensembl_id: str — Ensembl Gene ID
        output_dir: str — 出力ディレクトリ
    """
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1) 標的-疾患アソシエーション
    assoc = ot_target_disease_assoc(ensembl_id)
    assoc.to_csv(output_dir / "ot_associations.csv",
                 index=False)

    # 2) トップ疾患の薬剤エビデンス
    if not assoc.empty:
        top_disease = assoc.iloc[0]["disease_id"]
        drugs = ot_drug_evidence(ensembl_id, top_disease)
        drugs.to_csv(output_dir / "ot_drugs.csv",
                     index=False)

    print(f"OT pipeline: {gene_symbol} → {output_dir}")
    return {"associations": assoc}
```

---

## ToolUniverse 連携

| TU Key | ツール名 | 連携内容 |
|--------|---------|---------|
| `opentarget` | Open Targets | 標的-疾患アソシエーション GraphQL (~55 tools) |

## パイプライン統合

```
disease-research → opentargets-genetics → drug-target-profiling
  (疾患遺伝子)      (OT Platform API)     (標的プロファイリング)
        │                   │                    ↓
variant-interpretation ────┘          pharmacogenomics
  (ClinVar/VEP)      │              (薬理ゲノミクス)
                      ↓
            gnomad-variants
            (集団頻度)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/ot_associations.csv` | 標的-疾患スコア | → disease-research |
| `results/ot_drugs.csv` | 薬剤エビデンス | → drug-target-profiling |
