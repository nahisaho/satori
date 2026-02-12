---
name: scientific-chembl-assay-mining
description: |
  ChEMBL アッセイ・活性データマイニングスキル。ChEMBL REST API による
  アッセイ検索・バイオアクティビティデータ取得・IC50/Ki/EC50 SAR 解析・
  ターゲット-化合物マッピング・選択性プロファイリング・ATC 分類検索・
  構造アラート検出パイプライン。
---

# Scientific ChEMBL Assay Mining

ChEMBL REST API (EBI) を活用したバイオアクティビティデータマイニング
パイプラインを提供する。アッセイ検索、活性値解析、SAR (構造活性相関)、
ターゲット選択性プロファイリング、分子ドッキング前処理を統合。

## When to Use

- ChEMBL からターゲット特異的アッセイデータを取得するとき
- IC50/Ki/EC50 をバルク取得して SAR 解析するとき
- 特定標的に対する化合物選択性を評価するとき
- 分子構造類似性検索・サブ構造検索を行うとき
- ATC 分類から承認薬を特定し薬理プロファイルを構築するとき
- 構造アラート (PAINS, Dundee) をスクリーニングするとき

---

## Quick Start

## 1. ターゲット検索 & アッセイ取得

```python
import requests
import pandas as pd

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
HEADERS = {"Accept": "application/json"}


def search_target(query, organism="Homo sapiens", limit=10):
    """
    ChEMBL ターゲット検索。

    Parameters:
        query: str — ターゲット名 (例: "EGFR", "CDK4")
        organism: str — 生物種
        limit: int — 最大取得数

    ToolUniverse:
        ChEMBL_search_targets(pref_name__contains=query, organism=organism)
        ChEMBL_get_target(target_chembl_id=target_id)
    """
    url = f"{CHEMBL_API}/target.json"
    params = {
        "pref_name__icontains": query,
        "organism": organism,
        "limit": limit,
    }
    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    targets = resp.json().get("targets", [])

    rows = []
    for t in targets:
        rows.append({
            "target_chembl_id": t.get("target_chembl_id"),
            "pref_name": t.get("pref_name"),
            "target_type": t.get("target_type"),
            "organism": t.get("organism"),
        })

    df = pd.DataFrame(rows)
    print(f"ChEMBL targets matching '{query}': {len(df)}")
    return df
```

## 2. バイオアクティビティデータ取得

```python
def get_target_activities(target_chembl_id, standard_type="IC50",
                          max_value=10000, limit=500):
    """
    ターゲットに対するバイオアクティビティデータ取得。

    Parameters:
        target_chembl_id: str — ChEMBL ターゲットID
        standard_type: str — "IC50", "Ki", "EC50", "Kd" etc.
        max_value: float — nM 単位閾値
        limit: int — 最大取得数

    ToolUniverse:
        ChEMBL_search_activities(
            target_chembl_id=target_chembl_id,
            standard_type=standard_type,
            standard_value__lte=max_value
        )
        ChEMBL_get_target_activities(target_chembl_id__exact=target_chembl_id)
    """
    url = f"{CHEMBL_API}/activity.json"
    params = {
        "target_chembl_id": target_chembl_id,
        "standard_type": standard_type,
        "standard_value__lte": max_value,
        "standard_units": "nM",
        "limit": limit,
    }
    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    activities = resp.json().get("activities", [])

    rows = []
    for act in activities:
        rows.append({
            "molecule_chembl_id": act.get("molecule_chembl_id"),
            "canonical_smiles": act.get("canonical_smiles"),
            "standard_type": act.get("standard_type"),
            "standard_value": act.get("standard_value"),
            "standard_units": act.get("standard_units"),
            "pchembl_value": act.get("pchembl_value"),
            "assay_chembl_id": act.get("assay_chembl_id"),
            "assay_type": act.get("assay_type"),
            "target_chembl_id": act.get("target_chembl_id"),
        })

    df = pd.DataFrame(rows)
    if "standard_value" in df.columns:
        df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    print(f"Activities for {target_chembl_id} ({standard_type}): {len(df)}")
    return df
```

## 3. アッセイ詳細検索

```python
def search_assays(target_chembl_id=None, assay_type=None, limit=50):
    """
    ChEMBL アッセイ検索。

    Parameters:
        target_chembl_id: str — ターゲット ChEMBL ID
        assay_type: str — "B" (Binding), "F" (Functional), "A" (ADME)
        limit: int — 最大取得数

    ToolUniverse:
        ChEMBL_search_assays(
            target_chembl_id=target_chembl_id,
            assay_type=assay_type
        )
        ChEMBL_get_assay(assay_chembl_id=assay_id)
        ChEMBL_get_assay_activities(assay_chembl_id__exact=assay_id)
    """
    url = f"{CHEMBL_API}/assay.json"
    params = {"limit": limit}
    if target_chembl_id:
        params["target_chembl_id"] = target_chembl_id
    if assay_type:
        params["assay_type"] = assay_type

    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    assays = resp.json().get("assays", [])

    rows = []
    for a in assays:
        rows.append({
            "assay_chembl_id": a.get("assay_chembl_id"),
            "description": a.get("description", "")[:200],
            "assay_type": a.get("assay_type"),
            "assay_organism": a.get("assay_organism"),
            "confidence_score": a.get("confidence_score"),
            "target_chembl_id": a.get("target_chembl_id"),
        })

    df = pd.DataFrame(rows)
    print(f"Assays found: {len(df)}")
    return df
```

## 4. SAR (構造活性相関) 解析

```python
import numpy as np


def sar_analysis(activity_df, pchembl_col="pchembl_value"):
    """
    バイオアクティビティデータの SAR 解析。

    Parameters:
        activity_df: DataFrame — get_target_activities の出力
        pchembl_col: str — pChEMBL 値カラム名

    Returns:
        dict — SAR 統計サマリ
    """
    df = activity_df.copy()
    df[pchembl_col] = pd.to_numeric(df[pchembl_col], errors="coerce")
    df = df.dropna(subset=[pchembl_col])

    summary = {
        "n_compounds": len(df),
        "n_unique_molecules": df["molecule_chembl_id"].nunique(),
        "pchembl_mean": round(df[pchembl_col].mean(), 2),
        "pchembl_median": round(df[pchembl_col].median(), 2),
        "pchembl_std": round(df[pchembl_col].std(), 2),
        "pchembl_range": [
            round(df[pchembl_col].min(), 2),
            round(df[pchembl_col].max(), 2),
        ],
        "most_potent": df.loc[df[pchembl_col].idxmax()].to_dict()
            if len(df) > 0 else None,
        "potency_bins": {
            "high_potent_lt100nM": int((df[pchembl_col] >= 7.0).sum()),
            "moderate_100_1000nM": int(
                ((df[pchembl_col] >= 6.0) & (df[pchembl_col] < 7.0)).sum()
            ),
            "weak_gt1000nM": int((df[pchembl_col] < 6.0).sum()),
        },
    }

    print(f"SAR summary: {summary['n_unique_molecules']} molecules, "
          f"pChEMBL mean={summary['pchembl_mean']}")
    return summary
```

## 5. 選択性プロファイリング

```python
def selectivity_profile(molecule_chembl_id, limit=100):
    """
    化合物のマルチターゲット選択性評価。

    Parameters:
        molecule_chembl_id: str — 化合物 ChEMBL ID
        limit: int — 最大取得数

    ToolUniverse:
        ChEMBL_get_molecule_targets(
            molecule_chembl_id__exact=molecule_chembl_id
        )
        ChEMBL_search_activities(molecule_chembl_id=molecule_chembl_id)
    """
    url = f"{CHEMBL_API}/activity.json"
    params = {
        "molecule_chembl_id": molecule_chembl_id,
        "limit": limit,
    }
    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    activities = resp.json().get("activities", [])

    target_data = {}
    for act in activities:
        tid = act.get("target_chembl_id")
        pchembl = act.get("pchembl_value")
        if tid and pchembl:
            if tid not in target_data:
                target_data[tid] = {
                    "target_pref_name": act.get("target_pref_name", ""),
                    "pchembl_values": [],
                }
            target_data[tid]["pchembl_values"].append(float(pchembl))

    profile = []
    for tid, info in target_data.items():
        vals = info["pchembl_values"]
        profile.append({
            "target_chembl_id": tid,
            "target_name": info["target_pref_name"],
            "n_measurements": len(vals),
            "best_pchembl": round(max(vals), 2),
            "mean_pchembl": round(np.mean(vals), 2),
        })

    df = pd.DataFrame(profile).sort_values("best_pchembl", ascending=False)
    print(f"Selectivity: {molecule_chembl_id} tested on {len(df)} targets")
    return df
```

## 6. 分子類似性 & サブ構造検索

```python
def similarity_search(smiles, threshold=70, max_results=25):
    """
    SMILES 構造による類似性検索。

    Parameters:
        smiles: str — クエリ SMILES
        threshold: int — Tanimoto 類似性閾値 (%)
        max_results: int — 最大結果数

    ToolUniverse:
        ChEMBL_search_similar_molecules(
            query=smiles, similarity_threshold=threshold
        )
        ChEMBL_search_substructure(smiles=smiles)
    """
    url = f"{CHEMBL_API}/similarity/{smiles}/{threshold}.json"
    params = {"limit": max_results}
    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    molecules = resp.json().get("molecules", [])

    rows = []
    for mol in molecules:
        rows.append({
            "molecule_chembl_id": mol.get("molecule_chembl_id"),
            "pref_name": mol.get("pref_name"),
            "similarity": mol.get("similarity"),
            "canonical_smiles": mol.get("molecule_structures", {}).get(
                "canonical_smiles", ""
            ),
            "max_phase": mol.get("max_phase"),
        })

    df = pd.DataFrame(rows)
    print(f"Similar molecules (>{threshold}%): {len(df)}")
    return df
```

## 7. ATC 分類 & 承認薬検索

```python
def search_approved_drugs(target_chembl_id, limit=50):
    """
    ターゲットに対する承認薬を ATC 分類とともに検索。

    ToolUniverse:
        ChEMBL_search_drugs(max_phase=4)
        ChEMBL_search_mechanisms(target_chembl_id=target_chembl_id)
        ChEMBL_search_atc_classification()
    """
    url = f"{CHEMBL_API}/mechanism.json"
    params = {
        "target_chembl_id": target_chembl_id,
        "limit": limit,
    }
    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    mechanisms = resp.json().get("mechanisms", [])

    rows = []
    for mech in mechanisms:
        drug_id = mech.get("molecule_chembl_id")
        drug_resp = requests.get(
            f"{CHEMBL_API}/molecule/{drug_id}.json", headers=HEADERS
        )
        if drug_resp.ok:
            drug = drug_resp.json()
            rows.append({
                "molecule_chembl_id": drug_id,
                "pref_name": drug.get("pref_name"),
                "max_phase": drug.get("max_phase"),
                "mechanism": mech.get("mechanism_of_action"),
                "action_type": mech.get("action_type"),
                "first_approval": drug.get("first_approval"),
                "atc_classifications": drug.get("atc_classifications", []),
            })

    df = pd.DataFrame(rows)
    print(f"Drugs/mechanisms for {target_chembl_id}: {len(df)}")
    return df
```

## 8. 構造アラート検出

```python
def check_structural_alerts(molecule_chembl_id):
    """
    化合物の構造アラート (PAINS, Dundee) を検出。

    ToolUniverse:
        ChEMBL_search_compound_structural_alerts(
            molecule_chembl_id=molecule_chembl_id
        )
    """
    url = f"{CHEMBL_API}/compound_structural_alert.json"
    params = {"molecule_chembl_id": molecule_chembl_id, "limit": 100}
    resp = requests.get(url, params=params, headers=HEADERS)
    resp.raise_for_status()
    alerts = resp.json().get("compound_structural_alerts", [])

    rows = []
    for a in alerts:
        rows.append({
            "alert_set_name": a.get("alert", {}).get("alert_set", {}).get(
                "set_name", ""
            ),
            "smarts": a.get("alert", {}).get("smarts", ""),
            "alert_name": a.get("alert", {}).get("alert_name", ""),
        })

    df = pd.DataFrame(rows)
    if len(df) > 0:
        print(f"⚠ {molecule_chembl_id}: {len(df)} structural alerts found")
    else:
        print(f"✓ {molecule_chembl_id}: No structural alerts")
    return df
```

## 9. 統合 SAR マイニングパイプライン

```python
def chembl_sar_pipeline(target_query, organism="Homo sapiens",
                        standard_type="IC50", max_nm=10000):
    """
    ChEMBL SAR マイニング統合パイプライン。

    Pipeline:
        search_target → get_target_activities → sar_analysis →
        selectivity_profile (top hits) → check_structural_alerts

    Parameters:
        target_query: str — ターゲット名
        organism: str — 生物種
        standard_type: str — 活性値タイプ
        max_nm: float — nM 閾値
    """
    # Step 1: ターゲット検索
    targets = search_target(target_query, organism)
    if targets.empty:
        print(f"No targets found for '{target_query}'")
        return None

    target_id = targets.iloc[0]["target_chembl_id"]
    print(f"\nSelected target: {target_id}")

    # Step 2: バイオアクティビティ取得
    activities = get_target_activities(target_id, standard_type, max_nm)
    if activities.empty:
        print("No activities found")
        return None

    # Step 3: SAR 解析
    sar = sar_analysis(activities)

    # Step 4: トップ化合物の選択性
    top = activities.nlargest(3, "pchembl_value")
    selectivity_results = []
    for _, row in top.iterrows():
        mol_id = row["molecule_chembl_id"]
        sel = selectivity_profile(mol_id)
        selectivity_results.append({"molecule": mol_id, "profile": sel})

    # Step 5: 構造アラートチェック
    alert_results = {}
    for _, row in top.iterrows():
        mol_id = row["molecule_chembl_id"]
        alerts = check_structural_alerts(mol_id)
        alert_results[mol_id] = len(alerts)

    result = {
        "target": target_id,
        "sar_summary": sar,
        "top_compounds": top.to_dict("records"),
        "structural_alerts": alert_results,
    }

    print(f"\n=== ChEMBL SAR Pipeline Complete ===")
    print(f"Target: {target_id}")
    print(f"Compounds: {sar['n_unique_molecules']}")
    print(f"Top hit pChEMBL: {sar['pchembl_range'][1]}")

    return result
```

---

## パイプライン統合

```
drug-target-profiling → chembl-assay-mining → admet-pharmacokinetics
  (ターゲット候補)       (SAR データマイニング)   (ADMET/PK 評価)
        │                       │                      ↓
compound-screening       selectivity_profile    molecular-docking
  (ZINC ライブラリ)       (マルチターゲット)       (Vina/DiffDock)
```

## パイプライン出力

| ファイル | 説明 | 次スキル |
|---------|------|---------|
| `results/chembl_activities.csv` | バイオアクティビティデータ | → admet-pharmacokinetics |
| `results/sar_summary.json` | SAR 統計サマリ | → drug-target-profiling |
| `results/selectivity_profile.csv` | 選択性プロファイル | → compound-screening |
| `results/structural_alerts.json` | 構造アラート結果 | → molecular-docking |

## 利用可能ツール (ToolUniverse SMCP)

| ツール名 | 用途 |
|---------|------|
| `ChEMBL_search_targets` | ターゲット検索 |
| `ChEMBL_get_target` | ターゲット詳細 |
| `ChEMBL_search_assays` | アッセイ検索 |
| `ChEMBL_get_assay` | アッセイ詳細 |
| `ChEMBL_get_assay_activities` | アッセイ活性データ |
| `ChEMBL_search_activities` | 活性検索 (フィルタ付き) |
| `ChEMBL_get_activity` | 活性データ詳細 |
| `ChEMBL_get_molecule` | 化合物詳細 |
| `ChEMBL_get_molecule_targets` | 化合物-ターゲットマップ |
| `ChEMBL_search_similar_molecules` | 類似性検索 |
| `ChEMBL_search_substructure` | サブ構造検索 |
| `ChEMBL_search_drugs` | 承認薬検索 |
| `ChEMBL_search_mechanisms` | 作用機序検索 |
| `ChEMBL_search_atc_classification` | ATC 分類 |
| `ChEMBL_search_compound_structural_alerts` | 構造アラート |
| `ChEMBL_search_cell_lines` | 細胞株検索 |
| `ChEMBL_search_binding_sites` | 結合サイト |
| `ChEMBL_get_drug` | 薬剤情報 |
| `ChEMBL_get_drug_mechanisms` | 薬剤作用機序 |
